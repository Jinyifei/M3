cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Program to determine dust temperature distribution using the
c     algorithm and equations of Draine & Li 2001, except use H for
c     enthalpy not U (which is grain potential here)
c
c     Output is IRphot in phot s-1 cm-3 Hz-1 Sr-1 (4pi)
c
c     Brent Groves 5/11/01
c
      subroutine dustinit ()
c
c Call once dust setup is known to save repeat intialisation
c
      include 'cblocks.inc'
c
      real*8 massatom(2)
      integer*4 dtype, k,l,inl
c
      massatom(1)=12.d0*amu
      massatom(2)=60.0855d0*amu
      do l=1,dustbinmax
        grainvol(l)=ftpi*(grainrad(l))**3
        d_rad(l)=gradedge(l+1)-gradedge(l)
        do dtype=1,numtypes
c    number of atoms in the grain
          atom_no(dtype,l)=grainvol(l)*graindens(dtype)/massatom(dtype)
        enddo
      enddo
c
c     setup energy bns for quick calculation
c
      do inl=1,infph-1
        dustsigt(inl)=0.d0
        do dtype=1,numtypes
          dustsigt(inl)=dustsigt(inl)+dcrosec(inl,dtype)
        enddo
      enddo
c
      dele=photev(infph)/photev(1)
      dele=dlog10(dele)/dble(dinfph-1)
      edphot(1)=photev(1)
      do k=1,dinfph-1
        edphot(k+1)=photev(1)*10.d0**(k*dele)
        dengy(k)=0.5d0*(edphot(k+1)+edphot(k))*ev
        dde(k)=(edphot(k+1)-edphot(k))*ev
        v_e(k)=dsqrt(2.d0*dengy(k)/me)
      enddo
c
      return
      end
      subroutine dusttemp (t_e, hdens, n_e, fi, dr, ircount)
      include 'cblocks.inc'
c
      real*8 e_kt,wid,den
      real*8 delta_t,dt,invkt
      real*8 tmin,tmax, dtmin, dtmax
      real*8 dtmin0, dtmax0
c
      real*8 hdens,grarea,fi,dr
      real*8 n_e,t_e, s_f, phi
      real*8 newflx, oldflx, flxratio
      real*8 flxtest1, flxtest2, flxabs
      real*8 teq, tlimit
      real*8 totaldust,flxabstot
c
      real*8 t_av(mxdust,mxdtype),norm1,norm2
      real*8 tnorm(mxdtype), tsignorm(mxdtype)
c
      integer*4 dtype,ircount,nmax,nmax0,nlimit
      integer*4 absmax,dabsmax
      integer*4 i,j,k,l,inl
c
      character tempfile*32
c      character fn*32
c      character pfx*16,sfx*4
c
      logical tfine, quickir
c
c     old
c       real*8 dcool(dinfph, MxTempbin)
c
       real*8 mbdist
c
c      functions
c
       real*8 planck, a, b
c
      planck(a,b)=(kbb*b*b*b*2)/(dexp(b*a)-1.d0)
c
c  Dust parameters
c
c
c      massatom(1)=12.d0*amu
c      massatom(2)=60.0855d0*amu
c      do l=1,dustbinmax
c       Grainvol(l)=ftpi*(grainrad(l))**3
c       d_rad(l)=gradedge(l+1)-gradedge(l)
c       do dtype=1,numtypes
c    number of atoms in the grain
c        atom_no(dtype,l)=Grainvol(l) * graindens(dtype)/massatom(dtype)
c       enddo
c      enddo
c
c  Determine IR mode
c
      quickir=.true.
      nlimit=50
c
      if (irmode.eq.1) then
        quickir=.true.
        nlimit=50
      else if (irmode.eq.2) then
        quickir=.true.
        nlimit=100
      else if (irmode.eq.3) then
        quickir=.false.
      else
        write (*,*) 'IR program incorrect'
        return
      endif
c
c      if (IRtemp) then
c         pfx='temp'
c         sfx='dat'
c         call newfile(pfx,4,sfx,3,fn)
c         tempfile=fn(1:12)
c         open(15,file=tempfile,status='unknown',access='APPEND')
c      endif
c
c  Initialise Temperature bins
c
      tmin=1.d0
      tmax=1201.d0
c
c   force initial temp/grid calcs.
c
      dtmin0=0.0d0
      dtmax0=0.0d0
c
c     setup energy bns for quick calculation
c
      flxtest1=0.d0
      flxabstot=0.d0
      do inl=1,infph-1
        wid=(photev(inl+1)-photev(inl))*evplk
        phot(inl)=4.d0*pi*dustphot(inl)*wid
c  convert from Jnu (1/4pi) to Fnu.dnu (erg s-1 cm-2)
        engy(inl)=cphotev(inl)
        if (engy(inl).gt.bgrain) then
          phot(inl)=phot(inl)*(1.d0-yinf*(1.d0-bgrain/engy(inl)))
c  losses due to PE effect
        endif
        photabs(inl)=phot(inl)
        flxtest1=flxtest1+phot(inl)
        engy(inl)=engy(inl)*ev
      enddo
c      write(*,*)  ' dusttemp flux:',flxtest1
c
c      delE = photev(infph)/photev(1)
c      delE = dlog10(delE)/dble(dinfph-1)
c      edphot(1) = photev(1)
c      do k=1,dinfph-1
c        edphot(k+1) = photev(1)*10**(k*delE)
c        dengy(k)=0.5d0*(edphot(k+1)+edphot(k))*eV
c        ddE(k)=(edphot(k+1)-edphot(k))*eV
c      enddo
c
c
c  set_up collisional heating bins
c
      do k=1,dinfph-1
c       v_e=dsqrt(2.d0*dengy(k)/me)
        collheat(k)=mbdist(t_e,dengy(k))*v_e(k)*dengy(k)*dde(k)
      enddo
c
      if (t_e.ne.0.0) then
        e_kt=ev/(rkb*t_e)
      else
        e_kt=0.d0
      endif
c
c  Initialise Region Flux output
c
      do k=1,infph-1
        irflux(k)=0.d0
      enddo
c
c  Do for both Silicates and Graphite grains
c   Type=1=Silicates, Type=2=Graphite
c
      totaldust=0.d0
      do 60 dtype=1,numtypes
c
c  Set-up Temp range
        dtmax=tmax
c
c  determine temperature for grain sizes considered
c
        do l=mindust(dtype),maxdust(dtype)
          totaldust=totaldust+dustsig(l,dtype)*d_rad(l)
c
c  Calculate sticking factor
c
          phi=grainpot(l,dtype)*e_kt
          if (phi.ge.0.d0) then
            s_f=0.5d0*(1.d0+phi)
          else
            s_f=0.5d0*dexp(phi)
          endif
c
c Areas in sq. microns within Transmatrix
c
          grarea=dustsig(l,dtype)*d_rad(l)
          do inl=1,infph-1
            cabs(inl)=grarea*absorp(inl,l,dtype)
            if (cabs(inl).ne.0.d0) absmax=inl
          enddo
c
c     Set up absphot bins
c
          flxabs=0.d0
          j=1
          do k=1,dinfph-1
            dabsphot(k)=0.d0
   10       den=cphotev(j)
            if (dustsigt(j).eq.0.d0) goto 20
            if (den.lt.edphot(k+1)) then
c           write(*,*) Cabs(j),dustsigt(j),photabs(j)
c
c  determine energy in scaled dust bin
c  using fraction of absorption due to grain l (ie sigma_l/sigma_tot)
c
              dabsphot(k)=dabsphot(k)+(cabs(j)/dustsigt(j))*(photabs(j)/
     &         (dr*fi*hdens))
              flxabs=flxabs+(cabs(j)/dustsigt(j))*(photabs(j)/(dr*fi*
     &         hdens))
              if (dabsphot(k).ne.0.d0) dabsmax=k
              j=j+1
              if (j.gt.(infph-1)) goto 20
              goto 10
            endif
          enddo
   20     continue
          flxabstot=flxabstot+flxabs
          teq=(flxabs/(grarea*stefan))**(0.25)
c          if (IRtemp) then
c 105              format(i4,2(1pg15.7,1x))
c           write(15,105) dtype,grainrad(l),Teq
c          endif
c         write(*,2000) dtype,l,flxabs,Teq
c 2000   format(2(i4,1x),'flxabs= ',1pg12.4,'Teq= ',1pg12.4)
c
c  setup initial Temperature & flxmeasures
c     (maxmimum dTmax in smallest (first) grain)
c
          nmax=50
          nmax0=nmax
          dtmin=tmin
          oldflx=1.d2
          newflx=1.d0
          flxratio=1.d2
          tfine=.true.
c
c  set P(T) Tlimit
c
c        if (grainrad(l).lt.1d-5) then
          tlimit=1.0d-15
c        else
c          Tlimit=1.0d-10
c        endif
          do while (tfine)
c
c     Clear Flux distribution
c
            do k=1,infph-1
              flxdist(k,l)=0.d0
            enddo
c
c recompute Temp grids if either dTmax or dTmin change
c
            if ((dtmin.ne.dtmin0).or.(dtmax.ne.dtmax0)
     &       .or.(nmax.ne.nmax0)) then
c
c  setup Temp grid
c
              delta_t=(dtmax-dtmin)/nmax
              dt=delta_t*0.5d0
c          write(*,*) 'tmax tmin delt dt', dTmax, dTmin, delta_T, dT
              t_edge(1)=dtmin
              do i=1,nmax
                t_edge(i+1)=dtmin+i*delta_t
                t_grid(i)=t_edge(i+1)-dt
              enddo
c
c  Setup grid for Plank spectrum at all dust temp. and photon energies
c
c           open(69,file='dcool.dat')
              do i=1,nmax
                invkt=1.d0/(rkb*t_grid(i))
                do k=1,infph-1
                  bbemiss(k,i)=planck(invkt,engy(k))
                enddo
              enddo
c
c record Ts for most recent BB and grid calcs
c
              dtmin0=dtmin
              dtmax0=dtmax
c
            endif
c
            call initgrids (dtype, grainvol(l), atom_no(dtype,l),
     &       t_grid, t_edge, h_grid, hmin, hmax, deltah, nmax,
     &       mxtempbin)
            do j=1,nmax
              do i=1,nmax
                invh(i,j)=0.d0
                if (i.ne.j) invh(i,j)=1.d0/(h_grid(i)-h_grid(j))
              enddo
            enddo
c
c     solve for transition matrix
c
            if (dabsmax.ne.1) then
              call transmatrix (grarea, n_e, s_f, nmax, absmax, dabsmax,
     &          t_e)
c  With transition matrix obtained now solve for P vector
              call probsolve (tr_matrix, t_prob, nmax, mxtempbin)
            endif
c
c     Calculate Flux distribution and total flux
c
            newflx=0.d0
            do j=1,nmax
              if (t_prob(j).ne.0.d0) then
                do k=1,infph-1
                  if (bbemiss(k,j).lt.1.d-50) goto 30
                  flxdist(k,l)=flxdist(k,l)+t_prob(j)*bbemiss(k,j)*
     &             absorp(k,l,dtype)*grarea
                enddo
   30           continue
              endif
            enddo
            do k=1,infph-1
              wid=(photev(k+1)-photev(k))*ev
              newflx=newflx+flxdist(k,l)*4*pi*wid
            enddo
c
c
            if (irtemp.eq.1) then
              t_av(l,dtype)=0.d0
              do j=1,nmax
                t_av(l,dtype)=t_av(l,dtype)+t_prob(j)*t_grid(j)
              enddo
            endif
c
c     QuickIR using only 50 bins for energy conservation
c
            if (quickir) then
              k=nmax
              do while (t_prob(k).lt.tlimit)
                dtmax=t_edge(k)
                k=k-1
              enddo
              j=1
              do while (t_prob(j).lt.tlimit)
                dtmin=t_edge(j+1)
                j=j+1
              enddo
c  In case hotter than expected
              if ((k.eq.nmax).and.(dtmax.eq.tmax)) dtmax=tmax+500
              if ((j.lt.4).and.((nmax-k).lt.4)) then
                nmax=nmax*2
              endif
              if (nmax.gt.nlimit) tfine=.false.
c
c  Not quick do full IR for temp distributions
c
            else
              if (nmax.lt.200) then
                k=nmax
                do while (t_prob(k).lt.tlimit)
                  dtmax=t_edge(k)
                  k=k-1
                enddo
                j=1
                do while (t_prob(j).lt.tlimit)
                  dtmin=t_edge(j+1)
                  j=j+1
                enddo
c  In case hotter than expected
                if ((k.eq.nmax).and.(dtmax.eq.tmax)) dtmax=tmax+500
                if ((j.lt.4).and.((nmax-k).lt.4)) nmax=nmax*2
              else
                k=nmax
                do while (t_prob(k).lt.1.d-15)
                  dtmax=t_edge(k)
                  k=k-1
                enddo
                j=1
                do while (t_prob(j).lt.1.d-15)
                  dtmin=t_edge(j+1)
                  j=j+1
                enddo
                if ((j.lt.(nmax/20)).and.((nmax-k).lt.(nmax/20))) nmax=
     &           nmax*2
              endif
              if (nmax.gt.400) tfine=.false.
            endif
          enddo
c
c       flxtest1=flxtest1+flxabs
c
c  Absorbed & Emitted per Grain
c
c       flxtest2=0.d0
c       do k=1,infph-1
c          wid=(photev(k+1)-photev(k))*ev
c          flxtest2=flxtest2+flxdist(k,l)*wid
c       enddo
c       flxtest2=flxtest2*4*pi*d_rad(l)
c       write(*,*) grainrad(l),'E IRFlux Diff =',
c     &    (1.d0 - flxtest2/flxabs)
c
c
        enddo
c
c  Finally determine radiation field
c   (use trapezoidal integration over dust sizes)
c
        if (mindust(dtype).ne.maxdust(dtype)) then
          do l=mindust(dtype),maxdust(dtype)
            do k=1,infph-1
              if (flxdist(k,l).lt.1.0d-50) goto 40
c           IRFlux(k)=IRFlux(k)+0.5d0*d_rad(l)*(flxdist(k,l-1)
c                  +flxdist(k,l))
              irflux(k)=irflux(k)+flxdist(k,l)
            enddo
   40       continue
          enddo
        else
          l=mindust(dtype)
          do k=1,infph-1
            if (flxdist(k,l).lt.1.d-50) goto 50
            irflux(k)=irflux(k)+flxdist(k,l)
          enddo
   50     continue
        endif
   60 continue
c
c  Detemine IRphot (photons/s/cm^3/Hz)
c     flux in ergs s-1 Hatom-1 Hz-1 Sr-1 (due to dustsig)
c
c
      flxtest2=0.d0
      do k=1,infph-1
        wid=(photev(k+1)-photev(k))*ev
        flxtest2=flxtest2+irflux(k)*4*pi*wid
        irphot(k)=irflux(k)*hdens*plk/engy(k)
      enddo
      write (*,70) dr*fi*hdens*flxabstot,flxtest1,dr*fi*hdens*flxtest2
   70 format('IRin1=',1pg12.5,' IRin2=',1pg12.5,' IRout=',1pg12.5)
      teq=(flxtest2/(totaldust*stefan))**(0.25)
      write (*,80) teq
   80 format('Average Dust Temp: ',1pg14.7)
      if (irtemp.eq.1) then
        do dtype=1,numtypes
          norm1=0.d0
          norm2=0.d0
          tnorm(dtype)=0.d0
          tsignorm(dtype)=0.d0
          do l=mindust(dtype),maxdust(dtype)
            norm1=norm1+dustnum(l,dtype)*d_rad(l)
            norm2=norm2+dustsig(l,dtype)*d_rad(l)
          enddo
          do l=mindust(dtype),maxdust(dtype)
            tnorm(dtype)=tnorm(dtype)+t_av(l,dtype)*dustnum(l,dtype)*
     &       d_rad(l)
            tsignorm(dtype)=tsignorm(dtype)+t_av(l,dtype)*dustsig(l,
     &       dtype)*d_rad(l)
          enddo
          tnorm(dtype)=tnorm(dtype)/norm1
          tsignorm(dtype)=tsignorm(dtype)/norm2
        enddo
        tempfile='dusttemp.dat'
        open (15,file=tempfile,status='unknown',access='APPEND')
   90      format(i4,2x,2(1pg15.7,2x),2x,2(1pg15.7,2x))
        write (15,90) ircount,tnorm(1),tsignorm(1),tnorm(2),tsignorm(2)
        close (15)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine intialises Grids and grid values
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initgrids (dtype, v, atom_no, t_grid, t_edge, h_grid,
     &hmin, hmax, deltah, nmax, mxbin)
      implicit none
      integer*4 dtype,nmax, i, mxbin
      real*8 v,atom_no,t_grid(mxbin),t_edge(mxbin),
     & h_grid(mxbin),hmin(mxbin),hmax(mxbin),deltah(mxbin)
c     external functions
      real*8 sil_enth,gra_enth
c  calculate enthalpy for each temperature
      if (dtype.eq.1) then
        do i=1,nmax
          h_grid(i)=sil_enth(atom_no,v,t_grid(i))
          hmin(i)=sil_enth(atom_no,v,t_edge(i))
          hmax(i)=sil_enth(atom_no,v,t_edge(i+1))
          deltah(i)=hmax(i)-hmin(i)
        enddo
      else
        do i=1,nmax
          h_grid(i)=gra_enth(atom_no,t_grid(i))
          hmin(i)=gra_enth(atom_no,t_edge(i))
          hmax(i)=gra_enth(atom_no,t_edge(i+1))
          deltah(i)=hmax(i)-hmin(i)
        enddo
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Function integrates Silicate Heat Capacity to get Enthalpy (H-n)
c Heat capacity from experimental values (see G. Mark Voit 1991)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function sil_enth(n,v,t)
      implicit none
      real*8 n,v,t
      real*8 f, fint1, fint2, fint3, fint4
      real*8 c1, c2, c3, t2
      f=(1.d0-(2.d0/n))*v
      c1=5.83333333333333333d+7
      c2=8.76008564039770429d+8
      c3=8.54862846766113158d+9
      fint1=0.d0
      fint2=0.d0
      fint3=0.d0
      fint4=0.d0
      if (t.le.50.d0) then
        fint1=1.40d3*(t*t*t)/3.d0
      else if (t.le.150.d0) then
        fint1=c1
        t2=50.d0
        fint2=2.1647d4*((t**2.3d0)-(t2**2.3d0))/2.3d0
      else if (t.le.500.d0) then
        fint1=c1
        fint2=c2
        t2=150.d0
        fint3=4.8369d5*((t**1.68d0)-(t2**1.68d0))/1.68d0
      else
        fint1=c1
        fint2=c2
        fint3=c3
        t2=500.d0
        fint4=3.3103d7*(t-t2)
      endif
      sil_enth=f*(fint1+fint2+fint3+fint4)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Function integrates Graphite Heat Capacity to get Enthalpy (H-n)
c Heat capacity fitted (to within 3%) to experimental values (se G&D)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function gra_enth(n,t)
      implicit none
      real*8 n,t
      real*8 hatom
      hatom=(4.15d-22*t**3.3d0)/
     &        (1+6.51d-3*t+1.5d-6*t**2+8.3d-7*t**2.3d0)
      gra_enth=(1.d0-2.d0/n)*n*hatom
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Function evaluates Planck spectrum given Temp and photon energy
c  ergs s-1 cm-2 Hz-1 Sr-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      real*8 function Planck(invkT,E)
c invkT = 1/(kT)
c
c      include 'const.inc'
c      real*8 invkT, E, a, b
c
cc      Planck=(2.d0*E**3.d0/((cls**2)*(plk**3)))/(dexp(E/(rkb*T))-1.d0)
c      a = dexp(E*invkT)-1.d0
c      b = kbb*E*E*E*2
c      Planck=b/a
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Function evaluates Maxwell-Boltzmann distribution given
c  electron Temp and energy
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function mbdist(t,e)
      include 'const.inc'
      real*8 t,e, a, b, c, invpi
      parameter (invpi = 0.318309886183791d0)
c      MBdist=2.d0*dsqrt(E/(pi*(rkb*T)**3.d0))*dexp(-E/(rkb*T))
      a=1.d0/(rkb*t)
      b=dsqrt(e*(a*a*a)*invpi)
      c=dexp(-e*a)
      mbdist=2*b*c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine Makes Transisition matrix
c  |1 2 0 0 ... 0|
c  |3 1 2 0 ... 0|
c  |.   .       :|
c  |:     .     :|
c  |: ....3 1 2 0|
c  |3 ......3 1 2|
c  |3 ......3 3 1|
c  2 is Thermal continuous cooling approx. only from level above
c  3 is heating as in Draine & Li
c  1 is the sum of the column representing losses from that level
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine transmatrix (grarea, n_e, s_f, nmax, absmax, dabsmax,
     &t_e)
      include 'cblocks.inc'
      real*8 grarea
      real*8 n_e,s_f,stick,t_e
      integer*4 nmax
      integer*4 absmax,dabsmax
      real*8 den,wid
      real*8 cool,heat
      real*8 w1,w2,w3,w4,w5
      real*8 gfac,invfac
      real*8 dht1, invdht1
      real*8 hmn, hmx
      integer*4 i,j,k
c
c  Initialise
c
      do j=1,nmax
        do i=1,nmax
          tr_matrix(i,j)=0.d0
        enddo
      enddo
c
c  Cooling (#2)
c  Uses midpoint integration with Qabs energy points
c  Will improve later
c
      do i=2,nmax
        cool=0.d0
        do k=1,absmax
          den=cphote(k)
          wid=(photev(k+1)-photev(k))*ev
          if (den.gt.h_grid(i)) goto 10
          cool=cool+bbemiss(k,i)*cabs(k)*wid
        enddo
   10   tr_matrix(i-1,i)=(4.d0*pi)*invh(i,i-1)*cool
      enddo
      stick=n_e*grarea*s_f*(segrain*dexp(-t_e/2.0d5))
c      stick=0.d0
      do i=2,nmax-1
c
c  Heating (#3) from l->u (j->i)
c  Allows for finite width of H bins
c  Uses midpoint integration with Qabs energy points
c
        hmn=hmin(i)
        hmx=hmax(i)
        do j=1,i-1
          w1=hmn-hmax(j)
          w2=min((hmn-hmin(j)),(hmx-hmax(j)))
          w3=max((hmn-hmin(j)),(hmx-hmax(j)))
          w4=hmx-hmin(j)
          invfac=1.d0/(deltah(i)*deltah(j))
c
c  Photon & Collisional (currently turned off) heating
c
          heat=0.d0
          do k=1,dabsmax
            if (dengy(k).ge.w1) then
              if (dengy(k).lt.w2) then
                gfac=(dengy(k)-w1)*invfac
              else if (dengy(k).lt.w3) then
                gfac=min(deltah(i),deltah(j))*invfac
              else if (dengy(k).le.w4) then
                gfac=(w4-dengy(k))*invfac
              else
                goto 20
              endif
              heat=heat+gfac*(dabsphot(k)+stick*collheat(k))
            endif
          enddo
   20     continue
          tr_matrix(i,j)=deltah(i)*invh(i,j)*heat
        enddo
c  Allow for intrabin absorbtions
        heat=0.d0
        dht1=deltah(i-1)
        invdht1=1.d0/dht1
        do k=1,dabsmax
          if (dengy(k).gt.dht1) goto 30
          heat=heat+(1.d0-(dengy(k)*invdht1))*dabsphot(k)
        enddo
   30   continue
        tr_matrix(i,i-1)=tr_matrix(i,i-1)+invh(i,i-1)*heat
      enddo
      hmn=hmin(nmax)
      do j=1,nmax-1
        w1=hmn-hmax(j)
        w4=hmn-hmin(j)
        w5=1.d0/(w4-w1)
        heat=0.d0
c  Photon & collisional heating
        do k=1,dabsmax
          if ((dengy(k).gt.w1).and.(dengy(k).lt.w4)) then
            heat=heat+(dengy(k)-w1)*w5*(dabsphot(k)+stick*collheat(k))
          else if (dengy(k).ge.w4) then
            heat=heat+dabsphot(k)+stick*collheat(k)
          endif
        enddo
        tr_matrix(nmax,j)=1.d0*invh(nmax,j)*heat
      enddo
      heat=0.d0
      dht1=deltah(nmax-1)
      invdht1=1.d0/dht1
      do k=1,dabsmax
        if (dengy(k).le.dht1) then
          heat=heat+(1.d0-dengy(k)*invdht1)*dabsphot(k)
        else
          goto 40
        endif
      enddo
   40 continue
      tr_matrix(nmax,nmax-1)=tr_matrix(nmax,nmax-1)+invh(nmax,nmax-1)*
     &heat
c
c  Losses from i (Heating and Cooling) (#1)
c
c      open(69,file='Trmatrix.dat')
c      do i=1,Nmax
c       sum=0.d0
c       do j=1,Nmax
c         sum=sum+Tr_matrix(j,i)
c       enddo
c       Tr_matrix(i,i)=-sum
c       write(69,169) (Tr_matrix(j,i),j=1,Nmax)
c      enddo
c 169  format(400(1pg11.4,1x))
c      write(69,*)
c      close(69)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This subroutine solves for probability distribution given
c    Transition matrix Tr_matrix(i,j) (see G & D)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine probsolve (tr_matrix, t_prob, nmax, mxbin)
      implicit none
      integer*4 nmax,submax,mxbin
      parameter (submax=401)
      real*8 tr_matrix(mxbin,mxbin),t_prob(mxbin)
      real*8 bm(submax,submax), xo(submax)
      real*8 sum
      integer*4 i,j
      do i=1,nmax
        bm(nmax,i)=tr_matrix(nmax,i)
        do j=nmax-1,i,-1
          bm(j,i)=tr_matrix(j,i)+bm(j+1,i)
        enddo
      enddo
      xo(1)=1.d0
      do j=2,nmax
        sum=0.d0
        do i=1,j-1
          sum=sum+bm(j,i)*xo(i)
        enddo
        if (tr_matrix(j-1,j).ne.0.d0) then
          xo(j)=sum/tr_matrix(j-1,j)
        else
          print *,'T(j-1,j)=0. IN trouble'
          stop
        endif
      enddo
      sum=0.d0
      do j=1,nmax
        sum=sum+xo(j)
      enddo
      if (sum.eq.0.d0) then
        print *,'sum of X=0. IN trouble'
        stop
      endif
c solve for T_prob
      do j=1,nmax
        t_prob(j)=xo(j)/sum
        if (t_prob(j).lt.1.0d-20) t_prob(j)=0.d0
      enddo
      return
      end
