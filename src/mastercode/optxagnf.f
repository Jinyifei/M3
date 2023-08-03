cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Special FORTRAN MAPPINGS photsou version of
c the xspec 12 optxagnf model routine.
c Includes all optxagnf dependicies
c
c See:
c Done et al Mon. Not. R. Astron. Soc. 420, 1848–1860 (2012)
c
c ear(i), in keV , is the energy at the right of the bin:
c ear(i)-ear(i-1), from  i = 1 .. ne, so index 0 is added for
c first interval.  Strictly I think bin centre may be better
c than ear(i) for getting luminosity, but as most spectra are
c highly non-linear, even bin center is not correctly weighted...
c but if ec = 0.5d0*(ear(i)+ear(i-1))
c then approx:L(i) = photar(i) * ec to get luminosity of the bin,
c and sum L(i) to get total luminosity.
c
c optxagnf: creates integral photons/cm^2/s/bin
c This routine returns Inu = ergs/cm^2/s/Hz/pi, rebinned for
c current MAPPINGS photdat.
c
c Note, total luminosity may be a little less than expected
c if MAPPINGS PHOTDAT has a limited range and misses some
c of the optxagnf range.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine optxagnf(param, flux)
c
      include 'cblocks.inc'
c
      real*8 param(12)
      integer*4 i
      real*8 flux(mxinfph)
      real*8 phots(mxoptxbins)
      real*8 ec(mxoptxbins)
      real*8 optxhnu(mxoptxbins)
      real*8 dloge,ew,ewhz,ecergs
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Input paramters:
c
c     param(1) BH mass in solar
c     param(2) distance, default 1.5e18 cm
c     param(3) log mass accretion rate in L/LEdd
c     param(4) bhastar rotation : usually 0.0
c     param(5) rcorona/rg
c     param(6) log10 rout/rg
c     param(7) opt thick kte
c     param(8) opt thick tau
c     param(9) power law gamma
c     param(10) frac of coronal power in power law.
c     11,12 unused atm, 11 used to be redshift
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Now really simple :) RSS
c     always returns the total
c     no need for inefficient repeated -ve parameter calls
      do i = 1,10,1
      optxparam(i)=param(i)
      enddo
      optxparam(11) = 0.0d0
      optxparam(12) = 0.0d0
c
  100 format(/" OPTXAGNF Parameters: ",//
     &        "   M_BH: ",1pg12.4," M_0",/
     &        "   L   : ",1pg12.4," L/Ledd",/
     &        "   D   : ",1pg12.4," cm",/
     &        "   Rcor: ",1pg12.4," Rg      ",
     &        "   Rmax: ",1pg12.4," log[Rg] ",/
     &        "   C_T : ",1pg12.4," keV     ",
     &        "   Tau:  ",1pg12.4," Optical Depth ",/
     &        "   gam : ",1pg12.4," (-slope)",
     &        "   fpl : ",1pg12.4," Fraction ",/)
c
      write(*,100) optxparam( 1)
     & ,(10.d0**optxparam( 3))
     & ,optxparam( 2)
     & ,optxparam( 5), optxparam( 6)
     & ,optxparam( 7), optxparam( 8)
     & ,optxparam( 9), optxparam(10)
c
c set up energy bins, keV
c
      neoptx=5000
      optxear(0)=1.0d-4
      optxear(neoptx)=5.0d2
c
      dloge=dlog10(optxear(neoptx)/optxear(0))/dble(neoptx)
      do i=1,neoptx,1
        optxear  (i)=10.d0**(dlog10(optxear(0))+dloge*dble(i))
        optxtotal(i)=0.d0
        optxdisk (i)=0.d0
        optxcoro (i)=0.d0
        optxnont (i)=0.d0
        phots    (i)=0.d0
      enddo
c
c optxear(i), in keV , is the energy at the right of the bin:
c optxear(i)-optxear(i-1), from  i = 1 .. ne, so index 0 is added for
c first interval.  I think bin centre is better
c than optxear(i) for getting luminosity, but as most spectra are
c highly non-linear, even bin center is not correctly weighted...
c but if ec = 0.5d0*(optxear(i)+optxear(i-1))
c then approx:L(i) = optxtotal(i) * ec to get luminosity of the bin,
c and sum L(i) to get total luminosity.
c
c mydiskf Returns array optxtotal of bin integrals:
c         * as photons/cm^2/s per bin
c          not differential* phots/cm^2/s/keV etc etc
c that is : bin width is already integrated --like all XSPEC additive
c model components.
c
      call mydiskf (optxparam, phots)
c
c   convert photons/cm^2/s/bin to Inu at bin centres
c
      do i=1,neoptx,1
        ew=(optxear(i)-optxear(i-1))*1.0d3  ! eV
        ewhz=ew*ev/plk
        ec(i)=0.5d0*(optxear(i)+optxear(i-1))*1.0d3  ! eV
        ecergs=ec(i)*ev
c phots/cm2/s/bin.  *photon energy, / bin width / pi -> Hnu
        optxhnu(i)=phots(i)*ecergs/(pi*ewhz) ! bin centre mean Inu
      enddo
c
c     actually call with Inu in flux...
c     but rebin doesn't scale like readrebin does, so OK
c
      call rebin(ec,optxhnu,neoptx,flux)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     program to integrate the disk equations from shakura-sunyaev disk
c     as given by Novikov and Thorne
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mydiskf (param,phots)
      include 'cblocks.inc'
      real*8 param(12),phots(mxoptxbins)
c
      real*8 m,mdot,rsg
      real*8 bhastar,da,z1,z2,rms,fcor,fpl
      real*8 corona,cor,pow,rcor,mytemp,fcol
      real*8 rgcm,t,wt0,r,dr,dlogr,dflux
      real*8 el,ec,er,dfluxl,dfluxc,dfluxr
      real*8 dsk,tot,d,logrout
      integer*4 i,iin,imax,n,ne
      logical first
      real*8 flux(mxoptxbins)
      real*8 disc,discu,mdotedd,alpha,eff
      real*8 lfrac,hfrac,renorm,lonledd
      real*8 lpar(5),lphot(mxoptxbins),lphote(mxoptxbins)
      real*8 hpar(5),hphot(mxoptxbins),hphote(mxoptxbins)
c     system parameters
c
      ne=neoptx
c
      m=param(1)!in solar units
      d=param(2)!in cm now, not Mpc
      lonledd=10.d0**(param(3))
      mdotedd=10.d0**(param(3))
      bhastar=param(4)
      alpha=0.1d0
c
      rgcm=rgm*m
c     get rms
      z1=((1.d0-bhastar**2)**(1.d0/3.d0))
      z1=z1*(((1+bhastar)**(1.d0/3.d0))+((1-bhastar)**(1.d0/3.d0)))
      z1=1.d0+z1
      z2=dsqrt(3.d0*bhastar*bhastar+z1*z1)
      rms=3.d0+z2-dsqrt((3.d0-z1)*(3.d0+z1+2.0d0*z2))
      eff=1.0d0-dsqrt(1.0d0-2.0d0/(3.0d0*rms))
c
c      mdot in g/s Ledd=1.2572e38 m and c2=8.98755e20
      mdot=lonledd*ledd*m/(cls*cls*eff)
c     disc parameters
      rcor=dabs(param(5))!in rin/rg
      logrout=param(6)!log10rout/rg
c     iv -ve then calculate rout=rsg from laor & netzer 1989
      if (param(6).lt.0.0d0) then
        logrout=((m/1.d9)**(-2.d0/9.d0))*((mdotedd)**(4.d0/9.d0))
        logrout=2150.d0*logrout*(alpha**(2.d0/9.d0))
        logrout=dlog10(logrout)
      endif
      rsg=10.d0**logrout
      fpl=dabs(param(10))
c     initialise
      do n=1,ne,1
        optxtotal(n)=0.0d0
        optxdisk(n)=0.0d0
        optxcoro(n)=0.0d0
        optxnont(n)=0.0d0
        flux(n)=0.0d0
        phots(n)=0.d0
      enddo
      if (rcor.gt.rms) then
        iin=50
      else
        iin=0
        rcor=rms
      endif
      imax=1000
      first=.true.
      wt0=0.d0
      corona=0.0d0
      disc=0.0d0
      discu=0.0d0
      do i=1,iin+imax,1
        if (i.le.iin) then
          dlogr=dlog10(rcor/rms)/dble(iin-1)
          r=10.d0**(dlog10(rms)+dble(i-1)*dlogr+dlogr/2.d0)
        else
          dlogr=dlog10(rsg/rcor)/dble(imax-1)
          r=10.d0**(dlog10(rcor)+dble(i-iin-1)*dlogr+dlogr/2.d0)
        endif
        dr=10.d0**(dlog10(r)+dlogr/2.d0)-10.d0**(dlog10(r)-dlogr/2.d0)
c
        da=4.d0*pi*r*dr*rgcm*rgcm
c
        t=mytemp(m,bhastar,mdot,rms,r)
        if (t.gt.1.0d5) then
          fcol=(72.d0/(t/kkev))
          fcol=fcol**(1.d0/9.d0)
        else
          if (t.gt.3.0d4) then
c              linear from 1 to 2.7 between 4e4 and 1e5
            fcol=(t/3.0d4)**0.82d0
          else
            fcol=1.d0
          endif
        endif
        if (i.le.iin) then
c            discu=discu+4.d0*pi*r*dr*rgcm*rgcm*5.67d-5*(t**4)
          discu=discu+da*stefan*(t**4)
        else
          if (first) then
            wt0=t*fcol/kkev
c             write(*,*) 'seed photon',t,wt0*kkev
            first=.false.
          endif
        endif
c         disc=disc+4.d0*pi*r*dr*rgcm*rgcm*5.67d-5*(t**4)
        disc=disc+da*stefan*(t**4)
        t=t/kkev
c        go over each photon energy - 3rd order ppm integral average
        do n=1,ne,1
C          en=0.5d0*(dlog10(optxear(n))+dlog10(optxear(n-1)))
C          en=10.d0**en
          el = optxear(n-1)
          ec = 0.5d0*(optxear(n-1)+optxear(n))
          er = optxear(n)
c
c           do blackbody spectrum left edge
c
          if ((el.lt.30.d0*t*fcol).and.(r.gt.rcor)) then
c               dflux=pi*2.d0*plk*((en*kevhz)**3)/9.0d20
            dfluxl=pi*2.d0*plk*((el*kevhz)**3)/(cls*cls)
            dfluxl=dfluxl*da/(dexp(el/(t*fcol))-1.d0)
            dfluxl=dfluxl/(fcol**4)
          else
            dfluxl=0.0d0
          endif
c           do blackbody spectrum center
          if ((ec.lt.30.d0*t*fcol).and.(r.gt.rcor)) then
c               dflux=pi*2.d0*plk*((en*kevhz)**3)/9.0d20
            dfluxc=pi*2.d0*plk*((ec*kevhz)**3)/(cls*cls)
            dfluxc=dfluxc*da/(dexp(ec/(t*fcol))-1.d0)
            dfluxc=dfluxc/(fcol**4)
          else
            dfluxc=0.0d0
          endif
c           do blackbody spectrum right edge
          if ((er.lt.30.d0*t*fcol).and.(r.gt.rcor)) then
c               dflux=pi*2.d0*plk*((en*kevhz)**3)/9.0d20
            dfluxr=pi*2.d0*plk*((er*kevhz)**3)/(cls*cls)
            dfluxr=dfluxr*da/(dexp(er/(t*fcol))-1.d0)
            dfluxr=dfluxr/(fcol**4)
          else
            dfluxr=0.0d0
          endif
          dflux=(dfluxl+4.d0*dfluxc+dfluxr)/6.d0
          flux(n)=flux(n)+(dflux)
        enddo
      enddo
c
      fcor=discu/(disc-discu)
c      write(*,*) 'fcor', fcor, discu, disc
c
c      dsk=0.d0
c
      do n=1,ne,1
c         dsk=dsk+flux(n)*(optxear(n)-optxear(n-1))*kevhz
c        this is ergs cm^2 s-1 Hz^-1
        ec = 0.5d0*(optxear(n-1)+optxear(n))
        flux(n)=flux(n)/(4.d0*pi*d*d)
c        photons is flux/hv - photons cm^2 s-1 Hz^-1
        flux(n)=flux(n)/(plk*kevhz*ec)
c        now multiply by energy band in Hz
        optxdisk(n)=flux(n)*(optxear(n)-optxear(n-1))*kevhz
      enddo
c
c      write(*,*) 'd',dsk, dsk*fcor, dsk*(1.d0+fcor)
c
c     now add low temcomptonised emission with comptt
      lpar(1)=0.d0
      lpar(2)=wt0
      lpar(3)=param(7)
      lpar(4)=dabs(param(8))!if-veplotcomp
      lpar(5)=1
      call xstitg (optxear, neoptx, lpar, lphot, lphote)
c
c     now add high temp compton emission
      hpar(1)=param(9)
      hpar(2)=100.d0
      hpar(3)=wt0
      hpar(4)=0.d0
      hpar(5)=0.d0
      call donthcomp (optxear, neoptx, hpar, hphot, hphote)
c
      cor=0.d0
      pow=0.d0
      dsk=0.d0
      do n=1,ne,1
        ec =0.5d0*(optxear(n-1)+optxear(n))
        cor=cor+lphot(n)*ec
        pow=pow+hphot(n)*ec
        dsk=dsk+optxdisk(n)*ec
      enddo
c
c      write(*,*) 't2',dsk*kev*(4.d0*pi*d*d)
c
c     now add in the rest of the components
c
      lfrac=(dsk/cor)*fcor*(1.d0-fpl)
      hfrac=(dsk/pow)*fcor*fpl
c      write(*,*) lfrac,hfrac
      do n=1,ne,1
        optxcoro(n)=lphot(n)*lfrac
        optxnont(n)=hphot(n)*hfrac
      enddo
c
c renormalise to total spectrum to recover L/Ledd
c
      tot=0.d0
      do n=1,ne,1
        optxtotal(n)=(optxdisk(n)+optxcoro(n)+optxnont(n))
        ec=0.5d0*(optxear(n-1)+optxear(n))
        tot=tot+optxtotal(n)*ec
      enddo
      renorm=(lonledd*ledd*m)/(tot*(4.d0*pi*d*d)*kev)
c
c now safe to rescale even partial spectra
c
       do n=1,ne,1
         optxtotal(n)= optxtotal(n)*renorm
         optxdisk(n) = optxdisk(n)*renorm
         optxcoro(n) = optxcoro(n)*renorm
         optxnont(n) = optxnont(n)*renorm
         phots(n)    = optxtotal(n)
       enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find temperature as a function of mass, spin and r in units of rg
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function mytemp(m0,bhastar,mdot0,rms,r0)
      include 'const.inc'
      real*8 m0,bhastar,mdot0,r,r0,rgcm,y,yms
      real*8 y1,y2,y3,part1,part2,part3,rms,a,b
c
      rgcm=rgm*m0
c
      r=r0
      y=dsqrt(r0)
      yms=dsqrt(rms)
      y1=2.d0*dcos((dacos(bhastar)-pi)/3.d0)
      y2=2.d0*dcos((dacos(bhastar)+pi)/3.d0)
      y3=-2.d0*dcos((dacos(bhastar)/3.d0))
c
      part3=3.d0*((y3-bhastar)**2)*dlog((y-y3)/(yms-y3))
      part3=part3/(y*y3*(y3-y1)*(y3-y2))
      part2=3.d0*((y2-bhastar)**2)*dlog((y-y2)/(yms-y2))
      part2=part2/(y*y2*(y2-y1)*(y2-y3))
      part1=3.d0*((y1-bhastar)**2)*dlog((y-y1)/(yms-y1))
      part1=part1/(y*y1*(y1-y2)*(y1-y3))
      a=1.d0-yms/y-(3.d0*bhastar/(2.d0*y))*dlog(y/yms)-part1-part2-part3
      b=1.d0-3.d0/r+2.d0*bhastar/(r**1.5d0)
c
      mytemp=3.d0*grav*msun*m0*mdot0
      mytemp=mytemp/(8.d0*pi*stefan*(r*rgcm)**3)
      mytemp=(mytemp*a/b)**0.25d0
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     driver for the Comptonization code solving Kompaneets equation
c     seed photons - (disc) blackbody
c     reflection + Fe line with smearing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine donthcomp (ear, ne, param, photar, photer)
      include 'const.inc'
      integer*4 ne
      real*8 param(*) , ear(0:ne) , photar(ne), photer(ne)
c
c     Version optimized for a number of data files but with the same values
c     of parameters:
c
c     number of model parameters: 16
c     1: photon spectral index
c     2: plasma temperature in keV
c     3: (disc)blackbody temperature in keV
c     4: type of seed spectrum (0-blackbody, 1-diskbb)
c     5: redshift
c
      real*8  prim(0:mxoptxbins)
c
      integer*4 n , ierr, np
      real*8 xninv , normfac , normlum, spp, pa0(5)
      integer*4 nth
      real*8 xth(900) , spt(900)
      logical recalc
      real*8    zfactor
      integer*4 i,j, jl
      save pa0,normfac,normlum,xth,nth,spt
      data pa0/5*9999.d0/

c this model does not calculate errors
      do i=1,ne
        photer(i)=0.d0
      enddo
c     xtot is the energy array (units m_e c^2)
c     spnth is the nonthermal spectrum alone (E F_E)
c     sptot is the total spectrum array (E F_E), = spref if no reflection
      ierr=0
      zfactor=1.d0+param(5)
c  calculate internal source spectrum if input parameters have changed
      np=5
      recalc=.false.
      do n=1,np
        if (param(n).ne.pa0(n)) recalc=.true.
      enddo
      if (recalc) then
        if (param(4).lt.0.5d0) then
          call thcompton (param(3)/511.d0, param(2)/511.d0, param(1),
     &     xth, nth, spt)
        else
          call thdscompton (param(3)/511.d0, param(2)/511.d0, param(1),
     &     xth, nth, spt)
        endif
        xninv=511.d0/zfactor
        normfac=1.d0/spp(xninv,xth,nth,spt)
c Calculate luminosity normalization  (used in another model!)
        normlum=0.d0
        do i=2,nth-1
          normlum=normlum+0.5d0*(spt(i)/xth(i)+spt(i-1)/xth(i-1))*
     &     (xth(i)-xth(i-1))
        enddo
        normlum=normlum*normfac
      endif
c     zero arrays
      do i=1,ne
        photar(i)=0.0d0
        prim(i)=0.0d0
      enddo
      prim(0)=0.0d0
c     put primary into final array only if scale >= 0.
      j=1
      do i=0,ne
        do while (j.le.nth.and.511.d0*xth(j).lt.ear(i)*zfactor)
          j=j+1
        enddo
        if (j.le.nth) then
          if (j.gt.1) then
            jl=j-1
            prim(i)=spt(jl)+(ear(i)/511.d0*zfactor-xth(jl))*(spt(jl+1)-
     &       spt(jl))/(xth(jl+1)-xth(jl))
          else
            prim(i)=spt(1)
          endif
        endif
      enddo
      do i=1,ne
        photar(i)=0.5d0*(prim(i)/ear(i)**2+prim(i-1)/ear(i-1)**2)*
     &   (ear(i)-ear(i-1))*normfac
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Thermal Comptonization; solves Kompaneets eq. with some
c     relativistic corrections. See Lightman & Zdziarski (1987), ApJ
c     The seed spectrum is blackbody.
c  version: January 96
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine thcompton (tempbb, theta, gamma, x, jmax, sptot)
c
      implicit none
      real*8 delta,xmin,xmax,deltal,xnr,xr,taukn,arg,flz,planck,pi
      real*8 sptot(900),x(900),dphesc(900),dphdot(900)
     &  ,rel(900),bet(900),c2(900)
      integer*4 j,jmax,jmaxth,jnr,jrel
      double precision w,w1,z1,z2,z3,z4,z5,z6
      real*8 tautom
c  input parameters:
      real*8 tempbb,theta,gamma
c use internally Thomson optical depth
      tautom=dsqrt(2.25d0+3.d0/(theta*((gamma+.5d0)**2-2.25d0)))-1.5d0
c
      pi=4.d0*datan(1.0d0)
c clear arrays (important for repeated calls)
      do 10 j=1,900
        dphesc(j)=0.d0
        dphdot(j)=0.d0
        rel(j)=0.d0
        bet(j)=0.d0
        c2(j)=0.d0
        sptot(j)=0.d0
   10 continue
c
c JMAX - # OF PHOTON ENERGIES
c
c delta is the 10-log interval of the photon array.
      delta=0.02d0
      deltal=delta*dlog(10.d0)
      xmin=1.d-4*tempbb
      xmax=40.d0*theta
      jmax=min(899,int(dlog10(xmax/xmin)/delta)+1)
c
c X - ARRAY FOR PHOTON ENERGIES
c
      do 20 j=1,jmax+1
        x(j)=xmin*10.d0**(dble(j-1)*delta)
   20 continue
c
c compute c2(x), and rel(x) arrays
c c2(x) is the relativistic correction to Kompaneets equation
c rel(x) is the Klein-Nishina cross section divided by the Thomson crossection
      do 30 j=1,jmax
        w=x(j)
c c2 is the Cooper's coefficient calculated at w1
c w1 is x(j+1/2) (x(i) defined up to jmax+1)
        w1=dsqrt(x(j)*x(j+1))
        c2(j)=(w1**4/(1.d0+4.6d0*w1+1.1d0*w1*w1))
        if (w.le.0.05d0) then
c use asymptotic limit for rel(x) for x less than 0.05
          rel(j)=(1.d0-2.d0*w+26.d0*w*w*0.2d0)
        else
          z1=(1.d0+w)/w**3
          z2=1.d0+2.d0*w
          z3=dlog(z2)
          z4=2.d0*w*(1.d0+w)/z2
          z5=z3/2.d0/w
          z6=(1.d0+3.d0*w)/z2/z2
          rel(j)=(0.75d0*(z1*(z4-z3)+z5-z6))
        endif
   30 continue
c the thermal emission spectrum
      jmaxth=min(900,int(dlog10(50*tempbb/xmin)/delta))
      if (jmaxth.gt.jmax) then
c           print *,'thcomp: ',jmaxth,jmax
        jmaxth=jmax
      endif
      planck=15.d0/(pi*tempbb)**4
      do 40 j=1,jmaxth
        dphdot(j)=planck*x(j)**2/(dexp(x(j)/tempbb)-1)
   40 continue
c
c compute beta array, the probability of escape per Thomson time.
c bet evaluated for spherical geometry and nearly uniform sources.
c Between x=0.1 and 1.0, a function flz modifies beta to allow
c the increasingly large energy change per scattering to gradually
c eliminate spatial diffusion
c
      jnr=int(dlog10(0.1d0/xmin)/delta+1)
      jnr=min(jnr,jmax-1)
      jrel=int(log10(1/xmin)/delta+1)
      jrel=min(jrel,jmax)
      xnr=x(jnr)
      xr=x(jrel)
      do 50 j=1,jnr-1
        taukn=tautom*rel(j)
        bet(j)=1.d0/tautom/(1.d0+taukn/3.d0)
   50 continue
      do 60 j=jnr,jrel
        taukn=tautom*rel(j)
        arg=(x(j)-xnr)/(xr-xnr)
        flz=1-arg
        bet(j)=1.d0/tautom/(1.d0+taukn/3.d0*flz)
   60 continue
      do 70 j=jrel+1,jmax
        bet(j)=1.d0/tautom
   70 continue
c
      call thermlc (tautom, theta, deltal, x, jmax, dphesc, dphdot, bet,
     & c2)
c
c     the spectrum in E F_E
      do 80 j=1,jmax-1
        sptot(j)=dphesc(j)*x(j)**2
c          write(1,*) x(j), sptot(j)
   80 continue
c     the input spectrum
c      do 498 j=1,jmaxth
c         write(2,*) x(j), dphdot(j)*x(j)**2
c 498  continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This program computes the effects of Comptonization by
c nonrelativistic thermal electrons in a sphere including escape, and
c relativistic corrections up to photon energies of 1 MeV.
c the dimensionless photon energy is x=hv/(m*c*c)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine thermlc (tautom, theta, deltal, x, jmax, dphesc,
     &dphdot, bet, c2)
c
c The input parameters and functions are:
c dphdot(x), the photon production rate
c tautom, the Thomson scattering depth
c theta, the temperature in units of m*c*c
c c2(x), and bet(x), the coefficients in the K-equation and the
c   probability of photon escape per Thomson time, respectively,
c   including Klein-Nishina corrections
c The output parameters and functions are:
c dphesc(x), the escaping photon density
      implicit none
      integer*4 jmax
      real*8  tautom,theta,deltal
      real*8 x(900),dphesc(900), dphdot(900),bet(900),c2(900)
      integer*4 j,jj
      real*8  a(900),b(900),c(900),c20,w1,w2,t1,t2,t3,x32,aa
      real*8  d(900),alp(900),g(900),gam(900),u(900)
c u(x) is the dimensionless photon occupation number
      c20=tautom/deltal
c
c determine u
c define coefficients going into equation
c a(j)*u(j+1)+b(j)*u(j)+c(j)*u(j-1)=d(j)
      do 10 j=2,jmax-1
        w1=dsqrt(x(j)*x(j+1))
        w2=dsqrt(x(j-1)*x(j))
c  w1 is x(j+1/2)
c  w2 is x(j-1/2)
        a(j)=-c20*c2(j)*(theta/deltal/w1+0.5d0)
        t1=-c20*c2(j)*(0.5d0-theta/deltal/w1)
        t2=c20*c2(j-1)*(theta/deltal/w2+0.5d0)
        t3=x(j)**3*(tautom*bet(j))
        b(j)=t1+t2+t3
        c(j)=c20*c2(j-1)*(0.5d0-theta/deltal/w2)
        d(j)=x(j)*dphdot(j)
   10 continue
c define constants going into boundary terms
c u(1)=aa*u(2) (zero flux at lowest energy)
c u(jx2) given from region 2 above
      x32=dsqrt(x(1)*x(2))
      aa=(theta/deltal/x32+0.5d0)/(theta/deltal/x32-0.5d0)
c
c zero flux at the highest energy
      u(jmax)=0.d0
c
c invert tridiagonal matrix
      alp(2)=b(2)+c(2)*aa
      gam(2)=a(2)/alp(2)
      do 20 j=3,jmax-1
        alp(j)=b(j)-c(j)*gam(j-1)
        gam(j)=a(j)/alp(j)
   20 continue
      g(2)=d(2)/alp(2)
      do 30 j=3,jmax-2
        g(j)=(d(j)-c(j)*g(j-1))/alp(j)
   30 continue
      g(jmax-1)=(d(jmax-1)-a(jmax-1)*u(jmax)-c(jmax-1)*g(jmax-2))/
     &alp(jmax-1)
      u(jmax-1)=g(jmax-1)
      do 40 j=3,jmax-1
        jj=jmax+1-j
        u(jj)=g(jj)-gam(jj)*u(jj+1)
   40 continue
      u(1)=aa*u(2)
c compute new value of dph(x) and new value of dphesc(x)
      do 50 j=1,jmax
        dphesc(j)=x(j)*x(j)*u(j)*bet(j)*tautom
   50 continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Thermal Comptonization; solves Kompaneets eq. with some
c     relativistic corrections. See Lightman & Zdziarski (1987), ApJ
c     The seed spectrum is DISK blackbody.
c  version: January 96
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine thdscompton (tempbb, theta, gamma, x, jmax, sptot)
      implicit none
      real*8 delta,xmin,xmax,deltal,xnr,xr,taukn,arg,flz,pi
      real*8 sptot(900),x(900),dphesc(900),dphdot(900)
     &  ,rel(900),bet(900),c2(900)
      integer*4 j,jmax,jmaxth,jnr,jrel
      double precision w,w1,z1,z2,z3,z4,z5,z6
      real*8 tautom
c  input parameters:
      real*8 tempbb,theta,gamma
      real*8    ear(0:5000),photar(5000),photer(5000), parth(10)
      integer*4 ifl,ne
c use internally Thomson optical depth
      tautom=dsqrt(2.25d0+3.d0/(theta*((gamma+.5d0)**2-2.25d0)))-1.5d0
c
      pi=4.d0*datan(1.d0)
c clear arrays (important for repeated calls)
      do 10 j=1,900
        dphesc(j)=0.0d0
        dphdot(j)=0.0d0
        rel(j)=0.0d0
        bet(j)=0.0d0
        c2(j)=0.0d0
        sptot(j)=0.0d0
   10 continue
c
c JMAX - # OF PHOTON ENERGIES
c
c delta is the 10-log interval of the photon array.
      delta=0.02d0
      deltal=delta*dlog(10.0d0)
      xmin=1.d-4*tempbb
      xmax=40.0d0*theta
      jmax=min(899,int(dlog10(xmax/xmin)/delta)+1)
c
c X - ARRAY FOR PHOTON ENERGIES
c
      do 20 j=1,jmax+1
        x(j)=xmin*10.d0**((j-1)*delta)
   20 continue
c
c compute c2(x), and rel(x) arrays
c c2(x) is the relativistic correction to Kompaneets equation
c rel(x) is the Klein-Nishina cross section divided by the Thomson crossection
      do 30 j=1,jmax
        w=x(j)
c c2 is the Cooper's coefficient calculated at w1
c w1 is x(j+1/2) (x(i) defined up to jmax+1)
        w1=dsqrt(x(j)*x(j+1))
        c2(j)=(w1**4/(1.0d0+4.6d0*w1+1.1d0*w1*w1))
        if (w.le.0.05d0) then
c use asymptotic limit for rel(x) for x less than 0.05
          rel(j)=(1.d0-2.d0*w+26.d0*w*w*0.2d0)
        else
          z1=(1+w)/w**3
          z2=1.d0+2.d0*w
          z3=dlog(z2)
          z4=2.d0*w*(1.d0+w)/z2
          z5=z3/2.d0/w
          z6=(1.d0+3.d0*w)/z2/z2
          rel(j)=(0.75d0*(z1*(z4-z3)+z5-z6))
        endif
   30 continue
c the thermal emission spectrum
      jmaxth=min(900,int(dlog10(50*tempbb/xmin)/delta))
      if (jmaxth.gt.jmax) then
        print *,'thcomp: ',jmaxth,jmax
        jmaxth=jmax
      endif
c        planck=15/(pi*tempbb)**4
c        do 5 j=1,jmaxth
c          dphdot(j)=planck*x(j)**2/(exp(x(j)/tempbb)-1)
c  5     continue
      do j=1,jmaxth-1
        ear(j-1)=511.d0*dsqrt(x(j)*x(j+1))
      enddo
      parth(1)=tempbb*511
      ne=jmaxth-2
      call xsdskb (ear, ne, parth, ifl, photar, photer)
      do j=1,ne
        dphdot(j+1)=511.d0*photar(j)/(ear(j)-ear(j-1))
      enddo
      jmaxth=ne+1
      dphdot(1)=dphdot(2)
c
c compute beta array, the probability of escape per Thomson time.
c bet evaluated for spherical geometry and nearly uniform sources.
c Between x=0.1 and 1.0, a function flz modifies beta to allow
c the increasingly large energy change per scattering to gradually
c eliminate spatial diffusion
      jnr=int(dlog10(0.1/xmin)/delta+1)
      jnr=min(jnr,jmax-1)
      jrel=int(dlog10(1/xmin)/delta+1)
      jrel=min(jrel,jmax)
      xnr=x(jnr)
      xr=x(jrel)
      do 40 j=1,jnr-1
        taukn=tautom*rel(j)
        bet(j)=1.d0/tautom/(1.d0+taukn/3.d0)
   40 continue
      do 50 j=jnr,jrel
        taukn=tautom*rel(j)
        arg=(x(j)-xnr)/(xr-xnr)
        flz=1.d0-arg
        bet(j)=1.d0/tautom/(1.d0+taukn/3.d0*flz)
   50 continue
      do 60 j=jrel+1,jmax
        bet(j)=1/tautom
   60 continue
c
      call thermlc (tautom, theta, deltal, x, jmax, dphesc, dphdot, bet,
     & c2)
c
c     the spectrum in E F_E
      do 70 j=1,jmax-1
        sptot(j)=dphesc(j)*x(j)**2
c          write(1,*) x(j), sptot(j)
   70 continue
c      print *,'jmax: ',jmax,jmaxth
c      open(33,file='spec.dat')
cc     the input spectrum
c      do 498 j=1,min(jmaxth,jmax-1)
c         write(33,*) 511*x(j), dphdot(j)*x(j), dphesc(j)*x(j)
c 498  continue
c      close(33)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Auxillary function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function spp(y,xnonth,nnonth,spnth)
      implicit none
      integer*4 nnonth
      real*8 y , xnonth(nnonth) , spnth(nnonth)
      integer*4 il , ih
      real*8 xx
      save ih
      data ih/2/
      xx=1.d0/y
      if (xx.lt.xnonth(ih)) ih=2
      do while (ih.lt.nnonth.and.xx.gt.xnonth(ih))
        ih=ih+1
      enddo
      il=ih-1
      spp=spnth(il)+(spnth(ih)-spnth(il))*(xx-xnonth(il))/(xnonth(ih)-
     &xnonth(il))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Multicolour disk blackbody model used in ISAS, Japan.
c     See Mitsuda et al. PASJ, 36, 741 (1984)
c     & Makishima et al. ApJ 308, 635 1986)
c     Ken Ebisawa 1992/12/22
c     Modified to use real*8 and to make numerical
c     integration faster.
c     Ken Ebisawa 1993/07/29
c     Modified the algorithm. by Kazuhisa MITSUDA May 31, 1994
c     A numerical calculation was first done within an accruacy of 0.01 %.
c     Then an interpolation formula was found.
c     The interpolation is precise within 1e-5 level.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xsdskb (ear, ne, param, idt, photar, photer)
      integer*4 ne, idt
      real*8   ear(0:ne), param(1), photar(ne), photer(ne)
      integer*4 i, j
      real*8 xn, xh
      real*8 tin, e, photon
c     These coefficients are taken from the spectral fitting program SPFD
c     in ISAS.
c     These are used for Gaussian Integral in the given energy band
      real*8  gauss(5,2)
      save gauss
      data  gauss /
     &     0.236926885d0,  0.478628670d0, 0.568888888d0, 0.478628670d0,
     &     0.236926885d0, -0.906179846d0, -0.538469310d0, 0.0d0,
     &     0.538469310d0,  0.906179846d0 /
c suppress a warning message from the compiler
      i=idt
c this model has no errors
      do i=1,ne
        photer(i)=0.0d0
      enddo
      tin=dble(param(1))
      do 20 i=1,ne
        xn=(ear(i)-ear(i-1))/2.0d0
        photar(i)=0.0d0
        xh=xn+ear(i-1)
        do 10 j=1,5
          e=dble(xn)*gauss(j,2)+dble(xh)
          call mcdspc (e, tin, 1.0d0, photon)
          photar(i)=photar(i)+dble(gauss(j,1)*photon)
   10   continue
        photar(i)=photar(i)*xn
   20 continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Multi-Color Disk SPECTRUM
c     SEE MITSUDA ET AL. 1984 PASJ 36, 741.
c     & MAKISHIMA ET AL. APJ 308, 635 1986)
c     NORMALIZATION={RIN(KM)/(D/10KPC)}^2*COS(THETA)
c     TIN=     INNER TEMPERATURE OF THE DISK
c     E  =     ENERGY
c     Rin2 = Normalization factor in the unit of above normalization.
c     PHOTON = PHOTONS/S/KEV/CM2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mcdspc (e, tin, rin2, flux)
      real*8 e, tin, rin2, flux
c  E = X-ray energy (keV)
c  Tin = inner edge color-temperature (keV)
c  Rin2 = inner edge radius(km) ^2 * cos(inclination)/ [D (10kpc)]^2
c  Flux = photon flux, photons/sec/cm^2/keV
      real*8 normfact
      parameter (normfact=361.0d0)
      real*8 et
      real*8 value
      if (tin.eq.0.0d0) then
        flux=0.0d0
        return
      endif
      et=e/tin
      call mcdint (et, value)
      flux=value*tin*tin*rin2/normfact
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mcdint (et, value)
      implicit none
      integer*4 nres
      parameter(nres=98)
      real*8 value0, et0, step, a, b, beki
      parameter(value0=0.19321556d+03, et0=0.10000000d-02)
      parameter(step=0.60000000d-01, a=0.52876731d+00)
      parameter(b=0.16637530d+01, beki=-2.0d0/3.0d0)
      real*8 et, value
      real*8 gc(3),gw(3),gn(3),res(nres)
      real*8 loget, z, pos
      real*8 resfact, gaufact
      integer*4 j
      save gc, gw, gn, res
      data gc/  0.78196667d-01, -0.10662020d+01,  0.11924180d+01/
      data gw/  0.52078740d+00,  0.51345700d+00,  0.40779830d+00/
      data gn/  0.37286910d+00,  0.39775528d-01,  0.37766505d-01/
      data (res(j),j=1,40)/
     &   0.96198382d-03, 0.10901181d-02, 0.12310012d-02, 0.13841352d-02,
     &   0.15481583d-02, 0.17210036d-02, 0.18988943d-02, 0.20769390d-02,
     &   0.22484281d-02, 0.24049483d-02, 0.25366202d-02, 0.26316255d-02,
     &   0.26774985d-02, 0.26613059d-02, 0.25708784d-02, 0.23962965d-02,
     &   0.21306550d-02, 0.17725174d-02, 0.13268656d-02, 0.80657672d-03,
     &   0.23337584d-03,-0.36291778d-03,-0.94443569d-03,-0.14678875d-02,
     &  -0.18873741d-02,-0.21588493d-02,-0.22448371d-02,-0.21198179d-02,
     &  -0.17754602d-02,-0.12246034d-02,-0.50414167d-03, 0.32507078d-03,
     &   0.11811065d-02, 0.19673402d-02, 0.25827094d-02, 0.29342526d-02,
     &   0.29517083d-02, 0.26012166d-02, 0.18959062d-02, 0.90128649d-03/
      data (res(j),j=41,80)/
     &  -0.26757144d-03,-0.14567885d-02,-0.24928550d-02,-0.32079776d-02,
     &  -0.34678637d-02,-0.31988217d-02,-0.24080969d-02,-0.11936240d-02,
     &   0.26134145d-03, 0.17117758d-02, 0.28906898d-02, 0.35614435d-02,
     &   0.35711778d-02, 0.28921374d-02, 0.16385898d-02, 0.49857464d-04,
     &  -0.15572671d-02,-0.28578151d-02,-0.35924212d-02,-0.36253044d-02,
     &  -0.29750860d-02,-0.18044436d-02,-0.37796664d-03, 0.10076215d-02,
     &   0.20937327d-02, 0.27090854d-02, 0.28031667d-02, 0.24276576d-02,
     &   0.17175597d-02, 0.81030795d-03,-0.12592304d-03,-0.94888491d-03,
     &  -0.15544816d-02,-0.18831972d-02,-0.19203142d-02,-0.16905849d-02,
     &  -0.12487737d-02,-0.66789911d-03,-0.27079461d-04, 0.59931935d-03/
      data (res(j),j=81,98)/
     &   0.11499748d-02, 0.15816521d-02, 0.18709224d-02, 0.20129966d-02,
     &   0.20184702d-02, 0.19089181d-02, 0.17122289d-02, 0.14583770d-02,
     &   0.11760717d-02, 0.89046768d-03, 0.62190822d-03, 0.38553762d-03,
     &   0.19155022d-03, 0.45837109d-04,-0.49177834d-04,-0.93670762d-04,
     &  -0.89622968d-04,-0.401538532d-04/
      loget=dlog10(et)
      pos=(loget-dlog10(et0))/step+1
      j=int(pos)
      if (j.lt.1) then
        resfact=res(1)
      elseif (j.ge.nres) then
        resfact=res(nres)
      else
        pos=pos-j
        resfact=res(j)*(1.0d0-pos)+res(j+1)*pos
      endif
      gaufact=1.0d0
      do j=1,3
        z=(loget-gc(j))/gw(j)
        gaufact=gaufact+gn(j)*dexp(-z*z/2.0d0)
      enddo
      value=value0*(et/et0)**beki*(1.0d0+a*et**b)*dexp(-et)*gaufact*
     &(1.0d0+resfact)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c relativistic compton scattering of wien-law spectrum
c CONVERTED FROM TITARCHUK'S STAND-ALONE CODE. HANDLES OPTICALLY
c THIN AND THICK CASES. ADAPTED FOR XSPEC January 1995, T. Yaqoob.
c **       IMPORTANT **
c      You MUST consult Titarchuk 1994 (ApJ, 434, 313)
c       and HUA & Titarchuk 1994 (ApJ xxx,xxx) to fully understand
c      and appreciate the physical assumptions and approximations.
c       This routine will not warn you when the physical assumptions
c      break down.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xstitg (ear, npts, param, photar, photer)
      implicit none
      integer*4 npts
      real*8 ear(0:npts) , param(5) , photar(npts), photer(npts)
c---
      real*8 aa , z , pi , gam0 , yyit2 , argg , rdel , taue ,
     &                 xr , xt
      real*8 x , temp , t5 , beta , x0 , ro
      real*8 fx , alfax0 , algx0 , al3x0 , bol3
      real*8 taup , tb1 , tb2 , t23 , bol2i , bol7i , bol6
      real*8 compd0 , tal , tal0 , bola , betaint
      real*8 gammi , gamln , bol3g , arg1 , dumt
      real*8 xrfac , tfx , f0theta
      integer*4 i , ii
      real*8 zfac , t0 , apprx
      real*8 ens, oens, phot, ophot
c=====>
c 1. Z    ** REDSHIFT
c 2. T0   ** WIEN TEMPERATURE (keV)
c 3. TEMP ** PLASMA TEMPERATURE (KEV)
c 4. TAUP ** PLASMA OPTICAL DEPTH
c 5. APPRX * ABS(APPRX) <= 1.0 DISK
c           ABS(APPRX) >  1.0 SPHERE
c           APPRX < 0 : Get beta as a function of optical depth by
c                   by interpolation of accurately calculated points
c                   in Sunyaev & Titarchuk 1985 (A&A, 143, 374).
c           APPRX >= 0: Get beta as a function of optical depth from
c                   analytic approximation (e.g. Titarchuk 1994).
c suppress a warning message from the compiler
      pi=4.d0*datan(1.d0)
c this model has no errors
      do i=1,npts
        photer(i)=0.0d0
      enddo
      ro=1.d0
      bol2i=0.0d0
      temp=0.0d0
      taup=0.0d0
      zfac=1.d0+param(1)
      t0=param(2)
      temp=param(3)
      dumt=param(4)
      taup=dumt
      apprx=param(5)
      t5=temp/511.d0
c dimensionless soft photon energy
      x0=3.d0*t0/temp
      t23=taup+(2.d0/3.d0)
      taue=taup/(1.0d0+(temp/39.2d0)**0.86d0)
      rdel=0.0d0
      beta=0.0d0
      if (dabs(apprx).gt.1.0d0) then
c BETA FOR SPHERE
        if (taup.le.0.1d0) then
          beta=dlog(1.d0/0.75d0/taup)
        elseif (taup.ge.10.d0.or.apprx.ge.0.0d0) then
          tb1=pi*pi*(1.d0-dexp(-0.7d0*taup))/3.d0/t23/t23
          tb2=dexp(-1.4d0*taup)*dlog(4.d0/(3.d0*taup))
          beta=tb1+tb2
        elseif (taup.gt.0.1d0.and.taup.lt.10.d0.and.apprx.lt.0.0d0)
     &   then
          beta=betaint(taup,apprx)
        endif
c                 write(*,*) 'tau & beta = ',taup,beta
        if (taue.lt.0.01d0) then
          rdel=taue/2.d0
        else
          rdel=1.0-3.d0/taue*(1.d0-2.d0/taue+2.d0/taue/taue*(1d0-dexp(-
     &     taue)))
        endif
      elseif (dabs(apprx).le.1.0d0) then
c BETA FOR DISK
        if (taup.le.0.1d0) then
          beta=dlog(1.d0/taup/dlog(1.53d0/taup))
        elseif (taup.ge.10.d0.or.apprx.ge.0.0d0) then
          tb1=pi*pi*(1.d0-exp(-1.35d0*taup))/12.d0/t23/t23
          tb2=0.45d0*dexp(-3.7d0*taup)*dlog(10.d0/(3.d0*taup))
          beta=tb1+tb2
        elseif (taup.gt.0.1d0.and.taup.lt.10.d0.and.apprx.lt.0.0d0)
     &   then
          beta=betaint(taup,apprx)
        endif
c       write(*,*) 'tau beta = ',taup,beta
        if (taue.lt.0.01d0) then
          rdel=taue/2.d0
        else
          rdel=1.0d0-(1.0d0-dexp(-taue))/taue
        endif
      endif
      f0theta=2.5d0*t5+1.875d0*t5*t5*(1.0d0-t5)
      gam0=beta/t5
c        write(*,*) 'gam0 = ',gam0
c compute the power-law index
      tal0=dsqrt(2.25d0+gam0)-1.5d0
c     These changes by CM 22 Mar 2001 to improve the convergence of the
c     solution.  Now the loop goes for a maximum of 50 iterations, and
c     the convergence criterium is looser.  Finally, the corrected
c     expression using COMPD0 has been replaced.
      if (tal0.gt.10.d0) tal0=10.d0
      do ii=1,50
        bola=1.d0+(tal0+3.d0)*t5/(1.d0+t5)+4.d0*compd0(tal0)**(1.d0/
     &   tal0)*t5*t5
        tal=beta/dlog(bola)
        if (dabs(tal-tal0).le.1.0d-4) goto 10
        tal0=tal
      enddo
   10 continue
      alfax0=tal
      if (alfax0.gt.10.d0) alfax0=10.d0
c follwoing is outdated
c      IF (ithick.eq.1) then
c        gamx0=gam0/(1.+f0theta)
c      elseif (ithick.eq.0) then
c        gamx0=gam0*(8.+15.*t5)/(8.+19.*t5)
c      endif
      aa=3.d0/x0
c      alfax0 = SQRT(2.25+Gamx0) - 1.5
      argg=3.d0+alfax0
      arg1=dble(argg)
      xrfac=rdel**(1.d0/argg)
      algx0=2.0d0*alfax0+3.0d0
      al3x0=alfax0*(alfax0+3.d0)/2.d0/algx0
      ens=0.d0
      oens=0.d0
      phot=0.d0
      ophot=0.d0
      do 20 i=0,npts
        ens=ear(i)*zfac
        if (ens.gt.0.0d0) then
          xt=ens/temp
          x=xt
          z=xt*t5
c            write(*,*) 'z x =',z,x
          xr=x*xrfac
c            write(*,*) 'xr = ',xr
          bol3=x*aa
c            write(*,*) 'bol3 = ',bol3
          bol3g=dble(bol3)
c             write(*,*) 'bol3g = ',bol3g
c see if fx is going to overflow
          tfx=(-1.d0-alfax0)*dlog10(bol3)-x*0.4343d0+dlog10(algx0)
c      write(*,*) 'tfx = ',tfx
          if (dabs(tfx).lt.27.d0) then
            fx=bol3**(-alfax0-1.d0)*dexp(-x)*algx0
          else
            fx=0.0d0
          endif
c           write(*,*) 'arg1 bol3g gammi = ',arg1,bol3g,gammi(arg1,bol3g)
c         write(*,*) ' gamln exp(gamln)= ',gamln(arg1),exp(gamln(arg1))
          bol2i=gammi(arg1,bol3g)*dexp(gamln(arg1))
c         write(*,*) 'bol2i =',bol2i
c         write(*,*) 'fx yyit2',fx,yyit2(xr,alfax0,ro)
          bol7i=yyit2(xr,alfax0,ro)*fx*bol2i
c         write(*,*) 'bol7i = ',bol7i
          bol6=bol3**2.d0*dexp(-bol3)/(alfax0+3.d0)
          phot=dble(al3x0*(bol7i+bol6)*aa)
c      if (i.eq.0) then
c        write(*,*) 'beta = ',beta
c        write(*,*) 'f0theta = ',f0theta
c        write(*,*) 'gamx0 = ',gamx0
c        write(*,*) 'alfax0 = ',alfax0
c        write(*,*) 'argg = ',argg
c        write(*,*) 'bol2i = ',bol2i
c        write(*,*) 'x = ',x
c        write(*,*) 'bol3 = ',bol3
c        write(*,*) 'fx = ',fx
c        write(*,*) 'bol7i bol6 ',bol7i,bol6
c        write(*,*) 'gamma(arg) ',exp(gamln(arg1))
c        write(*,*) 'igamma   ',gammi(arg1,bol3g)
c      endif
        else
          phot=0.0d0
        endif
c            write(2,*) ens(i),bol7i,bol6,bol6+bol7i
        if (i.gt.0) then
          photar(i)=0.5d0*(phot+ophot)*(ens-oens)/zfac
        endif
        oens=ens
        ophot=phot
   20 continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function yyit2(x,alfa,ro)
      implicit none
c*** Start of declarations inserted by SPAG
      real*8 a2 , a3 , v , w , x , z , ro
      integer*4 i
c*** End of declarations inserted by SPAG
      real*8 alfa , db , a2n , v1
      real*8 gamln , argg
      dimension w(10) , z(10)
      data z/.1377934705d0, .7294545495d0, 1.8083429017d0,
     &     3.4014336978d0,5.55249614d0, 8.3301527467d0,
     &     11.8437858379d0, 16.2792578313d0,
     &     21.9965858119d0, 29.9206970122d0/
      data w/3.0844111576d-01 , 4.0111992915d-01 , 2.180682876d-01 ,
     &     6.208745609d-02 , 9.501516975d-03 , 7.530083885d-04 ,
     &     2.825923349d-05 , 4.249313984d-07 , 1.839564823d-09 ,
     &     9.911827219d-13/
      a2=alfa+3.d0
      a2n=alfa+2.d0
      a3=alfa-1.d0
      argg=dble(a2+a3+2.d0)
      db=gamln(argg)
      yyit2=0.d0
      do 10 i=1,10
        v=w(i)*dexp(a2n*dlog(ro*x+z(i))+alfa*dlog(z(i))-db)
c          write(*,*) 'alfa x z a2 ',alfa,x,z(i),a2
c          write(*,*) alfa*((x+z(i)) - a2)
        v1=v/alfa*((x+z(i))-a2)
c            v = w(i)*DEXP(a2*DLOG(ro*x+z(i))+a3*Dlog(z(i))-db)
        yyit2=yyit2+v1
   10 continue
      yyit2=yyit2
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SERIES REPRESENTATION OF INCOMPLETE GAMMA FUNCTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function gamln(az)
      implicit none
      real*8 z , a(7) , s , pi
      real*8 az
      integer*4 i
      pi=4.d0*datan(1.d0)
      a(1)=1.d0/12.d0
      a(2)=1.d0/30.d0
      a(3)=53.d0/210.d0
      a(4)=195.d0/371.d0
      a(5)=22999.d0/22737.d0
      a(6)=29944523.d0/19733142.d0
      a(7)=109535241009.d0/48264275462.d0
      z=0.d0
      z=az
c      write(*,*) 'gamln : az  z = ',az,z
      s=z
      do 10 i=1,6
        s=z+a(8-i)/s
   10 continue
      s=a(1)/s
      s=s-z+(z-0.5d0)*dlog(z)+0.5d0*dlog(2.d0*pi)
      gamln=s
      return
      end
c***********************************************************
      real*8 function gammi(a,x)
      implicit none
c CODED FROM NUMERICAL RECIPIES
c INCOMPLETE GAMMA FUNCTION
      real*8 a , x , gln , gammcf , gamser
      if (x.lt.0.d0.or.a.le.0.d0) then
      write(*,*) 'Inc. Gamma Fn. called with x < 0 or A <= 0'
      endif
c USE THE SERIES REPRESENTATION
      if (x.lt.(a+1.d0)) then
        call gserr (gamser, a, x, gln)
        gammi=gamser
      else
        call gcff (gammcf, a, x, gln)
        gammi=1.d0-gammcf
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gserr (gamser, a, x, gln)
      implicit none
      real*8 a , x , gln , gamser , ap , sum , del , eps
      real*8 gamln
      integer*4 itmax , n
      parameter (itmax=10000,eps=3.d-12)
c      write(*,*) 'gserr: a  x  logx',a,x,log(x)
      gln=gamln(a)
      if (x.le.0.0d0) then
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 10 n=1,itmax
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if (abs(del).lt.abs(sum)*eps) goto 20
   10 continue
      write (*,*)
     &'Warning: Inc. Gamma Fn. GSERR did not converge for A = ',a
   20 gamser=sum*dexp(-x+a*log(x)-gln)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c CONTINUED FRACTION METHOD FOR INCOMPLETE GAMMA FUNCTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gcff (gammcf, a, x, gln)
      implicit none
      real*8 a0 , b0 , a1 , b1 , gln , eps , gammcf , a , x , anf , g ,
     &     gold , fac , ana , an
      integer*4 itmax , n
      real*8 gamln
      parameter (itmax=10000,eps=3.d-12)
c      write(*,*) 'gcff: a  x  logx',a,x,log(x)
      gln=gamln(a)
      gold=0.0d0
      a0=1.d0
      a1=x
      b0=0.d0
      b1=1.d0
      fac=1.d0
      g=0.0d0
      do 10 n=1,itmax
        an=dble(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if (a1.ne.0.0d0) then
c            write(*,*) 'gcff: a1 = ',a1
          fac=1.d0/a1
          g=b1*fac
          if (dabs((g-gold)/g).lt.eps) goto 20
          gold=g
        endif
   10 continue
      write (*,*)
     &'Warning: Inc. Gamma Fn. GCFF did not converge for A = ',a
   20 gammcf=dexp(-x+a*dlog(x)-gln)*g
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function compd0(x)
      implicit none
      real*8 x , x1 , x2 , x3 , bol , db
      real*8 gamln , argg
      x1=x+1.d0
      x2=x+2.d0
      x3=x+3.d0
      bol=x3*x+4.d0
      argg=(2.d0*x1)
      db=gamln(argg)
      if (db.lt.-70.d0) then
        compd0=0.0d0
      elseif (db.gt.70.d0) then
        compd0=1.d34
      else
        compd0=3.d0*bol*dexp(db)/x3/x2/x2
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Get beta parameter as a function of optical depth for sphere or
c disk geometry by interpolating accurately calculated points in
c Sunyaev & Titrachuk (1985), A&A 143, 374. Optical depth must be
c inside the range 0.1-10.0.
c DISK: abs(apprx) le 1.0
c SPHERE: abs(apprx) gt 1.0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function betaint(tau,apprx)
      implicit none
      integer*4 npts , ix
      parameter (npts=99)
      real*8 apprx
      real*8 sbetal(npts) , sbetah(npts) , taul(npts) ,
     &                 tauh(npts)
      real*8 dbetal(npts) , dbetah(npts) , tau
      real*8 tau1 , tau2 , b1 , b2 , b3
c
      data taul/0.1d0 , 0.2d0 , 2*0.3d0 , 2*0.5d0 , 3*0.7d0 , 5*1.0d0 ,
     &     5*1.5d0 , 5*2.0d0 , 5*2.5d0 , 10*3.0d0 , 10*4.d0 , 20*5.0d0 ,
     &     30*7.d0/
      data tauh/0.2d0 , 0.3d0 , 2*0.5d0 , 2*0.7d0 , 3*1.0d0 , 5*1.5d0 ,
     &     5*2.0d0 , 5*2.5d0 , 5*3.0d0 , 10*4.0d0 , 10*5.d0 , 20*7.0d0 ,
     &     30*10.d0/
      data sbetal/2.97d0 , 1.97d0 , 2*1.60d0 , 2*1.18d0 , 3*0.925d0 ,
     &     5*0.688d0 , 5*0.462d0 , 5*0.334d0 , 5*0.252d0 , 10*0.197d0 ,
     &     10*0.130d0 , 20*0.092d0 , 30*0.0522d0/
      data sbetah/1.97d0 , 1.60d0 , 2*1.18d0 , 2*0.925d0 , 3*0.688d0 ,
     &     5*0.462d0 , 5*0.334d0 , 5*0.252d0 , 5*0.197d0 , 10*0.130d0 ,
     &     10*0.092d0 , 20*0.0522d0 , 30*0.0278d0/
      data dbetal/1.23d0 , 0.84d0 , 2*0.655d0 , 2*0.45d0 , 3*0.335d0 ,
     &     5*0.234d0 , 5*0.147d0 , 5*0.101d0 , 5*0.0735d0 ,
     &     10*0.0559d0 , 10*0.0355d0 , 20*0.0245d0 , 30*0.013393d0/
      data dbetah/0.84d0 , 0.655d0 , 2*0.45d0 , 2*0.335d0 , 3*0.234d0 ,
     &     5*0.147d0 , 5*0.101d0 , 5*0.0735d0 , 5*0.0559d0 ,
     &     10*0.0355d0 , 10*0.0245d0 , 20*0.013393d0 , 30*0.00706d0/
c get look-up index
      ix=int(tau*10.d0)
      tau1=taul(ix)
      tau2=tauh(ix)
      if (abs(apprx).gt.1.0d0) then
        b1=sbetal(ix)
        b2=sbetah(ix)
      elseif (abs(apprx).le.1.0d0) then
        b1=dbetal(ix)
        b2=dbetah(ix)
      endif
c      write(*,*) 'tau = ',tau
      call dinter (tau, tau1, tau2, b3, b1, b2)
      betaint=b3
c      write(*,*) tau1,tau2,b1,b2,b3
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Interpolate in log- space
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dinter (x0, x1, x2, y0, y1, y2)
      implicit none
      real*8 x0 , x1 , x2 , y0 , y1 , y2 , tn , tn1 , t0
      if (x1.ne.x2) then
        if (y1.gt.0.0d0) then
          tn=-dlog10(y1)
        else
          tn=36.0d0
        endif
        if (y2.gt.0.0d0) then
          tn1=-dlog10(y2)
        else
          tn1=36.0d0
        endif
        t0=(tn*(x2-x0)+tn1*(x0-x1))/(x2-x1)
        if (tn.ge.36.d0.and.tn1.ge.36.d0) then
          y0=0.0d0
        else
          y0=10.0d0**(-t0)
        endif
      else
        y0=y1
      endif
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
