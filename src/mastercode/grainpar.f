cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine solvegraink ()
c
      include 'cblocks.inc'
c
      real*8 volume,pow,pow1,const,fw
      real*8 dgrad,dustmax
      integer*4 i,j,atom
c
c     Solve for factor k in distribution (num dens wrt H)
c     and determine distribution
c
      if ((gdist.eq.'M').or.(gdist.eq.'P')) then
c
c      Power law distribution (incl MRN)
c
c
c     Graphite grains
c
c     solve for weight of dust
c
        atom=zmap(6)
        fw=atwei(atom)*zion0(atom)*deltazion(atom)*(1.d0-dion(atom))*
     &   (1.d0-pahcfrac)
        fw=fw*amu
        dustmass(1)=fw
c
c     Then determine 'volume and solve for k
c
        if (mindust(1).eq.maxdust(1)) then
          volume=ftpi*amax(1)**3.d0
        else
          if (galpha.ne.3.d0) then
            pow=galpha+3.d0
            pow1=pow+1.d0
            volume=ftpi*(((amax(1)**pow1)/pow1)-((amin(1)**pow1)/pow1))
          else
            volume=ftpi*(log(amax(1))-log(amin(1)))
          endif
        endif
        const=fw/(volume*graindens(1))
c
c     Evaluate graphite grain distribution per H atom
c     (1e4 is diff factor between cubic and square microns in cm)
c
        if (mindust(1).eq.maxdust(1)) then
          dgrad=gradedge(mindust(1)+1)-gradedge(mindust(1))
          dustnum(maxdust(1),1)=const/dgrad
          dustsig(maxdust(1),1)=pi*const*grainrad(maxdust(1))*
     &     grainrad(maxdust(1))/dgrad
        else
          do i=mindust(1),maxdust(1)
            dustnum(i,1)=const*grainrad(i)**galpha
            dustsig(i,1)=pi*const*grainrad(i)**(galpha+2.d0)
          enddo
        endif
c
c     Silicate grains - O, Mg, Al, Si, Fe
c
c     solve for weight of dust
c
        fw=0.d0
        do i=3,atypes
        if ((mapz(i).ne.6).and.(mapz(i).ne.10).and.(mapz(i).ne.18)) then
          fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
        endif
C         if (mapz(i).eq.8) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.12) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.13) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.14) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.26) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
        enddo
        fw=fw*amu
        dustmass(2)=fw
c
c     Then determine volume and solve for k
c
        if (mindust(2).eq.maxdust(2)) then
          volume=ftpi*amax(2)**3.d0
        else
          if (galpha.ne.3.d0) then
            pow=galpha+3.d0
            pow1=pow+1.d0
            volume=ftpi*(((amax(2)**pow1)/pow1)-((amin(2)**pow1)/pow1))
          else
            volume=ftpi*(log(amax(2))-log(amin(2)))
          endif
        endif
        const=fw/(volume*graindens(2))
c
c     Evaluate silicate grain distribution per H atom
c
        if (mindust(2).eq.maxdust(2)) then
          dgrad=gradedge(mindust(2)+1)-gradedge(mindust(2))
          dustnum(maxdust(2),2)=const/dgrad
          dustsig(maxdust(2),2)=pi*const*grainrad(maxdust(2))*
     &     grainrad(maxdust(2))/dgrad
        else
          do i=mindust(2),maxdust(2)
            dustnum(i,2)=const*grainrad(i)**galpha
            dustsig(i,2)=pi*const*grainrad(i)**(galpha+2.d0)
          enddo
        endif
      else if (gdist.eq.'S') then
c
c    Grain shattering distribution
c
c     Graphite grains
c
c     solve for weight of dust
c
        atom=zmap(6)
        fw=atwei(atom)*zion0(atom)*deltazion(atom)*(1.d0-dion(atom))*
     &   (1.d0-pahcfrac)
c        fw = atwei(zmap(6))*zion(zmap(6))*(1.d0-dion(zmap(6)))
        fw=fw*amu
c
c     Then determine 'volume and solve for k
c
        volume=0.d0
        do i=1,dustbinmax
          dgrad=gradedge(i+1)-gradedge(i)
          dustnum(i,1)=grainrad(i)**(galpha)*exp(-(grainrad(i)/amin(1))*
     &     *(-3)-(grainrad(i)/amax(1)))
          volume=volume+dustnum(i,1)*grainrad(i)**3*dgrad
        enddo
        const=fw/(ftpi*volume*graindens(1))
c
c     Evaluate graphite grain distribution per H atom
c     (1e4 is diff factor between cubic and square microns in cm)
c
        dustmax=0.d0
        do i=1,dustbinmax
          dustnum(i,1)=const*dustnum(i,1)
          dustsig(i,1)=pi*grainrad(i)**2*dustnum(i,1)
          if (dustnum(i,1).gt.dustmax) dustmax=dustnum(i,1)
        enddo
        i=1
        do while (dustnum(i,1).lt.1.0d-10*dustmax)
          i=i+1
        enddo
        mindust(1)=i
        i=i+1
        do while ((dustnum(i,1).gt.1.0d-10*dustmax)
     &   .and.(i.lt.dustbinmax))
          i=i+1
        enddo
        maxdust(1)=i
c
c     Silicate grains - O, Mg, Al, Si, Fe
c
c     solve for weight of dust
c
        fw=0.d0
        do i=3,atypes
        if ((mapz(i).ne.6).and.(mapz(i).ne.10).and.(mapz(i).ne.18)) then
            fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
        endif
C         if (mapz(i).eq.8) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.12) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.13) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.14) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
C         if (mapz(i).eq.26) then
C           fw=fw+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
C         endif
        enddo
        fw=fw*amu
c
c     Then determine volume and solve for k
c
        volume=0.d0
        do i=1,dustbinmax
          dgrad=gradedge(i+1)-gradedge(i)
          dustnum(i,2)=grainrad(i)**(galpha)*exp(-(grainrad(i)/amin(2))*
     &     *(-3)-(grainrad(i)/amax(2))**(3))
          volume=volume+dustnum(i,2)*grainrad(i)**3*dgrad
        enddo
        const=fw/(ftpi*volume*graindens(2))
c
c     Evaluate silicate grain distribution per H atom
c
        dustmax=0.d0
        do i=1,dustbinmax
          dustnum(i,2)=const*dustnum(i,2)
          dustsig(i,2)=pi*grainrad(i)**2*dustnum(i,2)
          if (dustnum(i,2).gt.dustmax) dustmax=dustnum(i,2)
        enddo
        i=1
        do while (dustnum(i,2).lt.1.0d-10*dustmax)
          i=i+1
        enddo
        mindust(2)=i
        i=i+1
        do while ((dustnum(i,2).gt.1.0d-10*dustmax)
     &   .and.(i.lt.dustbinmax))
          i=i+1
        enddo
        maxdust(2)=i
      endif
c
c     calc ratios
c
      graindgr=0.d0
c
      do i=3,atypes
        if ((mapz(i).ne.6).and.(mapz(i).ne.10).and.(mapz(i).ne.18)) then
C H, He, Ne and Ar not in dust
C Carbon treated differently
        graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
        endif
C       if (mapz(i).eq.8) then
C         graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i)
C    &     )
C       endif
C       if (mapz(i).eq.12) then
C         graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i)
C    &     )
C       endif
C       if (mapz(i).eq.14) then
C         graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i)
C    &     )
C       endif
C       if (mapz(i).eq.13) then
C         graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i)
C    &     )
C       endif
C       if (mapz(i).eq.26) then
C         graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i)
C    &     )
C       endif
        if (mapz(i).eq.6) then
        graindgr=graindgr+atwei(i)*zion0(i)*deltazion(i)*(1.d0-dion(i))
     &     *(1.d0-pahcfrac)
        endif
      enddo
c
      grainghr=0.d0
c
      do i=1,atypes
        grainghr=grainghr+atwei(i)*zion0(i)*deltazion(i)*dion(i)
      enddo
c
      graindgr=graindgr/grainghr
c
      atom=zmap(6)
      if ((pahmode.eq.1).and.(grainghr.gt.0.d0)) then
        grainpgr=(atwei(atom)*pahfrac*pahnc)/grainghr
      endif
c
      do i=1,numtypes
        siggrain(i)=0.d0
        do j=mindust(i),maxdust(i)
          dgrad=gradedge(j+1)-gradedge(j)
          siggrain(i)=siggrain(i)+dustsig(j,i)*dgrad
        enddo
      enddo
      return
c
      end
      subroutine graindepletegas (taper)
c
      include 'cblocks.inc'
c
      real*8 taper,delta
      integer*4 i
      do i=3,atypes
        delta=1.d0-dion0(i)
        dion(i)=dion0(i)+(delta*(1.d0-taper))
        invdion(i)=1.d0/dion(i)
        zion(i)=zion0(i)*deltazion(i)*dion(i)
        if ((clinpah.eq.0).and.(mapz(i).eq.6)) then
          zion(i)=zion0(i)*deltazion(i)*(dion(i)+(1.d0-dion(i))*
     &     pahcfrac)
        endif
      enddo
      return
      end
      subroutine grainfinalise
c
      include 'cblocks.inc'
c
      real*8 dgrad
      integer*4 i, k, dtype, inl
c
c     deplete...
c
      pahfrac=0.d0
      if (pahmode.eq.1) then
        i=zmap(6)
        if (clinpah.eq.0) zion(i)=zion0(i)*deltazion(i)*(dion(i)+(1.d0-
     &   dion(i))*pahcfrac)
        pahfrac=zion0(i)*deltazion(i)*(1.d0-dion(i))*pahcfrac/pahnc
      endif
c
c  Calculate extinction & scattering cross-sections for use in totphot
c     & newdif
c
      do dtype=1,numtypes
        do inl=1,infph-1
          dcrosec(inl,dtype)=0.d0
          dfscacros(inl,dtype)=0.d0
          dbscacros(inl,dtype)=0.d0
        enddo
      enddo
      do dtype=1,numtypes
        do k=mindust(dtype),maxdust(dtype)
          do inl=1,infph-1
            dgrad=gradedge(k+1)-gradedge(k)
            dcrosec(inl,dtype)=dcrosec(inl,dtype)+absorp(inl,k,dtype)*
     &       dustsig(k,dtype)*dgrad
c     &        (absorp(inl,k,dtype)+scatter(inl,k,dtype))
c              dfscacros(inl,dtype) = dfscacros(inl,dtype)
c     &             + 0.5d0*(1.d0+gcos(inl,k,dtype))*
c     &             scatter(inl,k,dtype)*dustsig(k,dtype)*dgrad
c              dbscacros(inl,dtype) = dbscacros(inl,dtype)
c     &             + 0.5d0*(1.d0-gcos(inl,k,dtype))*
c     &             scatter(inl,k,dtype)*dustsig(k,dtype)*dgrad
          enddo
        enddo
      enddo
c
c     Calculate visual extinction (per H column) (V=5470A=2.27eV)
c
      do inl=1,infph-1
        if ((photev(inl).lt.2.2666435).and.(photev(inl+1).gt.2.2666435))
     &    then
          dustav=2.5*dlog(dexp(1.d0))*(dcrosec(inl,1)+dcrosec(inl,2))
          goto 10
        endif
      enddo
   10 continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Recompute grains if composition changes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adjustgrains (t, qh, allowonandoff)
c
      include 'cblocks.inc'
c
      real*8 t, qh
      integer*4 allowonandoff
c
c      real*8 fzgas
c
      real*8 dustx, dustxq, x
c
      if (t.le.grainmaxtemp) then
        dustx=1.d0
      else
        x=dmax1((t/grainmaxtemp)-1.d0,0.d0)*(grainmaxtemp/
     &   grainscaletemp)
        x=dmin1(x,30.d0)
        dustx=dexp(-x)
        dustx=dmin1(dustx,1.d0)
        dustx=dmax1(dustx,0.d0)
        if (dustx.lt.1d-9) dustx=0.d0
      endif
c
      if (qh.le.grainmaxq) then
        dustxq=1.d0
      else
        x=dmax1((qh/grainmaxq)-1.d0,0.d0)*(grainmaxq/grainscaleq)
        x=dmin1(x,30.d0)
        dustxq=dexp(-x)
        dustxq=dmin1(dustxq,1.d0)
        dustxq=dmax1(dustxq,0.d0)
        if (dustxq.lt.1d-9) dustxq=0.d0
      endif
c
      dustx=dustx*dustxq
c
      if (dustx.ne.0.d0) then
        grainmode=1
      endif
      if ((allowonandoff.eq.1).and.(dustx.eq.0.d0)) then
        grainmode=0
      endif
c
      if (grainadjust.eq.1) then
        call graindepletegas (dustx)
        call solvegraink ()
        call grainfinalise ()
      endif
c
c  stop recomputing once fully dusty
c
      if ((dustx.eq.1.d0).and.(allowonandoff.eq.0)) grainadjust=0
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set grain parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine grainpar ()
c
c
      include 'cblocks.inc'
      include 'mpif.h'
c
c
      integer*4 i,inl,dtype
      character ilgg*4,qlimit*4
      real*8 qn,qh,uh,tlim
c
   10 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
   20  format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Dust Physics:'
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
   30  format(
     & '  Dust is enabled.')
   40  format(
     & '  Dust is currently disabled.')
   50  format(/'  Include dust calculations? (y/N) :',$)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     defaults
c
c     Standard Grain parameters
c
c
c     Bgrain: threshold parameter
c     Yinf: yield parameter
c     alh: accomodation parameter
c     segrain,spgrain: sticking for electrons and protons
c     siggrain: projected area of grains per H atom
c     graindens : grain  density in g/cm^-3
c
c
      siggrain(1)=1.5d-21
      siggrain(2)=1.5d-21
      bgrain=8.0d0
      yinf=0.05d0
      segrain=0.5d0
      spgrain=1.d0
      gheat=0.d0
      gcool=0.d0
      haccom=0.14d0
      graindens(1)=2.0d0
      amin(1)=0.005
      amax(1)=0.25
      amin(2)=0.01
      amax(2)=0.25
      galpha=-3.5d0
      pahnc=468
      mindust(1)=1
      maxdust(1)=dustbinmax
      mindust(2)=1
      maxdust(2)=dustbinmax
c
c if 1 then use the next two parameters to destroy dust
c
      graindestructmode=0
c allow grain mixture recalculations in adjustgrains routine
      grainadjust=0
c above this temp, grains are destroyed and gas undepleted
      grainmaxtemp=1.0d5
      grainscaletemp=1.0d4
      grainmaxq=cls
      grainscaleq=0.1d0*grainmaxq
c
      grainmode=0
      pahmode=0
      pahactive=0
      pahcfrac=0.d0
      pahfrac=0.d0
      irmode=0
      irtemp=0
      do dtype=1,numtypes
        do inl=1,infph-1
          dcrosec(inl,dtype)=0.d0
          dfscacros(inl,dtype)=0.d0
          dbscacros(inl,dtype)=0.d0
        enddo
      enddo
c
      grainghr=0.d0
      graindgr=0.d0
      grainpgr=0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write (*,20)
      if (grainmode.eq.1) then
        write (*,30)
      else
        write (*,40)
      endif
      write (*,50)
c
      if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.ne.'Y') ilgg='N'
c
      if (ilgg.eq.'Y') grainmode=1
c
      if (grainmode.eq.1) then
c
        call depcha
c
   60     format(/,' Allow grain destruction? (y/N) :',$)
        write (*,60)
        if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
        ilgg=ilgg(1:1)
        write (*,*)
c
        if (ilgg.eq.'y') ilgg='Y'
        if (ilgg.ne.'Y') ilgg='N'
        if (ilgg.eq.'Y') graindestructmode=1
c
        if (graindestructmode.eq.1) then
c
          grainadjust=1
c
   70     format(/' Choose the grain photoionisation limit '/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::'/
     &        '  Upper ionisation parameter in terms of:'/
     &        '    Q   :  QHDN [nion/nH =',0pf7.4,']',/
     &        '    H   :  QHDH',/
     &        '    U   :  U(H)',/
     &        ' :: ',$)
   80     write (*,70) zen
          if (taskid.eq.0) read (*,'(a)') qlimit
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(qlimit,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
          qlimit=qlimit(1:1)
          write (*,*)
c
          if (qlimit.eq.'q') qlimit='Q'
          if (qlimit.eq.'h') qlimit='H'
          if (qlimit.eq.'u') qlimit='U'
c
          if ((qlimit.ne.'Q').and.(qlimit.ne.'H').and.(qlimit.ne.'U'))
     &     goto 80
c
          if (qlimit.eq.'Q') then
   90    format(' Give QHDN dust limit (<=100 as log [10.500] ): ',$)
            write (*,90)
            if (taskid.eq.0) read (*,*) qn
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(qn,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
            write (*,*)
            if (qn.le.100.d0) qn=10.d0**qn
            qh=qn/zen
            grainmaxq=qh
            grainscaleq=0.1d0*grainmaxq
          endif
c
          if (qlimit.eq.'H') then
  100    format(' Give QHDH dust limit (<=100 as log : [10.500] ): ',$)
            write (*,100)
            if (taskid.eq.0) read (*,*) qh
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(qh,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
            write (*,*)
            if (qh.le.100.d0) qh=10.d0**qh
            grainmaxq=qh
            grainscaleq=0.1d0*grainmaxq
          endif
c
          if (qlimit.eq.'U') then
  110    format(' Give U(H) dust limit (<=100 as log : [0.00] ): ',$)
            write (*,110)
            if (taskid.eq.0) read (*,*) uh
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(uh,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)             
            write (*,*)
            if (uh.le.100.d0) uh=10.d0**uh
            grainmaxq=uh*cls
            grainscaleq=0.1d0*grainmaxq
          endif
c
  120    format(' Give dust temperature limit (<=10 as log : [5.00] ): '
     &     ,$)
          write (*,120)
          if (taskid.eq.0) read (*,*) tlim
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(tlim,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)           
          write (*,*)
          if (tlim.le.10.d0) tlim=10.d0**tlim
          grainmaxtemp=tlim
          grainscaletemp=0.1d0*grainmaxtemp
        endif
c
c
c     Initialize anyway for safety
c
        numtypes=2
c
c     Determine grain distribution
c
  130   write (*,140)
  140    format(/' Choose the grain distribution model '/
     &        ' ::::::::::::::::::::::::::::::::::::::'/
     &        '    P   :  Powerlaw     N(a) = k a^alpha'/
     &        '    M   :  MRN distribution'/
     &        '    S   :  Grain Shattering Profile'/
     &        '   (Note: need full DUSTDATA for P or S)'/
     &        /' :: ',$)
c     &        '    E   :  Exponential  N(a) = k exp(a/alpha)'/
      if (taskid.eq.0) read (*,'(a)') gdist
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(gdist,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)        
        gdist=gdist(1:1)
        write (*,*)
        if (gdist.eq.'p') gdist='P'
        if (gdist.eq.'m') gdist='M'
        if (gdist.eq.'s') gdist='S'
c
        if ((gdist.ne.'P').and.(gdist.ne.'M').and.(gdist.ne.'S')) goto
     &   130
c
        if (gdist.eq.'M') then
c
c     MRN dist N(a) = k a^-3.5 C amin=50A Sil amin=100A,amax=2500 A
c
          galpha=-3.5d0
          amin(1)=0.0050d0
          amin(2)=0.0100d0
          amax(1)=0.2500d0
          amax(2)=0.2500d0
c
c  convert to cm
c
          amin(1)=amin(1)*1.d-4
          amax(1)=amax(1)*1.d-4
          amin(2)=amin(2)*1.d-4
          amax(2)=amax(2)*1.d-4
        else if (gdist.eq.'P') then
c
c      Power - Law N(a) = A^alpha
c
          write (*,150)
  150     format(/,' Enter alpha (MRN=-3.5)::',$)
          if (taskid.eq.0) read (*,*) galpha
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(galpha,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
c     Enter min and max grain radii (within 1e-3, 10 mu)
c     then find corresponding grain radius bins
c
c     Graphite first (type=1)
c
  160     format(/' Give graphite grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  170     write (*,160)
          if (taskid.eq.0) read (*,*) amin(1)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amin(1),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amin(1).le.0.d0) goto 170
          if (amin(1).ge.10.d0) amin(1)=amin(1)/1d4
          if ((amin(1).lt.1.0d-3).or.(amin(1).gt.10.d0)) goto 170
c
  180     format(/' Give graphite grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  190     write (*,180)
          if (taskid.eq.0) read (*,*) amax(1)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amax(1),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amax(1).le.0.d0) goto 190
          if (amax(1).ge.10.d0) amax(1)=amax(1)/1d4
          if ((amax(1).lt.1.0d-3).or.(amax(1).gt.10.d0)) goto 190
          if (amax(1).lt.amin(1)) goto 170
c
c     Then Silicate (type=2)
c
  200        format(/' Give silicate grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  210     write (*,200)
          if (taskid.eq.0) read (*,*) amin(2)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amin(2),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amin(2).le.0.d0) goto 210
          if (amin(2).ge.10.d0) amin(2)=amin(2)/1d4
          if ((amin(2).lt.1.0d-3).or.(amin(2).gt.10.d0)) goto 210
c
  220        format(/' Give silicate grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  230     write (*,220)
          if (taskid.eq.0) read (*,*) amax(2)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amax(2),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amax(2).le.0.d0) goto 230
          if (amax(2).ge.10.d0) amax(2)=amax(2)/1d4
          if ((amax(2).lt.1.0d-3).or.(amax(2).gt.10.d0)) goto 230
          if (amax(2).lt.amin(2)) goto 210
c
c  convert to cm
c
          amin(1)=amin(1)*1.d-4
          amax(1)=amax(1)*1.d-4
          amin(2)=amin(2)*1.d-4
          amax(2)=amax(2)*1.d-4
c
        else if (gdist.eq.'S') then
c
c      Grain shattering distribution based on Jones et al 1996
c       with additional exponential cutoffs
c      N(a) = k * a^-3.3 * (exp[-(a/amin)^-3])/(exp[(a/amax)^3])
c
  240    format(/,' Shattered grain distribution:',/
     &   ,'N(a) = k * a^alpha * (exp[-(a/amin)^-3])/(exp[(a/amax)])',/)
          write (*,240)
c     Enter min and max grain radii (within 1e-3, 10 mu)
c     then find corresponding grain radius bins
c
          write (*,250)
  250     format(/,' Enter alpha (Jones=-3.3, MRN=-3.5)::',$)
          if (taskid.eq.0) read (*,*) galpha
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(galpha,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
c     Graphite first (type=1)
c
  260        format(/' Give graphite grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  270     write (*,260)
          if (taskid.eq.0) read (*,*) amin(1)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amin(1),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
c
          if (amin(1).le.0.d0) goto 270
          if (amin(1).ge.10.d0) amin(1)=amin(1)/1d4
          if ((amin(1).lt.1.0d-3).or.(amin(1).gt.10.d0)) goto 270
c
  280        format(/' Give graphite grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  290     write (*,280)
          if (taskid.eq.0) read (*,*) amax(1)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amax(1),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amax(1).le.0.d0) goto 290
          if (amax(1).ge.10.d0) amax(1)=amax(1)/1d4
          if ((amax(1).lt.1.0d-3).or.(amax(1).gt.10.d0)) goto 290
          if (amax(1).lt.amin(1)) goto 270
c
c     Then Silicate (type=2)
c
  300        format(/' Give silicate grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  310     write (*,300)
          if (taskid.eq.0) read (*,*) amin(2)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amin(2),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)            
          write (*,*)
c
          if (amin(2).le.0.d0) goto 310
          if (amin(2).ge.10.d0) amin(2)=amin(2)/1d4
          if ((amin(2).lt.1.0d-3).or.(amin(2).gt.10.d0)) goto 310
c
  320        format(/' Give silicate grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
  330     write (*,320)
          if (taskid.eq.0) read (*,*) amax(2)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(amax(2),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)  
          write (*,*)
c
          if (amax(2).le.0.d0) goto 330
          if (amax(2).ge.10.d0) amax(2)=amax(2)/1d4
          if ((amax(2).lt.1.0d-3).or.(amax(2).gt.10.d0)) goto 330
          if (amax(2).lt.amin(2)) goto 310
c
c  convert to cm
c
          amin(1)=amin(1)*1.d-4
          amax(1)=amax(1)*1.d-4
          amin(2)=amin(2)*1.d-4
          amax(2)=amax(2)*1.d-4
        else
c
c  Wrong Grain size distribution
c
          goto 130
        endif
c
c  get edges
c
        if ((gdist.eq.'M').or.(gdist.eq.'P')) then
          mindust(1)=1
          maxdust(1)=dustbinmax
          do i=1,dustbinmax
            if ((amin(1).ge.gradedge(i)).and.(amin(1).lt.gradedge(i+1)))
     &        then
              mindust(1)=i
              amin(1)=grainrad(i)
            endif
            if ((amax(1).gt.gradedge(i)).and.(amax(1).le.gradedge(i+1)))
     &        then
              maxdust(1)=i
              amax(1)=grainrad(i)
            endif
          enddo
          mindust(2)=1
          maxdust(2)=dustbinmax
          do i=1,dustbinmax
            if ((amin(2).ge.gradedge(i)).and.(amin(2).lt.gradedge(i+1)))
     &        then
              mindust(2)=i
              amin(2)=grainrad(i)
            endif
            if ((amax(2).gt.gradedge(i)).and.(amax(2).le.gradedge(i+1)))
     &        then
              maxdust(2)=i
              amax(2)=grainrad(i)
            endif
          enddo
        endif
        write (*,340) amin(1)*1d4
        write (*,350) amax(1)*1d4
  340    format(' ***************************',/,
     &      '  Graphite grain min radius (mu): ',e11.4)
  350    format('  Graphite grain max radius (mu): ',e11.4,/,
     &      ' ***************************')
        write (*,360) amin(2)*1d4
        write (*,370) amax(2)*1d4
  360    format(' ***************************',/,
     &      '  Silicate grain min radius (mu): ',e11.4)
  370    format('  Silicate grain max radius mu): ',e11.4,/,
     &      ' ***************************')
c
c      Photoelectric parameters
c
        bgrain=8.0d0
        yinf=0.05d0
c
        if ((gdist.eq.'M')) then
          graindens(1)=1.85d0
          graindens(2)=2.5d0
          goto 430
        endif
c
  380   write (*,390)
  390    format(/' Give graphite grain density (~1.85 g/cm^3) :',$)
c
        if (taskid.eq.0) read (*,*) graindens(1)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(graindens(1),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)          
        write (*,*)
        if (graindens(1).le.0.d0) goto 380
c
  400   write (*,410)
  410    format(/' Give silicate grain density (~2.5 g/cm^3) :',$)
        if (taskid.eq.0) read (*,*) graindens(2)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(graindens(2),1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr) 
        write (*,*)
        if (graindens(2).le.0.d0) goto 400
c
c     Brent: Still need to be changed to distribution
c
c 200        format(/' Give grain photoelectric ',
c     &           'threshold(~8 eV [5-13]) :',$)
c 190        write(*,200)
c         read(*,*) Bgrain
c
c         if (Bgrain.lt.5.d0) goto 190
c         if (Bgrain.gt.13.d0) goto 190
c
c 211        format(/' Give grain photoelectric ',
c     &           'yield (~0.1 at inf.) :',$)
c 210     write(*,211)
c         read(*,*) Yinf
c
c         if (Yinf.lt.0.d0) goto 210
c         if (Yinf.gt.1.d0) goto 210
c
c     get PAHs at this point so carbon depletion is consistent
c
        clinpah=0
        pahmode=0
        pahactive=0
c
  420    format(/' Include PAH molecules? (y/n) :',$)
  430   write (*,420)
        if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)        
        ilgg=ilgg(1:1)
        write (*,*)
c
        if (ilgg.eq.'y') ilgg='Y'
        if (ilgg.eq.'Y') pahmode=1
c
        if (pahmode.eq.1) then
c
          pahnc=468
c
  440      format(/' Give fraction of Carbon Dust',
     &            ' Depletion in PAHs :',$)
  450     write (*,440)
          if (taskid.eq.0) read (*,*) pahcfrac
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(pahcfrac,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)           
          write (*,*)
c
          if (pahcfrac.lt.0.d0) goto 450
          if (pahcfrac.gt.1.d0) goto 450
          if (pahcfrac.eq.0.d0) pahmode=0
          if (pahmode.eq.1) then
c
  460       write (*,470)
  470     format(/,' Choose PAH switch on;',/
     &     ,' :::::::::::::::::::::::::::::::::::::::::::::::::::',/
     &     ,'   H :  Habing photodissociation parameter < Value',/
     &     ,'   I :  PAH ionizing radiation/ local ISRF ratio < Value',/
     &     ,'   Q :  QHDH < Value',/
     &     ,' :: ',$)
        if (taskid.eq.0) read (*,*) ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)              
            ilgg=ilgg(1:1)
            write (*,*)
            if (ilgg.eq.'h') ilgg='H'
            if (ilgg.eq.'i') ilgg='I'
            if (ilgg.eq.'q') ilgg='Q'
c
            if ((ilgg.ne.'H').and.(ilgg.ne.'I').and.(ilgg.ne.'Q')) goto
     &       460
c
            pahend=ilgg(1:1)
c
  480       write (*,490)
  490       format(/
     &      ' Give PAH switch on Value (<=100 as log : [4.00] ):',$)
          if (taskid.eq.0) read (*,*) pahlimit
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(pahlimit,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)             
            if (pahlimit.le.100.d0) pahlimit=10.d0**pahlimit
            write (*,*)
c
            if (pahlimit.lt.0.d0) goto 480
c
c
c   Convert to #PAH grains per H atom
c
            pahfrac=zion(zmap(6))*(1.d0-dion(zmap(6)))*pahcfrac/pahnc
            write (*,'("PAH per H (ppm): ",1pg11.4)') pahfrac*1.d6
            write (*,500)
  500     format(' Do you wish graphite grains to be cospatial'
     &           ,' with PAHs?',/
     &           ,' (PAH destroyed = C grain destroyed):',$)
            if (taskid.eq.0) read (*,*) ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)            
            ilgg=ilgg(1:1)
            write (*,*)
            if (ilgg.eq.'y') ilgg='Y'
            if (ilgg.eq.'Y') clinpah=1
c
          endif
        endif
c
        write (*,510)
  510 format(/' Evaluate dust temperatures and IR flux? (P5/P6 only)',/,
     &  '  (WARNING: will slow down computing time) (y/N):',$)
        if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)        
        ilgg=ilgg(1:1)
        write (*,*)
        if (ilgg.eq.'y') ilgg='Y'
c
        if (ilgg.eq.'Y') then
  520     write (*,530)
  530    format(/,' Which IR model? ',/
     &      ,' ::::::::::::::::::::::::::::::::::::::',/
     &      ,'    Q :  Quick (50 bins, ~10% accurate IR)',/
     &      ,'    I :  Intermediate (100 bins, ~2% accurate',/
     &      ,'    S :  Slow (400 bins, accurate Temp. distributions)',/
     &      ,/,' :: ',$)
        if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)     
          ilgg=ilgg(1:1)
          if ((ilgg.eq.'q').or.(ilgg.eq.'Q')) then
            irmode=1
          else if ((ilgg.eq.'i').or.(ilgg.eq.'I')) then
            irmode=2
          else if ((ilgg.eq.'s').or.(ilgg.eq.'S')) then
            irmode=3
          else
            goto 520
          endif
c       if (expertmode.gt.0) then
          write (*,540)
  540    format(/,' Output dust temperature data (large file) (y/n):',$)
        if (taskid.eq.0) read (*,'(a)') ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)     
          ilgg=ilgg(1:1)
          if (ilgg.eq.'y') ilgg='Y'
          if (ilgg.eq.'Y') then
            irtemp=1
          else
            irtemp=0
          endif
c       endif
        endif
c
        call graindepletegas (1.d0)
        call solvegraink ()
        call grainfinalise ()
c
  550    format(//'  Projected dust area/H atom (graphite,silicate):'/
     &            ' ****************************'/
     &            '  ',1pg12.4,1x,' (cm^2)',1x,1pg12.4,1x,' (cm^2)' )
        write (*,550) siggrain(1),siggrain(2)
c
  560    format( /'  Composition Ratios:  '/
     &            ' ****************************'/
     &            '   Gas/H     :',1pg12.4/
     &            '   Dust/Gas  :',1pg12.4/
     &            '   PAH/Gas   :',1pg12.4/
     &            ' ****************************'/
     &           ' (ratios by mass) ')
c
        write (*,560) grainghr,graindgr,grainpgr
c
c
c   Init fixed parameters for new IR temp calcs
c   in dustemp.f
c
        call dustinit ()
c
      endif
c
      write (*,10)
c
      return
c
      end
