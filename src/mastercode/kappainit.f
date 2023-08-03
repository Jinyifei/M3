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
c     Subroutine to setup electron distributions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kappainit ()
c
      include 'cblocks.inc'
      include 'mpif.h'
c
      integer*4 k
      real*8 kap, x, delta
      character ilgg*4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Non-thermal Kappa excitation,will move to init routine when debugged
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kappaidx=9
      kappaa=1.d0
      kappab=0.d0
      usekappa=.false.
      usekappainterp=.false.
      kappa=1.0d99
c
   10 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
   20  format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Electron Energy Distributions:'
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
   30 format(/'  Maxwellian thermal distributions are being used.'/
     &       /'  Use Kappa electron distributions ? (y/N) : ',$)
   40 format (a)
   50 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Select a Kappa value'/
     & '   (2.0 - 1000.0, outside this range disabled)'/
     & '   (2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0, 100.0 are exact) :',$)
      write (*,20)
      write (*,30)
      if (taskid.eq.0) read (*,40) ilgg
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'Y') usekappa=.true.
      if (usekappa) then
        write (*,50)
        if (taskid.eq.0) read (*,*) kap
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(kap,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
        if ((kap.ge.2.d0).and.(kap.le.1000.d0)) then
          kappa=kap
          do k=1,(nkappas-1)
            if (kap.ge.kappas(k)) kappaidx=k
          enddo
          delta=(1.d0/kappas(kappaidx))-(1.d0/kappas(kappaidx+1))
c 2> 1so1/1 >1/2
          x=(1.d0/kappa)-(1.d0/kappas(kappaidx+1))
          kappaa=x/delta
          kappab=1.d0-kappaa
          if (kappaa.lt.0.995d0) usekappainterp=.true.
c
         write(*,*) kappa, kappaA, kappaB, kappaidx, useKappaInterp
c
c linear interp coeffs,
c         y = kappaA*spline(index1) + kappaB*spline(index2)
c
        else
          usekappa=.false.
        endif
      endif
c
      write (*,10)
c
      return
      end
