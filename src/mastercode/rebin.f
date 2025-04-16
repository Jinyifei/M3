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
c Conservative powerlaw interpolation of highly non-linear data
c rebinning method.
c
c First scales data from file units to ev and hnu, then
c renins and converts to Inu = 4*Hnu.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readrebin (fname, coordstype, rebinnedsrc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Use powerlaws between src points, alpha, and interval
c analytic integrals, sflux, to fill destination vector with
c integrals of the src, over bins defined by bin edges in the
c destination = then average by destination bin width to get
c final, conservatively averaged, bin data.
c
c find piecewise powerlaw integrals between src points.
c from low -> high energy
c
c Use log-log coordinates to get alpha = slope,
c and linear coordinates to do analytical
c integrals between source points.
c
c Input with flux densities at input energies, to average
c fluxes over bin with left edge energies.
c
c Use analytical powerlaw integrals, or linear interpolation
c if magnitude of alpha is too great
c  set some constants depending on coordstype
c
c  0 = Kurucz nm-Hnu two column atmosphere file
c  1 = SB99 spectrum 1 luminosity file
c  2 = AGN Model file NLS+BLS v1
c  3 = variable file, two columns.
c  4 = AGN Model file NLS+BLS v2
c  5 = MAPPINGS .sou file with different binning
c  6 = TLUSTY Hz-Hnu two column atmosphere file
c  7 = WMBASIC A-Hnu two column atmosphere file
c  8 =  Hz Hnu two column atmosphere file.
c  9 = CSPN Rauch TNMAP HNi/HCa grid lambda flambda
c      two column atmosphere file
c  11 = SLUG integrated galaxy spectrum  luminosity file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'mpif.h'
c
c arguments
c
      character*256 fname
      integer*4 coordstype
      real*8 rebinnedsrc(mxinfph)
c
c  locals
c
      logical iexi
      integer*4 luin
c
      character*256 caract
c
      integer*4 i,j,l,m
      integer*4 idx0, idx1, idx2, idx3
      integer*4 nfluxes, nheadlines, nskiplines
      integer*4 lowbin, hibin
      integer*4 nmodels, modelidx
      integer*4 nfluxes2,souid2
      character*4 xtype
      character*4 ytype
c
      real*8 agnkev, agnwhidth, agbflux, agnshift
      real*8 agncomp0, agncomp1, agncomp2
      real*8 xbump, xinter, xhigh
      real*8 ac,wid0,wid1,wev0,blum0,blum1
      real*8 wavescale
      real*8 freqscale
c
      real*8 sb99time, sb99total
      real*8 sb99wave, sb99star, sb99nebula
      real*8 slugtime, slugtotal, slugwave
      real*8 freq, lam, lamcm, time(999)
c
      real*8 ufrac, lfrac, mean, f
      real*8 delta, logdelta, swap
      real*8 a0, a1
      real*8 e0, e1
      real*8 w0, w1, scale, wid, sum, total
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Clear array first
c
      blum0=0.d0
      blum1=0.d0
      do i=1,infph
        rebinnedsrc(i)=0.d0
      enddo
c
      luin=17
      inquire (file=fname,exist=iexi)
c
      if (iexi) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fname exists
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set some constants depending on coordstype
c
c  0 = Kurucz nm-Hnu two column atmosphere file
c  1 = SB99 spectrum 1 luminosity file
c  2 = AGN Model file NLS+BLS v1
c  3 = variable file, two columns.
c  4 = AGN Model file NLS+BLS v2
c  5 = MAPPINGS .sou file with different binning
c  6 = TLUSTY Hz-Hnu two column atmosphere file
c  7 = WMBASIC A-Hnu two column atmosphere file
c  8 =  Hz Hnu two column atmosphere file.
c  9 = CSPN Rauch TNMAP HNi/HCa grid lambda flambda
c      two column atmosphere file
c  11 = SLUG integrated galaxy spectrum  luminosity file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        wavescale=1.d0
        freqscale=1.d0
        nfluxes=1
        nheadlines=0
c
        if (coordstype.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 0, Kurucz nm-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Kurucz standard nm - Hnu units
c nm  Hnu ergs/cm^2/s/Hz/sr in 1/4pi units
c  4.0*Hnu = Inu
c Columns 1 and 2
          wavescale=10.d0
          freqscale=1.d0
          nfluxes=1221
          nheadlines=4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 0, Kurucz nm-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.1) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 1, SB99 spectrum 1 luminosity file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c SB99 spectrum1 standard Angstroms - L_lambda units, in
c linear lambda, log L TIME [YR]    WAVELENGTH [A]   LOG
c TOTAL  LOG STELLAR  LOG NEBULAR  [ERG/SEC/A] Columns 2 and 3
c
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=1221
          nheadlines=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 1, SB99 spectrum 1 luminosity file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.6) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 6, TLUSTY Hz-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TLUSTY standard Hz - Hnu units
c Hz  Hnu ergs/cm^2/s/Hz/sr in 1/4pi units
c  4.0*Hnu = Inu
c Columns 1 and 2
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=19998
          nheadlines=4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 6, TLUSTY Hz-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.7) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 7, WMBASIC A-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c WMBASIC standard Hz - Hnu units
c A  Hnu ergs/cm^2/s/Hz/sr in 1/4pi units
c  4.0*Hnu = Inu
c Columns 1 and 2
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=0
          nheadlines=4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 7, WMBASIC A-Hnu two column atmosphere file..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.9) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN Rauch HNi/HCa grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Rauch  nonstandard Lammbda (A) - Flambda (ergs/cm^2/s/cm/pi)
c Flambda  ergs/cm^2/s/cm/pi in 1/pi units
c Columns 1 and 2
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=0
          nheadlines=37
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN Rauch HNi/HCa grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.11) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c coordstype = 11, SLUG itegrated_spec single/multiple step file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c SLUG spectrum  Angstroms - L_lambda units, in
c Three columnm file with extinct off and nebula off
c Time          Wavelength    L_lambda
c(yr)          (Angstrom)     (erg/s/A)
c
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=0
          nheadlines=3
          nskiplines=1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c coordstype = 11, SLUG itegrated spectrum single/multiple step file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.2) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 2, AGN Model file NLS+BLS v1 1000 bins
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Every file name shows the parameter values, e.g.
c 'M1e7_E0.1_NLS1.txt' means M=1.e7, Eddington ratio = 0.1, for BLS1;
c
c In every file:
c Column 1: center of energy bin in KeV;
c Column 2: half width of energy bin in KeV;
c Column 3: total model flux in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 4: flux of disc component in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 5: flux of low  temperature Comptonisation component
c            in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 6: flux of high temperature Comptonisation component
c           in keV*keV/cm^2/s/keV, i.e. E*F(E);
c
c The SEDs are all in the observer's frame with redshift = 0.1
c (i.e. distance = 413.5Mpc).
c scale energy by (1+z) = 1.1, fluxes by (1+z)^2 = 1.21
c
          agnshift=0.1d0
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=1000
          nheadlines=3
          xbump=1.d0
          xinter=1.d0
          xhigh=1.d0
          write (*,10)
   10 format(//
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  AGN Component Fractions :'/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    Thermal Bump        (0 - 1, default 1) :',$)
c          read (*,*) xbump
      if (taskid.eq.0) read (*,*) xbump
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xbump,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)    
c
          write (*,*)
          if (xbump.lt.0.0d0) xbump=0.0d0
          if (xbump.gt.1.0d0) xbump=1.0d0
          write (*,20)
   20 format(  '   Intermediate Compton (0 - 1, default 1) :',$)
c          read (*,*) xinter
      if (taskid.eq.0) read (*,*) xinter
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xinter,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c      
          write (*,*)
          if (xinter.lt.0.0d0) xinter=0.0d0
          if (xinter.gt.1.0d0) xinter=1.0d0
          write (*,30)
   30 format(  '   High Energy Nonthermal (0 - 1, default 1) :',$)
c          read (*,*) xhigh
      if (taskid.eq.0) read (*,*) xhigh
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xhigh,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)   
c      
          write (*,*)
          if (xhigh.lt.0.0d0) xhigh=0.0d0
          if (xhigh.gt.1.0d0) xhigh=1.0d0
c      write(*,*) xBump, xInter, xHigh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 2, AGN Model file NLS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.4) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 4, AGN Model file NLS+BLS v2 3000 bins
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Every file name shows the parameter values, e.g.
c 'M1e7_E0.1_BLS1.txt' means M=1.e7, Eddington ratio = 0.1, for BLS1;
c
c In every file:
c Column 1: center of energy bin in KeV;
c Column 2: half width of energy bin in KeV;
c Column 3: total model flux in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 4: flux of disc component in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 5: flux of low  temperature Comptonisation component
c            in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 6: flux of high temperature Comptonisation component
c           in keV*keV/cm^2/s/keV, i.e. E*F(E);
c
c The SEDs are all in the observer's frame with redshift = 0.0
c distance 1.5e18cm (1e4 Rg at 1e9)
c
          agnshift=0.1d0
          wavescale=1.d0
          freqscale=1.d0
          nfluxes=3000
          nheadlines=8
          xbump=1.d0
          xinter=1.d0
          xhigh=1.d0
          write (*,40)
   40 format(//
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  AGN Component Fractions :'/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '   Thermal Bump         (0 - 1, default 1) :',$)
c          read (*,*) xbump
      if (taskid.eq.0) read (*,*) xbump
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xbump,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c          
          write (*,*)
          if (xbump.lt.0.0d0) xbump=0.0d0
          if (xbump.gt.1.0d0) xbump=1.0d0
          write (*,50)
   50 format(  '   Intermediate Compton (0 - 1, default 1) :',$)
c          read (*,*) xinter
      if (taskid.eq.0) read (*,*) xinter
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xinter,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
          write (*,*)
          if (xinter.lt.0.0d0) xinter=0.0d0
          if (xinter.gt.1.0d0) xinter=1.0d0
          write (*,60)
   60 format(  '   High Energy Nonthermal (0 - 1, default 1) :',$)
c          read (*,*) xhigh
      if (taskid.eq.0) read (*,*) xhigh
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(xhigh,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c          
          write (*,*)
          if (xhigh.lt.0.0d0) xhigh=0.0d0
          if (xhigh.gt.1.0d0) xhigh=1.0d0
c      write(*,*) xBump, xInter, xHigh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 4, AGN Model file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Conservative powerlaw interpolation of highly non-linear data
c rebinning method.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Use powerlaws between src points, alpha, and interval
c analytic integrals, sflux, to fill destination vector with
c integrals of the src, over bins defined by bin edges in the
c destination = then average by destination bin width to get
c final, conservatively averaged, bin data.
c
c find piecewise powerlaw integrals between src points.
c from low -> high energy
c
c Use log-log coordinates to get alpha = slope,
c and linear coordinates to do analytical
c integrals between source points.
c
c Use analytical powerlaw integrals, or linear interpolation
c if magnitude of alpha is too great
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy, ( to make the integrals less
c confusing) Between each pair of points, use log-log coords
c to get the power law slope in the interval, and
c analytically integrate the fluxes over the interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (coordstype.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 0, Kurucz nm-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write (*,*) 'Reading ATLAS9 nm-Hnu Atmosphere File ',fname
c          write (*,*) 'N Header Lines:',nheadlines
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy, ( to make the integrals less
c confusing) Between each pair of points, use log-log coords
c to get the power law slope in the interval, and
c analytically integrate the fluxes over the interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do i=nfluxes,1,-1
            read (luin,*) stlwave(i),stlhnu(i)
c      cvt to eV from nm, LmeV is a std MAPPINGS constant
            stlwev(i)=lmev/(wavescale*stlwave(i))
          enddo
          close (luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 0, Kurucz nm-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (coordstype.eq.6) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 6, TLUSTY Hz-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write (*,*) 'Reading TLUSTY Hz-Hnu Atmosphere File ',fname
c          write (*,*) 'N Header Lines:',nheadlines
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy, ( to make the integrals less
c confusing) Between each pair of points, use log-log coords
c to get the power law slope in the interval, and
c analytically integrate the fluxes over the interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c reverse, as TLUSTY has high freq first
c
          do i=nfluxes,1,-1
            read (luin,*) stlwave(i),stlhnu(i)
c      cvt to eV from Hz, plk and ev are std MAPPINGS constant
            stlwev(i)=plk*stlwave(i)/ev
          enddo
          close (luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 6, TLUSTY Hz-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.1) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coords type = 1, SB99 spectrum 1 file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading SB99 spectrum1 File',fname
          open (luin,file=fname,status='OLD')
          do i=1,99999
            read (luin,'(a)',end=70) caract
            nheadlines=i
            idx0=index(caract,'TIME')
            idx1=index(caract,'[YR]')
            idx2=index(caract,'WAVELENGTH')
            idx3=index(caract,'[A]')
            if ((idx0.gt.0).and.(idx1.gt.idx0).and.(idx2.gt.idx1)
     &       .and.(idx3.gt.idx2)) goto 80
          enddo
   70     continue
          write (*,*) 'ERROR: SB99 file invalid header.'
          close (luin)
          return
   80     write (*,*) 'SB99 File: Found N Header Lines:',nheadlines
          close (luin)
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  scan SB99 spectrum 1 file for time steps
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nmodels=0
          do m=1,999
            do i=nfluxes,1,-1
              read (luin,*,end=90) sb99time,sb99wave,sb99total,sb99star,
     &         sb99nebula
            enddo
            nmodels=nmodels+1
            time(nmodels)=sb99time
          enddo
   90     continue
c
          close (luin)
c
          if (nmodels.lt.1) then
            write (*,*) ' ERROR:No complete spectra found'
            return
          endif
c
  100         format(/'  Available Model Time/s (yr), choose an index:'/
     & ' :::::::::::::::::::::::::::::::::::::s:::::::::::::::::::::::')
  110         format(4(2x,i2,':',1pg10.3))
  120         format(' :: ',$)
          write (*,100)
          write (*,110) (i,time(i),i=1,nmodels)
          write (*,120)
          if (taskid.eq.0) read (*,*) modelidx
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(modelidx,1,MPI_INTEGER4,0,MPI_COMM_WORLD,taskerr)          
          write (*,*)
c
          if (modelidx<1) modelidx=1
          if (modelidx>nmodels) modelidx=nmodels
c
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c skip models up to chosen one
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (modelidx.gt.1) then
            do m=1,(modelidx-1)
              do i=nfluxes,1,-1
                read (luin,*) sb99time,sb99wave,sb99total,sb99star,
     &           sb99nebula
              enddo
            enddo
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy, ( to make the integrals less confusing)
c Between each pair of points, use log-log coords to get the power law slope
c in the interval, and analytically integrate the fluxes over the interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do i=nfluxes,1,-1
            read (luin,*) sb99time,stlwave(i),sb99total,stlhnu(i),
     &       sb99nebula
            if (sb99time.eq.time(modelidx)) then
c
c wave in 1e-8 cm, A, and flux in log-L-lam erg/s/A
c
              stlhnu(i)=10.d0**stlhnu(i)
c
c convert L_lam ergs/s/A -> Hnu erg/cm^2/s/Hz/sr, assuming a
c std radius = 1.0d10 area = 4pir^2 = 1.256637061435917e21
c cm^2 1/area = 7.957747154594766e-22 sr = 1/4pi, also ifpi =
c invers 4 pi
c
              freq=cls/(stlwave(i)*1.d-8)
              stlhnu(i)=7.957747154594766d-22*((stlhnu(i)*stlwave(i))/
     &         freq)*ifpi
              stlwev(i)=lmev/(wavescale*stlwave(i))
            else
              write (*,*) 'ERROR: SB99 file Time Mismatch Error,',' corr
     &upt file?'
              return
            endif
          enddo
c
          close (luin)
c
  130  format(//
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Read SB99 Spectrum at t=', 1pg10.3, 'yr',/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
          write (*,130) time(modelidx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 1, SB99 spectrum 1 file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.11) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coords type = 1, SB99 spectrum 1 file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading SLUG integrated_spect 3 col File',fname
          open (luin,file=fname,status='OLD')
          read (luin,'(a)') caract
          nheadlines=1
          idx0=index(caract,'Time')
          idx1=index(caract,'Wavelength')
          idx2=index(caract,'L_lambda')
          if ((idx0.gt.0).and.(idx1.gt.idx0).and.(idx2.gt.idx1)) goto
     &     140
          write (*,*) 'ERROR: SLUG file invalid header.'
          close (luin)
          return
  140     read (luin,'(a)') caract
          nheadlines=2
          idx0=index(caract,'(yr)')
          idx1=index(caract,'(Angstrom)')
          idx2=index(caract,'(erg/s/A)')
          if ((idx0.gt.0).and.(idx1.gt.idx0).and.(idx2.gt.idx1)) goto
     &     150
          write (*,*) 'ERROR: SLUG file invalid header.'
          close (luin)
          return
  150     read (luin,'(a)') caract
          nheadlines=3
          idx0=index(caract,'===========')
          if (idx0.gt.0) goto 160
          write (*,*) 'ERROR: SLUG file invalid header.'
          close (luin)
          return
  160     write (*,*) 'SLUG File: Found Header Lines:',nheadlines
          close (luin)
c
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(" ",a80)') caract
          enddo
c     skip line
          do i=1,nskiplines
            read (luin,'(a)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  scan SLUG integrated_spec file for time steps
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nmodels=0
          do m=1,999
            nfluxes=0
            do i=1,99999
              read (luin,*,err=170,end=180) slugtime,slugwave,slugtotal
              nfluxes=nfluxes+1
            enddo
  170       continue
            nmodels=nmodels+1
            time(nmodels)=slugtime
          enddo
          goto 190
  180     continue
          nmodels=nmodels+1
          time(nmodels)=slugtime
  190     continue
          close (luin)
          write (*,*) 'SLUG File: Found Flux Lines:',nfluxes
c
          if (nmodels.lt.1) then
            write (*,*) ' ERROR:No complete spectra found'
            return
          endif
c
  200   format(/' Available Model Time/s (yr), choose an index:'/
     & '::::::::::::::::::::::::::::::::::::s:::::::::::::::::::::::')
  210   format(4(2x,i2,':',1pg10.3))
  220   format(' :: ',$)
          write (*,200)
          write (*,210) (i,time(i),i=1,nmodels)
          write (*,220)
          read (*,*) modelidx
          write (*,*)
c
          if (modelidx<1) modelidx=1
          if (modelidx>nmodels) modelidx=nmodels
c
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c skip models up to chosen one
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (modelidx.gt.1) then
            do m=1,(modelidx-1)
c     skip line
              do i=1,nskiplines
                read (luin,'(a)') caract
                write (*,'(" ",a80)') caract
              enddo
              do i=nfluxes,1,-1
                read (luin,*) slugtime,slugwave,slugtotal
              enddo
            enddo
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy, ( to make the integrals less confusing)
c Between each pair of points, use log-log coords to get the power law slope
c in the interval, and analytically integrate the fluxes over the interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do i=1,nskiplines
            read (luin,'(a)') caract
            write (*,'(" ",a80)') caract
          enddo
          do i=nfluxes,1,-1
            read (luin,*) slugtime,stlwave(i),stlhnu(i)
            if (slugtime.eq.time(modelidx)) then
c
c wave in 1e-8 cm, A, and flux in Linear-L-lam erg/s/A
c
c convert L_lam ergs/s/A -> Hnu erg/cm^2/s/Hz/sr, assuming a
c std radius = 1.0d10 area = 4pir^2 = 1.256637061435917e21
c cm^2 1/area = 7.957747154594766e-22 sr = 1/4pi, also ifpi =
c invers 4 pi
c
              freq=cls/(stlwave(i)*1.d-8)
              stlhnu(i)=7.957747154594766d-22*((stlhnu(i)*stlwave(i))/
     &         freq)*ifpi
              stlwev(i)=lmev/(wavescale*stlwave(i))
            else
              write (*,*) 'ERROR: SLUG file Time Mismatch Error,',' corr
     &upt file?'
              return
            endif
          enddo
c
          close (luin)
c
  230  format(//
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Read SLUG Spectrum at t=', 1pg10.3, 'yr',/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
          write (*,230) time(modelidx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 1, SB99 spectrum 1 file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (coordstype.eq.2) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 2, AGN Model file NLS+BLS v1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading AGN SED v1 File',fname
          write (*,*) 'N Header Lines:',nheadlines
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(" ",a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy Column 1: center of energy bin
c in KeV; Column 2: half width of energy bin in KeV; Column
c 3: total model flux in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 4: flux of disc component in keV*keV/cm^2/s/keV,
c i.e. E*F(E); Column 5: flux of low  temperature
c Comptonisation component in keV*keV/cm^2/s/keV, i.e.
c E*F(E); Column 6: flux of high temperature Comptonisation
c component in keV*keV/cm^2/s/keV, i.e. E*F(E);
c
c The SEDs are all in the observer's frame with redshift =
c 0.1 (i.e. distance = 413.5Mpc, h = 72.5km/sMpc).
c library approx correction applied to recover ~L/Ledd
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do i=1,nfluxes
            read (luin,*) agnkev,agnwhidth,agbflux,agncomp0,agncomp1,
     &       agncomp2
            agbflux=xbump*agncomp0+xinter*agncomp1+xhigh*agncomp2
c          cvt to keV to eV
            stlwev(i)=agnkev*1.0d3
c          -> Hz
            freq=stlwev(i)*evplk
c          E.F(E) -> F(E) , keV/cm^2/s/keV
            stlhnu(i)=(agbflux/agnkev)
c          keV/cm^2/s/keV -> keV/cm^2/s/Hz
            stlhnu(i)=stlhnu(i)*(agnkev/freq)
c          keV/cm^2/s/Hz  -> ergs/cm^2/s/Hz
            stlhnu(i)=stlhnu(i)*1.0e3*ev
c          de-redshift energy bins
            stlwev(i)=stlwev(i)*1.1d0
c
c Fix apparent inconsistencies in tables to
c recover approx correct Edd. luminosities
c
            stlhnu(i)=ifpi*stlhnu(i)/1.1d0
c                 413.5Mpc sphere to 1.0e16cm sphere
            stlhnu(i)=stlhnu(i)*1.627206670138892d22
          enddo
          close (luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 2, AGN Model file NLS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
        if (coordstype.eq.4) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 4, AGN Model files v2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading AGN SED v2 File',fname
          write (*,*) 'N Header Lines:',nheadlines
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(" ",a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data and convert to energy on the x-axis,
c sorted from low-high energy Column 1: center of energy bin
c in KeV; Column 2: half width of energy bin in KeV; Column
c 3: total model flux in keV*keV/cm^2/s/keV, i.e. E*F(E);
c Column 4: flux of disc component in keV*keV/cm^2/s/keV,
c i.e. E*F(E); Column 5: flux of low  temperature
c Comptonisation component in keV*keV/cm^2/s/keV, i.e.
c E*F(E); Column 6: flux of high temperature Comptonisation
c component in keV*keV/cm^2/s/keV, i.e. E*F(E);
c
c These SEDs are all already in the observer's frame with redshift =
c 0.0 at 1.5e18cm
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do i=1,nfluxes
            read (luin,*) agnkev,agnwhidth,agbflux,agncomp0,agncomp1,
     &       agncomp2
            agbflux=xbump*agncomp0+xinter*agncomp1+xhigh*agncomp2
            stlwev(i)=agnkev*1.0d3
            freq=stlwev(i)*evplk
            stlhnu(i)=(agbflux/agnkev)
            stlhnu(i)=stlhnu(i)*(agnkev/freq)
            stlhnu(i)=ifpi*stlhnu(i)*1.0e3*ev
          enddo
          close (luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 4, AGN Model file v2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.3) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 3, variable file, two columns.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading generic two column flux file',fname
  240      format(/' Give number of lines to skip : ',$)
          write (*,240)
          read (*,*) nheadlines
          write (*,*)
          if (nheadlines.lt.0) nheadlines=0
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
c
  250    format(/'  Choose x-coordinate units : '/
     &  ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &  '    A  : Wavelength Angstroms'/
     &  '    B  : Wavelength nm'/
     &  '    C  : Frequency Hz'/
     &  '    D  : Energy ergs'/
     &  '    E  : Energy eV'/
     &  '    F  : Energy keV'/
     &  ' :: ',$)
          write (*,250)
          read (*,'(a)') xtype
          xtype=xtype(1:1)
          write (*,*)
          if (xtype.eq.'a') xtype='A'
          if (xtype.eq.'b') xtype='B'
          if (xtype.eq.'c') xtype='C'
          if (xtype.eq.'d') xtype='D'
          if (xtype.eq.'e') xtype='E'
          if (xtype.eq.'f') xtype='F'
c
          if ((xtype.ne.'A').and.(xtype.ne.'B').and.(xtype.ne.'C')
     &     .and.(xtype.ne.'D').and.(xtype.ne.'E').and.(xtype.ne.'F'))
     &     xtype='A'
c
  260  format(/'  Choose y-coordinate units : '/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    A  : Inu  (1/pi ergs/cm^2/s/Hz) '//
     & '    B  : Flam (ergs/cm^2/s/A)'/
     & '    C  : Fnu  (ergs/cm^2/s/Hz)'//
     & '    D  : Eddington, Hlam (ergs/cm^2/s/A/sr) '/
     & '    E  : Eddington, Hnu  (ergs/cm^2/s/Hz/sr) '//
     & '    F  : Llam (ergs/s/A),  r = 1.0d10cm'/
     & '    L  : Log10 Llam (ergs/s/A),  r = 1.0d10cm'/
     & '    G  : Lnu  (ergs/s/Hz), r = 1.0d10cm'/
     & '    H  : Fnu  (keV/cm^2/s/keV)'/
     & '    I  : Flam (ergs/cm^2/s/cm/pi)'//
     & ' :: ',$)
          write (*,260)
          read (*,'(a)') ytype
          call toup (ytype)
          ytype=ytype(1:1)
          write (*,*)
          if ((ytype.ne.'A').and.(ytype.ne.'B').and.(ytype.ne.'C')
     &     .and.(ytype.ne.'D').and.(ytype.ne.'E').and.(ytype.ne.'F')
     &     .and.(ytype.ne.'G').and.(ytype.ne.'H').and.(ytype.ne.'I')
     &     .and.(ytype.ne.'L')) ytype='A'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data straight, and count fluxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nfluxes=0
          do i=1,mxflxbins
            read (luin,*,end=270) stlwave(i),stlhnu(i)
            nfluxes=nfluxes+1
          enddo
  270     continue
          close (luin)
c
          if (nfluxes.lt.2) then
            write (*,*) 'ERROR: no fluxes found in ',fname
            return
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c convert and sort the data into increasing eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Angstroms
          if (xtype.eq.'A') then
            if (stlwave(nfluxes).gt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=lmev/stlwave(i)
            enddo
          endif
c nm
          if (xtype.eq.'B') then
            if (stlwave(nfluxes).gt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=lmev/(10.d0*stlwave(i))
           enddo
          endif
c Hz
          if (xtype.eq.'C') then
            if (stlwave(nfluxes).lt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=plk*stlwave(i)/ev
            enddo
          endif
c ergs
          if (xtype.eq.'D') then
            if (stlwave(nfluxes).lt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=stlwave(i)/ev
            enddo
          endif
c eV
          if (xtype.eq.'E') then
            if (stlwave(nfluxes).lt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=stlwave(i)
            enddo
          endif
c keV
          if (xtype.eq.'F') then
            if (stlwave(nfluxes).lt.stlwave(1)) then
c reverse
              do i=1,nfluxes
                j=nfluxes-i+1
                if (i.lt.j) then
                  swap=stlwave(j)
                  stlwave(j)=stlwave(i)
                  stlwave(i)=swap
                  swap=stlhnu(j)
                  stlhnu(j)=stlhnu(i)
                  stlhnu(i)=swap
                endif
              enddo
            endif
            do i=1,nfluxes
              stlwev(i)=stlwave(i)*1.0d3
            enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Inu -> Hnu
          if (ytype.eq.'A') then
            do i=1,nfluxes
              stlhnu(i)=0.25*stlhnu(i)
            enddo
          endif
c Flam -> Hnu
          if (ytype.eq.'B') then
            do i=1,nfluxes
              lam=lmev/stlwev(i)
              lamcm=1.0d-8*lam
              freq=cls/lamcm
              stlhnu(i)=stlhnu(i)*lam/(freq*fpi)
            enddo
          endif
c special flam for cspn spectra
          if (ytype.eq.'I') then
            do i=1,nfluxes
              lam=lmev/stlwev(i)
              lamcm=1.0d-8*lam
              freq=cls/lamcm
              stlhnu(i)=0.25*stlhnu(i)*lamcm/freq
            enddo
          endif
c Fnu -> Hnu
          if (ytype.eq.'C') then
            do i=1,nfluxes
              stlhnu(i)=ifpi*stlhnu(i)
            enddo
          endif
c Eddington Hlam
          if (ytype.eq.'D') then
            do i=1,nfluxes
              lam=lmev/stlwev(i)
              lamcm=1.0d-8*lam
              freq=cls/lamcm
              stlhnu(i)=stlhnu(i)*lam/freq
            enddo
          endif
c Eddington Hnu
          if (ytype.eq.'E') then
            do i=1,nfluxes
              stlhnu(i)=stlhnu(i)
            enddo
          endif
c LLam, use r = 1.0d10
          if (ytype.eq.'F') then
            do i=1,nfluxes
              lam=lmev/stlwev(i)
              lamcm=1.0d-8*lam
              freq=cls/lamcm
              stlhnu(i)=ifpi*7.957747154594766e-22*stlhnu(i)*lam/freq
            enddo
          endif
c log10 LLam, use r = 1.0d10, like F
          if (ytype.eq.'L') then
            do i=1,nfluxes
              lam=lmev/stlwev(i)
              lamcm=1.0d-8*lam
              freq=cls/lamcm
              stlhnu(i)=10.d0**stlhnu(i)
              stlhnu(i)=ifpi*7.957747154594766e-22*stlhnu(i)*lam/freq
            enddo
          endif
c Lnu,  use r = 1.0d10
          if (ytype.eq.'G') then
            do i=1,nfluxes
              stlhnu(i)=ifpi*7.957747154594766e-22*stlhnu(i)
            enddo
          endif
c Xray Fnu: keV->ergs, keV->Hz
          if (ytype.eq.'H') then
            do i=1,nfluxes
              stlhnu(i)=ifpi*plk*stlhnu(i)
            enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 3, variable file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.7) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 7, WMBASIC A-Hnu two column atmosphere file,
c  variable lines,with a reversed of the energies
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data straight, and count fluxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nfluxes=0
          do i=1,mxflxbins
            read (luin,*,end=280) stlwave(i),stlhnu(i)
            nfluxes=nfluxes+1
          enddo
  280     continue
          close (luin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c convert and sort the data into increasing eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Angstroms
          if (stlwave(nfluxes).gt.stlwave(1)) then
c reverse if needed
            do i=1,nfluxes
              j=nfluxes-i+1
              if (i.lt.j) then
                swap=stlwave(j)
                stlwave(j)=stlwave(i)
                stlwave(i)=swap
                swap=stlhnu(j)
                stlhnu(j)=stlhnu(i)
                stlhnu(i)=swap
              endif
            enddo
          endif
          do i=1,nfluxes
            stlwev(i)=lmev/stlwave(i)
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 7, WMBASIC A-Hnu two column atmosphere file,
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
c
        if (coordstype.eq.9) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN Rauch HNi grid lambda flambda two
c  column atmosphere
c  Note curious ergs/cm^/s/cm/pi flambda units
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data straight, and count fluxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nfluxes=0
          do i=1,mxflxbins
            read (luin,*,end=290) stlwave(i),stlhnu(i)
            nfluxes=nfluxes+1
          enddo
  290     continue
          close (luin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c convert and sort the data into increasing eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Angstroms
          if (stlwave(nfluxes).gt.stlwave(1)) then
c reverse if needed
            do i=1,nfluxes
              j=nfluxes-i+1
              if (i.lt.j) then
                swap=stlwave(j)
                stlwave(j)=stlwave(i)
                stlwave(i)=swap
                swap=stlhnu(j)
                stlhnu(j)=stlhnu(i)
                stlhnu(i)=swap
              endif
            enddo
          endif
          do i=1,nfluxes
            stlwev(i)=lmev/stlwave(i)
          enddo
c
          do i=1,nfluxes
            lam=lmev/stlwev(i)
            lamcm=1.0d-8*lam
            freq=cls/lamcm
            stlhnu(i)=0.25*stlhnu(i)*lamcm/freq
          enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN HNi grid file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (coordstype.eq.5) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 5, MAPPINGS sou file with different PHOTDAT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write (*,*) 'Reading mismatched MAPPINGS .sou File: ',fname
          open (luin,file=fname,status='OLD')
c
c skip header lines and look for title, ID and nbins
c
          nheadlines=0
          do i=1,mxflxbins
            read (luin,'(a)') caract
            if (caract(1:1).ne.'%') goto 300
c            write (*,'(a80)') caract
            nheadlines=nheadlines+1
          enddo
  300     continue
c
c just use bin edges as given, without referring to an internal
c PHOTDAT initialised array since we are not trying to be binary
c exact.
c
          read (luin,*) souid2
          read (luin,*) nfluxes2
c          write(*,*) 'Found ID:',souID2,' and nfluxes:', nfluxes2
          do i=1,nfluxes2
            read (luin,*) stlwev(i),stlhnu(i)
          enddo
          close (luin)
c save total...
          blum0=0.d0
          do i=1,nfluxes2
            wid=(stlwev(i+1)-stlwev(i))
            blum0=blum0+wid*stlhnu(i)
          enddo
c
c use stlphnu(i) temporarily to use cell averages <a> to compute values
c at left edges, to get al, <a> and ar values which can conservatively
c convert to 3rd order cell centred values
c
c  ac = 1/4(6<a>-al-ar)  [ note: converse in PPM is <a>=1/6(4ac+al+ar) ]
c
c first add integral fluxes over two bins and get average left/right
c values in phnu
c al =  stlphnu(i), ar =  stlphnu(i+1)
c <a> = stlhnu(i)
c
          stlphnu(1)=stlhnu(1)
          stlphnu(nfluxes2)=stlhnu(nfluxes2)
          do i=2,nfluxes2-1
            wid0=stlwev(i)-stlwev(i-1)
            wid1=stlwev(i+1)-stlwev(i)
            stlphnu(i)=((stlhnu(i-1)*wid0)+(stlhnu(i)*wid1))/(stlwev(i+
     &       1)-stlwev(i-1))
          enddo
c
c use PPM transform from al, <a>, ar to al, ac, ar and
c set ac and mean wev to get 3rd order central points
c check for -ves
c
c  ac = 1/4(6<a>-al-ar)
c
          do i=1,nfluxes2-1
            ac=0.25d0*(6.d0*stlhnu(i)-stlphnu(i)-stlphnu(i+1))
            if ((stlhnu(i).gt.stlphnu(i)).and.(stlhnu(i).gt.stlphnu(i+1)
     &       )) then
              ac=stlhnu(i)
            endif
            ac=dmax1(ac,0.d0)
            stlhnu(i)=ac
          enddo
          stlhnu(nfluxes2)=stlhnu(nfluxes2-1)
c
c now set up central wevs
c
          wev0=stlwev(nfluxes2-1)
          do i=1,nfluxes2-1
            stlwev(i)=0.5d0*(stlwev(i)+stlwev(i+1))
          enddo
          stlwev(nfluxes2)=0.5d0*(wev0+stlwev(nfluxes2))
c
c          do i=1,nfluxes2
c            lghnu(i)=dlog10(stlhnu(i)+epsilon)
c          enddo
c
          nfluxes=nfluxes2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 5, MAPPINGS sou file with different PHOTDAT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  BEGIN rebinning
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Save the alphas and interval integrals; sflux
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,nfluxes-1
          stlsflux(i)=0.d0
          w0=stlwev(i)
          w1=stlwev(i+1)
          delta=w1-w0
          logdelta=dlog10(stlwev(i+1))-dlog10(stlwev(i))
          stlalpha(i)=dlog10(stlhnu(i+1)+epsilon)-dlog10(stlhnu(i)+
     &     epsilon)
          stlalpha(i)=stlalpha(i)/logdelta
          if (dabs(stlalpha(i)).lt.30.0d0) then
            a1=stlalpha(i)+1.d0
            scale=stlhnu(i+1)/(w1**stlalpha(i))
            sum=0.0d0
            if (a1.eq.0.0d0) then
              sum=dlog(w1/w0)
            else
              sum=((w1**(a1))-(w0**(a1)))/a1
            endif
            stlsflux(i)=(scale*sum)
          else
            stlsflux(i)=delta*0.5*(stlhnu(i+1)+stlhnu(i))
          endif
        enddo
        stlalpha(nfluxes)=0.d0
        stlsflux(nfluxes)=0.0d0
        rebinnedsrc(infph)=0.d0
        do i=1,infph-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          rebinnedsrc(i)=0.d0
          e0=photev(i)
          e1=photev(i+1)
          wid=e1-e0
          if ((e0.gt.stlwev(1)).and.(e1.lt.stlwev(nfluxes))) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c find source low and hi intervals that encompass the
c destination bin, each destination bin has edges of e0 - e1.
c  Each source data has points at energies wev.  Find source
c interval that contains e0 and source interval that contains
c e1, these may be the same or different source intervals.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            lowbin=0
            hibin=0
            do j=1,nfluxes
              if (e0.ge.stlwev(j)) then
                if (e0.lt.stlwev(j+1)) then
                  lowbin=j
                endif
              endif
              if (e1.ge.stlwev(j)) then
                if (e1.lt.stlwev(j+1)) then
                  hibin=j
                endif
              endif
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c get source integral between e0 and e1
c
c integrate from e0 to upper limit (w1) of the low interval,
c lobin integrate from the low limit (w0) to e1 in the high
c interval, hibin add the integrals (sflux) for any bins
c between the low and high interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ufrac=0.d0
            lfrac=0.d0
            sum=0.d0
            total=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c lfrac = integral of low interval from e0 to upper limit of lowbin,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(lowbin)
            w0=e0
            w1=stlwev(lowbin+1)
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=stlhnu(lowbin+1)/(w1**stlalpha(lowbin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              lfrac=(scale*sum)
            else
              mean=0.5*(stlhnu(lowbin)+stlhnu(lowbin+1))
              delta=(stlwev(lowbin+1)-stlwev(lowbin))
              f=(stlwev(lowbin+1)-e0)/delta
              lfrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c ufrac = integral of high interval from lower limit to e1,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(hibin)
            w0=stlwev(hibin)
            w1=e1
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=stlhnu(hibin+1)/(stlwev(hibin+1)**stlalpha(hibin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              ufrac=(scale*sum)
            else
              mean=0.5*(stlhnu(hibin)+stlhnu(hibin+1))
              delta=(stlwev(hibin+1)-stlwev(hibin))
              f=(e1-stlwev(hibin))/delta
              ufrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 1) if both e0 and e1 are in the same src interval, lowbin =
c hibin, then the integral between e0-e1 is
c lfrac+ufrac-stlsflux(lowbin), using saved interval integral in
c stlsflux(lowbin) from readin
c 2) if they are in adjacent bins, then the integral between
c e0 and e1 is: lfrac+ufrac
c 3) if in widely separated bins then the integral is
c lfrac+sum+ufrac where sum is the sum of stlsflux(i) bins
c between low and hibin
c BUT:  case 2 = 3 as sum = 0.0 in 2
c So only two cases  Same bin or different bins.
c ie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (lowbin.eq.hibin) then
              total=lfrac+ufrac-stlsflux(lowbin)
            else
              sum=0.d0
              if ((hibin-lowbin).gt.1) then
                do l=(lowbin+1),(hibin-1)
                  sum=sum+stlsflux(l)
                enddo
              endif
              total=lfrac+sum+ufrac
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Finally average integral over dest bin width, using same energy units
c as the integrals. Scale by 4 to get Inu for MAPPINGS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            rebinnedsrc(i)=4.0d0*(total/wid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (coordstype.eq.5) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 5, MAPPINGS sou file with different PHOTDAT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Check normalisation
c
          blum1=0.d0
          do i=1,infph-1
            wid=photev(i+1)-photev(i)
            blum1=blum1+wid*rebinnedsrc(i)
          enddo
c          write(*,*) 4.d0*blum0,blum1,4.d0*blum0/blum1
          blum0=4.d0*blum0/blum1
          do i=1,infph
            rebinnedsrc(i)=blum0*rebinnedsrc(i)
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 5, MAPPINGS sou file with different PHOTDAT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fname exists
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c
        write (*,*) ' PHOTSOU - READREBIN File Not Found: ',fname
        write (*,*) ' WARNING: No fluxes read.'
c
      endif
c
      return
      end
      subroutine readlinrebin (fname, coordstype, rebinnedsrc)
c
c Linear only interpolation for highly over sampled CMFGEN input spectra
c  8 =  Hz Hnu two column atmosphere file.
c  9 = CSPN Rauch TNMAP HNi/HCa grid lambda flambda
c      two column atmosphere file
c
      include 'cblocks.inc'
      include 'mpif.h'
c
c arguments
c
      character*256 fname
      integer*4 coordstype
      real*8 rebinnedsrc(mxinfph)
c
c  locals
c
      logical iexi
      integer*4 luin
c
      character*256 caract
c
      integer*4 i,j,l
      integer*4 nfluxes, nheadlines
      integer*4 lowbin, hibin
c
      real*8 ufrac, lfrac, mean, f
      real*8 delta, logdelta, swap
      real*8 a0, a1, freq, lam, lamcm
      real*8 e0, e1
      real*8 w0, w1
      real*8 scale, wid, sum, total
      real*8 wavescale
      real*8 freqscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Clear array first
c
      do i=1,mxinfph
        rebinnedsrc(i)=0.d0
      enddo
c
      if ((coordstype.ne.8).and.(coordstype.ne.9)) return
c
      luin=17
      inquire (file=fname,exist=iexi)
c
      if (iexi) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fname exists
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set some constants depending on coordstype
c
c  8 = Hz - Hnu two column atmosphere file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        wavescale=1.d0
        freqscale=1.d0
        nfluxes=0
c
        if (coordstype.eq.8) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 8, Hz Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        nheadlines=7
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          nfluxes=1
          do i=1,mxflxbins
            read (luin,*,end=10) stlwave(nfluxes),stlhnu(nfluxes)
            if (nfluxes.gt.1) then
              if (stlwave(nfluxes).ne.stlwave(nfluxes-1)) then
                nfluxes=nfluxes+1
              endif
            else
              nfluxes=nfluxes+1
            endif
          enddo
   10     continue
          nfluxes=nfluxes-1
c
c reverse, as CMFGEN has high freq first
c
          do i=1,nfluxes
            j=nfluxes-i+1
            if (i.lt.j) then
              swap=stlwave(j)
              stlwave(j)=stlwave(i)
              stlwave(i)=swap
              swap=stlhnu(j)
              stlhnu(j)=stlhnu(i)
              stlhnu(i)=swap
            endif
          enddo
c
          do i=1,nfluxes
            freq=stlwave(i)
            stlwev(i)=plk*freq/ev
          enddo
c
          close (luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 8, CMFGEN Hz-Hnu two column atmosphere file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (coordstype.eq.9) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN Rauch HNi grid lambda flambda two
c  column atmosphere
c  Note curious ergs/cm^/s/cm/pi flambda units
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nheadlines=37
          open (luin,file=fname,status='OLD')
          do i=1,nheadlines
            read (luin,'(a)') caract
            write (*,'(a80)') caract
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Scan the source data straight, and count fluxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nfluxes=0
          do i=1,mxflxbins
            read (luin,*,end=20) stlwave(i),stlhnu(i)
            nfluxes=nfluxes+1
          enddo
  20     continue
          close (luin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c convert and sort the data into increasing eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Angstroms
          if (stlwave(nfluxes).gt.stlwave(1)) then
c reverse if needed
            do i=1,nfluxes
              j=nfluxes-i+1
              if (i.lt.j) then
                swap=stlwave(j)
                stlwave(j)=stlwave(i)
                stlwave(i)=swap
                swap=stlhnu(j)
                stlhnu(j)=stlhnu(i)
                stlhnu(i)=swap
              endif
            enddo
          endif
          do i=1,nfluxes
            stlwev(i)=lmev/stlwave(i)
          enddo
c
          do i=1,nfluxes
            lam=lmev/stlwev(i)
            lamcm=1.0d-8*lam
            freq=cls/lamcm
            stlhnu(i)=0.25*stlhnu(i)*lamcm/freq
          enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  coordstype = 9, CSPN HNi grid file.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  BEGIN rebinning
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Save the alphas and interval integrals; sflux
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,nfluxes-1
          stlsflux(i)=0.d0
          w0=stlwev(i)
          w1=stlwev(i+1)
          delta=w1-w0
          logdelta=dlog10(stlwev(i+1))-dlog10(stlwev(i))
          stlalpha(i)=(dlog10(stlhnu(i+1))-dlog10(stlhnu(i)))/logdelta
          if (dabs(stlalpha(i)).lt.30.0d0) then
            a1=stlalpha(i)+1.d0
            scale=stlhnu(i+1)/(w1**stlalpha(i))
            sum=0.0d0
            if (a1.eq.0.0d0) then
              sum=dlog(w1/w0)
            else
              sum=((w1**(a1))-(w0**(a1)))/a1
            endif
            stlsflux(i)=(scale*sum)
          else
            stlsflux(i)=delta*0.5*(stlhnu(i+1)+stlhnu(i))
          endif
        enddo
        stlalpha(nfluxes)=0.d0
        stlsflux(nfluxes)=0.0d0
        rebinnedsrc(infph)=0.d0
        do i=1,infph-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          rebinnedsrc(i)=0.d0
          e0=photev(i)
          e1=photev(i+1)
          wid=e1-e0
          if ((e0.gt.stlwev(1)).and.(e1.lt.stlwev(nfluxes))) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c find source low and hi intervals that encompass the
c destination bin, each destination bin has edges of e0 - e1.
c  Each source data has points at energies wev.  Find source
c interval that contains e0 and source interval that contains
c e1, these may be the same or different source intervals.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            lowbin=0
            hibin=0
            do j=1,nfluxes
              if (e0.ge.stlwev(j)) then
                if (e0.lt.stlwev(j+1)) then
                  lowbin=j
                endif
              endif
              if (e1.ge.stlwev(j)) then
                if (e1.lt.stlwev(j+1)) then
                  hibin=j
                endif
              endif
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c get source integral between e0 and e1
c
c integrate from e0 to upper limit (w1) of the low interval,
c lobin integrate from the low limit (w0) to e1 in the high
c interval, hibin add the integrals (sflux) for any bins
c between the low and high interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ufrac=0.d0
            lfrac=0.d0
            sum=0.d0
            total=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c lfrac = integral of low interval from e0 to upper limit of lowbin,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(lowbin)
            w0=e0
            w1=stlwev(lowbin+1)
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=stlhnu(lowbin+1)/(w1**stlalpha(lowbin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              lfrac=(scale*sum)
            else
              mean=0.5*(stlhnu(lowbin)+stlhnu(lowbin+1))
              delta=(stlwev(lowbin+1)-stlwev(lowbin))
              f=(stlwev(lowbin+1)-e0)/delta
              lfrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c ufrac = integral of high interval from lower limit to e1,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(hibin)
            w0=stlwev(hibin)
            w1=e1
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=stlhnu(hibin+1)/(stlwev(hibin+1)**stlalpha(hibin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              ufrac=(scale*sum)
            else
              mean=0.5*(stlhnu(hibin)+stlhnu(hibin+1))
              delta=(stlwev(hibin+1)-stlwev(hibin))
              f=(e1-stlwev(hibin))/delta
              ufrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 1) if both e0 and e1 are in the same src interval, lowbin =
c hibin, then the integral between e0-e1 is
c lfrac+ufrac-stlsflux(lowbin), using saved interval integral in
c stlsflux(lowbin) from readin
c 2) if they are in adjacent bins, then the integral between
c e0 and e1 is: lfrac+ufrac
c 3) if in widely separated bins then the integral is
c lfrac+sum+ufrac where sum is the sum of stlsflux(i) bins
c between low and hibin
c BUT:  case 2 = 3 as sum = 0.0 in 2
c So only two cases  Same bin or different bins.
c ie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (lowbin.eq.hibin) then
              total=lfrac+ufrac-stlsflux(lowbin)
            else
              sum=0.d0
              if ((hibin-lowbin).gt.1) then
                do l=(lowbin+1),(hibin-1)
                  sum=sum+stlsflux(l)
                enddo
              endif
              total=lfrac+sum+ufrac
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Finally average integral over dest bin width, using same energy units
c as the integrals. Scale by 4 to get Inu for MAPPINGS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            rebinnedsrc(i)=4.0d0*(total/wid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fname exists
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c
        write (*,*) ' PHOTSOU - READLINREBIN File Not Found: ',fname
        write (*,*) ' WARNING: No fluxes read.'
c
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rebin (eev, hnu, nfluxes, rebinned)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Conservative powerlaw interpolation of highly non-linear data
c rebinning method.
c
c eev in eV, hnu in Hnu erg/cm^2/s/Hz/sr units, fluxes *at*
c eev points, output in <Hnu> over MV bins.
c
c Rebinned in current MV spectrum binning, average fluxes over
c bin with ephot at left edge of bin.
c
c Reuses stl common block temp vars to save stack space
c not intended for stellar libraries, readrebin is more
c efficient than making a local eev and hnu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Use power laws between src points, alpha, and interval
c analytic integrals, sflux, to fill destination vector with
c integrals of the src, over bins defined by bin edges in the
c destination = then average by destination bin width to get
c final, conservatively averaged, bin data.
c
c find piecewise powerlaw integrals between src points.
c from low -> high energy
c
c Use log-log coordinates to get alpha = slope,
c and linear coordinates to do analytical
c integrals between source points.
c
c Use analytical powerlaw integrals, or linear interpolation
c if magnitude of alpha is too great
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'mpif.h'
c
      integer*4 nfluxes
      real*8 eev(nfluxes),hnu(nfluxes)
      real*8 rebinned(mxinfph)
      integer*4 i,j,l
      integer*4 lowbin, hibin
      real*8 ufrac, lfrac, mean, f
      real*8 delta, logdelta, swap
      real*8 a0, a1
      real*8 e0, e1
      real*8 w0, w1, scale, wid, sum, total
c
      if (nfluxes.gt.mxflxbins) then
        write (*,*) ' Too many bins in rebin:',nfluxes,' > ',mxflxbins
      stop
      endif
c
c reverse if needed, low to high energy only...
c
      if (eev(nfluxes).lt.eev(1)) then
      do i=1,nfluxes
        j=nfluxes-i+1
        if (i.lt.j) then
          swap=eev(j)
          eev(j)=eev(i)
          eev(i)=swap
          swap=hnu(j)
          hnu(j)=hnu(i)
          hnu(i)=swap
        endif
      enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  BEGIN rebinning
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Save the alphas and interval integrals; sflux
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,nfluxes-1
          stlsflux(i)=0.d0
          w0=eev(i)
          w1=eev(i+1)
          delta=w1-w0
          logdelta=dlog10(eev(i+1))-dlog10(eev(i))
        stlalpha(i)=dlog10(hnu(i+1)+epsilon)-dlog10(hnu(i)+epsilon)
          stlalpha(i)=stlalpha(i)/logdelta
          if (dabs(stlalpha(i)).lt.30.0d0) then
            a1=stlalpha(i)+1.d0
            scale=hnu(i+1)/(w1**stlalpha(i))
            sum=0.0d0
            if (a1.eq.0.0d0) then
              sum=dlog(w1/w0)
            else
              sum=((w1**(a1))-(w0**(a1)))/a1
            endif
            stlsflux(i)=(scale*sum)
          else
            stlsflux(i)=delta*0.5*(hnu(i+1)+hnu(i))
          endif
        enddo
        stlalpha(nfluxes)=0.d0
        stlsflux(nfluxes)=0.0d0
        rebinned(infph)=0.d0
        do i=1,infph-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          rebinned(i)=0.d0
          e0=photev(i)
          e1=photev(i+1)
          wid=e1-e0
          if ((e0.gt.eev(1)).and.(e1.lt.eev(nfluxes))) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c find source low and hi intervals that encompass the
c destination bin, each destination bin has edges of e0 - e1.
c  Each source data has points at energies wev.  Find source
c interval that contains e0 and source interval that contains
c e1, these may be the same or different source intervals.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            lowbin=0
            hibin=0
            do j=1,nfluxes
              if (e0.ge.eev(j)) then
                if (e0.lt.eev(j+1)) then
                  lowbin=j
                endif
              endif
              if (e1.ge.eev(j)) then
                if (e1.lt.eev(j+1)) then
                  hibin=j
                endif
              endif
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c get source integral between e0 and e1
c
c integrate from e0 to upper limit (w1) of the low interval,
c lobin integrate from the low limit (w0) to e1 in the high
c interval, hibin add the integrals (sflux) for any bins
c between the low and high interval.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ufrac=0.d0
            lfrac=0.d0
            sum=0.d0
            total=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c lfrac = integral of low interval from e0 to upper limit of lowbin,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(lowbin)
            w0=e0
            w1=eev(lowbin+1)
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=hnu(lowbin+1)/(w1**stlalpha(lowbin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              lfrac=(scale*sum)
            else
              mean=0.5*(hnu(lowbin)+hnu(lowbin+1))
              delta=(eev(lowbin+1)-eev(lowbin))
              f=(eev(lowbin+1)-e0)/delta
              lfrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c ufrac = integral of high interval from lower limit to e1,
c using analytical powerlaw integrals, with saved alphas
c
c If alpha magnitude is too great, use linear integral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            a0=stlalpha(hibin)
            w0=eev(hibin)
            w1=e1
            if (dabs(a0).lt.30.0d0) then
              a1=a0+1.d0
              scale=hnu(hibin+1)/(eev(hibin+1)**stlalpha(hibin))
              sum=0.0d0
              if (a1.eq.0.0d0) then
                sum=dlog(w1/w0)
              else
                sum=((w1**(a1))-(w0**(a1)))/a1
              endif
              ufrac=(scale*sum)
            else
              mean=0.5*(hnu(hibin)+hnu(hibin+1))
              delta=(eev(hibin+1)-eev(hibin))
              f=(e1-eev(hibin))/delta
              ufrac=f*mean*delta
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 1) if both e0 and e1 are in the same src interval, lowbin =
c hibin, then the integral between e0-e1 is
c lfrac+ufrac-stlsflux(lowbin), using saved interval integral in
c stlsflux(lowbin) from readin
c 2) if they are in adjacent bins, then the integral between
c e0 and e1 is: lfrac+ufrac
c 3) if in widely separated bins then the integral is
c lfrac+sum+ufrac where sum is the sum of stlsflux(i) bins
c between low and hibin
c BUT:  case 2 = 3 as sum = 0.0 in 2
c So only two cases  Same bin or different bins.
c ie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (lowbin.eq.hibin) then
              total=lfrac+ufrac-stlsflux(lowbin)
            else
              sum=0.d0
              if ((hibin-lowbin).gt.1) then
                do l=(lowbin+1),(hibin-1)
                  sum=sum+stlsflux(l)
                enddo
              endif
              total=lfrac+sum+ufrac
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Finally average integral over dest bin width, using same energy units
c as the integrals. Keep as Hnu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            rebinned(i)=(total/wid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c src in range of destination
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop i over destination ephots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo
      return
      end
