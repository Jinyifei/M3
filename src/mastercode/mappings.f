cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c Modelling And Prediction in PhotoIonised Nebulae & Gasdynamical Shocks
c M         A   P             P    I       N         G            S
c
c Developed 1975-1994 mainly at Mt. Stromlo Stromlo and Siding Spring
c Observatories, Institute of Advanced Studies, The Australian National
c University.
c
c copyright 1994 Ralph S. Sutherland (1), Michael A. Dopita (1)
c                      Luc Binette (2), and Brent Groves (1)
c
c
c (1) Research School of Astronomy & Astrophysics
c     Mount Stromlo Observatory
c     The Australian National University
c     Canberra, Australia
c
c     email:
c     Ralph.Sutherland@anu.edu.au
c
c This code may be freely copied for scientific research and educational
c purposes.  It may not be copied for commercial purposes without prior
c consent from the authours.
c
c No portion of this source code may be used elsewhere without
c permission and appropriate acknowledgements.
c
c Private modifications may be made, but do not distribute modified
c copies.  If you have modifications that should be distributed, then
c contact (ralph@mso.anu.edu.au), so that
c new versions can be made for general distribution.
c
c When MAPPINGS V is used, please refer to the code and version number
c in all publications.  The references for MAPPINGS V itself are:
c
c Dopita, M.A., 1976, Ap. J., 209, 395.
c Binette, L., Dopita, M.A. & Tuohy, I.R., 1985, Ap.J., 297, 476.
c Sutherland, R.S., & Dopita, M.A. 1993,ApJS, 88, 253.
c Groves
c Allen
c
c If you intend to use MAPPINGS V often then send email to
c ralph@mso.anu.edu.au so that a mailing list for update
c notices and bug reports can be compiled.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Version v5.1.13 of MAPPINGS V
c
c     RSS 11/00
c     BAG edition 04/04 added dust properties
c     DN 2012 Kappa distributions
c     RSS 2013/4 Complete Overhaul
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      program mappings
c
      include 'cblocks.inc'
c
c
      character ilgg*4, modeltype*12
      character abundtitle*24
      logical initerr
c
      theversion='v5.1.13'
c
   10 format(/,
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & /
     & '  Welcome to MAPPINGS V ',a8,/
     & /
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/)
      write (*,10) theversion
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     call the initialsation routine to read datafiles etc...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      initerr=.false.
      call mapinit (initerr)
      if (initerr) goto 110
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (photonmode.eq.0) then
        write (*,*) ' ***********************************************'
        write (*,*) ' *                                             *'
        write (*,*) ' *   WARNING: PHOTON FIELD DISABLED.           *'
        write (*,*) ' *                                             *'
        write (*,*) ' ***********************************************'
c
      endif
c
      if (alphacoolmode.eq.1) then
        write (*,*) ' ***********************************************'
        write (*,*) ' *                                             *'
        write (*,*) ' *   WARNING: POWERLAW COOLING ENABLED.        *'
        write (*,*) ' *                                             *'
        write (*,*) ' ***********************************************'
c
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     begin interactive initialisation and then subroutine selector...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Elemental Abundances:'/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
   30 write (*,20)
c
c    ***CHANGE ELEMENTAL ABUNDANCES ?
c      FINDS ZGAS RELATIVE TO THE SUN
c
      call abecha
c
c read in an optional abundance offsets file
c
      call deltaabund
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     *** Setup electron distributions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call kappainit
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     *** Setup dust grain properties
c
c     Now available (averages over all grainsizes)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call grainpar
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***SELECT A PROGRAM
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   40  format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V: End of General Initialisation'/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      write (*,40)
      abundtitle=' Final Gas Abundances :'
      call dispabundances (6, zion, abundtitle)
      modeltype=' '
      if (usekappa) then
        if (grainmode.eq.1) then
          modeltype='Dust, Kappa'
        else
          modeltype='Kappa On'
        endif
      else
        if (grainmode.eq.1) then
          modeltype='Dust Enabled'
        else
          modeltype=' '
        endif
      endif
   50 if (usekappa) then
        write (*,70) modeltype,kappa
        write (*,90)
      else
        write (*,60) modeltype
        write (*,80)
      endif
   60 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V: Choose a Model: ',a,/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/)
   70 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V: Choose a Model: ',a,'(k=',f6.1,')',/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/)
   80 format(
     & '    MM  :  Multi-Level Ion Emissivity'/
     & '    CD  :  Multi-Level Ion Critical Densities'/
     & '    SS  :  Single slab models'/
     & '    CC  :  CIE Cooling curves'/
     & '    NC  :  NEQ Cooling curves'/
     & '    PP  :  PIE Photoionization curves'/
     & '        :'/
     & '    S2  :  Shock + diffuse field ( v < 150 kms/s )'/
     & '    S3  :  Shock, elect. + ion Temps ( v < 150 kms/s )'/
C     & '    S4  :  Shock, adaptive mesh, requires ext. preionisation'/
     & '    S5  :  Shock, adaptive mesh, requires ext. preionisation'/
     & '        :'/
     & '    P4  :  Photo. Abs. distance step, full cont.'/
     & '    P5  :  P4 + Radiation Pressure & Dust'/
     & '    P6  :  P5 + Robust Integrator'/
     & '    P7  :  P6 + Experimental'/
     & '        :'/
     & '    MP  :  P7 + Monte Carlo, (3D Geometry)'/
     & '    MA  :  Re-compute Monte Carlo'/
     & '        :'/
     & '    R   :  Reinitialise'/
     & '  E,X,Q :  Exit'/
     &       ' :: ',$)
   90 format(
     & '    MM  :  Multi-Level Ion Emissivity'/
     & '    SS  :  Single slab models'/
     & '    CC  :  CIE Cooling curves'/
     & '    NC  :  NEQ Cooling curves'/
     & '    PP  :  PIE Photoionization curves'/
     & '        :'/
     & '    P4  :  Photo. Abs. distance step, full cont.'/
     & '    P5  :  P4 + Radiation Pressure & Dust'/
     & '    P6  :  P5 + Robust Integrator'/
C     & '    P7  :  P6 + Experimental'/
     & '        :'/
     & '    R   :  Reinitialise'/
     & '  E,X,Q :  Exit'/
     &       ' :: ',$)
c
      read (*,100) ilgg
  100 format(a)
      ilgg=ilgg(1:2)
c
c clear up any gobbled line feeds...
c
      write (*,*)
c
c     fiddle lower case entries..
c
      if (ilgg(1:2).eq.'rm') ilgg='RM'
      if (ilgg(1:2).eq.'mm') ilgg='MM'
      if (ilgg(1:2).eq.'cd') ilgg='CD'
      if (ilgg(1:2).eq.'ss') ilgg='SS'
      if (ilgg(1:2).eq.'pp') ilgg='PP'
      if (ilgg(1:2).eq.'cc') ilgg='CC'
      if (ilgg(1:2).eq.'nc') ilgg='NC'
      if (ilgg(1:2).eq.'mp') ilgg='MP'
      if (ilgg(1:2).eq.'ma') ilgg='MA'
      if (ilgg(1:1).eq.'s') ilgg(1:1)='S'
      if (ilgg(1:1).eq.'p') ilgg(1:1)='P'
      if (ilgg(1:1).eq.'o') ilgg(1:1)='O'
c      if (ilgg(1:1) .eq.'t') ilgg(1:1) ='T'
      if (ilgg(1:1).eq.'r') ilgg(1:1)='R'
      if (ilgg(1:1).eq.'e') ilgg(1:1)='E'
c
c     X for exit as well
c
      if (ilgg(1:1).eq.'X') ilgg(1:1)='E'
      if (ilgg(1:1).eq.'x') ilgg(1:1)='E'
c
c     Q for exit as well
c
      if (ilgg(1:1).eq.'Q') ilgg(1:1)='E'
      if (ilgg(1:1).eq.'q') ilgg(1:1)='E'
c
c    ***ZERO BUFFER ARRAYS AND RESET COUNTERS
c
      call zer
c
      runname='Test Run'
c
      flgrc = .false.
c
      if (ilgg(1:2).eq.'S2') call shock2
      if (ilgg(1:2).eq.'S3') call shock3
      if (ilgg(1:2).eq.'S4') call shock4
      if (ilgg(1:2).eq.'S5') call shock5
      if (ilgg(1:2).eq.'P4') call photo4
      if (ilgg(1:2).eq.'P5') call photo5
      if (ilgg(1:2).eq.'P6') call photo6
      if (ilgg(1:2).eq.'P7') call photo7
      if (ilgg(1:1).eq.'E') goto 130
      if (ilgg(1:1).eq.'R') goto 30
      if (ilgg(1:2).eq.'CC') call coolc
      if (ilgg(1:2).eq.'NC') call neqc
      if (ilgg(1:2).eq.'MM') call ionemit
      if (ilgg(1:2).eq.'CD') call critdens
      if (ilgg(1:2).eq.'SS') call slab
      if (ilgg(1:2).eq.'PP') call phocrv
      if (ilgg(1:2).eq.'MP') call montph7
      if (ilgg(1:2).eq.'MA') then
         flgrc = .true.
         call montph7
         flgrc = .false.
      endif
c
c
      goto 50
c
c error message
c
  110 continue
      write (*,120)
  120 format(/,' ERROR: Failed to Initialise. Immediate Exit.',/)
  130 continue
c
  140 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V: Session Ended. ',/
     & ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',//)
      write (*,140)
      stop
c
      end
