cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE ALL ATOMIC RATES (RECOMBINATION,CHARGE EX.,
c                                    COLLIS. IONIS.,PHOTOIONISATION)
c     USES FLAGS TEM AND IPHOM TO DETERMINE IF PREVIOUS RATES
c     ARE STILL APPLICABLE.
c     (THE PHOTON FIELD IS CONSIDERED CHANGED WHENEVER SUBR. TOTPHOT
c      IS CALLED (IPHO IS RESET))
c
c     JJMOD = 'TEMP' TO CALCULATE ONLY TEMPERATURE DEPENDANT RATES
c     JJMOD = 'PHOT' TO CALCULATE ONLY PHOTOIONISATION/HEATING RATES
c     JJMOD = 'ALL'  TO CALCULATE BOTH TYPES
c
c     CALL SUBR. PHION,COLLION,CHAREX,RECOM
c
c
      subroutine allrates (t, jjmod)
c
      include 'cblocks.inc'
c
      real*8 t, telc
c
      character jjmod*4
c
      telc=dmax1(t,1.d-2)
      if (((jjmod.ne.'ALL').and.(jjmod.ne.'PHOT')).and.(jjmod.ne.'TEMP')
     &) then
        write (*,10) jjmod
   10 format('  MODE WRONGLY DEFINED FOR SUBR. ALLRATES :',a4)
        stop
      endif
c
c
c     always call cosmic first to get sec ion rate, heating unused here
c     so t de dh not relavent, so long as pop is up to date.
c
      call cosmic (0.d0, 0.d0, 0.d0)
c
      if ((iphom.ne.ipho).and.(jjmod.ne.'TEMP')) then
        call phion
        iphom=ipho
      endif
c
      if ((telc.ne.tem).and.(jjmod.ne.'PHOT')) then
        call recom (telc)
        call collion (telc)
        call charex (telc)
        tem=telc
      endif
c
      return
      end
