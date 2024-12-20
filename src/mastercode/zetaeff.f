cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE THE EFFECTIVE DENSITY OF PHOTONS PER PARTICLE
c     IT TAKES INTO ACCOUNT THE DIFFERENT CROSS SECTIONS OF
c     PHOTOIONISATION AND THE RELATIVE ABUNDANCES OF ELEMENTS.
c     GIVES ALSO TO THE FIRST ORDER THE IONISING FRONT VELOCITY VIOFR
c     CALL SUBR. ALLRATES
c
      subroutine zetaeff (fi, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 abr,dh,fi
      real*8 ph
      real*8 wei
c
      integer*4 iel,ion
c
      character jjmod*4
c
      if (photonmode.ne.0) then
c
c
c    ***CALCULATE NEW PHOTOIONISATION/HEATING RATES IF NECESSARY.
c
        jjmod='PHOT'
        call allrates (100.d0, jjmod)
c
c
c    ***COMPUTES ZETAE IN COMMON BLOCK /DENLIN/
c
        wei=zion(1)
        ph=zion(1)*rphot(1,1)
c
        do 20 iel=2,atypes
          abr=zion(iel)
          wei=wei+abr
          ion=maxion(iel)-1
          ph=ph+(abr*rphot(ion,iel))
          do 10 ion=1,maxion(iel)-2
            wei=wei+abr
            ph=ph+(abr*(rphot(ion,iel)+auphot(ion,iel)))
   10     continue
   20   continue
c
        zetae=ph/(dh*wei)
c
c
c    ***COMPUTES PHOTON IMPACT PARAMETER ( USING QTOSOH COMPUTED IN
c     SUBROUTINE PHION ) WHICH CORRESPONDS ALSO TO VELOCITY OF
c     IONISING FRONT .
c
c      rion = 0.0d0
c      do 222 j = 1, atypes
c      rion = rion+(zion(j)*(dble(maxion(j)-1)**0.15d0))
c  222 continue
c
c      qhdh = qtosoh/(dh*rion)
c      viofr = qtosoh/((fi*dh)*rion)
c
        qhdh=qtosoh/dh
        qhdn=qhdh/zen
        uhdh=qhdh/cls
        viofr=qtosoh/((fi*dh)*zen)
c
      endif
      return
      end
