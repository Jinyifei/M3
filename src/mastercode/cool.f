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
c*******COMPUTES TOTAL COOLING RATE OF PLASMA , TEMP : T
c     ELECTRON DENSITY : DE , HYDROGEN DENSITY : DH
c     RETURNS TLOSS : TOTAL COOLING RATE (ERG.CM-3.S-1)
c     CALL SUBR.  HYDRO,RESON,INTER,FORBID,FREFRE,COLOSS,
c                 PHEAT,NETGAIN,ALLRATES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine cool (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 tll, ell, egg, ett
      character jjmod*4
      real*8 t, de, dh
      real*8 q, v, alpha_cool, g, cooltime,csum
      real*8 frectim2, fpressu, frho, fcietim
      integer*4 i, j
c
c clear element/ion cooling totals
c
      do j=1,atypes
        coolz(j)=0.d0
        heatz(j)=0.d0
        do i=1,maxion(j)
          coolzion(i,j)=0.d0
          heatzion(i,j)=0.d0
        enddo
      enddo
      oiii5007loss=0.d0
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod='ALL'
      call allrates (t, jjmod)
c
c    ***New hydrogen + helium II line cooling
c    ***Also calls 2 photon and He I calcs
c
      hheloss=0.d0
      hloss=0.d0
      if (hhecollmode.eq.0) call hhecoll (t, de, dh)
      call hydro (t, de, dh)
c
      tll=hloss
c
c Heavy element recombination, negligible cooling
c
      call recom_cii (t, de, dh)
      call recom_nii (t, de, dh)
      call recom_oi (t, de, dh)
      call recom_oii (t, de, dh)
      call recom_neii (t, de, dh)
c
c
c    LEGACY RESONANCE,INTERCOMBINATION AND FORBIDDEN LINES
c
      call fine3 (t, de, dh)
      call inter (t, de, dh)
c      call inter2 (t, de, dh)
c
c multi-levels must be called before the reson routines
c
      call multilevel (t, de, dh)
      call multiiron (t, de, dh)
c
c      call reson (t, de, dh)
c      call reson2 (t, de, dh)
      rloss=0.0d0
      call reson3 (t, de, dh)
      call resonl (t, de, dh)
c
c      call helif (t, de, dh)
c

      tll=tll+rloss+fslos+fmloss+xr3loss+xrlloss
      tll=tll+f3loss+feloss
c
c
c    ***FREE-FREE COOLING
c
      call frefre (t, de, dh)
      tll=tll+fflos
c
      cmplos=0.d0
c
c    ***COLLISIONAL IONISATION LOSSES
c
      call coloss (de, dh) ! colos
      tll=tll+colos
c
C      tll=hloss+rloss+fslos+fmloss+xrloss+xr3loss+xrlloss+xiloss
C      tll=tll+fflos+colos
C      tll=tll+f3loss+feloss+gcool
c
c    ***RECOMB. COOLING (PLUS ESTIMATE OF ON THE SPOT HEATING IF APPL)
c
      call netgain (t, de, dh)
c      tll=tll-rngain
c
c    ***PHOTOIONISATION HEATING
c
      call pheat (de, dh) ! pgain
c      tll=tll-pgain
c
c    ***PHOTOELECTRIC GRAIN HEATING AND COLLISIONAL GRAIN LOSSES
c
      call hgrains (t, de, dh)
      call hpahs (t, de, dh)
c
c      tll=tll-gheat+gcool-paheat
c
c    ***CHARGE EXCHANGE HEATING
c
      call cheat (dh) ! chgain
c      tll=tll-chgain
c
c    ***COSMIC RAY HEATING
c
      call cosmic (t, de, dh)
c      tll=tll-cosgain
c
c    ***MICROTURBULENT DISSPATION HEATING
c
      q=0.0d0
      if (turbheatmode.gt.0) then
        g=1.66666666666667d0
        dens=frho(de,dh)
        pres=fpressu(t,dh,pop)
        v=sqrt(g*pres/dens)*admach
        alpha_cool=0.d0
        if (turbheatmode.eq.1) then
          alpha_cool=1.d0/frectim2(de)
        endif
        if (turbheatmode.eq.2) then
          alpha_cool=1.d0/fcietim(t,de,dh)
        endif
        if (turbheatmode.eq.3) then
          cooltime=1.5d0*pres/tll
          alpha_cool=1.d0/cooltime
        endif
        if (turbheatmode.eq.4) then
          alpha_cool=1.d0/alphaturbfixed
        endif
        q=0.5d0*dens*v*v*alpha_cool
        tll=tll-q
      endif
c
c    ***(SEPARATELY) EFFECTIVE LOSS AND GAIN  :  ELOSS,EGAIN
c
      ell=hloss+rloss+fslos+fmloss+xr3loss+xrlloss
      ell=ell+fflos+colos
      ell=ell+f3loss+feloss+gcool
c
      egg=pgain+cosgain+paheat+gheat+q
c
      if (rngain.lt.0.0d0) then
        ell=ell-rngain
      else
        egg=egg+rngain
      endif
c
      if (chgain.lt.0.0d0) then
        ell=ell-chgain
      else
        egg=egg+chgain
      endif
      ett=ell+egg
      if (ett.gt.0.d0) then
        dlos=(ell-egg)/ett
      else
        dlos=1.0d0
      endif
c
      csum=0.d0
      do j=1,atypes
        do i=1,maxion(j)
          if (dabs(coolzion(i,j)).lt.epsilon) coolzion(i,j)=0.d0
        enddo
        if (dabs(coolz(j)).lt.epsilon) coolz(j)=0.d0
        csum=csum+coolz(j)
      enddo
c
      tloss=tll
      egain=egg
      eloss=ell
c
      if (alphacoolmode.eq.1) then
        tloss=de*de*alphac0*((1.0d-6*t)**alphaclaw)
        egain=0.d0
        eloss=tloss
      endif
c
      if (dabs(tloss).lt.epsilon) tloss=0.d0
      if (dabs(egain).lt.epsilon) egain=0.d0
      if (dabs(eloss).lt.epsilon) eloss=0.d0
c
c       write(*,*) 'Cool:',tloss, egain, eloss, csum
c
      return
      end
