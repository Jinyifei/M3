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
c
c     NEQC; non-equilibrium cooling, simplified code based on Shock5
c     Calls compsh5
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine neqc ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i
c
c     set up logical unit numbers
c
      lual=20
      luop=21
      lurt=22
      ludy=23
      lusp=24
      lupf=25
      lupb=26
      lucl=27
c
      ieln=4
      do i=1,ieln
        luions(i)=28+i
      enddo
c
      fl=' '
      pfx='shckn'
      sfx='sh5'
      call newfile (pfx, 4, sfx, 3, fn)
      fl=fn(1:12)
c
      open (luop,file=fl,status='NEW')
c
c     zero arrays and define modes: all ions, no on-thespot,
c     full continuum
c     calculations..
c
      call zer
c
      jspot='NO'
      jcon='YES'
      jspec='NO'
      ispo='SII'
      allmod='NO'
      fclmod='NO'
      tsrmod='NO'
      ratmod='NO'
      dynmod='NO'
      vmod='NONE'
      mypfx='neqcl'
      specmode='DW'
      jgeo='P'
      jtrans='LODW'
      jden='B'
      subname='NEQ Cool'
      jnorm=0
c
      photonmode=1
c
      call neqcsetup ()
c
      call neqcheaders ()
c
      call compsh5 (0,1)
c
      close (luop)
c
      photonmode=1
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine neqcsetup ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 teinit,dhinit
      real*8 ptotal
      real*8 feldens,fpressu,frho
c
   10 format(a)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTRO - Define Input Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   20 format(///,
     & ' Non-equilibrium Cooling model selected:',/,
     & ' Based on v5.1 S5 model',/,
     & ' ****************************************************')
      write (*,20)
c
c     get ionisation balance to work with...
c
      subname='NEQ Cool'
c
   30 format(///,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Setting the Initial Conditions',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,30)
c
        teinit=1.0d9
        dhinit=1.0d-4
        magparam=1.0d0
        machnumber=1.0d0
   90   format(//,' Initial state:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '  Temperature, H density, Mach number, Magnetic Alpha',/,
     &  '  (T<10 taken as a log10, dh <0.0 as log10, )',/,
     &  '  (Mach<=1.0, Alpha>=0.0 Pmag/Pgas, < 0.0 as -microgauss)',/,
     &  ' : T (K), dh (N), M, Alpha',/,
     &  ' : ')
        write (*,90)
        read (*,*) teinit,dhinit,machnumber,magparam
c
        if (teinit.le.10.d0) teinit=10.d0**teinit
        if (dhinit.lt.0.d0) dhinit=10.d0**dhinit
c
c  overwrites t,te0,te1 by calling zer, needs fixing!
c  t and te0 are saved atm
c  need to change globals to fix properly
c
        call ciepops(teinit,dhinit)
c
        dh=dhinit
        te0=teinit
        te1=te0
c
        if (machnumber.gt.1.d0) machnumber=1.0d0
        if (machnumber.lt.0.d0) machnumber=0.0d0
c
        wdil=0.5d0
c
c     Magnetic fields expressed as Alpha = Pmag/Pgas
c     Alpha = 1.0 is equipartition, Alpha = 10.0 is strong magnetic
c     Alpha = 0.1 is weak magnetic field. Alpha = 0.0 is no magnetic
c     Alpha = 1.0/Beta where Beta is the usual magnetic parameter, but
c     Alpha allows Alpha = 0.0 for the non-magnetic limit
c
        pgas=fpressu(te0,dh,pop)
c
        if (magparam.lt.0.0d0) then
          bmag=-magparam
        endif
        if (dr.le.0.d0) dr=1.d0
        if (wdil.gt.0.5d0) wdil=0.5d0
        if (magparam.ge.0.d0) then
          pmag=magparam*pgas
          bmag=dsqrt(epi*pmag)*1.0d6
        else
          bmag=-magparam
        endif
        bmag=bmag*1.d-6
c
        bm0=bmag
c
c Now set constant values at steady velocity, no shock.
c
        dhpo=dh
        dh0=dh
        dh1=dh
c
        de=feldens(dh,pop)
        de0=feldens(dh,pop)
        de1=feldens(dh,pop)
c
        rho0=frho(de0,dh0)
        rho1=frho(de1,dh1)
        pr0=fpressu(te0,dh0,pop)
        pr1=fpressu(te1,dh1,pop)
c
        vpo=dsqrt(gammaEOS*pr1/rho1)*machnumber
        vel0=vpo
        vel1=vpo
c
        bm1=bm0*rho1/rho0
        tm00=te0
c
        bp0=(bm0*bm0)/epi
        bp1=(bm1*bm1)/epi
c
        machnumber=vel0/dsqrt(gammaeos*pr0/rho0)
        alfvennumber=vel0/dsqrt(2.d0*bp0/rho0)
        malpha=bp0/pr0
        pram=rho1*vel1*vel1
        gaseta=gammaeos*pr1/pram
        mageta=2.d0*pr1/pram
        ptotal=(pr1+bp1+pram)
c
  100   format(//' Shock Properties',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '   Velocity                  :',1pg12.5,' km/s',/,
     &  '   Flow Mach   Number        :',1pg12.5,/,
     &  '   Flow Alfven Number        :',1pg12.5,/,
     &  '   Flow Alpha (Pmag/Pgas)    :',1pg12.5,/,
     &  '   Flow Gas Eta (gPgas/Pram) :',1pg12.5,/,
     &  '   Flow Mag Eta (2Pmag/Pram) :',1pg12.5,/,
     &  '   Flow  v:',1pg12.5,' km/s',/,
     &  '   Flow ne:',1pg12.5,' cm^-3',/,
     &  '   Flow nH:',1pg12.5,' cm^-3',/,
     &  '   Flow d :',1pg12.5,' g/cm^-3',/,
     &  '   Flow P :',1pg12.5,' dyne/cm^2',/,
     &  '   Flow T :',1pg12.5,' K',/,
     &  '   Flow B :',1pg12.5,' microGauss',/,
     &  '   Flow BP:',1pg12.5,' dyne/cm^2',/,
     &  '   Flow Ram Pressure:',1pg12.5,' dyne/cm^2',/,
     &  '   Ptotal     :',1pg12.5,' dyne/cm^2',/,
     &  '   Pgas/Ptotal:',1pg12.5,' (Eta)',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,100) vel0*1d-5,machnumber,alfvennumber,
     & malpha,gaseta,mageta,
     & vel1*1d-5,de1,dh1,rho1,pr1,te1,bm1*1.d6,bp1,pram,
     & ptotal,pr1/ptotal
c
      if (alphacoolmode.eq.1) then
  110  format(//,
     &  ' Powerlaw Cooling enabled  Lambda T6 ^ alpha :',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Give index alpha, Lambda at 1e6K',/,
     &  '    (Lambda<0 as log)',/,
     &  ' :: ')
        write (*,110)
        read (*,*) alphaclaw,alphac0
        if (alphac0.lt.0.d0) alphac0=10.d0**alphac0
      endif
c
c     Diffuse field interaction
c
c
  120 format(//,
     & ' Choose Diffuse Field Option :',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    F  : Full diffuse field interaction (Default).',/,
     & '    Z  : Zero diffuse field interaction.',/,
     & ' :: ')
      write (*,120)
  130 read (*,10) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'z') ilgg='Z'
      if (ilgg.eq.'f') ilgg='F'
c
      if ((ilgg.ne.'Z').and.(ilgg.ne.'F')) goto 130
c
      if (ilgg.eq.'Z') photonmode=0
      if (ilgg.eq.'F') photonmode=1
c
c  140 format(//,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     & '  Calculation Settings ',/,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c      write (*,140)
c
c     Set up time step limits
c
c 140 format(//,' Choose a Time Step Behaviour :',/,
c    & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c    & '    A  : Auto. time steps, based on Atomic timescales.',/,
c    & '    M  : Choose a Maximum time step, Auto otherwise.',/,
c    & '    N  : Choose a miNimum time step, Auto otherwise.',/,
c    & '    S  : Strict preset timestep. (not recomended)',//,
c    & ' :: ')
cc
c 150 write (*,140)
c     read (*,10) tmod
c     tmod=tmod(1:1)
c     if (tmod.eq.'a') tmod='A'
c     if (tmod.eq.'m') tmod='M'
c     if (tmod.eq.'n') tmod='N'
c     if (tmod.eq.'s') tmod='S'
cc
c     if ((tmod.ne.'A').and.(tmod.ne.'M').and.(tmod.ne.'N')
c    &.and.(tmod.ne.'S')) goto 150
cc
c
c no long a useful option, set strict values here
c perhas enable for expert mode later
c
      tmod='A'
      atimefrac=0.1d0
      cltimefrac=0.05d0
      abtimefrac=1.0d0
c
c     if (tmod.ne.'S') then
c 160    format(//,
c    &   ' Give timescale fractions:',/,
c    &   ' ( 0 < dtau < 1 )',/,
c    &   ' ( recommend: 0.05 < dtau < 0.25 )',/,
c    &   ' : atomic, cooling, absorbsion',/,
c    &   ' : ')
c 170   write (*,160)
c       read (*,*) atimefrac,cltimefrac,abtimefrac
c
c       if (atimefrac.le.0.d0) goto 170
c       if (cltimefrac.le.0.d0) goto 170
c       if (abtimefrac.le.0.d0) goto 170
c
c       if (atimefrac.ge.1.d0) goto 170
c       if (cltimefrac.ge.1.d0) goto 170
c       if (abtimefrac.ge.1.d0) goto 170
c
c     endif
c
      utime=0.d0
c     if (tmod.ne.'A') then
c 180   format(//,
c    &  ' Give the time step:',/,
c    &  ' ( time < 100 taken as a log)',//,
c    &  ' : time (sec)',/,
c    &  ' : ')
c       write (*,180)
c       read (*,*) utime
c       if (utime.le.100.d0) utime=10.d0**utime
c     endif
c
c
  150 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Boundry Conditions ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,150)
c
c     Choose ending
c
  160 format(/,
     & '  Choose Ending :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard ending, 1% weighted ionisation.',/,
     & '    B  : Species Ionisation limit.',/,
     & '    C  : Temperature limit.',/,
     & '    D  : Distance limit.',/,
     & '    E  : Time limit.',/,
     & '    F  : Thermal Balance limit.',/,
     & '    G  : Heating limit.',/,
     & ' :: ')
c
  170 write (*,160)
      read (*,10) jend
      jend=jend(1:1)
c
      if (jend.eq.'a') jend='A'
      if (jend.eq.'b') jend='B'
      if (jend.eq.'c') jend='C'
      if (jend.eq.'d') jend='D'
      if (jend.eq.'e') jend='E'
      if (jend.eq.'f') jend='F'
      if (jend.eq.'g') jend='G'
c
      if ((jend.ne.'A').and.(jend.ne.'B').and.(jend.ne.'C')
     &.and.(jend.ne.'D').and.(jend.ne.'E').and.(jend.ne.'F')
     &.and.(jend.ne.'G')) goto 170
c
c     Secondary info:
c
      if (jend.eq.'B') then
  180   format(//,' Give Atom, Ion and limit fraction: ')
        write (*,180)
        read (*,*) ielen,jpoen,fren
      endif
      if (jend.eq.'C') then
  190   format(//,' Give final temperature (K > 10, log <= 10): ')
        write (*,190)
        read (*,*) tend
        if (tend.le.10.d0) tend=10.d0**tend
      endif
      if (jend.eq.'D') then
  200   format(//,' Give final distance (cm > 100, log<=100): ')
        write (*,200)
        read (*,*) diend
        if (diend.le.100.d0) diend=10.d0**diend
      endif
      if (jend.eq.'E') then
  210   format(//,' Give time limit (s > 100, log<=100): ')
        write (*,210)
        read (*,*) timend
        if (timend.le.100.d0) timend=10.d0**timend
      endif
c
  220 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Output Requirements ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,220)
c
c     default monitor elements
c
      iel1=zmap(1)
      iel2=zmap(2)
      iel3=zmap(6)
      iel4=zmap(8)
      mypfx='neqcl'
c
c  230 format(//,' Choose output settings : ',/,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     & '    A  : Standard output. ',/,
c     & '    B  : Standard plus ion balance files.',/,
c     & '    C  : Standard plus dynamics file.',/,
c     & '    D  : Standard plus rates file.',/,
c     & '    E  : Standard plus final down stream field.',/,
c     & '    F  : Standard plus upstream field at each step.',/,
c     & '    J  : Cooling Components by Elements file.',/,
c     & '       :',/,
c     & '    G  : B + E',/,
c     & '    H  : B + F',/,
c     & '       :',/,
c     & '    I  : Everything',//,
c     & ' :: ')
c  240 write (*,230)
c      read (*,10) ilgg
c      ilgg=ilgg(1:1)
c
      ilgg='J'
c
c     if (ilgg.eq.'a') ilgg='A'
c     if (ilgg.eq.'b') ilgg='B'
c     if (ilgg.eq.'c') ilgg='C'
c     if (ilgg.eq.'d') ilgg='D'
c     if (ilgg.eq.'e') ilgg='E'
c     if (ilgg.eq.'f') ilgg='F'
c     if (ilgg.eq.'g') ilgg='G'
c     if (ilgg.eq.'h') ilgg='H'
c     if (ilgg.eq.'i') ilgg='I'
c     if (ilgg.eq.'j') ilgg='J'
c
c      if ((ilgg.ne.'A').and.(ilgg.ne.'B').and.(ilgg.ne.'C')
c     &.and.(ilgg.ne.'D').and.(ilgg.ne.'E').and.(ilgg.ne.'F')
c     &.and.(ilgg.ne.'G').and.(ilgg.ne.'H').and.(ilgg.ne.'I')
c     &.and.(ilgg.ne.'J')) goto 240
c
c     Upstream fields
c
      lmod='Y'
c
c      if ((ilgg.eq.'F').or.(ilgg.eq.'H').or.(ilgg.eq.'I')) then
c        lmod='Y'
c      endif
c
c     F-Lambda plots, ionising and optical
c
      jspec='N'
c      if ((ilgg.eq.'E').or.(ilgg.eq.'G').or.(ilgg.eq.'I')) then
c        jspec='YES'
c      endif
c
c     timescales and rates file
c
      ratmod='N'
      if ((ilgg.eq.'D').or.(ilgg.eq.'I')) then
        ratmod='Y'
        pfx='rates'
        sfx='csv'
        call newfile (pfx, 5, sfx, 3, fn)
        fr=fn(1:13)
        open (lurt,file=fr,status='NEW')
      endif
c
c     Cooling Components File
c
      fclmod='N'
      if ((ilgg.eq.'J').or.(ilgg.eq.'I')) then
c      lmod='Y'
c      jspec='YES'
        jnorm=0
  250  format (' Cooling File Normalisation,',/
     &         ' (0=ne.nH, 1=nH^2, 2=ne.ni, 3=n^2, 4=ne^2):')
        write (*,250)
        read (*,*) jnorm
        if (jnorm.lt.0) jnorm=0
        if (jnorm.gt.4) jnorm=0
        fclmod='Y'
        pfx='cc'
        sfx='csv'
        fcl=' '
        call newfile (pfx, 2, sfx, 3, fn)
        fcl=fn(1:10)
        open (lucl,file=fcl,status='NEW')
      endif
c
c     flow dynamics file
c
      dynmod='N'
      if ((ilgg.eq.'C').or.(ilgg.eq.'I')) then
        dynmod='Y'
        pfx='dyn'
        sfx='csv'
        fd=' '
        call newfile (pfx, 3, sfx, 3, fn)
        fd=fn(1:11)
        open (ludy,file=fd,status='NEW')
      endif
c
      pfx='spec'
      sfx='csv'
      fsp=' '
      call newfile (pfx, 4, sfx, 3, fn)
      fsp=fn(1:12)
      open (lusp,file=fsp,status='NEW')
c
c     ionbalance files
c
      allmod='n'
      tsrmod='n'
c
      if ((ilgg.eq.'B').or.
     &    (ilgg.eq.'G').or.
     &    (ilgg.eq.'J').or.
     &    (ilgg.eq.'H').or.
     &    (ilgg.eq.'I')) then
c
        if ((ilgg.eq.'B').or.(ilgg.eq.'G').or.(ilgg.eq.'H')) then
c
  260    format(//,' Record all ions file (Y/N)?: ')
          write (*,260)
          read (*,10) allmod
          allmod=allmod(1:1)
          if (allmod.eq.'y') allmod='Y'
        endif
c
        if (ilgg.eq.'I') allmod='Y'
c
        if (allmod.eq.'Y') then
          pfx='allion'
          sfx='csv'
          call newfile (pfx, 6, sfx, 3, fn)
          fa=fn(1:14)
          open (lual,file=fa,status='NEW')
        endif
c
        if (ilgg.eq.'B') then
  270    format(//,' Record particular ion balances (Y/N)?: ')
          write (*,270)
          read (*,10) tsrmod
          tsrmod=tsrmod(1:1)
          if (tsrmod.eq.'y') tsrmod='Y'
        endif
c
        if (ilgg.eq.'I') tsrmod='Y'
        if (ilgg.eq.'J') tsrmod='Y'
c
        if (tsrmod.eq.'Y') then
  280       format(//,' Give 4 elements by atomic number: ')
          write (*,280)
          read (*,*) iel1,iel2,iel3,iel4
          iel1=zmap(iel1)
          iel2=zmap(iel2)
          iel3=zmap(iel3)
          iel4=zmap(iel4)
        endif
      endif
c
c     get final field file prefix
c
  290 format (a8)
  300 format(//,
     & ' Give a prefix for final field file (5chars): ')
      write (*,300)
      read (*,290) mypfx
c
c     get screen display mode
c
c  310 format(//,' Choose runtime screen display:',/,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     & '    A  : Standard display ',/,
c     & '    B  : Detailed Slab display.',/,
c     & '    C  : Full Display',/,
c     & '       : (Full slab display + timescales).',/,
c     & ' :: ')
c      write (*,310)
c      read (*,10) ilgg
c      ilgg=ilgg(1:1)
      vmod='MINI'
c
c      if (ilgg.eq.'a') ilgg='A'
c      if (ilgg.eq.'b') ilgg='B'
c      if (ilgg.eq.'c') ilgg='C'
c
c      if (ilgg.eq.'A') vmod='MINI'
c      if (ilgg.eq.'B') vmod='SLAB'
c      if (ilgg.eq.'C') vmod='FULL'
c
c     get runname
c
  320 format (a80)
  330 format(//,' Give a name/code for this run: ')
      write (*,330)
      read (*,320) runname
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine neqcheaders ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i,j
      character tab*1
      real*8 feldens
c
      tab=','
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Headers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   10   format(//' Flow Properties',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '   Velocity    :',1pg12.5,' km/s',/,
     &  '   Mach Number           :',1pg12.5,/,
     &  '   Mag. Alpha (Pmag/Pgas):',1pg12.5,//,
     &  '   ne  :',1pg12.5,' cm^-3',/,
     &  '   nH  :',1pg12.5,' cm^-3',/,
     &  '   d   :',1pg12.5,' g/cm^-3',/,
     &  '   T   :',1pg12.5,' K',/,
     &  '   B   :',1pg12.5,' microGauss',//,
     &  '   Pgas:',1pg12.5,' dyne/cm^2',/,
     &  '   Pmag:',1pg12.5,' dyne/cm^2',/,
     &  '   Pram:',1pg12.5,' dyne/cm^2',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
c     write ion balance files if requested
c
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (i.eq.1) ie=iel1
          if (i.eq.2) ie=iel2
          if (i.eq.3) ie=iel3
          if (i.eq.4) ie=iel4
          fn=' '
          pfx=elem(ie)
          sfx='csv'
          call newfile (pfx, elem_len(ie), sfx, 3, fn)
          filn(i)=fn(1:(elem_len(ie)+3+5))
        enddo
c
        do i=1,ieln
          open (luions(i),file=filn(i),status='NEW')
        enddo
c
        do i=1,ieln
          write (luions(i),*) fl,', ',runname
c
   20 format('Step [1], <X> [2], dX [3], <T> [4]',31(',',a6,'[',i2,']'))
          if (i.eq.1) then
            write (luions(i),20) (rom(j),j+4,j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),20) (rom(j),j+4,j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),20) (rom(j),j+4,j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),20) (rom(j),j+4,j=1,maxion(iel4))
          endif
c
        enddo
c
        do i=1,ieln
          close (luions(i))
        enddo
c
c     end trsmod = 'Y'
      endif
c
      de=feldens(dh,pop)
c
      fpb=' '
      pfx='bands'
      sfx='csv'
      call newfile (pfx, 5, sfx, 3, fn)
      fpb=fn(1:13)
      open (lupb,file=fpb,status='NEW')
c
c     Title
c
   30 format(' Nonequilibrium Cooling',/,
     &       ' ============================================',/,
     &       ' Diffuse Field, Full Continuum Calculations.',/,
     &       ' Calculated by MAPPINGS V ',a8)
      write (*,30) theversion
      write (luop,30) theversion
      write (lusp,30) theversion
      write (lupb,30) theversion
      if (ratmod.eq.'Y') write (lurt,30) theversion
      if (dynmod.eq.'Y') write (ludy,30) theversion
      if (allmod.eq.'Y') write (lual,30) theversion
      if (fclmod.eq.'Y') write (lucl,30) theversion
   40 format(//' Run  :',a,/
     &         ' File :',a)
      write (luop,40) runname,fl
      write (lusp,40) runname,fl
      write (lupb,40) runname,fl
      if (ratmod.eq.'Y') write (lurt,40) runname,fl
      if (dynmod.eq.'Y') write (ludy,40) runname,fl
      if (allmod.eq.'Y') write (lual,40) runname,fl
      if (fclmod.eq.'Y') write (lucl,40) runname,fl
c
c     Write Model Parameters
c
   50 format(//,' Model Parameters:',/,
     &        ' =================',//,
     &        ' Abundances     : ',a128,/,
     &        ' Pre-ionisation : ',a128,/,
     &        ' Photon Source  : ',a128)
   60 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12)
   70 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12,/,
     &        ' Electron Kappa : ',1pg11.4)
   80 format(//' ',t2,'Jden',t9,'Jgeo',t16,'Jtrans',
     &            t23,'Jend',t29,'Ielen',t35,
     &         'Jpoen',t43,'Fren',t51,'Tend',t57,'DIend',t67,
     &         'TAUen',t77,'Jeq',t84,'Teini')
   90 format('  ',4(a4,3x),2(i2,4x),0pf6.4,0pf6.0,
     &           2(1pg10.3),3x,a4,0pf7.1)
  100 format(//' Photon Source'/
     &         ' ============='/)
  110 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
  120 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
      write (luop,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
c
      write (lupb,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
c
      if (ratmod.eq.'Y') then
      write (lurt,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
      endif
c
      if (dynmod.eq.'Y') then
      write (ludy,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
      endif
c
      if (allmod.eq.'Y') then
      write (lual,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
      endif
c
      if (fclmod.eq.'Y') then
      write (lucl,10) vel0*1.d-5,machnumber,magparam,de,dh,rho0,te0,
     &bm0*1.d6,pr0,bp0,rho0*vel0*vel0
      endif
c
      write (luop,50) abnfile,ionsetup,srcfile
      write (lupb,50) abnfile,ionsetup,srcfile
      write (lusp,50) abnfile,ionsetup,srcfile
      if (ratmod.eq.'Y') write (lurt,50) abnfile,ionsetup,srcfile
      if (dynmod.eq.'Y') write (ludy,50) abnfile,ionsetup,srcfile
      if (allmod.eq.'Y') write (lual,50) abnfile,ionsetup,srcfile
      if (fclmod.eq.'Y') write (lucl,50) abnfile,ionsetup,srcfile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     abundances file header
c
c
      abundtitle=' Initial Abundances :'
      do i=1,atypes
        zi(i)=zion0(i)*deltazion(i)
      enddo
      call dispabundances (luop, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
c
      call wabund (lusp)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      cht='New Set'
      if (chargemode.eq.2) cht='Disabled'
      if (chargemode.eq.1) cht='Old Set'
      pht='Normal'
      if (photonmode.eq.0) pht='Zero Field, Decoupled'
      clt='Dere 2007 Collisions'
      if (collmode.eq.1) clt='Old A&R Collision Methods'
      if (collmode.eq.1) clt='LMS Collision Methods'
c
      if (usekappa) then
        write (luop,70) cht,pht,clt,kappa
      else
        write (luop,60) cht,pht,clt
      endif
      write (luop,80)
      write (luop,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,diend,
     &tauen,tmod,tm00
      write (luop,100)
      write (luop,110)
      write (luop,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      write (lupb,60) cht,pht,clt
      write (lupb,80)
      write (lupb,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,diend,
     &tauen,tmod,tm00
      write (lupb,100)
      write (lupb,110)
      write (lupb,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (ratmod.eq.'Y') then
        write (lurt,60) cht,pht,clt
        write (lurt,80)
        write (lurt,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lurt,100)
        write (lurt,110)
        write (lurt,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (dynmod.eq.'Y') then
        write (ludy,60) cht,pht,clt
        write (ludy,80)
        write (ludy,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (ludy,100)
        write (ludy,110)
        write (ludy,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (allmod.eq.'Y') then
c
        write (lual,60) cht,pht,clt
        write (lual,80)
        write (lual,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lual,100)
        write (lual,110)
        write (lual,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
      if (fclmod.eq.'Y') then
c
        write (lucl,60) cht,pht,clt
        write (lucl,80)
        write (lucl,90) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lucl,100)
        write (lucl,110)
        write (lucl,120) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
c     Shock parameters
c
  130 format(//,' Jump Conditions:',/,
     &         ' ================')
      write (luop,130)
      write (lupb,130)
      write (lusp,130)
      if (ratmod.eq.'Y') write (lurt,130)
      if (dynmod.eq.'Y') write (ludy,130)
      if (allmod.eq.'Y') write (lual,130)
      if (fclmod.eq.'Y') write (lucl,130)
c
  140 format(//' T1',1pg14.7,' V1',1pg14.7,
     &        ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7/
     &         ' T2',1pg14.7,' V2',1pg14.7,
     &        ' RH2',1pg14.7,' P2',1pg14.7,' B2',1pg14.7)
c
      dr=0.d0
      dv=vel1-vel0
c
      write (*,140) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (luop,140) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (lupb,140) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (lusp,140) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      if (ratmod.eq.'Y') write (lurt,140) te0,vel0,rho0,pr0,bm0,te1,
     &vel1,rho1,pr1,bm1
      if (dynmod.eq.'Y') write (ludy,140) te0,vel0,rho0,pr0,bm0,te1,
     &vel1,rho1,pr1,bm1
      if (fclmod.eq.'Y') write (lucl,140) te0,vel0,rho0,pr0,bm0,te1,
     &vel1,rho1,pr1,bm1
c
c
      close (lusp)
c
c
      wdilt0=te0
      t=te0
      ve=vel0
      dr=1.0
      dv=ve/100.d0
      fi=1.d0
      rad=1.d38
      if (wdil.eq.0.5) rad=0.d0
c
c     get the electrons...
c
      de=feldens(dh,pop)
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
c
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
      endif
      call cool (t, de, dh)
c
c
  150 format(//,' Precursor Conditions and Ionisation State',/
     &          ' =========================================',/)
      write (luop,150)
      write (lupb,150)
      if (fclmod.eq.'Y') write (lucl,150)
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
      wmod='FILE'
      call wmodel (luop, t, de, dh, dr, wmod)
      call wmodel (lupb, t, de, dh, dr, wmod)
      if (ratmod.eq.'Y') then
        call wmodel (lurt, t, de, dh, dr, wmod)
  160 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,'  en',t70,
     &     '  Dloss',t84,'  Eloss',t98,'  Egain',t112,'  Coll. H',t126,
     &     '  Chrg. Ex.',t140,'  Reson',t154,
     &     '  Xreson',t168,'  Inter/fine',t182,'  He int-forb',t196,
     &     '  Forbid',t210,'  Fe II',t224,'  2Photon',t238,
     &     '  Compton',t252,'  Free-Free',t266,'  Coll. Ion',t280,
     &     '  Photo.',t294,'  Recomb.',t308,'  Cosmic',t322,
     &     '  GrainsH-C',t336,'  G.Heat',t350,'  G.Cool',t364,
     &     '  PAHs')
        write (lurt,160)
      endif
c
      call wionabal (luop, pop)
      call wionabal (lupb, pop)
      if (allmod.eq.'Y') call wionabal (lual, pop)
      if (fclmod.eq.'Y') call wionabal (lucl, pop)
c
      close (luop)
  170 format(17(a12,a1))
      write (lupb,170) 'Te',tab,'de',tab,'dh',tab,'en',tab,'FHI',tab,'FH
     &II',tab,'mu',tab,'tloss',tab,'Lambda',tab,'ff/total',tab,'B0.0-0.1
     &keV',tab,'B0.1-0.5keV',tab,'B0.5-1.0keV',tab,'B1.0-2.0eV',tab,'B2.
     &0-10.0keV',tab,'Ball'
      close (lupb)
  180 format(//'Mean Zone Values'/
     &         '================'/)
  190 format(10(a12,a2),30(a12,a2),33(a12,a2))
  200 format(10(a12,a2),30(a12,a2),33(1pg12.5,a2))
c
  210 format( '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '======================================')
c
      if (fclmod.eq.'Y') then
        write (lucl,180)
c      write (lucl,560)
c      write (lucl,570)
        if (jnorm.eq.0) then
          write (lucl,190) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,
     & 'rho ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',
     & tab,'L/(ne.nH)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,
     & 'LHalpha',tab,'LLyalpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.1) then
          write (lucl,190) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,
     &'rho ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,
     &'L/(nH^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',
     &tab,'LLyalpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.2) then
          write (lucl,190) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,
     &'rho ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,
     &'L/(ne.ni)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',
     &tab,'LLyalpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.3) then
          write (lucl,190) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,
     &'rho ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,
     &'L/(n^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',
     &tab,'LLyalpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.4) then
          write (lucl,190) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,
     &'rho ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,
     &'L/(ne^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',
     &tab,'LLyalpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        write (lucl,200) '(K)',tab,'(/cm^3)',tab,'(/cm^3)',tab,
     &'(/cm^3)',tab,'(g/cm^3)',tab,' ',tab,' ',tab,'(amu)',tab,
     &'(erg/cm^3/s)',tab,'(erg cm^3/s)',tab,('(erg cm^3/s)',tab,
     &j=1,atypes),'(erg cm^3/s)',tab,'(erg cm^3/s)',tab,'(erg cm^3/s)'
c     &    (elem(j),tab,j=1,atypes)
        write (lucl,210)
      endif
      if (ratmod.eq.'Y') close (lurt)
      if (dynmod.eq.'Y') close (ludy)
      if (allmod.eq.'Y') close (lual)
      if (fclmod.eq.'Y') close (lucl)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
