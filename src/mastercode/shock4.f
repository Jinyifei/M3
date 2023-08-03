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
c     SHOCK4: Shock model with explicit Rankine-Hugoiot solution for
c     each step.  Time steps based on fastest atomic, cooling or
c     photon absorbion timescales.
c
c     Full continuum and diffuse field, requires precalculated ionisation
c     input.
c
c     RSS 1992
c     RSS 2010 minor bug fixes and fixed logical unit collisions
c     RSS 2013 improved setup for slow shocks, minor bug fixes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock4 ()
c
      include 'cblocks.inc'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           Variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 dr,drdw,drup,dvdw,dvup,frdw,ftime
      real*8 tdw,tlim,tup,xhii,bp0, bp1
      real*8 t, de, dh, dv, rad, fi, wdilt0
      real*8 tt0,l0,l11,l12,l1,cmpf,ue0,ue1
      real*8 r0,p0,v0, tnloss,cspd,wmol,mb
      real*8 hfree, bmag, ve, en,el,ue,rdis,dva
      real*8 pop1(mxion, mxelem),totn,ionbar,abstime
      real*8 press, mu, rhotot,cltime,utime
      real*8 ptime, ctime, rtime, htime,stime,eqtime
      real*8 atimefrac,cltimefrac,abtimefrac
      real*8 timend,teq,dheq,tl0,tl1,tt1
      real*8 band1,band2,band3,band4,band5,band6,bandall,wid
c
      real*8 meanion(mxelem),zi(mxelem)
c
      integer*4 ie,luop,lusp,lupb, lucl
      integer*4 count, step, np,i,j,lurt,ludy,lual,m, lupf
      integer*4 luions(4)
c
c     time steps - non-standard so commented out
c
c      integer*4 nt0,nt1,time
c      real tarray(2)
c      real dtime,dtarr
c
      character*24 abundtitle
      character fn*32
      character filn(4)*32
      character specmode*4, rmod*4, imod*4, tsrmod*4,spmod*4
      character fl*32, lmod*4, tab*1, ilgg*4, tmod*4, linemod*4
      character fd*32, fr*32, fa*32, fsp*32, allmod*4, ratmod*4,ispo*4
      character pfx*16, sfx*4, caller*4,wmod*4,pollfile*12,vmod*4
      character cht*32,pht*32,clt*32,dynmod*4, mypfx*16, model*32
      character fpf*64, fpb*64, fcl*64, fclmod*4
c
      logical iexi,pseudo
c
c
c           Functions
c
      real*8 fcolltim,feldens,fphotim,fpressu,frectim2,frho
      real*8 fmua
c
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
c      nt1 = time()
c      nt0 = nt1
c
      m=0
c
      ieln=4
      do i=1,ieln
        luions(i)=28+i
      enddo
c
      fl=' '
      pfx='shckn'
      sfx='sh4'
      call newfile (pfx, 4, sfx, 3, fn)
      fl=fn(1:12)
c
      open (luop,file=fl,status='NEW')
c
c
c     zero arrays and define modes: all ions, no on-thespot,
c     full continuum
c     calculations..
c
      call zer
      call zerbuf
c
      tab=','
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
      mypfx='psend'
      specmode='DW'
      jgeo='P'
      jtrans='LODW'
      jden='B'
      pseudo=.false.
      model='Shock 4'
c
   10 format(a)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTRO - Define Input Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   20 format(///,
     & ' SHOCK 4 model selected:',/,
     & ' Shock code with natural coordinates.',/,
     & ' ****************************************************')
      write (*,20)
c
c     get ionisation balance to work with...
c
      model='pre-ionisation'
      call popcha (model)
      model='Shock 4'
c
c     set up current field conditions
c
      call photsou (model)
c
c     Shock model preferences
c
c
   30 format(///,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Setting the shock conditions',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,30)
c
   40 format(//,
     & ' Choose Shock Jump Paramter:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    T : Shock jump in terms of Temperature jump',/,
     & '    V : Shock jump in terms of flow Velocity.',//,
     & ' :: ')
c
   50 write (*,40)
      read (*,10) ilgg
      ilgg=ilgg(1:1)
      if (ilgg.eq.'T') ilgg='t'
      if (ilgg.eq.'V') ilgg='v'
c
      if ((ilgg.ne.'t').and.(ilgg.ne.'v')) goto 50
c
c     Shock in term of flow velocity
c
      if (ilgg.eq.'v') then
   60   format(//,' Give Preshock Conditions:',/,
     &  '  T (K), dh (N), v (cm/s), B, Geo. dilution',/,
     &  ' (T<10 taken as a log, dh <= 0 taken as a log),',/,
     &  ' (B in micro Gauss, Dilution factor <= 0.5)',/,
     &  ' : ')
        write (*,60)
        read (*,*) t,dh,ve,bmag,wdil
c
        dr=1.d0
c
        if (t.le.10.d0) t=10.d0**t
        if (dh.le.0.d0) dh=10.d0**dh
        if (wdil.gt.0.5d0) wdil=0.5d0
        wdilt0=0.d0
c
        bmag=bmag*1.d-6
        vshoc=ve
        de=feldens(dh,pop)
c
c     fill in other default parameters
c
        ftime=0.d0
        tloss=0.d0
        rmod='NEAR'
c
        call rankhug (t, de, dh, ve, bmag, tloss, ftime)
        tm00=te0
c
      endif
c
c     Shock in terms of temp jump
c
      if (ilgg.eq.'t') then
c
   70   format(//,' Give Preshock Conditions:',/,
     &  '  (T<10 taken as a log, dh <= 0 taken as a log)',/,
     &  '  (B in micro Gauss, Dilution factor <= 0.5)',/,
     &  ' : T (K), dh (N), B, Geo. dilution',/,
     &  ' : ')
        write (*,70)
        read (*,*) te0,dh,bmag,wdil
c
        if (dr.le.0.d0) dr=1.d0
        if (wdil.gt.0.5d0) wdil=0.5d0
        bmag=bmag*1.d-6
c
   80   format(//,' Give Postshock Conditions:',/,
     &  '  (T<10 taken as a log)',//,
     &  ' : T (K)',/,
     &  ' : ')
        write (*,80)
        read (*,*) te1
c
        if (te0.le.10.d0) te0=10.d0**te0
        if (te1.le.10.d0) te1=10.d0**te1
c
        bm0=bmag
c
c        write(*,*) 'Calling velshock...'
        call velshock (dh, pop(2,1), te0, te1, bmag)
c        write(*,*) 'Back From velshock: te1 depo dhpo vel0'
c     &  , te1, depo, dhpo, vel0
c
        de0=feldens(dh,pop)
        de1=feldens(dhpo,pop)
        dh0=dh
        dh1=dhpo
        vel0=vshoc
        vel1=vpo
        rho0=frho(de0,dh0)
        rho1=frho(de1,dh1)
        pr0=fpressu(te0,dh0,pop)
        pr1=fpressu(te1,dh1,pop)
        bm1=bm0*rho1/rho0
        tm00=te0
      endif
c
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
c
   90   format(//' Shock Properties',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '   Velocity:',1pg12.5,' km/s',//,
     &  '   Preshock ne:',1pg12.5,' cm^-3',/,
     &  '   Preshock nH:',1pg12.5,' cm^-3',/,
     &  '   Preshock d :',1pg12.5,' g/cm^-3',/,
     &  '   Preshock P :',1pg12.5,' dyne/cm^2',/,
     &  '   Preshock T :',1pg12.5,' K',/,
     &  '   Preshock B :',1pg12.5,' microGauss',/,
     &  '   Preshock BP:',1pg12.5,' dyne/cm^2',//,
     &  '   Postshock ne:',1pg12.5,' cm^-3',/,
     &  '   Postshock nH:',1pg12.5,' cm^-3',/,
     &  '   Postshock d :',1pg12.5,' g/cm^-3',/,
     &  '   Postshock P :',1pg12.5,' dyne/cm^2',/,
     &  '   Postshock T :',1pg12.5,' K',/,
     &  '   Postshock B :',1pg12.5,' microGauss',/,
     &  '   Postshock BP:',1pg12.5,' dyne/cm^2',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,90) vel0*1d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,dhpo,
     &rho1,pr1,te1,bm1*1.d6,bp1
c
      if (alphacoolmode.eq.1) then
  100  format(//,
     &  ' Powerlaw Cooling enabled  Lambda T6 ^ alpha :',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Give index alpha, Lambda at 1e6K',/,
     &  '    (Lambda<0 as log)',/,
     &  ' :: ')
        write (*,100)
        read (*,*) alphaclaw,alphac0
        if (alphac0.lt.0.d0) alphac0=10.d0**alphac0
      endif
c
c     Diffuse field interaction
c
c
      photonmode=1
c
  110 format(//,
     & ' Choose Diffuse Field Option :',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    F  : Full diffuse field interaction (Default).',/,
     & '    Z  : Zero diffuse field interaction.',/,
     & ' :: ')
      write (*,110)
  120 read (*,10) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'z') ilgg='Z'
      if (ilgg.eq.'f') ilgg='F'
c
      if ((ilgg.ne.'Z').and.(ilgg.ne.'F')) goto 120
c
      if (ilgg.eq.'Z') photonmode=0
      if (ilgg.eq.'F') photonmode=1
c
  130 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Calculation Settings ',/,
     & ' ,::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,130)
c
c     Set up time step limits
c
  140 format(//,' Choose a Time Step Behaviour :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Auto. time steps, based on Atomic timescales.',/,
     & '    M  : Choose a Maximum time step, Auto otherwise.',/,
     & '    N  : Choose a miNimum time step, Auto otherwise.',/,
     & '    S  : Strict preset timestep. (not recomended)',//,
     & ' :: ')
c
  150 write (*,140)
      read (*,10) tmod
      tmod=tmod(1:1)
      if (tmod.eq.'a') tmod='A'
      if (tmod.eq.'m') tmod='M'
      if (tmod.eq.'n') tmod='N'
      if (tmod.eq.'s') tmod='S'
c
      if ((tmod.ne.'A').and.(tmod.ne.'M').and.(tmod.ne.'N')
     &.and.(tmod.ne.'S')) goto 150
c
c
      atimefrac=0.2d0
      cltimefrac=0.2d0
      abtimefrac=0.2d0
c
      if (tmod.ne.'S') then
  160    format(//,
     &   ' Give timescale fractions:',/,
     &   ' ( 0 < dtau < 1 )',/,
     &   ' ( recommend: 0.05 < dtau < 0.25 )',/,
     &   ' : atomic, cooling, absorbsion',/,
     &   ' : ')
  170   write (*,160)
        read (*,*) atimefrac,cltimefrac,abtimefrac
c
        if (atimefrac.le.0.d0) goto 170
        if (cltimefrac.le.0.d0) goto 170
        if (abtimefrac.le.0.d0) goto 170
c
        if (atimefrac.ge.1.d0) goto 170
        if (cltimefrac.ge.1.d0) goto 170
        if (abtimefrac.ge.1.d0) goto 170
c
      endif
c
      utime=0.d0
      if (tmod.ne.'A') then
  180   format(//,
     &  ' Give the time step:',/,
     &  ' ( time < 100 taken as a log)',//,
     &  ' : time (sec)',/,
     &  ' : ')
        write (*,180)
        read (*,*) utime
        if (utime.le.100.d0) utime=10.d0**utime
      endif
c
c
  190 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Boundry Conditions ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,190)
c
c     Choose ending
c
  200 format(/,
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
  210 write (*,200)
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
     &.and.(jend.ne.'G')) goto 210
c
c     Secondary info:
c
      if (jend.eq.'B') then
  220   format(//,' Give Atom, Ion and limit fraction: ')
        write (*,220)
        read (*,*) ielen,jpoen,fren
      endif
      if (jend.eq.'C') then
  230   format(//,' Give final temperature (K > 10, log <= 10): ')
        write (*,230)
        read (*,*) tend
        if (tend.le.10.d0) tend=10.d0**tend
      endif
      if (jend.eq.'D') then
  240   format(//,' Give final distance (cm > 100, log<=100): ')
        write (*,240)
        read (*,*) diend
        if (diend.le.100.d0) diend=10.d0**diend
      endif
      if (jend.eq.'E') then
  250   format(//,' Give time limit (s > 100, log<=100): ')
        write (*,250)
        read (*,*) timend
        if (timend.le.100.d0) timend=10.d0**timend
      endif
c
  260 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Output Requirements ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,260)
c
c     default monitor elements
c
      iel1=zmap(1)
      iel2=zmap(2)
      iel3=zmap(6)
      iel4=zmap(8)
      mypfx='psend'
c
  270 format(//,' Choose output settings : ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard output. ',/,
     & '    B  : Standard plus ion balance files.',/,
     & '    C  : Standard plus dynamics file.',/,
     & '    D  : Standard plus rates file.',/,
     & '    E  : Standard plus final down stream field.',/,
     & '    F  : Standard plus upstream field at each step.',/,
     & '    J  : Cooling Components by Elements file.',/,
     & '       :',/,
     & '    G  : B + E',/,
     & '    H  : B + F',/,
     & '       :',/,
     & '    I  : Everything',//,
     & ' :: ')
  280 write (*,270)
      read (*,10) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
      if (ilgg.eq.'d') ilgg='D'
      if (ilgg.eq.'e') ilgg='E'
      if (ilgg.eq.'f') ilgg='F'
      if (ilgg.eq.'g') ilgg='G'
      if (ilgg.eq.'h') ilgg='H'
      if (ilgg.eq.'i') ilgg='I'
      if (ilgg.eq.'j') ilgg='J'
c
      if ((ilgg.ne.'A').and.(ilgg.ne.'B').and.(ilgg.ne.'C')
     &.and.(ilgg.ne.'D').and.(ilgg.ne.'E').and.(ilgg.ne.'F')
     &.and.(ilgg.ne.'G').and.(ilgg.ne.'H').and.(ilgg.ne.'I')
     &.and.(ilgg.ne.'J')) goto 280
c
c     Upstream fields
c
      lmod='N'
      if ((ilgg.eq.'F').or.(ilgg.eq.'H').or.(ilgg.eq.'I')) then
        lmod='Y'
      endif
c
c     F-Lambda plots, ionising and optical
c
      jspec='NO'
      if ((ilgg.eq.'E').or.(ilgg.eq.'G').or.(ilgg.eq.'I')) then
        jspec='YES'
      endif
c
c     timescales and rates file
c
      ratmod='N'
      if ((ilgg.eq.'D').or.(ilgg.eq.'I')) then
        ratmod='Y'
        pfx='rates'
        sfx='sh4'
        call newfile (pfx, 5, sfx, 3, fn)
        fr=fn(1:13)
        open (lurt,file=fr,status='NEW')
      endif
c
c     Cooling Components File
c
      fclmod='N'
      if ((ilgg.eq.'J').or.(ilgg.eq.'I')) then
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
        sfx='sh4'
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
      if ((ilgg.eq.'B').or.(ilgg.eq.'G')
     &    .or.(ilgg.eq.'H').or.(ilgg.eq.'I')) then
c
        if ((ilgg.eq.'B').or.(ilgg.eq.'G').or.(ilgg.eq.'H')) then
c
  290    format(//,' Record all ions file (Y/N)?: ')
          write (*,290)
          read (*,10) allmod
          allmod=allmod(1:1)
          if (allmod.eq.'y') allmod='Y'
        endif
c
        if (ilgg.eq.'I') allmod='Y'
c
        if (allmod.eq.'Y') then
          pfx='allion'
          sfx='sh4'
          call newfile (pfx, 6, sfx, 3, fn)
          fa=fn(1:14)
          open (lual,file=fa,status='NEW')
        endif
c
        if (ilgg.eq.'B') then
  300    format(//,' Record particular ion balances (Y/N)?: ')
          write (*,300)
          read (*,10) tsrmod
          tsrmod=tsrmod(1:1)
          if (tsrmod.eq.'y') tsrmod='Y'
        endif
c
        if (ilgg.eq.'I') tsrmod='Y'
c
        if (tsrmod.eq.'Y') then
  310       format(//,' Give 4 elements by atomic number: ')
          write (*,310)
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
  320 format (a8)
  330 format(//,
     & ' Give a prefix for final field file (5chars): ')
      write (*,330)
      read (*,320) mypfx
c
c     get screen display mode
c
  340 format(//,' Choose runtime screen display:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard display ',/,
     & '    B  : Detailed Slab display.',/,
     & '    C  : Full Display',/,
     & '       : (Full slab display + timescales).',/,
     & ' :: ')
      write (*,340)
      read (*,10) ilgg
      ilgg=ilgg(1:1)
      vmod='MINI'
c
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
c
      if (ilgg.eq.'A') vmod='MINI'
      if (ilgg.eq.'B') vmod='SLAB'
      if (ilgg.eq.'C') vmod='FULL'
c
c     get runname
c
  350 format (a80)
  360 format(//,' Give a name/code for this run: ')
      write (*,360)
      read (*,350) runname
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Headers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
  370 format('Step [1], <X> [2], dX [3], <T> [4]',31(',',a6,'[',i2,']'))
          if (i.eq.1) then
            write (luions(i),370) (rom(j),j+4,j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),370) (rom(j),j+4,j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),370) (rom(j),j+4,j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),370) (rom(j),j+4,j=1,maxion(iel4))
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
      fpb=' '
      pfx='bands'
      sfx='sh4'
      call newfile (pfx, 5, sfx, 3, fn)
      fpb=fn(1:13)
      open (lupb,file=fpb,status='NEW')
c
c     Title
c
  380 format(' SHOCK4: Explicit Rankine-Hugoniot shock code: ',/,
     &       ' ============================================',/,
     &       ' Diffuse Field, Full Continuum Calculations.',/,
     &       ' Calculated by MAPPINGS V ',a8)
      write (*,380) theversion
      write (luop,380) theversion
      write (lusp,380) theversion
      write (lupb,380) theversion
      if (ratmod.eq.'Y') write (lurt,380) theversion
      if (dynmod.eq.'Y') write (ludy,380) theversion
      if (allmod.eq.'Y') write (lual,380) theversion
      if (fclmod.eq.'Y') write (lucl,380) theversion
  390 format(//' Run  :',a,/
     &         ' File :',a)
      write (luop,390) runname,fl
      write (lusp,390) runname,fl
      write (lupb,390) runname,fl
      if (ratmod.eq.'Y') write (lurt,390) runname,fl
      if (dynmod.eq.'Y') write (ludy,390) runname,fl
      if (allmod.eq.'Y') write (lual,390) runname,fl
      if (fclmod.eq.'Y') write (lucl,390) runname,fl
c
c     Write Model Parameters
c
  400 format(//,' Model Parameters:',/,
     &        ' =================',//,
     &        ' Abundances     : ',a128,/,
     &        ' Pre-ionisation : ',a128,/,
     &        ' Photon Source  : ',a128)
  410 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12)
  420 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12,/,
     &        ' Electron Kappa : ',1pg11.4)
  430 format(//' ',t2,'Jden',t9,'Jgeo',t16,'Jtrans',
     &            t23,'Jend',t29,'Ielen',t35,
     &         'Jpoen',t43,'Fren',t51,'Tend',t57,'DIend',t67,
     &         'TAUen',t77,'Jeq',t84,'Teini')
  440 format('  ',4(a4,3x),2(i2,4x),0pf6.4,0pf6.0,
     &           2(1pg10.3),3x,a4,0pf7.1)
  450 format(//' Photon Source'/
     &         ' ============='/)
  460 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
  470 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
      write (luop,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &dhpo,rho1,pr1,te1,bm1*1.d6,bp1
c
      write (lupb,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &dhpo,rho1,pr1,te1,bm1*1.d6,bp1
c
      if (ratmod.eq.'Y') then
        write (lurt,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &   dhpo,rho1,pr1,te1,bm1*1.d6,bp1
      endif
c
      if (dynmod.eq.'Y') then
        write (ludy,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &   dhpo,rho1,pr1,te1,bm1*1.d6,bp1
      endif
c
      if (allmod.eq.'Y') then
        write (lual,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &   dhpo,rho1,pr1,te1,bm1*1.d6,bp1
      endif
c
      if (fclmod.eq.'Y') then
        write (lucl,90) vel0*1.d-5,de,dh,rho0,pr0,te0,bm0*1.d6,bp0,de1,
     &   dhpo,rho1,pr1,te1,bm1*1.d6,bp1
      endif
c
      write (luop,400) abnfile,ionsetup,srcfile
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
      write (lupb,400) abnfile,ionsetup,srcfile
      if (ratmod.eq.'Y') write (lurt,400) abnfile,ionsetup,srcfile
      if (dynmod.eq.'Y') write (ludy,400) abnfile,ionsetup,srcfile
      if (allmod.eq.'Y') write (lual,400) abnfile,ionsetup,srcfile
      if (fclmod.eq.'Y') write (lucl,400) abnfile,ionsetup,srcfile
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
        write (luop,420) cht,pht,clt,kappa
      else
        write (luop,410) cht,pht,clt
      endif
      write (luop,430)
      write (luop,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &diend,tauen,tmod,tm00
      write (luop,450)
      write (luop,460)
      write (luop,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      write (lupb,410) cht,pht,clt
      write (lupb,430)
      write (lupb,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &diend,tauen,tmod,tm00
      write (lupb,450)
      write (lupb,460)
      write (lupb,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (ratmod.eq.'Y') then
        write (lurt,410) cht,pht,clt
        write (lurt,430)
        write (lurt,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lurt,450)
        write (lurt,460)
        write (lurt,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (dynmod.eq.'Y') then
        write (ludy,410) cht,pht,clt
        write (ludy,430)
        write (ludy,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (ludy,450)
        write (ludy,460)
        write (ludy,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (allmod.eq.'Y') then
c
        write (lual,410) cht,pht,clt
        write (lual,430)
        write (lual,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lual,450)
        write (lual,460)
        write (lual,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
      if (fclmod.eq.'Y') then
c
        write (lucl,410) cht,pht,clt
        write (lucl,430)
        write (lucl,440) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lucl,450)
        write (lucl,460)
        write (lucl,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
c     Shock parameters
c
  480 format(//,' Jump Conditions:',/,
     &         ' ================')
      write (luop,480)
      write (lupb,480)
      write (lusp,480)
      if (ratmod.eq.'Y') write (lurt,480)
      if (dynmod.eq.'Y') write (ludy,480)
      if (allmod.eq.'Y') write (lual,480)
      if (fclmod.eq.'Y') write (lucl,480)
c
  490 format(//' T1',1pg14.7,' V1',1pg14.7,
     &        ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7/
     &         ' T2',1pg14.7,' V2',1pg14.7,
     &        ' RH2',1pg14.7,' P2',1pg14.7,' B2',1pg14.7)
c
      dr=0.d0
      dv=vel1-vel0
c
      write (*,490) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (luop,490) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (lupb,490) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      write (lusp,490) te0,vel0,rho0,pr0,bm0,te1,vel1,rho1,pr1,bm1
      if (ratmod.eq.'Y') write (lurt,490) te0,vel0,rho0,pr0,bm0,te1,
     &vel1,rho1,pr1,bm1
      if (dynmod.eq.'Y') write (ludy,490) te0,vel0,rho0,pr0,bm0,te1,
     &vel1,rho1,pr1,bm1
      if (fclmod.eq.'Y') write (lucl,490) te0,vel0,rho0,pr0,bm0,te1,
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
  500 format(//,' Precursor Conditions and Ionisation State',/
     &          ' =========================================',/)
      write (luop,500)
      write (lupb,500)
      if (fclmod.eq.'Y') write (lucl,500)
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
  510 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,'  en',t70,
     &     '  Dloss',t84,'  Eloss',t98,'  Egain',t112,'  Coll. H',t126,
     &     '  Chrg. Ex.',t140,'  Reson',t154,
     &     '  Xreson',t168,'  Inter/fine',t182,'  He int-forb',t196,
     &     '  Forbid',t210,'  Fe II',t224,'  2Photon',t238,
     &     '  Compton',t252,'  Free-Free',t266,'  Coll. Ion',t280,
     &     '  Photo.',t294,'  Recomb.',t308,'  Cosmic',t322,
     &     '  GrainsH-C',t336,'  G.Heat',t350,'  G.Cool',t364,
     &     '  PAHs')
        write (lurt,510)
      endif
c
      call wionabal (luop, pop)
      call wionabal (lupb, pop)
      if (allmod.eq.'Y') call wionabal (lual, pop)
      if (fclmod.eq.'Y') call wionabal (lucl, pop)
c
      close (luop)
  520 format(17(a12,a1))
      write (lupb,520) 'Te',tab,'de',tab,'dh',tab,'en',tab,'FHI',tab,'FH
     &II',tab,'mu',tab,'Lambda',tab,'tloss',tab,'ff',tab,'B0.5-2.0keV',
     &tab,'B2.0-8.0keV',tab,'B0.5-8.0keV',tab,'B3.0-10.0eV',tab,'B1.0-2.
     &0keV',tab,'B2.0-3.0keV',tab,'Ball'
      close (lupb)
  530 format(//'Mean Zone Values'/
     &         '================'/)
  540 format(10(a12,a1),30(a12,a1))
  550 format(10(a12,a1),30(a12,a1),30(1pg12.5,a1))
  560 format( '                                                ',
     & '                             ',
     & '                                                ',
     & '                             ',
     & '        ===================== Lambda Cooling Function',
     & ' and Element Components,',
     & ' Components Scale Approximately Linearly with Component',
     & ' Abundances ==========',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '==============================')
  570 format( '==================================================',
     & '===========================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '=================================================',
     & '============================',
     & '============================================')
      if (fclmod.eq.'Y') then
        write (lucl,530)
        write (lucl,560)
        write (lucl,570)
        write (lucl,540) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho '
     &   ,tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(n
     &_H^2)',tab,(elem(j),tab,j=1,atypes)
        write (lucl,550) '(K)',tab,'(/cm^3)',tab,'(/cm^3)',tab,'(/cm^3)'
     &   ,tab,'(g/cm^3)',tab,' ',tab,' ',tab,'(amu)',tab,'(erg/cm^3/s)',
     &   tab,'(erg cm^3/s)',tab,('(erg cm^3/s)',tab,j=1,atypes),
     &   (atwei(j),tab,j=1,atypes)
        write (lucl,570)
      endif
      if (ratmod.eq.'Y') close (lurt)
      if (dynmod.eq.'Y') close (ludy)
      if (allmod.eq.'Y') close (lual)
      if (fclmod.eq.'Y') close (lucl)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Begin Main Calculation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      wdilt0=te1
      t=te1
      cmpf=rho1/rho0
      dh=dh*cmpf
      de=feldens(dh,pop)
      dr=1.d0
      ve=vel1
      vpo=vel1
      bmag=bm1
      bm0=bm1
      dv=ve*0.01d0
      fi=1.d0
      rad=1.d38
      if (wdil.eq.0.5d0) rad=0.d0
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
      endif
c
      call cool (t, de, dh)
c
      en=zen*dh
c
      ue=1.5d0*(en+de)*rkb*t
      tnloss=tloss/(en*de)
c
      cspd=dsqrt(1.666666666666d0*pr1/rho1)
      wmol=rho1/(en+de)
      mu=fmua(de,dh)
c
      step=0
c
  580 format(//,a4,a1,14(4x,a6,4x,a1))
  590 format(i4.3,a1,14(1pg13.6,a1))
c
  600 format(//' #',t9,'Te',t21,'ne',t33,'nH',t45,'ni',t57,
     &     'B',t69,'FHI'/t8,'Time',t22,'dt',t32,'Dist',t45,
     &     'dr',t57,'v',t67,'dlos'/)
  610 format (i4,6(1x,1pg11.4)/t4,6(1x,1pg11.4))
      open (luop,file=fl,status='OLD',access='APPEND')
      write (luop,580) 'Step',tab,'Dist.',tab,'Time',tab,'Te',tab,'de',
     &tab,'dh',tab,'en',tab,'mu',tab,'Lambda',tab,'Tloss',tab,'dlos',
     &tab,'FHI',tab,'FHII',tab,'Bbar',tab,'dr',tab
      write (luop,590) step,tab,dist(step+1),tab,timlps(step+1),tab,t,
     &tab,de,tab,dh,tab,en,tab,mu,tab,tloss/(dh*dh),tab,tloss,tab,dlos,
     &tab,pop(1,1),tab,pop(2,1),tab,bm0,tab,dr,tab
      close (luop)
c
      if (vmod.eq.'MINI') then
        write (*,600)
        write (*,610) step,t,de,dh,en,bm0,pop(1,1),timlps(step+1),ftime,
     &   dist(step+1),dr,veloc(step+1),dlos
      endif
c
      en=zen*dh
      press=(en+de)*rkb*t
      ue=(1.5d0)*(en+de)*rkb*t
c
c     effective cooling timescale, based on net loss
c
      cltime=dabs(ue/tloss)
c
      rhotot=frho(de,dh)
      wmol=rhotot/(en+de)
      mu=fmua(de,dh)
c
      hfree=bmag*bmag/(4.d0*pi*en)+0.5d0*(wmol*ve*ve)+2.5d0*(rkb*t)/mu
c
c
      htime=dabs(hfree/tloss)
c
      abstime=1.d0/epsilon
c
      ctime=fcolltim(de)
      rtime=frectim2(de)
      if (photonmode.ne.0) then
        ptime=fphotim()
        eqtime=(1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
        stime=(1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      else
        eqtime=(1.d0/ctime)+(1.d0/rtime)
        stime=(1.d0/ctime)-(1.d0/rtime)
      endif
c
      eqtime=dabs(1.d0/eqtime)
      stime=dabs(1.d0/stime)
c
      ftime=1.d0/((1.d0/(stime*atimefrac))+(1.d0/(cltime*cltimefrac)))
c     &             +(1.d0/(abstime*abtimefrac)))
c
      if ((tmod.eq.'M').and.(ftime.gt.utime)) ftime=utime
      if ((tmod.eq.'N').and.(ftime.lt.utime)) ftime=utime
      if (tmod.eq.'S') ftime=utime
c
  620 format ( //,' Tot. Internal :',1pg14.7,/,
     &     ' Tot.loss rate :',1pg14.7,/' Eff.loss rate :',1pg14.7,/,
     &     ' Enthaply time :',1pg14.7,/' Eff.Cool time :',1pg14.7,/,
     &     ' Collis.  time :',1pg14.7,/' Recomb.  time :',1pg14.7,/,
     &     ' Photo.   time :',1pg14.7,/' Atomic   time :',1pg14.7,/,
     &     ' Equilib  time :',1pg14.7,/' Absorb.  time :',1pg14.7,/,
     &     ' Choice   time :',1pg14.7/)
c
      if (vmod.eq.'FULL') then
        write (*,620) ue,tloss,tloss,htime,cltime,ctime,rtime,ptime,
     &   stime,eqtime,abstime,ftime
      endif
c
c      if (ratmod.eq.'Y') then
c        open (lurt,file=fr,status='OLD',access='APPEND')
c        write (lurt,610) ue,tloss,tloss,htime,cltime,ctime,rtime,ptime,
c     &   stime,eqtime,abstime,ftime
c        close (lurt)
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Setup Initial first step in post shock-front gas(#1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rmod='NEAR'
      call rankhug (t, de, dh, ve, bmag, tloss, ftime)
c
      dr=ftime*(vel0+vel1)/2.d0
      dv=vel1-vel0
c
c      write(*,200) tloss,dr,dv
c
      step=0
c
c     init arrays first step
c
      xh(1)=0.d0
      xh(2)=pop(2,1)
      veloc(1)=0.d0
      veloc(2)=vel0
      te(1)=0.d0
      te(2)=te0
      dhy(1)=0.d0
      dhy(2)=dh
      dist(1)=0.0d0
      dist(2)=0.0d0
      deel(1)=0.d0
      deel(2)=de
      timlps(1)=0.0d0
c
c     prepare for time step in ions
c
c
c     record initial conditions, t0 and l0 and pop
c
      call copypop (pop, pop0)
c
      tt0=te0
      l0=tloss
      l12=0.d0
      l11=0.d0
      r0=rho0
      p0=pr0
      v0=vel0
      de0=de
      dh0=dh
      count=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SHOCK Iteration reentry point:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  630 count=count+1
c
      call copypop (pop0, pop)
c
      cmpf=rho1/rho0
      dh=(dh0+dh1)*0.5d0
      de=feldens(dh,pop)
      t=(te1+te0)*0.5d0
      xhii=pop(2,1)
c
      call timion (t, de, dh, xhii, ftime)
c
      call copypop (pop, pop1)
c
c     after timion on the average t,de,dh, now
c     calculate new lossrate at te1,de1,dh1
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
      rad=dist(step+1)
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
      endif
c
      t=te1
      dh=dh1
      de=feldens(dh,pop)
c
      call cool (t, de, dh)
c
c     record first guess
c
      l12=l11
      l11=tloss
      l1=dabs(1.d0-dabs(l11)/(epsilon+dabs(l12)))
c
c
c    test for pseudo equilibrium
c
c      if ((eqtime.le.stime).and.
c     &   (eqtime.le.cltime).and.
c     &   (eqtime.le.abstime*abtimefrac).and.
      if ((rtime.le.ctime).and.(ptime.le.ctime).and.(dlos.lt.0.01d0))
     &then
        if (count.gt.1) then
          te1=te0
          dh1=dh0
          vel1=vel0
          bm1=bm0
        endif
        goto 640
      endif
c
      tloss=(l11+l0)*0.5d0
c
c     allow convergence to relax as dlos->0
c
      tlim=(1.d-5)/(dabs(dlos)+epsilon)
c
c     if within limits go and write step
c
      if (l1.le.tlim) then
        write (*,*) '    Converged',l1
        call copypop (pop, pop1)
        call averinto (0.5d0, pop0, pop1, pop)
        t=(te0+te1)*0.5d0
        dh=(dh0+dh1)*0.5d0
        de=feldens(dh,pop)
        call cool (t, de, dh)
        goto 650
      endif
c
c     try 8 times then give up
c
      if (count.gt.7) then
c         count = 0
c         ftime = ftime*0.5d0
        write (*,*) '    Failed to converge',l1
        call copypop (pop, pop1)
        call averinto (0.5d0, pop0, pop1, pop)
        t=(te0+te1)*0.5d0
        dh=(dh0+dh1)*0.5d0
        de=feldens(dh,pop)
        call cool (t, de, dh)
        goto 650
      endif
c
      en=zen*dh+de
      press=en*rkb*t
      ue=1.5d0*press
c
      cltime=dabs(ue/tloss)
c
      if (ftime.gt.cltime) then
        ftime=cltime
      endif
c
      t=tt0
      dh=dh0
      de=feldens(dh,pop0)
      ve=v0
      bmag=bm0
c
      call copypop (pop1, pop)
c
      rmod='NEAR'
      call rankhug (t, de, dh, ve, bmag, tloss, ftime)
c
      dr=ftime*(vel0+vel1)*0.5d0
      dv=vel1-vel0
c
c 200  format(t4,4(1pg15.7,x))
c      write(*,200) tloss,ve,dv,vel1
c
c     otherwise continue with the normal NEQ
c     mode
c
      goto 630
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pseudo equilibrium calculations for case where
c     heating and cooling very nearly match and extrapolation
c     based on net cooling rate is unstable....
c     = photoionized post-shock zone
c
c     pseudo equilibrium entry point
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  640 continue
c
      write (*,*) 'Equilibrium'
c
      pseudo=.true.
c
      count=count+1
c
      call copypop (pop0, pop)
c
      dh=dh1
      de=feldens(dh,pop)
c
c     initial internal energy /unit mass
c
      en=zen*dh
      ue0=(1.5d0)*(en+de)*rkb*te1
c
      tt0=te1
      t=te1
c
      imod='TIM'
      call teequi2 (tt0, t, de, dh, ftime, imod)
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
c     compute new internal energy
c     de is updated
c
      teq=t
      dheq=dh
c
      en=zen*dheq
      ue1=(1.5d0)*(en+de)*rkb*t
c
c     Net Loss
c
      tloss=(ue0-ue1)/ftime
c
c      if (ratmod.eq.'Y') then
c        open (lurt,file=fr,status='OLD',access='APPEND')
c        write (lurt,*) 'pse',ue0,ue1,tloss,tt0,t
c        close (lurt)
c      endif
c
      call copypop (pop, pop1)
      call copypop (pop0, pop)
      t=tt0
      dh=dh1
      de=feldens(dh,pop)
      ve=vel1
      bmag=bm1
      call copypop (pop1, pop)
c
      rmod='NEAR'
      call rankhug (t, de, dh, ve, bmag, tloss, ftime)
c
      if (dabs(1.d0-(teq/te1)).gt.1.d-3) then
        tl0=tloss*0.5
        bmag=bm1
        call rankhug (t, de, dh, ve, bmag, tl0, ftime)
        tt1=te1
        tl1=tloss*0.5
        bmag=bm1
        call rankhug (t, de, dh, ve, bmag, tl1, ftime)
        tl1=tl0+(((teq-tt1)/(te1-tt1))*(tl1-tl0))
        bmag=bm1
        call rankhug (t, de, dh, ve, bmag, tl1, ftime)
        tloss=tl1
      endif
c
      dr=ftime*(vel0+vel1)*0.5d0
      dv=vel1-vel0
c
      call averinto (0.5d0, pop0, pop1, pop)
      t=(te0+te1)*0.5d0
      dh=(dh0+dh1)*0.5d0
      de=feldens(dh,pop)
c
      call cool (t, de, dh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c End step iterations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  650 continue
c
      step=step+1
c
c     record step
c
      xh(step+1)=pop(2,1)
      veloc(step+1)=(vel1+vel0)*0.5d0
      te(step+1)=t
      dhy(step+1)=dh
      deel(step+1)=de
      dist(step+1)=dist(step)+dr
      timlps(step+1)=timlps(step)+ftime
c
c
  660 format(i4,1x,i4)
  670 format(1pg14.7,'K',1x,1pg14.7,'cm/s',1x,1pg14.7,'g/cm3',1x,
     &1pg14.7,'dyne/cm2',1x,1pg14.7,'Gauss',1x,1pg14.7,'ergs/cm3/s')
  680 format(2(1pg14.7,' cm',1x),2(1pg14.7,' s',1x))
  690 format(2(1pg14.7,' ergs/cm^3',1x),1pg14.7,' /cm^3',1x,
     &1pg14.7,' cm/s',1x,1pg14.7,' g ')
      en=zen*dh
c
      ue=1.5d0*(en+de)*rkb*t
      tnloss=tloss/(en*de)
c
      cspd=dsqrt(1.666666666666d0*(pr1+pr0)/(rho1+rho0))
      wmol=(rho1+rho0)/(2.d0*(en+de))
c
      mb=(bm0+bm1)*0.5
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
c
  700    format('     Case A-B (HI, HeII): ',2(1pg11.3))
      write (*,700) caseab(1),caseab(2)
c
      if (vmod.eq.'MINI') then
        write (*,610) step,t,de,dh,en,mb,pop(1,1),timlps(step+1),ftime,
     &   dist(step+1),dr,veloc(step+1),dlos
      endif
c
c         dtarr = dtime(tarray)
c         nt1   = time()
c
c
c         write(*,900) tarray(1),nt1-nt0
c
c         nt0 = nt1
c
c
      if (ratmod.eq.'Y') then
        open (lurt,file=fr,status='OLD',access='APPEND')
        wmod='LOSS'
        call wmodel (lurt, t, de, dh, dr, wmod)
        close (lurt)
      endif
c
      if (dynmod.eq.'Y') then
        open (ludy,file=fd,status='OLD',access='APPEND')
        write (ludy,660) step,count
        write (ludy,680) dist(step+1),dr,timlps(step+1),ftime
        write (ludy,670) te0,vel0,rho0,pr0,bm0,l0
        write (ludy,670) te1,vel1,rho1,pr1,bm1,l11
        write (ludy,690) ue,hfree,en,cspd,wmol
        close (ludy)
      endif
c
      open (luop,file=fl,status='OLD',access='APPEND')
      write (luop,590) step,tab,dist(step+1),tab,timlps(step+1),tab,t,
     &tab,de,tab,dh,tab,en,tab,mu,tab,tloss/(dh*dh),tab,tloss,tab,dlos,
     &tab,pop(1,1),tab,pop(2,1),tab,(0.5*(bm0+bm1)),tab,dr,tab
      close (luop)
  710 format(70(1pg12.5,a1))
      if (fclmod.eq.'Y') then
        do i=1,atypes
c
          meanion(i)=0.d0
c
          ionbar=0.d0
          do j=1,maxion(i)
            ionbar=ionbar+(pop(j,i)*j)
          enddo
c
          meanion(i)=dmax1(ionbar-1.0d0,0.0d0)
c
        enddo
c
        open (lucl,file=fcl,status='OLD',access='APPEND')
c
        write (lucl,710) t,tab,de,tab,dh,tab,en,tab,(0.5*(rho0+rho1)),
     &   tab,pop(1,1),tab,pop(2,1),tab,mu,tab,tloss,tab,tloss/(dh*dh),
     &   (tab,coolz(j)/(dh*dh),j=1,atypes),(tab,meanion(j),j=1,atypes)
        close (lucl)
c
      endif
c
      if (allmod.eq.'Y') then
        open (lual,file=fa,status='OLD',access='APPEND')
        write (lual,720)
        write (lual,730) step,t,de,dh,vel1,rho1,pr1,dist(step+1),
     &   timlps(step+1)
  720 format(//,' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  V1(cm/s)    ',
     &                ' Rho1(g/cm^3) ',
     &                ' Pr1(erg/cm^3)',
     &                '   Dist.(cm)  ',
     &                ' Elps. Time(s)')
  730 format(1x,i4,8(1pg14.7)/)
c
        call wionabal (lual, pop)
        call wionabal (lual, popint)
        close (lual)
      endif
c
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          open (luions(i),file=filn(i),status='OLD',access='APPEND')
c
  740 format(1x,i4,3(', ',1pg12.5),31(', ',1pg12.5))
          if (i.eq.1) then
            write (luions(i),740) step,dist(step)+(dr*0.5),dr,t,(pop(j,
     &       iel1),j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),740) step,dist(step)+(dr*0.5),dr,t,(pop(j,
     &       iel2),j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),740) step,dist(step)+(dr*0.5),dr,t,(pop(j,
     &       iel3),j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),740) step,dist(step)+(dr*0.5),dr,t,(pop(j,
     &       iel4),j=1,maxion(iel4))
          endif
c
          close (luions(i))
c
        enddo
c
c     end .tsr files
c
      endif
c
c
c     poll for photons file
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
c
      if ((lmod.eq.'Y').or.(iexi)) then
c
c     Local photon field
c
        specmode='LOCL'
        call localem (t, de, dh)
        rad=dist(step+1)
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
        caller='S4'
        pfx='plocl'
        np=5
        wmod='REAL'
        dva=vel1-vel0
        call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva, vshoc,
     &   dist(step+1), dr, timlps(step+1), ftime, 1.d0, tphot)
c
        pfx='slec'
        sfx='sh4'
        fpf=' '
        call newfile (pfx, 4, sfx, 3, fn)
        fpf=fn(1:12)
        write (*,*) 'fpf''',fa,''''
        open (lupf,file=fpf,status='NEW')
        spmod='ABS'
        linemod='LAMB'
        call speclocal (lupf, tloss, eloss, egain, dlos, t, de, dh,
     &   pop(1,1), dist(step+1), dr, linemod, spmod)
        close (lupf)
c
c     reset mode and tphot
c
        specmode='DW'
        rad=dist(step+1)
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
      endif
c
      open (lupb,file=fpb,status='OLD',access='APPEND')
c
c     total energy in Inu
c
c
      band1=0.d0
      band2=0.d0
      band3=0.d0
      band4=0.d0
      band5=0.d0
      band6=0.d0
      bandall=0.d0
c
      do i=1,infph-1
        wid=photev(i+1)-photev(i)
        bandall=bandall+tphot(i)*wid*evplk
        if ((photev(i).ge.5.0d2).and.(photev(i).lt.2.0d3)) then
          band1=band1+tphot(i)*wid*evplk
        endif
        if ((photev(i).ge.2.0d3).and.(photev(i).lt.8.0d3)) then
          band2=band2+tphot(i)*wid*evplk
        endif
        if ((photev(i).ge.5.0d2).and.(photev(i).lt.8.0d3)) then
          band3=band3+tphot(i)*wid*evplk
        endif
        if ((photev(i).ge.3.0d3).and.(photev(i).lt.10.0d3)) then
          band4=band4+tphot(i)*wid*evplk
        endif
        if ((photev(i).ge.1.0d3).and.(photev(i).lt.2.0d3)) then
          band5=band5+tphot(i)*wid*evplk
        endif
        if ((photev(i).ge.2.0d3).and.(photev(i).lt.3.0d3)) then
          band6=band6+tphot(i)*wid*evplk
        endif
      enddo
c
c     times 4 pi to cover whole sky
c
      band1=band1*fpi/dr
      band2=band2*fpi/dr
      band3=band3*fpi/dr
      band4=band4*fpi/dr
      band5=band5*fpi/dr
      band6=band6*fpi/dr
      bandall=bandall*fpi/dr
c
      if (band1.lt.epsilon) band1=0.d0
      if (band2.lt.epsilon) band2=0.d0
      if (band3.lt.epsilon) band3=0.d0
      if (band4.lt.epsilon) band4=0.d0
      if (band5.lt.epsilon) band5=0.d0
      if (band6.lt.epsilon) band6=0.d0
      if (bandall.lt.epsilon) bandall=0.d0
c
  750 format(17(1pg12.5,a1))
      write (lupb,750) t,tab,de,tab,dh,tab,en,tab,pop(1,1),tab,pop(2,1),
     &tab,mu,tab,tloss/(de*dh),tab,tloss,tab,fflos,tab,band1,tab,band2,
     &tab,band3,tab,band4,tab,band5,tab,band6,tab,bandall
c
      close (lupb)
c
c     get mean ionisation state for step
c
      rdis=(dist(step)+dist(step+1))*0.5d0
      cmpf=rho1/rho0
      fi=1.0d0
c
c     accumulate spectrum
c
      tdw=te0
      tup=te1
      drdw=dr
      dvdw=vel0-vel1
      drup=dist(step+1)
      dvup=dsqrt(vpo*vel0)
      frdw=0.5d0
c
      call localem (t, de, dh)
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
        call newdif2 (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     &   frdw, jtrans)
c
      endif
c
      imod='ALL'
      call sumdata (t, de, dh, fi, dr, dr, rdis, imod)
c
c     record line ratios
c
      if (ox3.ne.0) then
        hoiii(step+1)=(fluxm(7,ox3)+fluxm(10,ox3))/(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then
        hoii(step+1)=(fluxm(1,ox2)+fluxm(2,ox2))/(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then
        hnii(step+1)=(fluxm(7,ni2)+fluxm(10,ni2))/(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then
        hsii(step+1)=(fluxm(1,su2)+fluxm(1,su2))/(fluxh(2)+epsilon)
      endif
      if (ox1.ne.0) then
        if (ispo.eq.' OI') hsii(step+1)=fluxm(3,ox1)/(fluxh(1)+epsilon)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set up initial conditions for the next zone
c     (= end conditions in previous zone)
c     Calculate new cooling times and next step
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop1, pop0)
c
      t=te1
      tt0=te1
      l0=l11
      l12=0.d0
      l11=0.d0
      r0=rho1
      p0=pr1
      v0=vel1
      dh0=dh1
      de0=de1
      bm0=bm1
      bmag=bm0
      count=0
c
      call cool (t, de0, dh0)
c
      el=dmax1(0.d0,eloss-egain)
c
c      write(*,*) t,de0,dh0,tloss
c
c     now compute free enthalpy per particle
c
      en=zen*dh0
c
      press=(en+de0)*rkb*t
      ue=(1.5d0)*(en+de0)*rkb*t
c
      cltime=dabs(ue/tloss)
c
      wmol=r0/(en+de0)
      mu=fmua(de0,dh0)
c      mu = en*wmol/((en+de0)*(mp+me))
c
      hfree=bmag*bmag/(4.d0*pi*en)+0.5d0*(wmol*ve*ve)+2.5d0*(rkb*t)/mu
c
c     calculate timescales
c
      htime=dabs(hfree/tloss)
c
      abstime=1.d0/epsilon
c
      ctime=fcolltim(de)
      rtime=frectim2(de)
      ptime=fphotim()
c
      eqtime=(1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime=dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime=(1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime=dabs(1.d0/stime)
c
      ftime=1.d0/((1.d0/(stime*atimefrac))+(1.d0/(cltime*cltimefrac)))
c     &             +(1.d0/(abstime*abtimefrac)))
c
      if ((tmod.eq.'M').and.(ftime.gt.utime)) ftime=utime
      if ((tmod.eq.'N').and.(ftime.lt.utime)) ftime=utime
      if (tmod.eq.'S') ftime=utime
c
      if (vmod.eq.'FULL') then
        write (*,620) ue,tloss,tloss,htime,cltime,ctime,rtime,ptime,
     &   stime,eqtime,abstime,ftime
      endif
c
c      if (ratmod.eq.'Y') then
c        open (lurt,file=fr,status='OLD',access='APPEND')
c        write (lurt,610) ue,tloss,tloss,htime,cltime,ctime,rtime,ptime,
c     &   stime,eqtime,abstime,ftime
c        close (lurt)
c      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Program endings
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ionisation, finish when mean ionisation < 1%
c
      if (jend.eq.'A') then
c
        totn=0.d0
c
        totn=zen
c
        ionbar=0.d0
c
        do i=1,atypes
c
          do j=2,maxion(i)
            ionbar=ionbar+(pop(j,i)*zion(i))
          enddo
c
        enddo
c
        ionbar=ionbar/totn
c
        if (ionbar.lt.0.01d0) goto 760
c
      endif
c
c     Specific species ionisation limit
c
c
      if (jend.eq.'B') then
        if (pop(jpoen,ielen).lt.fren) goto 760
      endif
c
c
c     Temperature limit
c
      if ((jend.eq.'C').and.(t.lt.tend)) goto 760
c
c
c     Distance Limit
c
      if ((jend.eq.'D').and.(dist(step+1).ge.diend)) goto 760
c
c
c     Time Limit
c
      if ((jend.eq.'E').and.(timlps(step+1).ge.timend)) goto 760
c
c
c     thermal balance dlos<1e-2
c
      if ((jend.eq.'F').and.(dlos.lt.1.d-2)) goto 760
c
c
c     cooling function test, finish when tloss goes -ve
c
      if ((jend.eq.'G').and.(tloss.lt.0.d0)) goto 760
c
c     poll for terminate file
c
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 760
c
c     otherwise go to normal iteration loop, if interlocks
c     permit
c
      if ((t.gt.100.d0).and.(step.lt.(mxnsteps-1))) then
        goto 630
      endif
c
  760 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     End model
c     write out spectrum etc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jspec.eq.'YES') then
c
c
c     Final Downstream Field
c
c
        caller='S4'
        pfx='spdw'
        np=4
        dva=0
c
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva, vshoc,
     &   dist(step+1), dr, timlps(step+1), ftime, 0.5d0, tphot)
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva, vshoc,
     &   dist(step+1), dr, timlps(step+1), ftime, 0.5d0, tphot)
c
        wmod='NFNU'
        call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva, vshoc,
     &   dist(step+1), dr, timlps(step+1), ftime, 0.5d0, tphot)
c
      endif
c
c     Upstream photon field
c
      specmode='UP'
      rad=dist(step+1)
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
      caller='S4'
      pfx=mypfx
      np=5
      dva=vel1-vel0
c
      wmod='REAL'
      call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &vshoc, dist(step+1), dr, timlps(step+1), ftime, 0.5d0,
     &tphot)
c
      wmod='NFNU'
      call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &vshoc, dist(step+1), dr, timlps(step+1), ftime, 0.5d0,
     &tphot)
c
      wmod='LFLM'
      call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &vshoc, dist(step+1), dr, timlps(step+1), ftime, 0.5d0,
     &tphot)
c
c     dynamics
c
      open (luop,file=fl,status='OLD',access='APPEND')
c
  770 format(//,' Model ended.',a1,'Distance:',1pg14.7,a1,
     &'Time:',1pg14.7,a1,'Temp:',1pg12.5//)
      write (luop,770) tab,dist(step+1),tab,timlps(step+1),tab,t
c
      call avrdata
      call wrsppop (luop)
      close (luop)
c
      open (lusp,file=fsp,status='OLD',access='APPEND')
      spmod='REL'
      linemod='LAMB'
      call spec2 (lusp, linemod, spmod)
      close (lusp)
c
c     restore photonmode
c
      photonmode=1
c
      return
c
      end
