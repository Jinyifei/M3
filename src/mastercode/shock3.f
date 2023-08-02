cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c    TREATES THE ELECTRON GAS AND THE IONS AS A TWO-PHASE MEDIUM
c    INCLUDES DIFFUSE FIELD , USES CASE B
c
c      *COMPUTATIONS ARE PERFORMED IN SUBROUTINE COMPSH3
c
c    CALL SUBR.:VELSHOCK,COMPSH3,POPCHA,PHOTSOU
c
      subroutine shock3 ()
c
      include 'cblocks.inc'
c
c
c           Variables
c
      real*8 dhpr,epotmi,hmag
      real*8 qdl,tepo
      real*8 tepr,tmag,vshockm,xhpr
c
      integer*4 i,j,lucty,nstma2
c
      character carac*46,model*32
      character ilgg*4, caract*4 ,ibell*4
      character banfil*32
c
c
      lucty=6
c
      nstma2=mxnsteps-50
c
      jcon='YES'
      jspot='NO'
      model='Shock 3'
      epotmi=iphe
c
      ibell=char(7)
      write (*,10)
   10 format(///' SHOCK MODEL SELECTED : S3'/
     &' Diffuse photon field included'/
     &' Tion and Telec maintained separately')
c
c
c    ***SET UP INITIAL IONISATION CONDITIONS
c
      model='proto-ionisation'
      call popcha (model)
      model='Shock 3'
c
      call copypop (pop, propop)
c
   20 write (*,30)
   30 format(//' Do you want the pre-ionisation to be the same'/
     &' as the proto-ionisation (Y/N) ? ',$)
      read (*,60,err=20) ilgg
c
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
c
      if ((ilgg.ne.'Y').and.(ilgg.ne.'N')) goto 20
c
      model='pre-ionisation'
      if (ilgg.eq.'N') call popcha (model)
      model='Shock 3'
c
      carac(1:33)='   '
      j=1
      do 50 i=3,atypes
        j=1+index(carac(1:33),'   ')
        caract='   '
   40 format(a2)
        if ((ipote(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) then
          write (caract,40) elem(i)
        endif
        carac(j:j+2)=caract//'   '
   50 continue
      i=j+2
c
   60 format(a2)
   70 format(a)
c
   80 write (*,90) carac(1:i)
   90 format(//' Do you allow the elements : ',a/
     &' to recombine freely to neutral state (Y/N) ? ',$)
      read (*,70,err=80) ilgg
      ilgg=ilgg(1:1)
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
      if (ilgg.eq.'Y') goto 110
      if (ilgg.ne.'N') goto 80
      do 100 i=3,atypes
        arad(2,i)=dabs(arad(2,i))
        if ((ipote(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) arad(2,i)=-
     &   arad(2,i)
  100 continue
  110 continue
c
c
      write (*,120)
  120 format(//' select shock conditions :')
  130 write (*,140)
  140 format(//' Give initial hydrogen density and ',
     &'ionised fraction : ',$ )
      read (*,*,err=130) dhpr,xhpr
      if (((xhpr.lt.0.0d0).or.(dhpr.le.0.0d0)).or.(xhpr.gt.1.d0)) goto
     &130
c
c
  150 write (*,160)
  160 format(//' Give initial and postshock temperatures : ',$)
      read (*,*,err=150) tepr,tepo
      if (tepr.le.10) tepr=10**tepr
      if (tepo.le.10) tepo=10**tepo
c
      if ((tepo.le.tepr).or.(tepr.le.100.0d0)) then
        write (*,170) ibell
  170 format(' UNSUITABLE ANSWER(S) ****************************',a)
        goto 150
      endif
c
c
  180 write (*,190)
  190 format(//' Give magnetic field (micro-Gauss) : ',$)
      read (*,*,err=180) tmag
      if (tmag.lt.0.0d0) goto 180
      hmag=tmag*1.0d-6
c
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
c
      call velshock (dhpr, xhpr, tepr, tepo, hmag)
      vshockm=vshoc/1.d5
      write (*,200) vshockm
  200 format(//' Shock velocity :',f8.1)
c
c
c    ***PHOTOIONIZATION SOURCES
c
      iso='N'
  210 write (*,220)
  220 format(//' Any photon source at the shock front',
     &'(Y/N) ? ',$)
      read (*,60,err=210) ilgg
      ilgg=ilgg(1:1)
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
c
      if ((ilgg.ne.'Y').and.(ilgg.ne.'N')) goto 210
      if (ilgg.eq.'N') goto 250
c
      call photsou (model)
c
  230 write (*,240)
  240 format(/' Dilution factor (<=0.5) :',$ )
      read (*,*,err=230) diluf
      if ((diluf.lt.0.0).or.(diluf.gt.0.5)) goto 230
c
c
  250 continue
c
c
  260 write (*,270) nstma2
  270 format(//' Give number of steps (MAX:',i4,') '/
     &' Choose two elements to be followed'/
     &' and give number of global iterations : ',$)
      read (*,*,err=260) loop,iel1,iel2,iter
c
      if (iel1.gt.28) goto 260
      if (iel2.gt.28) goto 260
c
c
      if (iel1.le.0) goto 260
      if (iel2.le.0) goto 260
c
      iel1=zmap(iel1)
      iel2=zmap(iel2)
c
      if (iel1.eq.0) goto 260
      if (iel2.eq.0) goto 260
c
      if ((loop.gt.nstma2).or.(loop.lt.8)) goto 260
      if (iter.le.0) goto 260
c
c    ***SET TYPE OF ENDING FOR THE PROGRAM
c
      xhmin=0.01d0
      write (*,280) vshockm,qdl
  280 format(//' SHOCK VELOCITY :',f8.1,' KM/SEC',7x,'Q FINE-TUNING :'
     &,f7.2)
  290 write (*,300) xhmin,elem(iel1)
  300 format(//' CHOOSE TYPE OF ENDING FOR THE PROGRAM :'/t6,
     &' A-  NORMAL ENDING  (FINAL H+/H :',2pf4.1,'% )'/t6,
     &' B-  ACCORDING TO FINAL TEMPERATURE'/t6,
     &' C-  ACCORDING TO ELAPSED TIME'/t6,
     &' D-  ACCORDING TO THE IONISED FRACTION OF : ',a2/t6,
     &' E-  ACCORDING TO LINE INTENSITIES'/t6,
     &' F-  ACCORDING TO FINAL AVERAGE CHARGE'/t6,
     &' R-  RE-INITIALIZE'/' ',t55,':: ',$)
      read (*,310,err=290) jfin
  310 format(a)
      jfin=jfin(1:1)
c
      if (jfin.eq.'a') jfin='A'
      if (jfin.eq.'b') jfin='B'
      if (jfin.eq.'c') jfin='C'
      if (jfin.eq.'d') jfin='D'
      if (jfin.eq.'e') jfin='E'
      if (jfin.eq.'f') jfin='F'
      if (jfin.eq.'r') jfin='R'
c
      if ((jfin.lt.'A').or.((jfin.gt.'F').and.(jfin.ne.'R'))) goto 290
      if (jfin.eq.'R') return
      if (jfin.eq.'A') goto 370
      if (jfin.ne.'D') goto 330
c
      write (*,320)
  320 format(/' WHICH IONISING STAGE :',$ )
      read (*,*,err=290) istf
      if ((istf.lt.2).or.(istf.gt.maxion(iel1))) goto 290
c
  330 continue
c
      if (jfin.ne.'E') goto 350
      write (*,340)
  340 format(/' CHOOSE LINE NUMBER : ' ,
     &' [OIII]   2 [OII]   3 [NII]   4 [SII]  :: ',$)
      read (*,*,err=290) istf
      if ((istf.lt.1).or.(istf.gt.4)) goto 290
  350 continue
c
      write (*,360)
  360 format(/' GIVE FINAL VALUE : ',$)
      read (*,*,err=290) fval
      if (fval.lt.0.0d0) goto 290
      if (jfin.eq.'B') texi=fval
      if (jfin.eq.'F') xhmin=fval
  370 continue
c
      jpre='C'
      if (iter.le.1) goto 400
  380 write (*,390)
  390 format(/' DETERMINATION OF THE SUCCESSIVE ',
     &'PREIONISATION CONDITIONS (U/F/C) :'/t6,
     &' U-  USING UPSTREAM PHOTON FIELD'/t6,
     &' F-  USING CONDITIONS AT THE END OF SHOCK'/t6,
     &' C-  CONSTANT FROM ONE ITERATION TO THE OTHER'/t6,
     &' T-  THERMAL EQUILIBRIUM AT PRESHOCK TEMP.'/,t55,':: ',$)
      read (*,70,err=380) jpre
      jpre=jpre(1:1)
      if (jpre.eq.'u') jpre='U'
      if (jpre.eq.'f') jpre='F'
      if (jpre.eq.'c') jpre='C'
      if (jpre.eq.'t') jpre='T'
      if ((((jpre.ne.'U').and.(jpre.ne.'F')).and.(jpre.ne.'C'))
     &.and.(jpre.ne.'T')) goto 380
c
  400 continue
c
c     get runname
c
  410 format (a80)
  420 format(/'Give a name/code for this run: ',$)
      write (*,420)
      read (*,410) runname
c
c
c     remains of old vax batch system (not used)
c
c
      banfil='Interactive'
c
      call compsh3 (lucty, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADYFLOWPLANE PARALLEL RADIATING SHOCK
c    TREATES THE ELECTRON GAS AND THE IONS AS A TWO-PHASE MEDIUM
c    INCLUDES DIFFUSE FIELD , USES CASE B
c    CALL SUBROUTINES  ROOTS,COOL,TIMION,INTVEC,TAUDIST,EQUION,
c                      WRSPPOP,SUMDATA,AVRDATA,VELSHOCK,ZETAEFF,
c                      LOCALEM,NEWDIF,TOTPHOT,ZERBUF,PREION
c              AVERINTO,copypop
c
      subroutine compsh3 (luter, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      include 'cblocks.inc'
c
      character imod*4, jmod*4, lmod*4, mmod*4
c
      real*8 teel(mxnsteps), teion(mxnsteps)
      real*8 u, usq, qf, qs, dq, upr, gam, qfmo, qsmo
      real*8 mu2,eac
c
      real*8 adj,ag0,agdiv,agmu,agpro,agq,canor,cavx
      real*8 cfm,chasl,chate,cofcrit,convg
      real*8 d,de,de2,dedhma,dedhmi,def,dexi,dfa0,dfaa
      real*8 dfv1,dh,dh2,dhpr,distav,dl0,dlf,dli
      real*8 dll0,dlla,dlli,dlosa,dlosav,dqro
      real*8 dr,drdw,drta,drup,dstep
      real*8 dte,dtee,dtef,dtel0,dtema,dtemi,dteto
      real*8 dv,dva,dvdw,dvi,dvup
      real*8 egairr,el0,elf,eli,elosrr
      real*8 fa,fel0,fel1,felio,felnf,ff,ff0,fgae
      real*8 fi,floe,fmi,fmi0,fmui,fnorm,frdv,frdv0
      real*8 frdw,frta,frto,ftim,ga0,gaf,gai,hmag,humag
      real*8 poppro,prec1,prec2,prec3,pros,qdift,qdl,rad
      real*8 rdis,rfvl,rho,rla,rsl,t,taa,tdw,tee
      real*8 tefin,tei,teio,tel0,tel1,tel2,tela
      real*8 telf,telim,tepo,tepr,tioa,tiof,tl0,tlf,tli
      real*8 tlosa0,tlosrr,tmag,tmp107,tsa,tscal,tsfin
      real*8 tsmax,tsrsl,tst,tstep,tstep0,tstepi,tsterr
      real*8 tta,tup,vcons,velav,vshockm,vsi
      real*8 wdt,web,wec,wef,wei,xhav,xhpr,xhyf,zia,zsq
c
      integer*4 i,i107,ic107,icpl,iflag,iii,irsl
      integer*4 j,kk,l,m,mi,mma,ndi
      integer*4 luter,luop,lupf,lusl,np
      integer*4 nstma2
c
      character filna*32, filnb*32, filnc*32,banfil*32, fnam*32,ispo*4
      character caller*4, wmod*4, pfx*16
c
c     External Functions
c
      real*8 fcrit,feldens,felneur,fioncha,fmua,frho
c
c
      nstma2=mxnsteps-50
c
      luter=6
      mma=300
      jcon='YES'
      jspot='NO'
      ispo='SII'
c
      dfaa=0.d0
      dfa0=0.d0
      dfv1=0.d0
      dl0=0.d0
      dlf=0.d0
      dli=0.d0
      dll0=0.d0
      dte=0.d0
      el0=0.d0
      ga0=0.d0
      prec1=0.d0
      prec2=0.d0
      prec3=0.d0
      rfvl=0.d0
      tefin=0.d0
      tli=0.d0
      tsa=0.d0
      tscal=0.d0
      tstep0=0.d0
      tstepi=0.d0
c
      dedhmi=0.0d0
      dedhma=0.0d0
      do j=1,atypes
        dedhma=dedhma+zion(j)
        if (arad(2,j).le.0.d0) dedhmi=dedhmi+zion(j)
      enddo
c
      luop=21
      lupf=27
      lusl=45
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
c
      call s3fheader (fnam, filna, filnb, filnc, luop, lupf, lusl, dhpr,
     & xhpr, tepr, tepo, tmag, banfil)
c
c    ***SHOCK INITIALISATION FINISHED
c    COMPUTATION STARTS HERE ; ITERATION INDEX : MI
c
      do 410 mi=1,iter
        rad=1.d38
        wdil=diluf
        fi=1.0d0
c
        qdift=0.0d0
c
        if (mi.eq.1) goto 70
        dva=0.0d0
        tei=tepr
        dr=1.d0
c
c    ***SETTING UP OF PREIONISATION CONDITIONS
c
c
        if (jfin.eq.'A') tsmax=1.d2*tsmax
        open (luop,file=filna,status='OLD',access='APPEND')
        write (luop,10) jpre
c
   10 format(//' SETTING UP OF THE NEXT PREIONISATION CONDITIONS :'
     &,a)
        lmod='UP'
        call localem (tei, depr, dhpr)
c      call contcont(tei, depr, dhpr)
c
c
        call totphot (tei, dhpr, fi, rad, dr, dva, wdil, lmod)
        if (jpre.ne.'U') goto 30
        do 20 l=1,4
c
          vsi=vshoc
          call preion (luter, luop, tsmax, vshoc, fi, dhpr, depr, tepr,
     &     qdift, d)
          xhpr=pop(2,1)
c
          call velshock (dhpr, xhpr, tepr, tepo, hmag)
          if ((dabs(vsi-vshoc)/dmax1(vsi,vshoc)).lt.0.05) goto 30
   20   continue
   30   continue
c
        if (jpre.ne.'F') goto 40
        call copypop (popf, pop)
        xhpr=pop(2,1)
        tepr=tefin
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
c
   40   continue
c
        if (jpre.ne.'T') goto 50
        call equion (tepr, depr, dhpr)
        xhpr=pop(2,1)
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
c
   50   continue
        if (jpre.ne.'C') goto 60
c
        call copypop (popini, pop)
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
   60   continue
c
c
        close (luop)
c
        call zerbuf
   70   continue
c
c
        call copypop (pop, popi)
        call copypop (pop, popf)
        call copypop (pop, popini)
c
        xh(1)=xhpr
        xh(2)=xh(1)
        veloc(1)=vshoc
        veloc(2)=vpo
        vshockm=vshoc/1.d5
        te(1)=tepr
        te(2)=tepo
        teel(1)=te(1)
        teel(2)=teel(1)
        teion(1)=te(1)
        teion(2)=te(2)+(felneur(pop,1)*(te(2)-teel(2)))
        dhy(1)=dhpr
        dhy(2)=dhpo
        dist(1)=0.0d0
        dist(2)=0.0d0
        deel(1)=depr
        deel(2)=depo
c
        timlps(1)=0.0d0
        qdl=1.5+((0.2*(tepo/1.d5))*(nstma2/(loop)))
        fnorm=qdl/(1.0-dexp(-qdl))
        vcons=-dabs((vpo-vfin)/(loop-4.0))
        upr=vpo
        u=upr
        usq=u*u
        humag=hmago*(vpo/u)
        gam=(((pres*u)/fm)-usq)-humag
        qf=(2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+usq)
        qfmo=-(((5.d0*gam)+(4.d0*humag))+usq)
c
        open (luop,file=filna,status='OLD',access='APPEND')
        write (luop,80) loop,mi,tepr,qdift
   80 format(//i4,' STEPS',t40,' ITERATION #',i2,10x,'TEPR :'
     &     ,f8.0,10x,'PRESHOCK QTOT : ',1pg10.3/t40,' ==========='/)
c     write(*, 2027) mi, tepr, qdift, fnam
        write (luop,90) jpre
   90 format(' ',t9,'INITIAL STATE OF IONISATION ' ,
     &     'FOR ALL ELEMENTS  (JPRE : ' ,a1,' )  :'/)
c
        call wionpop (luop, popi)
c
        write (*,100)
  100 format(/' STEP NUMBERS AS THEY ARE COMPLETED :'/' L#',t8,'<Te>'
     &,t18,'<TLOSS>',t27,'<DLFR>',t36,'TSTEP',t45,'<DENS.>',t55
     &,'<EL.DENS>',t69,'DV')
c
        write (luop,110) elem(1),rom(1),elem(1),rom(2),(elem(iel1),
     &   rom(j),j=1,6),(elem(iel2),rom(j),j=1,6)
  110 format(//' STEP NUMBERS AS THEY ARE COMPLETED :'//' L#',t7,
     &'<LOSS.R.>',t18,'TSTEP',t26,14(1x,a2,a3,1x)/)
        close (luop)
c
c
c    DOMWNSTREAM COMPUTATION ; STEP INDEX : M
c
        drta=0.0d0
        frta=0.03d0
        eac=100.d0
        canor=3.5d0
        frdv=1.d-2
        ff=1.0d0
        frto=0.d0
        fmi=fnorm*dexp(-((qdl*frto)/loop))
        convg=0.004d0
        do 120 i=1,5
  120     cav(i)=canor
        cavx=canor
c
c
        m=1
  130   m=m+1
        iflag=1
        icpl=0
        eac=dsqrt(50.d0*dmax1(50.d0,eac))
        irsl=0
c
c    ***IT OVERWRITES STORED QUANTITIES (FROM 100 UP) WHEN NUMBER
c    OF STEPS IS BIGGER THEN 250
c
        if (m.eq.mma) then
          m=100
          do 140 i=0,1
            teel(m-i)=teel(mma-i)
            te(m-i)=te(mma-i)
            teion(m-i)=teion(mma-i)
            dhy(m-i)=dhy(mma-i)
            xh(m-i)=xh(mma-i)
            deel(m-i)=deel(mma-i)
            dist(m-i)=dist(mma-i)
            veloc(m-i)=veloc(mma-i)
            timlps(m-i)=timlps(mma-i)
            hoiii(m-i)=hoiii(mma-i)
            hoii(m-i)=hoii(mma-i)
            hnii(m-i)=hnii(mma-i)
            hsii(m-i)=hsii(mma-i)
  140     continue
        endif
        qs=qf
c
        qsmo=qfmo
        dtema=0.25d0
        dtemi=0.033d0
        if ((m.lt.(loop/3)).and.((((deel(m)/dhy(m))-dedhmi)/dedhma)
     &   .lt.0.2d0)) dtemi=dtemi/dabs(dlog(dmax1(((deel(m)/dhy(m))-
     &   dedhmi)/dedhma,1.d-10)))
        fmi0=fmi
        fmi=fnorm*dexp(-((qdl*frto)/loop))
        frdv=(frdv*(fmi0/fmi))*(((cavx+0.05)/canor)**(-0.5))
        if (m.gt.2) then
          if (((te(m).lt.200.0).and.(dexi.lt.7.d-2)).or.(dexi.lt.2.d-2))
     &      dtema=0.36
          if (wef.gt.0.75) ff=dmin1(1.d0,ff**dmin1(1.d0,0.5d0+((dteto/
     &     dtema)*(5.d0**((1.d0/wef)**3.d0)))))
          if (dteto.lt.(1.05*dtemi)) then
            fmui=0.97*dsqrt(dtemi/(dteto+1.d-9))
            frdv=frdv*dmin1(2.d0,fmui)
          else if ((dte.lt.(2.0*dtemi)).and.(dteto.gt.(2.5*dtemi)))
     &     then
            agq=1.60
            frdv=frdv*dmin1(agq-((agq-1.d0)*((dmin1(dteto,dtema)/dtema)*
     &       *0.5d0)),dsqrt((2.d0*dtemi)/dte))
          endif
        endif
        frdv=dmin1(frdv,1.d0)
        dvi=fmi*vcons
        if (dabs(dvi/veloc(m)).gt.(0.75d0*dtema)) dvi=-((veloc(m)*
     &   0.75d0)*dtema)
c
c    ***DERIVES GUESSED VALUES FOR : DE,DH,TLOSS,DR,TSTEP AND T
c
  150   ic107=iflag
c
c      *STEP INDEX M=2
        icpl=icpl+1
        if (m.eq.2) then
          xh(2)=popi(2,1)
          dh=dhy(2)
          de1=feldens(dh,popi)
          de=de1
          deel(2)=de
          t=te(2)
          tela=teel(2)
          tstep=0.0d0
          dr=1.d5
c
c
          lmod='DW'
          call localem (tela, de, dh)
c      call contcont(tela, de, dh)
          call totphot (tela, dh, fi, rad, dr, dva, wdil, lmod)
c
c
          call cool (tela, de, dh)
c
c
          el0=eloss
          ga0=egain
          dl0=dlos
          tl0=tloss
c
c      **STEP INDEX M=3
c
        else if (m.eq.3) then
          dfaa=dabs((fmi*frdv)/dfv1)
          tloss=tl0
          dh=dhy(3)
          de=deel(3)
          t=te(3)+(((0.3*(te(3)-te(2)))*te(3))/te(2))
          tela=teel(3)+(((0.3*(teel(3)-teel(2)))*teel(3))/teel(2))
          tela=dmax1((1.0-(dtema/1.9))*teel(m),dmin1((1.0+((dtema*(1.0+
     &     dtema))/1.9))*teel(m),tela))
          dlosa=dl0
c
c      *STEP INDEX M>3
c
        else
          if (dabs(dlf).gt.0.7d0) then
            dfaa=dabs((fmi*frdv)/dfv1)
            adj=dfaa/dfa0
            tsterr=(tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36))))*
     &       adj
c
c
          else if ((dabs(dlf).gt.0.3d0).and.(dabs(dl0).gt.0.3d0))
     &     then
            agmu=5.d0*(1.d0+((3.d0*dabs(dl0))**2.d0))
            dfaa=dabs((fmi*frdv)/dfv1)
            dfaa=dmin1(agmu,dmax1(1.d0/agmu,dfaa))
            adj=dfaa/dfa0
            adj=dmin1(agmu,dmax1(1.d0/agmu,adj))
            tsterr=(tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36))))*
     &       adj
            tsterr=dmin1(agmu*tstep0,dmax1(tstep0/agmu,tsterr))
          else
            adj=1.d0
            dfaa=1.d0
            agmu=5.d0
            tsterr=(tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36))))*
     &       adj
            tsterr=dmin1(agmu*tstep0,dmax1(tstep0/agmu,tsterr))
          endif
c
c
          if (irsl.eq.0) then
            dlosa=dl0+((0.5d0*(dl0-dli))*dfaa)
            if (dabs(dlosa).gt.1.d0) dlosa=dsign(1.d0,dlosa)
            if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.(ff*tta))) tsterr=
     &       ff*tta
            if (iflag.gt.1) tsterr=dmin1(tsterr,1.10*tsa)
            web=fcrit(dlosa,eac)
            web=fcrit(dlosa/(1.d0+(1.d1*web)),eac)
            tstep=dexp(((1.d0-web)*dlog(tsterr+1.d-30))+(web*((0.3d0*
     &       dlog(ff*tta))+(0.7d0*dlog(tsterr+1.d-30)))))
          else
            dlosa=dl0+(0.25d0*((dlf-dl0)+((dl0-dli)*dfaa)))
            if (dabs(dlosa).gt.1.d0) dlosa=dsign(1.d0,dlosa)
            if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.(ff*tta))) tsterr=
     &       ff*tta
            web=fcrit(dlosa,eac)
            web=fcrit(dlosa/(1.d0+(1.d1*web)),eac)
            tstep=dexp(((1.d0-web)*dlog(dmin1(tsterr,tsa)+1.d-30))+(web*
     &       dlog(tsa+1.d-30)))
          endif
c
c
          tloss=tl0+((0.5d0*(tl0-tli))*dfaa)
          tlf=tloss+(tloss-tl0)
c
c
          if ((tli/tl0).gt.0.d0) then
            dh=dhy(m)+(((0.5d0*(dhy(m)-dhy(m-1)))*(dhy(m)/dhy(m-1)))*
     &       dfaa)
            t=te(m)+(((0.5d0*(te(m)-te(m-1)))*(te(m)/te(m-1)))*dfaa)
            tel1=teel(m)+(((0.5*(teel(m)-teel(m-1)))*(teel(m)/teel(m-1))
     &       )*dfaa)
          else
            dh=dhy(m)-((0.5*(dhy(m)-dhy(m-1)))*dfaa)
            t=te(m)-((0.5*(te(m)-te(m-1)))*dfaa)
            tel1=teel(m)-((0.5*(teel(m)-teel(m-1)))*dfaa)
          endif
c
c
          t=dmax1((1.d0-(dtema/2.d0))*te(m),dmin1((1.d0+((dtema*(1.d0+
     &     dtema))/2.d0))*te(m),t))
          dlla=dsign(dmin1(dabs(t-teion(m)),dmin1(2.d0,dfaa)*dabs(dll0))
     &     ,dll0)
          tioa=teion(m)+(0.5d0*dlla)
          tel2=t+((t-tioa)/(felneur(popi,1)+1.d-10))
c
c
          if (((tstep.gt.(12.d0*tscal)).or.((dabs(te(m)-teel(m))/te(m))
     &     .lt.1.d-2)).or.((tstep.gt.(2.5d0*tscal)).and.((dabs(t-tel2)/
     &     t).lt.0.05d0))) then
            tela=tel1
          else
            wei=dmin1(0.6d0,(0.2d0*tstep)/tscal)
            tela=(wei*tel1)+((1.0-wei)*tel2)
          endif
c
c
          tela=dmax1((1.d0-(dtema/1.9d0))*teel(m),dmin1((1.d0+((dtema*
     &     (1.d0+dtema))/1.9d0))*teel(m),tela))
          de=(deel(m)*dh)/dhy(m)
        endif
c
c
c    ***DERIVES AVERAGES FOR DH,DE,MU,RHO AND TE,
c    COMPUTES CONDITIONS AT THE END OF TIMESTEP
c
  160   frdv0=frdv
        ff0=ff
        dv=dvi*dsign(frdv,tl0)
        dva=dabs(dv)
        u=upr+(dv)
        usq=u*u
        humag=hmago*(vpo/u)
        gam=(((pres*u)/fm)-usq)-humag
        qf=(2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+usq)
        qfmo=-(((5.d0*gam)+(4.d0*humag))+usq)
        dq=(qfmo-qsmo)/2.d0
        veloc(m+1)=u
        velav=veloc(m)+(dv/2.d0)
        dstep=((veloc(m+1)+veloc(m))*tstep)/2.d0
        dr=0.5d0*dstep
c
c
        if (iflag.gt.ic107) goto 190
        open (lupf,file=filnb,status='OLD',access='APPEND')
        write (*,170) m,tela,tloss,dlosa,tstep,dh,de,dv
  170 format(' ',i3,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
        write (lupf,180) m,iflag,ic107,teel(m),tl0,tlf,dlosa,dv,tstep,
     &   dh,de,frdv,ff
  180 format(' ',3(i3),0pf9.0,9x,7(1pg11.3),2(1pg9.2))
        close (lupf)
  190   continue
c
c
        lmod='DW'
        wei=0.5d0
        if (iflag.eq.ic107) then
          call localem (tela, de, dh)
c        call contcont(tela,de,dh)
        endif
c
        call totphot (tela, dh, fi, rad, dr, dva, wdil, lmod)
        call copypop (popi, pop)
        call timion (tela, def, dh, xhyf, tstep)
        call copypop (pop, popf)
c
c
        call averinto (wei, popi, popf, pop)
c
c
        xh(m+1)=xhyf
        xhav=(xh(m)+xh(m+1))/2.d0
        dhy(m+1)=(dhpo*vpo)/veloc(m+1)
        dh1=dhy(m)
        dh2=dhy(m+1)
        dh=(dh1+dh2)/2.d0
        de1=feldens(dh1,popi)
        de2=feldens(dh2,popf)
        deel(m+1)=de2
        de=(de1+de2)/2.d0
        mu2=fmua(de2,dh2)
        rho=frho(de,dh)
        te(m+1)=(mu2*gam)/rgas
        te(m+1)=dmax1(0.15d0*te(m),dmin1(te(m)/0.15d0,te(m+1)))
        t=(te(m+1)+te(m))/2.d0
        t=dmax1((1.d0-(dtema/1.8d0))*te(m),dmin1((1.d0+((dtema*(1.d0+
     &   dtema))/1.8d0))*te(m),t))
c
c
        ndi=15
        fel0=felneur(popi,1)
        fel1=felneur(popf,1)
        felio=felneur(pop,2)
        zsq=fioncha(pop,2)**2
        zia=fioncha(pop,1)
        tst=tstep/ndi
        do 210 iii=1,10
          telim=teel(m)
          taa=te(m)
          tiof=teion(m)
          do 200 kk=1,ndi
            rla=((9.42d0+(1.5d0*dlog(telim)))-(0.5d0*dlog(de+1.d-26)))-
     &       dlog(zia+1.d-20)
            tscal=((3.197d-3*fmua(0.d0,dh))*(((telim*1836.d0)+(tiof/
     &       fmua(0.d0,dh)))**1.5d0))/(((((1.d0+(1.d0/felio))*rla)*de)*
     &       zsq)+1.d-20)
            dlli=(taa-tiof)-((taa-tiof)*dexp(-(tst/tscal)))
            taa=te(m)+((kk*(te(m+1)-te(m)))/ndi)
            dlli=dsign(dmin1(dabs(taa-tiof),dabs(dlli)),dlli)
            tiof=tiof+dlli
            if (tiof.lt.((1.d0-dtema)*teion(m))) tiof=dmax1(tiof,0.2d0*
     &       teion(m))
            felnf=fel0+((kk*(fel1-fel0))/ndi)
            tee=taa+((taa-tiof)/felnf)
            telim=dmax1((1.d0-((4.d0*dtema)/ndi))*telim,dmin1((1.d0+
     &       (((4.0*dtema)*(1.d0+dtema))/ndi))*telim,tee))
            telim=dmax1((1.d0-(dtema/0.9d0))*teel(m),dmin1((1.d0+
     &       ((dtema*(1.d0+dtema))/0.9d0))*teel(m),telim))
  200     continue
c
c
          teion(m+1)=tiof
          if (teion(m+1).lt.((1.d0-dtema)*teion(m))) teion(m+1)=
     &     dmax1(teion(m+1),0.2d0*teion(m))
          tel0=tela
          teel(m+1)=te(m+1)+(((te(m+1))-teion(m+1))/(felneur(popf,1)+
     &     1.d-10))
          teel(m+1)=dmax1(0.1d0*teel(m),dmin1(teel(m)/0.1d0,teel(m+1)))
          telf=teel(m+1)
          telf=dmax1((1.d0-(dtema/0.9d0))*teel(m),dmin1((1.d0+((dtema*
     &     (1.d0+dtema))/0.9d0))*teel(m),telf))
          tela=(teel(m)+teel(m+1))/2.d0
          tela=dmax1((1.d0-(dtema/1.9d0))*teel(m),dmin1((1.d0+((dtema*
     &     (1.0+dtema))/1.9d0))*teel(m),tela))
          if (iii.gt.4) tela=(tela+tel0)/2.d0
          dtel0=dabs(tela-tel0)/dmax1(tela,tel0)
          if (dtel0.lt.4.d-5) goto 220
  210   continue
  220   continue
c
c
        wdt=dmax1(1.d0-(0.4d0*dabs(dlog10(1.d-30+(teel(m)/te(m))))),
     &   dmin1(0.75d0,0.5d0+(iflag-ic107)))
        dte=(dabs(te(m+1)-te(m))/dmax1(te(m),te(m+1)))+1.d-10
        dtef=dabs(teion(m)-teion(m+1))/dmax1(teion(m),teion(m+1))
        dtee=dabs(teel(m+1)-teel(m))/dmax1(teel(m),teel(m+1))
        dteto=dmax1(0.8d0*dte,wdt*dtee,0.25d0*dtef)
c
        dqro=dq*rho
        tlosrr=tloss
        elosrr=eloss
        egairr=egain
        tsterr=tstep
c
c    ***FINDS FINAL AND AVERAGE COOLING RATE
c
        wei=0.5d0
        lmod='DW'
        mmod='DIS'
        call localem (tela, de, dh)
c      call contcont(tela,de,dh)
c
        if ((dabs(dlosa).lt.0.99d0).or.(dabs(dl0).lt.0.99d0)) then
          call totphot (tela, dh, fi, rad, dstep, dva, wdil, lmod)
        endif
        call taudist (dh, fi, frta, drta, rad, pop, mmod)
        call copypop (popf, pop)
        call cool (telf, de2, dh2)
        call averinto (wei, popi, popf, pop)
c
c
        elf=eloss
        eloss=(el0+elf)/2.d0
        gaf=egain
        egain=(ga0+gaf)/2.d0
        dlf=dlos
        dlosa=(dl0+dlf)/2.d0
        tlf=tloss
        tloss=(tl0+tlf)/2.d0
        rsl=tl0/dsign(dmax1(dabs(tlf),0.01d0*dabs(tl0),1.d-38),tlf)
        wef=fcrit(dlosa,eac)
        wef=fcrit(dlosa/(1.d0+(1.d1*wef)),eac)
c
c
        tsfin=dabs((2.d0*dqro)/(tl0+tlf))
        tsrsl=dmin1((tsfin),dabs((2.d0*dqro)/dmax1(dabs(tl0),dabs(tlf)))
     &   )
        tta=(drta/((veloc(m)+veloc(m+1))/2.d0))+1.d-20
        iflag=iflag+1
        i107=iflag-ic107
c
c    test convergence
c
        pros=dmax1(1.d0,((i107*fcrit(dlosa,1.d1))+1.d-4)**0.35d0)
        cofcrit=dmin1(1.001d0,convg+(1.5d0*fcrit(dlosa/pros,3.d1)))
        if ((((rsl.lt.0.0d0).and.(dabs(rsl).lt.0.5d0)).and.(irsl.lt.2))
     &   .and.(dteto.gt.((1.d0+((0.1d0*(irsl+1))*(irsl+1)))*dtemi)))
     &   goto 230
c
c
        if (dteto.gt.dtema) goto 230
        floe=dabs((eloss-elosrr)/dmax1(dabs(eloss),dabs(elosrr)))
        fgae=dabs((egain-egairr)/dmax1(egain,egairr,1.d-33))
        ftim=dabs((tsfin-tsterr)/dmax1(tsfin,tsterr))
c
c    ***ITERATES IF NECESSARY ON TIME STEP TO IMPROVE CONVERGENCE
c
        if ((ftim.lt.convg).or.((i107.gt.1).and.(ftim.lt.cofcrit)))
     &   goto 300
  230   continue
        tmp107=dble(i107)
        eac=dmin1(2.d3,eac*dmax1(1.0d0,2.0d0*dsqrt(tmp107/8.d0)))
        chate=(0.7d0*dtemi)/dteto
        if (rsl.gt.0.d0) then
          if (((((rsl.le.3.0d0).or.(irsl.ge.2)).or.(dteto.le.dtemi))
     &     .or.(i107.ge.2)).or.(dteto.gt.dtema)) goto 240
          chasl=0.5d0
          irsl=irsl+1
          frdv=frdv0*dmax1(chate,chasl)
  240     continue
c
c
          tsa=dexp(((1.d0-wef)*dlog(tsfin+1.d-30))+((wef*0.5d0)*
     &     (dlog(ff*tta)+dlog(dmin1(tsterr,ff*tta)+1.d-30))))
          agdiv=dmax1(10.d0,((tsterr/(tsfin+1.d-30))+0.01d0)**(0.25d0+
     &     (0.025d0*dble(i107))))
          tsa=dmax1(tsa,(fcrit(dlosav,1.d0)*tsterr)/agdiv)
          if ((i107.gt.4).and.(dabs(dlosa).gt.0.35d0)) then
            wei=dmin1(0.5d0,1.d0*((i107/1.d1)**2))*(1.d0-fcrit(dlosa,
     &       1.d0))
            tsa=((1.d0-wei)*tsa)+(wei*dmin1(1.45d0*tsa,dmax1(0.7d0*tsa,
     &       tsterr)))
          endif
c
c
        else
          if (((dabs(rsl).ge.0.5d0).or.(irsl.ge.2)).or.(dteto.le.dtemi))
     &      goto 250
          chasl=dabs(rsl)/5.d0
          irsl=irsl+1
          frdv=frdv0*dmax1(chate,chasl)
  250     continue
c
c
          fa=dabs(tsterr-tsfin)/dmax1(tsterr,tsfin)
          tsa=(tsfin+tsterr)/2.d0
          cfm=0.5d0-(0.35d0*fcrit(dlosa,2.5d1))
          if (fa.gt.cfm) tsa=tsterr*(1.0d0+dsign(cfm,dlog(1.d-38+(tsfin/
     &     (tsterr+1.d-30)))))
          tsa=dmin1(ff*tta,dmax1(tsa,tsrsl/2.0))
        endif
c
c
        tstep=tsa
        open (lupf,file=filnb,status='OLD',access='APPEND')
        write (*,260) tela,tloss,dlosa,tsterr,tsfin,tta,rsl
  260 format('    ',f10.0,1pg11.3,0pf7.3,4(1pg11.3))
        write (lupf,270) m,iflag,ic107,teel(m),teel(m+1),tl0,tlf,dlosa,
     &   dv,tsterr,tsfin,tta,tstep,rsl
  270 format(' ',3(i3),2(0pf9.0),8(1pg11.3),0pf8.3)
        close (lupf)
c
c
        ag0=1.d0
        if ((dteto.gt.dtema).and.((i107.gt.1).or.(dabs(dlosa).gt.0.3d0))
     &   ) then
          ag0=dmin1(0.6d0*(1.d0-dtema),0.7d0*((dtema/dteto)**2.d0))
          eac=dmin1(2.d3,1.5*eac)
        else if ((i107.gt.7).and.((((deel(m)/dhy(m))-dedhmi)/dedhma)
     &   .lt.0.2d0)) then
          ag0=0.125d0
        else if (i107.gt.7) then
          ag0=0.3d0
        endif
c
c
        if (ag0.lt.1.d0) then
          wec=wef
          frdv=frdv*dexp((1.0-wec)*dlog(ag0))
          ff=ff*dexp(wec*dlog(ag0))
        endif
c
c
        if (iflag.gt.30) goto 280
        if ((frdv.ne.frdv0).or.(ff.ne.ff0)) goto 150
        goto 160
c
c    ***IF CONVERGENCE UNACHIEVABLE : EXIT GRACEFULLY
c
  280   write (*,290) m
  290 format(/' <<<<<<<<<<<< NO CONVERGENCE AT STEP#',i4/)
  300   continue
c
c    ***WHEN CONVERGENCE ACHIEVED , PREPARES NEXT TIME STEP
c
        do 310 j=1,4
  310     cav(j)=cav(j+1)
        agpro=i107+dmin1(1.5d0,icpl-1.d0)
        cav(5)=(fcrit(dlosa,8.d0)*canor)+((1.d0-fcrit(dlosa,8.d0))*
     &   agpro)
        cavx=0.d0
        i=5
        do 320 j=5,6-i,-1
  320     cavx=cavx+(j*cav(j))
        cavx=cavx/((dble(i)+1.d0)*(dble(i)/2.d0))
c
c
        dist(m+1)=dist(m)+dstep
        timlps(m)=timlps(m-1)+tsterr
        tsmax=5.d0*timlps(m)
        tefin=te(m+1)
c
c    ***INTEGRATES DIFFUSE FIELD,COMPUTES AVERAGE SPECTRUM AND
c    SUMS UP RELEVANT QUANTITIES
c
        dexi=((de/dh)-dedhmi)/dedhma
        rdis=(dist(m)+dist(m+1))/2.d0
        tdw=tela
        tup=((tela*te(2))*teel(2))**(1.d0/3.d0)
        drdw=dstep
        dvdw=dva
        drup=dist(m)
        dvup=dsqrt(vpo*u)
        frdw=0.5d0
        jmod='LODW'
        imod='ALL'
c
c
        call cool (tela, de, dh)
        call zetaeff (fi, dh)
        call newdif (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     &   frdw, jmod)
c
c
c    ***OUTPUT PART OF THE DATA AND TEST ENDING
c
        call sumdata (tela, de, dh, fi, drdw, drdw, rdis, imod)
c
        if (ox3.ne.0) then
          hoiii(m)=(fluxm(7,ox3)+fluxm(10,ox3))/(fluxh(2)+epsilon)
        endif
        if (ox2.ne.0) then
          hoii(m)=(fluxm(1,ox2)+fluxm(2,ox2))/(fluxh(2)+epsilon)
        endif
        if (ni2.ne.0) then
          hnii(m)=(fluxm(7,ni2)+fluxm(10,ni2))/(fluxh(2)+epsilon)
        endif
        if (su2.ne.0) then
          hsii(m)=(fluxm(1,su2)+fluxm(2,su2))/(fluxh(2)+epsilon)
        endif
        if (ox1.ne.0) then
          if (ispo.eq.' OI') hsii(m)=fluxm(3,ox1)/(fluxh(1)+epsilon)
        endif
        if (jfin.ne.'E') goto 330
        if (istf.eq.1) rfvl=hoiii(m)
        if (istf.eq.2) rfvl=hoii(m)
        if (istf.eq.3) rfvl=hnii(m)
        if (istf.eq.4) rfvl=hsii(m)
  330   continue
c
c
        open (luop,file=filna,status='OLD',access='APPEND')
        open (lupf,file=filnb,status='OLD',access='APPEND')
        open (lusl,file=filnc,status='OLD',access='APPEND')
c
c
        write (*,260) tela,tloss,dlosa,tstep,tsfin,tta,rsl
        write (luop,340) m,tloss,tstep,pop(1,1),pop(2,1),(pop(j,iel1),j=
     &   1,6),(pop(i,iel2),i=1,6)
  340 format(' ',i3,1pg10.2,1pg9.2,1x,14(0pf7.4))
        write (lupf,270) m,iflag,ic107,teel(m),teel(m+1),tl0,tlf,dlosa,
     &   dv,tsterr,tsfin,tta,tstep,rsl
        write (lupf,350) cavx,caseab(1),dexi,fgae,floe,ftim,cofcrit,t,
     &   tela,dteto,wef
  350 format(' CAVX :',0pf6.3,3x,'CASE A,B :' ,0pf4.2,3x,'FEXIT :'
     &,1pg10.3,3x,4(0pf7.4),2(0pf10.0),2(0pf7.4))
c
        write (lusl,*) 'It#   distance         dr            Te    ',' L
     &oss rate        de            dh   '
        write (lusl,360) m,dist(m),dstep,tela,tloss,de,dh
  360 format(i3.3,4(1pg14.7))
        call wionabal (lusl, pop)
c
        close (luop)
        close (lupf)
        close (lusl)
c
        if (iflag.gt.30) goto 370
        if ((t.lt.3000.0d0).and.(dexi.lt.xhmin)) goto 370
        if ((jfin.eq.'B').and.(te(m+1).le.texi)) goto 370
        if (((jfin.eq.'C').and.(m.gt.5)).and.(fval.lt.timlps(m))) goto
     &   370
        if (((((jfin.eq.'D').and.(m.gt.5)).and.(pop(istf+1,iel1)
     &   .lt.prec2)).and.(pop(istf,iel1).lt.prec1)).and.(pop(istf,iel1)
     &   .lt.fval)) goto 370
        if ((((jfin.eq.'E').and.(m.gt.(loop/5.0))).and.(rfvl.lt.fval))
     &   .and.(prec3.gt.fval)) goto 370
        prec2=pop(istf+1,iel1)
        prec1=pop(istf,iel1)
        prec3=rfvl
c
c
        call copypop (popf, popi)
        call copypop (popf, pop)
c
c
        upr=u
        tli=tl0
        tl0=tlf
        eli=el0
        el0=elf
        gai=ga0
        ga0=gaf
        dli=dl0
        dl0=dlf
        tstepi=tstep0
        tstep0=tsterr
        if (dlosa.gt.0.5d0) tstep0=tsfin
        tlosa0=tloss
        dfv1=fmi*frdv
        dfa0=dfaa
        dll0=dlf
        frto=frto+frdv
        goto 130
c
c    ***OUTPUT QUANTITIES STORED IN ARRAYS
c
  370   continue
c
c
        if (mi.eq.iter) then
          open (luop,file=filna,status='OLD',access='APPEND')
          write (luop,380)
  380 format(//' LIST OF STORED ' ,'QUANTITIES : ' //' L*',t10,'<TE>'
     &,t18,'<TEEL>',t27,'<TEION>',t37,'<DENS.>',t46,'<EL.DENS>',t56,
     &'<DIST>',t66,'<V>',t73,'<ELPS.TM>',t83,'OIII/H-B',t92,'OII/H-B'
     &,t101,'NII/H-A',t110,'SII/H-A'/)
c
c
          do 400 l=2,m
            distav=(dist(l)+dist(l+1))/2.0
            velav=(veloc(l+1)+veloc(l))/2.0
            t=(te(l)+te(l+1))/2.0
            tela=(teel(l)+teel(l+1))/2.0
            teio=(teion(l)+teion(l+1))/2.0
            dh=(dhy(l)+dhy(l+1))/2.0
            de=(deel(l)+deel(l+1))/2.0
            write (luop,390) l,t,tela,teio,dh,de,distav,velav,timlps(l),
     &       hoiii(l),hoii(l),hnii(l),hsii(l)
  390 format(' ',i3,3(0pf10.0),2(1pg10.3),7(1pg9.2))
  400     continue
          close (luop)
c
c
        endif
c
c
        dva=dsqrt(vpo*u)
        dr=1.d0
        lmod='UP'
c
        call totphot (tepr, dhpr, fi, rad, dr, dva, wdil, lmod)
c
        caller='S3'
        pfx='psoup'
        np=5
        wmod='REAL'
c
        call wpsou (caller, pfx, np, wmod, tepr, depr, dhpr, veloc(m),
     &   vpo, vshoc, dist(m), dstep, timlps(m), tsterr, (2.d0/4.d0),
     &   tphot)
c
c***    MODEL TERMINATED
c
  410 continue
c
c
      open (luop,file=filna,status='OLD',access='APPEND')
c
c
      if (iter.le.1) goto 500
      write (luop,420)
  420 format(///' INITIAL PARAMETERS CHARACTERISING THE ',
     &'LAST ITERATION :'/)
c
c
  430 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' ====================',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write (luop,440) zgas
  440  format(' ',t9,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun)'/
     &' and the proto-ionisation used.')
c
      call wionabal (luop, poppro)
c
  450  format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
      write (luop,460)
  460  format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')
c
      write (luop,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii,
     &diluf
c
  470  format(' ',a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
c
      vshockm=vshoc/1.d5
      write (luop,480) tepo,vshockm
c
  480  format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' =====================')
c
      write (luop,430) dhpr,xhpr,tepr,tmag
      write (luop,90) jpre
c
      call wionpop (luop, popini)
c
      write (luop,450)
      write (luop,460)
      write (luop,470) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii,
     &diluf
      write (luop,490) qdift
  490 format(/' ',t10,' PRESHOCK UPSTREAM IONISING FLUX : ' ,1pg10.3/)
      write (luop,480) tepo,vshockm
c
c    ***OUTPUT SPECTRUM AND AVERAGE IONIC POP. AND TEMP.
c
c
  500 continue
c
c
      call avrdata
c
      call wrsppop (luop)
      close (luop)
c
c
      write (*,510) fnam
  510 format(/' OUTPUT CREATED IN &&&&&& FILE : ',a/)
c
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Subroutine to get the file header buisness out where
c    it can be worked on, and/or modified for other headers.
c
c    first used for photo3 and photo4
c
c    RSS 8/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine s3fheader (fnam, filna, filnb, filnc, luop, lupf, lusl,
     & dhpr, xhpr, tepr, tepo, tmag, banfil)
c
      include 'cblocks.inc'
c
c
      real*8 dhpr, xhpr, tepr, tepo,tmag, vshockm
c
      integer*4 luop,lupf,lusl
c
      character filna*32, filnb*32,filnc*32
      character banfil*32,fnam*32
      character fn*32
      character pfx*16,sfx*4
c
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
c
c    ***DERIVE FILENAME : shckn**.sh3
c
      fn=' '
      pfx='shckn'
      sfx='sh3'
      call newfile (pfx, 5, sfx, 3, fn)
      filna=fn(1:13)
      fnam=filna
c
c
      fn=' '
      pfx='shapn'
      sfx='sh3'
      call newfile (pfx, 5, sfx, 3, fn)
      filnb=fn(1:13)
c
c
      fn=' '
      pfx='allion'
      sfx='sh3'
      call newfile (pfx, 6, sfx, 3, fn)
      filnc=fn(1:14)
c
c
      open (luop,file=filna,status='NEW')
      open (lupf,file=filnb,status='NEW')
      open (lusl,file=filnc,status='NEW')
c
c     main file
c
      write (luop,10) fnam,runname,banfil
   10 format('PLANE PARALLEL SHOCK MODEL: S3',3x,
     & '(separate el. and ion temps)',/,' FILE  ',a16,/,
     & ' RUN: ',a80,/,' INPUT: ',a20,//)
      if (usekappa) then
   20    format(/' Kappa Electron Distribution Enabled :'/
     &           ' Electron Kappa : ',1pg11.4)
        write (luop,20) kappa
      endif
c
c     write main header
c
      write (luop,30) dhpr,xhpr,tepr,tmag
   30 format(/' PRE-SHOCK CONDITIONS :',t28,
     & 'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     & 'FRACTIONAL IONISATION :',1pg10.3/' ====================',t33,
     & 'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     & 'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write (luop,40) zgas
   40 format(9x,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun) ',
     &' and the proto-ionisation used.')
c
      call wionabal (luop, propop)
c
      write (luop,50)
   50 format(/' MOD',t6,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t46,'FQHI',t55,'FQHEI',t65,'FQHEII',t75,'DILUF')
      write (luop,60) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii,
     &diluf
c
   60 format(1x,a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
      vshockm=vshoc/1.d5
      write (luop,70) tepo,vshockm
   70 format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' =====================')
c
c     apn file
c
      write (lupf,80) fnam,runname
   80 format(
     &'PLANE PARALLEL SHOCK MODEL S3 (Separate el. and ion temps)',
     & /,' REF.FILE: ',a16,/,' RUN: ',a80,//)
c
      write (lupf,30) dhpr,xhpr,tepr,tmag
      write (lupf,40) zgas
      call wionabal (lupf, propop)
c
c     allion file
c
      write (lusl,80) fnam,runname
      write (lusl,30) dhpr,xhpr,tepr,tmag
      write (lusl,40) zgas
      call wionabal (lusl, propop)
c
c
      close (luop)
      close (lupf)
      close (lusl)
c
      return
c
      end
