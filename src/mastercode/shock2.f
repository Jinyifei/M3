cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c    INCLUDES DIFFUSE FIELD , USES CASE B
c
c      *COMPUTATIONS ARE PERFORMED IN SUBROUTINE COMPSH2
c
c    CALL SUBR.:VELSHOCK,COMPSH2,POPCHA,PHOTSOU
c
c
      subroutine shock2 ()
c
      include 'cblocks.inc'
c
c
c           Variables
c
      real*8 dhpr,epotmi,hmag
      real*8 qdl,tepo,tepr
      real*8 tmag,vshockm,xhpr
c
      integer*4 i,j
      integer*4 lucty,nstma2
c
      character carac*46
      character ilgg*4, caract*4,ibell*4
      character banfil*32, model*32
c
      lucty=5
c
      nstma2=mxnsteps-50
      jcon='YES'
      jspot='NO'
      model='Shock 2'
      epotmi=iphe
c
      ibell(1:4)=char(7)
c
c
      write (*,10)
   10 format(///' SHOCK MODEL SELECTED : S2'/
     &' PHOTON DIFFUSE FIELD INCLUDED')
c
c     set up ionisation condiditons
c
      model='proto-ionisation'
      call popcha (model)
      model='Shock 2'
c
      call copypop (pop, propop)
c
      carac(1:33)='   '
      j=1
      do 30 i=3,atypes
        j=1+index(carac(1:33),'   ')
        caract='   '
   20 format(a2)
        if ((ipote(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) then
          write (caract,20) elem(i)
        endif
        carac(j:j+2)=caract//'   '
   30 continue
      i=j+2
c
   40 format(a2)
   50 format(a)
c
   60 write (*,70) carac(1:i)
   70 format(//' Do you allow the elements : ',a/
     &' to recombine freely to neutral state (Y/N) ? ',$)
      read (*,50,err=60) ilgg
      ilgg=ilgg(1:1)
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
      if (ilgg.eq.'Y') goto 90
      if (ilgg.ne.'N') goto 60
      do 80 i=3,atypes
        arad(2,i)=dabs(arad(2,i))
        if ((ipote(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) arad(2,i)=-
     &   arad(2,i)
   80 continue
   90 continue
c
c
      write (*,100)
  100 format(//' Select shock conditions :')
  110 write (*,120) (propop(2,1))
  120 format(//' Give initial hydrogen density and ',
     &'ionised fraction (',1pg14.7,').'/
     &' Fractions less than zero are taken as logs'/
     &' Fractions greater than one cause default to be used'/
     &' :: ',$)
      read (*,*) dhpr,xhpr
      if (dhpr.le.0.d0) dhpr=10.d0**dhpr
      if (xhpr.gt.1.d0) xhpr=propop(2,1)
      if (xhpr.lt.0.d0) xhpr=10.d0**xhpr
      if (((xhpr.lt.0.0d0).or.(dhpr.le.0.0d0)).or.(xhpr.gt.1.d0)) goto
     &110
c
c
  130 write (*,140)
  140 format(//' Give initial and postshock temperatures : ',$)
      read (*,*) tepr,tepo
      if (tepr.le.10.d0) tepr=10.d0**tepr
      if (tepo.le.10.d0) tepo=10.d0**tepo
c
      if ((tepo.le.tepr).or.(tepr.le.100.0d0)) then
        write (*,150) ibell
  150 format(' UNSUITABLE ANSWER(S) ****************************',a)
        goto 130
      endif
c
c
  160 write (*,170)
  170 format(//' Give magnetic field (micro-Gauss) : ',$)
      read (*,*,err=160) tmag
      if (tmag.lt.0.0d0) goto 160
      hmag=tmag*1.0d-6
c
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
c
      call velshock (dhpr, xhpr, tepr, tepo, hmag)
      vshockm=vshoc/1.d5
      write (*,180) vshockm
  180 format(//' Shock velocity :',f8.1)
c
c
c    ***PHOTOIONIZATION SOURCES
c
  190 write (*,200)
  200 format(//' Any photon source at the shock front',
     &'(Y/N) ? ',$)
      read (*,40,err=190) ilgg
      ilgg=ilgg(1:1)
      ilgg=ilgg(1:1)
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
c
      if ((ilgg.ne.'Y').and.(ilgg.ne.'N')) goto 190
      if (ilgg.eq.'N') goto 230
c
      call photsou (model)
c
  210 write (*,220)
  220 format(/' Dilution factor (<=0.5) : ',$)
      read (*,*,err=210) diluf
      if ((diluf.lt.0.0).or.(diluf.gt.0.5)) goto 210
c
c
  230 continue
c
c
  240 write (*,250) nstma2
  250 format(//' Give number of steps (MAX:',i4,') '/
     &' and give number of global iterations : ',$)
      read (*,*,err=240) loop,iter
c
c
      if ((loop.gt.nstma2).or.(loop.lt.8)) goto 240
      if (iter.le.0) goto 240
c
      iel1=0
      iel2=0
      iel3=0
      iel4=0
      ieln=4
  260 write (*,270)
  270 format(//' Choose four elements (atomic numbers)'/
     &' to be followed : ',$)
c
      read (*,*) iel1,iel2,iel3,iel4
      if (iel1.gt.28) goto 260
      if (iel2.gt.28) goto 260
      if (iel3.gt.28) goto 260
      if (iel4.gt.28) goto 260
c
c
      if (iel1.le.0) goto 260
      if (iel2.le.0) goto 260
      if (iel3.le.0) goto 260
      if (iel4.le.0) goto 260
c
      iel1=zmap(iel1)
      iel2=zmap(iel2)
      iel3=zmap(iel3)
      iel4=zmap(iel4)
c
      if (iel1.eq.0) goto 260
      if (iel2.eq.0) goto 260
      if (iel3.eq.0) goto 260
      if (iel4.eq.0) goto 260
c
c    ***SET TYPE OF ENDING FOR THE PROGRAM
c
      xhmin=0.01d0
      qdl=1.5d0+((0.2d0*(tepo/1.d5))*(nstma2/(loop)))
      write (*,280) vshockm,qdl
  280 format(//' SHOCK VELOCITY :',f8.1,' KM/SEC',7x,'Q FINE-TUNING :'
     &,f7.2)
  290 write (*,300) xhmin,elem(iel1)
  300 format(//' CHOOSE TYPE OF ENDING FOR THE PROGRAM :'/t6,
     &' A-  NORMAL ENDING  (FINAL H+/H :',2pf4.1,'% )'/t6,
     &' B-  ACCORDING TO FINAL TEMPERATURE'/t6,
     &' C-  ACCORDING TO ELAPSED TIME'/t6,
     &' D-  ACCORDING TO THE IONISED FRACTION OF : ',a2/t6,
     &' E-  ACCORDING TO LINE INTENSITIES'/t6,' R-  RE-INITIALIZE'/' '
     &,t55,':: ',$)
      read (*,310,err=290) jfin
  310 format(a)
      jfin=jfin(1:1)
c
      if (jfin.eq.'a') jfin='A'
      if (jfin.eq.'b') jfin='B'
      if (jfin.eq.'c') jfin='C'
      if (jfin.eq.'d') jfin='D'
      if (jfin.eq.'e') jfin='E'
      if (jfin.eq.'r') jfin='R'
c
      if ((jfin.lt.'A').or.((jfin.gt.'E').and.(jfin.ne.'R'))) goto 290
      if (jfin.eq.'R') return
      if (jfin.eq.'A') goto 370
      if (jfin.ne.'D') goto 330
c
      write (*,320)
  320 format(/' WHICH IONISING STAGE : ',$)
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
  360 format(/' GIVE FINAL VALUE : ')
      read (*,*,err=290) fval
      if (fval.lt.0.d0) fval=10**fval
c
      if (jfin.eq.'B') texi=fval
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
      read (*,50,err=380) jpre
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
  420 format(//' Give a name/code for this run: ',$)
      write (*,420)
      read (*,410) runname
c
c
      banfil='Interactive'
c
      call compsh2 (lucty, dhpr, xhpr, tepr, tepo, hmag, banfil)
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
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c     INCLUDES DIFFUSE FIELD , USES CASE B
c     CALL SUBROUTINES  ROOTS,COOL,TIMION,INTVEC,TAUDIST,EQUION,
c                       WRSPPOP,SUMDATA,AVRDATA,VELSHOCK,ZETAEFF,
c                       LOCALEM,NEWDIF,TOTPHOT,ZERBUF,PREION
c         AVERINTO,copypop
c
c
      subroutine compsh2 (luter, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      include 'cblocks.inc'
c
      real*8 mu2,dhpr, xhpr, tepr, tmag,t,de,dh,dr
      real*8 cspd,tnloss,press,en,ue,rhotot,wmol,fmui
c
      real*8 adj,canor,cavx,cfm,chasl,chate,cofcrit,convg
      real*8 d,de2,dedhma,dedhmi,def,dexi,dfa0,dfaa,dfv1,dh2
      real*8 distav,dl0,dlf,dli,dlosa,dlosav,dq,dqro
      real*8 drdw,drta,drup,dstep,dte,dtema,dtemi
      real*8 dv,dva,dvdw,dvi,dvup,eac,egairr
      real*8 el0,elf,eli,elosrr,fa,ff,fgae,fi,floe,fmi,fmi0
      real*8 fnorm,frdv,frdv0,frdw,frta,frto,ftim
      real*8 ga0,gaf,gai,gam,hmag,humag
      real*8 prec1,prec2,prec3,pros,qdift,qdl,qf,qs
      real*8 rad,rdis,rfvl,rho,rsl
      real*8 tdw,tefin,tei,tepo,tl0,tlf,tli
      real*8 tlosa0,tlosrr,tsa,tsfin,tsmax,tsrsl
      real*8 tstau,tstep,tstep0,tstepi,tsterr
      real*8 tup,u,u2,vcons,velav,vshockm,vsi
      real*8 wei,xhav,xhy,xhyf
c
      integer*4 i,i107,ic107,iflag,irsl,j,l,m,mi
      integer*4 luop,lupf,lusl,np
      integer*4 luter,luions(4),lsp,nstma2
c
      character imod*4, jmod*4, lmod*4, mmod*4
      character carac*46, pfx*16, tab*4
      character filna*32, filnb*32, filnc*32,filn(4)*32
      character banfil*32, fnam*32,caller*4,wmod*4,ispo*4
c
c     External Functions
c
      real*8 fcrit,feldens,fmua,fpressu,frho
c
      jspot='NO'
      jcon='YES'
      ispo='SII'
c init
      tstepi=0.d0
      tstep0=0.d0
      tli=0.d0
      tefin=0.d0
      rfvl=0.d0
      prec1=0.d0
      prec2=0.d0
      prec3=0.d0
      ga0=0.d0
      el0=0.d0
      dli=0.d0
      dl0=0.d0
      dfv1=0.d0
      dfaa=0.d0
      dfa0=1.d0
c
      nstma2=mxnsteps-50
c
      if ((iel1.eq.5).and.(iel2.eq.5)) ispo=' OI'
c
      dedhmi=0.0d0
      dedhma=0.0d0
      do 10 j=1,atypes
        dedhma=dedhma+zion(j)
        if (arad(2,j).le.0.d0) dedhmi=dedhmi+zion(j)
c
   10 continue
c
      lsp=26
      luop=21
      lupf=27
      lusl=45
      luter=6
c
      do i=1,ieln
        luions(i)=21+i
      enddo
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
      call s2fheader (fnam, filna, filnb, filnc, filn, luop, lupf, lusl,
     & luions, dhpr, xhpr, tepr, tepo, tmag, banfil)
c
c
      do 400 mi=1,iter
        rad=1.d38
        wdil=diluf
        fi=1.0d0
c
        qdift=0.0d0
        if (mi.eq.1) goto 80
        dva=0.0d0
        tei=tepr
        dr=1.d0
c
c    ***SETTING UP OF PREIONISATION CONDITIONS
c
c
        if (jfin.eq.'A') tsmax=1.d2*tsmax
c
        open (luop,file=filna,status='OLD',access='APPEND')
        write (luop,20) jpre
c
   20 format(//' SETTING UP OF THE NEXT PREIONISATION CONDITIONS : '
     &,a1)
c
        lmod='UP'
        call localem (tei, depr, dhpr)
c      call contcont(tei, depr, dhpr)
c
c
        call totphot (tei, dhpr, fi, rad, dr, dva, wdil, lmod)
c
        if (jpre.ne.'U') goto 40
        do 30 l=1,4
c
          vsi=vshoc
c
          call preion (luter, luop, tsmax, vshoc, fi, dhpr, depr, tepr,
     &     qdift, d)
c
          xhpr=pop(2,1)
c
          call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
          if ((dabs(vsi-vshoc)/dmax1(vsi,vshoc)).lt.0.05) goto 40
   30   continue
c
   40   continue
c
        if (jpre.ne.'F') goto 50
        call copypop (popf, pop)
        xhpr=pop(2,1)
c
        tepr=tefin
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
   50   continue
c
        if (jpre.ne.'T') goto 60
c
        call equion (tepr, depr, dhpr)
        xhpr=pop(2,1)
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
   60   continue
c
        if (jpre.ne.'C') goto 70
c
        call copypop (popini, pop)
c
        call velshock (dhpr, xhpr, tepr, tepo, hmag)
c
   70   continue
c
c
        close (luop)
c
        call zerbuf
c
c
   80   continue
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
        dhy(1)=dhpr
        dhy(2)=dhpo
        dist(1)=0.0d0
        dist(2)=0.0d0
        deel(1)=depr
        deel(2)=depo
        timlps(1)=0.0d0
c
        qdl=1.5d0+((0.2d0*(tepo/1.d5))*(nstma2/(loop)))
        fnorm=qdl/(1.0d0-dexp(-qdl))
        vcons=-dabs((vpo-vfin)/(loop-4.0d0))
c
c
        qf=qst
        open (luop,file=filna,status='OLD',access='APPEND')
c
        write (luop,90) loop,mi,tepr,qdift
   90 format(//i4,' STEPS',t40,' ITERATION #' ,i2,10x,'TEPR : '
     &,f8.0,10x,'PRESHOCK QTOT : ' ,1pg10.3/t40,' ===========' /)
c
        write (*,100) mi,tepr,qdift,fnam
  100 format(/' ITER. #',i2,5x,'TEPR :',f8.0,'K',5x,'PRESHOCK QTOT : '
     &,1pg10.3,5x,a12)
        write (luop,110) jpre
  110 format(' ',t9,'INITIAL STATE OF IONISATION ' ,
     &'FOR ALL ELEMENTS  (JPRE : ' ,a1,' )  :'/)
c
        call wionabal (luop, popi)
c
        carac(1:33)='I  II IIIIV V  VI '
        if (luter.gt.0) write (luter,120)
  120 format(/' STEP NUMBERS AS THEY ARE COMPLETED :' /' L*',t8,'<TE>'
     &,t18,'<TLOSS>',t27,'<DLFR>',t36,'TSTEP',t45,'<DENS.>',t55,
     &'<EL.DENS>',t69,'DV')
        write (luop,130) elem(1),carac(1:3),elem(1),carac(4:7),
     &   (elem(iel1),carac(j:j+3),j=1,16,3),(elem(iel2),carac(j:j+3),j=
     &   1,16,3)
  130 format(//' STEP NUMBERS AS THEY ARE COMPLETED :' //' L*',t7,
     &'<LOSS.R.>',t18,'TSTEP',t26,14(1x,a2,a3,1x)/)
c
        close (luop)
c
c     DOMWNSTREAM COMPUTATION ; STEP INDEX : M
c
        drta=0.0d0
        frta=0.007d0
        eac=100.0d0
        canor=4.0d0
        frdv=1.0d0
        ff=1.0d0
        frto=0.0d0
        fmi=fnorm*dexp(-((qdl*frto)/loop))
        convg=0.004d0
        do 140 i=1,5
  140     cav(i)=canor
        cavx=canor
c
c
        do 350 m=2,299
c
          iflag=1
          irsl=0
c
          if (luter.gt.0) write (luter,*)
c
          qs=qf
          dtema=0.25d0
          dtemi=0.033d0
          if ((m.lt.(loop/3.0)).and.((((deel(m)/dhy(m))-dedhmi)/dedhma)
     &     .lt.0.2)) dtemi=dtemi/dabs(dlog(dmax1(((deel(m)/dhy(m))-
     &     dedhmi)/dedhma,1.d-10)))
          fmi0=fmi
          fmi=fnorm*dexp(-((qdl*frto)/loop))
          dte=dabs(te(m)-te(m-1))/dmax1(te(m-1),te(m))
          frdv=(frdv*(fmi0/fmi))*(((cavx+0.05d0)/canor)**(-0.5d0))
          fmui=0.97d0*dsqrt(dtemi/(dte+1.d-9))
          if (dte.lt.(1.05d0*dtemi)) frdv=frdv*dmin1(2.0d0,fmui)
          if (frdv.gt.1.0d0) frdv=1.0d0
          dvi=fmi*vcons
c
c    ***DERIVES GUESSED VALUES FOR : DE,DH,TLOSS,DR,TSTEP AND T
c
          if (dabs(dvi/veloc(m)).gt.0.2d0) dvi=-(veloc(m)*0.2d0)
c
  150     ic107=iflag
c
c      *STEP INDEX M=2
c
          if (m.eq.2) then
            xh(2)=popi(2,1)
            dh=dhy(2)
            de1=feldens(dh,popi)
            de=de1
            deel(2)=de
            t=tepo
            tstep=0.0d0
c
c
            dr=1.d5
            lmod='DW'
            call localem (t, de, dh)
c      call contcont (t,de,dh)
            call totphot (t, dh, fi, rad, dr, dva, wdil, lmod)
            call cool (t, de, dh)
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
c
            dfaa=dabs((fmi*frdv)/dfv1)
            tloss=tl0
            dh=dhy(3)
            de=deel(3)
            t=te(3)
            dlosa=dl0
            te(4)=0.0
c
          else
c
c      *STEP INDEX M>3
c
            if (((te(m).lt.200.d0).and.(dexi.lt.7.d-2)).or.(dexi.lt.2.d-
     &       2)) dtema=0.4d0
            dfaa=dabs((fmi*frdv)/dfv1)
            adj=dfaa/dfa0
c
            if ((dl0.lt.0.4d0).and.(dli.lt.0.4d0)) adj=1.0d0
c
            if (irsl.eq.0) then
              dlosa=dl0+((0.5d0*(dl0-dli))*dfaa)
              if (dabs(dlosa).gt.1.d0) dlosa=dsign(1.d0,dlosa)
              tsterr=(tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36))))
     &         *adj
c
              if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.tstau)) tsterr=
     &         tstau
c
              wei=fcrit(dlosa,eac)
              wei=fcrit(dlosa/(1.d0+(10.d0*wei)),eac)
              tstep=dexp(((1.d0-wei)*dlog(tsterr+1.d-30))+(wei*((0.2d0*
     &         dlog(tstau))+(0.8d0*dlog(tsterr+1.d-30)))))
            else
              dlosa=dl0+(0.25d0*((dlf-dl0)+((dl0-dli)*dfaa)))
              if (dabs(dlosa).gt.1.d0) dlosa=dsign(1.d0,dlosa)
              tsterr=(tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36))))
     &         *adj
c
              if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.tstau)) tsterr=
     &         tstau
c
              wei=fcrit(dlosa,eac)
              wei=fcrit(dlosa/(1.d0+(10.d0*wei)),eac)
              tstep=dexp(((1.d0-wei)*dlog(tsterr+1.d-30))+(wei*dlog(tsa+
     &         1.d-30)))
            endif
c
c
            tloss=tl0+((0.5d0*(tl0-tli))*dfaa)
            tlf=tloss+(tloss-tl0)
            if ((tli/tl0).gt.0.d0) then
              dh=dhy(m)+(((0.5d0*(dhy(m)-dhy(m-1)))*(dhy(m)/dhy(m-1)))*
     &         dfaa)
              t=te(m)+(((0.5d0*(te(m)-te(m-1)))*(te(m)/te(m-1)))*dfaa)
c
            else
c
              dh=dhy(m)-((0.5d0*(dhy(m)-dhy(m-1)))*dfaa)
              t=te(m)-((0.5d0*(te(m)-te(m-1)))*dfaa)
            endif
            t=dmax1((1.d0-(dtema/2.d0))*te(m),dmin1((1.d0+((dtema*(1.d0+
     &       dtema))/2.d0))*te(m),t))
            de=(deel(m)*dh)/dhy(m)
            te(m+1)=0.0d0
c
          endif
c
c    ***DERIVES AVERAGES FOR DH,DE,MU,RHO AND TE,
c     COMPUTES CONDITIONS AT THE END OF TIMESTEP
c
  160     frdv0=frdv
          dv=dvi*dsign(frdv,tl0)
          dva=dabs(dv)
          veloc(m+1)=veloc(m)+dv
          velav=veloc(m)+(dv/2.d0)
          u=veloc(m+1)
          u2=u*u
          dstep=((veloc(m+1)+veloc(m))*tstep)/2.d0
          dr=0.5d0*dstep
c
c
          if (iflag.gt.ic107) goto 190
          open (lupf,file=filnb,status='OLD',access='APPEND')
c
          if (luter.gt.0) write (luter,170) m,t,tloss,dlosa,tstep,dh,de,
     &     dv
  170 format(' ',i3,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
          write (lupf,180) m,iflag,ic107,te(m),tl0,tlf,dlosa,dv,tstep,
     &     dh,de,frdv,ff
  180 format(' ',3(i3),0pf9.0,9x,7(1pg11.3),2(1pg9.2))
          close (lupf)
c
  190     continue
c
          lmod='DW'
          if (iflag.eq.ic107) then
            call localem (t, de, dh)
c         call contcont(t,de,dh)
          endif
c
          call totphot (t, dh, fi, rad, dr, dva, wdil, lmod)
          call copypop (popi, pop)
          call timion (t, def, dh, xhyf, tstep)
          call copypop (pop, popf)
c
c
          xh(m+1)=xhyf
          xhav=(xh(m)+xh(m+1))/2.d0
          humag=(hmago*vpo)/u
          gam=(((pres*u)/fm)-u2)-humag
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
          if (te(m+1).lt.((1.d0-dtema)*te(m))) te(m+1)=dmax1(te(m+1),
     &     0.15d0*te(m))
          t=(te(m+1)+te(m))/2.d0
          t=dmax1((1.d0-(dtema/1.8d0))*te(m),dmin1((1.d0+((dtema*(1.d0+
     &     dtema))/1.8d0))*te(m),t))
          dte=dabs(te(m+1)-te(m))/dmax1(te(m),te(m+1))
          qf=(2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+u2)
          dq=(qf-qs)*0.5d0
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
          call averinto (wei, popi, popf, pop)
          call localem (t, de, dh)
c      call contcont(t,de,dh)
c
          if ((dabs(dlosa).lt.0.99d0).or.(dabs(dl0).lt.0.99d0)) then
            call totphot (t, dh, fi, rad, dstep, dva, wdil, lmod)
          endif
c
          call taudist (dh, fi, frta, drta, rad, pop, mmod)
          call copypop (popf, pop)
          call cool (te(m+1), de2, dh2)
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
          rsl=tl0/tlf
c
          tsfin=dabs((2.d0*dqro)/(tl0+tlf))
          tsrsl=dmin1(tsfin,dabs((2.d0*dqro)/dmax1(dabs(tl0),dabs(tlf)))
     &     )
c
c      if photon field is very small, then skipbin may be true for all
c     bins, in this case drta will be exactly zero, so tstau will be
c     set equal to the dynamic timescale from the Rankine-Hugenoit
c     eqns.  Note; once the step size increases then the field will
c     build up and drta should return a non zero result.
c
c      if (drta.gt.0.d0) then
c         tstau = (drta/((veloc(m)+veloc(m+1))/2.d0))+1.d-20
c      else
          tstau=tsfin
c      endif
c
          iflag=iflag+1
          i107=iflag-ic107
c
c    ***TEST CONVERGENCE
c
          pros=dmax1(1.d0,((i107*fcrit(dlosa,10.d0))+1.d-4)**0.35d0)
          cofcrit=dmin1(1.001d0,convg+(1.5d0*fcrit(dlosa/pros,3.d1)))
          if ((((rsl.lt.0.d0).and.(dabs(rsl).lt.0.5d0))
     &     .and.(dte.gt.dtemi)).and.(irsl.lt.2)) goto 200
          if (dte.gt.dtema) goto 200
          floe=dabs((eloss-elosrr)/dmax1(dabs(eloss),dabs(elosrr)))
          fgae=dabs((egain-egairr)/dmax1(egain,egairr,1.d-33))
          ftim=dabs((tsfin-tsterr)/dmax1(tsfin,tsterr))
          if ((ftim.lt.convg).or.((i107.gt.1).and.(ftim.lt.cofcrit)))
     &     goto 270
c
c    ***ITERATES IF NECESSARY ON TIME STEP TO IMPROVE CONVERGENCE
c
  200     continue
          chate=(0.8d0*dtemi)/dte
          if (rsl.gt.0.d0) then
            if (((((rsl.le.3.d0).or.(irsl.ge.2)).or.(dte.le.dtemi))
     &       .or.(i107.ge.2)).or.(dte.gt.dtema)) goto 210
            chasl=0.5d0
            irsl=irsl+1
            frdv=frdv0*dmax1(chate,chasl)
  210       wei=fcrit(dlosa,eac)
            wei=fcrit(dlosa/(1.d0+(1.d1*wei)),eac)
            tsa=dexp(((1.d0-wei)*dlog(tsfin+1.d-30))+((wei*0.5d0)*
     &       (dlog(tstau)+dlog(tsterr+1.d-30))))
            tsa=dmax1(tsa,(fcrit(dlosav,1.d0)*tsterr)/1.d1)
            if ((i107.gt.3).and.(dabs(dlosa).gt.0.25d0)) then
              wei=dmin1(0.5d0,1.6d0*((i107/10.d0)**2))*(1.d0-
     &         fcrit(dlosa,eac/4.d0))
c
              tsa=((1.d0-wei)*tsa)+(wei*dmin1(1.45d0*tsa,dmax1(0.7d0*
     &         tsa,tsterr)))
c
            endif
          else
            if (((dabs(rsl).ge.0.5d0).or.(irsl.ge.2)).or.(dte.le.dtemi))
     &        goto 220
            chasl=dabs(rsl)/5.0d0
            irsl=irsl+1
            frdv=frdv0*dmax1(chate,chasl)
  220       continue
            fa=dabs(tsterr-tsfin)/dmax1(tsterr,tsfin)
            tsa=(tsfin+tsterr)/2.d0
            cfm=0.5d0-(0.35d0*fcrit(dlosa,2.5d1))
            if (fa.gt.cfm) tsa=tsterr*(1.d0+dsign(cfm,dlog(epsilon+
     &       (tsfin/(tsterr+1.d-30)))))
            tsa=dmin1(tstau,dmax1(tsa,tsrsl/2.d0))
c
          endif
c
c
          tstep=tsa
          if (dte.gt.dtema) then
            frdv=frdv0*dmin1(0.6d0*(1.d0-dtema),0.7d0*((dtema/dte)**
     &       1.5d0))
          else if ((i107.gt.7).and.((((deel(m)/dhy(m))-dedhmi)/dedhma)
     &     .lt.0.2d0)) then
            frdv=frdv0/8.0d0
          else if (i107.gt.7) then
            frdv=frdv0/3.0d0
c
          endif
c
c
          open (lupf,file=filnb,status='OLD',access='APPEND')
c
          if (luter.gt.0) write (luter,230) t,tloss,dlosa,tsterr,tsfin,
     &     tstau,rsl
  230 format('    ',f10.0,1pg11.3,0pf7.3,4(1pg11.3))
          write (lupf,240) m,iflag,ic107,te(m),te(m+1),tl0,tlf,dlosa,dv,
     &     tsterr,tsfin,tstau,tstep,rsl
  240 format(' ',3(i3),2(0pf9.0),8(1pg11.3),0pf8.3)
          close (lupf)
          if (iflag.gt.30) goto 250
          if (frdv.ne.frdv0) goto 150
c
c    ***IF CONVERGENCE UNACHIEVABLE : EXIT GRACEFULLY
c
          goto 160
  250     continue
          write (*,260) m
  260 format(/' <<<<<<<<<<<< NO CONVERGENCE AT STEP#',i4/)
c
  270     continue
c
c    ***WHEN CONVERGENCE ACHIEVED , PREPARES NEXT TIME STEP
c
          do 280 j=1,4
  280       cav(j)=cav(j+1)
          cav(5)=(fcrit(dlosa,8.d0)*canor)+((1.d0-fcrit(dlosa,8.d0))*
     &     i107)
          cavx=0.0d0
          i=5
          do 290 j=5,6-i,-1
  290       cavx=cavx+(j*cav(j))
c
          cavx=cavx/((i+1.d0)*(i/2.d0))
          dist(m+1)=dist(m)+dstep
          timlps(m)=timlps(m-1)+tsterr
          tsmax=5.d0*timlps(m)
          tefin=te(m+1)
c
c    ***INTEGRATES DIFFUSE FIELD,COMPUTES AVERAGE SPECTRUM AND
c     SUMS UP RELEVANT QUANTITIES
c
          dexi=((de/dh)-dedhmi)/dedhma
          rdis=(dist(m)+dist(m+1))/2.d0
          tdw=t
          tup=dsqrt(tepo*t)
          drdw=dstep
          dvdw=dva
          drup=dist(m)
          dvup=dsqrt(vpo*u)
          frdw=0.5d0
c
c
          jmod='LODW'
          imod='ALL'
          call cool (t, de, dh)
          call zetaeff (fi, dh)
          call newdif (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     &     frdw, jmod)
c
c
c    ***OUTPUT PART OF THE DATA AND TEST ENDING
c
          call sumdata (t, de, dh, fi, drdw, drdw, rdis, imod)
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
            hsii(m)=(fluxm(1,su2)+fluxm(2,su2))/(fluxh(1)+epsilon)
          endif
          if (ox1.ne.0) then
            if (ispo.eq.' OI') hsii(m)=fluxm(3,ox1)/(fluxh(1)+epsilon)
          endif
          if (jfin.ne.'E') goto 300
          if (istf.eq.1) rfvl=hoiii(m)
          if (istf.eq.2) rfvl=hoii(m)
          if (istf.eq.3) rfvl=hnii(m)
          if (istf.eq.4) rfvl=hsii(m)
  300     continue
c
          open (luop,file=filna,status='OLD',access='APPEND')
          open (lupf,file=filnb,status='OLD',access='APPEND')
          open (lusl,file=filnc,status='OLD',access='APPEND')
c
          tab=char(9)
          do i=1,ieln
            open (luions(i),file=filn(i),status='OLD',access='APPEND')
c
  310 format(i3.3,a1,1pg12.5,a1,31(1pg11.4,a1))
            if (i.eq.1) then
              write (luions(i),310) m,tab,t,tab,(pop(j,iel1),tab,j=1,
     &         maxion(iel1))
            endif
            if (i.eq.2) then
              write (luions(i),310) m,tab,t,tab,(pop(j,iel2),tab,j=1,
     &         maxion(iel2))
            endif
            if (i.eq.3) then
              write (luions(i),310) m,tab,t,tab,(pop(j,iel3),tab,j=1,
     &         maxion(iel3))
            endif
            if (i.eq.4) then
              write (luions(i),310) m,tab,t,tab,(pop(j,iel4),tab,j=1,
     &         maxion(iel4))
            endif
c
            close (luions(i))
c
          enddo
c
          if (luter.gt.0) write (luter,230) t,tloss,dlosa,tstep,tsfin,
     &     tstau,rsl
          write (luop,320) m,tloss,tstep,pop(1,1),pop(2,1),(pop(j,iel1),
     &     j=1,6),(pop(i,iel2),i=1,6)
  320 format(' ',i3,1pg10.2,1pg9.2,1x,14(0pf7.4))
          write (lupf,240) m,iflag,ic107,te(m),te(m+1),tl0,tlf,dlosa,dv,
     &     tsterr,tsfin,tstau,tstep,rsl
          write (lupf,330) cavx,caseab(1),dexi,fgae,floe,ftim,cofcrit,t,
     &     dte
  330 format(' CAVX :',0pf6.3,3x,'CASE A,B :' ,0pf4.2,3x,'FEXIT :'
     &,1pg10.3,3x,4(0pf7.4),1(0pf10.0),1(0pf7.4))
c
c
          press=fpressu(t,dh,pop)
c
          en=zen*dh
c
          ue=3/2*(en+de)*rkb*t
          tnloss=tloss/(en*de)
c
          rhotot=frho(de,dh)
          cspd=dsqrt(5/3*press/rhotot)
          wmol=rhotot/(en+de)
c
          tab=char(9)
          write (lusl,340) m,tab,dist(m),tab,dstep,tab,t,tab,de,tab,dh,
     &     tab,en,tab,tloss,tab,tnloss,tab,rhotot,tab,press,tab,ue,tab,
     &     cspd,tab,wmol,tab
c
  340 format(i3.3,a1,13(1pg12.5,a1))
c
c      wmod = 'FILE'
c      call wmodel(luter,t,de,dh,dstep,wmod)
c
c      call wionabal(luter,pop)
c
          close (luop)
          close (lupf)
          close (lusl)
c
c
c      wmod = 'SCRN'
c      call wmodel(luter,t,de,dh,dstep,wmod)
c
c first upstream integrated field
c
c     diagnostics turned off
c
          if (.false.) then
            dva=dsqrt(vpo*u)
            lmod='UP'
c
            call totphot (t, dh, fi, rad, dstep, dva, wdil, lmod)
c
            caller='S2'
            pfx='psoup'
            np=5
            wmod='REAL'
c
            call wpsou (caller, pfx, np, wmod, t, de, dh, veloc(m+1),
     &       dva, vshoc, dist(m+1), dstep, timlps(m), tsterr, (2.d0/
     &       4.d0), tphot)
c
c     then local downstream
c
            dva=0.0d0
            lmod='DW'
c
            call totphot (t, dh, fi, rad, dstep, dva, wdil, lmod)
c
            caller='S2'
            pfx='psodw'
            np=5
            wmod='REAL'
c
            dv=dvi*dsign(frdv,tl0)
c
            call wpsou (caller, pfx, np, wmod, t, de, dh, veloc(m+1),
     &       dv, vshoc, dist(m+1), dstep, timlps(m), tsterr, (2.d0/4.d0)
     &       , tphot)
          endif
c
c
          if (iflag.gt.30) goto 360
          if ((t.lt.1500.0).and.(dexi.lt.xhmin)) goto 360
          if ((t.lt.200.0)) goto 360
          if ((jfin.eq.'B').and.(te(m+1).le.texi)) goto 360
          if (((jfin.eq.'C').and.(m.gt.5)).and.(fval.lt.timlps(m)))
     &     goto 360
          if (((((jfin.eq.'D').and.(m.gt.5)).and.(pop(istf+1,iel1)
     &     .lt.prec2)).and.(pop(istf,iel1).lt.prec1)).and.(pop(istf,
     &     iel1).lt.fval)) goto 360
          if ((((jfin.eq.'E').and.(m.gt.(loop/5.0))).and.(rfvl.lt.fval))
     &     .and.(prec3.gt.fval)) goto 360
          prec2=pop(istf+1,iel1)
          prec1=pop(istf,iel1)
          prec3=rfvl
c
          call copypop (popf, popi)
          call copypop (popf, pop)
c
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
          frto=frto+frdv
c
c
c     END MAIN LOOP IN M
c
c
  350   continue
c
        m=m-1
  360   continue
c
c    ***OUTPUT QUANTITIES STORED IN ARRAYS
c
        open (luop,file=filna,status='OLD',access='APPEND')
        write (luop,370) ispo
  370 format(//' LIST OF STORED ' ,'QUANTITIES : ' //' L*',t9,'<TE>'
     &,t16,'<FHII>',t24,'<DENS.>',t34,'<EL.DENS>',t43,'<DIST>',t53,'<V>'
     &,t63,'<ELPS.TM>',t74,'OIII/H-B',t84,'OII/H-B',t94,'NII/H-A'
     &,t104,a3,'/H-A'/)
        do 390 l=2,m
          distav=(dist(l)+dist(l+1))/2.0
          velav=(veloc(l+1)+veloc(l))/2.0
          t=(te(l)+te(l+1))/2.0
          dh=(dhy(l)+dhy(l+1))/2.0
          de=(deel(l)+deel(l+1))/2.0
          xhy=(xh(l)+xh(l+1))/2.0
          write (luop,380) l,t,xhy,dh,de,distav,velav,timlps(l),hoiii(l)
     &     ,hoii(l),hnii(l),hsii(l)
  380 format(' ',i3,f10.0,f7.4,9(1pg10.3))
  390   continue
c
c
        close (luop)
c
c     output upstream field photon source file
c
        dva=dsqrt(vpo*u)
        dr=1.d0
        lmod='UP'
c
        call totphot (tepr, dhpr, fi, rad, dr, dva, wdil, lmod)
c
        caller='S2'
        pfx='psoup'
        np=5
        wmod='REAL'
c
        call wpsou (caller, pfx, np, wmod, tepr, depr, dhpr, veloc(m),
     &   vpo, vshoc, dist(m), dstep, timlps(m), tsterr, (2.d0/4.d0),
     &   tphot)
c
c     main loop
c
  400 continue
c
      open (luop,file=filna,status='OLD',access='APPEND')
c
c     formats
c
  410 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' ====================',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
c
      if (iter.le.1) goto 480
      write (luop,420)
  420 format(///' INITIAL PARAMETERS CARACTERISING THE ' ,
     &'LAST ITERATION :' /)
      write (luop,410) dhpr,xhpr,tepr,tmag
      write (luop,110) jpre
c
      call wionpop (luop, popini)
c
c
  430 format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
  440 format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')
  450 format(' ',a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
      write (luop,430)
      write (luop,440)
      write (luop,450) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii,
     &diluf
      write (luop,460) qdift
  460 format(/' ',t10,' PRESHOCK UPSTREAM IONISING FLUX : ' ,1pg10.3/)
  470  format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' =====================')
      write (luop,470) tepo,vshockm
  480 continue
c
c    ***OUTPUT SPECTRUM AND AVERAGE IONIC POP. AND TEMP.
c
c
      call avrdata
c
c
      call wrsppop (luop)
c
      close (luop)
      write (*,490) fnam
c
  490 format(/' OUTPUT PRINTED &&&&&& FILE : ' ,a12/)
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
      subroutine s2fheader (fnam, filna, filnb, filnc, filn, luop, lupf,
     & lusl, luions, dhpr, xhpr, tepr, tepo, tmag, banfil)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 vshockm
      real*8 dhpr, xhpr, tepr,tepo, tmag
c
      integer*4 i,ie,j
      integer*4 luop,lupf,lusl,luions(4)
c
      character filna*32, filnb*32,filnc*32,filn(4)*32
      character banfil*32,fnam*32
      character fn*32
      character pfx*16,sfx*4,tab*4
c
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
c
c    ***DERIVE FILENAME : shckn**.sh2
c
      fn=' '
      pfx='shckn'
      sfx='sh2'
      call newfile (pfx, 5, sfx, 3, fn)
      filna=fn(1:13)
      fnam=filna
c
c
      fn=' '
      pfx='shapn'
      sfx='sh2'
      call newfile (pfx, 5, sfx, 3, fn)
      filnb=fn(1:13)
c
c
      fn=' '
      pfx='allion'
      sfx='sh2'
      call newfile (pfx, 6, sfx, 3, fn)
      filnc=fn(1:14)
c
c
      do i=1,ieln
        if (i.eq.1) ie=iel1
        if (i.eq.2) ie=iel2
        if (i.eq.3) ie=iel3
        if (i.eq.4) ie=iel4
        fn=' '
        pfx=elem(ie)
        sfx='sh2'
        call newfile (pfx, elem_len(ie), sfx, 3, fn)
        filn(i)=fn(1:(elem_len(ie)+3+5))
      enddo
c
c
      open (luop,file=filna,status='NEW')
      open (lupf,file=filnb,status='NEW')
      open (lusl,file=filnc,status='NEW')
c
      do i=1,ieln
        open (luions(i),file=filn(i),status='NEW')
      enddo
c
c     main file
c
      write (luop,10) fnam,runname,banfil
   10 format(' ',t40,'PLANE PARALLEL SHOCK MODEL: S2',3x,
     &'(DIFFUSE FIELD INCLUDED)',t100,'FILE  ',a16/t40,
     &'RUN:',a80,t100,'INPUT ',a12//)
c
      if (usekappa) then
   20    format(/' Kappa Electron Distribution Enabled :'/
     &           ' Electron Kappa : ',1pg11.4)
        write (luop,20) kappa
      endif
c     write main header
c
      write (luop,30) dhpr,xhpr,tepr,tmag
   30 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' ====================',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write (luop,40) zgas
   40 format(' ',t9,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun) '/
     &' and the proto-ionisation used.')
c
      call wionabal (luop, propop)
c
      write (luop,50)
   50 format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')
      write (luop,60) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii,
     &diluf
c
   60 format(' ',a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
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
   80 format(' ',t40,
     &'PLANE PARALLEL SHOCK MODEL S2 (DIFFUSE FIELD INCLUDED)',10x,
     &'REF.FILE : ',a16/'RUN: ',a80/
     &t40,'=========================='/)
c
      write (lupf,30) dhpr,xhpr,tepr,tmag
      write (lupf,70) tepo,vshockm
      write (lupf,40) zgas
      call wionabal (lupf, propop)
c
c     allion file
c
      write (lusl,80) fnam,runname
      write (lusl,30) dhpr,xhpr,tepr,tmag
      write (lusl,70) tepo,vshockm
      write (lusl,40) zgas
      call wionabal (lusl, propop)
c
      tab=char(9)
c
      write (lusl,*) 'It#',tab,'distance',tab,'dr',tab,'Te',tab,'de',
     &tab,'dh',tab,'nt',tab,'Tloss',tab,'Nloss',tab,'density',tab,'press
     &ure',tab,'internal energy',tab,'sound speed',tab,'mol. wgt'
c
c
      do i=1,ieln
c
        write (luions(i),80) fnam,runname
        write (luions(i),30) dhpr,xhpr,tepr,tmag
        write (luions(i),70) tepo,vshockm
        write (luions(i),40) zgas
c
   90 format(2a4,29a6)
        if (i.eq.1) then
          write (luions(i),90) 'It#','Te ',(rom(j),j=1,maxion(iel1))
        endif
        if (i.eq.2) then
          write (luions(i),90) 'It#','Te ',(rom(j),j=1,maxion(iel2))
        endif
        if (i.eq.3) then
          write (luions(i),90) 'It#','Te ',(rom(j),j=1,maxion(iel3))
        endif
        if (i.eq.4) then
          write (luions(i),90) 'It#','Te ',(rom(j),j=1,maxion(iel4))
        endif
c
      enddo
c
c
      close (luop)
      close (lupf)
      close (lusl)
c
      do i=1,ieln
        close (luions(i))
      enddo
c
      return
c
      end
