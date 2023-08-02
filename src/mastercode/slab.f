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
c     A single slab model
c     full diffuse field, includes slab depth.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine slab
c
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 cab,de,dh,dr,dv,epotmi,exl,fi
      real*8 fron,prescc,rad,t,tex,xhyf
      real*8 tf,trea,trec,tsl,tstep1,tstep2,tstep3,wd00
      real*8 popz(mxion, mxelem),dift,tst,tstep
c
      integer*4 l,luop,np
      integer*4 i,j
c
      character imod*4, lmod*4, kmod*4, nmod*4,wmod*4
      character ilgg*4, jprin*4, ndee*4, ill*4
      character fn*32
      character fl*32,pfx*16,sfx*4,caller*4,model*32
c
c           Functions
c
      real*8 feldens,fpressu,frectim
c
      luop=23
      cab=0.0d0
      epotmi=iphe
      model='Single slab model'
c
      fn=' '
      pfx='slab'
      sfx='txt'
      call newfile (pfx, 4, sfx, 3, fn)
      fl=fn(1:12)
c
      trea=epsilon
c
      open (luop,file=fl,status='NEW')
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      jspot='NO'
      jcon='YES'
c
c     jspec is used here only in mode F, and in evoltem to
c     write out a full spectrum for each step of evoltem.
c     default to 'NO'
c
      jspec='NO'
c
   10 write (*,20)
   20 format(//' Choose a single slab model type :'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     A :  Fixed degree of ionisation at a given temp.'/
     &'     B :  Equilibrium ionisation at a fixed temp.'/
     &'     C :  Equilibrium ionisation and equilibrium temp.'/
     &'     D :  Time dependent ionisation at equilibrium temp.'/
     &'     E :  Time dependent ionisation at a fixed temp.'/
     &'     F :  Time dependent ionisation and temperature.'/
     &'  :: ',$ )
      read (*,30) ilgg
   30 format(a)
c
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
      if (ilgg.eq.'d') ilgg='D'
      if (ilgg.eq.'e') ilgg='E'
      if (ilgg.eq.'f') ilgg='F'
c
      if ((ilgg.lt.'A').or.(ilgg.gt.'F')) goto 10
c
      l=2
      do 50 i=1,atypes
        do 40 j=1,maxion(i)
   40     pop(j,i)=0.0d0
        if ((ilgg.eq.'C').or.(ilgg.eq.'B')) pop(l,i)=1.0d0
        if (((ilgg.eq.'D').or.(ilgg.eq.'E')).or.(ilgg.eq.'F')) then
          l=1
          if (ipote(1,i).lt.epotmi) l=2
          pop(l,i)=1.0d0
        endif
   50 continue
c
      call popcha (model)
      call photsou (model)
      call copypop (pop, pop0)
c
c
   60 format(//' Give initial conditions:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    t (K), dh (N), dr (cm), Geo. dilution, Case A-B:'/
     &'    (t<10 taken as a log, dh<=0 taken as a log),'/
     &'    (dr<100 taken as log, Dilution factor <= 1.0)'/
     &'    (Case A-B (H) taken as 0(A) <-> 1 (B), <0 auto)',//
     &' :: ',$)
   70 write (*,60)
      read (*,*) t,dh,dr,wdil,cab
c
      if (t.le.10) t=10**t
c
      if (((wdil.gt.1.0d0).or.(wdil.lt.0.0d0)).or.(dh.le.0.0d0)) goto
     &70
c
      if ((t.le.0.0).and.((ilgg.lt.'C').or.(ilgg.eq.'E'))) goto 70
      tf=t
c
      de=dh*zion(1)
      if ((ilgg.eq.'A').and.(ndee.eq.'N')) de=feldens(dh,pop)
      de=feldens(dh,pop)
c
      tstep1=0.0d0
      tstep2=0.0d0
      tstep3=0.0d0
c
   80 write (*,90)
   90 format(//' Spectrum printout required? (y/n) ',$)
      read (*,30) jprin
c
      jprin=jprin(1:1)
      if (jprin.eq.'y') jprin='Y'
      if (jprin.eq.'n') jprin='N'
c
      if ((jprin.ne.'N').and.(jprin.ne.'Y')) goto 80
c
      if (jprin.eq.'Y') then
c
        call zerbuf
c
        if (ilgg.eq.'F') then
  100     write (*,110)
  110    format(//' Spectrum for each iteration? (y/n) ',$)
          read (*,30) jspec
c
          jspec=jspec(1:1)
          if (jspec.eq.'y') jspec='Y'
          if (jspec.eq.'n') jspec='N'
c
          if ((jspec.ne.'N').and.(jspec.ne.'Y')) goto 100
c
          if (jspec.eq.'Y') jspec='YES'
          if (jspec.eq.'N') jspec='NO'
c
        endif
c
      endif
c
      lmod='ALL'
      imod='ALL'
c
c     get runname
c
  120 format (a80)
  130 format(//' Give a name/code for this run: ',$)
      write (*,130)
      read (*,120) runname
c
      call copypop (pop0, pop)
c
      rad=1.d38
      if (wdil.ge.0.5d0) rad=0.d0
      if (wdil.gt.0.5d0) wdil=0.5d0
c
      dv=0.0d0
      fi=1.d0
c
      if (cab.ge.0) caseab(1)=cab
      if (cab.ge.0) caseab(2)=cab
c
c ilgg = 'A' fixed ionisation
c
c      write(*,*) 'Initial localem'
      call localem (t, de, dh)
      call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff (fi, dh)
c
      if (ilgg.eq.'B') then
c
c eq ionisation
c
        call equion (t, de, dh)
c
        i=0
c
        trea=1.d-12
        dift=0.d0
c
c
  140   call copypop (pop, popz)
c
        if (cab.ge.0) caseab(1)=cab
        if (cab.ge.0) caseab(2)=cab
c
c        write(*,*) 'Loop localem'
        call localem (t, de, dh)
        call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
        call zetaeff (fi, dh)
        call equion (t, de, dh)
c
        call difpop (pop, popz, trea, atypes, dift)
        i=i+1
c
        if ((dift.ge.1.d-4).or.(i.le.4)) goto 140
c
        call copypop (pop, popz)
c
        call cool (t, de, dh)
        call localem (t, de, dh)
        call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
        call zetaeff (fi, dh)
c
      endif
c
      trec=frectim(t,de,dh)
c
      if (ilgg.gt.'B') then
        if (ilgg.eq.'C') then
          nmod='EQUI'
        else
          if (ilgg.ne.'F') then
            nmod='TIM'
            write (*,150) zetae,qhdh
  150          format(/' Zetae :',1pg10.3,'  QHDH :',1pg10.3/
     &              ' Give time step (dt<100 as log) :',$ )
            read (*,*) tstep1
c
            if (tstep1.lt.100.d0) tstep=10**tstep
c
          else
c
            nmod='TIM'
  160       write (*,170) zetae,qhdh
  170          format(//' Zetae :',1pg10.3,'  QHDH :',1pg10.3/
     &              ' Give final time and source lifetime:',$ )
            read (*,*,err=160) tst,tsl
            if (tst.le.0.d0) goto 160
            fron=0.d0
            tstep3=dmax1(0.d0,tst-tsl)
            tstep1=fron*dmin1(tsl,tst)
            tstep2=dmax1(0.0d0,dmin1(tsl,tst)-tstep1)
c
  180       write (*,190)
  190 format(/' Density: Isochoric or Isobaric (c/b) :',$ )
            read (*,30) jden
c
            jden=jden(1:1)
            if (jden.eq.'c') jden='C'
            if (jden.eq.'b') jden='B'
c
            if ((jden.ne.'C').and.(jden.ne.'B')) goto 180
c
          endif
          if (ill.eq.'N') then
            call copypop (pop0, pop)
          else
            call copypop (pop, pop0)
          endif
        endif
        if (ilgg.lt.'E') then
          call teequi (t, tf, de, dh, tstep1, nmod)
          trec=frectim(tf,de,dh)
        else if (ilgg.eq.'E') then
c
          tst=0.0d0
          call cool (t, de, dh)
          write (*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
          pfx='sltd'
          np=4
          caller='TE'
          wmod='REAL'
c
          call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0,
     &     0.d0, 0.d0, dr, 0.d0, 0.d0, 1.d0, tphot)
c
          call timion (t, de, dh, xhyf, tstep1)
          tst=tst+tstep1
          call cool (t, de, dh)
          write (*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
          pfx='sltd'
          np=4
          caller='TE'
          wmod='REAL'
c
          call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0,
     &     0.d0, 0.d0, dr, 0.d0, 0.d0, 1.d0, tphot)
c
          dift=0.d0
c
  200     call copypop (pop, popz)
c
          if (cab.ge.0) caseab(1)=cab
          if (cab.ge.0) caseab(2)=cab
c
          call cool (t, de, dh)
          call localem (t, de, dh)
          call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
          call zetaeff (fi, dh)
          call timion (t, de, dh, xhyf, tstep1)
          call cool (t, de, dh)
          tst=tst+tstep1
          write (*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
          pfx='sltd'
          np=4
          caller='TE'
          wmod='REAL'
c
          call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0,
     &     0.d0, 0.d0, dr, 0.d0, 0.d0, 1.d0, tphot)
c
          call difpop (pop, popz, trea, atypes, dift)
c
          tstep1=tstep1*1.05
c
          i=i+1
c
          if ((dift.ge.0.001)) goto 200
c
          trec=frectim(t,de,dh)
        else
          wd00=0.0d0
          tex=100.d0
          exl=1.d-20
          if (tstep1.gt.0.d0) call teequi (t, tf, de, dh, tstep1, nmod)
          prescc=fpressu(t,dh,pop)
          call evoltem (tf, tf, de, dh, prescc, tstep2, 0.0d0, 0.0d0,
     &     luop)
          trec=frectim(tf,de,dh)
c
          if (cab.ge.0) caseab(1)=cab
          if (cab.ge.0) caseab(2)=cab
c
          call localem (t, de, dh)
          call totphot (t, dh, fi, rad, dr, dv, wd00, lmod)
          call evoltem (tf, tf, de, dh, prescc, tstep3, exl, tex, luop)
        endif
        call difpop (pop, popz, trea, atypes, dift)
        t=tf
      endif
c
c Sum up fluxes etc
c
      if (ilgg.lt.'E') then
        call sumdata (t, de, dh, fi, dr, dr, dr, imod)
      endif
      call avrdata
c
c Output..
c
      wmod='SCRN'
      call wmodel (luop, t, de, dh, dr, wmod)
c
  210  format(' Plasma Slab model diffuse field included.'/
     &' produced by MAPPINGS V ',a8,'     Run:',a80/)
      write (luop,210) theversion,runname
c
      wmod='FILE'
      call wmodel (luop, t, de, dh, dr, wmod)
c
      write (luop,220)
  220 format(/' ABUNDANCES OF THE ELEMENTS RELATIVE' ,
     &' TO HYDROGEN AND THEIR STATE OF IONISATION :' )
c
      call wionabal (luop, pop)
c
      write (luop,230) zetae,qhdh,tstep1,tstep2,tstep3,dift,trec
c
  230 format(/' ZETAE:',1pg10.3,4x,'QHDH:',1pg10.3,4x,'TSTEP:'
     &,3(1pg10.3),4x,'DIFT:',0pf7.3,4x,'TREC:',1pg9.2/)
      if (jprin.ne.'N') then
        write (luop,240) caseab(1),caseab(2)
  240 format(//' CASE A,B (H, He) : ',f4.2,f4.2)
c
        kmod='REL'
C        linemod='LAMB'
C        call spec2 (luop, linemod, kmod)
c
        call spectrum(luop, kmod)
c
      endif
c
      close (luop)
c
c     output diffuse field photon source file
c
      pfx='slso'
      np=4
      caller='TE'
      wmod='REAL'
c
      call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, 0.d0,
     &0.d0, dr, 0.d0, 0.d0, 1.d0, tphot)
c
      pfx='slfn'
      np=4
      caller='TE'
      wmod='NFNU'
c
      call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, 0.d0,
     &0.d0, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c      write out balance file
c
      pfx='slbal'
      np=5
      call wbal (caller, pfx, np, pop)
c
      write (*,250)
c
  250 format(//' Model ended, output created.' )
c
      return
      end
