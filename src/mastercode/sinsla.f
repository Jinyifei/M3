cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       single slab ionisation balance routine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sinsla (model)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 cab,de,dh,dr,dv,epotmi,exl,fi,fron
      real*8 prescc,rad
      real*8 tex,tf,trea,trec,tsl,tst,tstep1,tstep2
      real*8 tstep3,wd00,xhyf,t
      real*8 popz(mxion, mxelem),dift
c
      integer*4 i,j,l,luop
c
      character imod*4, lmod*4,  nmod*4
      character ilgg*4, ill*4,model*32
c
c           Functions
c
      real*8 fpressu,frectim
c
c
      luop=23
      cab=0.0d0
      epotmi=iphe
c
c
      trea=epsilon
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      jspot='NO'
      jcon='YES'
c
c
   10 write (*,20)
   20 format(//' Choose a single slab model:'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     B :  Equilibrium ionisation at a fixed temp.'/
     &'     C :  Equilibrium ionisation and equilibrium temp.'/
     &'     D :  Time dependent ionisation at equilibrium temp.'/
     &'     E :  Time dependent ionisation at a fixed temp.'/
     &'     F :  Time dependent ionisation and temperature.'/
     &'     X :  Exit.'/
     &'  :: ',$ )
      read (*,30) ilgg
c
   30 format(a)
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
      if (ilgg.eq.'d') ilgg='D'
      if (ilgg.eq.'e') ilgg='E'
      if (ilgg.eq.'f') ilgg='F'
      if (ilgg.eq.'x') ilgg='X'
c
      l=1
c      if ((ilgg .eq. 'CC').or.(ilgg .eq. 'BB')) l = 2
      if (ilgg.eq.'X') goto 150
c
      if ((ilgg.lt.'B').or.(ilgg.gt.'F')) goto 10
c
      do 50 i=1,atypes
        do 40 j=1,maxion(i)
   40     pop(j,i)=0.0d0
c
        if ((ilgg.eq.'C').or.(ilgg.eq.'B')) pop(l,i)=1.0d0
        if (((ilgg.eq.'D').or.(ilgg.eq.'E')).or.(ilgg.eq.'F')) then
          l=1
          if (ipote(1,i).lt.epotmi) l=2
          pop(l,i)=1.0d0
        endif
   50 continue
c
c
      call photsou (model)
c
      call copypop (pop, pop0)
c
c
      if (ilgg.gt.'B') ill='N'
c
   60 write (*,70)
   70 format(/' Give plasma slab conditions:'/
     &' ( T (K) T<10 as log, dr in cm, Dilution <= 1.0 )'/
     &' : T, H density, slab thickness (dr), Geo. Dilution'/
     &' : ',$)
      read (*,*) t,dh,dr,wdil
      if (t.le.10) t=10**t
      if (((wdil.gt.1.0d0).or.(wdil.lt.0.0d0)).or.(dh.le.0.0d0)) goto
     &60
      if ((t.le.0.0).and.((ilgg.lt.'C').or.(ilgg.eq.'E'))) goto 60
      tf=t
      t=dabs(t)
      de=dh*zion(1)
c
      tstep1=0.0d0
      tstep2=0.0d0
      tstep3=0.0d0
c
      lmod='SO'
      imod='ALL'
      rad=1.d38
      dr=1.0d0
      dv=0.0d0
      fi=1.d0
c
c
      caseab(1)=cab
      caseab(2)=cab
c
      call localem (t, de, dh)
      call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff (fi, dh)
c
      if (ilgg.eq.'B') then
        call equion (t, de, dh)
c
        i=0
c
        trea=1.d-12
        dift=0.d0
c
c
   80   call copypop (pop, popz)
c
        call localem (t, de, dh)
        call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
        call zetaeff (fi, dh)
        call equion (t, de, dh)
c
        call difpop (pop, popz, trea, atypes, dift)
c
        call copypop (pop, popz)
c
        i=i+1
c
        if ((dift.ge.1.0d-4).or.(i.le.4)) goto 80
c
c         write (*,*) i,dift
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
   90       write (*,100) zetae,qhdh
  100          format(/' ZETAE :',1pg10.3,4x,'QHDH :',1pg10.3/
     &   'Give time step (s): ',$)
            read (*,*) tstep1
c
            if (tstep1.lt.0.d0) goto 90
          else
            nmod='TIM'
  110       write (*,120) zetae,qhdh
  120          format(/' ZETAE :',1pg10.3,4x,'QHDH :',1pg10.3/
     &   ' Give age of nebula and source life-time (s, s): ',$)
            read (*,*) tst,tsl
            if (tst.le.0.d0) goto 110
            fron=0.d0
            tstep3=dmax1(0.d0,tst-tsl)
            tstep1=fron*dmin1(tsl,tst)
            tstep2=dmax1(0.0d0,dmin1(tsl,tst)-tstep1)
c
  130       write (*,140)
  140    format(' Density behaviour, Isochoric or Isobaric (C/B): ',$)
            read (*,30) jden
            if ((jden.ne.'C').and.(jden.ne.'B')) goto 130
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
          call timion (t, de, dh, xhyf, tstep1)
          call cool (t, de, dh)
          trec=frectim(t,de,dh)
        else
          write (*,*) ' Time Dependent Cooling and Recombination'
          wd00=0.0d0
          tex=100.d0
          exl=1.d-20
          if (tstep1.gt.0.d0) call teequi (t, tf, de, dh, tstep1, nmod)
          prescc=fpressu(t,dh,pop)
          call evoltem (tf, tf, de, dh, prescc, tstep2, 0.0d0, 0.0d0,
     &     luop)
          trec=frectim(tf,de,dh)
          call localem (t, de, dh)
          call totphot (t, dh, fi, rad, dr, dv, wd00, lmod)
          call evoltem (tf, tf, de, dh, prescc, tstep3, exl, tex, luop)
        endif
        call difpop (pop, pop0, trea, atypes, dift)
        t=tf
      endif
c
  150 continue
c
c     new population is now in pop
c
      write (*,*) ' T :',t
      write (*,*) ' dh:',dh
      write (*,*) ' de :',de
c
      return
      end
