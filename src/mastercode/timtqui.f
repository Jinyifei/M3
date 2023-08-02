cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******FINDS THE EQUILIBRIUM TEMPERATURE AND THE
c     CORRESPONDING IONISING STATE OF THE GAS AFTER TIME TSTEP
c     TEI : INITIAL GUESS FOR TEMPERATURE
c     TEF,EDENS: FINAL TEMPERATURE AND ELECTRONIC DENSITY
c     CALL SUBR. COOL,TIMION
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine timtqui (tei, tef, edens, hdens, tstep)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a,b,c1,c2,cc,cof,csub,delmin
      real*8 dimax,dinc,dl1,emul,ex
      real*8 t2,tl1,tl2,tm1,tmin
      real*8 tsub,tti,u,xhyf
      real*8 popp0(mxion, mxelem), tei, tef, edens, hdens
     &                 , tstep
c
      integer*4 iel,ion,km
      integer*4 luop,n,nf
c
      character ntest*4
c
      real*8 arctanh
c
      arctanh(u)=dlog((1.d0+u)/(1.d0-u))/2.d0
      ntest='NO'
      if (tef.lt.0.0d0) ntest='Y'
      luop=23
      delmin=1.d-3
      nf=12
      dimax=1.4d0
      ex=1.0d0
      tmin=500.d0
      emul=4.0d0
      cof=0.8d0
c
c    ***USING 2 TEMP. (TM1,T2)  AND 2 FRACT. RESID. (DLOS1,DLOS2)
c     IT FITS THE FUNCTION : ARCTANH(DLOS)=B+A*LN(T)
c     NEXT VALUE OF T=DEXP(-B/A)
c
      km=0
      te0=tei
      if (te0.le.0.d0) te0=2.d4
      if (te0.lt.(3.d0*tmin)) te0=3.d0*tmin
      tm1=te0
      tl1=dlog(tm1)
      do 20 iel=1,atypes
        do 10 ion=1,maxion(iel)
          popp0(ion,iel)=pop(ion,iel)
   10   continue
   20 continue
c
      call timion (tm1, edens, hdens, xhyf, tstep)
c
      call cool (tm1, edens, hdens)
c     if (dabs(dlos).lt.(delmin/4.d0)) goto 200
      if (dabs(dabs(dlos)-1.d0).lt.1.d-6) dlos=(1.d0-1.d-6)*dsign(1.0d0,
     &dlos)
      c1=arctanh(dlos)
      dinc=dimax*(dabs(dlos)**ex)
      tl2=tl1-(dinc*dsign(1.d0,dlos))
      t2=dexp(tl2)
      if (t2.lt.tmin) t2=tmin
      tl2=dlog(t2)
      tti=t2
c
c
      dl1=dlos
      do 110 n=1,nf
        do 40 iel=1,atypes
          do 30 ion=1,maxion(iel)
            pop(ion,iel)=popp0(ion,iel)
   30     continue
   40   continue
c
        call timion (t2, edens, hdens, xhyf, tstep)
c
        call cool (t2, edens, hdens)
        if (dabs(dabs(dlos)-1.d0).lt.1.d-6) dlos=(1.0d0-1.d-6)*
     &   dsign(1.d0,dlos)
        c2=arctanh(dlos)
        csub=c2-c1
        tsub=tl2-tl1
        if (dabs(tsub).lt.1.d-36) tsub=1.d-36*dsign(1.0d0,tsub)
        if (tsub.eq.0.d0) tsub=1.d-36*dsign(1.0d0,csub)
        a=csub/tsub
        b=c1-(a*tl1)
        if ((80.d0*dabs(a)).lt.dabs(b)) a=(dabs(b)/80.d0)*dsign(1.0d0,a)
        b=c1-(a*tl1)
        cc=c1
        c1=c2
        tl1=tl2
        tm1=t2
        tl2=-(b/a)
        t2=dexp(tl2)
        if ((t2.le.(tm1*emul)).and.(t2.ge.(tm1/emul))) goto 50
        if (t2.gt.(tm1*emul)) t2=tm1*emul
        if (t2.lt.(tm1/emul)) t2=tm1/emul
        km=km+1
        tl2=dlog(t2)
        emul=1.d0+(cof*(emul-1.d0))
   50   continue
c
c      *PRINT OUT WHEN TESTED BY SUBR. TESTI
c
        if (t2.lt.tmin) goto 130
        if (ntest.ne.'Y') goto 100
        if (n.ne.1) goto 80
        write (luop,60)
        write (*,60)
   60 format(' ',t4,'T1',t14,'DLOS',t27,'T2',t39,'A',t49,'B',t58,'C1'
     &,t68,'C2')
        write (*,70) te0,dl1,tti,dinc
        write (luop,70) te0,dl1,tti,dinc
   70 format(' ',1pg12.5,g10.3,g12.5,5x,'(DINC:',g10.3,')')
   80   write (*,90) tm1,dlos,t2,a,b,cc,c2,km
        write (luop,90) tm1,dlos,t2,a,b,cc,c2,km,emul
c
   90 format(' ',1pg12.5,g10.3,g12.5,4g10.3,i4,1pg10.3)
  100   if ((dabs(dlos).lt.(delmin*5.d0)).or.(km.ge.3)) goto 130
  110 continue
c
      write (*,120) dlos,t2
  120 format(' CONVERGENCE FOR EQUIL. TEMP. IS TOO SLOW   ' ,'DL:'
     &,1pg9.2,'  TE:',1pg9.3)
  130 continue
      if (t2.lt.tmin) t2=tmin
      tef=t2
      do 150 iel=1,atypes
        do 140 ion=1,maxion(iel)
          pop(ion,iel)=popp0(ion,iel)
  140   continue
  150 continue
c
      call timion (tef, edens, hdens, xhyf, tstep)
c
      call cool (tef, edens, hdens)
c
      return
      end
