cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO SOLVE THE DIFFERENIAL EQUATION THAT
c     DETERMINES THE IONIC ABUNDANCE OF HII OR HI
c     FR : INITIAL AND (OUTPUT) FINAL FRACTION
c
c
      subroutine sdifeq (a, b, c, fr, tstep)
c
      include 'const.inc'
c
      real*8 a, b, c, fr, frf, do, du,tstep
      real*8 dq, di, dt, ds, dm, d1, d2, dl, dr, tim, bb
c
      frf=fr
      tim=tstep
      bb=dabs(b)
      dm=1.d-38
      if (tim.le.0.d0) goto 80
      if (((a.eq.0.d0).and.(b.eq.0.d0)).and.(c.eq.0.d0)) goto 80
c
c*******CASE WHERE C<>0
      if (c.eq.0.d0) goto 70
      if ((b.eq.0.d0).and.(a.eq.0.d0)) goto 60
      di=(4.d0*c)*a
      ds=dlog(bb)
      dq=bb-(di/bb)
      if (dq.le.0.d0) goto 10
      dr=dlog(dq)
      dq=dexp((ds+dr)/2.d0)
      if (dq.gt.0.d0) goto 30
   10 if (((tim*dm)*dabs(a)).gt.1.d0) tim=1.d0/(dm*dabs(a))
      di=b/(2.d0*a)
      ds=fr+di
      if (dabs(ds).lt.dm) goto 20
      dr=(1.d0/ds)-(a*tim)
      frf=(1.d0/dr)-di
      goto 80
   20 dr=1.d0-((a*ds)*tim)
      frf=(ds/dr)-di
      goto 80
   30 continue
      if (((dabs(dq-bb).lt.(1.d-8*bb)).and.(dabs((2.d0*di)/b).lt.1.d-8))
     &.and.(dabs(c).lt.(1.d-4*bb))) goto 70
      if (((tim*dm)*dq).gt.1.d0) tim=1.d0/(dq*dm)
      di=(2.d0*fr)*c
      d1=di+(b-dq)
      d2=di+(b+dq)
      if (dabs(d1).lt.dm) d1=-dm
      if ((d2.lt.dm).and.(d2.ge.0.d0)) d2=dm
      if ((d2.gt.(-dm)).and.(d2.lt.0.d0)) d2=-dm
      ds=dsign(1.d0,d1)*dsign(1.d0,d2)
      if ((dabs(b-dq)+dabs(b+dq)).lt.(1.d-7*dabs(di))) goto 40
      du=dlog(dabs(d1))
      do=dlog(dabs(d2))
      goto 50
   40 continue
      du=(b-dq)/di
      do=(b+dq)/di
   50 continue
      dl=(du-do)+(dq*tim)
      dm=-dlog(dm)
      if (dl.gt.dm) dl=dm
      dt=ds*dexp(dl)
      du=1.d0+dt
      do=1.d0-dt
      if (dabs(do).lt.1.d-6) do=-(ds*dl)
      if (dabs(du).lt.1.d-6) du=-(ds*dl)
      dr=du/do
      frf=((dq*dr)-b)/(2.d0*c)
      goto 80
c
   60 if (((tim*dm)*dabs(c)).gt.1.d0) tim=1.d0/(dm*dabs(c))
      di=1.d0-((c*tim)*fr)
      frf=fr/di
      goto 80
   70 di=a+(b*fr)
      if ((di.lt.dm).and.(di.ge.0.d0)) di=dm
      if ((di.gt.(-dm)).and.(di.lt.0.d0)) di=-dm
      ds=dsign(1.d0,di)
      if (((tim*dm)*bb).gt.1.d0) tim=1.d0/(dm*bb)
      dl=(dlog(dabs(di)+dm)+(tim*b))-dlog(bb)
      dm=-dlog(dm)
      if (dl.gt.dm) dl=dm
      dt=(ds*dexp(dl))*dsign(1.d0,b)
      frf=dt-(a/b)
   80 if (frf.lt.0.d0) frf=0.d0
      if (frf.gt.1.d0) frf=1.d0
      fr=frf
c
      return
      end
