cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO FIND INTERPOLATION VALUES USING NEIGHBORING POINTS TO
c     POINT*INL IN VECTOR BUFPHO .
c     BOUNDARIES OF BINS IN VECTOR EPHOT (EV) .
c     ASSUMES POWER LAW FOR INTERPOLATION :
c      DLOG(Y)=CFLO+a*DLOG(NU)
c
c
      subroutine interpol (inl, cflo, a, bufpho)
c
      include 'cblocks.inc'
c
c
      real*8 xx(5), yy(5), al(5), cc(5), bufpho(mxinfph)
      real*8 cflo,a, he, hemi, yc, yg, yd
      integer*4 inl, n, n1, modd, modg, k, ini, ifi
c
c
c   ***FIND INTERPOL. VALUES IN THE INTERVAL EPHOT(INL),EPHOT(INL+1)
c     FOR THE MEAN INTENSITY IN BUFPHO(INL) ASSUMING A POWER LAW
c     FUNCTION OF ENERGY .
c
      he=bufpho(inl)
      k=1
      al(1)=0.0d0
      al(2)=0.0d0
      al(3)=0.0d0
      al(4)=0.0d0
      al(5)=0.0d0
      cc(1)=-1.d20
      cc(2)=-1.d20
      cc(3)=-1.d20
      cc(4)=-1.d20
      cc(5)=-1.d20
      if (he.le.0.0d0) goto 50
      k=5
c
      cc(5)=dlog(he)
      if (inl.gt.2) then
        ini=1
        modg=1
      else
        ini=3
        modg=3
      endif
      if (inl.lt.(infph-2)) then
        ifi=5
        modd=4
      else
        ifi=3
        modd=2
      endif
c
      if ((inl.le.1)) goto 50
      if (((ifi-ini).ge.4).and.((bufpho(inl-1)+bufpho(inl+1)).lt.(he/
     &7.0d0))) goto 50
c
      do 10 n=ini,ifi
        n1=(inl+n)-3
        hemi=bufpho(n1)
        if ((hemi.gt.(he/200.0d0)).and.(hemi.lt.(he*200.0d0))) then
          yy(n)=dlog(hemi)
          xx(n)=dlog(cphotev(n1))
        else
          if (n.lt.3) modg=3
          if (n.gt.3) modd=2
        endif
   10 continue
c
      if (modg.gt.modd) goto 50
      do 20 n=modg,modd
        al(n)=(yy(n+1)-yy(n))/(xx(n+1)-xx(n))
        cc(n)=((yy(n)*xx(n+1))-(yy(n+1)*xx(n)))/(xx(n+1)-xx(n))
   20 continue
c
      if ((modg.eq.3).or.(modd.eq.2)) goto 30
      yc=cc(5)
      yg=cc(1)+(al(1)*xx(3))
      yd=cc(4)+(al(4)*xx(3))
      if (dabs(yg-yc).gt.dabs(yd-yc)) then
        modg=3
      else
        modd=2
      endif
c
   30 if (modg.eq.3) goto 40
      k=2
      goto 50
   40 if (modd.eq.2) goto 50
      k=3
   50 continue
c
c     CONVERTS INTO FUNCTION OF FREQUENCY INSTEAD OF 2XENERGY (EV)
c
      a=al(k)
c
      cflo=cc(k)-(33.1191d0*a)
c
      return
      end
