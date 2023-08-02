cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag (jl, alph, x)
c
      include 'const.inc'
c
      real*8 aa,aaa,dx
      real*8 ra,rab
      real*8 alph(7, 6), b(7, 6), a(7, 6), x(6), r(6)
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      integer*4 jl
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag3 (jl, alph, x)
c
      include 'const.inc'
c
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 alph(5, 4), b(5, 4), a(5, 4), x(4), r(4)
      integer*4 jl
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag6 (alph, x)
c
      include 'const.inc'
c
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 alph(8, 7), b(8, 7), a(8, 7), x(7), r(7)
      integer*4 jl
c
      jl=6
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag9 (alph, x)
c
      include 'const.inc'
c
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 alph(11, 10), b(11, 10), a(11, 10), x(10), r(10)
      integer*4 jl
c
      jl=9
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c
c
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mdiag16 (alph, x)
c
      include 'const.inc'
c
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 alph(18,17),b(18,17),a(18,17),x(17),r(17)
      integer*4 jl
c
      jl=16
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiagn (alph, x, n)
c
      include 'const.inc'
c
      real*8 alph((mxnl+2), (mxnl+1)), x(mxnl+1)
      integer*4 n
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 b((mxnl+2), (mxnl+1))
      real*8 a((mxnl+2), (mxnl+1)), r(mxnl+1)
      integer*4 jl
c
      jl=n
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c     BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c     BACK SUBSTITUTION
c
c     JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c     ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c     ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c     R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c     X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiagfe (alph, x, n)
c
      include 'const.inc'
c
      real*8 alph((mxfenl+2), (mxfenl+1)), x(mxfenl+1)
      integer*4 n
      real*8 aa,aaa,dx
      real*8 ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      real*8 b((mxfenl+2), (mxfenl+1))
      real*8 a((mxfenl+2), (mxfenl+1)), r(mxfenl+1)
      integer*4 jl
c
      jl=n
c
      jc=jl+1
      jca=jc+1
      do 20 il=1,jl
        do 10 ic=1,jca
   10     a(ic,il)=alph(ic,il)
   20 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa=0.0d0
      do 30 il=1,jl
        aa=aa+dabs(a(il,il))
        x(il)=0.0d0
   30 continue
      do 50 il=1,jl
        do 40 ic=1,jc
          a(ic,il)=a(ic,il)/aa
          b(ic,il)=a(ic,il)
   40   continue
        a(jca,il)=a(jc,il)
        b(jca,il)=a(jc,il)
   50 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 110 il=1,jl
        do 70 kl=il,jl
          rab=a(il,kl)
          ra=dabs(rab)
          if (ra.lt.epsilon) goto 70
          do 60 kc=1,jc
            a(kc,kl)=a(kc,kl)/rab
   60     continue
   70   continue
c
        do 80 kl=il,jl
          r(kl)=a(il,kl)
   80   continue
        do 100 kl=1,jl
          if (kl.eq.il) goto 100
          ra=dabs(r(kl))
          if (ra.lt.epsilon) goto 100
          do 90 kc=1,jc
   90       a(kc,kl)=a(kc,kl)-a(kc,il)
  100   continue
  110 continue
c
c    ***BACK SUBSTITUTION
c
      do 140 ila=1,jl
        il=(jl+1)-ila
        ics=il+1
        aaa=0.0d0
        if (ics.gt.jl) goto 130
        do 120 ic=ics,jl
  120     aaa=aaa+(a(ic,il)*x(ic))
  130   dx=(a(jc,il)-aaa)/a(il,il)
        x(il)=x(il)+dx
  140 continue
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
c
      subroutine matsolve (alph, x, n, np)
      include 'const.inc'
c
c general call for n levels in arrays of up to MAXN
c
      integer*4 maxn
      parameter (maxn=100)
c
      real*8 alph((np+2), (np+1)), x(np+1)
      integer*4 n,np
c
      integer*4 i,j,indx(maxn)
      real*8 a(maxn,maxn),b(maxn),d
c
      if (n.gt.maxn) then
        write (*,*) 'ERROR in matsolve, n exceed current MAXN'
        write (*,*) n,maxn
        write (*,*) 'edit matsolve to increase the maximum'
        stop
      endif
c
      do i=1,n
        do j=1,n
          a(i,j)=alph(j,i)
        enddo
      enddo
      do j=1,n
        b(j)=alph(n+1,j)
      enddo
c
      call ludcmp (a, n, maxn, indx, d)
      call lubksb (a, n, maxn, indx, b)
c
      do i=1,n
        x(i)=b(i)
      enddo
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
c
      subroutine matsolvemulti (alph, x, n)
      include 'const.inc'
c
c for nlevels up to mxnl, as used by multidat
c
      real*8 alph((mxnl+2), (mxnl+1)), x(mxnl+1)
      integer*4 n,np
      integer*4 i,j,indx(mxnl)
      real*8 a(mxnl,mxnl),b(mxnl),d
c
      np=mxnl
      do i=1,n
        do j=1,n
          a(i,j)=alph(j,i)
        enddo
      enddo
      do j=1,n
        b(j)=alph(n+1,j)
      enddo
c
      call ludcmp (a, n, np, indx, d)
      call lubksb (a, n, np, indx, b)
c
      do i=1,n
        x(i)=b(i)
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matsolvefe (alph, x, n)
      include 'const.inc'
c
c for nlevels up to mxfenl, as used by fedat
c with additional interative step to help large matrices
c stay stable
c
      real*8 alph((mxfenl+2), (mxfenl+1)), x(mxfenl+1)
      integer*4 n,np
      integer*4 i,j,indx(mxfenl)
      real*8 alu(mxfenl,mxfenl),b(mxfenl),d
      real*8 a(mxfenl,mxfenl)
c
      np=mxfenl
      do i=1,n
        do j=1,n
          alu(i,j)=alph(j,i)
        enddo
      enddo
      do j=1,n
        b(j)=alph(n+1,j)
      enddo
c
      call ludcmp (alu, n, np, indx, d)
      call lubksb (alu, n, np, indx, b)
c save x
      do i=1,n
        x(i)=b(i)
      enddo
c redo a and b
      do i=1,n
        do j=1,n
          a(i,j)=alph(j,i)
        enddo
      enddo
      do j=1,n
        b(j)=alph(n+1,j)
      enddo
c iterate solution, result directly in x
      call mprove (a, alu, n, np, indx, b, x)
c
c     do i=1,n
c     write(*,*) i, x(i)
c     enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NR fortran routines for LU decomp and back sub
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ludcmp (a, n, np, indx, d)
      include 'const.inc'
      integer*4  n,np,indx(n),nmax
      real*8 d,a(np,np),tiny
      parameter (nmax=500,tiny=1.0d-20)
      integer*4  i,imax,j,k
      real*8 aamax,dum,sum,vv(nmax)
      d=1.d0
      imax=0
      do i=1,n
        aamax=0.d0
        do j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        enddo
        if (aamax.lt.epsilon) then
          write (*,*) 'ERROR: singular matrix in ludcmp, in mdiag.f'
          stop
        endif
        vv(i)=1.d0/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j).eq.0.d0) a(j,j)=tiny
        if (j.ne.n) then
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lubksb (a, n, np, indx, b)
      implicit none
      integer*4 n,np,indx(n)
      real*8 a(np,np),b(n)
      integer*4 i,ii,j,ll
      real*8 sum
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        elseif (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mprove (a, alud, n, np, indx, b, x)
      implicit none
      integer*4 n,np,indx(n),nmax
      real*8 a(np,np),alud(np,np),b(n),x(n)
      parameter (nmax=500)
cu    uses lubksb
      integer*4 i,j
      real*8 r(nmax)
      real*8 sdp
      do i=1,n
        sdp=-b(i)
        do j=1,n
          sdp=sdp+a(i,j)*x(j)
        enddo
        r(i)=sdp
      enddo
      call lubksb (alud, n, np, indx, r)
      do i=1,n
        x(i)=x(i)-r(i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
