cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******INTEGRATES ARRAY BUFPHO THAT CONTAINS THE LOCAL MEAN INTENS.
c   JNU TO DERIVE THE AVAILABLE NUMBER OF PHOTONS TO IONISE H,HE
c   (HENCE INTEGRATES 4PI*JNU/HNU)
c   NB.  BUFPHO MUST BE DEXPRESSED IN  ERGS.CM-2.SEC-1.HZ-1.STER-1 (1/4pi)
c   CALL SUBR. INTERPOL
c
c
      subroutine intvec (bufpho, qahi, qahei, qaheii, qatot)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 bufpho(mxinfph)
      real*8 qahi, qahei, qaheii, qatot
      real*8 phq(2, 2)
      real*8 wid,q
c
      integer*4 i, j, inl
c
      do i=1,2
        do j=1,2
          phq(j,i)=0.d0
        enddo
      enddo
c
      do 30 inl=1,infph-1
        if (bufpho(inl).le.epsilon) goto 30
        wid=(photev(inl+1)-photev(inl))*evplk
        q=fpi*(wid*bufpho(inl)/cphote(inl))
        do 20 i=1,2
          do 10 j=1,maxion(i)-1
            if (cphote(inl).lt.ipote(j,i)) goto 10
            phq(j,i)=phq(j,i)+q
   10     continue
   20   continue
   30 continue
c
      qaheii=phq(2,2)
      qahei=dmax1(0.d0,phq(1,2)-phq(2,2))
      qahi=dmax1(0.d0,phq(1,1)-phq(1,2))
c
      qatot=phq(1,1)
c
      return
      end
