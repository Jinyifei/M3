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
c*******CALCULATES CHARGE EXCHANGE RATES WITH HYDROGEN
c     RETURNS RATES IN CHARTE(3,X) AND CHARTE(4,X)
c     ELEMENT NUMBER IS IN CHARTE(1,X)
c     AND THE LOWEST IONIC SPECIE FOR THAT ELEMENT (1<=N<=5)
c     IN CHARTE(2,X)   (EX:N <> N+ CORRESPONDS TO SPECIE 1 )
c     NB. IF CHARCO(1,X) IS PRESET NEGATIVE , REACTION IS DISABLED
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine charex (t)
c
      include 'cblocks.inc'
c
      real*8 t,tt,tl,tlim
      real*8 alph,sigma,tl0,tl1
      real*8 ktev,lowt
      real*8 x,t4,invt4
      integer*4 i,j,jn,nft,id
      real*8 a,b,c,d,e
      real*8 tlo4,thi4
c
c  internal functions
c
      real*8 fonct1, fonct2
      fonct1(x,a)=sigma*((x*1.0d-4)**a)
      fonct2(x,a)=sigma*dexp(-(a/x))
c
c
      do i=1,nchxr
        chxr(i)=0.d0
      enddo
      do i=1,nchxi
        chxi(i)=0.d0
      enddo
c
      if (chargemode.ge.3) return
c
      tt=dmax1(t,1.d0)
      tl=dlog10(tt)
      t4=tt*1.d-4
      invt4=1.d0/t4
c
      if (chargemode.eq.2) then
c
        tlim=6.0d4
        tl=dmin1(tt,tlim)
c
c     original calcs
c
        do 30 j=1,nchxold
          do 10 jn=1,5
   10       charte(jn,j)=0.0d0
c
          if (charco(1,j).le.0.0d0) goto 30
          charte(5,j)=charco(7,j)
          do 20 i=1,2
            charte(i,j)=charco(i,j)
            sigma=charco(2+i,j)
            alph=charco(4+i,j)
            if (alph.lt.20.0d0) charte(2+i,j)=fonct1(tl,alph)
            if (alph.ge.20.0d0) charte(2+i,j)=fonct2(tt,alph)
   20     continue
   30   continue
c
c     end old mappings exchange reactions
c
      endif
c
      if (chargemode.eq.1) then
c
        do i=1,nlegacychxr
          chxrlegacy(i)=0.d0
        enddo
        do i=1,nlegacychxi
          chxilegacy(i)=0.d0
        enddo
c
c     legacy A&R 1985 calcs,
c     returns rates in chxrlegacy and chxilegacy
c
        tt=t+1.d0
        t4=tt*1.d-4
c
c     first recombination rates
c
        do i=1,nlegacychxr
c
          tl0=dlog10(chxrlegacytemp(1,i))-0.5d0
          tl1=dlog10(chxrlegacytemp(2,i))+0.5d0
          tl=dlog10(t)
c
          if ((tl.gt.tl0).and.(tl.lt.tl1)) then
            chxrlegacy(i)=(chxrlegacycos(1,i)*(t4)**chxrlegacycos(2,i))*
     &       (1.d0+chxrlegacycos(3,i)*dexp(t4*chxrlegacycos(4,i)))
c
c Allow for hard fields by entering constant adiabatic limit
c
          else if (tl.lt.tl0) then
            lowt=10**(tl0-4)
            chxrlegacy(i)=(chxrlegacycos(1,i)*(lowt)**chxrlegacycos(2,i)
     &       )*(1.d0+chxrlegacycos(3,i)*dexp(lowt*chxrlegacycos(4,i)))
c
          else
            chxrlegacy(i)=0.0d0
          endif
c
        enddo
c
c     then ionising rates
c
        do i=1,nlegacychxi
c
          tl0=dlog10(chxilegacytemp(1,i))-0.5d0
          tl1=dlog10(chxilegacytemp(2,i))+0.5d0
          tl=dlog10(tt)
c
          if ((tl.gt.tl0).and.(tl.lt.tl1)) then
c
            ktev=rkb*tt/ev
c
c     special case for OI-HII reaction
c
            if (chxilegacyat(i).eq.zmap(8)) then
              chxilegacy(i)=chxilegacycos(1,i)*dexp(-1.d0*
     &         chxilegacycos(4,i)/ktev)*(1.d0-(0.93d0*dexp(-1.d0*
     &         chxilegacycos(3,i)*t4)))
c
            else
c
c     the rest of the reactions
c
              chxilegacy(i)=(chxilegacycos(1,i)*((t4)**chxilegacycos(2,
     &         i)))*dexp(-1.d0*chxilegacycos(3,i)*t4)*dexp(-1.d0*
     &         chxilegacycos(4,i)/ktev)
c
            endif
          else
            chxilegacy(i)=0.0d0
          endif
c
        enddo
c
      endif
c
      if (chargemode.eq.0) then
cc
cc      testloop logT = 0-10
cc
c      do j = 1,22
c      do idx = 0,200
c      tt   = 10.d0**(0.d0+idx*0.05d0)
cc
c
c     new calcs, returns rates in chxr and chxi
c
c     first recombination rates
c
        do i=1,nchxr
          chxr(i)=0.d0
c
          nft=chxrnfit(i)
          id=chxrfitid(i)
c
          a=chxrcos(1,i)
          b=chxrcos(2,i)
          c=chxrcos(3,i)
          d=chxrcos(4,i)
          e=0.0d0
c
          tl0=chxrtemp(1,i)
          tl1=chxrtemp(2,i)
c
          if ((tl.ge.tl0).and.(tl.lt.tl1)) then
            chxr(i)=a*(t4**b)*(1.d0+c*dexp(d*t4))
          endif
c
          if ((tl.lt.tl0).and.(id.eq.1)) then
            tlo4=10.d0**(tl0-4.d0)
            chxr(i)=a*(tlo4**b)*(1.d0+c*dexp(d*tlo4))
          endif
          if ((tl.ge.tl1).and.(id.eq.nft)) then
            thi4=10.d0**(tl1-4.d0)
            chxr(i)=a*(thi4**b)*(1.d0+c*dexp(d*thi4))
            chxr(i)=chxr(i)*(t4/thi4)**(-0.5d0)
          endif
c
          chxr(i)=dmax1(chxr(i),0.d0)
          if (chxr(i).le.1.d-17) chxr(i)=0.d0
c
ccc
cc test loop output
ccc
c      if (i .eq. j) then
c
c         if (idx .eq. 0) then
cc
c            if (elem_len(chxrat(i)).eq.2) then
c            write(*,'(a2,"^+",i1,x,$)')
c     &      elem(chxrat(i)), chxrio(i)-1
c            endif
c            if (elem_len(chxrat(i)).eq.1) then
c            write(*,'(a1,"^+",i1,x,$)')
c     &      elem(chxrat(i)), chxrio(i)-1
c            endif
cc
c         endif
cc
c            write(*,'(x,1pe9.2,$)')
c     &      chxr(i)*1e9
cc
c      endif
ccc
cc end test loop output
ccc
c
        enddo
c
c     then ionising rates
c
        do i=1,nchxi
c
          chxi(i)=0.d0
c
          nft=chxinfit(i)
          id=chxifitid(i)
c
          a=chxicos(1,i)
          b=chxicos(2,i)
          c=chxicos(3,i)
          d=chxicos(4,i)
          e=dmax1(0.d0,chxicos(5,i))
c
          tl0=chxitemp(1,i)
          tl1=chxitemp(2,i)
c
          if (chxix(i).eq.zmap(1)) then
            if ((tl.ge.tl0).and.(tl.lt.tl1)) then
              chxi(i)=a*(t4**b)*(1.d0+c*dexp(d*t4))*dexp(-e*invt4)
            endif
          endif
c
          if (chxix(i).eq.zmap(2)) then
            if ((tl.ge.tl0).and.(tl.lt.tl1)) then
c A&R1985 formula note -c, and we use E4 not eV
              chxi(i)=a*(t4**b)*dexp(-c*t4)*dexp(-e*invt4)
            endif
          endif
c
          if (chxix(i).eq.zmap(1)) then
            if ((tl.lt.tl0).and.(id.eq.1)) then
              tlo4=10.d0**(tl0-4.d0)
              chxi(i)=a*(tlo4**b)*(1.d0+c*dexp(tlo4))*dexp(-e/tlo4)
            endif
            if ((tl.ge.tl1).and.(id.eq.nft)) then
              thi4=10.d0**(tl1-4.d0)
              chxi(i)=a*(thi4**b)*(1.d0+c*dexp(d*thi4))*dexp(-e/thi4)
              chxi(i)=chxi(i)*(t4/thi4)**(-0.5d0)
            endif
          endif
c
          if (chxix(i).eq.zmap(2)) then
            if ((tl.lt.tl0).and.(id.eq.1)) then
              tlo4=10.d0**(tl0-4.d0)
              chxi(i)=a*(tlo4**b)*dexp(-c*tlo4)*dexp(-e/tlo4)
            endif
            if ((tl.ge.tl1).and.(id.eq.nft)) then
              thi4=10.d0**(tl1-4.d0)
              chxi(i)=a*(thi4**b)*dexp(-c*thi4)*dexp(-e/thi4)
              chxi(i)=chxi(i)*(t4/thi4)**(-0.5d0)
            endif
          endif
c
          chxi(i)=dmax1(chxi(i),0.d0)
          if (chxi(i).le.1.d-17) chxi(i)=0.d0
c
        enddo
cc
cc  testloops
cc
c      enddo
c      write(*,*)
c      enddo
c      stop
ccc
cc  end testloops
ccc
      endif
c
      return
      end
