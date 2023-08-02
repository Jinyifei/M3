cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE PHOTOIONISATION RATES IN RPHOT(mxion, mxelem) ,
c     PHOTOHEATING RATES IN HEAPH(5,28,16) AND
c     BASIS FOR SECONDARY IONISATION IN ANR(2,28,16),WNR(2,28,16)
c
c     USES LOCAL MEAN INTENSITY OF RADIATION JNU CONTAINED
c     IN VECTOR : TPHOT     ( NUMBER OF CROSS SECTIONS : IONUM )
c     IT MULTIPLIES JNU BY 4PI AND INTEGRATES THE NUMBER OF PHOTONS
c
      subroutine phion ()
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 augen,augerg,augt,dain
      real*8 daiau,daip,dait,e0,e1,e2,eau2,eau2ev,ed,etr
      real*8 cedge,nedge,oedge,hiedge
      real*8 wid,q,rr1,rr2
c
      integer*4 i,ie,inl,j,n
c
c   **  EFF(A,B,ED)=B+1.0/(ED/A+1.0/(1.0-B+1.0D-12))
c   **  COEFFICIENTS A,B FROM SHULL,M.J. APJ.234,P761 (1979)
c
      real*8 eff1, eff2, eff3, eff4, eff5
c
c
      eff1(ed)=1.0d0
      eff2(ed)=0.708d0+(1.0d0/((ed*0.03225806451612903d0)+
     &3.424657534246576d0))
      eff3(ed)=0.383d0+(1.0d0/((ed*0.0617283950617284d0)+
     &1.620745542949757d0))
      eff4(ed)=0.194d0+(1.0d0/((ed*0.1388888888888889d0)+
     &1.240694789081886d0))
      eff5(ed)=0.113d0+(1.0d0/((ed*0.3076923076923077d0)+
     &1.127395715896280d0))
c
      cedge=2.80d+02
      nedge=3.95d+02
      oedge=5.33d+02
      hiedge=1.32d+03
c
      eau2=30.0d0
      eau2ev=eau2*ev
      qtosoh=0.0d0
c
      do i=1,atypes
        do j=1,maxion(i)
          rasec(j,i)=0.0d0
          auphot(j,i)=0.0d0
          rphot(j,i)=0.0d0
          do n=1,2
            anr(n,j,i)=0.0d0
            wnr(n,j,i)=0.0d0
            heaph(n,j,i)=0.0d0
          enddo
          do n=3,5
            heaph(n,j,i)=0.0d0
          enddo
        enddo
      enddo
cc
c in zero field mode there is no photoionsation
c
      if (photonmode.eq.0) return
c
c    If there is dust include atoms contained in dust above Auger limit
c
c
      do 20 inl=ionstartbin,infph-1
c
        if (tphot(inl).le.epsilon) goto 20
c
        e0=cphotev(inl)
        e1=photev(inl)
        e2=photev(inl+1)
        wid=(e2-e1)*evplk
        q=fpi*(wid*tphot(inl)/cphote(inl))
        if (e0.ge.iph) then
          qtosoh=qtosoh+q
        endif
c
c          if (e1.lt.1.d0) goto 30
c
        do 10 i=1,ionum
          if (inl.ge.photbinstart(i)) then
c
            ie=atpho(i)
            j=ionpho(i)
c          pz=zion(ie)*pop(j,ie)
cc
c          if (pz.gt.pzlimit) then
            augen=augpho(i)
c
            dain=q*photxsec(i,inl)
c
            if (grainmode.ne.0) then
              if ((mapz(ie).eq.6).and.(e2.gt.cedge)) then
                dain=dain*invdion(ie)
              else if ((mapz(ie).eq.7).and.(e2.gt.nedge)) then
                dain=dain*invdion(ie)
              else if ((mapz(ie).eq.8).and.(e2.gt.oedge)) then
                dain=dain*invdion(ie)
              else if ((mapz(ie).gt.8).and.(e2.gt.hiedge)) then
                dain=dain*invdion(ie)
              endif
            endif
c
            if (dain.lt.1.d-35) dain=0.d0
            if (dain.le.0.d0) goto 10
            if (augen.le.0.0d0) then
              rphot(j,ie)=rphot(j,ie)+dain
            else
              auphot(j,ie)=auphot(j,ie)+dain
            endif
c
            daip=dain*(cphotev(inl)-ipotpho(i))*ev
c
            if (daip.gt.0.0d0) then
              etr=ipotpho(i)
              rr2=e2/etr
              rr1=e1/etr
              ed=dmax1(1.d-3,(((rr1*etr)+e2)*0.5d0)-(etr+10.2d0))
              heaph(1,j,ie)=heaph(1,j,ie)+(eff1(ed)*daip)
              heaph(2,j,ie)=heaph(2,j,ie)+(eff2(ed)*daip)
              heaph(3,j,ie)=heaph(3,j,ie)+(eff3(ed)*daip)
              heaph(4,j,ie)=heaph(4,j,ie)+(eff4(ed)*daip)
              heaph(5,j,ie)=heaph(5,j,ie)+(eff5(ed)*daip)
              dait=daip-(dain*iphe)
              if (dait.ge.0.0d0) then
                anr(1,j,ie)=anr(1,j,ie)+dait
                wnr(1,j,ie)=wnr(1,j,ie)+dain
                dait=daip-(dain*eau2ev)
                if (dait.ge.0.0d0) then
                  anr(2,j,ie)=anr(2,j,ie)+dait
                  wnr(2,j,ie)=wnr(2,j,ie)+dain
                endif
              endif
            endif
c
            if (augen.gt.iph) then
              augt=augen-10.2d0
              augerg=augen*ev
              daiau=dain*augerg
              heaph(1,j,ie)=heaph(1,j,ie)+(eff1(augt)*daiau)
              heaph(2,j,ie)=heaph(2,j,ie)+(eff2(augt)*daiau)
              heaph(3,j,ie)=heaph(3,j,ie)+(eff3(augt)*daiau)
              heaph(4,j,ie)=heaph(4,j,ie)+(eff4(augt)*daiau)
              heaph(5,j,ie)=heaph(5,j,ie)+(eff5(augt)*daiau)
              dait=dain*(augerg-iphe)
              if (dait.ge.0.0d0) then
                anr(1,j,ie)=anr(1,j,ie)+dait
                wnr(1,j,ie)=wnr(1,j,ie)+dain
                dait=dain*(augerg-eau2ev)
                if (dait.ge.0.0d0) then
                  anr(2,j,ie)=anr(2,j,ie)+dait
                  wnr(2,j,ie)=wnr(2,j,ie)+dain
                endif
              endif
            endif
          endif
c          endif
   10   continue
c
   20 continue
c
c     add cosmic ray ionization rate, crate, to neutral hydrogen
c     and Helium
c     crate of order 10^-17, cosphi 0-0.75
c
c      if (crate.gt.0.d0) then
c      rphot(1,1)=rphot(1,1)+crate*(1.d0+cosphi)
c      rphot(1,2)=rphot(1,2)+crate*(1.d0+cosphi)
c      endif
c
      do i=1,atypes
        do j=1,maxion(i)
          if (rphot(j,i).lt.1.d-35) rphot(j,i)=0.d0
        enddo
      enddo
c
      return
      end
