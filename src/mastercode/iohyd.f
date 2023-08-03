cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine iohyd (dh, xhy, t, tstep, def, xhyf, mod)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 co,rc,phi,rkn,eps,abu
      real*8 dh, xhy, t, tstep, def, xhyf
      real*8 a1,a2,aa,b1,b2,c1,c2
      real*8 de,de2,del,dep,dey,fhi
      real*8 fhii,ph,ph2,spo,yh,dehyd
      real*8 maxrate, eqhi, eqhii, rateratio
      real*8 rectotal, iontotal
      real*8 rce, coe, phe
c
c     Old charge exchange
c      real*8 chex
c
      integer*4 at,ion,j,iel,ies
c
      character jjmod*4, mod*4
c
c           Functions
c
      real*8 feldens
c
      if (tstep.ge.0.0d0) goto 20
      write (*,10)
   10 format(' STOP: NEGATIVE TIME STEP'//)
      stop
   20 continue
c
      if ((jspot.eq.'YES').or.(jspot.eq.'NO')) goto 40
      write (*,30) jspot
   30 format(//' FLAG "JSPOT" IMPROPERLY DEFINED : ',a4)
      stop
   40 continue
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod='ALL'
c
c    ***INITIAL CONDITIONS
c
      call allrates (t, jjmod)
c
      fhii=pop(2,1)
      fhi=pop(1,1)
      if ((xhy.lt.0.0d0).or.(xhy.gt.1.0d0)) goto 50
      fhii=xhy
      fhi=1.d0-xhy
   50 continue
      dehyd=feldens(dh,pop)
c
c     **COL- , PHOTO- IONISATION AND RECOMBINATION RATES
c
      del=dmax1(0.d0,dehyd-((dh*zion(1))*fhii))/dh
      co=dh*col(1,1)
      rc=dh*rec(2,1)
      phi=rphot(1,1)
c
c    **IONISATION RATE CORRECTION DUE TO THE ON THE SPOT APPROX.
c
      spo=0.0d0
c
      if (jspot.ne.'YES') goto 60
c
      call spotap (de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
      rc=dh*rec(3,1)
      aa=(dh*zion(2))*pop(2,2)
      spo=aa*((yh*(rec(2,2)-rec(4,2)))+(ph*rec(4,2)))
      aa=(dh*zion(2))*pop(3,2)
      spo=spo+((aa*ph2)*rec(5,2))
   60 continue
c
c     **IONISATION RATE CORRECTION DUE TO NON-THERMAL ELECTRONS
c
      call ionsec
      phi=phi+rasec(1,1)
c
c    ***ADD CHARGE EXCHANGE RATES BETWEEN HYDROGEN AND HEAVIER ELEMENTS
c
      rkn=0.0d0
      eps=0.0d0
c
      if (chargemode.eq.0) then
c
c     Recomb reactions
c     with neutral H -> ionise H
c
        do j=1,nchxr
          if (chxr(j).gt.0.d0) then
            at=chxrat(j)
            ies=chxrx(j)
            if (ies.eq.zmap(1)) then
              ion=chxrio(j)
              abu=zion(at)*pop(ion,at)
              eps=eps+(chxr(j)*abu*dh)
            endif
          endif
        enddo
c
c     Ionising reactions
c     with ionised H -> recomb H
c
        do j=1,nchxi
          if (chxi(j).gt.0.d0) then
            at=chxiat(j)
            ies=chxix(j)
            if (ies.eq.zmap(1)) then
              ion=chxiio(j)
              abu=zion(at)*pop(ion,at)
              rkn=rkn+(chxi(j)*abu*dh)
            endif
          endif
        enddo
      endif
c
      if (chargemode.eq.1) then
c
c     old mappings charge exchanges
c
        do 70 j=1,nchxold
          iel=idint(charte(1,j))
          if (iel.lt.2) goto 70
          ion=idint(charte(2,j))
          if (rec(ion+1,iel).le.0.0d0) goto 70
          ies=idint(charte(5,j))
          if (ies.ne.1) goto 70
          abu=dh*zion(iel)
          if (abu*pop(ion,iel).lt.pzlimit) goto 70
          if (abu*pop(ion+1,iel).lt.pzlimit) goto 70
          rkn=rkn+((charte(3,j)*abu)*pop(ion,iel))
          eps=eps+((charte(4,j)*abu)*pop(ion+1,iel))
   70   continue
c
c     end old charge reactions
c
      endif
c
      coe=dehyd*col(1,1)
      rce=dehyd*rec(2,1)
      phe=phi
c
c harmonic dampening of CHX rates for stability in cold dark gas
c
      if (eps.gt.0.d0) then
        maxrate=dmax1(co,phi,spo)
        eps=1.d0/((2.d0/maxrate)+(1.d0/eps))
      endif
      if (rkn.gt.0.d0) then
        rkn=1.d0/((2.d0/rce)+(1.d0/rkn))
      endif
c
      if (mod.eq.'EQUI') then
        rectotal=dmax1((rce+rkn),0.d0)
        iontotal=dmax1((coe+phe-spo+eps),0.d0)
c      write(*,*) 'iohyd:  ', rectotal,iontotal,rc,rkn,coe,phe,eps
        if (iontotal.gt.rectotal) then
          rateratio=rectotal/(iontotal+epsilon)
          eqhi=rateratio/(rateratio+1.d0)
          eqhi=dmax1(eqhi,0.d0)
          eqhii=1.d0-eqhi
        else
          rateratio=iontotal/(rectotal+epsilon)
          eqhii=rateratio/(rateratio+1.d0)
          eqhii=dmax1(eqhii,0.d0)
          eqhi=1.d0-eqhii
        endif
        fhi=eqhi
        fhii=eqhii
c      write(*,*) 'iohyd: Simple EQ ', fhi, fhii
      endif
c
      a1=(del*(co+spo))+phi+eps
      b1=((1.d0-del)*co)-phi+spo-(del*rc)-(eps+rkn)
      c1=-(rc+co)
c
      a2=(1.d0+del)*(rc-spo)+rkn
      b2=(((-(co*(del+1.d0)))-phi)+spo)-rc*(del+2.d0)-(eps+rkn)
      c2=-c1
c
c      if (mod.ne.'EQUI') then
c
      call sdifeq (a1, b1, c1, fhii, tstep)
      call sdifeq (a2, b2, c2, fhi, tstep)
c
c      endif
c
      if (fhi.lt.pzlimit) then
        fhi=0.d0
        fhii=1.d0
      endif
c
      if (fhii.lt.pzlimit) then
        fhii=0.d0
        fhi=1.d0
      endif
c
c c
c c    ***WARNING MESSAGE AND NORMALISATION IF ABUNDANCES
c c   DO NOT CONSERVE
c c
c       da = fhi+fhii
c       dma = 1.d2*dabs(1.d0-da)
c       if (dma.ge.1.d0) then
c        write(*, 940) dma
c   940  format(' WARNING:'/
c      & ' ABUND. OF H DID NOT CONSERVE BEFORE NORMALISATION :',g9.3,2h %
c      & //)
c        write(*,*) rc,co, phi
c        write(*,*) rphot(1,1),rasec(1,1)
c        write(*,*) del, dehyd, dh
c        write(*,*) pop(2,1),pop(1,1)
c        fhi = fhi/da
c        fhii = fhii/da
c       endif
c
c    ***SET OUTPUT VARIABLES AND DETERMINE TIME DERIVATIVE
c
      pop(1,1)=fhi
      pop(2,1)=fhii
      xhyf=fhii
      def=feldens(dh,pop)
      dndt(1,1)=dh*((((c2*fhi)+b2)*fhi)+a2)
c
      dndt(2,1)=dh*((((c1*fhii)+b1)*fhii)+a1)
      return
      end
