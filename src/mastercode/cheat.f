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
c
c     Charge exchange heating.  Uses the energy defect from
c     charge exchange reactions.
c
c     chgain is a heating term so -ve is cooling
c     rates setup in allrates
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine cheat (dh)
c
      include 'cblocks.inc'
c
      real*8 dh, delta, ab
      real*8 chxheat,e4e
      integer*4 i,ion,at,rat
c
      chgain=0.d0
c
      if (chargemode.ge.2) return
      if (chargeheatmode.eq.0) return
c
      if (chargemode.eq.0) then
c
c     new reactions being used
c
c
        e4e=(rkb*1.0d4)
c
        do i=1,nchxr
          if ((chxr(i).gt.0.d0).and.(chxrde(i).ne.0.d0)) then
            rat=chxrx(i)
            at=chxrat(i)
            ion=chxrio(i)
            delta=chxrde(i)*e4e
            ab=dh*zion(at)*pop(ion,at)*dh*zion(rat)*pop(1,rat)
            if (ab.ge.pzlimit) then
              chxheat=chxr(i)*ab*delta
              chgain=chgain+chxheat
              heatz(at)=heatz(at)+chxheat
              heatzion(ion,at)=heatzion(ion,at)+chxheat
            endif
          endif
        enddo
        do i=1,nchxi
          if ((chxi(i).gt.0.d0).and.(chxide(i).ne.0.d0)) then
            at=chxiat(i)
            rat=chxix(i)
            ion=chxiio(i)
            delta=chxide(i)*e4e
            ab=dh*zion(at)*pop(ion,at)*dh*zion(rat)*pop(2,rat)
            if (ab.ge.pzlimit) then
              chxheat=chxi(i)*ab*delta
              chgain=chgain+chxheat
              heatz(at)=heatz(at)+chxheat
              heatzion(ion,at)=heatzion(ion,at)+chxheat
            endif
          endif
        enddo
      endif
c
      if (chargemode.eq.1) then
c
c     legacy A&R reactions being used
c
        do i=1,nlegacychxi
          if (chxilegacy(i).gt.0.d0) then
            at=chxilegacyat(i)
            rat=chxilegacyx(i)
            ion=chxilegacyio(i)
            delta=chxilegacycos(4,i)
            ab=dh*zion(at)*pop(ion,at)*dh*pop(2,rat)
            if (ab.ge.pzlimit) then
              chgain=chgain+chxilegacy(i)*ab*ev*delta
              heatz(at)=heatz(at)+chxilegacy(i)*ab*ev*delta
              heatzion(ion,at)=heatzion(ion,at)+chxilegacy(i)*ab*ev*
     &         delta
            endif
          endif
        enddo
c
      endif
c
      if (dabs(chgain).lt.epsilon) chgain=0.d0
c
      return
c
      end
