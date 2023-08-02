cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******OUTPUT EMISSION SPECTRUM ON DEVICE*LUOP
c     USES SUMMED LINE DATA IN /SLINE/
c     ASSUMED INTENSITY IS RELATIVE TO H-BETA
c     WHEN MODE='REL'
c
c     Internal Ritz (vac) wavelengths converted to std Air on output.
c
c
      subroutine spectrum (luop, mode)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 chcksum,chklim
c
      integer*4 i,j,jl,k,max,atom,ion
      integer*4 luop,ionidx,lj, line, series
      integer*4 nl, nr, nd, ns, ne,nt
      integer*4 idx,jdx,nlo,nhi,trdx
      integer*4 l0,l1
      real*8 lineflx(5)
      real*8 linelam(5)
      integer*4 linelo(5),linehi(5)
      character*24 lineid(5)
      character*16 lineion(5)
      character*32 pad,ionstr
c
      character mode*4
c
c function
c
      real*8 fnair
c
c     Formats:
c
   10 format(t6,  'Lambda (A)',t16,5(0pf18.3))
   20 format(t6,  'Lambda (m)',t16,5(0pf18.3))
   30 format(t5, 'Int./H-beta',t16,5(1pe18.4))
   40 format(t4,'CaseA/H-beta',t16,5(1pe18.4))
   50 format(t4,'CaseB/H-beta',t16,5(1pe18.4))
   60 format(t4,'CaseC/H-beta',t16,5(1pe18.4))
   70 format(t6,  '    Levels',t16,5(15x,i1,'-',i1))
   80 format(t6,  '    Levels',t16,5(11x,i3,'-',i3))
   90 format(t6,  'Transition',t16,5(9x,a9))
  100 format(t6,  'Transition',t16,5(8x,a10))
  110 format(t6,  'Transition',t16,5(2x,a16))
  120 format(t6,  'Transition',t16,5a18)
  130 format('')
  140 format('')
  150 format(t2,'Ion: ',a3,a6)
  160 format(t4,'         Ion',t16,5(9x,a9))
  170 format(t2,'Ion: ',a3,a6,x,i3,' levels')
c     chklim is a lower limit for line fluxes to print out.
c     1e-4 saves a lot of paper
c
      chklim=epsilon
      if (mode.eq.'REL') then
        chklim=1.d-5
      endif
c
c     Top Twenty Five Lines - ie rough spectrum, include HBeta printout now too.
c
      call spec2 (luop, 'TTWN', mode)
c
c     The full spectrum....
c
      write (luop,180)
  180 format(// ' Detailed Emission Spectrum: Hbeta = 1.000:'/
     &          ' ==============================================')
c
c      write (luop,190) h2qav,h2qflux(1)
      write (luop,190) h2qflux(1)
  190 format(//,' Hydrogen Spectrum, includes collisions:'/
     &          ' =============================================='/
     &          ' HI Two-Photon emission :',1pg10.3)
c
  200 format(/,' Lyman n = 1 Series'/
     &         ' ======================')
      write (luop,200)
c
      do i=1,nhlines
        if (hydroflux(i,1).lt.chklim) hydroflux(i,1)=0.d0
      enddo
c
      write (luop,150) elem(1),rom(1)
      write (luop,130)
      write (luop,10) (hlambda(j,1)/fnair(hlambda(j,1)),j=1,5)
      write (luop,30) (hydroflux(j,1),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,1)/fnair(hlambda(j,1)),j=6,10)
      write (luop,30) (hydroflux(j,1),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,1)/fnair(hlambda(j,1)),j=11,15)
      write (luop,30) (hydroflux(j,1),j=11,15)
c
  210 format(/,' Balmer n = 2 Series'/
     &         ' ======================')
      write (luop,210)
c
      write (luop,150) elem(1),rom(1)
      write (luop,130)
      write (luop,10) (hlambda(j,2)/fnair(hlambda(j,2)),j=1,5)
      write (luop,30) (hydroflux(j,2),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,2)/fnair(hlambda(j,2)),j=6,10)
      write (luop,30) (hydroflux(j,2),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,2)/fnair(hlambda(j,2)),j=11,15)
      write (luop,30) (hydroflux(j,2),j=11,15)
c
  220 format(/,' Paschen n = 3 Series'/
     &         ' ======================')
      write (luop,220)
c
      write (luop,150) elem(1),rom(1)
      write (luop,130)
      write (luop,10) (hlambda(j,3)/fnair(hlambda(j,3)),j=1,5)
      write (luop,30) (hydroflux(j,3),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,3)/fnair(hlambda(j,3)),j=6,10)
      write (luop,30) (hydroflux(j,3),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,3)/fnair(hlambda(j,3)),j=11,15)
      write (luop,30) (hydroflux(j,3),j=11,15)
c
  230 format(/,' Brackett n = 4 Series'/
     &         ' ======================')
      write (luop,230)
c
      write (luop,130)
      write (luop,150) elem(1),rom(1)
      write (luop,10) (hlambda(j,4)/fnair(hlambda(j,4)),j=1,5)
      write (luop,30) (hydroflux(j,4),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,4)/fnair(hlambda(j,4)),j=6,10)
      write (luop,30) (hydroflux(j,4),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,4)/fnair(hlambda(j,4)),j=11,15)
      write (luop,30) (hydroflux(j,4),j=11,15)
c
  240 format(/,' Pfund n = 5 Series'/
     &         ' ======================')
      write (luop,240)
c
      write (luop,130)
      write (luop,150) elem(1),rom(1)
      write (luop,10) (hlambda(j,5)/fnair(hlambda(j,5)),j=1,5)
      write (luop,30) (hydroflux(j,5),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,5)/fnair(hlambda(j,5)),j=6,10)
      write (luop,30) (hydroflux(j,5),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,5)/fnair(hlambda(j,5)),j=11,15)
      write (luop,30) (hydroflux(j,5),j=11,15)
c
  250 format(/,' Humphreys n = 6 Series'/
     &         ' ======================')
      write (luop,250)
c
      write (luop,130)
      write (luop,150) elem(1),rom(1)
      write (luop,10) (hlambda(j,6)/fnair(hlambda(j,6)),j=1,5)
      write (luop,30) (hydroflux(j,6),j=1,5)
c
      write (luop,130)
      write (luop,10) (hlambda(j,6)/fnair(hlambda(j,6)),j=6,10)
      write (luop,30) (hydroflux(j,6),j=6,10)
c
      write (luop,130)
      write (luop,10) (hlambda(j,6)/fnair(hlambda(j,6)),j=11,15)
      write (luop,30) (hydroflux(j,6),j=11,15)
c
      if (zmap(2).ne.0) then
c
c        write (luop,260) heii2qa,h2qflux(2)
        write (luop,260) h2qflux(2)
  260 format(//,' Helium II Spectrum, includes collisons:'/
     &          ' ==============================================='/
     &          ' HeII Two-Photon emission :',1pg10.3)
c
  270 format(/,' n = 1 Series'/
     &         ' ============='/
     &         ' Bowen Total :', 1pg10.3)
        write (luop,270) heiioiiibfsum
c
        do i=1,10
          if (heliflux(i,1).lt.chklim) heliflux(i,1)=0.d0
        enddo
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
        write (luop,10) (helambda(j,1)/fnair(helambda(j,1)),j=1,5)
        write (luop,30) (heliflux(j,1),j=1,5)
c
        write (luop,130)
        write (luop,10) (helambda(j,1)/fnair(helambda(j,1)),j=6,10)
        write (luop,30) (heliflux(j,1),j=6,10)
c
        write (luop,130)
        write (luop,10) (helambda(j,1)/fnair(helambda(j,1)),j=11,15)
        write (luop,30) (heliflux(j,1),j=11,15)
c
  280 format(/,' n = 2 Series'/
     &         ' ============')
        write (luop,280)
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
        write (luop,10) (helambda(j,2)/fnair(helambda(j,2)),j=1,5)
        write (luop,30) (heliflux(j,2),j=1,5)
c
        write (luop,130)
        write (luop,10) (helambda(j,2)/fnair(helambda(j,2)),j=6,10)
        write (luop,30) (heliflux(j,2),j=6,10)
c
        write (luop,130)
        write (luop,10) (helambda(j,2)/fnair(helambda(j,2)),j=11,15)
        write (luop,30) (heliflux(j,2),j=11,15)
c
  290 format(/,' n = 3 Series'/
     &         ' ============')
        write (luop,290)
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
        write (luop,10) (helambda(j,3)/fnair(helambda(j,3)),j=1,5)
        write (luop,30) (heliflux(j,3),j=1,5)
c
        write (luop,130)
        write (luop,10) (helambda(j,3)/fnair(helambda(j,3)),j=6,10)
        write (luop,30) (heliflux(j,3),j=6,10)
c
        write (luop,130)
        write (luop,10) (helambda(j,3)/fnair(helambda(j,3)),j=11,15)
        write (luop,30) (heliflux(j,3),j=11,15)
c
  300 format(/,' Pickering Series'/
     &         ' ================')
        write (luop,300)
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
        write (luop,10) (helambda(j,4)/fnair(helambda(j,4)),j=1,5)
        write (luop,30) (heliflux(j,4),j=1,5)
c
        write (luop,130)
        write (luop,10) (helambda(j,4)/fnair(helambda(j,4)),j=6,10)
        write (luop,30) (heliflux(j,4),j=6,10)
c
        write (luop,130)
        write (luop,10) (helambda(j,4)/fnair(helambda(j,4)),j=11,15)
        write (luop,30) (heliflux(j,4),j=11,15)
c
  310 format(/,' n = 5 Series'/
     &         ' ============')
        write (luop,310)
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
        write (luop,10) (helambda(j,5)/fnair(helambda(j,5)),j=1,5)
        write (luop,30) (heliflux(j,5),j=1,5)
c
        write (luop,130)
        write (luop,10) (helambda(j,5)/fnair(helambda(j,5)),j=6,10)
        write (luop,30) (heliflux(j,5),j=6,10)
c
        write (luop,130)
        write (luop,10) (helambda(j,5)/fnair(helambda(j,5)),j=11,15)
        write (luop,30) (heliflux(j,5),j=11,15)
c
  320 format(/,' n = 6 Series'/
     &         ' ============')
        write (luop,320)
c
        write (luop,150) elem(2),rom(2)
        write (luop,130)
c
        write (luop,10) (helambda(j,6)/fnair(helambda(j,6)),j=1,5)
        write (luop,30) (heliflux(j,6),j=1,5)
        write (luop,130)
c
        write (luop,10) (helambda(j,6)/fnair(helambda(j,6)),j=6,10)
        write (luop,30) (heliflux(j,6),j=6,10)
        write (luop,130)
c
        write (luop,130)
        write (luop,10) (helambda(j,6)/fnair(helambda(j,6)),j=11,15)
        write (luop,30) (heliflux(j,6),j=11,15)
c
        write (luop,330)
  330 format(//,' Helium I Singlet Spectrum: HBeta = 1.000'/
     &          ' ==============================================='/)
        write (luop,90) (heislamid(j),j=1,5)
        write (luop,10) (heislam(j)/fnair(heislam(j)),j=1,5)
        write (luop,30) (fluxheis(j),j=1,5)
        write (luop,130)
        write (luop,90) (heislamid(j),j=6,10)
        write (luop,10) (heislam(j)/fnair(heislam(j)),j=6,10)
        write (luop,30) (fluxheis(j),j=6,10)
        write (luop,130)
        write (luop,90) (heislamid(j),j=11,15)
        write (luop,10) (heislam(j)/fnair(heislam(j)),j=11,15)
        write (luop,30) (fluxheis(j),j=11,15)
        write (luop,130)
c
        write (luop,340)
  340 format(//,' Helium I Triplet Spectrum: HBeta = 1.000'/
     &          ' ==============================================='/)
        write (luop,90) (heitlamid(j),j=1,5)
        write (luop,10) (heitlam(j)/fnair(heitlam(j)),j=1,5)
        write (luop,30) (fluxheit(j),j=1,5)
        write (luop,130)
        write (luop,90) (heitlamid(j),j=6,10)
        write (luop,10) (heitlam(j)/fnair(heitlam(j)),j=6,10)
        write (luop,30) (fluxheit(j),j=6,10)
        write (luop,130)
        write (luop,90) (heitlamid(j),j=11,15)
        write (luop,10) (heitlam(j)/fnair(heitlam(j)),j=11,15)
        write (luop,30) (fluxheit(j),j=11,15)
        write (luop,130)
        write (luop,90) (heitlamid(j),j=16,18)
        write (luop,10) (heitlam(j)/fnair(heitlam(j)),j=16,18)
        write (luop,30) (fluxheit(j),j=16,18)
c
        write (luop,350) hei2qa
  350 format(//,' Legacy Helium I Spectrum: HBeta = 1.000'/
     &          ' ==============================================='/
     &          ' HeI Two-Photon emission :',1pg10.3)
c
        write (luop,150) elem(2),rom(1)
        write (luop,130)
        write (luop,10) ((heilam(j)*1.d8)/fnair(heilam(j)*1.d8),j=1,5)
        write (luop,30) (fluxhei(j),j=1,5)
        write (luop,130)
        write (luop,10) ((heilam(j)*1.d8)/fnair(heilam(j)*1.d8),j=6,10)
        write (luop,30) (fluxhei(j),j=6,10)
        write (luop,130)
        write (luop,10) ((heilam(j)*1.d8)/fnair(heilam(j)*1.d8),j=11,15)
        write (luop,30) (fluxhei(j),j=11,15)
c
      endif
c
      if (atypes.gt.2) then
        write (luop,360)
  360    format(//,' Heavy Atom Two Photon Continuum: HBeta = 1.000'/
     &             ' ===============================================')
        idx=0
        do atom=3,atypes
          if (h2qflux(atom).gt.chklim) idx=idx+1
        enddo
        if (idx.gt.0) then
          do atom=3,atypes
            if (h2qflux(atom).gt.chklim) then
  370          format(' ',a2,' Two-Photon emission :',1pg10.3)
              write (luop,370) elem(atom),h2qflux(atom)
            endif
          enddo
        endif
c
        idx=0
        do atom=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              if (xhydroflux(line,series,atom).gt.chklim) idx=idx+1
            enddo
          enddo
        enddo
        if (idx.gt.0) then
          write (luop,380)
  380 format(//,' Heavy Atom Hydrogenic Spectrum, pure recomb.',/
     &          ' (Collision components in Collisonal Cascade below)',/
     &          ' =============================',
     &          '==============================')
c
  390 format(/,a2,' n = 1 Series'/
     &         ' =============')
  400 format(/,a2,' n = 2 Series'/
     &         ' ============')
c
          do atom=3,atypes
c
            if (xhydroflux(1,1,atom).gt.chklim) then
c
              write (luop,390) elem(atom)
c
              write (luop,130)
              write (luop,150) elem(atom),rom(mapz(atom))
              write (luop,130)
              write (luop,10)
     &        (xhlambda(j,1,atom)/fnair(xhlambda(j,1,atom)),j=1,5)
              write (luop,30) (xhydroflux(j,1,atom),j=1,5)
c
              write (luop,130)
              write (luop,10)
     &        (xhlambda(j,1,atom)/fnair(xhlambda(j,1,atom)),j=6,10)
              write (luop,30) (xhydroflux(j,1,atom),j=6,10)
c
              write (luop,130)
              write (luop,400) elem(atom)
c
              write (luop,150) elem(atom),rom(mapz(atom))
              write (luop,130)
              write (luop,10)
     &        (xhlambda(j,2,atom)/fnair(xhlambda(j,2,atom)),j=1,5)
              write (luop,30) (xhydroflux(j,2,atom),j=1,5)
c
              write (luop,130)
              write (luop,10)
     &        (xhlambda(j,2,atom)/fnair(xhlambda(j,2,atom)),j=6,10)
              write (luop,30) (xhydroflux(j,2,atom),j=6,10)
c
            endif
c
c     end atom loop
c
          enddo
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c CII   CIII->CII recombination spectrum
c
      idx=0
      do j=1,nrccii
        if (fluxrccii_a(j).gt.chklim) idx=idx+1
        if (fluxrccii_b(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,410)
  410 format(//,' C II Doublet Recomb. Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrccii,5
          if ((idx+nr).gt.nrccii) nr=nrccii-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrccii_a(j).gt.chklim) lj=lj+1
              if (fluxrccii_b(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,100) (rccii_tr(idx+j),j=1,nr)
              write (luop,10)
     &        (rccii_lam(idx+j)/fnair(rccii_lam(idx+j)),j=1,nr)
              write (luop,40) (fluxrccii_a(idx+j),j=1,nr)
              write (luop,50) (fluxrccii_b(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NII   NIII->NII recombination spectrum
c
      idx=0
      do j=1,nrcnii
        if (fluxrcnii_a(j).gt.chklim) idx=idx+1
        if (fluxrcnii_b(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,420)
  420 format(//,' N II Recombination Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrcnii,5
          if ((idx+nr).gt.nrcnii) nr=nrcnii-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrcnii_a(j).gt.chklim) lj=lj+1
              if (fluxrcnii_b(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,110) (rcnii_tr(idx+j),j=1,nr)
              write (luop,10)
     &        (rcnii_lam(idx+j)/fnair(rcnii_lam(idx+j)),j=1,nr)
              write (luop,40) (fluxrcnii_a(idx+j),j=1,nr)
              write (luop,50) (fluxrcnii_b(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c OI   OII->OI recombination spectrum, pure quintets and triplets
c
c
      idx=0
      do j=1,nrcoi_q
        if (fluxrcoi_qa(j).gt.chklim) idx=idx+1
        if (fluxrcoi_qa(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,430)
  430 format(//,' O I Quintet Recomb. Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrcoi_q,5
          if ((idx+nr).gt.nrcoi_q) nr=nrcoi_q-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrcoi_qa(j).gt.chklim) lj=lj+1
              if (fluxrcoi_qb(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,110) (rcoi_qtr(idx+j),j=1,nr)
              write (luop,10)
     &        (rcoi_qlam(idx+j)/fnair(rcoi_qlam(idx+j)),j=1,nr)
              write (luop,40) (fluxrcoi_qa(idx+j),j=1,nr)
              write (luop,50) (fluxrcoi_qa(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
c
      idx=0
      do j=1,nrcoi_t
        if (fluxrcoi_qa(j).gt.chklim) idx=idx+1
        if (fluxrcoi_qb(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,440)
  440 format(//,' O I Triplet Recomb. Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrcoi_t,5
          if ((idx+nr).gt.nrcoi_t) nr=nrcoi_t-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrcoi_ta(j).gt.chklim) lj=lj+1
              if (fluxrcoi_tb(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,110) (rcoi_ttr(idx+j),j=1,nr)
              write (luop,10)
     &        (rcoi_tlam(idx+j)/fnair(rcoi_tlam(idx+j)),j=1,nr)
              write (luop,40) (fluxrcoi_ta(idx+j),j=1,nr)
              write (luop,50) (fluxrcoi_qb(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c OII   OIII->OII recombination spectrum
c
      idx=0
      do j=1,nrcoii
        if (fluxrcoii_a(j).gt.chklim) idx=idx+1
        if (fluxrcoii_b(j).gt.chklim) idx=idx+1
        if (fluxrcoii_c(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,450)
  450 format(//,' O II Recomb. Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrcoii,5
          if ((idx+nr).gt.nrcoii) nr=nrcoii-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrcoii_a(j).gt.chklim) lj=lj+1
              if (fluxrcoii_b(j).gt.chklim) lj=lj+1
              if (fluxrcoii_c(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,120) (rcoii_tr(idx+j),j=1,nr)
              write (luop,10)
     &        (rcoii_lam(idx+j)/fnair(rcoii_lam(idx+j)),j=1,nr)
              write (luop,40) (fluxrcoii_a(idx+j),j=1,nr)
              write (luop,50) (fluxrcoii_b(idx+j),j=1,nr)
              write (luop,60) (fluxrcoii_c(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NeII   NeIII->NeII recombination spectrum, 6 brightest multiplets
c                    resolved into lines
c
      idx=0
      do j=1,nrcneii
        if (fluxrcneii_a(j).gt.chklim) idx=idx+1
        if (fluxrcneii_b(j).gt.chklim) idx=idx+1
      enddo
      if (idx.gt.0) then
        write (luop,460)
  460 format(//,' Ne II Recomb. Spectrum: HBeta = 1.000'/
     &          ' ============================================'/)
        nr=5
        do idx=0,nrcneii,5
          if ((idx+nr).gt.nrcneii) nr=nrcneii-idx
          if (nr.ge.1) then
            lj=0
            do j=1,nr
              if (fluxrcneii_a(j).gt.chklim) lj=lj+1
              if (fluxrcneii_b(j).gt.chklim) lj=lj+1
            enddo
            if (lj.gt.0) then
              write (luop,110) (rcneii_tr(idx+j),j=1,nr)
              write (luop,10)
     &        (rcneii_lam(idx+j)/fnair(rcneii_lam(idx+j)),j=1,nr)
              write (luop,40) (fluxrcneii_a(idx+j),j=1,nr)
              write (luop,50) (fluxrcneii_b(idx+j),j=1,nr)
              write (luop,130)
            endif
          endif
        enddo
      endif
c
      if (nfmions.ne.0) then
c
        write (luop,470)
  470 format(//,' Multi-Level Atom Models: HBeta = 1.000'/
     &          ' ============================================'/)
        chcksum=0.0d0
        do i=1,nfmions
          do j=1,nfmtrans(i)
            chcksum=chcksum+fluxm(j,i)
          enddo
        enddo
  480 format(/,' Total Multi-Level Line Int/H-Beta:',1pg14.6)
        write (luop,480) chcksum
        chcksum=0.0d0
        do i=1,nfmions
          do j=1,nfmtrans(i)
            if (fmlam(j,i).lt.(iphlam*1.d-8)) then
              chcksum=chcksum+fluxm(j,i)
            endif
          enddo
        enddo
  490 format(' Ionizing Multi-Level Lines Int :',1pg14.6)
        write (luop,490) chcksum
c
        do ionidx=1,nfmions
          idx=0
          do j=1,nfmtrans(ionidx)
            if (fluxm(j,ionidx).gt.0.d0) idx=idx+1
          enddo
          if (idx.gt.0) then
            nl=(idx/5)
            nr=idx-(nl*5)
            write (luop,130)
            write (luop,170) elem(fmatom(ionidx)),rom(fmion(ionidx)),
     &       fmnl(ionidx)
            jdx=1
            do nd=1,nl
              ns=((nd-1)*5)+1
              ne=ns+4
              chcksum=0.0d0
c
c get next 5 lines with flux
c
              idx=1
  500         if (fluxm(jdx,ionidx).gt.0.d0) then
                lineid(idx)=nfmid(jdx,ionidx)
                lineflx(idx)=fluxm(jdx,ionidx)
                linelam(idx)=(fmlam(jdx,ionidx)*1.d8)
                linelam(idx)=linelam(idx)/fnair(linelam(idx))
                linelo(idx)=nfmlower(jdx,ionidx)
                linehi(idx)=nfmupper(jdx,ionidx)
                idx=idx+1
              endif
              jdx=jdx+1
              if (idx.gt.5) goto 510
              goto 500
  510         continue
c
              write (luop,130)
              write (luop,70) (linelo(j),linehi(j),j=1,5)
              write (luop,110) (lineid(j),j=1,5)
              write (luop,10) (linelam(j),j=1,5)
              write (luop,30) (lineflx(j),j=1,5)
c
            enddo
            if (nr.gt.0) then
              idx=1
  520         if (fluxm(jdx,ionidx).gt.0.d0) then
                lineid(idx)=nfmid(jdx,ionidx)
                lineflx(idx)=fluxm(jdx,ionidx)
                linelam(idx)=(fmlam(jdx,ionidx)*1.d8)
                linelam(idx)=linelam(idx)/fnair(linelam(idx))
                linelo(idx)=nfmlower(jdx,ionidx)
                linehi(idx)=nfmupper(jdx,ionidx)
                idx=idx+1
              endif
              jdx=jdx+1
              if (idx.gt.nr) goto 530
              goto 520
  530         continue
c
              write (luop,130)
              write (luop,70) (linelo(j),linehi(j),j=1,nr)
              write (luop,110) (lineid(j),j=1,nr)
              write (luop,10) (linelam(j),j=1,nr)
              write (luop,30) (lineflx(j),j=1,nr)
c
            endif
          endif
        enddo
c
      endif
c
      if (nfeions.ne.0) then
c
        write (luop,540)
  540 format(//,' Multi-Level Iron-Like Ions: HBeta = 1.000'/
     &          ' ============================================'/)
        chcksum=0.0d0
        do i=1,nfeions
          do j=1,nfetrans(i)
            chcksum=chcksum+fluxfe(j,i)
          enddo
        enddo
  550 format(' Total Iron Lines Int/H-Beta:',1pg14.6)
        write (luop,550) chcksum
        chcksum=0.0d0
        do i=1,nfeions
          do j=1,nfetrans(i)
            if (felam(j,i).lt.3.d-5) then
              chcksum=chcksum+fluxfe(j,i)
            endif
          enddo
        enddo
  560 format(' UV Total Iron Lines <3000A Int/H-Beta:',1pg14.6)
        write (luop,560) chcksum
        chcksum=0.0d0
        do i=1,nfeions
          do j=1,nfetrans(i)
            if (felam(j,i).ge.7.d-5) then
              chcksum=chcksum+fluxfe(j,i)
            endif
          enddo
        enddo
  570 format(' IR Total Iron Lines >7000A Int/H-Beta:',1pg14.6)
        write (luop,570) chcksum
        do ionidx=1,nfeions
          idx=0
          do j=1,nfetrans(ionidx)
            if (fluxfe(j,ionidx).gt.chklim) idx=idx+1
          enddo
          if (idx.gt.0) then
            nl=(idx/5)
            nr=idx-(nl*5)
            write (luop,140)
            write (luop,170) elem(featom(ionidx)),rom(feion(ionidx)),
     &       fenl(ionidx)
c
            chcksum=0.0d0
            do j=1,nfetrans(ionidx)
              chcksum=chcksum+fluxfe(j,ionidx)
            enddo
  580 format(/,' Total ',a3,a6,' Lines Int/H-Beta:',1pg14.6)
            write (luop,580) elem(featom(ionidx)),rom(feion(ionidx)),
     &       chcksum
c
            chcksum=0.0d0
            do j=1,nfetrans(ionidx)
              if (felam(j,ionidx).lt.3.d-5) then
                chcksum=chcksum+fluxfe(j,ionidx)
              endif
            enddo
  590 format(' UV ',a3,a6,' Lines <3000A Int/H-Beta:',1pg14.6)
            write (luop,590) elem(featom(ionidx)),rom(feion(ionidx)),
     &       chcksum
c
            chcksum=0.0d0
            do j=1,nfetrans(ionidx)
              if (felam(j,ionidx).ge.7.d-5) then
                chcksum=chcksum+fluxfe(j,ionidx)
              endif
            enddo
  600 format(' IR ',a3,a6,' Lines >7000A Int/H-Beta:',1pg14.6)
            write (luop,600) elem(featom(ionidx)),rom(feion(ionidx)),
     &       chcksum
c
            idx=0
            nr=0
            nt=nfetrans(ionidx)
c
            do jdx=1,nt
              nlo=nfelower(jdx,ionidx)
              nhi=nfeupper(jdx,ionidx)
              trdx=nfetridx(nlo,nhi,ionidx)
              if (fluxfe(jdx,ionidx).gt.chklim) then
c
c get next 5 lines with flux
c
                idx=idx+1
                nr=idx
                lineid(idx)=nfeid(jdx,ionidx)
                lineflx(idx)=fluxfe(jdx,ionidx)
                linelam(idx)=(felam(jdx,ionidx)*1.d8)
                linelam(idx)=linelam(idx)/fnair(linelam(idx))
                linelo(idx)=nfelower(jdx,ionidx)
                linehi(idx)=nfeupper(jdx,ionidx)
                if ((idx.ge.5).or.((jdx.eq.nt).and.(nr.gt.0))) then
                  write (luop,140)
                  write (luop,80) (linelo(j),linehi(j),j=1,nr)
                  write (luop,110) (lineid(j),j=1,nr)
                  write (luop,10) (linelam(j),j=1,nr)
                  write (luop,30) (lineflx(j),j=1,nr)
                  idx=0
                  nr=0
                endif
              endif
              if ((jdx.eq.nt).and.(nr.gt.0)) then
                write (luop,140)
                write (luop,80) (linelo(j),linehi(j),j=1,nr)
                write (luop,110) (lineid(j),j=1,nr)
                write (luop,10) (linelam(j),j=1,nr)
                write (luop,30) (lineflx(j),j=1,nr)
                idx=0
                nr=0
              endif
            enddo
c
          endif
        enddo
c
      endif
c
      if (nf3ions.ne.0) then
c
        write (luop,610)
  610 format(//,' Three-Level H/Proton Excitation: HBeta = 1.000'/
     &          ' ================================================'/)
c
        chcksum=0.0d0
        do i=1,nf3ions
          do j=1,nf3trans
            chcksum=chcksum+fluxf3(j,i)
          enddo
        enddo
  620 format(/,' Total 3LA Lines Int/H-Beta:',1pg14.6)
        write (luop,620) chcksum
c
        do i=1,nf3ions
          chcksum=0.0d0
          do j=1,nf3trans
            chcksum=chcksum+fluxf3(j,i)
          enddo
          if (chcksum.ge.chklim) then
            write (luop,130)
            write (luop,150) elem(f3atom(i)),rom(f3ion(i))
            write (luop,20) (f3lam(j,i),j=1,nf3trans)!nocorrection
            write (luop,30) (fluxf3(j,i),j=1,nf3trans)
          endif
        enddo
c
      endif
c
c
c      if ((mlines+xilines+xhelines).ne.0) then
      if ((mlines).ne.0) then
c
        max=(mlines+4)/5
        write (luop,630)
  630 format(//,' Two Level FS and Semi-Forbidden: HBeta = 1.000'/
     &          ' ================================================'/)
        chcksum=0.0d0
        do i=1,mlines
          chcksum=chcksum+fluxi(i)
        enddo
c       do i=1,xilines
c         chcksum=chcksum+fluxxi(i)
c       enddo
  640 format(/,' Total Semi-Forbidden Line Int/H-Beta:',1pg14.6)
        write (luop,640) chcksum
        chcksum=0.0d0
        do i=1,mlines
          if (fslam(i).lt.(iphlam*1.d-8)) then
            chcksum=chcksum+fluxi(i)
          endif
        enddo
c       do i=1,xilines
c         if (xilam(i).lt.iphlam) then
c           chcksum=chcksum+fluxxi(i)
c         endif
c       enddo
  650 format(/,' Semi-Forbidden Lines < 911A Int/H-Beta:',1pg14.6)
        write (luop,650) chcksum
        do i=1,max
          jl=5
          k=(i-1)*5
          if ((i*5).gt.mlines) jl=mlines-k
          chcksum=0.0d0
          do j=1,jl
            chcksum=chcksum+fluxi(j+k)
            atom=ielfs(k+j)
            ion=ionfs(k+j)
            pad=elem(atom)
            ionstr=pad(1:3)//rom(ion)
            l0=rom_len(ion)
            l1=9-(l0+3)
            pad='         '
            if (l1.gt.0) then
              ionstr=pad(1:l1)//ionstr
            endif
            lineion(j)=ionstr(1:9)
            linelam(j)=fslam(k+j)*1.d8
            linelam(j)=linelam(j)/fnair(linelam(j))
          enddo
          if (chcksum.ge.chklim) then
            write (luop,130)
            write (luop,160) (lineion(j),j=1,jl)
            write (luop,10) (linelam(j),j=1,jl)
            write (luop,30) (fluxi(j+k),j=1,jl)
          endif
        enddo
c
c
c       max=(xilines+4)/5
c       write (luop,670)
c 670 format(//)
c       do 680 i=1,max
c         jl=5
c         k=(i-1)*5
c         if ((i*5).gt.xilines) jl=xilines-k
c         chcksum=0.0d0
c         do j=1,jl
c           chcksum=chcksum+fluxxi(j+k)
c           atom=xiat(k+j)
c           ion=xiion(k+j)
c           pad=elem(atom)
c           ionstr=pad(1:3)//rom(ion)
c           l0=rom_len(ion)
c           l1=9-(l0+3)
c           pad='         '
c           if (l1.gt.0) then
c             ionstr=pad(1:l1)//ionstr
c           endif
c           lineion(j)=ionstr(1:9)
c           linelam(j)=xilam(k+j)
c           linelam(j)=linelam(j)/fnair(linelam(j))
c         enddo
c         if (chcksum.ge.chklim) then
c           write (luop,130)
c           write (luop,160) (lineion(j),j=1,jl)
c           write (luop,10) (linelam(j),j=1,jl)
c           write (luop,30) (fluxxi(j+k),j=1,jl)
c         endif
c 680   continue
c
c       max=(xhelines+4)/5
c       write (luop,675)
c 675 format(//)
c       do 685 i=1,max
c         jl=5
c         k=(i-1)*5
c         if ((i*5).gt.xhelines) jl=xhelines-k
c         chcksum=0.0d0
c         do j=1,jl
c           chcksum=chcksum+fheif(j+k)
c           atom=xheat(k+j)
c           ion=xheion(k+j)
c           pad=elem(atom)
c           ionstr=pad(1:3)//rom(ion)
c           l0=rom_len(ion)
c           l1=9-(l0+3)
c           pad='         '
c           if (l1.gt.0) then
c             ionstr=pad(1:l1)//ionstr
c           endif
c           lineion(j)=ionstr(1:9)
c           linelam(j)=xhelam(k+j)
c           linelam(j)=linelam(j)/fnair(linelam(j))
c         enddo
c         if (chcksum.ge.chklim) then
c           write (luop,130)
c           write (luop,160) (lineion(j),j=1,jl)
c           write (luop,10) (linelam(j),j=1,jl)
c           write (luop,30) (fheif(j+k),j=1,jl)
c         endif
c 685   continue
c
      endif
c
      if ((nxr3lines).ne.0) then
        max=(nxr3lines+4)/5
        write (luop,660)
  660 format(// ' Cascade Resonance Lines: HBeta = 1.000'/
     &          ' ================================================'/)
        chcksum=0.0d0
        do i=1,nxr3lines
          chcksum=chcksum+xr3lines_flux(i)
        enddo
  670 format(' Total Cascade Resonance Line   Int/H-Beta:',1pg14.6)
        write (luop,670) chcksum
        chcksum=0.0d0
        do i=1,nxr3lines
          if (xr3lines_lam(i).lt.iphlam) then
            chcksum=chcksum+xr3lines_flux(i)
          endif
        enddo
  680 format(' Cascade Resonance Lines < 911A Int/H-Beta:',1pg14.6)
        write (luop,680) chcksum
        do i=1,max
          jl=5
          k=(i-1)*5
          if ((i*5).gt.nxr3lines) jl=nxr3lines-k
          chcksum=0.0d0
          do j=1,jl
            chcksum=chcksum+xr3lines_flux(j+k)
            atom=xr3lines_at(k+j)
            ion=xr3lines_ion(k+j)
            pad=elem(atom)
            ionstr=pad(1:3)//rom(ion)
            l0=rom_len(ion)
            l1=9-(l0+3)
            pad='         '
            if (l1.gt.0) then
              ionstr=pad(1:l1)//ionstr
            endif
            lineion(j)=ionstr(1:9)
            ionstr=xr3lines_tran(j+k)
            lineid(j)=ionstr(1:18)
            linelo(j)=xr3lines_imap_i(j+k)
            linehi(j)=xr3lines_imap_j(j+k)
            linelam(j)=xr3lines_lam(j+k)
            linelam(j)=linelam(j)/fnair(linelam(j))
            lineflx(j)=xr3lines_flux(j+k)
          enddo
          if (chcksum.ge.chklim) then
            write (luop,130)
            write (luop,160) (lineion(j),j=1,jl)
            write (luop,80)  (linelo(j),linehi(j),j=1,jl)
            write (luop,120) (lineid(j),j=1,jl)
            write (luop,10) (linelam(j),j=1,jl)
            write (luop,30) (lineflx(j),j=1,jl)
          endif
        enddo
      endif
c
      if ((nxrllines).ne.0) then
        max=(nxrllines+4)/5
        write (luop,690)
  690 format(// ' Large Atom Cascade Resonance Lines: HBeta = 1.000'/
     &          ' =================================================='/)
        chcksum=0.0d0
        do i=1,nxrllines
          chcksum=chcksum+xrllines_flux(i)
        enddo
  700 format(' Total Cascade Resonance Line   Int/H-Beta:',1pg14.6)
        write (luop,700) chcksum
        chcksum=0.0d0
        do i=1,nxrllines
          if (xrllines_lam(i).lt.iphlam) then
            chcksum=chcksum+xrllines_flux(i)
          endif
        enddo
  710 format(' Cascade Resonance Lines < 911A Int/H-Beta:',1pg14.6)
        write (luop,710) chcksum
        do i=1,max
          jl=5
          k=(i-1)*5
          if ((i*5).gt.nxrllines) jl=nxrllines-k
          chcksum=0.0d0
          do j=1,jl
            chcksum=chcksum+xrllines_flux(j+k)
            atom=xrllines_at(k+j)
            ion=xrllines_ion(k+j)
            pad=elem(atom)
            ionstr=pad(1:3)//rom(ion)
            l0=rom_len(ion)
            l1=9-(l0+3)
            pad='         '
            if (l1.gt.0) then
              ionstr=pad(1:l1)//ionstr
            endif
            lineion(j)=ionstr(1:9)
            ionstr=xrllines_tran(j+k)
            lineid(j)=ionstr(1:18)
            linelo(j)=xrllines_imap_i(j+k)
            linehi(j)=xrllines_imap_j(j+k)
            linelam(j)=xrllines_lam(j+k)
            linelam(j)=linelam(j)/fnair(linelam(j))
            lineflx(j)=xrllines_flux(j+k)
          enddo
          if (chcksum.ge.chklim) then
            write (luop,130)
            write (luop,160) (lineion(j),j=1,jl)
            write (luop,80) (linelo(j),linehi(j),j=1,jl)
            write (luop,120) (lineid(j),j=1,jl)
            write (luop,10) (linelam(j),j=1,jl)
            write (luop,30) (lineflx(j),j=1,jl)
          endif
        enddo
      endif
c
      if ((nlines).ne.0) then
c
        max=(nlines+4)/5
        write (luop,720)
  720 format(//' Two Level Resonance Lines: HBeta = 1.000'/
     &          ' ================================================'/)
        chcksum=0.0d0
        do i=1,nlines
          chcksum=chcksum+fluxr(i)
        enddo
c       do i=1,xlines
c         chcksum=chcksum+fluxx(i)
c       enddo
  730 format(/,' Total Resonance Line Int/H-Beta:',1pg14.6)
        write (luop,730) chcksum
        chcksum=0.0d0
        do i=1,nlines
          if (rlam(i).lt.(iphlam*1.d-8)) then
            chcksum=chcksum+fluxr(i)
          endif
        enddo
c       do i=1,xlines
c         if (xrlam(i).lt.iphlam) then
c           chcksum=chcksum+fluxx(i)
c         endif
c       enddo
  740 format(/,' Resonance Lines < 911A Int/H-Beta:',1pg14.6)
        write (luop,740) chcksum
        do i=1,max
          jl=5
          k=(i-1)*5
          if ((i*5).gt.nlines) jl=nlines-k
          chcksum=0.0d0
          do j=1,jl
            chcksum=chcksum+fluxr(j+k)
            atom=ielr(k+j)
            ion=ionr(k+j)
            pad=elem(atom)
            ionstr=pad(1:3)//rom(ion)
            l0=rom_len(ion)
            l1=9-(l0+3)
            pad='         '
            if (l1.gt.0) then
              ionstr=pad(1:l1)//ionstr
            endif
            lineion(j)=ionstr(1:9)
            linelam(j)=rlam(k+j)*1.d8
            linelam(j)=linelam(j)/fnair(linelam(j))
          enddo
          if (chcksum.ge.chklim) then
            write (luop,130)
            write (luop,160) (lineion(j),j=1,jl)
            write (luop,10) (linelam(j),j=1,jl)
            write (luop,30) (fluxr(j+k),j=1,jl)
          endif
        enddo
c
c
c       max=(xlines+4)/5
c       write (luop,750)
c 750 format(//)
c       do 760 i=1,max
c         jl=5
c         k=(i-1)*5
c         if ((i*5).gt.xlines) jl=xlines-k
c         chcksum=0.0d0
c         do j=1,jl
c           chcksum=chcksum+fluxx(j+k)
c           atom=xrat(k+j)
c           ion=xion(k+j)
c           pad=elem(atom)
c           ionstr=pad(1:3)//rom(ion)
c           l0=rom_len(ion)
c           l1=9-(l0+3)
c           pad='         '
c           if (l1.gt.0) then
c             ionstr=pad(1:l1)//ionstr
c           endif
c           lineion(j)=ionstr(1:9)
c           linelam(j)=xrlam(k+j)
c           linelam(j)=linelam(j)/fnair(linelam(j))
c         enddo
c         if (chcksum.ge.chklim) then
c           write (luop,130)
c           write (luop,160) (lineion(j),j=1,jl)
c           write (luop,10)  (linelam(j),j=1,jl)
c           write (luop,30)  (fluxx(j+k),j=1,jl)
c         endif
c 760   continue
c
      endif
c
      write (luop,750) h2qav
  750 format(// ' Legacy Hydrogen Lines: HBeta = 1.000'/
     &          ' ================================================'/
     &          ' HI Two-Photon emission :',1pg10.3)
      write (luop,150) elem(1),rom(1)
      write (luop,130)
      write (luop,10) ((hlam(j)*1.d8),j=1,5)
      write (luop,30) (fluxh(j),j=1,5)
c
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1979-2012 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     this subroutine is meant to be
c     used to write out various spectral components in a more
c     machine useable form (for plotting packages etc)
c     In general it will make looooonnnnggg lists
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine spec2 (luop, list, mode)
c
      include 'cblocks.inc'
c
c     Variables
c
      real*8 eloslog,fhbtlog,forblog
      real*8 chklim,totallines,total2p
c
c     Storage for sorted output.
c
      real*8 linelam(mxspeclines)
      real*8 linespec(mxspeclines)
      integer*4 lineid(mxspeclines,2)
      integer*4 lineidx(mxspeclines)
      integer*4 lineacc(mxspeclines)
      character*4 linekind(mxspeclines)
      character*4 kind1, kind2
      integer*4 linecount,idx,count
c
      integer*4 i,j,line,series
      integer*4 ionidx, nt, omtype
      integer*4 luop
c
      character mode*4,list*4, llist*4
c
c function
c
      real*8 fnair
c
c     Formats:
c
c
   10 format(' ',f12.3,', ',1pe11.4,', ',1pe11.4,
     &       ', ',a3,a6,', ',a4,', ',i4)
   20 format(' ',f12.3,', ',1pe11.4,', ',1pe11.4,
     &       ', ',a3,a6,', ',2a4,', ',i4,', *')
   30 format(/' Flux Scale (H Beta = 1.000) :'/
     &        ' ============================='/
     &        ' H-beta  :,',1pe11.4, ', (log(ergs/cm^2/s)), ',
     &        1pe11.4,', (ergs/cm^2/s)'/)
   40 format(/' Flux Scale (H Beta = 1.000) :'/
     &        ' ============================='/
     &        ' H-beta  :,',1pe11.4, ', (log(ergs/s)), ',
     &        1pe11.4,', (ergs/s)'/)
   50 format(/' Total Emission Line Flux Hbeta = 1.000  :'/
     &        ' ========================================='/
     &' F/Hbeta:,',1pe11.4, ' ,:, ', 1pe11.4,', (ergs/cm^2/s)')
   60 format(/' Total Two-Photon Continuum Hbeta = 1.000 :'/
     &        ' =========================================='/
     &' F/Hbeta:,'1pe11.4, ' ,:, ', 1pe11.4,', (ergs/cm^2/s)')
   70 format(/' Flux Scale (Absolute units ergs/s) :'/
     &        ' ============================='/
     &        ' Units (Log) :,',1pe11.4, ', (log(ergs/s)), ',
     &         1pe11.4,', (ergs/s)'/)
   80 format(//' Line List >,',1pe11.4,', HBeta = 1.000 :,'/
     &         ' Number of lines:, ',i5//
     &         ' ======================================'//
     &         '     Lambda(A),  E (eV) , Flux (HB=1.0),',
     &         ' Species  , Kind, Accuracy (1-5)'/
     &         ' ====================================',
     &         '=============================')
   90 format(//' Line List >,',1pe11.4,', HBeta = 1.000 :,'/
     &         ' Number of lines:, ',i5,
     &         ' , Duplicates Combined:, ',i5//
     &         ' ======================================'//
     &         '     Lambda(A),  E (eV) , Flux (HB=1.0),',
     &         ' Species  , Kind, Accuracy (1-5)'/
     &         ' ====================================',
     &         '=============================')
  100 format(//' Flux Sorted Line List >',1pe11.4,' HBeta = 1.000 :'/
     &         ' =================================================='//
     &         '     Lambda(A),  E (eV) , Flux (HB=1.0),',
     &         ' Species  , Kind, Accuracy (1-5)'/
     &         ' ===================================',
     &         '====================================')
c
c     chklim is a lower limit for individual line fluxes to print out.
c     1e-5 saves a lot of paper
c
c     chklim = 0.d0
c
      chklim=1.0d-5
      if ((mode.eq.'ABS').or.(mode.eq.'TOTL')) then
        chklim=epsilon
      endif
      if (list.eq.'TTWN') then
        chklim=1.0d-3
      endif
c
c     Copy list into llist (local list) so it can be modified if
c     necessary.
c
      llist=list
      if ((llist.ne.'TTWN').and.(llist.ne.'TOTL').and.(llist.ne.'FLUX')
     &.and.(llist.ne.'LAMB').and.(llist.ne.'ALL')) then
        write (*,*) ' Mode error in SPEC2 : ',list
        write (*,*) ' Using lambda sorting.'
        llist='LAMB'
      endif
c
      linecount=0
c
      do series=1,nhseries
        do line=1,nhlines
          if (hydroflux(line,series).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=hlambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=hydroflux(line,series)
            lineid(linecount,1)=1
            lineid(linecount,2)=1
            lineacc(linecount)=1
            linekind(linecount)='RCAB'
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
      do series=1,nheseries
        do line=1,nhelines
          if (heliflux(line,series).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=helambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=heliflux(line,series)
            lineid(linecount,1)=2
            lineid(linecount,2)=2
            lineacc(linecount)=1
            linekind(linecount)='RCAB'
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do i=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            if (xhydroflux(line,series,i).gt.chklim) then
              linecount=linecount+1
              linelam(linecount)=xhlambda(line,series,i)
              linelam(linecount)=linelam(linecount)/
     &         fnair(linelam(linecount))
              linespec(linecount)=xhydroflux(line,series,i)
              lineid(linecount,1)=i
              lineid(linecount,2)=mapz(i)
              lineacc(linecount)=3
c no collisons in these components
              linekind(linecount)='RA  '
            endif
          enddo
        enddo
      enddo
c
      do ionidx=1,nfmions
        do nt=1,nfmtrans(ionidx)
          if (fluxm(nt,ionidx).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(fmlam(nt,ionidx)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fluxm(nt,ionidx)
            lineid(linecount,1)=fmatom(ionidx)
            lineid(linecount,2)=fmion(ionidx)
            lineacc(linecount)=3
            linekind(linecount)='CM  '
            omtype=fmomdatatype(nt,ionidx)
            if (omtype.eq.0) then
              j=nfmupper(nt,ionidx)
              i=nfmlower(nt,ionidx)
              if (tdepm(j,i,ionidx).eq.0.d0) lineacc(linecount)=4
            endif
            if ((omtype.gt.0).and.(fmombtn(nt,ionidx).ge.9)) then
              lineacc(linecount)=2
            endif
          endif
        enddo
      enddo
c
      do ionidx=1,nfeions
        do nt=1,nfetrans(ionidx)
          if (fluxfe(nt,ionidx).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(felam(nt,ionidx)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fluxfe(nt,ionidx)
            lineid(linecount,1)=featom(ionidx)
            lineid(linecount,2)=feion(ionidx)
            lineacc(linecount)=3
            linekind(linecount)='CM  '
            omtype=feomdatatype(nt,ionidx)
            if (omtype.eq.0) then
              j=nfeupper(nt,ionidx)
              i=nfelower(nt,ionidx)
              if (tdepfe(j,i,ionidx).eq.0.d0) lineacc(linecount)=4
            endif
            if ((omtype.gt.0).and.(feombtn(nt,ionidx).gt.9)) then
              lineacc(linecount)=2
            endif
          endif
        enddo
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          if (fluxf3(j,i).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(f3lam(j,i)*1.d4)!nocorrection
            linespec(linecount)=fluxf3(j,i)
            lineid(linecount,1)=f3atom(i)
            lineid(linecount,2)=f3ion(i)
            lineacc(linecount)=4
            linekind(linecount)='CHP '
          endif
        enddo
      enddo
c
      do i=1,mlines
        if (fluxi(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(fslam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxi(i)
          lineid(linecount,1)=ielfs(i)
          lineid(linecount,2)=ionfs(i)
          lineacc(linecount)=5
          linekind(linecount)='C   '
        endif
      enddo
c
c     do i=1,xilines
c       if (fluxxi(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xilam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=fluxxi(i)
c         lineid(linecount,1)=xiat(i)
c         lineid(linecount,2)=xiion(i)
c         lineacc(linecount)=5
c         linekind(linecount)='C   '
c       endif
c     enddo
c
c     do i=1,xhelines
c       if (fheif(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xhelam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=fheif(i)
c         lineid(linecount,1)=xheat(i)
c         lineid(linecount,2)=xheion(i)
c         lineacc(linecount)=5
c         linekind(linecount)='CHe '
c       endif
c     enddo
c
      do i=1,nlines
        if (fluxr(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rlam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxr(i)
          lineid(linecount,1)=ielr(i)
          lineid(linecount,2)=ionr(i)
          lineacc(linecount)=4
          linekind(linecount)='C   '
        endif
      enddo
c
      do i=1,nxr3lines
        if (xr3lines_flux(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=xr3lines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xr3lines_flux(i)
          lineid(linecount,1)=xr3lines_at(i)
          lineid(linecount,2)=xr3lines_ion(i)
          lineacc(linecount)=3
          linekind(linecount)='CC  '
        endif
      enddo
c
      do i=1,nxrllines
        if (xrllines_flux(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=xrllines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xrllines_flux(i)
          lineid(linecount,1)=xrllines_at(i)
          lineid(linecount,2)=xrllines_ion(i)
          lineacc(linecount)=3
          linekind(linecount)='CCL '
        endif
      enddo
c
c     do i=1,xlines
c       if (fluxx(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xrlam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=fluxx(i)
c         lineid(linecount,1)=xrat(i)
c         lineid(linecount,2)=xion(i)
c         lineacc(linecount)=4
c         linekind(linecount)='C   '
c       endif
c     enddo
c
c Original He I lines, just the first two
c
      do i=1,2
        if (fluxhei(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heilam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxhei(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=5
          linekind(linecount)='RCB '
        endif
      enddo
c
c New He I singlet lines
c
      do i=1,nheislines
        if (fluxheis(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heislam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxheis(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
          linekind(linecount)='RCBS'
        endif
      enddo
c
c New He I triplet lines
c
      do i=1,nheitlines
        if (fluxheit(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heitlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxheit(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
          linekind(linecount)='RCBT'
        endif
      enddo
c
      do i=1,nrccii
        if (fluxrccii_b(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rccii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrccii_b(i)
          lineid(linecount,1)=zmap(6)
          lineid(linecount,2)=2
          lineacc(linecount)=4
          linekind(linecount)='RB  '
        endif
      enddo
c
      do i=1,nrcnii
        if (fluxrcnii_b(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcnii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrcnii_b(i)
          lineid(linecount,1)=zmap(7)
          lineid(linecount,2)=2
          lineacc(linecount)=4
          linekind(linecount)='RCB '
        endif
      enddo
c
      do i=1,nrcoi_q
        if (fluxrcoi_qb(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_qlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrcoi_qb(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=4
          linekind(linecount)='RBQ '
        endif
      enddo
c
      do i=1,nrcoi_t
        if (fluxrcoi_tb(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_tlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrcoi_tb(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=5
          linekind(linecount)='RBT '
        endif
      enddo
c
      do i=1,nrcoii
        if (fluxrcoii_b(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrcoii_b(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=2
          lineacc(linecount)=4
          linekind(linecount)='RCB '
        endif
      enddo
c
      do i=1,nrcneii
        if (fluxrcneii_b(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcneii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fluxrcneii_b(i)
          lineid(linecount,1)=zmap(10)
          lineid(linecount,2)=2
          lineacc(linecount)=4
          linekind(linecount)='RB  '
        endif
      enddo
c
      if ((llist.eq.'TOTL')) then
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10((fhbeta*(fpi))+epsilon)+vunilog
          forblog=dlog10(tforbi+epsilon)+vunilog
          eloslog=dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,40) fhbtlog,10.d0**fhbtlog
          else
            write (luop,30) fhbtlog,10.d0**fhbtlog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(fhbeta*fpi+epsilon)+vunilog
          write (luop,70) fhbtlog,10.d0**fhbtlog
c
        endif
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Spec2 output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
        call heapindexsort (linecount, linelam, lineidx)
c
        totallines=0.d0
        do i=1,linecount
          totallines=totallines+linespec(lineidx(i))
        enddo
        write (luop,50) totallines,totallines*(fhbeta*fpi)
c
        total2p=hei2qa
        do i=1,atypes
          total2p=total2p+h2qflux(i)
        enddo
        write (luop,60) total2p,total2p*(fhbeta*fpi)
c
        write (luop,80) chklim,linecount
c
      endif
c
      if ((llist.eq.'LAMB').or.(llist.eq.'ALL')) then
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10((fhbeta*(fpi))+epsilon)+vunilog
          forblog=dlog10(tforbi+epsilon)+vunilog
          eloslog=dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,40) fhbtlog,10.d0**fhbtlog
          else
            write (luop,30) fhbtlog,10.d0**fhbtlog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(fhbeta*fpi+epsilon)+vunilog
          write (luop,70) fhbtlog,10.d0**fhbtlog
c
        endif
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Spec2 output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
        call heapindexsort (linecount, linelam, lineidx)
c
        totallines=0.d0
        do i=1,linecount
          totallines=totallines+linespec(lineidx(i))
        enddo
        write (luop,50) totallines,totallines*(fhbeta*fpi)
c
        total2p=hei2qa
        do i=1,atypes
          total2p=total2p+h2qflux(i)
        enddo
        write (luop,60) total2p,total2p*(fhbeta*fpi)
c        write (luop,60) heiioiiibfsum
c  Combine dublicates
c
        idx=0
        count=0
  110   idx=idx+1
        if (idx.lt.linecount) then
          if ((lineid(lineidx(idx),1).eq.lineid(lineidx(idx+1),1))
     &     .and.(lineid(lineidx(idx),2).eq.lineid(lineidx(idx+1),2))
     &     .and.(linelam(lineidx(idx)).eq.linelam(lineidx(idx+1))))
     &     then
            count=count+1
            idx=idx+1
            goto 110
          endif
        endif
        count=count+1
        if (idx.lt.linecount) goto 110
c
        write (luop,90) chklim,linecount,count
c
        idx=0
c
  120   idx=idx+1
        if (idx.lt.linecount) then
          if ((lineid(lineidx(idx),1).eq.lineid(lineidx(idx+1),1))
     &     .and.(lineid(lineidx(idx),2).eq.lineid(lineidx(idx+1),2))
     &     .and.(linelam(lineidx(idx)).eq.linelam(lineidx(idx+1))))
     &     then
            linespec(lineidx(idx))=linespec(lineidx(idx))+
     &       linespec(lineidx(idx+1))
            kind1=linekind(lineidx(idx))
            kind2=linekind(lineidx(idx+1))
            write (luop,20) linelam(lineidx(idx)),lmev/
     &       (linelam(lineidx(idx))),linespec(lineidx(idx)),
     &       elem(lineid(lineidx(idx),1)),rom(lineid(lineidx(idx),2)),
     &       kind1,kind2,lineacc(lineidx(idx))
            idx=idx+1
            goto 120
          endif
        endif
c
        write (luop,10) linelam(lineidx(idx)),lmev/(linelam(lineidx(idx)
     &   )),linespec(lineidx(idx)),elem(lineid(lineidx(idx),1)),
     &   rom(lineid(lineidx(idx),2)),linekind(lineidx(idx)),
     &   lineacc(lineidx(idx))
c
        if (idx.lt.linecount) goto 120
c
c       do i=1,linecount
cC         write (luop,10) linelam(lineidx(i)),lmev/(linelam(lineidx(i)))
cC    &     ,linespec(lineidx(i)),elem(lineid(lineidx(i),1)),
cC    &     rom(lineid(lineidx(i),2)),linekind(lineidx(i)),
cC    &     lineacc(lineidx(i))
c       enddo
c
      endif
c
      if ((llist.eq.'FLUX').or.(llist.eq.'ALL')) then
c
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10((fhbeta*(fpi))+epsilon)+vunilog
          forblog=dlog10(tforbi+epsilon)+vunilog
          eloslog=dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,40) fhbtlog,10.d0**fhbtlog
          else
            write (luop,30) fhbtlog,10.d0**fhbtlog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(fhbeta*fpi+epsilon)+vunilog
          write (luop,70) fhbtlog
c
        endif
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Spec2 output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
        call heapindexsort (linecount, linespec, lineidx)
c
        totallines=0.d0
        do i=1,linecount
          totallines=totallines+linespec(lineidx(i))
        enddo
        write (luop,50) totallines,totallines*(fhbeta*fpi)
c
        total2p=hei2qa
        do i=1,atypes
          total2p=total2p+h2qflux(i)
        enddo
        write (luop,60) total2p,total2p*(fhbeta*fpi)
c        write (luop,60) heiioiiibfsum
c
        write (luop,100) chklim
c
        do i=linecount,1,-1
         write (luop,10) linelam(lineidx(i)),lmev/(linelam(lineidx(i)))
     &   ,linespec(lineidx(i)),elem(lineid(lineidx(i),1))
     &   ,rom(lineid(lineidx(i),2)),linekind(lineidx(i))
     &   ,lineacc(lineidx(i))
        enddo
c
      endif
c
c
      if (llist.eq.'TTWN') then
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10(fhbeta*fpi+epsilon)+vunilog
          forblog=dlog10(tforbi+epsilon)+vunilog
          eloslog=dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,40) fhbtlog,10.d0**fhbtlog
          else
            write (luop,30) fhbtlog,10.d0**fhbtlog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(fhbeta*fpi+epsilon)+vunilog
          write (luop,70) fhbtlog
c
        endif
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Spec2 output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
c sort on line flux:
c
c        call heapindexsort (linecount, linespec, lineidx)
        call heapindexsort (linecount, linelam, lineidx)
c
        totallines=0.d0
        do i=1,linecount
          totallines=totallines+linespec(lineidx(i))
        enddo
        write (luop,50) totallines,totallines*(fhbeta*fpi)
c
        total2p=hei2qa
        do i=1,atypes
          total2p=total2p+h2qflux(i)
        enddo
        write (luop,60) total2p,total2p*(fhbeta*fpi)
c
        write (luop,80) chklim,linecount
c
        do i=1,linecount
          write (luop,10) linelam(lineidx(i)),lmev/(linelam(lineidx(i)))
     &     ,linespec(lineidx(i)),elem(lineid(lineidx(i),1)),
     &     rom(lineid(lineidx(i),2)),linekind(lineidx(i)),
     &     lineacc(lineidx(i))
        enddo
c
      endif
c
      return
      end
      subroutine heapindexsort (n, ra, idx)
c
c     heapsort (because of partial ordering already present)
c     return order in arrau idx, and ra is untouched.
c
      include 'const.inc'
c
      integer*4 n,idx(mxspeclines)
      real*8 ra(mxspeclines)
c
      integer*4 i,ir,j,l,itmp
      real*8  rra
c
      do i=1,n
        idx(i)=i
      enddo
c
      if (n.lt.2) return
c
      l=n/2+1
      ir=n
c
   10 if (l.gt.1) then
        l=l-1
        itmp=idx(l)
        rra=ra(itmp)
      else
        itmp=idx(ir)
        rra=ra(itmp)
        idx(ir)=idx(1)
        ir=ir-1
        if (ir.eq.1) then
          idx(1)=itmp
          return
        endif
      endif
c
      i=l
      j=l+l
   20 if (j.le.ir) then
        if (j.lt.ir) then
          if (ra(idx(j)).lt.ra(idx(j+1))) j=j+1
        endif
        if (rra.lt.ra(idx(j))) then
          idx(i)=idx(j)
          i=j
          j=j+j
        else
          j=ir+1
        endif
        goto 20
      endif
c
      idx(i)=itmp
c
      goto 10
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1979-2012 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     this subroutine is meant to be
c     used to write out various spectral components in a more
c     machine useable form (for plotting packages etc)
c     In general it will make looooonnnnggg lists
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine speclocal (luop, tl, el, eg, dl, t, dh, de, fh1, di,
     &dr, list, mode)
c
      include 'cblocks.inc'
c
c     Variables
c
      real*8 eloslog,fhbtlog,forblog
      real*8 tl,el,eg,dl
      real*8 t,dh,de,fh1
      real*8 di,dr
      real*8 chklim
c
c     Storage for sorted output.
c
      real*8 linelam(mxspeclines)
      real*8 linespec(mxspeclines)
      integer*4 lineid(mxspeclines,2)
      integer*4 lineidx(mxspeclines)
      integer*4 lineacc(mxspeclines)
      integer*4 linecount
c
      real*8 shortlam(mxspeclines)
      real*8 shortspec(mxspeclines)
      integer*4 shortlineid(mxspeclines,2)
      integer*4 shortidx(mxspeclines)
      integer*4 shortacc(mxspeclines)
      integer*4 shortcount
c
      integer*4 i,j,line,series
      integer*4 luop
c
      character mode*4,list*4, llist*4
c
c function
c
      real*8 fnair
c
c     Formats:
c
c
   10 format(' ',1pg12.6,1x,1pg12.6,1x,1pg12.4,1x,a3,a6,1x,i4)
   20 format(' ',a3,a6,1x,1pg12.6,1x,1pg12.6,1x,1pg12.4,1x,i4)
   30 format(/' Local Emission Line Spectrum File:'/
     &        ' =================================='/)
   40 format( ' MAPPINGS V version :',a8,/
     &        ' RUN:',a80,/)
   50 format( ' Local Slab Properties:'/
     &' ==================================================='/
     &'  TOTAL LOSS  :',1pg11.3/'  EFF. LOSS   :',1pg11.3/
     &'  EFF.  GAIN  :',1pg11.3/'  FRAC. RESID.:',1pg11.3/
     &'  TEMP.       :',1pg11.4/'  H DENSITY   :',1pg11.3/
     &'  ELECTR. DENS:',1pg11.3/'  FR. NEUT. H :',1pg11.3/
     &'  DISTANCE.   :',1pg11.3/'  SLAB DEPTH  :',1pg11.3/
     &' ===================================================')
   60 format(/' Local Flux Scale (H Beta = 1.000) :'/
     &        ' ============================='/
     &' H-beta  :',f8.4,
     &' (log(ergs/cm^2/s))'/
     &' TOTLOSS :',f8.4)
   70 format(/' Local Flux Scale (H Beta = 1.000) :'/
     &        ' ============================='/
     &' H-beta :',f8.4,
     &' (log(ergs/s))'/
     &' TOTLOSS :',f8.4)
   80 format(/' Local Flux Scale (Absolute units) :'/
     &        ' ============================='/
     &        ' Units (Log) :',f8.4)
   90 format(//' Local Line List >',1pg9.2,' :'/
     &  ' Number of lines:',i4,/
     &  ' ==========================='//
     &  ' Lambda(A)     E(eV)         Flux        Species    Accuracy'/
     &  ' ===========================================================')
  100 format(//' Flux Sorted Local Line List >',1pg9.2,' :'/
     &  ' Number of lines:',i4,/
     &  ' ======================================='//
     &  ' Lambda(A)     E(eV)         Flux        Species    Accuracy'/
     &  ' ===========================================================')
  110 format(//' Strongest Local Lines :'/
     &  ' Number of lines:',i4,/
     &  ' ======================='//
     &  ' Species     Lambda(A)      E(eV)        Flux       Accuracy'/
     &  ' ===========================================================')
c     Header
c
      write (luop,30)
      write (luop,40) runname,theversion
      write (luop,50) tl,el,eg,dl,t,dh,de,fh1,di,dr
c
c     chklim is a lower limit for individual line fluxes to print out.
c     1e-4 saves a lot of paper
c
c      chklim = epsilon
c
      chklim=1.0d-4*(hbeta+epsilon)
      if (mode.eq.'ABS') then
        chklim=epsilon
      endif
c
c     Copy list into llist (local list) so it can be modified if
c     necessary.
c
      llist=list
      if ((llist.ne.'TTWN').and.(llist.ne.'FLUX').and.(llist.ne.'LAMB')
     &.and.(llist.ne.'ALL')) then
        write (*,*) ' Mode error in SPEC2 : ',list
        write (*,*) ' Using lambda sorting.'
        llist='LAMB'
      endif
c
      linecount=0
c
      do series=1,nhseries
        do line=1,nhlines
          if (hydrobri(line,series).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=hlambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=hydrobri(line,series)/(hbeta+epsilon)
            lineid(linecount,1)=1
            lineid(linecount,2)=1
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          if (helibri(line,series).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=helambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=helibri(line,series)/(hbeta+epsilon)
            lineid(linecount,1)=2
            lineid(linecount,2)=2
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do i=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            if (xhydrobri(line,series,i).gt.chklim) then
              linecount=linecount+1
              linelam(linecount)=xhlambda(line,series,i)
              linelam(linecount)=linelam(linecount)/
     &         fnair(linelam(linecount))
              linespec(linecount)=xhydrobri(line,series,i)/(hbeta+
     &         epsilon)
              lineid(linecount,1)=i
              lineid(linecount,2)=mapz(i)
              lineacc(linecount)=2
            endif
          enddo
        enddo
      enddo
c
      do i=1,nfmions
        do j=1,nfmtrans(i)
          if (fluxm(j,i).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(fmlam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fluxm(j,i)
            lineid(linecount,1)=fmatom(i)
            lineid(linecount,2)=fmion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nfeions
        do j=1,nfetrans(i)
          if (fluxfe(j,i).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(felam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fluxfe(j,i)
            lineid(linecount,1)=featom(i)
            lineid(linecount,2)=feion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          if (f3bri(j,i).gt.chklim) then
            linecount=linecount+1
            linelam(linecount)=(f3lam(j,i)*1.d4)!nocorrection
            linespec(linecount)=f3bri(j,i)/(hbeta+epsilon)
            lineid(linecount,1)=f3atom(i)
            lineid(linecount,2)=f3ion(i)
            lineacc(linecount)=3
          endif
        enddo
      enddo
c
      do i=1,mlines
        if (fsbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(fslam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fsbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=ielfs(i)
          lineid(linecount,2)=ionfs(i)
          lineacc(linecount)=3
c
c   Special case more accurate lines: CII, NeII, SIV
c
c         if ((xrat(i).eq.zmap(6)).and.(xion(i).eq.2)) then
c           lineacc(linecount)=2
c         endif
c         if ((xrat(i).eq.zmap(10)).and.(xion(i).eq.2)) then
c           lineacc(linecount)=2
c         endif
c         if ((xrat(i).eq.zmap(16)).and.(xion(i).eq.4)) then
c           lineacc(linecount)=2
c         endif
        endif
      enddo
c
c     do i=1,xilines
c       if (xibri(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xilam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=xibri(i)/(hbeta+epsilon)
c         lineid(linecount,1)=xiat(i)
c         lineid(linecount,2)=xiion(i)
c         lineacc(linecount)=5
c       endif
c     enddo
c
c     do i=1,xhelines
c       if (xhebri(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xhelam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=xhebri(i)/(hbeta+epsilon)
c         lineid(linecount,1)=xheat(i)
c         lineid(linecount,2)=xheion(i)
c         lineacc(linecount)=5
c       endif
c     enddo
c
      do i=1,nlines
        if (rbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rlam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=ielr(i)
          lineid(linecount,2)=ionr(i)
          lineacc(linecount)=4
        endif
      enddo
c
c     do i=1,xlines
c       if (xrbri(i).gt.chklim) then
c         linecount=linecount+1
c         linelam(linecount)=xrlam(i)
c         linelam(linecount)=linelam(linecount)
c    &                        /fnair(linelam(linecount))
c         linespec(linecount)=xrbri(i)/(hbeta+epsilon)
c         lineid(linecount,1)=xrat(i)
c         lineid(linecount,2)=xion(i)
c         lineacc(linecount)=4
c       endif
c     enddo
c
      do i=1,nxr3lines
        if (xr3lines_bri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=xr3lines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xr3lines_bri(i)/(hbeta+epsilon)
          lineid(linecount,1)=xr3lines_at(i)
          lineid(linecount,2)=xr3lines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nxrllines
        if (xrllines_bri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=xrllines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xrllines_bri(i)/(hbeta+epsilon)
          lineid(linecount,1)=xrllines_at(i)
          lineid(linecount,2)=xrllines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
c Original He I lines, just the first two
c
      do i=1,2
        if (heibri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heilam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heibri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
c New He I singlet lines
c
      do i=1,nheislines
        if (heisbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heislam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heisbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
c New He I triplet lines
c
      do i=1,nheitlines
        if (heitbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heitlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heitbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nrccii
        if (rccii_bbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rccii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rccii_bbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(6)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcnii
        if (rcnii_bbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcnii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcnii_bbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(7)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_q
        if (rcoi_qbbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_qlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_qbbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_t
        if (rcoi_tbbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_tlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_tbbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
      do i=1,nrcoii
        if (rcoii_bbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoii_bbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcneii
        if (rcneii_bbri(i).gt.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcneii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcneii_bbri(i)/(hbeta+epsilon)
          lineid(linecount,1)=zmap(10)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      if ((llist.eq.'LAMB').or.(llist.eq.'ALL')) then
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10((hbeta*fpi)+epsilon)+vunilog
          forblog=dlog10(fmloss+feloss+epsilon)+vunilog
          eloslog=dlog10(eloss+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,70) fhbtlog,eloslog
          else
            write (luop,60) fhbtlog,eloslog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(hbeta+epsilon)+vunilog
          write (luop,80) fhbtlog
c
        endif
c
        write (luop,90) chklim,linecount
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Spectrum output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
        call heapindexsort (linecount, linelam, lineidx)
c
        do i=1,linecount
          write (luop,10) linelam(lineidx(i)),lmev/(linelam(lineidx(i)))
     &     ,linespec(lineidx(i)),elem(lineid(lineidx(i),1)),
     &     rom(lineid(lineidx(i),2)),lineacc(lineidx(i))
        enddo
c
      endif
c
      if ((llist.eq.'FLUX').or.(llist.eq.'ALL')) then
c
c
        if (mode.eq.'REL') then
c
          fhbtlog=dlog10((hbeta*(4.d0*pi))+epsilon)+vunilog
          forblog=dlog10(fmloss+feloss+epsilon)+vunilog
          eloslog=dlog10(eloss+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
            write (luop,70) fhbtlog,eloslog
          else
            write (luop,60) fhbtlog,eloslog
          endif
c
        endif
c
        if (mode.eq.'ABS') then
c
          fhbtlog=dlog10(hbeta+epsilon)+vunilog
          write (luop,80) fhbtlog
c
        endif
c
        write (luop,100) chklim,linecount
c
        if (linecount.gt.mxspeclines) then
          write (*,*) 'ERROR: Speclocal output has too many lines'
          write (*,*) linecount,' exceeds',mxspeclines,' limit.'
          write (*,*) 'Edit const.inc, increase mxspeclines parameter'
          write (*,*) 'and rebuild mappings.'
          stop
        endif
c
        call heapindexsort (linecount, linespec, lineidx)
c
        do i=linecount,1,-1
          write (luop,10) linelam(lineidx(i)),lmev/(linelam(lineidx(i)))
     &     ,linespec(lineidx(i)),elem(lineid(lineidx(i),1)),
     &     rom(lineid(lineidx(i),2)),lineacc(lineidx(i))
        enddo
c
      endif
c
c
      if (llist.eq.'TTWN') then
c
        write (luop,110) linecount
c
        call heapindexsort (linecount, linespec, lineidx)
c
        j=1
        if (linecount.gt.50) j=linecount-49
c
c     Copy top 50 lines into short array for lambda sorting
c
        shortcount=linecount-j+1
        do i=j,linecount
          shortlineid(i-j+1,1)=lineid(lineidx(i),1)
          shortlineid(i-j+1,2)=lineid(lineidx(i),2)
          shortlam(i-j+1)=linelam(lineidx(i))
          shortspec(i-j+1)=linespec(lineidx(i))
          shortacc(i-j+1)=lineacc(lineidx(i))
        enddo
c
        call heapindexsort (shortcount, shortlam, shortidx)
c
        do i=1,shortcount
          write (luop,20) elem(shortlineid(shortidx(i),1)),
     &     rom(shortlineid(shortidx(i),2)),shortlam(shortidx(i)),lmev/
     &     (shortlam(shortidx(i))),shortspec(shortidx(i)),
     &     shortacc(shortidx(i))
        enddo
c
      endif
c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1979-2012 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine speclocallines (fluxes)
c
      include 'cblocks.inc'
c
c
      real*8 chklim
c
c     Storage for sorted output.
c
      real*8 fluxes(mxmonlines)
      real*8 delta
c
      real*8 linelam(mxspeclines)
      real*8 linespec(mxspeclines)
      integer*4 lineid(mxspeclines,2)
      integer*4 lineidx(mxspeclines)
      integer*4 lineacc(mxspeclines)
      integer*4 linecount
c
      integer*4 i,j,line,series
c      integer*4 luop
c
c      character mode*4,list*4, llist*4
c
c function
c
      real*8 fnair
c
      chklim = epsilon
c
c       chklim=0.d0
c
      do i=1,mxmonlines
          fluxes(i)=0.d0
      enddo
c
      if (njlines.le.0) return
c
      linecount=0
c
      do series=1,nhseries
        do line=1,nhlines
          if (hydrobri(line,series).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=hlambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=hydrobri(line,series)
            lineid(linecount,1)=1
            lineid(linecount,2)=1
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          if (helibri(line,series).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=helambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=helibri(line,series)
            lineid(linecount,1)=2
            lineid(linecount,2)=2
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do i=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            if (xhydrobri(line,series,i).ge.chklim) then
              linecount=linecount+1
              linelam(linecount)=xhlambda(line,series,i)
              linelam(linecount)=linelam(linecount)/
     &         fnair(linelam(linecount))
              linespec(linecount)=xhydrobri(line,series,i)
              lineid(linecount,1)=i
              lineid(linecount,2)=mapz(i)
              lineacc(linecount)=2
            endif
          enddo
        enddo
      enddo
c
      do i=1,nfmions
        do j=1,nfmtrans(i)
          if (fmbri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(fmlam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fmbri(j,i)
            lineid(linecount,1)=fmatom(i)
            lineid(linecount,2)=fmion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nfeions
        do j=1,nfetrans(i)
          if (febri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(felam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=febri(j,i)
            lineid(linecount,1)=featom(i)
            lineid(linecount,2)=feion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          if (f3bri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(f3lam(j,i)*1.d4)!nocorrection
            linespec(linecount)=f3bri(j,i)
            lineid(linecount,1)=f3atom(i)
            lineid(linecount,2)=f3ion(i)
            lineacc(linecount)=3
          endif
        enddo
      enddo
c
      do i=1,mlines
        if (fsbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(fslam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fsbri(i)
          lineid(linecount,1)=ielfs(i)
          lineid(linecount,2)=ionfs(i)
          lineacc(linecount)=3

        endif
      enddo
c
      do i=1,nlines
        if (rbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rlam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rbri(i)
          lineid(linecount,1)=ielr(i)
          lineid(linecount,2)=ionr(i)
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nxr3lines
        if (xr3lines_bri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=xr3lines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xr3lines_bri(i)
          lineid(linecount,1)=xr3lines_at(i)
          lineid(linecount,2)=xr3lines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nxrllines
        if (xrllines_bri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=xrllines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xrllines_bri(i)
          lineid(linecount,1)=xrllines_at(i)
          lineid(linecount,2)=xrllines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
c Original He I lines, just the first two
c
      do i=1,2
        if (heibri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heilam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heibri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
c New He I singlet lines
c
      do i=1,nheislines
        if (heisbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heislam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heisbri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
c New He I triplet lines
c
      do i=1,nheitlines
        if (heitbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heitlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heitbri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nrccii
        if (rccii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rccii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rccii_bbri(i)
          lineid(linecount,1)=zmap(6)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcnii
        if (rcnii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcnii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcnii_bbri(i)
          lineid(linecount,1)=zmap(7)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_q
        if (rcoi_qbbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_qlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_qbbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_t
        if (rcoi_tbbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_tlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_tbbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
      do i=1,nrcoii
        if (rcoii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoii_bbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcneii
        if (rcneii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcneii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcneii_bbri(i)
          lineid(linecount,1)=zmap(10)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
c      write (*,*) chklim,linecount
c
      if (linecount.gt.mxspeclines) then
        write (*,*) 'ERROR: Spectrum output has too many lines'
        write (*,*) linecount,' exceeds',mxspeclines,' limit.'
        write (*,*) 'Edit const.inc, increase mxspeclines parameter'
        write (*,*) 'and rebuild mappings.'
        stop
      endif
c
      call heapindexsort (linecount, linelam, lineidx)
c
      do i=1,linecount
        do j=1,njlines
           delta=dabs(linelam(lineidx(i))-emlinlist(j))
           if (delta.le.emlindeltas(j)) then
               fluxes(j)=fluxes(j)+linespec(lineidx(i))
           endif
        enddo
      enddo
c
      return
      end

      subroutine speclocallineids (lineat,lineion)
c
      include 'cblocks.inc'
c
c
      real*8 chklim
c
c     Storage for sorted output.
c
      integer*4 lineat(mxmonlines)
      integer*4 lineion(mxmonlines)
      real*8 delta
c
      real*8 linelam(mxspeclines)
      real*8 linespec(mxspeclines)
      integer*4 lineid(mxspeclines,2)
      integer*4 lineidx(mxspeclines)
      integer*4 lineacc(mxspeclines)
      integer*4 linecount
c
      integer*4 i,j,line,series
c      integer*4 luop
c
c      character mode*4,list*4, llist*4
c
c function
c
      real*8 fnair
c
      chklim = 0.0d0
c
c       chklim=0.d0
c
      do i=1,mxmonlines
          lineat(i)=1
          lineion(i)=1
      enddo
c
      if (njlines.le.0) return
c
      linecount=0
c
      do series=1,nhseries
        do line=1,nhlines
          if (hydrobri(line,series).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=hlambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=hydrobri(line,series)
            lineid(linecount,1)=1
            lineid(linecount,2)=1
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          if (helibri(line,series).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=helambda(line,series)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=helibri(line,series)
            lineid(linecount,1)=2
            lineid(linecount,2)=2
            lineacc(linecount)=1
            if (series.eq.1) then
              lineacc(linecount)=3
            endif
          endif
        enddo
      enddo
c
      do i=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            if (xhydrobri(line,series,i).ge.chklim) then
              linecount=linecount+1
              linelam(linecount)=xhlambda(line,series,i)
              linelam(linecount)=linelam(linecount)/
     &         fnair(linelam(linecount))
              linespec(linecount)=xhydrobri(line,series,i)
              lineid(linecount,1)=i
              lineid(linecount,2)=mapz(i)
              lineacc(linecount)=2
            endif
          enddo
        enddo
      enddo
c
      do i=1,nfmions
        do j=1,nfmtrans(i)
          if (fmbri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(fmlam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=fmbri(j,i)
            lineid(linecount,1)=fmatom(i)
            lineid(linecount,2)=fmion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nfeions
        do j=1,nfetrans(i)
          if (febri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(felam(j,i)*1.d8)
            linelam(linecount)=linelam(linecount)/
     &       fnair(linelam(linecount))
            linespec(linecount)=febri(j,i)
            lineid(linecount,1)=featom(i)
            lineid(linecount,2)=feion(i)
            lineacc(linecount)=2
          endif
        enddo
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          if (f3bri(j,i).ge.chklim) then
            linecount=linecount+1
            linelam(linecount)=(f3lam(j,i)*1.d4)!nocorrection
            linespec(linecount)=f3bri(j,i)
            lineid(linecount,1)=f3atom(i)
            lineid(linecount,2)=f3ion(i)
            lineacc(linecount)=3
          endif
        enddo
      enddo
c
      do i=1,mlines
        if (fsbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(fslam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=fsbri(i)
          lineid(linecount,1)=ielfs(i)
          lineid(linecount,2)=ionfs(i)
          lineacc(linecount)=3

        endif
      enddo
c
      do i=1,nlines
        if (rbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rlam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rbri(i)
          lineid(linecount,1)=ielr(i)
          lineid(linecount,2)=ionr(i)
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nxr3lines
        if (xr3lines_bri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=xr3lines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xr3lines_bri(i)
          lineid(linecount,1)=xr3lines_at(i)
          lineid(linecount,2)=xr3lines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nxrllines
        if (xrllines_bri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=xrllines_lam(i)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=xrllines_bri(i)
          lineid(linecount,1)=xrllines_at(i)
          lineid(linecount,2)=xrllines_ion(i)
          lineacc(linecount)=3
        endif
      enddo
c
c Original He I lines, just the first two
c
      do i=1,2
        if (heibri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heilam(i)*1.d8)
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heibri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
c New He I singlet lines
c
      do i=1,nheislines
        if (heisbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heislam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heisbri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
c New He I triplet lines
c
      do i=1,nheitlines
        if (heitbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(heitlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=heitbri(i)
          lineid(linecount,1)=zmap(2)
          lineid(linecount,2)=1
          lineacc(linecount)=3
        endif
      enddo
c
      do i=1,nrccii
        if (rccii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rccii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rccii_bbri(i)
          lineid(linecount,1)=zmap(6)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcnii
        if (rcnii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcnii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcnii_bbri(i)
          lineid(linecount,1)=zmap(7)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_q
        if (rcoi_qbbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_qlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_qbbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcoi_t
        if (rcoi_tbbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoi_tlam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoi_tbbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=1
          lineacc(linecount)=5
        endif
      enddo
c
      do i=1,nrcoii
        if (rcoii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcoii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcoii_bbri(i)
          lineid(linecount,1)=zmap(8)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
      do i=1,nrcneii
        if (rcneii_bbri(i).ge.chklim) then
          linecount=linecount+1
          linelam(linecount)=(rcneii_lam(i))
          linelam(linecount)=linelam(linecount)/fnair(linelam(linecount)
     &     )
          linespec(linecount)=rcneii_bbri(i)
          lineid(linecount,1)=zmap(10)
          lineid(linecount,2)=2
          lineacc(linecount)=4
        endif
      enddo
c
c      write (*,*) chklim,linecount
c
      if (linecount.gt.mxspeclines) then
        write (*,*) 'ERROR: Spectrum output has too many lines'
        write (*,*) linecount,' exceeds',mxspeclines,' limit.'
        write (*,*) 'Edit const.inc, increase mxspeclines parameter'
        write (*,*) 'and rebuild mappings.'
        stop
      endif
c
      call heapindexsort (linecount, linelam, lineidx)
c
      do i=1,linecount
        do j=1,njlines
           delta=dabs(linelam(lineidx(i))-emlinlist(j))
           if (delta.le.emlindeltas(j)) then
               lineat(j)=lineid(lineidx(i),1)
               lineion(j)=lineid(lineidx(i),2)
           endif
        enddo
      enddo
c
      return
      end
