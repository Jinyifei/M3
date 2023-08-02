cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recom_neii (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      integer*4 idx,atom,j,nte
      real*8 t4,tm,logt4,abde,emiss,tr
      real*8 x(8),y(8),y2(8)
      real*8 neiite(8)
c
      real*8 fsplint
c
      data neiite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
c
      tr=t
      if (usekappa) then
        tr=t*(kappa-1.5d0)/kappa
      endif
      t4=tr*1.0d-4
      if (t4.le.0.001d0) t4=0.001d0
      logt4=dlog10(t4)
      tm=1.d0-t4
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'NeII recombination', nrcNeII
c effective recombination rates per line
c Case A and CaseB, 6 brightest multiplets
c
      atom=zmap(10)
c
      do idx=1,nrcneii
        rcneii_abri(idx)=0.d0
        rcneii_bbri(idx)=0.d0
      enddo
c
c check NeIII population, that is doing the recomb
c
      if (zion(atom)*pop(3,atom).gt.pzlimit) then
c
        nte=8
        do j=1,nte
          x(j)=dlog10(neiite(j))-4.d0
        enddo
c
        abde=de*dh*zion(atom)*pop(3,atom)
c
c All quartets and doublets together
c
        do idx=1,nrcneii
          emiss=0.d0
          do j=1,nte
            y(j)=rcneii_y(1,j,idx)
            y2(j)=rcneii_y2(1,j,idx)
          enddo
          emiss=10.d0**(fsplint(x,y,y2,nte,logt4)-14.d0)
          emiss=emiss*rcneii_eij(idx)
          rcneii_abri(idx)=abde*emiss*ifpi
          emiss=0.d0
          do j=1,nte
            y(j)=rcneii_y(2,j,idx)
            y2(j)=rcneii_y2(2,j,idx)
          enddo
          emiss=10.d0**(fsplint(x,y,y2,nte,logt4)-14.d0)
          emiss=emiss*rcneii_eij(idx)
          rcneii_bbri(idx)=abde*emiss*ifpi
        enddo
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recom_cii (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      integer*4 idx,atom,j,nte
      real*8 logte,abde,emiss,tr
      real*8 x(15),y(15),y2(15)
c
      real*8 fsplint
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'CII recombination', nrcCII
c
      atom=zmap(6)
      nte=15
c
      do idx=1,nrccii
        rccii_abri(idx)=0.d0
        rccii_bbri(idx)=0.d0
      enddo
c
c check CIII population, that is doing the recomb
c
      if (zion(atom)*pop(3,atom).gt.pzlimit) then
c
        tr=t
        if (usekappa) then
          tr=t*(kappa-1.5d0)/kappa
        endif
        logte=dlog10(tr)
c
        abde=de*dh*zion(atom)*pop(3,atom)
c
c All lines together
c
        do j=1,nte
          x(j)=rccii_logte(j)
        enddo
c
        do idx=1,nrccii
c
c spline interpolate in log electron density - log emiss space
c
c Case A
c
          do j=1,nte
            y(j)=rccii_coeffy(1,j,idx)
            y2(j)=rccii_coeffy2(1,j,idx)
          enddo
          emiss=10.d0**(fsplint(x,y,y2,nte,logte)-14.d0)
          emiss=emiss*rccii_eij(idx)
          rccii_abri(idx)=abde*emiss*ifpi
c
c Case B
c
          do j=1,nte
            y(j)=rccii_coeffy(2,j,idx)
            y2(j)=rccii_coeffy2(2,j,idx)
          enddo
          emiss=10.d0**(fsplint(x,y,y2,nte,logte)-14.d0)
          emiss=emiss*rccii_eij(idx)
          rccii_bbri(idx)=abde*emiss*ifpi
c
        enddo
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recom_nii (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      integer*4 idx,atom,ni,j,nte
      real*8 tr,t4,logt4,abde,emiss,logde
      real*8 x(4),y(4),y2(4)
      real*8 niite(8)
      real*8 xdens(8),ydens(8),ydens2(8)
c
      real*8 fsplint
c
      data x/2.0d0,3.0d0,4.0d0,5.0d0/
      data niite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'NII recombination', nrcNII
c
c     log({alpha}_eff_)+15 = a+b*t+c*t^2^+(d+e*t+f*t^2^)log(t)+g*(log(t))^2^+h/t
c
      atom=zmap(7)
c
      do idx=1,nrcnii
        rcnii_abri(idx)=0.d0
        rcnii_bbri(idx)=0.d0
      enddo
c
c check NIII population, that is doing the recomb
c
      if (zion(atom)*pop(3,atom).gt.pzlimit) then
c
        nte=8
        do j=1,nte
          xdens(j)=dlog10(niite(j))-4.d0
        enddo
        tr=t
        if (usekappa) then
          tr=t*(kappa-1.5d0)/kappa
        endif
        t4=tr*1.0d-4
        if (t4.le.0.001d0) t4=0.001d0
        logt4=dlog10(t4)
        logde=dlog10(de)
c
        abde=de*dh*zion(atom)*pop(3,atom)
c
c All lines together
c
        do idx=1,nrcnii
c
c spline interpolate in log electron density - log emiss space
c
          do j=1,nte
            ydens(j)=rcniin20_y(1,j,idx)
            ydens2(j)=rcniin20_y2(1,j,idx)
          enddo
          y(1)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin30_y(1,j,idx)
            ydens2(j)=rcniin30_y2(1,j,idx)
          enddo
          y(2)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin40_y(1,j,idx)
            ydens2(j)=rcniin40_y2(1,j,idx)
          enddo
          y(3)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin50_y(1,j,idx)
            ydens2(j)=rcniin50_y2(1,j,idx)
          enddo
          y(4)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
c
          ni=4
c
c        call endgrads (x,y,ni,dydx0,dydxn)
          call spline (x, y, ni, 1.0d30, 1.0d30, y2)
          emiss=10.d0**(fsplint(x,y,y2,ni,logde)-15.d0)
          emiss=emiss*rcnii_eij(idx)
          rcnii_abri(idx)=abde*emiss*ifpi
c        write(*,*)  rcNII_Abri(idx)
c
          do j=1,nte
            ydens(j)=rcniin20_y(2,j,idx)
            ydens2(j)=rcniin20_y2(2,j,idx)
          enddo
          y(1)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin30_y(2,j,idx)
            ydens2(j)=rcniin30_y2(2,j,idx)
          enddo
          y(2)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin40_y(2,j,idx)
            ydens2(j)=rcniin40_y2(2,j,idx)
          enddo
          y(3)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcniin50_y(2,j,idx)
            ydens2(j)=rcniin50_y2(2,j,idx)
          enddo
          y(4)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
c        call endgrads (x,y,ni,dydx0,dydxn)
          call spline (x, y, ni, 1.0d30, 1.0d30, y2)
          emiss=10.d0**(fsplint(x,y,y2,ni,logde)-15.d0)
          emiss=emiss*rcnii_eij(idx)
          rcnii_bbri(idx)=abde*emiss*ifpi
c
c        if ( idx .eq. 8 ) then
c        write(*,*) 'NII Recomb: B'
c        write(*,'(4(1pg14.7))') zion(atom),pop(3,atom),rcNII_Bbri(idx)
c        write(*,'(4(1pg14.7))') logDe,t4,invT4,logT4
c        write(*,'(4(1pg14.7))') (x(i), i=1,4)
c        write(*,'(4(1pg14.7))') (y(i), i=1,4)
c        write(*,'(4(1pg14.7))') fsplint(x,y,y2,ni,logDe),rcNII_eij(idx),
c     &         abde, emiss
c        endif
c
        enddo
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
c      real*8 function fniiemiss(t4,invT4,logT4,coef,idx,caseAB)
c      implicit none
c      real*8 t4,invT4,logT4
c      real*8 coef(2,8,55)
c      integer*4 idx,caseAB
c      real*8 a,b,c,d,e,f,g,h
c      real*8 emiss1,emiss2,emiss3,emiss
c      a=coef(caseAB,1,idx)
c      b=coef(caseAB,2,idx)
c      c=coef(caseAB,3,idx)
c      d=coef(caseAB,4,idx)
c      e=coef(caseAB,5,idx)
c      f=coef(caseAB,6,idx)
c      g=coef(caseAB,7,idx)
c      h=coef(caseAB,8,idx)
c      emiss1=a+t4*(b+t4*c)
c      emiss2=(d+t4*(e+t4*f))*logT4
c      emiss3=(g*logT4*logT4)+h*invT4
c      emiss=(emiss1+emiss2+emiss3)-15.d0
c      if (emiss.lt.-20.d0) emiss=-20.d0
c      if (emiss.gt.-5.d0) emiss=-5.d0
c      fniiemiss=emiss
c      return
c      end
c
      subroutine recom_oi (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      integer*4 idx,atom
      real*8 tr,t4,lnt4,abde,emiss
      real*8 a,b,c
c
c
      tr=t
      if (usekappa) then
        tr=t*(kappa-1.5d0)/kappa
      endif
      t4=tr*1.0d-4
      if (t4.le.0.001d0) t4=0.001d0
      lnt4=dlog(t4)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'OI recombination', nrcOI_Q, nrcOI_T
c effective recombination rates per line
c Case A and CaseB, Quintets and Triplets
c
      atom=zmap(8)
c
      do idx=1,nrcoi_q
        rcoi_qabri(idx)=0.d0
        rcoi_qbbri(idx)=0.d0
      enddo
c
      do idx=1,nrcoi_t
        rcoi_tabri(idx)=0.d0
        rcoi_tbbri(idx)=0.d0
      enddo
c
c check OII population, that is doing the recomb
c
      if (zion(atom)*pop(2,atom).gt.pzlimit) then
c
        abde=de*dh*zion(atom)*pop(2,atom)
c
c Quintets....
c
        do idx=1,nrcoi_q
          a=rcoi_qcoef(1,1,idx)
          b=rcoi_qcoef(1,2,idx)
          c=rcoi_qcoef(1,3,idx)
          emiss=(a*(t4**(-b*(1.d0+c*lnt4))))
          emiss=emiss*rcoi_qeij(idx)
          rcoi_qabri(idx)=abde*emiss*ifpi
          a=rcoi_qcoef(2,1,idx)
          b=rcoi_qcoef(2,2,idx)
          c=rcoi_qcoef(2,3,idx)
          emiss=(a*(t4**(-b*(1.d0+c*lnt4))))
          emiss=dmax1(0.d0,emiss)
          emiss=emiss*rcoi_qeij(idx)
          rcoi_qbbri(idx)=abde*emiss*ifpi
        enddo
c
c Triplets....
c
        do idx=1,nrcoi_t
          a=rcoi_tcoef(1,1,idx)
          b=rcoi_tcoef(1,2,idx)
          c=rcoi_tcoef(1,3,idx)
          emiss=(a*(t4**(-b*(1.d0+c*lnt4))))
          emiss=emiss*rcoi_teij(idx)
          rcoi_tabri(idx)=abde*emiss*ifpi
          a=rcoi_tcoef(2,1,idx)
          b=rcoi_tcoef(2,2,idx)
          c=rcoi_tcoef(2,3,idx)
          emiss=(a*(t4**(-b*(1.d0+c*lnt4))))
          emiss=dmax1(0.d0,emiss)
          emiss=emiss*rcoi_teij(idx)
          rcoi_tbbri(idx)=abde*emiss*ifpi
        enddo
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recom_oii (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      integer*4 idx,atom,ni,j,nte
      real*8 tr,t4,tm,logt4,abde,emiss,logde
      real*8 x(4),y(4),y2(4)
      real*8 oiite(8)
      real*8 xdens(8),ydens(8),ydens2(8)
c
      real*8 fsplint
c
      data x/2.0d0,4.0d0,5.0d0,6.0d0/
      data oiite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
c
      tr=t
      if (usekappa) then
        tr=t*(kappa-1.5d0)/kappa
      endif
      t4=tr*1.0d-4
      if (t4.le.0.001d0) t4=0.001d0
      logt4=dlog10(t4)
      logde=dlog10(de)
      tm=1.d0-t4
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c effective recombination rates per line
c Case A, Case B, Case C, 35 Mutiplets
c
      atom=zmap(8)
c
      do idx=1,nrcoii
        rcoii_abri(idx)=0.d0
        rcoii_bbri(idx)=0.d0
        rcoii_cbri(idx)=0.d0
      enddo
c
c check OIII population, that is doing the recomb
c
      if (zion(atom)*pop(3,atom).gt.pzlimit) then
        nte=8
        do j=1,nte
          xdens(j)=dlog10(oiite(j))-4.d0
        enddo
        abde=de*dh*zion(atom)*pop(3,atom)
        do idx=1,nrcoii
          emiss=0.d0
          do j=1,nte
            ydens(j)=rcoiin20_y(1,j,idx)
            ydens2(j)=rcoiin20_y2(1,j,idx)
          enddo
          y(1)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin40_y(1,j,idx)
            ydens2(j)=rcoiin40_y2(1,j,idx)
          enddo
          y(2)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin50_y(1,j,idx)
            ydens2(j)=rcoiin50_y2(1,j,idx)
          enddo
          y(3)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin60_y(1,j,idx)
            ydens2(j)=rcoiin60_y2(1,j,idx)
          enddo
          y(4)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          ni=4
c
          call spline (x, y, ni, 1.0d30, 1.0d30, y2)
          emiss=10.d0**(fsplint(x,y,y2,ni,logde)-14.d0)
          emiss=emiss*rcoii_eij(idx)
          rcoii_abri(idx)=abde*emiss*ifpi
c caseB
          emiss=0.d0
          do j=1,nte
            ydens(j)=rcoiin20_y(2,j,idx)
            ydens2(j)=rcoiin20_y2(2,j,idx)
          enddo
          y(1)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin40_y(2,j,idx)
            ydens2(j)=rcoiin40_y2(2,j,idx)
          enddo
          y(2)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin50_y(2,j,idx)
            ydens2(j)=rcoiin50_y2(2,j,idx)
          enddo
          y(3)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin60_y(2,j,idx)
            ydens2(j)=rcoiin60_y2(2,j,idx)
          enddo
          y(4)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          ni=4
c
          call spline (x, y, ni, 1.0d30, 1.0d30, y2)
          emiss=10.d0**(fsplint(x,y,y2,ni,logde)-14.d0)
          emiss=emiss*rcoii_eij(idx)
          rcoii_bbri(idx)=abde*emiss*ifpi
c case C
          emiss=0.d0
          do j=1,nte
            ydens(j)=rcoiin20_y(3,j,idx)
            ydens2(j)=rcoiin20_y2(3,j,idx)
          enddo
          y(1)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin40_y(3,j,idx)
            ydens2(j)=rcoiin40_y2(3,j,idx)
          enddo
          y(2)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin50_y(3,j,idx)
            ydens2(j)=rcoiin50_y2(3,j,idx)
          enddo
          y(3)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          do j=1,nte
            ydens(j)=rcoiin60_y(3,j,idx)
            ydens2(j)=rcoiin60_y2(3,j,idx)
          enddo
          y(4)=fsplint(xdens,ydens,ydens2,nte,logt4)
c
          ni=4
c
          call spline (x, y, ni, 1.0d30, 1.0d30, y2)
          emiss=10.d0**(fsplint(x,y,y2,ni,logde)-14.d0)
          emiss=emiss*rcoii_eij(idx)
          rcoii_cbri(idx)=abde*emiss*ifpi
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
