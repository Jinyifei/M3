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
c*******TO ZERO BUFFER ARRAYS,RESET COUNTERS AND
c     SET DEFAULT VALUES
c
c
      subroutine zer ()
c
      include 'cblocks.inc'
c
      integer*4 j, i, k, atom
c
      call zeroemiss
c
c    ***MAKE RECOMB. COEFFICIENTS EFFECTIVE FOR ALL ELEMENTS
c
      fhbeta=0.0d0
      do i=1,atypes
        do j=1,maxion(i)
          arad(j,i)=dabs(arad(j,i))
        enddo
      enddo
c
c    ***SET DEFAULT POP. AND ZERO THE WEIGHTED IONIC QUANTITIES
c
      tem=-1
      do i=1,atypes
        do j=1,mxion
          pop(j,i)=0.0d0
          teav(j,i)=0.0d0
          popint(j,i)=0.0d0
          deam(j,i)=0.0d0
          rdisa(j,i)=0.0d0
          pam(j,i)=0.0d0
        enddo
        pop(1,i)=1.0d0
      enddo
      deav=0.0d0
      dhav=0.0d0
      dhavv=0.0d0
      vunilog=0.0d0
      qtosoh=0.0d0
      qhdha=0.0d0
      zetaeav=0.0d0
      tlosac=0.0d0
      tforbi=0.0d0
      fheiilr=0.0d0
      weoiii=0.0d0
      wenii=0.0d0
      wesii=0.0d0
      weoii=0.0d0
      roiii=0.0d0
      rnii=0.0d0
      rsii=0.0d0
      roii=0.0d0
      deoiii=0.0d0
      denii=0.0d0
      desii=0.0d0
      deoii=0.0d0
      toiii=0.0d0
      tnii=0.0d0
      tsii=0.0d0
      toii=0.0d0
c
c
c     get ion entry numbers for forbidden OI,OIII,NII,SII,OII
c
      ox3=0
      ni2=0
      su2=0
      ox2=0
      ox1=0
c
      do i=1,nfmions
        if ((fmion(i).eq.3).and.(fmatom(i).eq.zmap(8))) ox3=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(7))) ni2=i
        if ((fmion(i).eq.1).and.(fmatom(i).eq.zmap(8))) ox1=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(16))) su2=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(8))) ox2=i
      enddo
c
c    ***ZERO ALL PHOTON FIELD VECTORS
c
      do 10 i=1,infph
        skipbin(i)=.true.
        soupho(i)=0.0d0
        dwdif(i)=0.0d0
        updif(i)=0.0d0
        emidif(i)=0.0d0
        dwdifcont(i)=0.0d0
        updifcont(i)=0.0d0
        emidifcont(i)=0.0d0
        tphot(i)=0.0d0
        irphot(i)=0.d0
   10 continue
c
c     zero grain heating
c
      gheat=0.d0
      gcool=0.d0
      paheat=0.d0
c
c     Dust Column
c
      if ((grainmode.eq.1).or.(graindestructmode.eq.1)) then
        dustint=0.d0
        do i=1,atypes
          dion(i)=dion0(i)
          zion(i)=zion0(i)*deltazion(i)*dion(i)
          invdion(i)=1.d0/dion(i)
        enddo
        if (pahmode.eq.1) then
          do i=1,pahi
            pahint(i)=0.d0
          enddo
          pahactive=0
          i=zmap(6)
          if (clinpah.eq.1) then
            zion(i)=zion0(i)*deltazion(i)*dion(i)
          else
            zion(i)=zion0(i)*deltazion(i)*(dion(i)+(1.d0-dion(i))*
     &       pahcfrac)
          endif
        endif
      endif
c
c
      qhi=0.0d0
      qhei=0.0d0
      qheii=0.0d0
      zstar=0.0d0
      alnth=0.0d0
      teff=0.0d0
      crate=0.d0
      cosphi=0.d0
      do k=1,numtypes
        avgpot(k)=0.d0
      enddo
c
c     flow globals
c
      vel0=0.d0
      vel1=0.d0
      rho0=0.d0
      rho1=0.d0
      pr0=0.d0
      pr1=0.d0
      dh0=0.d0
      dh1=0.d0
      de0=0.d0
      de1=0.d0
      te0=0.d0
      te1=0.d0
c
      cut=0.0d0
      iphom=-1
      ipho=0
c
c    ***SET DEFAULT VALUES FOR "ON THE SPOT APPROX" : JSPOT
c    ***FOR CASE A,B (0=CASE A  1=CASE B)  :  CASEAB
c    ***NUMBER OF ABSORBING ATOMIC ELEMENTS : LIMPH
c
      jspot='NO'
      jcon='YES'
      jspec='NO'
      do atom=1,atypes
        caseab(atom)=0.0d0
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
c
c*******TO ZERO BUFFER VECTORS ONLY
c
c
      subroutine zerbuf ()
c
      include 'cblocks.inc'
c
c
      integer*4 j, i
c
      call zeroemiss
c
      deav=0.0d0
      dhav=0.0d0
      dhavv=0.0d0
      tlosac=0.0d0
      tforbi=0.0d0
      fheiilr=0.0d0
c
c    ***ZERO IONIC WEIGHTED QUANTITIES
c
      do i=1,atypes
        do j=1,mxion
          pop(j,i)=0.0d0
          teav(j,i)=0.0d0
          popint(j,i)=0.0d0
          deam(j,i)=0.0d0
          rdisa(j,i)=0.0d0
          pam(j,i)=0.0d0
        enddo
        pop(1,i)=1.0d0
c         popint(maxion(i)+1,i) = 0.0d0
      enddo
c
      weoiii=0.0d0
      wenii=0.0d0
      wesii=0.0d0
      weoii=0.0d0
      roiii=0.0d0
      rnii=0.0d0
      rsii=0.0d0
      roii=0.0d0
      deoiii=0.0d0
      denii=0.0d0
      desii=0.0d0
      deoii=0.0d0
      toiii=0.0d0
      tnii=0.0d0
      tsii=0.0d0
      toii=0.0d0
c
c
c     get ion entry numbers for forbidden OI,OII,OIII,NII,SII
c
      ni2=0
      su2=0
      ox3=0
      ox2=0
      ox1=0
c
      do i=1,nfmions
        if ((fmion(i).eq.3).and.(fmatom(i).eq.zmap(8))) ox3=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(8))) ox2=i
        if ((fmion(i).eq.1).and.(fmatom(i).eq.zmap(8))) ox1=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(7))) ni2=i
        if ((fmion(i).eq.2).and.(fmatom(i).eq.zmap(16))) su2=i
      enddo
c
      return
c
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
      subroutine zeroemiss ()
c
      include 'cblocks.inc'
c
c
      integer*4 j, i, series, line
c
c    Clear all line and emission data only
c
      do i=1,nfmions
        do j=1,nfmtrans(i)
          fluxm(j,i)=0.0d0
          fmbri(j,i)=0.0d0
        enddo
      enddo
      do i=1,nfeions
        do j=1,nfetrans(i)
          fluxfe(j,i)=0.0d0
          febri(j,i)=0.0d0
        enddo
      enddo
      do i=1,10
        fluxh(i)=0.0d0
        hbri(i)=0.0d0
      enddo
      do i=1,nheilines
        fluxhei(i)=0.0d0
        heibri(i)=0.0d0
      enddo
      do i=1,nheislines
        fluxheis(i)=0.0d0
        heisbri(i)=0.0d0
      enddo
      do i=1,nheitlines
        fluxheit(i)=0.0d0
        heitbri(i)=0.0d0
      enddo
      heiioiiibfsum=0.d0
      heiioiiibf=0.d0
      do series=1,nhseries
        do line=1,nhlines
          hydroflux(line,series)=0.0d0
          hydrobri(line,series)=0.0d0
          hydlin(1,line,series)=0.0d0
          hydlin(2,line,series)=1.0d0
          hyduplin(1,line,series)=0.0d0
          hyduplin(2,line,series)=1.0d0
          hyddwlin(1,line,series)=0.0d0
          hyddwlin(2,line,series)=1.0d0
        enddo
      enddo
      do series=1,nheseries
        do line=1,nhelines
          heliflux(line,series)=0.0d0
          helibri(line,series)=0.0d0
          hellin(1,line,series)=0.0d0
          hellin(2,line,series)=1.0d0
          heluplin(1,line,series)=0.0d0
          heluplin(2,line,series)=1.0d0
          heldwlin(1,line,series)=0.0d0
          heldwlin(2,line,series)=1.0d0
        enddo
      enddo
c
      do i=1,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            xhydroflux(line,series,i)=0.d0
            xhydrobri(line,series,i)=0.0d0
            xhydrobri(line,series,i)=0.0d0
            xhydlin(1,line,series,i)=0.0d0
            xhydlin(2,line,series,i)=1.d0
            xhyduplin(1,line,series,i)=0.0d0
            xhyduplin(2,line,series,i)=1.0d0
            xhyddwlin(1,line,series,i)=0.0d0
            xhyddwlin(2,line,series,i)=1.0d0
          enddo
        enddo
      enddo
c
      do i=1,mlines
        fluxi(i)=0.0d0
        fsbri(i)=0.0d0
      enddo
c
      do i=1,nlines
        fluxr(i)=0.0d0
        rbri(i)=0.0d0
      enddo
c
      do i=1,nxr3lines
        xr3lines_flux(i)=0.0d0
        xr3lines_bri(i)=0.0d0
      enddo
c
      do i=1,nxrllines
        xrllines_flux(i)=0.0d0
        xrllines_bri(i)=0.0d0
      enddo
c
      do j=1,nrccii
        fluxrccii_a(j)=0.0d0
        fluxrccii_b(j)=0.0d0
        rccii_abri(j)=0.0d0
        rccii_bbri(j)=0.0d0
      enddo
c
      do j=1,nrcnii
        fluxrcnii_a(j)=0.d0
        fluxrcnii_b(j)=0.d0
        rcnii_abri(j)=0.0d0
        rcnii_bbri(j)=0.0d0
      enddo
c
      do j=1,nrcoi_q
        fluxrcoi_qa(j)=0.0d0
        fluxrcoi_qb(j)=0.0d0
        rcoi_qabri(j)=0.0d0
        rcoi_qbbri(j)=0.0d0
      enddo
      do j=1,nrcoi_t
        fluxrcoi_ta(j)=0.0d0
        fluxrcoi_tb(j)=0.0d0
        rcoi_tabri(j)=0.0d0
        rcoi_tbbri(j)=0.0d0
      enddo
c
      do j=1,nrcoii
        fluxrcoii_a(j)=0.0d0
        fluxrcoii_b(j)=0.0d0
        fluxrcoii_c(j)=0.0d0
        rcoii_abri(j)=0.0d0
        rcoii_bbri(j)=0.0d0
        rcoii_cbri(j)=0.0d0
      enddo
c
      do j=1,nrcneii
        fluxrcneii_a(j)=0.0d0
        fluxrcneii_b(j)=0.0d0
        rcneii_abri(j)=0.0d0
        rcneii_bbri(j)=0.0d0
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          fluxf3(j,i)=0.0d0
          f3bri(j,i)=0.0d0
        enddo
      enddo
c
      h2qav=0.0d0
      hei2qa=0.0d0
      heii2qa=0.0d0
c
      do i=1,atypes
        h2qflux(i)=0.0d0
      enddo
c
      fhbeta=0.0d0
      qtosoh=0.0d0
      qhdha=0.0d0
      zetaeav=0.0d0
c
c    ***ZERO PHOTON FIELD ARRAYS : UPDIF,DWDIF,UPLIN,DWLIN....
c     BUT NOT THE SOURCE VECTOR : SOUPHO
c
      do i=1,infph
        dwdif(i)=0.0d0
        updif(i)=0.0d0
        emidif(i)=0.0d0
        dwdifcont(i)=0.0d0
        updifcont(i)=0.0d0
        emidifcont(i)=0.0d0
        tphot(i)=0.0d0
        irphot(i)=0.d0
      enddo
c
      do i=1,nxr3lines
        xr3lines_emilin(1,i)=0.0d0
        xr3lines_emilin(2,i)=1.0d0
        xr3lines_uplin(1,i)=0.0d0
        xr3lines_uplin(2,i)=1.0d0
        xr3lines_dwlin(1,i)=0.0d0
        xr3lines_dwlin(2,i)=1.0d0
      enddo
c
      do i=1,nxrllines
        xrllines_emilin(1,i)=0.0d0
        xrllines_emilin(2,i)=1.0d0
        xrllines_uplin(1,i)=0.0d0
        xrllines_uplin(2,i)=1.0d0
        xrllines_dwlin(1,i)=0.0d0
        xrllines_dwlin(2,i)=1.0d0
      enddo
c
      h2qav=0.0d0
      hei2qa=0.0d0
      heii2qa=0.0d0
c
      do i=1,atypes
        h2qflux(i)=0.0d0
      enddo
c
      return
c
      end
