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
c*******TO INTEGRATE THE LINE AND TEMPERATURE DATA AND THE
c       COLUMN DENSITY
c       USES TEMPERATURE (T), EL. DENS. (DE), HYDR. DENS (DH),
c       FILLING FACT. (FI), THE VOLUME OF THE SLAB CONSIDERED (DVOLUNI)
c       (IN ARBITRARY UNITS (VUNIT)) AND FINALLY THE RADIUS STEP
c       SIZE (DR) (IN CM)
c       THE EQUIVALENT IN CC OF THE VOLUME UNIT : VUNIT  ARE IN
c       VUNILOG (LOG10)
c
c       IMOD = 'ALL'  : ADS ALL QUANTITIES IN BUFFER ARRAYS
c       IMOD = 'COLD' : DERIVES ONLY COLUMN DENSITY AND IMPACT PAR.
c       IMOD = 'REST' : ADS ALL THE OTHER QUANTITIES
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sumdata (t, de, dh, fi, dvoluni, dr, rdis, imod)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 rno,wei,wei2,dfb
      real*8 t, de, dh, fi, dvoluni, dr, rdis
c      real*8 fl,ft
      integer*4 series, line,i,j,atom
      integer*4 luty
c
      character imod*4
c
      rno=1.d17
      luty=6
c
      if ((imod.ne.'ALL').and.(imod.ne.'COLD').and.(imod.ne.'REST'))
     &then
        write (luty,10) imod
   10    format(/' MODE IMPROPERLY SET FOR SUBR. SUMDATA  : ',a4)
        stop
      endif
c
c     ***INTEGRATION FOR THE AVERAGE IONISATION : PAM(6,11) ,
c     WEIGHT IN : DEAV  ;  VOLUME STEP : DVOLUNI
c     AVERAGE ELECTRONIC DENSITIES IN DEAM(6,11)
c     AVERAGE DISTANCES IN RDISA(6,11)
c
      if (imod.ne.'COLD') then
        wei=((de*dh)*fi)*dvoluni
        wei2=wei*de
        dhav=dhav+wei2
        deav=deav+wei
        do 30 j=1,atypes
          do 20 i=1,maxion(j)
            deam(i,j)=deam(i,j)+(wei2*pop(i,j))
            rdisa(i,j)=rdisa(i,j)+((wei*pop(i,j))*(rdis/rno))
            pam(i,j)=pam(i,j)+(wei*pop(i,j))
   20     continue
   30   continue
c
c     ***INTEGRATION FOR THE AVERAGE TEMPERATURES : TEAV(6,11) ,
c     WEIGHT IN : PAM(6,11)  ;  VOLUME STEP : DVOLUNI
c
        do 50 j=1,atypes
          do 40 i=1,maxion(j)
            teav(i,j)=teav(i,j)+((t*wei)*pop(i,j))
   40     continue
   50   continue
c
c     ***INTEGRATION FOR THE AVERAGE TEMPERATURES GIVEN BY LINE
c     RATIOS : TOIII,TNII  (INTEG. EL. DENS. IN : DEOIII,DENII)
c     WEIGHT IN : WEOIII,WENII  ;  VOLUME STEP : DVOLUNI
c
        if (ox3.ne.0) then
          wei=(fmbri(6,ox3)+fmbri(8,ox3))
          weoiii=weoiii+((wei*dvoluni)*fi)
          roiii=roiii+(((wei*dvoluni)*fi)*(wei/(epsilon+fmbri(10,ox3))))
          deoiii=deoiii+(((de*wei)*dvoluni)*fi)
        endif
        if (ni2.ne.0) then
          wei=(fmbri(7,ni2)+fmbri(10,ni2))
          wenii=wenii+((wei*dvoluni)*fi)
          rnii=rnii+(((wei*dvoluni)*fi)*(wei/(epsilon+fmbri(13,ni2))))
          denii=denii+(((de*wei)*dvoluni)*fi)
        endif
c
c     ***INTEGRATION FOR THE AVERAGE DENSITIES GIVEN BY LINE
c     RATIOS : TOII,TSII  (INTEG. TEMP. IN : TOII,TSII)
c     WEIGHT IN : WEOII,WESII  ;  VOLUME STEP : DVOLUNI
c
        if (ox2.ne.0) then
          wei=fmbri(1,ox2)
          weoii=weoii+((wei*dvoluni)*fi)
          roii=roii+(((wei*dvoluni)*fi)*(wei/(epsilon+fmbri(2,ox2))))
          toii=toii+(((t*wei)*dvoluni)*fi)
        endif
        if (su2.ne.0) then
          wei=fmbri(2,su2)
          wesii=wesii+((wei*dvoluni)*fi)
          rsii=rsii+(((wei*dvoluni)*fi)*(wei/(epsilon+fmbri(1,su2))))
          tsii=tsii+(((t*wei)*dvoluni)*fi)
        endif
c
c     ***INTEGRATION OF LINE INTENSITIES IN /SLINE/
c     WEIGHT IN : FHBETA  ;  VOLUME STEP : DVOLUNI
c
        wei=dvoluni*fi
c
        do j=1,nxr3lines
          xr3lines_flux(j)=xr3lines_flux(j)+(xr3lines_bri(j)*wei)
        enddo
c
        do j=1,nxrllines
          xrllines_flux(j)=xrllines_flux(j)+(xrllines_bri(j)*wei)
        enddo
C
C       do j=1,xlines
C         fluxx(j)=fluxx(j)+(xrbri(j)*wei)
C       enddo
C
C       do j=1,xilines
C         fluxxi(j)=fluxxi(j)+(xibri(j)*wei)
C       enddo
C
C       do j=1,xhelines
C         fheif(j)=fheif(j)+(xhebri(j)*wei)
C       enddo
c
        do 60 j=1,nlines
          fluxr(j)=fluxr(j)+(rbri(j)*wei)
   60   continue
c
        do 70 j=1,mlines
          fluxi(j)=fluxi(j)+(fsbri(j)*wei)
   70   continue
c
        do j=1,nfmions
          do i=1,nfmtrans(j)
            fluxm(i,j)=fluxm(i,j)+(fmbri(i,j)*wei)
          enddo
        enddo
c
        do j=1,nfeions
          do i=1,nfetrans(j)
            fluxfe(i,j)=fluxfe(i,j)+(febri(i,j)*wei)
          enddo
        enddo
c
        do j=1,nf3ions
          do i=1,nf3trans
            fluxf3(i,j)=fluxf3(i,j)+(f3bri(i,j)*wei)
          enddo
        enddo
c
c
c CII recomb - all in one loop
c
        do j=1,nrccii
          fluxrccii_a(j)=fluxrccii_a(j)+(rccii_abri(j)*wei)
          fluxrccii_b(j)=fluxrccii_b(j)+(rccii_bbri(j)*wei)
        enddo
c
c NII recomb - all in one loop
c
        do j=1,nrcnii
          fluxrcnii_a(j)=fluxrcnii_a(j)+(rcnii_abri(j)*wei)
          fluxrcnii_b(j)=fluxrcnii_b(j)+(rcnii_bbri(j)*wei)
        enddo
c
c  OI rec: Separate Quintets and Triplets
c
        do j=1,nrcoi_q
          fluxrcoi_qa(j)=fluxrcoi_qa(j)+(rcoi_qabri(j)*wei)
          fluxrcoi_qb(j)=fluxrcoi_qb(j)+(rcoi_qbbri(j)*wei)
        enddo
        do j=1,nrcoi_t
          fluxrcoi_ta(j)=fluxrcoi_ta(j)+(rcoi_tabri(j)*wei)
          fluxrcoi_tb(j)=fluxrcoi_tb(j)+(rcoi_tbbri(j)*wei)
        enddo
c
c OII recomb - all in one loop
c
        do j=1,nrcoii
          fluxrcoii_a(j)=fluxrcoii_a(j)+(rcoii_abri(j)*wei)
          fluxrcoii_b(j)=fluxrcoii_b(j)+(rcoii_bbri(j)*wei)
          fluxrcoii_c(j)=fluxrcoii_c(j)+(rcoii_cbri(j)*wei)
        enddo
c
c NeII recomb - all in one loop
c
        do j=1,nrcneii
          fluxrcneii_a(j)=fluxrcneii_a(j)+(rcneii_abri(j)*wei)
          fluxrcneii_b(j)=fluxrcneii_b(j)+(rcneii_bbri(j)*wei)
        enddo
c
        do j=1,10
          fluxh(j)=fluxh(j)+(hbri(j)*wei)
        enddo
c
        do j=1,nheilines
          fluxhei(j)=fluxhei(j)+(heibri(j)*wei)
        enddo
c
        do j=1,nheislines
          fluxheis(j)=fluxheis(j)+(heisbri(j)*wei)
        enddo
c
        do j=1,nheitlines
          fluxheit(j)=fluxheit(j)+(heitbri(j)*wei)
        enddo
c
        heiioiiibfsum=heiioiiibfsum+(heiioiiibf*helibri(1,1)*wei)
c
        do series=1,nhseries
          do line=1,nhlines
            hydroflux(line,series)=hydroflux(line,series)+
     &       (hydrobri(line,series)*wei)
          enddo
        enddo
c
        do series=1,nheseries
          do line=1,nhelines
            if ((series.eq.1).and.(line.eq.1)) then
              dfb=1.d0-heiioiiibf
              heliflux(line,series)=heliflux(line,series)+dfb*
     &         (helibri(line,series)*wei)
            else
              heliflux(line,series)=heliflux(line,series)+(helibri(line,
     &         series)*wei)
            endif
          enddo
        enddo
c
        do i=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              xhydroflux(line,series,i)=xhydroflux(line,series,i)+
     &         (xhydrobri(line,series,i)*wei)
            enddo
          enddo
        enddo
c
c
        fheiilr=fheiilr+(heiilr*wei)
        h2qav=h2qav+(h2ql*wei)
        hei2qa=hei2qa+(hein(2)*wei)
        heii2qa=heii2qa+(heii2ql*wei)
        do atom=1,atypes
          h2qflux(atom)=h2qflux(atom)+(h2qbri(atom)*wei)
        enddo
        fhbeta=fhbeta+(hbeta*wei)
        tlosac=tlosac+(eloss*wei)
        tforbi=tforbi+(fmloss*wei)
c
      endif
c
c     ***INTEGRATION OF THE COLUMN DENSITY IN POPINT(6,11)
c     INTEGRATION OF PHOTON IMPACT PARAMETERS
c     INTEGRATION OF Hydrogen COLUMN (for dust)
c     INTEGRATION OF PAH COLUMN
c     RADIUS STEP SIZE : DR
c
      if (imod.ne.'REST') then
        wei2=((de*dh)*fi)*dvoluni
        dhavv=dhavv+wei2
        qhdha=qhdha+(wei2*qhdh)
        zetaeav=zetaeav+(wei2*zetae)
        wei=(dh*fi)*dr
        do j=1,atypes
          do i=1,maxion(j)
            popint(i,j)=popint(i,j)+((wei*pop(i,j))*zion(j))
          enddo
c            popint(maxion(j)+1,j) = popint(maxion(j)+1,j)+wei*Zion(j)
        enddo
        if (grainmode.eq.1) then
          dustint=dustint+wei
          if ((pahmode.eq.1).and.(pahactive.eq.1)) then
            do i=1,pahi
              pahint(i)=pahint(i)+wei*pahz(i)
            enddo
          endif
        endif
      endif
c
      return
      end
