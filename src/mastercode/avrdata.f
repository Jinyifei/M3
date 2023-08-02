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
c*******TO OBTAIN THE AVERAGE TEMP. , IONIC POP.,DENSITIES AND
c       LINE INTENSITIES USING SUM DATA IN  : /SLINE/,/TONLIN/
c       AND /DELIN/
c       CALL SUBR. FINDTDE
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine avrdata ()
c
      include 'cblocks.inc'
c
c
      real*8 rno
      integer*4 i, j, line,series
c
c    ***AVERAGES IONIC POP,TEMP.,DISTANCES AND DENSITIES
c
      rno=1.d17
      do 20 j=1,atypes
        do 10 i=1,maxion(j)
          teav(i,j)=teav(i,j)/(epsilon+pam(i,j))
          deam(i,j)=deam(i,j)/(epsilon+pam(i,j))
          rdisa(i,j)=(rdisa(i,j)/(epsilon+pam(i,j)))*rno
          pam(i,j)=pam(i,j)/(epsilon+deav)
   10   continue
   20 continue
c
      deav=dhav/deav
      qhdha=qhdha/(epsilon+dhavv)
c
      zetaeav=zetaeav/(epsilon+dhavv)
      roiii=roiii/(epsilon+weoiii)
      deoiii=deoiii/(epsilon+weoiii)
      rnii=rnii/(epsilon+wenii)
      denii=denii/(epsilon+wenii)
      roii=roii/(epsilon+weoii)
      toii=toii/(epsilon+weoii)
      rsii=rsii/(epsilon+wesii)
      tsii=tsii/(epsilon+wesii)
c
c    ***FIND TNII,TOIII,DESII,DEOII
c
      call findtde
c
c    ***FORM RATIO OF EMISSION LINES RELATIVE TO  H-BETA
c
      do i=1,nfmions
        do j=1,nfmtrans(i)
          fluxm(j,i)=fluxm(j,i)/(epsilon+fhbeta)
        enddo
      enddo
c
      do i=1,nfeions
        do j=1,nfetrans(i)
          fluxfe(j,i)=fluxfe(j,i)/(epsilon+fhbeta)
        enddo
      enddo
c
      do i=1,nf3ions
        do j=1,nf3trans
          fluxf3(j,i)=fluxf3(j,i)/(epsilon+fhbeta)
        enddo
      enddo
c
      do j=1,10
        fluxh(j)=fluxh(j)/(epsilon+fhbeta)
      enddo
      do j=1,nheilines
        fluxhei(j)=fluxhei(j)/(epsilon+fhbeta)
      enddo
      do j=1,nheislines
        fluxheis(j)=fluxheis(j)/(epsilon+fhbeta)
      enddo
      do j=1,nheitlines
        fluxheit(j)=fluxheit(j)/(epsilon+fhbeta)
      enddo
c CII
      do j=1,nrccii
        fluxrccii_a(j)=fluxrccii_a(j)/(epsilon+fhbeta)
        fluxrccii_b(j)=fluxrccii_b(j)/(epsilon+fhbeta)
      enddo
c NII
      do j=1,nrcnii
        fluxrcnii_a(j)=fluxrcnii_a(j)/(epsilon+fhbeta)
        fluxrcnii_b(j)=fluxrcnii_b(j)/(epsilon+fhbeta)
      enddo
c OI
      do j=1,nrcoi_q
        fluxrcoi_qa(j)=fluxrcoi_qa(j)/(epsilon+fhbeta)
        fluxrcoi_qb(j)=fluxrcoi_qb(j)/(epsilon+fhbeta)
      enddo
      do j=1,nrcoi_t
        fluxrcoi_ta(j)=fluxrcoi_ta(j)/(epsilon+fhbeta)
        fluxrcoi_tb(j)=fluxrcoi_tb(j)/(epsilon+fhbeta)
      enddo
c OII
      do j=1,nrcoii
        fluxrcoii_a(j)=fluxrcoii_a(j)/(epsilon+fhbeta)
        fluxrcoii_b(j)=fluxrcoii_b(j)/(epsilon+fhbeta)
        fluxrcoii_c(j)=fluxrcoii_c(j)/(epsilon+fhbeta)
      enddo
c
c NeII
c
      do j=1,nrcneii
        fluxrcneii_a(j)=fluxrcneii_a(j)/(epsilon+fhbeta)
        fluxrcneii_b(j)=fluxrcneii_b(j)/(epsilon+fhbeta)
      enddo
c
      heiioiiibfsum=heiioiiibfsum/(epsilon+fhbeta)
c
      do series=1,nhseries
        do line=1,nhlines
          hydroflux(line,series)=hydroflux(line,series)/(epsilon+fhbeta)
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          heliflux(line,series)=heliflux(line,series)/(epsilon+fhbeta)
        enddo
      enddo
c
      do i=1,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            xhydroflux(line,series,i)=xhydroflux(line,series,i)/
     &       (epsilon+fhbeta)
          enddo
        enddo
      enddo
c
      do j=1,nlines
        fluxr(j)=fluxr(j)/(epsilon+fhbeta)
      enddo
c
      do j=1,mlines
        fluxi(j)=fluxi(j)/(epsilon+fhbeta)
      enddo
c
      do j=1,nxr3lines
        xr3lines_flux(j)=xr3lines_flux(j)/(epsilon+fhbeta)
      enddo
c
      do j=1,nxrllines
        xrllines_flux(j)=xrllines_flux(j)/(epsilon+fhbeta)
      enddo
c
C     do j=1,xlines
C       fluxx(j)=fluxx(j)/(epsilon+fhbeta)
C     enddo
c
C     do j=1,xilines
C       fluxxi(j)=fluxxi(j)/(epsilon+fhbeta)
C     enddo
c
C     do j=1,xhelines
C       fheif(j)=fheif(j)/(epsilon+fhbeta)
C     enddo
c
      fheiilr=fheiilr/(epsilon+fhbeta)
      h2qav=h2qav/(epsilon+fhbeta)
      hei2qa=hei2qa/(epsilon+fhbeta)
      heii2qa=heii2qa/(epsilon+fhbeta)
c
      do i=1,atypes
        h2qflux(i)=h2qflux(i)/(epsilon+fhbeta)
      enddo
c
      return
      end
