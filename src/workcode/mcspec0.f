      subroutine mcspec(ix,iy,iz,luop)

      include 'cblocks.inc'
c
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
      integer*4 i,j,line,series,atom
      integer*4 ionidx, nt, omtype
      integer*4 luop,ix,iy,iz
c
      real*8 fnair,wei,dfb
c
c
   10 format(' ',i3,' ',i3,' ',i3,
     &       ' ',f12.3,', ',1pe11.4,', ',1pe11.4,
     &       ', ',a3,a6,', ',a4,', ',i4)
c
c
c      chklim=epsilon
      chklim=1.0d-30
      linecount=0
c
c
        wei=1.d0
c
        do j=1,nxr3lines
          xr3lines_flux(j)=xr3lines_bri(j)*wei
        enddo
c
        do j=1,nxrllines
          xrllines_flux(j)=xrllines_bri(j)*wei
        enddo
C
        do 60 j=1,nlines
          fluxr(j)=rbri(j)*wei
   60   continue
c
        do 70 j=1,mlines
          fluxi(j)=fsbri(j)*wei
   70   continue
c
        do j=1,nfmions
          do i=1,nfmtrans(j)
            fluxm(i,j)=fmbri(i,j)*wei
          enddo
        enddo
c
        do j=1,nfeions
          do i=1,nfetrans(j)
            fluxfe(i,j)=febri(i,j)*wei
          enddo
        enddo
c
        do j=1,nf3ions
          do i=1,nf3trans
            fluxf3(i,j)=f3bri(i,j)*wei
          enddo
        enddo
c
c
c CII recomb - all in one loop
c
        do j=1,nrccii
          fluxrccii_a(j)=rccii_abri(j)*wei
          fluxrccii_b(j)=rccii_bbri(j)*wei
        enddo
c
c NII recomb - all in one loop
c
        do j=1,nrcnii
          fluxrcnii_a(j)=rcnii_abri(j)*wei
          fluxrcnii_b(j)=rcnii_bbri(j)*wei
        enddo
c
c  OI rec: Separate Quintets and Triplets
c
        do j=1,nrcoi_q
          fluxrcoi_qa(j)=rcoi_qabri(j)*wei
          fluxrcoi_qb(j)=rcoi_qbbri(j)*wei
        enddo
        do j=1,nrcoi_t
          fluxrcoi_ta(j)=rcoi_tabri(j)*wei
          fluxrcoi_tb(j)=rcoi_tbbri(j)*wei
        enddo
c
c OII recomb - all in one loop
c
        do j=1,nrcoii
          fluxrcoii_a(j)=rcoii_abri(j)*wei
          fluxrcoii_b(j)=rcoii_bbri(j)*wei
          fluxrcoii_c(j)=rcoii_cbri(j)*wei
        enddo
c
c NeII recomb - all in one loop
c
        do j=1,nrcneii
          fluxrcneii_a(j)=rcneii_abri(j)*wei
          fluxrcneii_b(j)=rcneii_bbri(j)*wei
        enddo
c
        do j=1,10
          fluxh(j)=hbri(j)*wei
        enddo
c
        do j=1,nheilines
          fluxhei(j)=heibri(j)*wei
        enddo
c
        do j=1,nheislines
          fluxheis(j)=heisbri(j)*wei
        enddo
c
        do j=1,nheitlines
          fluxheit(j)=heitbri(j)*wei
        enddo
c
        heiioiiibfsum=heiioiiibf*helibri(1,1)*wei
c
        do series=1,nhseries
          do line=1,nhlines
            hydroflux(line,series)=hydrobri(line,series)*wei
          enddo
        enddo
c
        do series=1,nheseries
          do line=1,nhelines
            if ((series.eq.1).and.(line.eq.1)) then
              dfb=1.d0-heiioiiibf
              heliflux(line,series)=helibri(line,series)*wei
            else
              heliflux(line,series)=helibri(line,series)*wei
            endif
          enddo
        enddo
c
        do i=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              xhydroflux(line,series,i)=xhydrobri(line,series,i)*wei
            enddo
          enddo
        enddo
c
c
        fheiilr=heiilr*wei
        h2qav=h2ql*wei
        hei2qa=hein(2)*wei
        heii2qa=heii2ql*wei
        do atom=1,atypes
          h2qflux(atom)=h2qbri(atom)*wei
        enddo
        fhbeta=hbeta*wei
        tlosac=eloss*wei
        tforbi=fmloss*wei
c
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
      call heapindexsort (linecount, linelam, lineidx)

      do idx=1,linecount
c
        if (linelam(lineidx(idx)).le.8000) then
           write (luop,10) ix,iy,iz,
     &           linelam(lineidx(idx)),lmev/(linelam(lineidx(idx)
     &           )),linespec(lineidx(idx)),elem(lineid(lineidx(idx),1)),
     &           rom(lineid(lineidx(idx),2)),linekind(lineidx(idx)),
     &           lineacc(lineidx(idx))
        endif
c
      enddo
c
      return
      end
