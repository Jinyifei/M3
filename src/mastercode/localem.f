cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES LOCAL EMISSIVITY IN VECTORS EMIDIF,EMILIN,EMIDIFCONT
c     TO BE USED IN SUBROUTINE TOTPHOT AND NEWDIF
c
c     NOTE : UNITS FOR THE VECTORS EMIDIF,EMILIN ARE IN NUMBERS
c            OF PHOTONS (INSTEAD OF ERGS LIKE ALL OTHER VECTORS
c          (Units are photons/s/cm^3/Hz/sr (1/4pi))
c
c     CALL SUBR. RESON,INTER,HYDRO,ALLRATES
c     NOTE : CALLING THESE SUBR. COULD CHANGE THE INTENSITIES
c            OF THE LINES IN THE BUFFER VECTORS,SO IT IS NECESSARY
c            TO CALL SUBR. COOL BEFORE SUMMING LINE INTENS. (SUMDATA)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine localem (t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t, de, dh, telc
      real*8 xp(mxion, mxelem)
      real*8 energ,rkt
      real*8 wid
c
c      integer*4 ie,ii2,jn
      integer i,j,inl,nz
      integer line,series,atom,ion,trans
c
      telc=dmax1(t,1.d-2)
      rkt=rkb*telc
c
      do j=1,atypes
        do i=1,maxion(j)-1
          xp(i,j)=1.0d0
        enddo
      enddo
c
c     ***FIND LINES BRIGHTNESSES
c
      call cool (t, de, dh)
c
      jcon='YES'
c
      do inl=1,infph
        emidifcont(inl)=0.d0
      enddo
c
      call freebound (t, de, dh)
      call freefree (t, de, dh)
      call twophoton (t, de, dh)
c
      do inl=1,infph
c
        emidifcont(inl)=ffph(inl)+fbph(inl)+p2ph(inl)
        emidif(inl)=emidifcont(inl)
c
      enddo
cc
c      endif
cc
c  If dust and IR included add dust continuum to continuum emission
c
      if (grainmode.eq.1) then
        if (irmode.ne.0) then
          do inl=1,infph
            emidifcont(inl)=emidifcont(inl)+irphot(inl)
            emidif(inl)=emidifcont(inl)
          enddo
        endif
        if ((pahmode.eq.1).and.(pahactive.eq.1).and.(irmode.ne.0)) then
          do inl=1,infph
            emidifcont(inl)=emidifcont(inl)
     &                     +paheng*pahflux(inl)*pahfrac*dh
            emidif(inl)=emidifcont(inl)
          enddo
        endif
      endif
c
c     Now add weak lines directly to vector emidif...but not
c     emidifcont
c
c    Multi-level atoms
c
      do ion=1,nfmions
        do trans=1,nfmtrans(ion)
          j=fmbin(trans,ion)
          if ((j.ne.0).and.(fmbri(trans,ion).gt.epsilon)) then
            energ=fmeij(trans,ion)
            wid=evplk*(photev(j+1)-photev(j))
            emidif(j)=emidif(j)+(fmbri(trans,ion)/(wid*energ))
          endif
        enddo
      enddo
cc
cc    Fe II-VII lines
c
      do ion=1,nfeions
        do trans=1,nfetrans(ion)
          j=febin(trans,ion)
          if ((j.ne.0).and.(febri(trans,ion).gt.epsilon)) then
            energ=feeij(trans,ion)
            wid=evplk*(photev(j+1)-photev(j))
            emidif(j)=emidif(j)+(febri(trans,ion)/(wid*energ))
          endif
        enddo
      enddo
cc
c
c    Three level atoms
c
      do ion=1,nf3ions
        do trans=1,nf3trans
          j=f3bin(trans,ion)
          if ((j.ne.0).and.(f3bri(trans,ion).gt.epsilon)) then
            energ=(plk*cls)/f3lam(trans,ion)
            wid=evplk*(photev(j+1)-photev(j))
            emidif(j)=emidif(j)+(f3bri(trans,ion)/(wid*energ))
          endif
        enddo
      enddo
c
c     Add old intercombination lines to vector EMDIF.
c
      do i=1,mlines
        j=lcbin(i)
        if ((j.ne.0).and.(fsbri(i).gt.epsilon)) then
          energ=e12fs(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+(fsbri(i)/(wid*energ))
        endif
      enddo
c
C     do i=1,xilines
C       j=xibin(i)
C       if ((j.gt.0).and.(xibri(i).gt.epsilon)) then
C         energ=ev*xiejk(i)
C         wid=evplk*(photev(j+1)-photev(j))
C         emidif(j)=emidif(j)+(xibri(i)/(wid*energ))
C       endif
C     enddo
c
C     do i=1,xhelines
C       j=xhebin(i)
C       if ((j.gt.0).and.(xhebri(i).gt.epsilon)) then
C         energ=ev*xiejk(i)
C         wid=evplk*(photev(j+1)-photev(j))
C         emidif(j)=emidif(j)+(xhebri(i)/(wid*energ))
C       endif
C     enddo
c
c old intercombination lines
c
      do i=1,nlines
        j=lrbin(i)
        if ((j.ne.0).and.(rbri(i).gt.epsilon)) then
          energ=e12r(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rbri(i)/(wid*energ)
        endif
      enddo
c Original He I lines, just the first three
      do i=1,3
        j=heibin(i)
        if ((j.ne.0).and.(heibri(i).gt.epsilon)) then
          energ=ev*(lmev/(heilam(i)*1.d8))
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+heibri(i)/(wid*energ)
        endif
      enddo
c
c New He I singlet lines
c
      do i=1,nheislines
        j=heisbin(i)
        if ((j.ne.0).and.(heisbri(i).gt.epsilon)) then
          energ=heiseij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+heisbri(i)/(wid*energ)
        endif
      enddo
c
c New He I triplet lines
c
      do i=1,nheitlines
        j=heitbin(i)
        if ((j.ne.0).and.(heitbri(i).gt.epsilon)) then
          energ=heiteij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+heitbri(i)/(wid*energ)
        endif
      enddo
c
c heavy element recomb lines
c
      do i=1,nrccii
        j=rccii_bin(i)
        if ((j.ne.0).and.(rccii_bbri(i).gt.epsilon)) then
          energ=rccii_eij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rccii_bbri(i)/(wid*energ)
        endif
      enddo
c
      do i=1,nrcnii
        j=rcnii_bin(i)
        if ((j.ne.0).and.(rcnii_bbri(i).gt.epsilon)) then
          energ=rcnii_eij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rcnii_bbri(i)/(wid*energ)
        endif
      enddo
c
      do i=1,nrcoi_q
        j=rcoi_qbin(i)
        if ((j.ne.0).and.(rcoi_qbbri(i).gt.epsilon)) then
          energ=rcoi_qeij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rcoi_qbbri(i)/(wid*energ)
        endif
      enddo
c
      do i=1,nrcoi_t
        j=rcoi_tbin(i)
        if ((j.ne.0).and.(rcoi_tbbri(i).gt.epsilon)) then
          energ=rcoi_teij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rcoi_tbbri(i)/(wid*energ)
        endif
      enddo
c
      do i=1,nrcoii
        j=rcoii_bin(i)
        if ((j.ne.0).and.(rcoii_bbri(i).gt.epsilon)) then
          energ=rcoii_eij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rcoii_bbri(i)/(wid*energ)
        endif
      enddo
c
      do i=1,nrcneii
        j=rcneii_bin(i)
        if ((j.ne.0).and.(rcneii_bbri(i).gt.epsilon)) then
          energ=rcneii_eij(i)
          wid=evplk*(photev(j+1)-photev(j))
          emidif(j)=emidif(j)+rcneii_bbri(i)/(wid*energ)
        endif
      enddo
c
c     adds xr3lines to xr3lines_emilin(1,line)
c
      do line=1,nxr3lines
        xr3lines_emilin(1,line)=0.d0
        if (xr3lines_bri(line).gt.epsilon) then
          j=xr3lines_bin(line)
          if (j.ne.0) then
          energ=xr3lines_egij(line)
          wid=evplk*(photev(j+1)-photev(j))
          xr3lines_emilin(1,line)=xr3lines_bri(line)/(wid*energ)
        endif
        endif
      enddo
c
c     adds xrllines to xrllines_emilin(1,line)
c
      do line=1,nxrllines
        xrllines_emilin(1,line)=0.d0
        if (xrllines_bri(line).gt.epsilon) then
          j=xrllines_bin(line)
          if (j.ne.0) then
          energ=xrllines_egij(line)
          wid=evplk*(photev(j+1)-photev(j))
          xrllines_emilin(1,line)=xrllines_bri(line)/(wid*energ)
        endif
        endif
      enddo
c
c     adds xlines to emilin(1,line)
c
C     do line=1,xlines
C       emilin(1,line)=0.d0
C       if (xrbri(line).gt.epsilon) then
C         j=xbin(line)
C         if (j.ne.0) then
C         energ=ev*xejk(line)
C         wid=evplk*(photev(j+1)-photev(j))
C         emilin(1,line)=xrbri(line)/(wid*energ)
C       endif
C       endif
C     enddo
c
c
c     adds hydrogen and helium lines to hydlin(1,line,series) and
c      hellin(1,line,series)
c
      do series=1,nhseries
        do line=1,nhlines
          hydlin(1,line,series)=0.d0
          if (hydrobri(line,series).gt.epsilon) then
          j=hbin(line,series)
          if (j.ne.0) then
            energ=ev*lmev/hlambda(line,series)
            wid=evplk*(photev(j+1)-photev(j))
            hydlin(1,line,series)=hydrobri(line,series)/(wid*energ)
          endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          hellin(1,line,series)=0.d0
          if (helibri(line,series).gt.epsilon) then
          j=hebin(line,series)
          if (j.ne.0) then
            energ=ev*lmev/helambda(line,series)
            wid=evplk*(photev(j+1)-photev(j))
            hellin(1,line,series)=helibri(line,series)/(wid*energ)
          endif
          endif
        enddo
      enddo
c
c     Adds heavy hydrogenic series
c
      if (atypes.gt.3) then
        do atom=3,atypes
c
          nz=mapz(atom)
c
          do series=1,nxhseries
            do line=1,nxhlines
              xhydlin(1,line,series,atom)=0.d0
              if (xhydrobri(line,series,atom).gt.epsilon) then
              j=xhbin(line,series,atom)
              if (j.ne.0) then
c            energ = 0.5d0*(ev*(photev(j+1)+photev(j)))
                energ=ev*lmev/xhlambda(line,series,atom)
                wid=evplk*(photev(j+1)-photev(j))
                xhydlin(1,line,series,atom)=xhydrobri(line,series,atom)/
     &           (wid*energ)
              endif
              endif
            enddo
          enddo
c      endif
        enddo
      endif
c
      return
c
      end
