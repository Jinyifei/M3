cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine to initialise mappings for non interactive data
c     ie data files and so on....
c
c     if the file map.prefs is found in the current working directory
c     then it will be used in place of data/ATDAT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine mapinit (error)
c
c
      include 'cblocks.inc'
c
c
c           Variables
c
      real*8 en,x,ez,xe,efit,lz
      integer*4 i,ion,j,k,kpre,line,luin,atom,ionindex
      integer*4 luop,min,npre,series,trans,m,nz
      integer*4 checksum
c
      character ibell*4,ilgg*4, ibuf(19)*4
cc
cc Test io variables
c      character sfx*4, fl*20
c      real*8 rate,e,px,t
c      real*8 fkramer,vernerphoto
c      integer*4 ir,l
cc
c
      logical error,iexi
c
c externals
c
      integer*4 lenv
      real*8 fgfflog
c
      error=.false.
c
c Default kappa > 1000 -> boltzmann excitation
c
      kappa=1.0d99
      kappaidx=9
      usekappa=.false.
      usekappainterp=.false.
      kappaa=1.d0
      kappab=0.d0
c
      luin=20
      luop=23
c
      ibell(1:4)=char(7)
c
      inquire (file='data/ATDAT.txt',exist=iexi)
      if (iexi) then
        datadir='./'
      else
        datadir='/usr/local/share/mappings/'
      endif
      dtlen=lenv(datadir)
      inquire (file=datadir(1:dtlen)//'data/ATDAT.txt',exist=iexi)
      if (iexi.eqv..false.) then
        filename=datadir(1:dtlen)//'data/ATDAT.txt'
        m=lenv(filename)
        write (*,*) 'ERROR in mapinit: ',filename(1:m),' NOT FOUND.'
        write (*,*) '   MV requires a valid local data/ directory or'
        write (*,*) '   a valid shared /usr/local/share/mappings/data/'
        write (*,*) '   directory.  '
        stop
      endif
c
      filename=''
c
c    ***START READ-IN OF DATA
c
   10 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V : Begin Data Initialisation.',/
     & '     DATA DIR.: ',a,
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      if ((dtlen.eq.2).and.(datadir(1:dtlen).eq.'./')) then
        write (*,10) 'Local = data/'
      else
        write (*,10) 'Shared = '//datadir(1:dtlen)//'data/'
      endif
c
c     formats for reading files
c
   20 format(19a4)
   30 format(' ',19a4)
c
c     set up roman numerals array: 1-30
c
      rom(1)='I'
      rom_len(1)=1
      rom(2)='II'
      rom_len(2)=2
      rom(3)='III'
      rom_len(3)=3
      rom(4)='IV'
      rom_len(4)=2
      rom(5)='V'
      rom_len(5)=1
      rom(6)='VI'
      rom_len(6)=2
      rom(7)='VII'
      rom_len(7)=3
      rom(8)='VIII'
      rom_len(8)=4
      rom(9)='IX'
      rom_len(9)=2
      rom(10)='X'
      rom_len(10)=1
      rom(11)='XI'
      rom_len(11)=2
      rom(12)='XII'
      rom_len(12)=3
      rom(13)='XIII'
      rom_len(13)=4
      rom(14)='XIV'
      rom_len(14)=3
      rom(15)='XV'
      rom_len(15)=2
      rom(16)='XVI'
      rom_len(16)=3
      rom(17)='XVII'
      rom_len(17)=4
      rom(18)='XVIII'
      rom_len(18)=5
      rom(19)='XIX'
      rom_len(19)=3
      rom(20)='XX'
      rom_len(20)=2
      rom(21)='XXI'
      rom_len(21)=3
      rom(22)='XXII'
      rom_len(22)=4
      rom(23)='XXIII'
      rom_len(23)=5
      rom(24)='XXIV'
      rom_len(24)=4
      rom(25)='XXV'
      rom_len(25)=3
      rom(26)='XXVI'
      rom_len(26)=4
      rom(27)='XXVII'
      rom_len(27)=5
      rom(28)='XXVIII'
      rom_len(28)=6
      rom(29)='XXIX'
      rom_len(29)=4
      rom(30)='XXX'
      rom_len(30)=3
      rom(31)='XXXI'
      rom_len(31)=4
c
      do 40 i=1,30
        zmap(i)=0
        mapz(i)=0
   40 continue
c
      jiel='NO'
      jiem='NO'
c
      inquire (file='map.prefs',exist=iexi)
      if (iexi) then
        open (luin,file='map.prefs',status='OLD')
      else
        filename=datadir(1:dtlen)//'data/ATDAT.txt'
        open (luin,file=filename,status='OLD')
      endif
c
   50 read (luin,fmt=20) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 50
      write (*,30) (ibuf(j),j=1,19)
      zen=0.d0
      read (luin,fmt=*) atypes
      do 60 i=1,atypes
        read (luin,*) j,elem(i),k,atwei(i),zion0(i),dion0(i),maxion(i)
        zion(i)=zion0(i)
        deltazion(i)=1.d0
        zen=zen+zion0(i)
        dion(i)=dion0(i)
        invdion(i)=1.d0/dion(i)
        zmap(k)=j
        mapz(j)=k
        if (i.ne.j) goto 380
   60 continue
c
      do i=1,atypes
        j=mapz(i)
        elem_len(i)=2
        if (j.eq.1) elem_len(i)=1
        if (j.eq.5) elem_len(i)=1
        if (j.eq.6) elem_len(i)=1
        if (j.eq.7) elem_len(i)=1
        if (j.eq.8) elem_len(i)=1
        if (j.eq.9) elem_len(i)=1
        if (j.eq.15) elem_len(i)=1
        if (j.eq.16) elem_len(i)=1
        if (j.eq.19) elem_len(i)=1
        if (j.eq.23) elem_len(i)=1
      enddo
c
      close (luin)
c
c
c default switch settings
c
      collmode=0
      phionmode=0
      ph2mode=0
      recommode=0
      drrecmode=0
      chargemode=0
      chargeheatmode=1
      freefreemode=1
      x=fgfflog(-1,0.0d0,0.0d0,1,1)
      photonmode=1
      crate=0.d0
      hhecollmode=0
      linecoolmode=0
      alphacoolmode=0
      turbheatmode=0
      admach=0.0d0
      radpressmode=0
      expertmode=0
c
      filename=datadir(1:dtlen)//'data/switches.txt'
      inquire (file=filename,exist=iexi)
      if (iexi) then
        open (luin,file=filename,status='OLD')
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) collmode
        if (collmode.lt.0) collmode=0
        if (collmode.gt.2) collmode=0
        write (*,*) '*** collmode     :',collmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) phionmode
        if (phionmode.lt.0) phionmode=0
        if (phionmode.gt.1) phionmode=0
        write (*,*) '*** phionmode    :',phionmode
c
        ph2mode=0
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) ph2mode
        if (ph2mode.lt.0) ph2mode=0
        if (ph2mode.gt.1) ph2mode=0
        write (*,*) '*** ph2mode      :',ph2mode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) recommode
        if (recommode.lt.0) recommode=0
        if (recommode.gt.2) recommode=0
        write (*,*) '*** recommode    :',recommode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) drrecmode
        if (drrecmode.lt.0) recommode=0
        if (drrecmode.gt.2) drrecmode=0
        write (*,*) '*** drrecmode    :',drrecmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) chargemode
        if (chargemode.lt.0) chargemode=0
        if (chargemode.gt.3) chargemode=0
        write (*,*) '*** chargexfermode   :',chargemode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) chargeheatmode
        if (chargeheatmode.lt.0) chargeheatmode=1
        if (chargeheatmode.gt.1) chargeheatmode=1
        write (*,*) '*** chargexfer heating   :',chargeheatmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) freefreemode
        if (freefreemode.lt.0) freefreemode=1
        if (freefreemode.gt.1) freefreemode=1
        if (freefreemode.eq.1) x=fgfflog(-1,0.0d0,0.0d0,1,1)
        write (*,*) '*** freefreemode :',freefreemode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) photonmode
        if (photonmode.lt.0) photonmode=1
        if (photonmode.gt.1) photonmode=1
        crate=0.d0
        write (*,*) '*** photonmode   :',photonmode,crate
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) hhecollmode
        if (hhecollmode.lt.0) hhecollmode=0
        if (hhecollmode.gt.1) hhecollmode=0
        write (*,*) '*** hhecollmode  :',hhecollmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) linecoolmode
        if (linecoolmode.lt.0) linecoolmode=0
        if (linecoolmode.gt.1) linecoolmode=0
        write (*,*) '*** linecoolmode:',linecoolmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) alphacoolmode
        if (alphacoolmode.lt.0) alphacoolmode=0
        if (alphacoolmode.gt.1) alphacoolmode=0
        write (*,*) '*** alphacoolmode:',alphacoolmode
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) turbheatmode
        if (turbheatmode.lt.0) turbheatmode=0
        if (turbheatmode.gt.1) turbheatmode=0
        admach=0.0d0
        write (*,*) '*** turbheatmode :',turbheatmode,admach
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) radpressmode
        if (radpressmode.lt.0) radpressmode=0
        if (radpressmode.gt.1) radpressmode=0
        write (*,*) '*** radiation pressure mode :',radpressmode
c
c     expert mode triggers many more detailed questions,
c     useful for testing and program design.  NOT DOCUMENTED
c     SO THERE.  Contact Ralph Sutherland (Ralph.Sutherland@anu.edu.au)
c     for details if you are programming and modifying mappings.
c
        read (luin,20) (ibuf(j),j=1,19)
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) expertmode
        if (expertmode.lt.0) expertmode=0
        if (expertmode.gt.1) expertmode=0
        write (*,*) '*** expertmode   :',expertmode
c
      endif
c
      close (luin)
c
      do atom=1,mxelem
        do ion=1,mxion
          nialines(ion,atom)=0
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Read ionisation and recombination data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readiondata (luin, error)
      call readphiondata (luin, error)
      call readhhecoldata (luin, error)
      call readhhedata (luin, error)
      call readcolldata (luin, error)
      call readrecomdata (luin, error)
      call readdirecomdata (luin, error)
      call readcont (luin, error)
      call readchx (luin, error)
c energy vector
      call readphotdat (luin, error)
c older and misc atomic data
      call readion2 (luin, error)
      call readcoll2 (luin, error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LINES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Read line data
c
c      call readxdata (luin, error)
      call readxr3data (luin, error)
      call readxrldata (luin, error)
c
      call readheavyrec (luin, error)
      call read2level (luin, error)
      call read3level (luin, error)
c
c     readkappadat must be called before readmultilevel
c
      call readkappadat (luin, error)
      call readmultilevel (luin, error)
      call readmultife (luin, error)
c
      if (expertmode.gt.0) then
   70 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V : Line Counts by Ion and Element',
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
        write (*,70)
   80 format(/t9,30(2x,a2,1x))
        write (*,80) (elem(i),i=1,atypes)
   90 format(' ',a6,30(i5))
        do j=1,mxion
          checksum=0
          do i=1,atypes
            checksum=checksum+nialines(j,i)
          enddo
          if (checksum.gt.0) then
            write (*,90) rom(j),(nialines(j,i),i=1,atypes)
          endif
        enddo
  100 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/)
        write (*,100)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DUST
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readdust (luin, error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readstardat (luin, error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  110 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  MAPPINGS V: Data read successfully '
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      write (*,110)
  120 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     &  '  MAPPINGS V: Begin General Initialisation'
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      write (*,120)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Tidy Up and Precalc some tables, assign bins to lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      xe=ee2phe(zmap(2))/(0.75d0*ipote(1,zmap(2)))
c
c fill in arrays of two photon upper energy levels, interpolating
c for ions without primary atomic data. Potential scale for H
c poly fit correction for He.
c
      do atom=1,atypes
        nz=mapz(atom)
        ez=0.75d0*ipote(nz,atom)
        if (ee2p(atom).le.0.d0) ee2p(atom)=ez
        nz=mapz(atom)-1
        if (nz.gt.0) then
          ez=xe*(0.75d0*ipote(nz,atom))
          lz=dlog10(dble(nz))
          efit=(1.d0+lz*(1.5688838763d-01+lz*(-5.8259412168d-02+lz*(-
     &     8.1946404393d-03+lz*7.0718088682d-03))))
          if (ee2phe(atom).le.0.d0) ee2phe(atom)=ez/efit
        endif
c         write(*,*) mapz(atom),ee2p(atom),ee2phe(atom)
      enddo
c
c     Find bin to start phion loops
c
      do i=1,ionum
        photbinstart(i)=1
        en=ipotpho(i)
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            if (i.eq.1) ionstartbin=j
            photbinstart(i)=j
            goto 130
          endif
        enddo
  130   continue
      enddo
c
c      *FIND CORRESPONDING ENERGY BINS IN EPHOT(N) FOR
c      V3 resonance lines
c
      xr3minxbin=infph-1
      do i=1,nxr3lines
        en=xr3lines_egij(i)/ev
        xr3lines_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            xr3lines_bin(i)=j
            if (j.le.xr3minxbin) xr3minxbin=j
            goto 140
          endif
        enddo
  140   continue
      enddo
c
c
      xrlminxbin=infph-1
      do i=1,nxrllines
        en=xrllines_egij(i)/ev
        xrllines_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            xrllines_bin(i)=j
            if (j.le.xrlminxbin) xrlminxbin=j
            goto 150
          endif
        enddo
  150   continue
      enddo
c
c      minxbin=infph-1
c     do i=1,xlines
c       en=xejk(i)
c       xbin(i)=0
c       do j=1,infph-1
c         if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
c           xbin(i)=j
c           if (j.le.minxbin) minxbin=j
c           goto 110
c         endif
c       enddo
c 110   continue
c     enddo
c
c     Hydrogen
c
      do line=1,nhlines
        do series=1,nhseries
          en=lmev/hlambda(line,series)
          hbin(line,series)=0
          do j=1,infph-1
            if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
              hbin(line,series)=j
              goto 160
            endif
          enddo
  160     continue
        enddo
      enddo
c
c helium ii
c
      do line=1,nhelines
        do series=1,nheseries
          en=lmev/helambda(line,series)
          hebin(line,series)=0
          do j=1,infph-1
            if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
              hebin(line,series)=j
              goto 170
            endif
          enddo
  170     continue
        enddo
      enddo
c
c     Heavy Hydrogenic ions
c
      do atom=1,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            xhlambda(line,series,atom)=hlambda(line,series)/(mapz(atom)*
     &       mapz(atom))
            en=lmev/xhlambda(line,series,atom)
            xhbin(line,series,atom)=0
            do j=1,infph-1
              if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
                xhbin(line,series,atom)=j
                goto 180
              endif
            enddo
  180       continue
          enddo
        enddo
      enddo
c
c     Three Level Ions
c
      do ion=1,nf3ions
        do trans=1,nf3trans
c             en = 1.985940d-16/f3lam(trans,ion)
          en=(plk*cls/ev)/f3lam(trans,ion)
          f3bin(trans,ion)=0
          do j=1,infph-1
            if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
              f3bin(trans,ion)=j
              goto 190
            endif
          enddo
  190     continue
        enddo
      enddo
c
c     Multi-Level Ions
c
      do ionindex=1,nfmions
        do trans=1,nfmtrans(ionindex)
          en=(plk*cls/ev)/fmlam(trans,ionindex)
          fmbin(trans,ionindex)=0
          do j=1,infph-1
            if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
              fmbin(trans,ionindex)=j
              goto 200
            endif
          enddo
  200     continue
        enddo
      enddo
cc
c
c     Multi-Level Fe Ions
c
      do ionindex=1,nfeions
        do trans=1,nfetrans(ionindex)
          en=(plk*cls/ev)/felam(trans,ionindex)
          febin(trans,ionindex)=0
          do j=1,infph-1
            if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
              febin(trans,ionindex)=j
              goto 210
            endif
          enddo
  210     continue
        enddo
      enddo
cc
c     do i=1,xilines
c       en=xiejk(i)
c       xibin(i)=0
c       do j=1,infph-1
c         if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
c           xibin(i)=j
c           goto 180
c         endif
c       enddo
c 180   continue
c     enddo
c
c     do i=1,xhelines
c       en=xhejk(i)
c       xhebin(i)=0
c       do j=1,infph-1
c         if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
c           xhebin(i)=j
c           goto 185
c         endif
c       enddo
c 185   continue
c     enddo
c
      do i=1,nlines
        en=e12r(i)/ev
        lrbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            lrbin(i)=j
            goto 220
          endif
        enddo
  220   continue
      enddo
c
      do i=1,nheilines
        en=(lmev/(heilam(i)*1.d8))
        heibin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            heibin(i)=j
            goto 230
          endif
        enddo
  230   continue
      enddo
c New HeI lines
      do i=1,nheislines
        en=(lmev/heislam(i))
        heisbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            heisbin(i)=j
            goto 240
          endif
        enddo
  240   continue
      enddo
      do i=1,nheitlines
        en=(lmev/heitlam(i))
        heitbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            heitbin(i)=j
            goto 250
          endif
        enddo
  250   continue
      enddo
c heavy element recomb lines
      do i=1,nrccii
        en=(lmev/rccii_lam(i))
        rccii_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rccii_bin(i)=j
            goto 260
          endif
        enddo
  260   continue
      enddo
      do i=1,nrcnii
        en=(lmev/rcnii_lam(i))
        rcnii_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rcnii_bin(i)=j
            goto 270
          endif
        enddo
  270   continue
      enddo
      do i=1,nrcoi_q
        en=(lmev/rcoi_qlam(i))
        rcoi_qbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rcoi_qbin(i)=j
            goto 280
          endif
        enddo
  280   continue
      enddo
      do i=1,nrcoi_t
        en=(lmev/rcoi_tlam(i))
        rcoi_tbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rcoi_tbin(i)=j
            goto 290
          endif
        enddo
  290   continue
      enddo
      do i=1,nrcoii
        en=(lmev/rcoii_lam(i))
        rcoii_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rcoii_bin(i)=j
            goto 300
          endif
        enddo
  300   continue
      enddo
      do i=1,nrcneii
        en=(lmev/rcneii_lam(i))
        rcoii_bin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            rcneii_bin(i)=j
            goto 310
          endif
        enddo
  310   continue
      enddo
      do i=1,mlines
        en=e12fs(i)/ev
        lcbin(i)=0
        do j=1,infph-1
          if ((photev(j+1).gt.en).and.(photev(j).le.en)) then
            lcbin(i)=j
            goto 320
          endif
        enddo
  320   continue
      enddo
c
      npre=0
      kpre=-1
      do 350 i=1,nlines
        min=infph
        do 330 j=1,nlines
          if ((lrbin(j).eq.npre).and.(j.gt.kpre)) goto 340
          if ((lrbin(j).ge.min).or.(lrbin(j).le.npre)) goto 330
          min=lrbin(j)
          k=j
  330   continue
c
        linpos(1,i)=min
        linpos(2,i)=k
        kpre=k
        npre=min
        goto 350
  340   linpos(1,i)=lrbin(j)
        linpos(2,i)=j
        kpre=j
        npre=lrbin(j)
  350 continue
c
c     put branching ratios in xbr, plus xdismul init
c
c     do i=1,xlines
c       ni=xiso(i)
c       ntr=idint(xtrans(i))
c       if (dble(ntr).eq.xtrans(i)) then
c         xbr(i)=brr(ni,ntr)
c       else
c         xbr(i)=1.d0
c       endif
c     enddo
c
c      *FINDS IONISATIONS CROSS SECTION NUMBER FOR HEII
c
      do 360 i=1,ionum
        jhe2p=i
        if ((atpho(i).eq.2).and.(ionpho(i).eq.2)) goto 370
  360 continue
  370 continue
c
      goto 390
c
  380 error=.true.
c
  390 continue
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readiondata (luin, error)
c
      include 'cblocks.inc'
c
      integer*4 luin
      logical error
c
c           Variables
c
      real*8 ep
      integer*4 at,atom,i,ion,io,j,nions
      character ilgg*4, ibuf(19)*4
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
c
      filename=datadir(1:dtlen)//'data/ionisation/IONDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
c     read primary ionisation potentials
c     and convert to ergs
c
      read (luin,fmt=*) nions
      do i=1,nions
        read (luin,fmt=*) at,io,ep
        if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            ion=io
            ipotev(ion,atom)=ep
            ipote(ion,atom)=ep*ev
          endif
        endif
      enddo
c
      close (luin)
c
      do atom=1,atypes
        ee2p(atom)=0.0d0
        ee2phe(atom)=0.0d0
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readphiondata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 s,si
      real*8 p1(5), p2(7)
      real*8 be,den,ep,g,rt,au,pl,edge
c
      integer*4 z,atom,io,ion,luin,nentries,ni
      integer*4 ph2,shell,ne,fb,wr,wi
      integer*4 i,j,k
c shell angular momentum quantum numbers
      integer*4 lnumber(7)
c
      character ilgg*4, ibuf(19)*4
      logical error,found(1696)
c
      real*8 vernerphoto
c
      data lnumber/0,0,1,0,1,2,0/
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      if (phionmode.eq.1) then
c
        filename=datadir(1:dtlen)//'data/ionisation/PHIONDAT.txt'
        open (luin,file=filename,status='OLD')
   30   read (luin,fmt=10) (ibuf(j),j=1,19)
        ilgg=ibuf(1)
        if (ilgg(1:1).eq.'%') goto 30
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) nentries
        ionum=0
        den=0.0d0
        do 50 i=1,nentries
          atom=0
          ion=0
          read (luin,*) z,io,ni,ep,au,si,be,s,g
c
c     Check to see if the atom type is one chosen in ATDAT
c
          if (zmap(z).ne.0) then
            atom=zmap(z)
c
c     Check to see if the ion is within the limits set in ATDAT
c
            if (io.le.maxion(atom)) then
c
              ion=io
              ionum=ionum+1
c
              atpho(ionum)=atom
              ionpho(ionum)=ion
              ipotpho(ionum)=ep
              augpho(ionum)=au
              sigpho(ionum)=si*1.0d-18
              betpho(ionum)=be
              spho(ionum)=s
              stwtpho(ionum)=g
              if (ipotpho(ionum).ge.den) goto 40
              write (*,*) ' INCOMPATIBLE data/ionisation/PHIONDAT'
              write (*,*) ionum,ipotpho(ionum),den
              error=.true.
              stop
   40         den=ipotpho(ionum)
            endif
          endif
   50   continue
c
        close (luin)
c
      endif
c
      if (phionmode.eq.0) then
c new verner photo fits and new auger data
        filename=datadir(1:dtlen)//'data/ionisation/PHIONDAT2.txt'
        open (luin,file=filename,status='OLD')
   60   read (luin,fmt=10) (ibuf(j),j=1,19)
        ilgg=ibuf(1)
        if (ilgg(1:1).eq.'%') goto 60
        write (*,20) (ibuf(j),j=1,19)
        read (luin,*) nentries
        ionum=0
        den=0.0d0
        do i=1,nentries
c
          atom=0
          ion=0
c
          read (luin,*) z,io,ni,ph2,shell,ne,fb,wr,wi,rt,ep,pl,au,(p1(j)
     &     ,j=1,5),(p2(k),k=1,7)
c
c     Check to see if the atom type is one chosen in ATDAT
c
          if (zmap(z).ne.0) then
c
            atom=zmap(z)
c
c     Check to see if the ion is within the limits set in ATDAT
c
            if (io.le.maxion(atom)) then
c
              ion=io
              ionum=ionum+1
              if (ionum.gt.mxnegdes) then
                write (*,*) 'ERROR in PHIONDAT2.txt'
                write (*,*) 'Max number of edges exceeded.'
                write (*,*) 'Change parameter mxnegdes in const.inc'
                write (*,*) 'eg: 696 for 16 elements'
                write (*,*) '    1696 for 30 elements'
                write (*,*) 'rebuild and try again.'
                write (*,*) ionum,' Edges.'
                stop
              endif
c
              atpho(ionum)=atom
              ionpho(ionum)=ion
              ipotpho(ionum)=ep
              ph2limitpho(ionum)=pl
              augpho(ionum)=au
c
              wrpho(ionum)=wr
              wipho(ionum)=wi
              rt=dble(wr)/dble(wi)
              stwtpho(ionum)=rt
c
              shellpho(ionum)=shell
              lquantpho(ionum)=lnumber(shell)
c            write(*,*)shell,lnumber(shell)
              nelecpho(ionum)=ne
              hasph2(ionum)=ph2
              isfbpho(ionum)=fb
c
              do j=1,5
                ph1pho(j,ionum)=p1(j)
              enddo
c
              do j=1,7
                ph2pho(j,ionum)=p2(j)
              enddo
c fill in old sigma with new crossection for tau end calculations
              sigpho(ionum)=vernerphoto(ep,ionum)
c              betpho(ionum)=0.d0
c              spho(ionum)=0.d0
c
              if (ipotpho(ionum).ge.den) goto 70
              write (*,*) ' INCOMPATIBLE data/ionisation/PHIONDAT2'
              write (*,*) ' Not in ascending energy order:'
              write (*,*) ionum,ipotpho(ionum),den
              error=.true.
              stop
   70         den=ipotpho(ionum)
            endif
          endif
        enddo
c
        close (luin)
c
      endif
c
c      write(*,*) ' PHIONDAT ',ionum,' edges read.'
c
c Synch phion data with ipotev
c
      do i=1,ionum
        found(i)=.false.
      enddo
c
      do i=1,ionum
c
        atom=atpho(i)
        ion=ionpho(i)
        ep=ipotpho(i)
c
        edge=dabs(ep-ipotev(ion,atom))/ipotev(ion,atom)
c
        if ((edge.lt.5.0d-3).and.(found(i).eqv..false.)) then
          ipotpho(i)=ipotev(ion,atom)
          found(i)=.true.
        endif
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readcolldata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 luin
      integer*4 at,io,ni,is
      integer*4 nions,atom,ion
      integer*4 i,j,k,l
      character ilgg*4, ibuf(19)*4
      logical error
c
      real*8 diffe,ep,tp,a(20)
      real*8 tmin,tmax,xmin,xmax
      real*8 f,dx,x0,x1
      real*8 x(20),y(20),y2(20)
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/ionisation/COLLDAT2.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,10) (ibuf(i),i=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(i),i=1,19)
      read (luin,*) nions
c
      ncol=0
c
      do 50 i=1,nions
c
        atom=0
        ion=0
c
        read (luin,*) at,io,is,ni,ep,tp,(a(l),l=1,20)
c
        if (zmap(at).ne.0) then
c
          atom=zmap(at)
c
          if (io.lt.maxion(atom)) then
c
            ncol=ncol+1
c
            if (ncol.gt.mxncol) then
              write (*,*) 'ERROR in COLLDAT2.txt'
              write (*,*) 'Max number of ions exceeded.'
              write (*,*) 'Change parameter mxncol in const.inc'
              write (*,*) 'eg: 209 for 16 elements'
              write (*,*) '    465 for 20 elements'
              write (*,*) 'rebuild and try again.'
              write (*,*) ncol,' ions.'
              stop
            endif
c
c ionisation rates belong to the lower ion, un like recombinations
c
            ion=io
c
            atcol(ncol)=atom
            ioncol(ncol)=ion
            diffe=(2.d0*(ep-ipotev(ion,atom))/(ep+ipotev(ion,atom)))
            if (dabs(diffe).gt.1e-3) then
              write (*,*) 'ERROR; bad potential in COLLDAT2.txt'
              write (*,*) 'LINE:',i
              write (*,*) 'Atom Ion:',at,io,atom,ion
              write (*,*) 'EP ipotev:',ep,ipotev(ion,atom)
              write (*,*) 'Make it consistent with IONDAT.txt'
              error=.true.
              stop
            endif
c
            ipotcol(ncol)=ipotev(ion,atom)
c
c  convert ipotcol from eV to scaled inverse version for computing x
c  and speeding up main code
c
            tmin=tp*rkb/(ipotcol(ncol)*ev)
            tmax=1.0d9*rkb/(ipotcol(ncol)*ev)
            dx=1.d0/20.d0
            f=2.d0
c
            xmin=1.d0-(dlog(f)/dlog(tmin+f))
            xmax=1.d0-(dlog(f)/dlog(tmax+f))
c
c create x coords array with ni entries
c
            k=1
c 0.d0 - 1.d0 loop, in 0.05 steps, incl 0 and 1
            do j=1,21
              x0=dx*(j-1)
              x1=dx*j
              if (x0.lt.xmin) then
                x(k)=xmin
              elseif (x0.gt.xmax) then
                x(k)=xmax
                goto 40
              else
                x(k)=x0
              endif
              if (x1.ge.xmin) k=k+1
            enddo
   40       continue
c
c debug check on x(i) construction
c
            if (k.ne.ni) then
              write (*,*) 'ERROR; bad x vector '
              write (*,*) 'LINE:',i
              write (*,*) 'Atom Ion:',at,io
              write (*,*) k,ni
              do j=1,k
                write (*,*) j,x(j)
              enddo
              error=.true.
              stop
            endif
c
c get y values and make free splines
c
            nsplcol(ncol)=ni
            do k=1,ni
              y(k)=a(k)
            enddo
c            call endgrads (x, y, ni, dydx0, dydxn)
            call spline (x, y, ni, 1.0d30, 1.0d30, y2)
c save..
            do k=1,ni
              xcol(k,ncol)=x(k)
              ycol(k,ncol)=y(k)
              y2col(k,ncol)=y2(k)
            enddo
c
          endif
        endif
   50 continue
c
      close (luin)
c
c   Precompute splines
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readrecomdata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a1,a2,a3,a4,a5,a6
      integer*4 luin
      integer*4 at,io,ni,tp
      integer*4 nions,i,j,atom,ion
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
c
      filename=datadir(1:dtlen)//'data/ionisation/RRECOMDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,10) (ibuf(i),i=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(i),i=1,19)
      read (luin,*) nions
c
      nrec=0
c
      do 40 i=1,nions
c
        atom=0
        ion=0
        read (luin,*) at,io,ni,tp,a1,a2,a3,a4,a5,a6
        if (zmap(at).ne.0) then
c
          atom=zmap(at)
c
          if (io.lt.maxion(atom)) then
c
            nrec=nrec+1
c
            ion=io+1
            if (nrec.gt.mxnrec) then
              write (*,*) 'ERROR in RRECOMDAT.txt'
              write (*,*) 'Max number of ions exceeded.'
              write (*,*) 'Change parameter mxnrec in const.inc'
              write (*,*) 'eg: 209 for 16 elements'
              write (*,*) '    465 for 20 elements'
              write (*,*) 'rebuild and try again.'
              write (*,*) nrec,' ions.'
              stop
            endif
c
            atrec(nrec)=atom
            ionrec(nrec)=ion
            typerrec(nrec)=tp
c
            arrec(1,nrec)=a1
            arrec(2,nrec)=a2
            arrec(3,nrec)=a3
            arrec(4,nrec)=a4
            arrec(5,nrec)=a5
            arrec(6,nrec)=a6
c
          endif
        endif
   40 continue
c
      close (luin)
c
      nnorad=0
      do i=1,mxion
        do j=1,mxelem
          noradid(i,j)=0
        enddo
      enddo
c
c      if (recommode.eq.2) then
        filename=datadir(1:dtlen)//'data/ionisation/NORADREC.txt'
        open (luin,file=filename,status='OLD')
   50   read (luin,10) (ibuf(i),i=1,19)
        ilgg=ibuf(1)
        if (ilgg(1:1).eq.'%') goto 50
        write (*,20) (ibuf(i),i=1,19)
        read (luin,*) nions
c
        do 60 i=1,nions
c
          atom=0
          ion=0
          read (luin,*) at,io,tp,ni
          if (zmap(at).ne.0) then
c
            atom=zmap(at)
c
            if (io.lt.maxion(atom)) then
c
              nnorad=nnorad+1
c
              ion=io+1
c
              atndrec(nnorad)=atom
              ionndrec(nnorad)=ion
              noradid(ion,atom)=nnorad
              typendrec(nnorad)=tp
              nodesndrec(nnorad)=ni
c
              read (luin,*) (tendrec(j,nnorad),j=1,ni)
              read (luin,*) (andrec(j,nnorad),j=1,ni)
              read (luin,*) (a2ndrec(j,nnorad),j=1,ni)
c
            else
              read (luin,10) (ibuf(j),j=1,19)
              read (luin,10) (ibuf(j),j=1,19)
              read (luin,10) (ibuf(j),j=1,19)
            endif

          else
              read (luin,10) (ibuf(j),j=1,19)
              read (luin,10) (ibuf(j),j=1,19)
              read (luin,10) (ibuf(j),j=1,19)
          endif
   60   continue
c
        close (luin)
c      endif
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readdirecomdata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 en(8),c(8)
      integer*4 luin,atom,ion
      integer*4 at,io,ni,tp
      integer*4 nions,i,j,l,m,ndr
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/ionisation/DRECOMDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,10) (ibuf(i),i=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(i),i=1,19)
      read (luin,*) nions
c
      ndr=0
      do 40 i=1,nions
c
        atom=0
        ion=0
        read (luin,*) at,io,ni,tp,(c(l),l=1,8),(en(m),m=1,8)
c
        if (zmap(at).ne.0) then
c
          atom=zmap(at)
c
          if (io.lt.maxion(atom)) then
c
            ion=io+1
            ndr=ndr+1
c
            if (ndr.gt.mxnrec) then
              write (*,*) 'ERROR in DRECOMDAT.txt'
              write (*,*) 'Max number of ions exceeded.'
              write (*,*) 'Change parameter mxnrec in const.inc'
              write (*,*) 'eg: 209 for 16 elements'
              write (*,*) '    465 for 30 elements'
              write (*,*) 'rebuild and try again.'
              write (*,*) ndr,' ions.'
              stop
            endif
c
            if (atrec(ndr).ne.atom) then
              write (*,*) 'ERROR: data/ionisation/RRECOMDAT and DRECOMDA
     &T'
              write (*,*) '      inconsistent ion in DRECOMDAT'
              write (*,*) 'Line :',i,' atom',at,' ion',io
              stop
            endif
            if (ionrec(ndr).ne.ion) then
              write (*,*) 'ERROR: data/ionisation/RRECOMDAT and DRECOMDA
     &T'
              write (*,*) '      inconsistent ion in DRECOMDAT'
              write (*,*) 'Line :',i,' atom',at,' ion',io
              stop
            endif
            atrec(ndr)=atom
            ionrec(ndr)=ion
            typedrec(ndr)=tp
c
            do j=1,8
              adrec(j,ndr)=c(j)
              bdrec(j,ndr)=en(j)
            enddo
c
          endif
        endif
   40 continue
c
      close (luin)
c
      if (ndr.ne.nrec) then
        write (*,*) 'ERROR: data/ionisation/RRECOMDAT and DRECOMDAT'
        stop
      endif
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readhhedata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 frac,ti,lge,lgs
      real*8 coeffs(18)
      real*8 lam,a,b,c,d
      real*8 y(13), y2(13)
      integer*4 loid,upid,id,n
c      real*8 cevbin,enevlo,enevhi
c      real*8 ip,erow0,erow1
c      integer*4 atom,ion,z,n,inl
      integer*4 i,j,k,l,m,series,at
      integer*4 nu,nl,nuref,nlref
      integer*4 level,line,luin,nentries,ni,nd,nj
      character ilgg*4, ibuf(19)*4
      integer*4 l1, l2,idx
      character s1*16, s2*16
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
   30 format(4(1x,i2,1x),8(1pg10.3,1x))
c
      error=.false.
      filename=datadir(1:dtlen)//'data/hydrogenic/HPXDAT.txt'
      open (luin,file=filename,status='OLD')
   40 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 40
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nhipx
      do i=1,nhipx
        read (luin,*) lge,lgs
        hipxlge(i)=lge
        hipxlgs(i)=lgs
      enddo
c
c    create spline coeffs once
c
c      call endgrads (hipxlge, hipxlgs, nhipx, dydx0, dydxn)
      call spline (hipxlge, hipxlgs, nhipx, 1.0d30, 1.0d30, hipxlgs2)
c
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HRECDAT.txt'
      open (luin,file=filename,status='OLD')
   50 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 50
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nhrec
      do i=1,nhrec
        read (luin,*) hreclgt(i),hreclga(i),hreclga1(i)
      enddo
c
c    create spline coeffs once
c
c      call endgrads (hreclgt, hreclga, nhrec, dydx0, dydxn)
      call spline (hreclgt, hreclga, nhrec, 1.0d30, 1.0d30, hreclga2)
c
c      call endgrads (hreclgt, hreclga1, nhrec, dydx0, dydxn)
      call spline (hreclgt, hreclga1, nhrec, 1.0d30, 1.0d30, hreclga12)
c
c      do  i = 1, nhrec
c          write(*,*) hrecomt(i),hreclga(i),hreclga2(i)
c      enddo
c
c   Read emissivity data
c
      read (luin,10) (ibuf(k),k=1,19)
c
c   Densities
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nd
      read (luin,*) (rhhe(k),k=1,nd)
      do k=1,nd
        rhhe(k)=dlog10(rhhe(k))
      enddo
c
c   Case A Hydrogen Hbeta
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) th42(i),(h42a(i,k),k=1,nd)
        do k=1,nd
          h42a(i,k)=dlog10(h42a(i,k)+epsilon)
        enddo
      enddo
c
c   Case B Hydrogen Hbeta
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) th42(i),(h42b(i,k),k=1,nd)
        do k=1,nd
          h42b(i,k)=dlog10(h42b(i,k)+epsilon)
        enddo
      enddo
c
      do i=1,ni
        th42(i)=dlog10(th42(i))
      enddo
c
      do k=1,nd
        do i=1,ni
          y(i)=h42a(i,k)
        enddo
        call spline (th42, y, ni, 1.0d30, 1.0d30, y2)
        do i=1,ni
          h42a2(i,k)=y2(i)
        enddo
        do i=1,ni
          y(i)=h42b(i,k)
        enddo
        call spline (th42, y, ni, 1.0d30, 1.0d30, y2)
        do i=1,ni
          h42b2(i,k)=y2(i)
        enddo
      enddo
c
c   Case A Helium II 4686
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) the43(i),(he43a(i,k),k=1,nd)
        do k=1,nd
          he43a(i,k)=dlog10(he43a(i,k)+epsilon)
        enddo
      enddo
c
c   Case B Helium II
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) the43(i),(he43b(i,k),k=1,nd)
        do k=1,nd
          he43b(i,k)=dlog10(he43b(i,k)+epsilon)
        enddo
      enddo
c
      do i=1,ni
        the43(i)=dlog10(the43(i))
      enddo
c
      do k=1,nd
        do i=1,ni
          y(i)=he43a(i,k)
        enddo
        call spline (the43, y, ni, 1.0d30, 1.0d30, y2)
        do i=1,ni
          he43a2(i,k)=y2(i)
        enddo
        do i=1,ni
          y(i)=he43b(i,k)
        enddo
        call spline (the43, y, ni, 1.0d30, 1.0d30, y2)
        do i=1,ni
          he43b2(i,k)=y2(i)
        enddo
      enddo
c
c   Hydrogen wavelengths by line, series: Hbeta = 2,2
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nhlines,nhseries
      if (nlines.gt.mxnhlines) then
        write (*,*) 'ERROR: too many hydrogen spectral lines',nlines
        stop
      endif
      if (nhseries.gt.mxnhs) then
        write (*,*) 'ERROR: too many hydrogen spectral series',nhseries
        stop
      endif
      do i=1,nhlines
        read (luin,*) (hlambda(i,k),k=1,nhseries)
      enddo
      nialines(1,zmap(1))=nialines(1,zmap(1))+(nhlines*nhseries)
c
c   Helium II wavelengths by line, series: H4686 = 1,3
c
c   Series 1 is hydrogen ionizing
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nhelines,nheseries
      if (nhelines.gt.mxnhelines) then
        write (*,*) 'ERROR: too many HeII spectral lines',nhelines
        stop
      endif
      if (nheseries.gt.mxnhes) then
        write (*,*) 'ERROR: too many HeII spectral series',nheseries
        stop
      endif
      do i=1,nhelines
        read (luin,*) (helambda(i,k),k=1,nheseries)
      enddo
      nialines(2,zmap(2))=nialines(2,zmap(2))+(nhelines*nheseries)
c
c
c   Hydrogenic gf values (same for H and He)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni,nj
      do i=1,ni
        read (luin,*) (hydrogf(i,k),k=1,nj)
      enddo
c
      close (luin)
c
      nxhlines=mxnxhlines
      nxhseries=mxnxhs
      do at=1,atypes
        nialines(at,at)=nialines(at,at)+nxhlines*nxhseries
      enddo
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HELIORECDAT.txt'
      open (luin,file=filename,status='OLD')
   60 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 60
      write (*,20) (ibuf(j),j=1,19)
c
c read level energies and level names
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nheis
      do j=1,nheis
        read (luin,*) id,heisej(j),heislvl(j)
        if (id.ne.j) then
          write (*,*) 'ERROR reading HELIORECDAT.txt singlet levels',i,
     &     j
          stop
        endif
      enddo
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nheit
      do j=1,nheit
        read (luin,*) id,heitej(j),heitlvl(j)
        if (id.ne.j) then
          write (*,*) 'ERROR reading HELIORECDAT.txt triplet levels',i,
     &     j
          stop
        endif
      enddo
c
c effective recombination rate coeffs
c
c singlet lines, case B
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nheislines
      do j=1,nheislines
        read (luin,*) upid,loid,lam,a,b,c,d
c
        heiseij(j)=(heisej(upid)-heisej(loid))*cls*plk
        heislam(j)=1.0d8/(heisej(upid)-heisej(loid))!vacritz
        heisupid(j)=upid
        heisloid(j)=loid
        heisreccoef(1,j)=a
        heisreccoef(2,j)=b
        heisreccoef(3,j)=c
        heisreccoef(4,j)=d
        s1=heislvl(upid)
        l1=0
        do idx=1,16
          if (s1(idx:idx).ne.' ') then
            l1=l1+1
          endif
        enddo
        if (l1.lt.1) l1=1
        s2=heislvl(loid)
        l2=0
        do idx=1,16
          if (s2(idx:idx).ne.' ') then
            l2=l2+1
          endif
        enddo
        if (l2.lt.1) l2=1
        heislamid(j)=s1(1:l1)//'-'//s2(1:l2)
c        write(*,*) trim(heislamid(j))
      enddo
      nialines(1,zmap(2))=nialines(1,zmap(2))+nheislines
c
c triplet line recomb coeffs
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nheitlines
      do j=1,nheitlines
        read (luin,*) upid,loid,lam,a,b,c,d
        heiteij(j)=(heitej(upid)-heitej(loid))*cls*plk
        heitlam(j)=1.0d8/(heitej(upid)-heitej(loid))!vacritz
        heitupid(j)=upid
        heitloid(j)=loid
        heitreccoef(1,j)=a
        heitreccoef(2,j)=b
        heitreccoef(3,j)=c
        heitreccoef(4,j)=d
        s1=heitlvl(upid)
        l1=0
        do idx=1,16
          if (s1(idx:idx).ne.' ') then
            l1=l1+1
          endif
        enddo
        if (l1.lt.1) l1=1
        s2=heitlvl(loid)
        l2=0
        do idx=1,16
          if (s2(idx:idx).ne.' ') then
            l2=l2+1
          endif
        enddo
        if (l2.lt.1) l2=1
        heitlamid(j)=s1(1:l1)//'-'//s2(1:l2)
c        write(*,*) trim(heitlamid(j))
      enddo
      nialines(1,zmap(2))=nialines(1,zmap(2))+nheitlines
c
c singlet levels, collisional contribution to upper levels
c nheiscr(14) heiscrA(6 14)  heiscrB(6 14)  heiscrC(6 14)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nentries
      if (nentries.ne.nheis) then
        write (*,*) 'ERROR: HELIORECDAT.txt singlet coll. data',
     &   nentries
        stop
      endif
      do j=1,nheis
        read (luin,*) id,n,(coeffs(i),i=1,18)
        nheitcr(j)=n
        do i=1,6
          heiscra(i,j)=coeffs(i)
          heiscrb(i,j)=coeffs(i+6)
          heiscrc(i,j)=coeffs(i+12)
        enddo
        if (id.ne.j) then
          write (*,*) 'ERROR: HELIORECDAT.txt singlet coll. data',id,j
          stop
        endif
      enddo
c
c triplet levels, collisional contribution to upper levels
c
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nentries
      if (nentries.ne.nheit) then
        write (*,*) 'ERROR: HELIORECDAT.txt triplet coll. data',
     &   nentries
        stop
      endif
      do j=1,nheit
        read (luin,*) id,n,(coeffs(i),i=1,18)
        nheitcr(j)=n
        do i=1,6
          heitcra(i,j)=coeffs(i)
          heitcrb(i,j)=coeffs(i+6)
          heitcrc(i,j)=coeffs(i+12)
        enddo
        if (id.ne.j) then
          write (*,*) 'ERROR: HELIORECDAT.txt triplet coll. data',id,j
          stop
        endif
      enddo
c
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HRECRATA.txt'
      open (luin,file=filename,status='OLD')
   70 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 70
      read (luin,*) nentries
c
      line=1
      series=1
      do j=1,nentries
        read (luin,*) nu,nl,nuref,nlref
        line=nu-nl
        series=nl
        read (luin,*) ni,nj
        do i=1,ni
          read (luin,*) ti,(hylratsa(line,series,i,k),k=1,nj)
          do k=1,nj
            hylratsa(line,series,i,k)=
     &     dlog10(hylratsa(line,series,i,k)+epsilon)
          enddo
        enddo
c precompute orthogonal splines
        do k=1,nj
          do i=1,ni
            y(i)=hylratsa(line,series,i,k)
          enddo
          call spline (th42, y, ni, 1.0d30, 1.0d30, y2)
          do i=1,ni
            hylratsa2(line,series,i,k)=y2(i)
          enddo
        enddo
      enddo
cc
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HRECRATB.txt'
      open (luin,file=filename,status='OLD')
   80 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 80
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
c
      line=1
      series=1
      do j=1,nentries
c        read (luin,10) (ibuf(k),k=1,19)
        read (luin,*) nu,nl,nuref,nlref
        line=nu-nl
        series=nl
        read (luin,*) ni,nj
        do i=1,ni
          read (luin,*) ti,(hylratsb(line,series,i,k),k=1,nj)
          do k=1,nj
            hylratsb(line,series,i,k)=
     &     dlog10(hylratsb(line,series,i,k)+epsilon)
          enddo
        enddo
c precompute orthogonal splines
        do k=1,nj
          do i=1,ni
            y(i)=hylratsb(line,series,i,k)
          enddo
          call spline (th42, y, ni, 1.0d30, 1.0d30, y2)
          do i=1,ni
            hylratsb2(line,series,i,k)=y2(i)
          enddo
        enddo
      enddo
c
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HERECDAT.txt'
      open (luin,file=filename,status='OLD')
   90 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 90
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nherec
      do i=1,nherec
        read (luin,*) hereclgt(i),hereclga(i),hereclga1(i)
      enddo
c
c    create spline coeffs once
c
c      call endgrads (hereclgt, hereclga, nherec, dydx0, dydxn)
      call spline (hereclgt, hereclga, nherec, 1.0d30, 1.0d30,
     &hereclga2)
c
c      call endgrads (hereclgt, hreclga1, nherec, dydx0, dydxn)
      call spline (hereclgt, hreclga1, nherec, 1.0d30, 1.0d30,
     &hereclga12)
c
c      do  i = 1, nhrec
c          write(*,*) herecomt(i),hereclga(i),hereclga2(i)
c      enddo
c
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HERECRATA.txt'
      open (luin,file=filename,status='OLD')
  100 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 100
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
c
      line=1
      series=1
      do j=1,nentries
c        read (luin,10) (ibuf(k),k=1,19)
        read (luin,*) nu,nl,nuref,nlref
        line=nu-nl
        series=nl
        read (luin,*) ni,nj
        do i=1,ni
          read (luin,*) ti,(helratsa(line,series,i,k),k=1,nj)
          do k=1,nj
            helratsa(line,series,i,k)=
     &        dlog10(helratsa(line,series,i,k)+epsilon)
          enddo
        enddo
c precompute orthogonal splines
        do k=1,nj
          do i=1,ni
            y(i)=helratsa(line,series,i,k)
          enddo
          call spline (the43, y, ni, 1.0d30, 1.0d30, y2)
          do i=1,ni
            helratsa2(line,series,i,k)=y2(i)
          enddo
        enddo
      enddo
c
      close (luin)
c
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HERECRATB.txt'
      open (luin,file=filename,status='OLD')
  110 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 110
      read (luin,*) nentries
c
      line=1
      series=1
      do j=1,nentries
c        read (luin,10) (ibuf(k),k=1,19)
        read (luin,*) nu,nl,nuref,nlref
        line=nu-nl
        series=nl
        read (luin,*) ni,nj
        do i=1,ni
          read (luin,*) ti,(helratsb(line,series,i,k),k=1,nj)
          do k=1,nj
            helratsb(line,series,i,k)=
     &          dlog10(helratsb(line,series,i,k)+epsilon)
          enddo
        enddo
c precompute orthogonal splines
        do k=1,nj
          do i=1,ni
            y(i)=helratsb(line,series,i,k)
          enddo
          call spline (the43, y, ni, 1.0d30, 1.0d30, y2)
          do i=1,ni
            helratsb2(line,series,i,k)=y2(i)
          enddo
        enddo
      enddo
c
      close (luin)
c
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HCOLLGNP.txt'
      open (luin,file=filename,status='OLD')
  120 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 120
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,30) (colid(j,i),i=1,4),(ccoln(j,k),k=1,4),(dcoln(j,l)
     &   ,l=1,4)
      enddo
c
c
c   collisional excitation branching ratios
c
c
      read (luin,10) (ibuf(k),k=1,19)
c
c     case A
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) line,series,level,frac
        collha(line,series,level)=frac
      enddo
c
c     case B
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) line,series,level,frac
        collhb(line,series,level)=frac
      enddo
c
      close (luin)
c
c **** INLCLUDES legacy HEI data ****
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HYDRODAT.txt'
      open (luin,file=filename,status='OLD')
  130 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 130
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) (hlam(j),j=1,5)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 140 j=1,5
        read (luin,*) m,(balcoe(j,i),i=1,5)
        if (m.eq.j) goto 140
        write (*,*) 'Error Reading HYDRODAT'
        error=.true.
        stop
  140 continue
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 150 j=1,5
        read (luin,*) m,(balcoe(j,i),i=6,11)
        if (m.eq.j) goto 150
        write (*,*) 'Error Reading HYDRODAT'
        error=.true.
        stop
  150 continue
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 160 j=1,5
        read (luin,*) m,(hbema(j,i),i=1,6)
        if (m.eq.j) goto 160
        write (*,*) 'Error Reading HYDRODAT'
        error=.true.
        stop
  160 continue
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 170 j=1,5
        read (luin,*) m,(hbemb(j,i),i=1,6)
        if (m.eq.j) goto 170
        write (*,*) 'Error Reading HYDRODAT'
        error=.true.
        stop
  170 continue
c
c Legacy He I recob data, input changed to cm^-1 and A, reconvert here
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) (heien(i),i=1,4)!newreadincm^-1
      heilam(1)=1.0d0/heien(1)!vacritzlambda,cm
      heilam(2)=1.0d0/heien(3)
      heilam(3)=1.0d0/heien(4)
      heien(1)=heien(1)*plk*cls!nowinbackinergs
      heien(2)=heien(2)*plk*cls
      heien(3)=heien(3)*plk*cls
      heien(4)=heien(4)*plk*cls
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) (qren(i),i=1,4)!newreadincm^-1
      qren(1)=qren(1)*plk*cls!nowinbackinergs
      qren(2)=qren(2)*plk*cls
      qren(3)=qren(3)*plk*cls
      qren(4)=qren(4)*plk*cls
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do j=1,12
        read (luin,*) (hel0(i,j),i=1,5)
        hel0(1,j)=hel0(1,j)*1.0d-8!nowvacritzcm
        heilam(j+3)=hel0(1,j)!nowvacritzcm
      enddo
      nheilines=15
      close (luin)
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HYDROFBDAT.txt'
      open (luin,file=filename,status='OLD')
  180 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 180
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nfbe
      do j=1,nfbe
        read (luin,*) nfbenu(j),fbenu(j),fben(j)
      enddo
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      read (luin,*) (fblte(i),i=1,ni)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni,nj
c
      do j=1,nj
        read (luin,*) (fbloggamma(i,j),i=1,ni)
      enddo
c
      do j=1,nj
        do i=1,ni
          fbloggamma(i,j)=dlog10(fbloggamma(i,j))
        enddo
      enddo
c
      close (luin)
c
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HELIOFBDAT.txt'
      open (luin,file=filename,status='OLD')
  190 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 190
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) nhefbe
      do j=1,nhefbe
        read (luin,*) nhefbenu(j),hefbenu(j),hefben(j)
      enddo
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      read (luin,*) (hefblte(i),i=1,ni)
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni,nj
c
      do j=1,nj
        read (luin,*) (hefbloggamma(i,j),i=1,ni)
      enddo
c
      do j=1,nj
        do i=1,ni
          hefbloggamma(i,j)=dlog10(hefbloggamma(i,j))
        enddo
      enddo
c
      close (luin)
c
      j108m=0
      k=6
      do j=1,12
        if ((dabs(hel0(1,j)-10830d-8)/hel0(1,j)).lt.3.d-2) j108m=j
        if ((hel0(2,j).gt.0.d0).and.(j108m.ne.j)) then
          hlam(k)=hel0(1,j)
          k=min0(10,k+1)
        endif
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readhhecoldata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 luin
      logical error
c
      integer*4 nentries,idx
      character ilgg*4, ibuf(19)*4
c
      integer*4 id,at,io,ni
      integer*4 nl,i,j
      real*8 ei,gi,tc,aij,gfij
      character*24 termi
      integer*4 nc,idc,tp,nspl,ids
c
c      real*8 hheeij (mxhhelvls,mxhhelvls)  ! cm^-1
c      real*8 hheaji (mxhhelvls,mxhhelvls)
c      real*8 hhegfij(mxhhelvls,mxhhelvls)
c
      real*8 btx (mxhhenspl)
      real*8 bty (mxhhenspl)
      real*8 bty2(mxhhenspl)
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/hydrogenic/HCOLLDATA.txt'
      open (luin,file=filename,status='OLD')
c
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
      do idx=1,mxelem
        ishheion(idx)=0
      enddo
c
      nhheions=0
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) id
        read (luin,*) at,io
        if (id.ne.idx) then
          write (*,*) 'ERROR READING HCOLLDATA ION: ',idx,id
          stop
        endif
c       atom=zmap(na)
c       if (atom.ne.0) then
c         if (ion.le.maxion(atom)) then
        nhheions=nentries
        hheat(idx)=zmap(at)
        hheion(idx)=io
        ishheion(zmap(at))=1
c
c get level energies, gs and names
c
        read (luin,*) nl
        if (nl.gt.mxhhelvls) then
          write (*,*) 'ERROR IN HCOLLDATA Level Limit Exceeded: '
          write (*,*) idx,nl,mxhhelvls
          write (*,*) mapz(hheat(idx)),hheion(idx)
          stop
        endif
        hheni(idx)=nl
        do i=1,nl
          read (luin,*) j,ni,ei,gi,termi
          if (i.ne.j) then
            write (*,*) 'ERROR IN HCOLLDATA LVLS: ',idx,i,j
            stop
          endif
          hheei(i,idx)=ei
          hhen(i,idx)=ni
          hhegi(i,idx)=gi
          hheinvgi(i,idx)=1.d0/gi
          hheid(i,idx)=termi(1:16)
c
          if (i.eq.2) then
            ee2p(hheat(idx))=ei*plk*cls
          endif
c
        enddo
c
c i lower, j upper
c       do i=1,nl
c         do j=1,nl
c           hheegij(i,j,idx)=plk*cls*dabs(hheei(j,idx)-hheei(i,idx))
c         enddo
c       enddo
c
c get collisions
c
        read (luin,*) nc
        if (nc.gt.mxhhecols) then
          write (*,*) 'ERROR IN HCOLLDATA Collisions Limit Exceeded: '
          write (*,*) idx,nc,mxhhecols
          write (*,*) mapz(hheat(idx)),hheion(idx)
          stop
        endif
        nhheioncol(idx)=nc
        do idc=1,nc
c          read (luin,*) id,i,j,tp,nspl,tc
          read (luin,*) id,i,j,tp,nspl,aij,gfij,tc
          if (id.ne.idc) then
            write (*,*) 'ERROR IN HCOLLDATA COLL: ',idx,id,idc
            stop
          endif
          hhecol_i(idc,idx)=i
          hhecol_j(idc,idx)=j
          hhecol_tc(idc,idx)=tc
          hhecol_typespl(idc,idx)=tp
          hhecol_nspl(idc,idx)=nspl
          read (luin,*) (bty(ids),ids=1,nspl)
          do ids=1,nspl
            btx(ids)=(1.d0/(nspl-1))*(ids-1)
            bty2(ids)=0.d0
          enddo
c
          call spline (btx, bty, nspl, 1.d30, 1.d30, bty2)
c
          do ids=1,nspl
            hhecol_x(ids)=btx(ids)
            hhecol_y(ids,idc,idx)=bty(ids)
            hhecol_y2(ids,idc,idx)=bty2(ids)
          enddo
        enddo
      enddo
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readphotdat (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 edge1,edge2,ep
      integer*4 i,j,luin,nentries
      integer*4 atom, ion, is, z
      character ilgg*4, ibuf(19)*4
      logical error, found(1696)
c
      real*8 eleft,e,eright
      real*8 c0,c1,c2
      real*8 x,lge
c
      integer inl
      real*8 ip,erow1,erow0,cevbin,enevlo,enevhi
c
      real*8 be, si, s, invpot
c
      real*8 acrs,eph,at,bet,se
c
c     internal function
c
      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c     external functions
c
      real*8 vernerphoto,fsplint
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
      do i=1,ionum
        found(i)=.false.
      enddo
      filename=datadir(1:dtlen)//'data/PHOTDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
      fieldversion=0
      ionstartbin=1
      read (luin,*) nentries
      if (nentries.lt.0) then
        fieldversion=nentries
        read (luin,*) nentries
      endif
c
      if (nentries.gt.mxinfph) then
        write (*,*) ' INCOMPATIBLE data/PHOTDAT.txt, need',mxinfph,' bin
     &s or less.'
        error=.true.
        stop
      endif
c
      infph=nentries
c
      do j=1,infph
        read (luin,*) photev(j)
      enddo
c
      close (luin)
c
c Synchronize any older PHOTDATs with any newer potential data.
c photev is synchronized with ipotpho, which was synched with ipotev
c the reference
c
      do j=1,infph
        do i=1,ionum
c
c          atom=atpho(i)
c          ion=ionpho(i)
          ep=ipotpho(i)
c
          if ((ep.ge.photev(j)).and.(ep.lt.photev(j+1))) then
            edge1=dabs(ep-photev(j))/photev(j)
            edge2=dabs(ep-photev(j+1))/photev(j+1)
c
            if (edge1.le.edge2) then
              if ((edge1.lt.5.0d-6).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j),ep,edge1
                photev(j)=ep
                phote(j)=ep*ev
                found(i)=.true.
              endif
            endif
c
            if (edge2.lt.edge1) then
              if ((edge2.lt.5.0d-6).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j+1),ep,edge2
                photev(j+1)=ep
                phote(j+1)=ep*ev
                found(i)=.true.
              endif
            endif
c
          endif
c
        enddo
      enddo
c
      do j=1,infph
        do i=1,ionum
c          atom=atpho(i)
c          ion=ionpho(i)
          ep=ipotpho(i)
          if ((ep.ge.photev(j)).and.(ep.lt.photev(j+1))) then
            edge1=dabs(ep-photev(j))/photev(j)
            edge2=dabs(ep-photev(j+1))/photev(j+1)
c
            if (edge1.le.edge2) then
              if ((edge1.lt.5.0d-5).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j),ep,edge1
                photev(j)=ep
                phote(j)=ep*ev
                found(i)=.true.
              endif
            endif
c
            if (edge2.lt.edge1) then
              if ((edge2.lt.5.0d-5).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j+1),ep,edge2
                photev(j+1)=ep
                phote(j+1)=ep*ev
                found(i)=.true.
              endif
            endif
c
          endif
c
        enddo
      enddo
c
      do j=1,infph
        do i=1,ionum
c          atom=atpho(i)
c          ion=ionpho(i)
          ep=ipotpho(i)
          if ((ep.ge.photev(j)).and.(ep.lt.photev(j+1))) then
            edge1=dabs(ep-photev(j))/photev(j)
            edge2=dabs(ep-photev(j+1))/photev(j+1)
c
            if (edge1.le.edge2) then
              if ((edge1.lt.5.0d-4).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j),ep,edge1
                photev(j)=ep
                phote(j)=ep*ev
                found(i)=.true.
              endif
            endif
c
            if (edge2.lt.edge1) then
              if ((edge2.lt.5.0d-4).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j+1),ep,edge2
                photev(j+1)=ep
                phote(j+1)=ep*ev
                found(i)=.true.
              endif
            endif
c
          endif
c
        enddo
      enddo
c
      do j=1,infph
        do i=1,ionum
c          atom=atpho(i)
c          ion=ionpho(i)
          ep=ipotpho(i)
          if ((ep.ge.photev(j)).and.(ep.lt.photev(j+1))) then
            edge1=dabs(ep-photev(j))/photev(j)
            edge2=dabs(ep-photev(j+1))/photev(j+1)
c
            if (edge1.le.edge2) then
              if ((edge1.lt.5.0d-3).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j),ep,edge1
                photev(j)=ep
                phote(j)=ep*ev
                found(i)=.true.
              endif
            endif
c
            if (edge2.lt.edge1) then
              if ((edge2.lt.5.0d-3).and.(found(i).eqv..false.)) then
c              write(*,*) mapz(atom),ion,photev(j+1),ep,edge2
                photev(j+1)=ep
                phote(j+1)=ep*ev
                found(i)=.true.
              endif
            endif
c
          endif
c
        enddo
      enddo
c
      do i=1,ionum
        atom=atpho(i)
        ion=ionpho(i)
        ep=ipotpho(i)
        if (found(i).eqv..false.) then
          write (*,*) 'WARNING IN FILE : PHOTDAT'
          write (*,*) 'Potential Not found',mapz(atom),ion,ep
          write (*,*) 'Using nearest bin edge...'
        endif
      enddo
c
      do 50 j=1,infph-1
        if (photev(j).lt.photev(j+1)) goto 50
   40 format(//'ERROR IN FILE : PHOTDAT @#',i4,2(x,1pg14.7))
        write (*,40) j,photev(j-1),photev(j),photev(j+1)
        error=.true.
        stop
   50 continue
c
      do j=1,infph-1
        lgphotev(j)=dlog10(photev(j))
        lgphote(j)=dlog10(phote(j))
        cphotev(j)=0.5d0*(photev(j)+photev(j+1))
        cphote(j)=cphotev(j)*ev
        lgcphotev(j)=dlog10(cphotev(j))
        lgcphote(j)=dlog10(cphote(j))
        widbinev(j)=photev(j+1)-photev(j)
      enddo
      lgphotev(infph)=dlog10(photev(infph))
      lgphote(infph)=dlog10(phote(infph))
      widbinev(infph)=0.d0
c
c     initialise skipbin
c
      do i=1,infph
        skipbin(i)=.true.
      enddo
c
c     precompute photo cross sections
c
      if (phionmode.eq.0) then
c
        do i=1,ionum
          ep=ipotpho(i)
          invpot=1.d0/ipotpho(i)
          atom=atpho(i)
          ion=ionpho(i)
          is=(mapz(atom)-ion+1)
          do j=1,infph-1
            eleft=photev(j)
            eright=photev(j+1)
            e=cphotev(j)
            x=0.d0
            if (eright.gt.ep) then
              if (eleft.lt.ep) then
                eleft=ep
                e=0.5d0*(eleft+eright)
              endif
c
              c0=vernerphoto(eleft,i)
              c1=vernerphoto(e,i)
              c2=vernerphoto(eright,i)
c  3rd order PPM/Simpsons integral average over bin
              x=(c0+(4.d0*c1)+c2)/6.d0
              if (is.eq.1) then
c
c                Hydrogenic Cross Sections from spline
c                cubic pline in log s - logE/E0 space
c
                lge=dlog10(eleft/ipotpho(i))
                c0=fsplint(hipxlge,hipxlgs,hipxlgs2,nhipx,lge)
                c0=10.d0**(c0)
                lge=dlog10(e/ipotpho(i))
                c1=fsplint(hipxlge,hipxlgs,hipxlgs2,nhipx,lge)
                c1=10.d0**(c1)
                lge=dlog10(eright/ipotpho(i))
                c2=fsplint(hipxlge,hipxlgs,hipxlgs2,nhipx,lge)
                c2=10.d0**(c2)
c  3rd order PPM/Simpsons integral average over bin
                x=1.0d-18*(c0+(4.d0*c1)+c2)/(ion*ion*6.d0)
              endif
c
            endif
            photxsec(i,j)=x
          enddo
          photxsec(i,infph)=photxsec(i,infph-1)
        enddo
      endif
c
      if (phionmode.eq.1) then
c
        do 70 i=1,ionum
          be=betpho(i)
          si=sigpho(i)
          s=spho(i)
          atom=atpho(i)
          ion=ionpho(i)
          invpot=1.d0/ipotpho(i)
          do 60 j=1,infph-1
            eleft=photev(j)*invpot
            eright=photev(j+1)*invpot
            e=cphotev(j)*invpot
            x=0.d0
            if (eright.gt.1.d0) then
              if (eleft.lt.1.d0) then
                eleft=1.d0
                e=0.5d0*(eleft+eright)
              endif
              c0=acrs(eleft,si,be,s)
              c1=acrs(e,si,be,s)
              c2=acrs(eright,si,be,s)
c  3rd order PPM/Simpsons integral average over bin
              x=(c0+(4.d0*c1)+c2)/6.d0
            endif
            photxsec(i,j)=x
   60     continue
          photxsec(i,infph)=photxsec(i,infph-1)
   70   continue
      endif
c
c prepare free bound n>= 2 edges
c
c
c set up energy vector pre-calculated arrays
c
      do atom=1,atypes
        z=mapz(atom)
        ion=z
c
c fbee = electron free energy = for exp[-ee/kT] T in K
c nfb  = n for dominant edge
c jfb  = lower index for row in fbloggamma table
c
        ip=ipotev(ion,atom)
        enevlo=ip*fbenu(1)
        enevhi=ip*fbenu(nfbe)
        do 100 inl=1,infph
          cevbin=cphotev(inl)
          jfb(atom,inl)=0
          if ((cevbin.gt.enevlo).and.(cevbin.lt.enevhi)) then
            do 80 j=1,(nfbe-1)
              erow0=ip*fbenu(j)
              erow1=ip*fbenu(j+1)
              if ((cevbin.ge.erow0).and.(cevbin.lt.erow1)) then
c                n=nfbenu(j)
                jfb(atom,inl)=j
                fbee(atom,inl)=(cevbin-(ip*fben(j)))*ev
c      if (atom.eq.1) then
c      write(*,*) 'HI',inl,cevbin,n,j,(ip*fben(j)),(cevbin-(ip*fben(j)))
c      write(*,*) 'HI',inl,cevbin,n,j,(ip/dble(n*n)),(cevbin-(ip/dble(n*n)))
c      endif
                goto 90
              endif
   80       continue
          endif
   90     continue
  100   continue
c
      enddo
c
      do j=1,nfbe
        fbenu(j)=dlog10(fbenu(j))
      enddo
c
c
c set up energy vector pre-calculated arrays
c
      do atom=1,atypes
        z=mapz(atom)
        ion=z-1
        if (ion.ge.1) then
c
c fbee = electron free energy = for exp[-ee/kT] T in K
c nfb  = n for dominant edge
c jfb  = lower index for row in fbloggamma table
c
          ip=ipotev(ion,atom)
          enevlo=ip*hefbenu(1)
          enevhi=ip*hefbenu(nhefbe)
          do 130 inl=1,infph
            cevbin=cphotev(inl)
            jhefb(atom,inl)=0
            if ((cevbin.gt.enevlo).and.(cevbin.lt.enevhi)) then
              do 110 j=1,(nhefbe-1)
                erow0=ip*hefbenu(j)
                erow1=ip*hefbenu(j+1)
                if ((cevbin.ge.erow0).and.(cevbin.lt.erow1)) then
                  jhefb(atom,inl)=j
                  hefbee(atom,inl)=(cevbin-(ip*hefben(j)))*ev
c         if (atom.eq.2) then
c         write(*,*) 'HeI',inl,cevbin,n,j,(ip*hefben(j)),(cevbin-(ip*hefben(j)))
c         write(*,*) 'HeI',inl,cevbin,n,j,(ip/dble(n*n)),(cevbin-(ip/dble(n*n)))
c         endif
                  goto 120
                endif
  110         continue
            endif
  120       continue
  130     continue
        endif
      enddo
c
      do j=1,nhefbe
        hefbenu(j)=dlog10(hefbenu(j))
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readchx (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 deltae
      real*8 tl0,tl1,ax,bx,cx,dx
      integer*4 at,atom,i,io,ion,j,luin,nentries
      integer*4 ni,nsequ,nft,id
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/ionisation/CHXDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
c     read new charge exchange reactions
c
      nchxr=0
      nchxi=0
c
      read (luin,*) nsequ
      do i=1,nsequ
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do j=1,nentries
c
          atom=0
          ion=0
c
          read (luin,*) ni,at,io,nft,id,tl0,tl1,ax,bx,cx,dx,deltae
c
          if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
            atom=zmap(at)
            if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
              if (i.lt.3) then
c
c     recombination rates
c
                nchxr=nchxr+1
                chxrx(nchxr)=ni
                chxrat(nchxr)=atom
                chxrio(nchxr)=io
                chxrnfit(nchxr)=nft
                chxrfitid(nchxr)=id
                chxrtemp(1,nchxr)=tl0
                chxrtemp(2,nchxr)=tl1
                chxrcos(1,nchxr)=ax*1.0d-9
                chxrcos(2,nchxr)=bx
                chxrcos(3,nchxr)=cx
                chxrcos(4,nchxr)=dx
c ev4
                chxrde(nchxr)=deltae
                chxr(nchxr)=0.d0
              else
c
c     ionisation rates
c
                nchxi=nchxi+1
                chxix(nchxi)=ni
                chxiat(nchxi)=atom
                chxiio(nchxi)=io
                chxinfit(nchxi)=nft
                chxifitid(nchxi)=id
                chxitemp(1,nchxi)=tl0
                chxitemp(2,nchxi)=tl1
                chxicos(1,nchxi)=ax*1.0d-9
                chxicos(2,nchxi)=bx
                chxicos(3,nchxi)=cx
                chxicos(4,nchxi)=dx
c ev4
                chxicos(5,nchxi)=deltae
                chxide(nchxi)=deltae
                chxi(nchxi)=0.d0
              endif
c
c     end inclusion ifs
c
            endif
          endif
c
c     end nsequ,nentries
c
        enddo
      enddo
c
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readion2 (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 aa,ac,ad,ar,bd
      real*8 ta,tc,ti,to,xr
      integer*4 at,atom,io,ion,j,luin
      integer*4 ni,nions,nsequ,nsecel
      character ilgg*4, ibuf(19)*4
      logical error
c
      real*8 lt(8)
      real*8 alp1,alp2,hhe
      integer*4 i,m,nentries
      real*8 sg1,sg2
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/ionisation/IONDAT2.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nions
c
      do 40 j=1,nions
        atom=0
        ion=0
        read (luin,*) at,io,ni,ac,tc,ar,xr,ad,bd,to,ti,aa,ta
        if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
          atom=zmap(at)
          if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
            ion=io
            acol(ion,atom)=ac
            tcol(ion,atom)=tc
c
            ion=io+1
            arad(ion,atom)=ar
            xrad(ion,atom)=xr
c
c     watch where you're coming from, not where you're going...
c     recom rates belong to the n+1 ionisation stage
c
            adi(ion,atom)=ad
            bdi(ion,atom)=bd
            t0(ion,atom)=to
            t1(ion,atom)=ti
c
          endif
        endif
   40 continue
c
c     read in misc ionisation rate stuff
c
c     H & He rates plus entries for on the spot approx
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(arad(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(arad(j,zmap(2)),j=2,6)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(xrad(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(xrad(j,zmap(2)),j=2,6)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(adi(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(adi(j,zmap(2)),j=2,6)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(t0(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(t0(j,zmap(2)),j=2,6)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(bdi(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(bdi(j,zmap(2)),j=2,6)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      if (zmap(1).ne.0) read (luin,*) m,(t1(j,zmap(1)),j=2,6)
      if (zmap(2).ne.0) read (luin,*) m,(t1(j,zmap(2)),j=2,6)
c
c     low T dielec correction
c      check that non set entries are zero
c
      do atom=1,atypes
        do ion=1,maxion(atom)
          adilt(ion,atom)=0.d0
          bdilt(ion,atom)=0.d0
          cdilt(ion,atom)=0.d0
          ddilt(ion,atom)=0.d0
          fdilt(ion,atom)=0.d0
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) at,(lt(j),j=2,8)
        if (zmap(at).ne.0) then
          atom=zmap(at)
          do ion=2,8
            adilt(ion,atom)=lt(ion)
          enddo
        endif
      enddo
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) at,(lt(j),j=2,8)
        if (zmap(at).ne.0) then
          atom=zmap(at)
          do ion=2,8
            bdilt(ion,atom)=lt(ion)
          enddo
        endif
      enddo
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) at,(lt(j),j=2,8)
        if (zmap(at).ne.0) then
          atom=zmap(at)
          do ion=2,8
            cdilt(ion,atom)=lt(ion)
          enddo
        endif
      enddo
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) at,(lt(j),j=2,8)
        if (zmap(at).ne.0) then
          atom=zmap(at)
          do ion=2,8
            ddilt(ion,atom)=lt(ion)
          enddo
        endif
      enddo
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) at,(lt(j),j=2,8)
        if (zmap(at).ne.0) then
          atom=zmap(at)
          do ion=2,8
            fdilt(ion,atom)=lt(ion)
          enddo
        endif
      enddo
c
c     Ancient and Old charge exchange and secondary ionisation data
c     read Ancient charge exchange reactions, chargemode=2
c
      nchxold=0
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) m,at,io,sg1,sg2,alp1,alp2,hhe
        if (zmap(at).ne.0) then
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            nchxold=nchxold+1
            charco(1,nchxold)=dble(atom)
            charco(2,nchxold)=dble(io)
            charco(3,nchxold)=sg1
            charco(4,nchxold)=sg2
            charco(5,nchxold)=alp1
            charco(6,nchxold)=alp2
            charco(7,nchxold)=dble(hhe)
          endif
        endif
      enddo
c
c     read legacy charge exchange reactions, chargemode=1
c
      nlegacychxr=0
      nlegacychxi=0
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nsequ
      do i=1,nsequ
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do j=1,nentries
          atom=0
          ion=0
c
c     just lazy, reusing some spare locals...
c
          read (luin,*) ni,at,io,to,ti,ac,tc,ar,xr
          if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
            atom=zmap(at)
            if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
              if (i.lt.3) then
c
c     recombination rates
c
                nlegacychxr=nlegacychxr+1
                chxrlegacyx(nlegacychxr)=ni
                chxrlegacyat(nlegacychxr)=atom
                chxrlegacyio(nlegacychxr)=io
                chxrlegacytemp(1,nlegacychxr)=to
                chxrlegacytemp(2,nlegacychxr)=ti
                chxrlegacycos(1,nlegacychxr)=ac*1.0d-9
                chxrlegacycos(2,nlegacychxr)=tc
                chxrlegacycos(3,nlegacychxr)=ar
                chxrlegacycos(4,nlegacychxr)=xr
                chxrlegacy(nlegacychxr)=0.d0
              else
c
c     ionisation rates
c
                nlegacychxi=nlegacychxi+1
                chxilegacyx(nlegacychxi)=ni
                chxilegacyat(nlegacychxi)=atom
                chxilegacyio(nlegacychxi)=io
                chxilegacytemp(1,nlegacychxi)=to
                chxilegacytemp(2,nlegacychxi)=ti
                chxilegacycos(1,nlegacychxi)=ac*1.0d-9
                chxilegacycos(2,nlegacychxi)=tc
                chxilegacycos(3,nlegacychxi)=ar
                chxilegacycos(4,nlegacychxi)=xr
                chxilegacy(nlegacychxi)=0.d0
              endif
c
c     end inclusion ifs
c
            endif
          endif
c
c     end nsequ,nentries
c
        enddo
      enddo
c
c Legacy sec electron rates
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) m,at,io,nsecel
        if (zmap(at).ne.0) then
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            secat(i)=atom
            secio(i)=io
            secel(i)=nsecel
          endif
        endif
      enddo
c
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readcoll2 (luin, error)
c
      include 'cblocks.inc'
c
c Variables
c
      real*8 aa,bb,cc,dd,el,po
      integer*4 i,j,luin,atom,ion,is,at,io
      integer*4 ni,nions,ns,nsequ,sh,shell
c
      character ilgg*4, ibuf(19)*4
      logical error
c
c formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
   30 format(a2,1x,a2,1x,a3,1x,3(i2,1x),g9.3,g9.2,1x,g9.2,g9.2,1x,g9.2)
c
      error=.false.
      atom=0
      ion=0
c
      filename=datadir(1:dtlen)//'data/ionisation/COLLDAT.txt'
      open (luin,file=filename,status='OLD')
   40 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 40
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=*) nsequ
c
c     read in each sequence
c
      do 70 i=1,nsequ
        read (luin,fmt=10) (ibuf(j),j=1,19)
c         write(*, 3) (ibuf(j),j = 1, 11)
        read (luin,*) nions,ns
c
c     for each ion in the sequence
c
        do 60 j=1,nions
          do 50 shell=1,ns
            atom=0
            ion=0
            read (luin,30) el,is,sh,at,io,ni,po,aa,bb,cc,dd
            if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
              atom=zmap(at)
              if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
                ion=io
                collpot(shell,ion,atom)=po
                aar(shell,ion,atom)=aa
                bar(shell,ion,atom)=bb
                car(shell,ion,atom)=cc
                dar(shell,ion,atom)=dd
              endif
            endif
   50     continue
          if ((atom.ne.0).and.(ion.ne.0)) nshells(ion,atom)=ns
          if ((at.eq.11).and.(ion.eq.1)) nshells(ion,atom)=2
   60   continue
   70 continue
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readcont (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 dg20,dgu0,g20,gu0,e2
      real*8 zo
      real*8 tmpgffx(81),tmpgffy(81),tmpgffy2(81)
      real*8 tmpgfbx(41),tmpgfby(41),tmpgfby2(41)
      integer*4 at,atom,i,io,ion,isos,j,ll,pqn
      integer*4 k,luin,nentries,ni,no,zto
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/CONTDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=*) nentries
      do i=1,nentries
        read (luin,*) at,io,isos,no,zto,zo,e2
        if (zmap(at).ne.0) then
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            ion=io
            fn0(ion,atom)=no
            zt0(ion,atom)=zto
            zn0(ion,atom)=zo
            ipot2ev(ion,atom)=e2
            ipot2e(ion,atom)=e2*ev
          endif
        endif
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=*) nklgfb,nklfbe
      read (luin,*) (klpxe(j),j=1,nklfbe)
      do i=1,nklgfb
        do j=1,i
          read (luin,*) pqn,ll,(klgfb(k,j,i),k=1,nklfbe)
        enddo
        do j=i+1,nklgfb
          do k=1,nklfbe
            klgfb(k,j,i)=0.0d0
          enddo
        enddo
      enddo
c
c     precalculate row splines
c
c
      do i=1,nklfbe
        tmpgfbx(i)=klpxe(i)
      enddo
      do pqn=1,nklgfb
        do ll=1,nklgfb
          do i=1,nklfbe
            tmpgfby(i)=klgfb(i,ll,pqn)
            tmpgfby2(i)=0.d0
          enddo
c          call endgrads (tmpgfbx, tmpgfby, nklfbe, dydx0, dydxn)
          call spline (tmpgfbx, tmpgfby, nklfbe, 1.0d30, 1.0d30,
     &     tmpgfby2)
          do i=1,nklfbe
            klgfby2(i,ll,pqn)=tmpgffy2(i)
          enddo
c          write(*,*) 'n l:', pqn, ll-1
c          write(*,'(10(1pg9.2))') (tmpgfbx(i),i=1,nklfbe)
c          write(*,'(10(1pg9.2))') (tmpgfby(i),i=1,nklfbe)
c          write(*,'(10(1pg9.2))') (tmpgfby2(i),i=1,nklfbe)
        enddo
      enddo
c
c Two photon, 0.75*Ipot:
c
      do atom=1,mxelem
        do ion=1,mxion
          ezz(ion,atom)=0.d0
        enddo
      enddo
c
      do atom=1,atypes
        ion=mapz(atom)
        ezz(ion,atom)=ipote(ion,atom)*0.75d0
        if (mapz(atom).gt.1) then
          ion=mapz(atom)-1
          ezz(ion,atom)=ipote(ion,atom)*0.75d0
        endif
      enddo
c
c   Read 2S0 effective recombination data
c
c    Convert to logs for interpolation
c
      read (luin,10) (ibuf(k),k=1,19)
c
c   Case A Hydrogen
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) th42(i),(r2s1a(i,k),k=1,9)
        do k=1,9
          r2s1a(i,k)=dlog10(r2s1a(i,k)+epsilon)
        enddo
      enddo
c
c   Case B Hydrogen
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) th42(i),(r2s1b(i,k),k=1,9)
        do k=1,9
          r2s1b(i,k)=dlog10(r2s1b(i,k)+epsilon)
        enddo
      enddo
      do i=1,ni
        th42(i)=dlog10(th42(i))
      enddo
c
c   Case A Helium II
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) the43(i),(r2s2a(i,k),k=1,9)
        do k=1,9
          r2s2a(i,k)=dlog10(r2s2a(i,k)+epsilon)
        enddo
      enddo
c
c   Case B Helium II
c
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ni
      do i=1,ni
        read (luin,*) the43(i),(r2s2b(i,k),k=1,9)
        do k=1,9
          r2s2b(i,k)=dlog10(r2s2b(i,k)+epsilon)
        enddo
      enddo
      do i=1,ni
        the43(i)=dlog10(the43(i))
      enddo
c
c
c   Read free free gaunt factors
c
      read (luin,10) (ibuf(k),k=1,19)
      write (*,20) (ibuf(j),j=1,19)
c
c   gamma2 and u (kT/IPH), (hnu/kT)
c
      read (luin,*) g20,dg20
      read (luin,*) gu0,dgu0
c
c   gffs & compute bins
c
      read (luin,*) ngffg2,ngffu
c      write(*,'(i4,1x,i4)') ngffg2, ngffu
c
      do i=1,ngffg2
        gffg2(i)=g20+dg20*dble(i-1)
      enddo
c
      do j=1,ngffu
        read (luin,*) (gff(i,j),i=1,ngffg2)
        gffu(j)=gu0+dgu0*dble(j-1)
      enddo
c
c    Convert to logs for interpolation
c
      do j=1,ngffu
        do i=1,ngffg2
          gff(i,j)=dlog10(gff(i,j))
        enddo
      enddo
c
      do i=1,ngffg2
        do j=1,ngffu
          tmpgffx(j)=gffu(j)
          tmpgffy(j)=gff(i,j)
          tmpgffy2(j)=0.d0
          gffy2(i,j)=dlog10(gff(i,j))
        enddo
c        call endgrads (tmpgffx, tmpgffy, ngffu, dydx0, dydxn)
        call spline (tmpgffx, tmpgffy, ngffu, 1.0d30, 1.0d30, tmpgffy2)
        do j=1,ngffu
          gffy2(i,j)=tmpgffy2(j)
        enddo
      enddo
cc
c      read(luin, 1) (ibuf(k),k = 1,19)
c      read(luin,*) ni
c      do  i = 1, ni
c          read(luin,*) m,(gffintspl(i,k),k = 1,5)
c      enddo
cc
      read (luin,10) (ibuf(k),k=1,19)
      read (luin,*) ngffint
      do i=1,ngffint
        read (luin,*) gffintg2(i),gffinty(i)
      enddo
c
c    create spline coeffs once
c
c      call endgrads (gffintg2, gffinty, ngffint, dydx0, dydxn)
      call spline (gffintg2, gffinty, ngffint, 1.0d30, 1.0d30, gffinty2)
c
c      do  i = 1, ngffint
c          write(*,*)gffintg2(i),gffinty(i),gffinty2(i)
c      enddo
c
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine readxdata (luin, error)
cc
c     include 'cblocks.inc'
cc
cc           Variables
cc
c     real*8 br,eg
c     real*8 rntc,xev
c     real*8 xf,xl,xt
c     integer*4 at,atom,i,io,j,luin,nct,nentries,ni
c     character ilgg*4, ibuf(19)*4
c     logical error
cc
cc     formats for reading files
cc
c  10 format(19a4)
c  20 format(' ',19a4)
cc
c     error=.false.
c
c     xlines=0
c     xilines=0
c
c     filename=datadir(1:dtlen)//'data/lines/XLINDAT.txt'
c     open (luin,file=filename,status='OLD')
c
c  30 read (luin,fmt=10) (ibuf(j),j=1,19)
c     ilgg=ibuf(1)
c     if (ilgg(1:1).eq.'%') goto 30
c     write (*,20) (ibuf(j),j=1,19)
c     xlines=0
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) at,io,ni,xl,xev,xf,xt,eg
c       if (zmap(at).ne.0) then
c         atom=zmap(at)
c         if (io.le.maxion(atom)) then
c           xlines=xlines+1
c           xrat(xlines)=atom
c           xion(xlines)=io
c           xiso(xlines)=ni
c           xrlam(xlines)=xl
c           xejk(xlines)=xev
c           xegj(xlines)=xev*eg
c           xfef(xlines)=xf
c           xtrans(xlines)=xt
c         endif
c       endif
c     enddo
c
c     close (luin)
c
c     filename=datadir(1:dtlen)//'data/lines/XINTERDAT.txt'
c     open (luin,file=filename,status='OLD')
cc
c  40 read (luin,fmt=10) (ibuf(j),j=1,19)
c     ilgg=ibuf(1)
c     if (ilgg(1:1).eq.'%') goto 40
c     write (*,20) (ibuf(j),j=1,19)
c     xilines=0
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) at,io,ni,xl,xev,xf,xt
c       if (zmap(at).ne.0) then
c         atom=zmap(at)
c         if (io.le.maxion(atom)) then
c           xilines=xilines+1
c           xiat(xilines)=atom
c           xiion(xilines)=io
c           xiiso(xilines)=ni
c           xilam(xilines)=xl
c           xiejk(xilines)=xev
c           xomeg(xilines)=xf
c           xitr(xilines)=xt
c         endif
c       endif
c     enddo
cc
c     close (luin)
c
c     filename=datadir(1:dtlen)//'data/lines/XHEIFDAT.txt'
c     open (luin,file=filename,status='OLD')
c
c  45 read (luin,fmt=10) (ibuf(j),j=1,19)
c     ilgg=ibuf(1)
c     if (ilgg(1:1).eq.'%') goto 45
c     write (*,20) (ibuf(j),j=1,19)
c     xhelines=0
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) at,io,ni,xl,xev,xf,xt
c       if (zmap(at).ne.0) then
c         atom=zmap(at)
c         if (io.le.maxion(atom)) then
c           xhelines=xhelines+1
c           xheat(xhelines)=atom
c           xheion(xhelines)=io
c           xhelam(xhelines)=xl
c           xhejk(xhelines)=xev
c           xhef(xhelines)=xf
c         endif
c       endif
c     enddo
c     read (luin,fmt=10) (ibuf(j),j=1,19)
c     write (*,20) (ibuf(j),j=1,19)
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) at,xf
c       if (zmap(at).ne.0) then
c         atom=zmap(at)
c         xhebra(atom)=xf
c       endif
c     enddo
c
c     close (luin)
c
c     filename=datadir(1:dtlen)//'data/lines/LMTDAT.txt'
c     open (luin,file=filename,status='OLD')
c  50 read (luin,fmt=10) (ibuf(j),j=1,19)
c     ilgg=ibuf(1)
c     if (ilgg(1:1).eq.'%') goto 50
c     write (*,20) (ibuf(j),j=1,19)
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) ni,rntc,i,br
c       nct=idint(rntc)
c       brr(ni,nct)=1.d0
c       if (dble(nct).eq.rntc) then
c         clkup(ni,nct)=i
c         brr(ni,nct)=br
c       endif
c     enddo
c     read (luin,fmt=10) (ibuf(j),j=1,19)
c     read (luin,*) nentries
c     do j=1,nentries
c       read (luin,*) i,arg(j),brg(j),crg(j),drg(j),erg(j)
c       if (i.ne.j) then
c         write (*,*) 'error reading LMTDAT'
c         write (*,*) i,j,arg(j),brg(j),crg(j),drg(j),erg(j)
c       endif
c     enddo
cc
c     close (luin)
cc
c     return
c     end
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
      subroutine readxr3data (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 luin
      logical error
c
      integer*4 nentries,idx
      character ilgg*4, ibuf(19)*4
c
      integer*4 id,at,io,atom,ion,is
      integer*4 nl,i,j,n
      real*8 ei,gi,aji,wl,gf,tc,br,frac
      character*24 termi
      integer*4 nc,idc,tp,nspl,ids,idcs,ncasc
      integer*4 jj,kk,nt
      integer*4 foundh12s2s12,foundhe12s2s12
c
      integer*4 maxcascade, maxlevels, maxcolls
c
      real*8 btx (mxxr3nspl)
      real*8 bty (mxxr3nspl)
      real*8 bty2(mxxr3nspl)
c local 2d matrix arrays, converted to line list indices
      real*8 xr3eij (mxxr3lvls,mxxr3lvls)  ! cm^-1
      real*8 xr3aji (mxxr3lvls,mxxr3lvls)  ! ^-1
      real*8 xr3gfij(mxxr3lvls,mxxr3lvls)
      real*8 xr3frij(mxxr3lvls,mxxr3lvls)
c
      integer*4 l0,l1,l2
      character*512 s1
      character*512 s2
      character*512 s3
      character*512 s4
c
      integer*4 lenv
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
   40 format('  - Loaded:',i4,' ions, ', i5,' lines.')
c
      error=.false.
      maxcascade=0
      maxlevels=0
      maxcolls=0
c
      filename=datadir(1:dtlen)//'data/lines/XR3DATA.txt'
      open (luin,file=filename,status='OLD')
c
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
      nxr3ions=0
      nxr3lines=0
      read (luin,*) nentries
c
      idx=0
      do n=1,nentries
        read (luin,*) id
        read (luin,*) at,io
c        write (*,*) 'READING XR3DATA ION: ',n,id
        if (id.ne.n) then
          write (*,*) 'ERROR READING XR3DATA ION: ',n,id
          stop
        endif
c
        atom=zmap(at)
        ion=io
        is=at-io+1
c
c Only add lines for species present
c
        if (atom.ne.0) then
          if (ion.le.maxion(atom)) then
c
            idx=idx+1
            xr3at(idx)=atom
            xr3ion(idx)=ion
            nxr3ions=nxr3ions+1
            xr32s12j(idx)=0
c
c get level energies, gs and names
c
            read (luin,*) nl
c        maxlevels=maxlevels+nl
            if (maxlevels.lt.nl) maxlevels=nl
            if (nl.gt.mxxr3lvls) then
              write (*,*) 'ERROR IN  XR3DATA Level Limit Exceeded: '
              write (*,*) idx,nl,mxxr3lvls
              write (*,*) mapz(xr3at(idx)),xr3ion(idx)
              stop
            endif
c
            if (expertmode.gt.0) then
              write (*,*) ' Atom: ',elem(atom),' Ion: ',rom(ion),' nl:',
     &         nl
            endif
c
c clear local matrices
c
            do i=1,nl
              do j=1,nl
                xr3eij(i,j)=0.d0
                xr3aji(i,j)=0.d0
                xr3gfij(i,j)=0.d0
                xr3frij(i,j)=0.d0
              enddo
            enddo
c
            xr3ni(idx)=nl
            foundh12s2s12=0
            foundhe12s2s12=0
            do i=1,nl
              read (luin,*) j,ei,gi,termi
              if (i.ne.j) then
                write (*,*) 'ERROR IN  XR3DATA LVLS: ',idx,i,j
                stop
              endif
              xr3ei(i,idx)=ei
              xr3gi(i,idx)=gi
              xr3invgi(i,idx)=1.d0/gi
              xr3id(i,idx)=termi(1:16)
              if (is.eq.1) then
c
c find H I like 2s 2S1/2 level for 2photon collisions
c first level abound ground with gi=02 termi=2S1|2
                if ((termi(1:5).eq.'2S1|2').and.(gi.eq.2).and.(i.gt.1)
     &           .and.(foundh12s2s12.eq.0)) then
                  xr32s12j(idx)=i
                  foundh12s2s12=1
                  ee2p(atom)=ei*plk*cls
c               write(*,*) at, io, is, idx, xr32S12j(idx),termi(1:5)
                endif
              endif
              if (is.eq.2) then
c
c find He I like  2s 2S1/2 level for 2photon collisions
c
                if ((termi(1:3).eq.'1S0').and.(gi.eq.1).and.(i.gt.1)
     &           .and.(foundhe12s2s12.eq.0)) then
                  xr32s12j(idx)=i
                  foundhe12s2s12=1
                  ee2phe(atom)=ei*plk*cls
c               write(*,*) at, io, is, idx, xr32S12j(idx),termi(1:5)
                endif
              endif
            enddo
            do i=1,nl
              do j=1,nl
                xr3eij(i,j)=dabs(xr3ei(j,idx)-xr3ei(i,idx))
c            xr3egij(i,j,idx)=plk*cls*xr3eij(i,j) ! ergs
              enddo
            enddo
c
c get collisions, and tag lines as they appear in the cascades
c
            read (luin,*) nc
c        maxcolls=maxcolls+nc
            if (maxcolls.lt.nc) maxcolls=nc
            if (nc.gt.mxxr3cols) then
            write (*,*) 'ERROR IN  XR3DATA Collisions Limit Exceeded: '
              write (*,*) idx,nc,mxxr3cols
              write (*,*) mapz(xr3at(idx)),xr3ion(idx)
              stop
            endif
            nxr3ioncol(idx)=nc
            do idc=1,nc
              read (luin,*) id,i,j,tp,nspl,tc,frac
              if (id.ne.idc) then
                write (*,*) 'ERROR IN  XR3DATA COLL: ',idx,id,idc
                stop
              endif
              xr3col_i(idc,idx)=i
              xr3col_j(idc,idx)=j
              xr3col_fr(idc,idx)=frac
              xr3col_tc(idc,idx)=tc
              read (luin,*) (bty(ids),ids=1,nspl)
              do ids=1,nspl
                btx(ids)=(1.d0/(nspl-1))*(ids-1)
                bty2(ids)=0.d0
              enddo
c
              call spline (btx, bty, nspl, 1.d30, 1.d30, bty2)
c
              do ids=1,nspl
                xr3col_x(ids)=btx(ids)
                xr3col_y(ids,idc,idx)=bty(ids)
                xr3col_y2(ids,idc,idx)=bty2(ids)
              enddo
              xr3col_typespl(idc,idx)=tp
              xr3col_nspl(idc,idx)=nspl
              read (luin,*) ncasc
c          maxcascade=maxcascade+ncasc
              if (maxcascade.lt.ncasc) maxcascade=ncasc
              if (ncasc.gt.mxxr3casc) then
                write (*,*) 'ERROR IN  XR3DATA Cascade Limit Exceeded: '
                write (*,*) idx,ncasc,mxxr3casc
                write (*,*) mapz(xr3at(idx)),xr3ion(idx)
                stop
              endif
              xr3col_nl(idc,idx)=ncasc
              do idcs=1,ncasc
                read (luin,*) jj,kk,br,wl,aji,frac
                xr3col_jj(idcs,idc,idx)=jj
                xr3col_kk(idcs,idc,idx)=kk
                xr3col_br(idcs,idc,idx)=br
                xr3aji(kk,jj)=aji
                wl=1.0d8/xr3eij(kk,jj)
                gi=xr3gi(kk,idx)
                gf=wl*wl*aji*1.49920542e-16*gi
                xr3gfij(kk,jj)=gf
                xr3frij(kk,jj)=frac
              enddo
            enddo
c
c matrix-linelist maps:
c
            nt=nxr3lines
            do i=1,(nl-1)
              s1='                        '
              s1=xr3id(i,idx)
              do j=(i+1),nl
                if (xr3gfij(i,j).gt.0.0) then
c
                  nt=nt+1
                  if (nt.gt.mxxr3lines) then
              write (*,*) 'ERROR IN XR3DATA line list Limit Exceeded: '
                    write (*,*) idx,nt,mxxr3lines
                    write (*,*) mapz(xr3at(idx)),xr3ion(idx)
                    stop
                  endif
c
c map idx,i,j to line number
c
                  xr3lines_map(i,j,idx)=nt
                  xr3lines_map(j,i,idx)=nt
c
c inverse map line number to idx,i,j
c
                  xr3lines_imapidx(nt)=idx
                  xr3lines_imap_i(nt)=i
                  xr3lines_imap_j(nt)=j
c
                  xr3lines_at(nt)=atom
                  xr3lines_ion(nt)=ion
                  nialines(ion,atom)=nialines(ion,atom)+1
c
c                  xr3lines_eij(nt)=xr3eij(i,j,idx)
                  xr3lines_egij(nt)=plk*cls*xr3eij(i,j)
c                  xr3lines_evij(nt)=xr3lines_egij(nt)/ev
c
                  xr3lines_lam(nt)=1.0d8/xr3eij(i,j)!vacritz
                  xr3lines_aji(nt)=xr3aji(i,j)
                  xr3lines_gf(nt)=xr3gfij(i,j)
                  xr3lines_frac(nt)=xr3frij(i,j)
c
                  s2='                        '
                  s2=xr3id(j,idx)
                  l1=lenv(s1)
                  l2=lenv(s2)
                  s3=s1(1:l1)//'-'//s2(1:l2)
                  l0=lenv(s3)
                  l0=18-l0
                  s4='                  '
                  if (l0.gt.0) then
                    s3=s4(1:l0)//s3
                  endif
                  xr3lines_tran(nt)=s3(1:18)
                endif
              enddo
            enddo
            nxr3lines=nt
          else
c
c skip ion: level energies, gs and names
c
            read (luin,*) nl
            do i=1,nl
              read (luin,*) j,ei,gi,termi
            enddo
c
c skip collisions
c
            read (luin,*) nc
            do idc=1,nc
              read (luin,*) id,i,j,tp,nspl,tc,frac
              read (luin,*) (bty(ids),ids=1,nspl)
              read (luin,*) ncasc
              do idcs=1,ncasc
                read (luin,*) jj,kk,br,wl,aji,frac
              enddo
            enddo
          endif
        else
c
c skip atom: level energies, gs and names
c
          read (luin,*) nl
          do i=1,nl
            read (luin,*) j,ei,gi,termi
          enddo
c
c skip collisions
c
          read (luin,*) nc
          do idc=1,nc
            read (luin,*) id,i,j,tp,nspl,tc,frac
            read (luin,*) (bty(ids),ids=1,nspl)
            read (luin,*) ncasc
            do idcs=1,ncasc
              read (luin,*) jj,kk,br,wl,aji,frac
            enddo
          enddo
        endif
      enddo
      close (luin)
c
      if (expertmode.gt.0) then
        write (*,*) 'MaxCascade: ',maxcascade
        write (*,*) 'MaxLevels : ',maxlevels
        write (*,*) 'MaxColls  : ',maxcolls
      endif
c
      write (*,40) nxr3ions,nxr3lines
      xr3minxbin=1

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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readxrldata (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 luin
      logical error
c
      integer*4 nentries,idx
      character ilgg*4, ibuf(19)*4
c
      integer*4 id,at,io,atom,ion
      integer*4 nl,i,j,n
      real*8 ei,gi,aji,wl,gf,tc,br,frac
      character*24 termi
      integer*4 nc,idc,tp,nspl,ids,idcs,ncasc
      integer*4 jj,kk,nt
c
      integer*4 maxcascade, maxlevels, maxcolls
c
      real*8 btx (mxxr3nspl)
      real*8 bty (mxxr3nspl)
      real*8 bty2(mxxr3nspl)
c same for large ions
      real*8 xrleij (mxxrllvls,mxxrllvls)  ! cm^-1
      real*8 xrlaji (mxxrllvls,mxxrllvls)  ! ^-1
      real*8 xrlgfij(mxxrllvls,mxxrllvls)
      real*8 xrlfrij(mxxrllvls,mxxrllvls)
c
      integer*4 l0,l1,l2
      character*512 s1
      character*512 s2
      character*512 s3
      character*512 s4
c
      integer*4 lenv
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
   40 format('  - Loaded:',i4,' ions, ', i5,' lines.')
c
      error=.false.
      maxcascade=0
      maxlevels=0
      maxcolls=0
c
      nxrlions=0
      nxrllines=0
c
      filename=datadir(1:dtlen)//'data/lines/XRLDATA.txt'
      open (luin,file=filename,status='OLD')
c
   50 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 50
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,*) nentries
c
      idx=0
      do n=1,nentries
        read (luin,*) id,at,io
c        write (*,*) 'READING XRLDATA ION: ',n,id
        if (id.ne.n) then
          write (*,*) 'ERROR READING XRLDATA ION: ',n,id
          stop
        endif
c
        atom=zmap(at)
        ion=io
c
c Only add lines for species present
c
        if (atom.ne.0) then
          if (ion.le.maxion(atom)) then
c
            idx=idx+1
            xrlat(idx)=atom
            xrlion(idx)=ion
            nxrlions=nxrlions+1
c
c get level energies, gs and names
c
            read (luin,*) nl
c        maxlevels=maxlevels+nl
            if (maxlevels.lt.nl) maxlevels=nl
            if (nl.gt.mxxrllvls) then
              write (*,*) 'ERROR IN  XRLDATA Level Limit Exceeded: '
              write (*,*) idx,nl,mxxrllvls
              write (*,*) mapz(xrlat(idx)),xrlion(idx)
              stop
            endif
            if (expertmode.gt.0) then
              write (*,*) ' Atom: ',elem(atom),' Ion: ',rom(ion),' nl:',
     &         nl
            endif
c
c clear local matrices
c
            do i=1,nl
              do j=1,nl
                xrleij(i,j)=0.d0
                xrlaji(i,j)=0.d0
                xrlgfij(i,j)=0.d0
                xrlfrij(i,j)=0.d0
              enddo
            enddo
c
            xrlni(idx)=nl
            do i=1,nl
              read (luin,*) j,ei,gi,termi
              if (i.ne.j) then
                write (*,*) 'ERROR IN  XRLDATA LVLS: ',idx,i,j
                stop
              endif
              xrlei(i,idx)=ei
              xrlgi(i,idx)=gi
              xrlinvgi(i,idx)=1.d0/gi
              xrlid(i,idx)=termi(1:16)
            enddo
            do i=1,nl
              do j=1,nl
                xrleij(i,j)=dabs(xrlei(j,idx)-xrlei(i,idx))
c            xrlegij(i,j,idx)=plk*cls*xrleij(i,j) ! ergs
              enddo
            enddo
c
c get collisions, and tag lines as they appear in the cascades
c
            read (luin,*) nc
c        maxcolls=maxcolls+nc
            if (maxcolls.lt.nc) maxcolls=nc
            if (nc.gt.mxxrlcols) then
            write (*,*) 'ERROR IN  XRLDATA Collisions Limit Exceeded: '
              write (*,*) idx,nc,mxxrlcols
              write (*,*) mapz(xrlat(idx)),xrlion(idx)
              stop
            endif
            nxrlioncol(idx)=nc
            do idc=1,nc
              read (luin,*) id,i,j,tp,nspl,tc,frac
              if (id.ne.idc) then
                write (*,*) 'ERROR IN  XRLDATA COLL: ',idx,id,idc
                stop
              endif
              xrlcol_i(idc,idx)=i
              xrlcol_j(idc,idx)=j
              xrlcol_fr(idc,idx)=frac
              xrlcol_tc(idc,idx)=tc
              read (luin,*) (bty(ids),ids=1,nspl)
              do ids=1,nspl
                btx(ids)=(1.d0/(nspl-1))*(ids-1)
                bty2(ids)=0.d0
              enddo
c
              call spline (btx, bty, nspl, 1.d30, 1.d30, bty2)
c
              do ids=1,nspl
                xrlcol_x(ids)=btx(ids)
                xrlcol_y(ids,idc,idx)=bty(ids)
                xrlcol_y2(ids,idc,idx)=bty2(ids)
              enddo
              xrlcol_typespl(idc,idx)=tp
              xrlcol_nspl(idc,idx)=nspl
              read (luin,*) ncasc
cc          maxcascade=maxcascade+ncasc
              if (maxcascade.lt.ncasc) maxcascade=ncasc
              if (ncasc.gt.mxxrlcasc) then
                write (*,*) 'ERROR IN  XRLDATA Cascade Limit Exceeded: '
                write (*,*) idx,ncasc,mxxrlcasc
                write (*,*) mapz(xrlat(idx)),xrlion(idx)
                stop
              endif
              xrlcol_nl(idc,idx)=ncasc
              do idcs=1,ncasc
                read (luin,*) jj,kk,br,wl,aji,frac
                xrlcol_jj(idcs,idc,idx)=jj
                xrlcol_kk(idcs,idc,idx)=kk
                xrlcol_br(idcs,idc,idx)=br
                xrlaji(kk,jj)=aji
                wl=1.0d8/xrleij(kk,jj)
                gi=xrlgi(kk,idx)
                gf=wl*wl*aji*1.49920542e-16*gi
                xrlgfij(kk,jj)=gf
                xrlfrij(kk,jj)=frac
              enddo
            enddo
c
c matrix-linelist maps:
c
            nt=nxrllines
            do i=1,(nl-1)
              s1='                        '
              s1=xrlid(i,idx)
              do j=(i+1),nl
                if (xrlgfij(i,j).gt.0.0) then
c
                  nt=nt+1
                  if (nt.gt.mxxrllines) then
              write (*,*) 'ERROR IN XRLDATA line list Limit Exceeded:'
                    write (*,*) idx,nt,mxxrllines
                    write (*,*) mapz(xrlat(idx)),xrlion(idx)
                    stop
                  endif
c
c map idx,i,j to line number
c
                  xrllines_map(i,j,idx)=nt
                  xrllines_map(j,i,idx)=nt
c
c inverse map line number to idx,i,j
c
                  xrllines_imapidx(nt)=idx
                  xrllines_imap_i(nt)=i
                  xrllines_imap_j(nt)=j
c
                  xrllines_at(nt)=atom
                  xrllines_ion(nt)=ion
                  nialines(ion,atom)=nialines(ion,atom)+1
c
c                  xrllines_eij(nt)=xrleij(i,j,idx)
                  xrllines_egij(nt)=plk*cls*xrleij(i,j)
c                  xrllines_evij(nt)=xrllines_egij(nt)/ev
c
                  xrllines_lam(nt)=1.0d8/xrleij(i,j)!vacritz
                  xrllines_aji(nt)=xrlaji(i,j)
                  xrllines_gf(nt)=xrlgfij(i,j)
                  xrllines_frac(nt)=xrlfrij(i,j)
c
                  s2='                        '
                  s2=xrlid(j,idx)
                  l1=lenv(s1)
                  l2=lenv(s2)
                  s3=s1(1:l1)//'-'//s2(1:l2)
                  l0=lenv(s3)
                  l0=18-l0
                  s4='                  '
                  if (l0.gt.0) then
                    s3=s4(1:l0)//s3
                  endif
                  xrllines_tran(nt)=s3(1:18)
                endif
              enddo
            enddo
            nxrllines=nt
          else
c
c skip ion: level energies, gs and names
c
            read (luin,*) nl
            do i=1,nl
              read (luin,*) j,ei,gi,termi
            enddo
c
c skip collisions
c
            read (luin,*) nc
            do idc=1,nc
              read (luin,*) id,i,j,tp,nspl,tc,frac
              read (luin,*) (bty(ids),ids=1,nspl)
              read (luin,*) ncasc
              do idcs=1,ncasc
                read (luin,*) jj,kk,br,wl,aji,frac
              enddo
            enddo
          endif
        else
c
c skip atom: level energies, gs and names
c
          read (luin,*) nl
          do i=1,nl
            read (luin,*) j,ei,gi,termi
          enddo
c
c skip collisions
c
          read (luin,*) nc
          do idc=1,nc
            read (luin,*) id,i,j,tp,nspl,tc,frac
            read (luin,*) (bty(ids),ids=1,nspl)
            read (luin,*) ncasc
            do idcs=1,ncasc
              read (luin,*) jj,kk,br,wl,aji,frac
            enddo
          enddo
        endif
      enddo
      close (luin)
c
      if (expertmode.gt.0) then
        write (*,*) 'MaxCascade: ',maxcascade
        write (*,*) 'MaxLevels : ',maxlevels
        write (*,*) 'MaxColls  : ',maxcolls
      endif
c
      write (*,40) nxrlions,nxrllines
      xrlminxbin=1
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fnii(t4,coef)
      implicit none
      real*8 t4,invt4,logt4
      real*8 coef(8)
      real*8 a,b,c,d,e,f,g,h
      real*8 emiss1,emiss2,emiss3,emiss
      invt4=1.d0/t4
      logt4=dlog10(t4)
      a=coef(1)
      b=coef(2)
      c=coef(3)
      d=coef(4)
      e=coef(5)
      f=coef(6)
      g=coef(7)
      h=coef(8)
      emiss1=a+t4*(b+t4*c)
      emiss2=(d+t4*(e+t4*f))*logt4
      emiss3=(g*logt4*logt4)+h*invt4
      emiss=(emiss1+emiss2+emiss3)
      fnii=emiss
      return
      end
      real*8 function foii(t4,abcd)
      implicit none
      real*8 t4,tm
      real*8 abcd(4)
      real*8 a,b,c,d
      real*8 emiss
      tm=1.d0-t4
      a=abcd(1)
      b=abcd(2)
      c=abcd(3)
      d=abcd(4)
      emiss=(a*(t4**b))*(1.d0+tm*(c+tm*d))
      foii=dlog10(emiss)
      return
      end
      real*8 function fneii(t4,abcdf)
      implicit none
      real*8 t4,tm
      real*8 abcdf(5)
      real*8 a,b,c,d,f
      real*8 emiss
      tm=1.d0-t4
      a=abcdf(1)
      b=abcdf(2)
      c=abcdf(3)
      d=abcdf(4)
      f=abcdf(5)
      emiss=(a*(t4**f))*(1.d0+tm*(b+tm*(c+tm*d)))
      fneii=dlog10(emiss)
      return
      end
      subroutine readheavyrec (luin, error)
c
      include 'cblocks.inc'
c
      integer*4 luin
      logical error
c
      integer*4 nentries,j,idx,ntemps,ndens
      integer*4 nte,ite
      integer*4 id,ida,idb,mult,idm,ncoefs,l0
      real*8 eij,wav,lm,t4
      real*8 aa,ba,ca,da,fa
      real*8 ab,bb,cb,db,fb
      real*8 abcd(4),abcdf(5)
      character*512 trans
      character*512 s3
      character ilgg*4,ibuf(19)*4
      real*8 lgtem(15)
      real*8 coef(15),coef2(15)
      real*8 niicoef(8)
      real*8 splinex(8),spliney(8),spliney2(8)
      real*8 niite(8),niilgte(8)
      real*8 oiite(8),oiilgte(8)
      real*8 neiite(8),neiilgte(8)
c
c Read data for heavy element effective recombination rates
c on a case by case basis.  Each ion treated separately with
c a mix of input file formats depending on each ion's data
c source.
c
c Current ions include:
c    OI : OII -> OI recombination, Escalante and Victor 1992
c    assuming case B for the Triplets intially
c
c           Functions
c
      real*8 fnii,foii,fneii
      integer*4 lenv
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      data niite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
      data oiite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
      data neiite/125.d0,500.d0,1000.d0,5000.d0,1.0d4,1.5d4,2.0d4,2.0d5/
c
      error=.false.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Carbon
c CIII->CII CII recombination
c
      nrccii=0
c
c CII recomb, 16 brightest multiplets resolved into 19 lines
c
      filename=datadir(1:dtlen)//'data/lines/RECOMB_CIIDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,*) nentries
      nrccii=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,eij,wav,trans
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_CIIDAT ENERGIES: ',id,j
          stop
        endif
        lm=1.0d8/eij!vacritz
        rccii_eij(j)=eij*plk*cls
        rccii_lam(j)=lm
        rccii_tr(j)=trans(1:10)
        nialines(2,zmap(6))=nialines(2,zmap(6))+1
      enddo
c
c Temperature grid for recomb coefs, cvt to log
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) ntemps
      read (luin,*) (lgtem(idx),idx=1,ntemps)
      do j=1,ntemps
        rccii_logte(j)=dlog10(lgtem(j))
      enddo
c
c read coeffs,  cvt to log and make splines
c
c Case A first
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) id,(coef(idx),idx=1,ntemps)
        do idx=1,ntemps
          coef(idx)=dlog10(coef(idx))
          rccii_coeffy(1,idx,j)=coef(idx)
        enddo
c
c compute spline coeffs
c
        call spline (rccii_logte, coef, ntemps, 1.0d30, 1.0d30, coef2)
        do idx=1,ntemps
          rccii_coeffy2(1,idx,j)=coef2(idx)
        enddo
c
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_CIIDAT CASEA: ',id,j,ntemps
          stop
        endif
      enddo
c
c Case B
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) id,(coef(idx),idx=1,ntemps)
        do idx=1,ntemps
          coef(idx)=dlog10(coef(idx))
          rccii_coeffy(2,idx,j)=coef(idx)
        enddo
c
c compute spline coeffs
c
        call spline (rccii_logte, coef, ntemps, 1.0d30, 1.0d30, coef2)
        do idx=1,ntemps
          rccii_coeffy2(2,idx,j)=coef2(idx)
        enddo
c
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_CIIDAT CASEB: ',id,j,ntemps
          stop
        endif
      enddo
c end CII
      close (luin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Nitrogen
c NIII->NII NII recombination
c
      nte=8
      do ite=1,nte
        niilgte(ite)=dlog10(niite(ite))-4.d0
      enddo
c
      nrcnii=0
c NII recomb, 7 brightest multiplets resolved into 55 lines
      filename=datadir(1:dtlen)//'data/lines/RECOMB_NIIDAT.txt'
      open (luin,file=filename,status='OLD')
   40 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 40
      write (*,20) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
      read (luin,*) nentries
      nrcnii=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,mult,eij,wav,trans
        l0=lenv(trans)
        l0=16-l0
        if (l0.lt.0) l0=0
        s3='                '
        if (l0.gt.0) then
          trans=s3(1:l0)//trans
        endif
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_NIIDAT: ',id,j,trans
          stop
        endif
        lm=1.0d8/eij!vacritz
        rcnii_eij(j)=eij*plk*cls
        rcnii_lam(j)=lm
        rcnii_tr(j)=trans(1:16)
        nialines(2,zmap(7))=nialines(2,zmap(7))+1
      enddo
c
c fits to recomb coeffs, one set for each density and case A + B
c
      do ite=1,nte
        splinex(ite)=niilgte(ite)
      enddo
c
c Case A first, must find 4 densities
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      if (nentries.ne.4) then
        write (*,*) 'ERROR READING RECOMB_NIIDAT: coeffs A',nentries
        stop
      endif
c for each density
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseA, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin20_y(1,ite,idx)=spliney(ite)
          rcniin20_y2(1,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseA, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin30_y(1,ite,idx)=spliney(ite)
          rcniin30_y2(1,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseA, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin40_y(1,ite,idx)=spliney(ite)
          rcniin40_y2(1,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseA, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin50_y(1,ite,idx)=spliney(ite)
          rcniin50_y2(1,ite,idx)=spliney2(ite)
        enddo
      enddo
c
c Case B second, must find 4 densities
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      if (nentries.ne.4) then
        write (*,*) 'ERROR READING RECOMB_NIIDAT: coeffs B',nentries
        stop
      endif
c for each density
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseB, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin20_y(2,ite,idx)=spliney(ite)
          rcniin20_y2(2,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseB, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin30_y(2,ite,idx)=spliney(ite)
          rcniin30_y2(2,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseB, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin40_y(2,ite,idx)=spliney(ite)
          rcniin40_y2(2,ite,idx)=spliney2(ite)
        enddo
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) (niicoef(ite),ite=1,8)
cc
cc caseB, coef, line,  case=1 A, case=2 B
cc
        do ite=1,nte-1
          t4=niite(ite)*1.0d-4
          spliney(ite)=fnii(t4,niicoef)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcniin50_y(2,ite,idx)=spliney(ite)
          rcniin50_y2(2,ite,idx)=spliney2(ite)
        enddo
      enddo
c
c end NII
c
      close (luin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c OII->OI OI recombination
c
      nrcoi_q=0
      nrcoi_t=0
c
c O I recomb, Triplet and Quintet systems
c
      filename=datadir(1:dtlen)//'data/lines/RECOMB_OIDAT.txt'
      open (luin,file=filename,status='OLD')
   50 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 50
      write (*,20) (ibuf(j),j=1,19)
c
c Quintets first..
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,*) nentries
      nrcoi_q=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,idm,ida,idb,eij,wav,trans
        l0=lenv(trans)
        l0=16-l0
        if (l0.lt.0) l0=0
        s3='                '
        if (l0.gt.0) then
          trans=s3(1:l0)//trans
        endif
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_OIDAT: Q',id,j,trans
          stop
        endif
        lm=1.0d8/eij!vacritz
        rcoi_qeij(j)=eij*plk*cls
        rcoi_qlam(j)=lm
        rcoi_qtr(j)=trans(1:16)
        nialines(1,zmap(8))=nialines(1,zmap(8))+1
      enddo
c
c second order powerlaw fits to recomb coeffs
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) id,aa,ba,ca,ab,bb,cb
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_OIDAT: Q',id,j
          stop
        endif
c
c case, coef, line,  case=1 A, case=2 B
c
        rcoi_qcoef(1,1,j)=aa
        rcoi_qcoef(1,2,j)=ba
        rcoi_qcoef(1,3,j)=ca
        rcoi_qcoef(2,1,j)=ab
        rcoi_qcoef(2,2,j)=bb
        rcoi_qcoef(2,3,j)=cb
      enddo
c
c Triplets...
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,*) nentries
      nrcoi_t=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,idm,ida,idb,eij,wav,trans
        l0=lenv(trans)
        l0=16-l0
        if (l0.lt.0) l0=0
        s3='                '
        if (l0.gt.0) then
          trans=s3(1:l0)//trans
        endif
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_OIDAT: T',id,j,trans
          stop
        endif
        lm=1.0d8/eij!vacritz
        rcoi_teij(j)=eij*plk*cls
        rcoi_tlam(j)=lm
        rcoi_ttr(j)=trans(1:16)
        nialines(1,zmap(8))=nialines(1,zmap(8))+1
      enddo
c
c second order powerlaw fits to recomb coeffs
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) id,aa,ba,ca,ab,bb,cb
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_OIDAT: T',id,j
          stop
        endif
c
c case, coef, line,  case=1 A, case=2 B
c
        rcoi_tcoef(1,1,j)=aa
        rcoi_tcoef(1,2,j)=ba
        rcoi_tcoef(1,3,j)=ca
        rcoi_tcoef(2,1,j)=ab
        rcoi_tcoef(2,2,j)=bb
        rcoi_tcoef(2,3,j)=cb
      enddo
c end OI
      close (luin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Oxygen
c OIII->OII OII recombination
c
      nrcoii=0
      nte=8
      do ite=1,nte
        oiilgte(ite)=dlog10(oiite(ite))-4.d0
      enddo
c
c OII recomb, 35 brightest multiplets resolved into 200 lines
c
      filename=datadir(1:dtlen)//'data/lines/RECOMB_OIIDAT.txt'
      open (luin,file=filename,status='OLD')
   60 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 60
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,*) nentries
      nrcoii=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,mult,eij,wav,trans
        l0=lenv(trans)
        l0=18-l0
        if (l0.lt.0) l0=0
        s3='                  '
        if (l0.gt.0) then
          trans=s3(1:l0)//trans
        endif
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_OIIDAT ENERGIES: ',id,j,
     &     trans
          stop
        endif
        lm=1.0d8/eij!vacritz
        rcoii_eij(j)=eij*plk*cls
        rcoii_lam(j)=lm
        rcoii_mltid(j)=mult
        rcoii_tr(j)=trans(1:18)
        nialines(2,zmap(8))=nialines(2,zmap(8))+1
      enddo
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
c
c read coeffs,  cvt to log and make splines
c
      ncoefs=7
      ndens=4
      do ite=1,nte
        splinex(ite)=oiilgte(ite)
      enddo
c
c Case A first
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) id,(coef(j),j=1,ncoefs)
        abcd(1)=coef(1)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin20_y(1,ite,idx)=spliney(ite)
          rcoiin20_y2(1,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(2)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin40_y(1,ite,idx)=spliney(ite)
          rcoiin40_y2(1,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(3)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin50_y(1,ite,idx)=spliney(ite)
          rcoiin50_y2(1,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(4)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin60_y(1,ite,idx)=spliney(ite)
          rcoiin60_y2(1,ite,idx)=spliney2(ite)
        enddo
c
        if (id.ne.idx) then
          write (*,*) 'ERROR READING RECOMB_OIIDAT CASEA: ',id,idx
          stop
        endif
c
      enddo
c
c Case B
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) id,(coef(j),j=1,ncoefs)
        abcd(1)=coef(1)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin20_y(2,ite,idx)=spliney(ite)
          rcoiin20_y2(2,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(2)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin40_y(2,ite,idx)=spliney(ite)
          rcoiin40_y2(2,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(3)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin50_y(2,ite,idx)=spliney(ite)
          rcoiin50_y2(2,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(4)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin60_y(2,ite,idx)=spliney(ite)
          rcoiin60_y2(2,ite,idx)=spliney2(ite)
        enddo
c
        if (id.ne.idx) then
          write (*,*) 'ERROR READING RECOMB_OIIDAT CASEB: ',id,idx
          stop
        endif
c
      enddo
c
c Case C
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do idx=1,nentries
        read (luin,*) id,(coef(j),j=1,ncoefs)
        abcd(1)=coef(1)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin20_y(3,ite,idx)=spliney(ite)
          rcoiin20_y2(3,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(2)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin40_y(3,ite,idx)=spliney(ite)
          rcoiin40_y2(3,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(3)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin50_y(3,ite,idx)=spliney(ite)
          rcoiin50_y2(3,ite,idx)=spliney2(ite)
        enddo
c
        abcd(1)=coef(4)
        abcd(2)=coef(5)
        abcd(3)=coef(6)
        abcd(4)=coef(7)
        do ite=1,nte-1
          t4=oiite(ite)*1.0d-4
          spliney(ite)=foii(t4,abcd)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcoiin60_y(3,ite,idx)=spliney(ite)
          rcoiin60_y2(3,ite,idx)=spliney2(ite)
        enddo
c
        if (id.ne.idx) then
          write (*,*) 'ERROR READING RECOMB_OIIDAT CASEC: ',id,idx
          stop
        endif
c
      enddo
c end OII
      close (luin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Neon
c NeIII->NeII NeII recombination
c
      nrcneii=0
      nte=8
      do ite=1,nte
        neiilgte(ite)=dlog10(neiite(ite))-4.d0
      enddo
c
c NeII recomb, 6 brightest multiplets resolved into 38 lines
c
      filename=datadir(1:dtlen)//'data/lines/RECOMB_NEIIDAT.txt'
      open (luin,file=filename,status='OLD')
   70 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 70
      write (*,20) (ibuf(j),j=1,19)
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
c
c transition energies and ids
c
      read (luin,*) nentries
      nrcneii=nentries
      do j=1,nentries
        trans='                        '
        read (luin,*) id,mult,eij,wav,trans
        l0=lenv(trans)
        l0=16-l0
        if (l0.lt.0) l0=0
        s3='                '
        if (l0.gt.0) then
          trans=s3(1:l0)//trans
        endif
        if (id.ne.j) then
          write (*,*) 'ERROR READING RECOMB_NEIIDAT: ',id,j,trans
          stop
        endif
        lm=1.0d8/eij!vacritz
        rcneii_eij(j)=eij*plk*cls
        rcneii_lam(j)=lm
        rcneii_tr(j)=trans(1:16)
        nialines(2,zmap(10))=nialines(2,zmap(10))+1
      enddo
c
c third order + powerlaw fits to recomb coeffs
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do ite=1,nte
        splinex(ite)=neiilgte(ite)
      enddo
      do idx=1,nentries
        read (luin,*) id,mult,aa,ba,ca,da,fa,ab,bb,cb,db,fb
        if (id.ne.idx) then
          write (*,*) 'ERROR READING RECOMB_NeIIDAT: ',id,idx
          stop
        endif
c
c case, coef, line,  case=1 A, case=2 B
c
        abcdf(1)=aa
        abcdf(2)=ba
        abcdf(3)=ca
        abcdf(4)=da
        abcdf(5)=fa
        do ite=1,nte-1
          t4=neiite(ite)*1.0d-4
          spliney(ite)=fneii(t4,abcdf)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcneii_y(1,ite,idx)=spliney(ite)
          rcneii_y2(1,ite,idx)=spliney2(ite)
        enddo
        abcdf(1)=ab
        abcdf(2)=bb
        abcdf(3)=cb
        abcdf(4)=db
        abcdf(5)=fb
        do ite=1,nte-1
          t4=neiite(ite)*1.0d-4
          spliney(ite)=fneii(t4,abcdf)
        enddo
        spliney(nte)=spliney(nte-1)-1.d0
        call spline (splinex, spliney, nte, 1.0d30, 1.0d30, spliney2)
        do ite=1,nte
          rcneii_y(2,ite,idx)=spliney(ite)
          rcneii_y2(2,ite,idx)=spliney2(ite)
        enddo
      enddo
c end neII
      close (luin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine read2level (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a21,e12,fab
      real*8 om,w1,w2
      integer*4 at,atom,io,j,luin
      integer*4 ml,nentries,nl
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/lines/RESONDAT2.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
      nlines=0
      read (luin,*) nentries
      do 40 j=1,nentries
        read (luin,*) at,io,e12,om,fab
        if (zmap(at).ne.0) then
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            nlines=nlines+1
            nl=nlines
            ielr(nl)=atom
            ionr(nl)=io
            e12r(nl)=e12
            omres(nl)=om
            fabsr(nl)=fab
            tkexres(nl)=e12/rkb
            rlam(nl)=plk*cls/e12r(nl)!vacritz
            nialines(io,atom)=nialines(io,atom)+1
          endif
        endif
   40 continue
      close (luin)
c
      filename=datadir(1:dtlen)//'data/lines/INTERDAT2.txt'
      open (luin,file=filename,status='OLD')
c
   50 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 50
      write (*,20) (ibuf(j),j=1,19)
      mlines=0
      read (luin,*) nentries
      do 60 j=1,nentries
        read (luin,*) at,io,e12,w1,w2,a21,om
        if (zmap(at).ne.0) then
          atom=zmap(at)
          if (io.le.maxion(atom)) then
            mlines=mlines+1
            ml=mlines
            ielfs(ml)=atom
            ionfs(ml)=io
            e12fs(ml)=e12
            w1fs(ml)=w1
            w2fs(ml)=w2
            a21fs(ml)=a21
            omfs(ml)=om
            tkexfs(ml)=e12/rkb
            fslam(ml)=plk*cls/e12fs(ml)!vacritz
            nialines(io,atom)=nialines(io,atom)+1
          endif
        endif
   60 continue
c
      close (luin)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine read3level (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a21,e12,g1,g2
      integer*4 at,atom,i,io,j,luin,m,nentries,ni
      real*8 p1,p2,tlam,wl
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/lines/THREEDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
c
      nf3ions=0
      nf3trans=3
      write (*,20) (ibuf(j),j=1,19)
c
c   Type 1 fine structure transitions
c
      read (luin,10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      nentries=nentries/3
      do i=1,nentries
        do j=1,3
          read (luin,*) ni,at,io,wl,e12,tlam,a21,g1,p1,g2,p2,m
          if (zmap(at).ne.0) then
            atom=zmap(at)
            if (io.le.maxion(atom)) then
              if (j.eq.1) then
                nf3ions=nf3ions+1
                f3type(nf3ions)=ni
                f3ion(nf3ions)=io
                f3atom(nf3ions)=atom
              endif
              f3lam(j,nf3ions)=tlam!vacritz
              nialines(io,atom)=nialines(io,atom)+1
              wi3(j,nf3ions)=wl
              tkex3(j,nf3ions)=e12
              ei3(j,nf3ions)=e12
              ai3(j,nf3ions)=a21
              gam3(1,j,nf3ions)=g1
              gam3(2,j,nf3ions)=g2
              pow3(1,j,nf3ions)=p1
              pow3(2,j,nf3ions)=p2
            endif
          endif
        enddo
      enddo
c
c   Now type 2 transitions
c
      read (luin,10) (ibuf(j),j=1,19)
      read (luin,*) nentries
      nentries=nentries/3
      do i=1,nentries
        do j=1,3
          read (luin,*) ni,at,io,wl,e12,tlam,a21,g1,p1,g2,p2,m
          if (zmap(at).ne.0) then
            atom=zmap(at)
            if (io.le.maxion(atom)) then
              if (j.eq.1) then
                nf3ions=nf3ions+1
                f3type(nf3ions)=ni
                f3ion(nf3ions)=io
                f3atom(nf3ions)=atom
              endif
              f3lam(j,nf3ions)=tlam!vacritzmicrons
              nialines(io,atom)=nialines(io,atom)+1
              wi3(j,nf3ions)=wl
              ei3(j,nf3ions)=e12
              tkex3(j,nf3ions)=e12
              ai3(j,nf3ions)=a21
              gam3(1,j,nf3ions)=g1
              gam3(2,j,nf3ions)=g2
              pow3(1,j,nf3ions)=p1
              pow3(2,j,nf3ions)=p2
            endif
          endif
        enddo
      enddo
c
      close (luin)
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Read multi-level data, mostly 5 levels, up to 9x9
c Replaces  read5level, read6level, read9level,
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readmultilevel (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 atom,ion,na,nl,ndata,maxlevels
      integer*4 i,id,j,l,la,m,n,ni,nt,nentries,ntrans
      integer*4 luin, l0, l1, l2, omtype, idx
      character ilgg*4, ibuf(19)*4
      character s1*16, s2*16, s3*16
      logical error
      real*8 a, om, td, tc, f1, f2, f3
C     integer*4 itr
C     real*8 t, logt
C     real*8 ups(mxtr)
      integer*4 btn
      real*8 btx (mxnspl)
      real*8 bty (mxnspl)
      real*8 bty2(mxnspl)
c
      real*8 eijpercm(mxnl)
c
c      real*8 fupsilontr
c
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
      filename=datadir(1:dtlen)//'data/lines/MULTIDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
c
      nfmions=0
      maxlevels=0
      do i=1,mxion
        do j=1,mxelem
          fmspindex(i,j)=0
        enddo
      enddo
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do ndata=1,nentries
        read (luin,*) id
        read (luin,*) na,ion
        read (luin,*) nl
c          read (luin,*) nl,omtype,btn
        ntrans=(nl*(nl-1))/2!expectednumber
c
        if (id.ne.ndata) then
          write (*,*) ' INCOMPATIBLE data/lines/MULTIDAT, index:',id,
     &     ndata
          write (*,*) ' Ion:',na,ion
          stop
        endif
c
        atom=zmap(na)
        if (atom.ne.0) then
          if (ion.le.maxion(atom)) then
c
            if (nl.gt.maxlevels) maxlevels=nl
c
            nfmions=nfmions+1
            ni=nfmions
            fmion(ni)=ion
            fmatom(ni)=atom
c            write (*,*) ' MM Ion:',nfmions,mapz(atom),ion,nl
            fmspindex(ion,atom)=ni
            fmnl(ni)=nl
            fmkappaidx(ni)=0
c  find index of kappa correction data if it exists
            do i=1,nkappaions
              if ((kapatom(i).eq.fmatom(ni)).and.(kapion(i).eq.fmion(ni)
     &         )) then
                fmkappaidx(ni)=i
c          write(*,*)  " Atom: ", elem(atom),
c     &                " Ion: ",rom(ion)," Has kappa data idx:", id
              endif
            enddo
c
            do i=1,mxnl
              do j=1,mxnl
                nfmtridx(i,j,ni)=0
              enddo
            enddo
c
c matrix-linelist maps:
c
            nt=0
            do i=1,(nl-1)
              do j=(i+1),nl
                nt=nt+1
                nfmtridx(i,j,ni)=nt
                nfmtridx(j,i,ni)=nt
                nfmlower(nt,ni)=i
                nfmupper(nt,ni)=j
              enddo
            enddo
c
            nfmtrans(ni)=nt
c            write(*,*)  " Atom: ", elem(atom),
c     &                " Ion: ",rom(ion)," Lvls:", nl
c
            do i=1,nl
              read (luin,*) j,eijpercm(i),wim(i,ni),nfmterm(i,ni)
            enddo
c
            do i=1,nl
              do j=1,nl
                eim(i,j,ni)=plk*cls*dabs(eijpercm(i)-eijpercm(j))
              enddo
            enddo
            do i=1,nl
              do j=1,nl
                aim(i,j,ni)=0.d0
                omim(i,j,ni)=0.d0
                tdepm(i,j,ni)=0.d0
              enddo
            enddo
c
c default values
c
            omtype=3
            btn=5
            do i=1,nt
              fmombtc(i,ni)=1.0d4
              fmomdatatype(i,ni)=omtype
              fmombtn(i,ni)=btn
              do m=1,btn
                fmombsplx(m,i,ni)=(1.d0/(btn-1))*(m-1)
                fmombsply(m,i,ni)=0.d0
                fmombsply2(m,i,ni)=0.d0
              enddo
            enddo
c
            read (luin,*) ntrans
c
            if (ntrans.gt.nt) then
              write (*,*) ' INCOMPATIBLE data/lines/MULTIDAT, index:',
     &         ndata
              write (*,*) ' Ion:',na,ion
              write (*,*) ' ntrans> nt:',ntrans,nt
              stop
            endif
c
            do l=1,ntrans
c
              read (luin,*) i,j,omtype,btn,f1,f2,f3
c              if (omtype.eq.0) then
c              read (luin,*) i,j,f1,f2,f3
c              else
c              read (luin,*) i,j,f1,f2
c              f3=0.d0
c              endif
c
              fmombtn(l,ni)=btn
              la=nfmtridx(i,j,ni)
              fmombtn(la,ni)=btn
              fmomdatatype(la,ni)=omtype
c
c  f1 and f2 have different meanings depending on omtype
c  f3 is only used by omtype 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 0: MAPPINGS III const or powerlaw at t4 data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if (omtype.eq.0) then
                a=f1
                om=f2
                td=f3
c               write(*,'(i2,i2,3(x,1pe11.4))')  i, j, a, om, td
c  type 0 , old mappings const or power law om data
c fill sym and non-sym arrays
                aim(j,i,ni)=a
                omim(j,i,ni)=om
                omim(i,j,ni)=om
                tdepm(j,i,ni)=td
                tdepm(i,j,ni)=td
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 1: n-point linear CHIANTI (Type 2) BT92 spline for omega
c Type 2: n-point log10  CHIANTI (Type 6) BT92 spline for omega
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.1).or.(omtype.eq.2)) then
c                read (luin,*) i,j,a,bec
                aife(j,i,ni)=f1
                tc=f2*eim(i,j,ni)/rkb
c convert dc to Tc in K
                fmombtc(l,ni)=tc
                do m=1,btn
                  btx(m)=(1.d0/(btn-1))*(m-1)
                  bty2(m)=0.d0
                enddo
                read (luin,*) (bty(m),m=1,btn)
c    write(*,'(i2,i2,7(x,1pe11.4))')  i, j, a, tc, (btx(m), m = 1, btn)
c
c natural spline with zero endpoint second derivs,
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  fmombsplx(m,l,ni)=btx(m)
                  fmombsply(m,l,ni)=bty(m)
                  fmombsply2(m,l,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Type 3..5: up to 17 point generalised BT92 spline for omega
c precomputed splines, arbitrary scaled x, c in K
c 3: n points, y values only, compute even spaced 0-1 xs, compute y2
c 4: n points, x & y values only, compute y2
c 5: n points, x, y and y2 values precomputed
c 6: as 3, but log10(y)
c 7: as 4, but log10(y)
c 8: as 5, but log10(y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13))
     &         then
c
c n point BT92 spline with TC scaling in K
c
c               read (luin,*) i,j,a,btc
                aim(j,i,ni)=f1
c   in K
                fmombtc(la,ni)=f2
c
                do m=1,btn
                  btx(m)=(1.d0/(btn-1))*(m-1)
                  bty2(m)=0.d0
                enddo
c               read(luin,*) (btx (m), m = 1, btn)
                read (luin,*) (bty(m),m=1,btn)
c               read(luin,*) (bty2(m), m = 1, btn)
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  fmombsplx(m,l,ni)=btx(m)
                  fmombsply(m,l,ni)=bty(m)
                  fmombsply2(m,l,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.4).or.(omtype.eq.7)) then
c
c have C in K, x and y values, allows x to be uneven
c
c               read (luin,*) i,j,a,btc
c
                aim(j,i,ni)=f1
                fmombtc(la,ni)=f2
c
                do m=1,btn
                  bty2(m)=0.d0
                enddo
c
                read (luin,*) (btx(m),m=1,btn)
                read (luin,*) (bty(m),m=1,btn)
c               read(luin,*) (bty2(m), m = 1, btn)
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  fmombsplx(m,l,ni)=btx(m)
                  fmombsply(m,l,ni)=bty(m)
                  fmombsply2(m,l,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.5).or.(omtype.eq.8)) then
c have C in K, x, y and y2 values, allows full control
c
c               read (luin,*) i,j,a,btc
c
                aim(j,i,ni)=f1
                fmombtc(la,ni)=f2
c
                do m=1,btn
                  bty2(m)=0.d0
                enddo
                read (luin,*) (btx(m),m=1,btn)
                read (luin,*) (bty(m),m=1,btn)
                read (luin,*) (bty2(m),m=1,btn)
c
c               call spline(btx,bty,btn,1.d30,1.d30,bty2)
c
                do m=1,btn
                  fmombsplx(m,l,ni)=btx(m)
                  fmombsply(m,l,ni)=bty(m)
                  fmombsply2(m,l,ni)=bty2(m)
                enddo
              endif
c
            enddo
c
            do i=1,nl
              do j=1,nl
                tkexm(j,i,ni)=eim(j,i,ni)/rkb
              enddo
            enddo
            n=0
            do i=1,(nl-1)
              s1=nfmterm(i,ni)
              l1=0
              do idx=1,16
                if (s1(idx:idx).ne.' ') then
                  l1=l1+1
                endif
              enddo
              if (l1.lt.1) l1=1
              do j=(i+1),nl
                n=n+1
                s2=nfmterm(j,ni)
                l2=0
                do idx=1,16
                  if (s2(idx:idx).ne.' ') then
                    l2=l2+1
                  endif
                enddo
                if (l2.lt.1) l2=1
                l0=16-(l1+1+l2)
                if (l0.lt.0) l0=0
                s3='                '
                if (l0.gt.0) then
                  nfmid(n,ni)=s3(1:l0)//s1(1:l1)//'-'//s2(1:l2)
                else
                  nfmid(n,ni)=s1(1:l1)//'-'//s2(1:l2)
                endif
                fmeij(n,ni)=eim(j,i,ni)
                fmlam(n,ni)=plk*cls/eim(j,i,ni)!vacritzcm
              enddo
            enddo
c
          else
c
c            write(*,*) " 1 - Skipping Ion: ", na,ion
c            write(*,*) " 1 - Skipping Transitions: ", ntrans
c
c skip the levels
c
            do l=1,(nl)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            enddo
c
c skip the transitions
c
c               ntrans=(nl*(nl-1))/2
            read (luin,*) ntrans
            do l=1,ntrans
c                 if (omtype.eq.0) then
c                   read (luin,*) i,j,f1,f2,f3
c                 else
c                   read (luin,*) i,j,f1,f2
c                   f3=0.d0
c                 endif
              read (luin,*) i,j,omtype,btn,f1,f2,f3
              if ((omtype.eq.1).or.(omtype.eq.2)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13))
     &         then
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.4).or.(omtype.eq.7)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.5).or.(omtype.eq.8)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
            enddo
c
          endif
        else
c
c            write(*,*) " 2 - Skipping Ion: ", na,ion,omtype
c            write(*,*) " 2 - Skipping Transitions: ", ntrans
c
c skip the levels
c
          do l=1,(nl)
            read (luin,fmt=10) (ibuf(j),j=1,19)
          enddo
c
c skip the colls
c
c               ntrans=(nl*(nl-1))/2
          read (luin,*) ntrans
          do l=1,ntrans
            read (luin,*) i,j,omtype,btn,f1,f2,f3
            if ((omtype.eq.1).or.(omtype.eq.2)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.4).or.(omtype.eq.7)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.5).or.(omtype.eq.8)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
          enddo
c
        endif
      enddo
c
      close (luin)
      nt=0
      do ni=1,nfmions
        ion=fmion(ni)
        atom=fmatom(ni)
        nl=fmnl(ni)
        if (expertmode.gt.0) then
          write (*,*) ' Atom: ',elem(atom),' Ion: ',rom(ion),' nl:',nl
        endif
        do i=1,nl-1
          do j=i+1,nl
            if (aim(j,i,ni).gt.0.0d0) then
              nt=nt+1
              nialines(ion,atom)=nialines(ion,atom)+1
            endif
          enddo
        enddo
      enddo
c
      if (nt.gt.0) then
   40 format(' - Multi-level Data: ',i4,' ions, ',i4,' lines')
        write (*,40) nfmions,nt
c   60 format('                  with  : ',i4,' Max levels')
c      write (*,60) maxlevels
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TEST Upsilons
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      usekappa=.false.
      usekappainterp=.false.
      kappaa=1.d0
      kappab=0.d0
c     kappa    = 2.d0
c     kappaidx = 2
c
c
c      do idx = 1, nfmions
c
c      atom = fmatom(idx)
c      ion  = fmion(idx)
c
c     if ( mapz(atom) .eq. 8) then
c     if ( ion .eq. 1) then
cc
c     nt   = nfmtrans(idx)
cc
c     write(*,*) 'M5 O I Kappa 10 Upsilonij:', nt, 'transitions'
c     write(*,'(17(x,a))') "     LogT  ", (nfmid(itr,idx), itr= 1,nt)
cc
c     do i = 1, 121
c     logt = 1.0d0 + 0.05d0*(i-1)
c     t = 10.0d0**(logt)
c     do itr = 1, nt
c     ups(itr) = fupsilontr(t,itr,idx)
c     enddo
c     write(*,'(17(x,1pe14.6))') logt, (ups(itr), itr= 1,nt)
c     enddo
cc
c     endif
c     endif
c
c     if ( mapz(atom) .eq. 8) then
c     if ( ion .eq. 2) then
cc
c     nt   = nfmtrans(idx)
cc
c     write(*,*) 'M3 O II Kappa 10 Upsilonij:', nt, 'transitions'
c     write(*,'(17(x,a))') "     LogT  ", (nfmid(itr,idx), itr= 1,nt)
cc
c     do i = 1, 121
c     logt = 1.0d0 + 0.05d0*(i-1)
c     t = 10.0d0**(logt)
c     do itr = 1, nt
c     ups(itr) = fupsilontr(t,itr,idx)
c     enddo
c     write(*,'(17(x,1pe14.6))') logt, (ups(itr), itr= 1,nt)
c     enddo
cc
c     endif
c     endif
cc
c       if ( mapz(atom) .eq. 8) then
c       if ( ion .eq. 3) then
cc
c       nt   = nfmtrans(idx)
cc
c       write(*,*) 'MV O III Upsilonij:', nt, 'transitions'
c       if (UseKappa) write(*,*) 'Kappa:', nt
c       write(*,'(17(x,a))') "     LogT  ", (nfmid(itr,idx), itr= 1,nt)
cc
c       do i = 1, 30
c       logt = 3.0d0 + 0.05d0*(i-1)
c       t = 10.d0**logt
c       do itr = 1, nt
c       ups(itr) = fupsilontr(t,itr,idx)
c       enddo
c       write(*,'(17(x,1pe14.6))') logt, (ups(itr), itr= 1,nt)
c       enddo
cc
c       endif
c       endif
cc
c     if ( mapz(atom) .eq. 16) then
c     if ( ion .eq. 2) then
cc
c     nt   = nfmtrans(idx)
cc
c     write(*,*) 'M3 S II Kappa 10 Upsilonij:', nt, 'transitions'
c     write(*,'(17(x,a))') "     LogT  ", (nfmid(itr,idx), itr= 1,nt)
cc
c     do i = 1, 121
c     logt = 1.0d0 + 0.05d0*(i-1)
c     t = 10.0d0**(logt)
c     do itr = 1, nt
c     ups(itr) = fupsilontr(t,itr,idx)
c     enddo
c     write(*,'(17(x,1pe14.6))') logt, (ups(itr), itr= 1,nt)
c     enddo
cc
c     endif
c     endif
cc
c       if ( mapz(atom) .eq. 16) then
c       if ( ion .eq. 3) then
cc
c       nt   = nfmtrans(idx)
cc
c       write(*,*) 'M3 S III Kappa 10 Upsilonij:', nt, 'transitions'
c       write(*,'(17(x,a))') "     LogT  ", (nfmid(itr,idx), itr= 1,nt)
cc
c       do i = 1, 121
c       logt = 1.0d0 + 0.05d0*(i-1)
c       t = 10.0d0**(logt)
c       do itr = 1, nt
c       ups(itr) = fupsilontr(t,itr,idx)
c       enddo
c       write(*,'(17(x,1pe14.6))') logt, (ups(itr), itr= 1,nt)
c       enddo
cc
c       endif
c       endif
c
c       enddo
c
c      stop
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Read multi-level data, for up to 213 level iron III data
c extended to other Fe-like large n ions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readmultife (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 atom,ion,na,nl,ndata,maxlevels
      integer*4 i,j,l,la,m,n,ni,nt,nentries,ntrans
      integer*4 luin, l0, l1, l2, omtype, idx
      character ilgg*4, ibuf(19)*4
      character s1*16, s2*16, s3*16
      logical error
      real*8 a, om, td, tc, f1, f2, f3
c      integer*4 itr
c      real*8 t, logt
c      real*8 ups(mxfetr)
      integer*4 btn
      real*8 btx (mxfenspl)
      real*8 bty (mxfenspl)
      real*8 bty2(mxfenspl)
c
      real*8 eijpercm(mxfenl)
c
c      real*8 ffeupsilontr
c
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
      filename=datadir(1:dtlen)//'data/lines/FEDAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
c
      nfeions=0
      maxlevels=0
      do i=1,mxion
        do j=1,mxelem
          fespindex(i,j)=0
        enddo
      enddo
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do ndata=1,nentries
        read (luin,*) i
        read (luin,*) na,ion
        read (luin,*) nl
c        read (luin,*) nl,omtype,btn
        ntrans=(nl*(nl-1))/2
c
        if (i.ne.ndata) then
          write (*,*) ' INCOMPATIBLE data/lines/FEDAT, index:',i,ndata
          write (*,*) ' Ion:',na,ion
          stop
        endif
c
        atom=zmap(na)
        if (atom.ne.0) then
          if (ion.le.maxion(atom)) then
c
            if (nl.gt.maxlevels) maxlevels=nl
c
            nfeions=nfeions+1
            ni=nfeions
            feion(ni)=ion
            featom(ni)=atom
            fespindex(ion,atom)=ni
            fenl(ni)=nl
c            fekappaidx(ni)=0
c
            do i=1,mxfenl
              do j=1,mxfenl
                nfetridx(i,j,ni)=0
              enddo
            enddo
c
c matrix-linelist maps:
c
            nt=0
            do i=1,(nl-1)
              do j=(i+1),nl
                nt=nt+1
                nfetridx(i,j,ni)=nt
                nfetridx(j,i,ni)=nt
                nfelower(nt,ni)=i
                nfeupper(nt,ni)=j
              enddo
            enddo
c
            nfetrans(ni)=nt
c
            if (nt.ne.ntrans) then
              write (*,*) ' INCOMPATIBLE data/lines/FEDAT, index:',
     &         ndata
              write (*,*) ' Ion:',na,ion
              write (*,*) ' Unexected nt:',nt,((nl-1)*nl)/2
              stop
            endif
c
            do i=1,nl
              read (luin,*) j,eijpercm(i),wife(i,ni),nfeterm(i,ni)
            enddo
c
            do i=1,nl
              do j=1,nl
                eife(i,j,ni)=plk*cls*dabs(eijpercm(i)-eijpercm(j))
              enddo
            enddo
            do i=1,nl
              do j=1,nl
                aife(i,j,ni)=0.d0
                omife(i,j,ni)=0.d0
                tdepfe(i,j,ni)=0.d0
              enddo
            enddo
c
c default values
c
            omtype=3
            btn=11
            do i=1,nt
              feombtc(i,ni)=1.0d4
              feomdatatype(i,ni)=omtype
              feombtn(i,ni)=btn
              do m=1,btn
                feombsplx(m,i,ni)=(1.d0/(btn-1))*(m-1)
                feombsply(m,i,ni)=0.d0
                feombsply2(m,i,ni)=0.d0
              enddo
            enddo
c
c  Diff from normal multi-level data, missing coll str are skipped, so
c  read number of coll present and loop over that.
c  above sets unused colls to 0.0 type 0 with zero splines.
            read (luin,*) ntrans
c             ntrans = nt
c
            if (ntrans.gt.nt) then
              write (*,*) ' INCOMPATIBLE data/lines/FEDAT, index:',
     &         ndata
              write (*,*) ' Ion:',na,ion
              write (*,*) ' ntrans> nt:',ntrans,nt
              stop
            endif
c
            do l=1,ntrans
c
c             if (omtype.eq.0) then
c             read (luin,*) i,j,f1,f2,f3
c             else
c             read (luin,*) i,j,f1,f2
c             f3=0.d0
c             endif
              read (luin,*) i,j,omtype,btn,f1,f2,f3
c
              la=nfetridx(i,j,ni)
              feombtn(la,ni)=btn
              feomdatatype(la,ni)=omtype
c
c  f1 and f2 have different meanings depending on omtype
c  f3 is only used by omtype 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Type 0: MAPPINGS III const or powerlaw at t4 data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if (omtype.eq.0) then
                a=f1
                om=f2
                td=f3
c     get actual index even if is or js are being skipped
c                write(*,'(i2,i2,3(x,1pe11.4))')  i, j, a, om, td
c   type 0 , old mappings const or power law om data
c  fill sym and non-sym arrays
                aife(j,i,ni)=a
                omife(j,i,ni)=om
                omife(i,j,ni)=om
                tdepfe(j,i,ni)=td
                tdepfe(i,j,ni)=td
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 1: n-point linear CHIANTI (Type 2) BT92 spline for omega
c Type 2: n-point log10  CHIANTI (Type 6) BT92 spline for omega
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.1).or.(omtype.eq.2)) then
c               read (luin,*) i,j,a,bec
                aife(j,i,ni)=f1
                tc=f2*eim(i,j,ni)/rkb
c convert dc to Tc in K
                feombtc(la,ni)=tc
                do m=1,btn
                  btx(m)=(1.d0/(btn-1))*(m-1)
                  bty2(m)=0.d0
                enddo
                read (luin,*) (bty(m),m=1,btn)
c
c write(*,'(i2,i2,7(x,1pe11.4))')  i, j, a, tc, (btx(m), m = 1, btn)
c natural spline with zero endpoint second derivs,
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  feombsplx(m,la,ni)=btx(m)
                  feombsply(m,la,ni)=bty(m)
                  feombsply2(m,la,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Type 3..5: up to 17 point generalised BT92 spline for omega
c precomputed splines, arbitrary scaled x, c in K
c 3: n points, y values only, compute even spaced 0-1 xs, compute y2
c 4: n points, x & y values only, compute y2
c 5: n points, x, y and y2 values precomputed
c 6: as 3, but log10(y)
c 7: as 4, but log10(y)
c 8: as 5, but log10(y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13))
     &         then
c
c n point BT92 spline with TC scaling in K
c
c               read (luin,*) i,j,a,btc
                aife(j,i,ni)=f1
c   in K
                feombtc(la,ni)=f2
c
                do m=1,btn
                  btx(m)=(1.d0/(btn-1))*(m-1)
                  bty2(m)=0.d0
                enddo
c               read(luin,*) (btx (m), m = 1, btn)
                read (luin,*) (bty(m),m=1,btn)
c               read(luin,*) (bty2(m), m = 1, btn)
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  feombsplx(m,la,ni)=btx(m)
                  feombsply(m,la,ni)=bty(m)
                  feombsply2(m,la,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.4).or.(omtype.eq.7)) then
c
c have C in K, x and y values, allows x to be uneven
c
c               read (luin,*) i,j,a,btc
c
                aife(j,i,ni)=f1
                feombtc(la,ni)=f2
c
                do m=1,btn
                  bty2(m)=0.d0
                enddo
c
                read (luin,*) (btx(m),m=1,btn)
                read (luin,*) (bty(m),m=1,btn)
c               read(luin,*) (bty2(m), m = 1, btn)
c
                call spline (btx, bty, btn, 1.d30, 1.d30, bty2)
c
                do m=1,btn
                  feombsplx(m,la,ni)=btx(m)
                  feombsply(m,la,ni)=bty(m)
                  feombsply2(m,la,ni)=bty2(m)
                enddo
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if ((omtype.eq.5).or.(omtype.eq.8)) then
c have C in K, x, y and y2 values, allows full control
c
c                read (luin,*) i,j,a,btc
c
                aife(j,i,ni)=f1
                feombtc(la,ni)=f2
c
                do m=1,btn
                  bty2(m)=0.d0
                enddo
                read (luin,*) (btx(m),m=1,btn)
                read (luin,*) (bty(m),m=1,btn)
                read (luin,*) (bty2(m),m=1,btn)
c
c               call spline(btx,bty,btn,1.d30,1.d30,bty2)
c
                do m=1,btn
                  feombsplx(m,la,ni)=btx(m)
                  feombsply(m,la,ni)=bty(m)
                  feombsply2(m,la,ni)=bty2(m)
                enddo
              endif
c
            enddo
c
            do i=1,nl
              do j=1,nl
                tkexfe(j,i,ni)=eife(j,i,ni)/rkb
              enddo
            enddo
            n=0
            do i=1,(nl-1)
              s1=nfeterm(i,ni)
              l1=0
              do idx=1,16
                if (s1(idx:idx).ne.' ') then
                  l1=l1+1
                endif
              enddo
              if (l1.lt.1) l1=1
              do j=(i+1),nl
                n=n+1
                s2=nfeterm(j,ni)
                l2=0
                do idx=1,16
                  if (s2(idx:idx).ne.' ') then
                    l2=l2+1
                  endif
                enddo
                if (l2.lt.1) l2=1
                l0=16-(l1+1+l2)
                if (l0.lt.0) l0=0
                s3='                '
                if (l0.gt.0) then
                  nfeid(n,ni)=s3(1:l0)//s1(1:l1)//'-'//s2(1:l2)
                else
                  nfeid(n,ni)=s1(1:l1)//'-'//s2(1:l2)
                endif
                feeij(n,ni)=eife(j,i,ni)
                felam(n,ni)=plk*cls/eife(j,i,ni)!vacritz
c
              enddo
            enddo
c
          else
c
c            write(*,*) " 1 - Skipping Ion: ", na,ion
c            write(*,*) " 1 - Skipping Transitions: ", ntrans
c
c skip the levels
c
            do l=1,(nl)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            enddo
c
c skip the colls
c
c              ntrans=(nl*(nl-1))/2
            read (luin,*) ntrans
            do l=1,ntrans
c                if (omtype.eq.0) then
c                   read (luin,*) i,j,f1,f2,f3
c                else
c                  read (luin,*) i,j,f1,f2
c                  f3=0.d0
c                endif
              read (luin,*) i,j,omtype,btn,f1,f2,f3
              if ((omtype.eq.1).or.(omtype.eq.2)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13))
     &         then
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.4).or.(omtype.eq.7)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
              if ((omtype.eq.5).or.(omtype.eq.8)) then
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
                read (luin,fmt=10) (ibuf(j),j=1,19)
              endif
            enddo
c
          endif
        else
c
c            write(*,*) " 2 - Skipping Ion: ", na,ion
c            write(*,*) " 2 - Skipping Transitions: ", ntrans
c
c skip the levels
c
          do l=1,(nl)
            read (luin,fmt=10) (ibuf(j),j=1,19)
          enddo
c
c skip the colls
c
c              ntrans=(nl*(nl-1))/2
          read (luin,*) ntrans
          do l=1,ntrans
c                 if (omtype.eq.0) then
c                    read (luin,*) i,j,f1,f2,f3
c                 else
c                   read (luin,*) i,j,f1,f2
c                   f3=0.d0
c                 endif
            read (luin,*) i,j,omtype,btn,f1,f2,f3
            if ((omtype.eq.1).or.(omtype.eq.2)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.3).or.(omtype.eq.6).or.(omtype.eq.13)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.4).or.(omtype.eq.7)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
            if ((omtype.eq.5).or.(omtype.eq.8)) then
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            endif
          enddo
c
        endif
      enddo
c
      close (luin)
c
c count actual number of lines with non-zero Aij
c
      nt=0
      do ni=1,nfeions
        ion=feion(ni)
        atom=featom(ni)
        nl=fenl(ni)
        if (expertmode.gt.0) then
          write (*,*) ' Atom: ',elem(atom),' Ion: ',rom(ion),' nl:',nl
        endif
        do i=1,nl-1
          do j=i+1,nl
            if (aife(j,i,ni).gt.0.0) then
              nt=nt+1
              nialines(ion,atom)=nialines(ion,atom)+1
            endif
          enddo
        enddo
      enddo
      if (nt.gt.0) then
   40 format(' - Large Multi-level Data: ',i2,' ions, ',i4,' lines')
        write (*,40) nfeions,nt
c   60 format('                  Max Levels : ',i4,' Levels')
c      write (*,60) maxlevels
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TEST Fe Upsilons
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
c      UseKappa = .false.
c      useKappaInterp = .false.
c      kappaA   = 1.d0
c      kappaB   = 0.d0
c      kappa    = 10.d0
c      kappaidx = 5
cc
c      do idx = 1, nfeions
cc
c       atom = featom(idx)
c       ion  = feion(idx)
cc
c      if ( mapz(atom) .eq. 26) then
c      if ( ion .eq. 3) then
cc
c      nt   = nfetrans(idx)
c
c      m    = nt/20
c      do j = 1, m+1
c      l = ((j-1)*20)+1
c
c       write(*,*) 'M3 Fe III Upsilonij:', nt, 'transitions'
c       write(*,'(21(x,a14))') "  LogT  ", (nfeid(itr,idx), itr= l,l+19)
cc
c       do i = 1, 41
c       logt = 2.0d0 + 0.1d0*(i-1)
c       t = 10.0d0**(logt)
c       do itr = l, l+19
c       ups(itr) = ffeupsilontr(t,itr,idx)
c       enddo
c       write(*,'(21(x,1pe14.6))') logt, (ups(itr), itr= l,l+19)
c       enddo
c
c      enddo
cc
c      endif
c      endif
cc
c      enddo
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readkappadat (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 atom,ion,na,nl,ndata,k,nk
      integer*4 i,j,l,m,ni,nt,nentries
      integer*4 luin, ktype
      character ilgg*4, ibuf(19)*4
      logical error
      integer*4 btn
      real*8 btc
      real*8 kaps (mxkappas)
      real*8 btx (mxnspl)
      real*8 bty (mxnspl)
      real*8 bty2(mxnspl)
c
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
      filename=datadir(1:dtlen)//'data/lines/KAPPADAT.txt'
      open (luin,file=filename,status='OLD')
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
c
      nkappaions=0
      write (*,20) (ibuf(j),j=1,19)
c     read kappas
      read (luin,*) nk
      nkappas=nk
      read (luin,*) (kaps(m),m=1,nk)
      do k=1,nk
        kappas(k)=kaps(k)
      enddo
c      loop over ions
      nkappaions=0
      read (luin,*) nentries
c      write(*,*) " - Kappa Upsilon Data for : ", nentries, " Species."
      do ndata=1,nentries
        read (luin,*) i
        if (i.ne.ndata) then
          write (*,*) ' INCOMPATIBLE data/lines/KAPPADAT, ion error',i
          error=.true.
          stop
        endif
        read (luin,fmt=*) na,ion
        read (luin,fmt=*) nl
        if (zmap(na).ne.0) then
          atom=zmap(na)
          if (ion.le.maxion(atom)) then
c      write(*,*) "   Kappa Upsilon Data for : ", elem(atom), rom(ion)
            nkappaions=nkappaions+1
            ni=nkappaions
            kapion(ni)=ion
            kapatom(ni)=atom
            kapnl(ni)=nl
            do i=1,nl
              do j=1,nl
                nkaptridx(i,j,ni)=0
              enddo
            enddo
c
c matirix-linelist maps:
c
            nt=0
            do i=1,(nl-1)
              do j=(i+1),nl
                nt=nt+1
                nkaptridx(i,j,ni)=nt
                nkaptridx(j,i,ni)=nt
                nkaplower(nt,ni)=i
                nkapupper(nt,ni)=j
              enddo
            enddo
            nkaptrans(ni)=nt
c
c               write(*,*)  nl," Level Atom: ", elem(atom),
c     &                        " Ion: ",rom(ion),
c     &                        " Transitions: ",nt
c
            do l=1,nt
              read (luin,*) i,j
              read (luin,*) ktype
              kupstype(l,ni)=ktype
              if (ktype.eq.0) then
                read (luin,*) btn,btc
                kupsbtn(l,ni)=btn
                kupsbtc(l,ni)=btc
                read (luin,*) (btx(m),m=1,btn)
                do m=1,btn
                  kupssplx(m,l,ni)=btx(m)
                enddo
                do k=1,nkappas
                  read (luin,*) (bty(m),m=1,btn)
                  read (luin,*) (bty2(m),m=1,btn)
                  do m=1,btn
                    kupssply(m,k,l,ni)=bty(m)
                    kupssply2(m,k,l,ni)=bty2(m)
                  enddo
                enddo
              else
                write (*,*) ' INCOMPATIBLE data/lines/KAPPADAT, fit type
     & error',i
                error=.true.
                stop
              endif
            enddo
c
          else
            nt=(nl*(nl-1))/2
            do i=1,nt*(nkappas*2+4)
              read (luin,fmt=10) (ibuf(j),j=1,19)
            enddo
          endif
        else
          nt=(nl*(nl-1))/2
          do i=1,nt*(nkappas*2+4)
            read (luin,fmt=10) (ibuf(j),j=1,19)
          enddo
        endif
      enddo
c
      close (luin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        and Luc Binette
c       Version 1.0.0q
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readstardat (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 i,j,luin,m
      character ilgg*4, ibuf(20)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
c
      error=.false.
c
      filename=datadir(1:dtlen)//'data/STARDAT.txt'
      open (luin,file=filename,status='OLD')
   20 read (luin,fmt=10) (ibuf(j),j=1,11)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 20
c      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 30 j=1,8
        read (luin,*) m
        read (luin,*) (tmet(j,i),i=1,4)
        read (luin,*) (tmet(j,i),i=5,8)
        read (luin,*) (tmet(j,i),i=9,11)
        if (m.eq.j) goto 30
        write (*,*) 'Error Reading STARDAT.txt'
        error=.true.
        stop
   30 continue
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 40 j=1,16
        read (luin,*) m
        read (luin,*) (tno4(j,i),i=1,4)
        read (luin,*) (tno4(j,i),i=5,8)
        read (luin,*) (tno4(j,i),i=9,11)
        if (m.eq.j) goto 40
        write (*,*) 'Error Reading STARDAT'
        error=.true.
        stop
   40 continue
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      do 50 j=1,16
        read (luin,*) m
        read (luin,*) (tno5(j,i),i=1,4)
        read (luin,*) (tno5(j,i),i=5,8)
        read (luin,*) (tno5(j,i),i=9,11)
        if (m.eq.j) goto 50
        write (*,*) 'Error Reading STARDAT'
        error=.true.
        stop
   50 continue
      close (luin)
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readdust (luin, error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer spmx
      parameter(spmx=1201)
      real*8 x0(spmx),aa(spmx),bb(spmx),cc(spmx),dd(spmx)
      real*8 x, ebinc,dq,dx0,dqdx, pahsigma,pahsum
      real*8 yield(5,5),yesc,y1,y2,ynu
      real*8 englow, engmax, beta,ec
      integer*4 i,ii,inl,j,k,luin,nentries,drad,npahy
      character ilgg*4, ibuf(19)*4
      logical error
c
c     formats for reading files
c
   10 format(19a4)
   20 format(' ',19a4)
c
      error=.false.
c
c     PAH Data
c
c     PAH crossection
c     assumes Coronene C(24)H(12), area=4*pi*a^2 from Draine
c     a=3.7A for Coronene from Draine
c     new uses 10A Draine & Li Grain
c
c
c     scattering set to Zero for moment
c     have dropped 4pi from sigma as it seems to fit.
c
      pahsigma=(100.d0)*1.d-16
c
c   Capacitance
      ec=((eesu*eesu)*pi)/(2*6.6136223d-08)
c
      filename=datadir(1:dtlen)//'data/dust/DUSTDATpah.txt'
      open (luin,file=filename,status='OLD')
c
   30 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 30
      write (*,20) (ibuf(j),j=1,19)
c
      npahy=0
c
c   PAH ionization potentials
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) pahion(j),pahip(j)
      enddo
c
c   PAH Yield
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      npahy=nentries
      do j=1,nentries
        read (luin,*) yield(1,j),yield(2,j)
      enddo
c
c     Negative PAHs
c
      do i=1,2
        do j=1,infph-1
          if (photev(j).gt.pahip(i)) then
            ebinc=cphotev(j)
c          yesc=1.d0
            pahyield(j,i)=1.d0
            do k=1,npahy-1
              if (ebinc.lt.(pahip(i)+yield(1,k+1))) then
                y1=pahip(i)+yield(1,k)
                y2=pahip(i)+yield(1,k+1)
                ynu=(ebinc-y1)/(y2-y1)
                y1=yield(2,k)
                y2=yield(2,k+1)
                pahyield(j,i)=(y1+(ynu*(y2-y1)))
                goto 40
              endif
            enddo
   40       continue
          endif
        enddo
      enddo
c
c     Neutral & positive PAHs
c
      do i=3,4
        englow=-(pahion(i)+1)*ec
        do j=1,infph-1
          if (photev(j).gt.pahip(i)) then
            ebinc=cphotev(j)
            engmax=ebinc-pahip(i)
            beta=1.d0/(engmax-englow)**2
            yesc=1.d0-2.d0*beta*englow**2
            pahyield(j,i)=1.d0
            do k=1,npahy-1
              if (ebinc.lt.(pahip(i)+yield(1,k+1))) then
                y1=pahip(i)+yield(1,k)
                y2=pahip(i)+yield(1,k+1)
                ynu=(ebinc-y1)/(y2-y1)
                y1=yield(2,k)
                y2=yield(2,k+1)
                pahyield(j,i)=(y1+(ynu*(y2-y1)))*yesc
                goto 50
              endif
            enddo
   50       continue
          endif
        enddo
      enddo
c
c  Read PAH emission data
c  Listed as lorentz profiles,with x_0,peak (aa) & sigma (bb)
c  where F_PAH=peak*sigma^2/(sigma^2+(x-x_0)^2)
c  Also uses fixed exponential cut-offs of 1.5eV and 0.01 eV
c   to prevent overflow
c
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) x0(j),aa(j),bb(j)
      enddo
      do inl=1,infph-1
        ebinc=cphotev(j)
        do j=1,nentries
          pahflux(inl)=pahflux(inl)+aa(j)*bb(j)*bb(j)/(bb(j)*bb(j)+
     &     (ebinc-x0(j))**2)
        enddo
        pahflux(inl)=pahflux(inl)*exp(-(ebinc/1.5)**3)*exp(-(0.01/ebinc)
     &   **3)
      enddo
c
c     Normalise energy flux to 1 erg s-1
c
      pahsum=0.d0
      do i=1,infph-1
        ebinc=(photev(i+1)-photev(i))*evplk
        pahsum=pahsum+pahflux(i)*ebinc
      enddo
c
c  convert to photons/s/Hz/Sr
c
      do i=1,infph-1
        ebinc=cphote(j)
        pahflux(i)=pahflux(i)/(fpi*pahsum*ebinc)
      enddo
c
c      open(1011,file='pahflux.sou',status='UNKNOWN')
c      do i=1,infph
c         write(1011,865) photev(i),pahflux(i)
c      enddo
c 865  format(1pg14.7,' ',1pg14.7)
c      close(1011)
c
c  Read data for neutral PAH
c
   60 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 60
      write (*,20) (ibuf(j),j=1,19)
c
c PAHabs
      read (luin,fmt=10) (ibuf(j),j=1,19)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c     Remap to current binning,use logarithmic extrapolation for
c        ebin < min x
c
      dx0=dlog10(x0(2))-dlog10(x0(1))
      dq=dlog10(aa(2))-dlog10(aa(1))
      dqdx=dq/dx0
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahnabs(j)=pahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))*
     &     dqdx+log10(aa(1)))
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahnabs(j)=pahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii))))
          else
            pahnabs(j)=0.d0
          endif
        endif
      enddo
c
c PAHsca
      read (luin,*)
      read (luin,fmt=10) (ibuf(j),j=1,19)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c     Remap to current binning,use logarithmic extrapolation for
c        ebin < min x
c
      dx0=dlog10(x0(2))-dlog10(x0(1))
      dq=dlog10(aa(2))-dlog10(aa(1))
      dqdx=dq/dx0
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahnsca(j)=0.d0
cpahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
c     &                *dqdx + log10(aa(1)))
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahnsca(j)=0.d0
cpahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
c     &                       +x*dd(ii))))
          else
            pahnsca(j)=0.d0
          endif
        endif
      enddo
c
c    PAHcos
      read (luin,*)
      read (luin,fmt=10) (ibuf(j),j=1,19)
c        write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
      enddo
c
c     Remap to current binning
c
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahncos(j)=aa(1)
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahncos(j)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
          else
            pahncos(j)=0.d0
          endif
        endif
      enddo
c
c  Read data for ionised PAH
c
   70 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 70
      write (*,20) (ibuf(j),j=1,19)
c
c PAHabs
      read (luin,fmt=10) (ibuf(j),j=1,19)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c     Remap to current binning,use logarithmic extrapolation for
c        ebin < min x
c
      dx0=dlog10(x0(2))-dlog10(x0(1))
      dq=dlog10(aa(2))-dlog10(aa(1))
      dqdx=dq/dx0
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahiabs(j)=pahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))*
     &     dqdx+log10(aa(1)))
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahiabs(j)=pahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii))))
          else
            pahiabs(j)=0.d0
          endif
        endif
      enddo
c
c PAHsca
      read (luin,*)
      read (luin,fmt=10) (ibuf(j),j=1,19)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do j=1,nentries
        read (luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c     Remap to current binning,use logarithmic extrapolation for
c        ebin < min x
c
      dx0=dlog10(x0(2))-dlog10(x0(1))
      dq=dlog10(aa(2))-dlog10(aa(1))
      dqdx=dq/dx0
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahisca(j)=0.d0
cpahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
c     &                 *dqdx + log10(aa(1)))
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahisca(j)=0.d0
cpahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
c     &                       +x*dd(ii))))
          else
            pahisca(j)=0.d0
          endif
        endif
      enddo
c
c    PAHcos
      read (luin,*)
      read (luin,fmt=10) (ibuf(j),j=1,19)
c        write(*, 3) (ibuf(j),j = 1, 11)
      read (luin,*) nentries
      do i=1,nentries
        read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
      enddo
c
c     Remap to current binning
c
      do j=1,infph-1
        ebinc=cphotev(j)
        if (ebinc.lt.x0(1)) then
          pahicos(j)=aa(1)
        else
          ii=0
          do i=1,nentries-1
            if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
          enddo
          if (ii.gt.0) then
            x=ebinc-x0(ii)
            pahicos(j)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
          else
            pahicos(j)=0.d0
          endif
        endif
      enddo
c
c     get PAH extinction
c
      do i=1,infph-1
        pahnext(i)=pahnabs(i)+pahnsca(i)
        pahiext(i)=pahiabs(i)+pahisca(i)
      enddo
      close (luin)
c
c     Dust grain absorption/extinct curves
c
c
c    Splines for graphite grains type=1
c
      filename=datadir(1:dtlen)//'data/dust/DUSTDATgra.txt'
      open (luin,file=filename,status='OLD')
c
   80 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 80
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) dustbinmax
c      write(*,*) dustbinmax
      do 90 drad=1,dustbinmax
        read (luin,*)
        read (luin,*) grainrad(drad)
c        write(*,*) drad,grainrad(drad)
c
c    Spline for Qabs
c
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning,use logarithmic extrapolation for
c        ebin < min x
c
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dq=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dq/dx0
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            absorp(j,drad,1)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))*dqdx+
     &       log10(aa(1)))
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              absorp(j,drad,1)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              absorp(j,drad,1)=0.d0
            endif
          endif
        enddo
c
c    Spline for Qsca
c
        read (luin,*)
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning
c
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dq=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dq/dx0
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            scatter(j,drad,1)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))*
     &       dqdx+log10(aa(1)))
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              scatter(j,drad,1)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              scatter(j,drad,1)=0.d0
            endif
          endif
        enddo
c
c    Spline for Gcos
c
        read (luin,*)
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning
c
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            gcos(j,drad,1)=aa(1)
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              gcos(j,drad,1)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              gcos(j,drad,1)=0.d0
            endif
          endif
        enddo
c
c     Determine extinction
c
        do j=1,infph-1
          extinct(j,drad,1)=absorp(j,drad,1)+scatter(j,drad,1)
        enddo
   90 continue
      close (luin)
c
c
c    Splines for silicate grains type=2
c
      filename=datadir(1:dtlen)//'data/dust/DUSTDATsil.txt'
      open (luin,file=filename,status='OLD')
c
  100 read (luin,fmt=10) (ibuf(j),j=1,19)
      ilgg=ibuf(1)
      if (ilgg(1:1).eq.'%') goto 100
      write (*,20) (ibuf(j),j=1,19)
      read (luin,fmt=10) (ibuf(j),j=1,19)
      write (*,20) (ibuf(j),j=1,19)
      read (luin,*) dustbinmax
      do 110 drad=1,dustbinmax
        read (luin,*)
        read (luin,*) grainrad(drad)
c
c    Spline for Qabs
c
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning
c
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dq=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dq/dx0
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            absorp(j,drad,2)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))*dqdx+
     &       log10(aa(1)))
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              absorp(j,drad,2)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              absorp(j,drad,2)=0.d0
            endif
          endif
        enddo
c
c    Spline for Qsca
c
        read (luin,*)
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning
c
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dq=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dq/dx0
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            scatter(j,drad,2)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))*
     &       dqdx+log10(aa(1)))
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              scatter(j,drad,2)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              scatter(j,drad,2)=0.d0
            endif
          endif
        enddo
c
c    Spline for Gcos
c
        read (luin,*)
        read (luin,fmt=10) (ibuf(j),j=1,19)
        read (luin,*) nentries
        do i=1,nentries
          read (luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c     Remap to current binning
c
        do j=1,infph-1
          ebinc=cphotev(j)
          if (ebinc.lt.x0(1)) then
            gcos(j,drad,2)=aa(1)
          else
            ii=0
            do i=1,nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii=i
            enddo
            if (ii.gt.0) then
              x=ebinc-x0(ii)
              gcos(j,drad,2)=aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
            else
              gcos(j,drad,2)=0.d0
            endif
          endif
        enddo
c
c     Determine extinction
c
        do j=1,infph-1
          extinct(j,drad,2)=absorp(j,drad,2)+scatter(j,drad,2)
        enddo
  110 continue
c
c     Calculate grain radius bin edges (+-0.025 in dex)
c
      do drad=1,dustbinmax
c        write(*,*) grainrad(drad)
        grainrad(drad)=grainrad(drad)*1.0d-04
        gradedge(drad)=grainrad(drad)*9.440608763d-01
      enddo
      gradedge(dustbinmax+1)=grainrad(dustbinmax)*1.059253725d0
      close (luin)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
