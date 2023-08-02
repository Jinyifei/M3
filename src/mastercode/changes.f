cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***CHANGE ELEMENTAL ABUNDANCES ?
c     FINDS ZGAS RELATIVE TO THE SUN
c     Using Asplund 2010 Solar and not including H or He for Zgas.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine dispabundances (luop, zelem, title)
c
      include 'cblocks.inc'
c
      integer*4 luop
      character*24 title
      real*8 zelem(mxelem)
c
      integer*4 i,j
      real*8 xi(mxelem)
      real*8 zmt,rm
      real*8 mu_neu, mu_ion, mu_hi
      real*8 msum, msum2, nsum, nsum2, zsum
c
      zmt=0.d0
c
      do i=4,atypes
        zmt=zmt+zelem(i)/zsol(i)
      enddo
c
      if (atypes.gt.3) then
        zmt=zmt/(atypes-3)
      else
        zmt=zsol(2)/zelem(2)
      endif
c
      zgas=zmt
c
      rm=1.d0
c
      mu_neu=0.0d0
      mu_ion=0.0d0
      mu_hi=0.0d0
c
      msum=0.0d0
      msum2=0.0d0
      nsum=0.0d0
      nsum2=0.0d0
      zsum=0.0d0
c
      do i=1,atypes
        nsum=nsum+zelem(i)
        nsum2=nsum2+zelem(i)*mapz(i)
        msum=msum+zelem(i)*atwei(i)
        msum2=msum2+zelem(i)*mapz(i)*5.485799e-04
      enddo
c
      mu_neu=msum/nsum
      mu_ion=(msum+msum2)/(nsum2+nsum)
      mu_hi=msum
      zen=nsum
c
      do i=1,atypes
        xi(i)=zelem(i)*atwei(i)/msum
      enddo
c
      write (luop,10) title
   10 format( /a24/
     &         ' ==================================',
     &         '================================')
   20 format(  ' ==================================',
     &         '================================')
   30 format(  ' Linear by number H =  1.000      :')
      write (luop,30)
      write (luop,40) (elem(j),zelem(j),j=1,atypes)
   40 format(4(4x,a2,1pe10.3))
   50 format(  ' Logarithmic by number H =  0.000 :')
      write (luop,50)
      write (luop,60) (elem(j),dlog10(zelem(j)),j=1,atypes)
   60 format(4(4x,a2,0pf8.3,2x))
   70 format(  ' Logarithmic by number H = 12.000 :')
      write (luop,70)
      write (luop,80) (elem(j),dlog10(zelem(j))+12.d0,j=1,atypes)
   80 format(4(4x,a2,f8.3,2x))
      write (luop,20)
      write (luop,90) zgas,zen
      write (luop,100) mu_neu,mu_ion,mu_hi
      write (luop,110) xi(1),xi(2),(1-xi(1)-xi(2))
   90    format('  [ Metallicity (Zgas): ',1pg11.5,
     &             'xSolar,    ni/nh: ',1pg11.5,' ]')
  100    format('  [ mu_neu: ',1pg11.5, '   mu_ion: ',
     &              1pg11.5,'  mu_h: ',1pg11.5,' ]')
  110    format('  [ Mass Fractions  X: ',1pg11.5,
     &           ' Y: ',1pg11.5,' Z: ',1pg11.5,' ]'/)
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
c     ***CHANGE ELEMENTAL ABUNDANCES ?
c     FINDS ZGAS RELATIVE TO THE SUN
c     Using Asplund 2010 Solar and not including H or He for Zgas.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine abecha ()
c
      include 'cblocks.inc'
c
      integer*4 lenv
c
      character*4 ilgg,ibuf(19)
      character abundtitle*24
      character*512 abdir
      character*512 fnam
      character*80 numstring
      integer*4 atom,at,luf,i,j,nentries
      integer*4 l,m,n
      integer*4 numtype
      real*8 zi
      real*8 zmt
c
      logical iexi
c
c
   10 format(19a4)
   20 format(' ',19a4)
c     read local solar abundance file first
c
      iexi=.false.
      fnam='abund/solar.txt'
      m=lenv(fnam)
      inquire (file=fnam(1:m),exist=iexi)
      if (iexi) then
        abdir='abund/'
      else
        abdir='/usr/local/share/mappings/abund/'
      endif
c
      l=lenv(abdir)
      fnam=abdir(1:l)//'solar.txt'
      m=lenv(fnam)
      inquire (file=fnam(1:m),exist=iexi)
      if (iexi) then
c
c     found the file...
c
c     if the file is not present the
c     just put zion0 in zsol
c
        luf=99
        open (unit=luf,file=fnam(1:m),status='OLD')
   40   read (unit=luf,fmt=10) (ibuf(j),j=1,19)
        ilgg=ibuf(1)
        if (ilgg(1:1).eq.'%') goto 40
        numtype=0
        read (luf,fmt='(a80)') numstring
c
c read up to two integers from the single line, using type = 0
c if numtype is *not found*
c
        read (unit=numstring,fmt='(I2,x,I2)') nentries,numtype
        if (numtype.lt.0) numtype=0
c
        if (numtype.eq.0) write (*,*) ' Reading mixed values...'
        if (numtype.eq.1) write (*,*) ' Reading log values...'
        if (numtype.eq.2) write (*,*) ' Reading linear values...'
        if (numtype.eq.3) write (*,*) ' Reading base 12 log values...'
c
        do i=1,nentries
          read (luf,*) at,zi
          if (zmap(at).ne.0) then
            atom=zmap(at)
            if (numtype.eq.0) then
              if (zi.le.0.d0) zi=10.d0**zi
            endif
            if (numtype.eq.1) then
              zi=10.d0**zi
            endif
            if (numtype.eq.3) then
              zi=10.d0**(zi-12.d0)
            endif
            zsol(atom)=zi
          endif
        enddo
        close (unit=luf)
      else
        do i=1,atypes
          zsol(i)=zion0(i)
        enddo
      endif
c
cc
      abnfile='Solar/Asplund 2009 Abundances'
c
      do i=1,atypes
        deltazion(i)=1.0d0
        zion(i)=zion0(i)*deltazion(i)
      enddo
c
      abundtitle=' Total Abundances :'
   50 call dispabundances (6, zion, abundtitle)
c
      write (*,60)
   60 format(/' Change abundances (y/N) : ',$)
c
      read (*,70) ilgg
   70 format(a)
c
      ilgg=ilgg(1:1)
      write (*,*)
c
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
      if (ilgg.ne.'Y') ilgg='N'
c
      if (ilgg.eq.'N') goto 150
c
c     Change abundances....
c
      ilgg='F'
c
c     end of manual input
c
      if (ilgg.eq.'F') then
c
c     Read an abundance file
c
   80   fnam=' '
        write (*,90)
   90    format(/' Enter abundance file name : ',$)
c
        read (*,100) fnam
  100   format(a)
        write (*,*)
c
c look locally and in share
c
        m=lenv(fnam)
        abnfile=fnam(1:m)
c abnfile
        iexi=.false.
        inquire (file=abnfile(1:m),exist=iexi)
        if (iexi.eqv..false.) then
          abdir='abund/'
          l=lenv(abdir)
c look in local abund/
          write (*,*) abnfile(1:m),' NOT FOUND.'
          write (*,*) ' Looking in abund/...'
          abnfile=abdir(1:l)//fnam(1:m)
          n=lenv(abnfile)
          inquire (file=abnfile(1:n),exist=iexi)
        endif
        if (iexi.eqv..false.) then
c look in shared abndir
          write (*,*) abnfile(1:m),' NOT FOUND.'
          abdir='/usr/local/share/mappings/abund/'
          l=lenv(abdir)
          write (*,*) ' Looking in ',abdir(1:l),'...'
          abnfile=abdir(1:l)//fnam(1:m)
          n=lenv(abnfile)
          inquire (file=abnfile(1:n),exist=iexi)
        endif
c
        if (iexi) then
c
c     found the file...
c
          m=lenv(abnfile)
c
          write (*,*) 'FOUND: ',abnfile(1:m)
c
c     FORMAT REQUIRED:
c
c     see example file allen.abn
c
c     otherwise, simply:
c
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max ~80 chars)
c     3) # of entries in the file (integer)  number type
c        (0 or none mixed, 1 all log, 2 all linear)
c     4) two columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Abundance, if mixed -ve: log abundance by number
c                            +ve: mean number abundance
c
          luf=99
          open (unit=luf,file=abnfile(1:m),status='OLD')
  110     read (unit=luf,fmt=10) (ibuf(j),j=1,19)
          ilgg=ibuf(1)
          if (ilgg(1:1).eq.'%') goto 110
          write (*,*) ' Read abundances from :'
          write (*,20) (ibuf(j),j=1,19)
          numtype=0
          read (luf,fmt='(a80)') numstring
c
c read up to two integers from the single line, using type = 0
c is numtype is not found
c
          numtype=-1
          read (unit=numstring,fmt='(I3,I3)') nentries,numtype
          if (numtype.lt.0) numtype=0

          if(nentries.lt.atypes) then
          write(*,*)'   WARNING: abundance file with fewer entries'
          write(*,*)'   than current atom types, some elements will'
          write(*,*)'   not be set as expected.'
          write(*,*)'   Entries:',nentries,' Atoms:',atypes
          write(*,*)'   If entries is unexpectedly small there may be'
          write(*,*)'   a file format error, leading spaces are not'
          write(*,*)'   permitted.'
          endif
c
          if (numtype.eq.0) then
          write (*,*) ' Reading ',nentries,' mixed values...'
          endif
          if (numtype.eq.1) then
          write (*,*) ' Reading ',nentries,' log values...'
          endif
          if (numtype.eq.2) then
          write (*,*) ' Reading ',nentries,' linear values...'
          endif
          if (numtype.eq.3) then
          write (*,*) ' Reading ',nentries,' base 12 log values...'
          endif
c
          do i=1,nentries
            read (luf,*) at,zi
            if (zmap(at).ne.0) then
              atom=zmap(at)
              if (numtype.eq.0) then
                if (zi.le.0.d0) zi=10.d0**zi
              endif
              if (numtype.eq.1) then
                zi=10.d0**zi
              endif
              if (numtype.eq.3) then
                zi=10.d0**(zi-12.d0)
              endif
              zion(atom)=zi
            endif
          enddo
          close (unit=luf)
c
        else
          write (*,120) fnam(1:m)
  120      format(' File: ',a,' NOT FOUND...'/
     &              'Try again? or cancel? (a/c) : ',$)
          read (*,130) ilgg
  130       format(a)
          ilgg=ilgg(1:1)
          write (*,*)
c
          if (ilgg.eq.'a') ilgg='A'
c
          if (ilgg.eq.'A') goto 80
        endif
      endif
c
c     loop back to abundance display
c
      goto 50
c
c     Happy....
c
  150 continue
c
      zmt=0.d0
      zen=0.0d0
      do i=1,atypes
        zion0(i)=zion(i)
        zen=zen+zion(i)
        if (i.gt.3) then
          zmt=zmt+zion(i)/zsol(i)
        endif
      enddo
c
      if (atypes.gt.3) then
        zmt=zmt/(atypes-3)
      else
        zmt=zsol(2)/zion(2)
      endif
c
      zgas=zmt
  160 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      write (*,160)
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
c     ***Apply an optoinal abundance adjustment file, after abecha,
c      before depletion
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deltaabund ()
c
      include 'cblocks.inc'
c
      integer*4 lenv
c
      character*4 ilgg,ibuf(19)
      character*256 abdir
      character*512 fnam
      character*80 numstring
      integer*4 atom,at,luf,i,j,nentries
      integer*4 l,m,n
      integer*4 numtype
      real*8 del
      real*8 zmt,rm, rmf
c
      logical iexi
c
c
   10 format(19a4)
   20 format(' ',19a4)
      do i=1,atypes
        deltazion(i)=1.0d0
        zion(i)=zion0(i)*deltazion(i)
      enddo
c
      deltafile='Default Zero Offsets'
      rmf=1.0d0
      rm=1.d0
c
   30 write (*,40)
   40 format(//'  Abundance Scalings/Offsets are:'/
     & '  ======================================',
     & '============================')
   50 format( '  ===============================',
     &        '===================================')
   60 format(  '  Linear Scalings     :')
      write (*,60)
      write (*,70) (elem(j),deltazion(j),j=1,atypes)
   70 format(4(4x,a2,1pe10.3))
      write (*,80)
   80 format(  '  Logarithmic Offsets :')
      write (*,90) (elem(j),dlog10(deltazion(j)),j=1,atypes)
   90 format(4(4x,a2,0pf8.3,2x))
      write (*,50)
c
      write (*,100)
  100  format(/' Change Abundance Offsets (y/N)? : ',$)
      read (*,110) ilgg
  110 format(a)
      ilgg=ilgg(1:1)
      write (*,*)
c
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
      if (ilgg.ne.'Y') ilgg='N'
c
      if (ilgg.eq.'N') goto 180
c
c     Apply abundance offsets, only from a file....
c
      ilgg='F'
c
      if (ilgg.eq.'F') then
c
        fnam='abund/solar.txt'
        m=lenv(fnam)
        inquire (file=fnam(1:m),exist=iexi)
        if (iexi) then
          abdir='abund/'
        else
          abdir='/usr/local/share/mappings/abund/'
        endif
        l=lenv(abdir)
c
c     Read an abundance offsets file
c
  120   fnam=' '
        write (*,130)
  130    format(' Enter abundance offsets file name : ',$)
c
        read (*,140) fnam
  140   format(a)
        write (*,*)
c
c look locally
c
        iexi=.false.
        m=lenv(fnam)
        deltafile=fnam(1:m)
c deltafile
        inquire (file=deltafile(1:m),exist=iexi)
c
c look locally and in share
c
        m=lenv(fnam)
        abnfile=fnam(1:m)
c abnfile
        iexi=.false.
        inquire (file=abnfile(1:m),exist=iexi)
        if (iexi.eqv..false.) then
          abdir='abund/'
          l=lenv(abdir)
c look in local abund/
          write (*,*) abnfile(1:m),' NOT FOUND.'
          write (*,*) ' Looking in abund/...'
          deltafile=abdir(1:l)//fnam(1:m)
          n=lenv(deltafile)
          inquire (file=deltafile(1:n),exist=iexi)
        endif
        if (iexi.eqv..false.) then
c look in shared abndir
          write (*,*) abnfile(1:m),' NOT FOUND.'
          abdir='/usr/local/share/mappings/abund/'
          l=lenv(abdir)
          write (*,*) ' Looking in ',abdir(1:l),'...'
          abnfile=abdir(1:l)//fnam(1:m)
          n=lenv(deltafile)
          inquire (file=deltafile(1:n),exist=iexi)
        endif
c
        if (iexi) then
c
c     found the file...
c
          n=lenv(deltafile)
c
          write (*,*) 'FOUND: ',deltafile(1:n)
c
c     FORMAT REQUIRED:
c
c     see example file offsets.txt
c
c     otherwise, simply:
c
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max ~80 chars)
c     3) # of entries in the file (integer)  number type
c        (0 or none mixed, 1 all log, 2 all linear)
c     4) two columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Abundance, if mixed -ve: log abundance by number
c                            +ve: mean number abundance
c
          luf=99
          open (unit=luf,file=deltafile(1:n),status='OLD')
  150     read (unit=luf,fmt=10) (ibuf(j),j=1,19)
          ilgg=ibuf(1)
          if (ilgg(1:1).eq.'%') goto 150
          write (*,*) ' Read offsets from :'
          write (*,20) (ibuf(j),j=1,19)
          numtype=0
          read (luf,fmt='(a80)') numstring
c
c read up to two integers from the single line, using type = 0
c if numtype is not found
c
          read (unit=numstring,fmt='(I2,x,I2)') nentries,numtype
          if (numtype.lt.0) numtype=0
          if (numtype.gt.2) numtype=0
c
          if (numtype.eq.0) write (*,*) ' Reading mixed values...'
          if (numtype.eq.1) write (*,*) ' Reading log values...'
          if (numtype.eq.2) write (*,*) ' Reading linear values...'
          if (numtype.eq.3) write (*,*) ' Reading base 12 log values...'
c
          do i=1,nentries
            read (luf,*) at,del
            if (zmap(at).ne.0) then
              atom=zmap(at)
              if (numtype.eq.0) then
                if (del.le.0.d0) del=10.d0**del
              endif
              if (numtype.eq.1) then
                del=10.d0**del
              endif
              if (numtype.eq.3) then
                del=10.d0**(del-12.d0)
              endif
              deltazion(atom)=del
            endif
          enddo
          close (unit=luf)
c
        else
          write (*,160) fnam(1:m)
  160      format(' File: ',a,' NOT FOUND...'/
     &              'Try again? or cancel? (a/c) : ',$)
          read (*,170) ilgg
  170       format(a)
          ilgg=ilgg(1:1)
          write (*,*)
c
          if (ilgg.eq.'a') ilgg='A'
c
          if (ilgg.eq.'A') goto 120
        endif
      endif
c
c     loop back to offsets display
c
      goto 30
c
c     Happy....
c
  180 continue
c
      zmt=0.d0
      zen=0.0d0
      do i=1,atypes
        zion(i)=zion0(i)*deltazion(i)
        zen=zen+zion(i)
        if (i.gt.3) then
          zmt=zmt+zion(i)/zsol(i)
        endif
      enddo
c
      if (atypes.gt.3) then
        zmt=zmt/(atypes-3)
      else
        zmt=zsol(2)/zion(2)
      endif
c
      zgas=zmt
  190 format(
     & /' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::')
      write (*,190)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Set Depletion Pattern
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine depcha ()
c
      include 'cblocks.inc'
c
      real*8 zi,rm,rmf
c
      integer*4 i,j,nentries
      integer*4 atom,at,luf,m,numtype,l,n
c
      logical iexi
c
      character*4 ilgg,ibuf(19)
      character*80 numstring
      character*256 abdir
      character*512 fnam
c
      integer*4 lenv
c
   10 format(19a4)
   20 format(' ',19a4)
c
      depfile='Default Depletions'
      rmf=1.0d0
      rm=1.d0
c
   30 write (*,40)
   40 format(//'  Present dust depletion factors are:'/
     & '  ======================================',
     & '============================')
   50 format( '  ===============================',
     &        '===================================')
   60 format(  '  Linear Depletion     :')
      write (*,60)
      write (*,70) (elem(j),dion0(j),j=1,atypes)
   70 format(4(4x,a2,1pe10.3))
      write (*,80)
   80 format(  '  Logarithmic Depletion :')
      write (*,90) (elem(j),dlog10(dion0(j)),j=1,atypes)
   90 format(4(4x,a2,0pf8.3,2x))
      write (*,50)
      write (*,100)
  100 format(/' Change dust depletion (y/N) : ',$)
      read (*,110) ilgg
  110 format(a)
c
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'y') ilgg='Y'
      if (ilgg.eq.'n') ilgg='N'
c
      if (ilgg.eq.'N') goto 180
      if (ilgg.ne.'Y') goto 180
c
c     Change depletions....
c
      ilgg='F'
c
      if (ilgg.eq.'F') then
c
        fnam='abund/solar.txt'
        m=lenv(fnam)
        inquire (file=fnam(1:m),exist=iexi)
        if (iexi) then
          abdir='abund/'
        else
          abdir='/usr/local/share/mappings/abund/'
        endif
        l=lenv(abdir)
c
c     read depletion file
c
  120   fnam=' '
        write (*,130)
  130   format(/' Enter dust depletion file name : ',$)
c
  140   format(a)
        read (*,140) fnam
        m=lenv(fnam)
        depfile=fnam(1:m)
c deltafile
        inquire (file=depfile(1:m),exist=iexi)
        if (iexi.eqv..false.) then
          abdir='abund/'
          l=lenv(abdir)
c look in local abund/
          write (*,*) abnfile(1:m),' NOT FOUND.'
          write (*,*) ' Looking in ',abdir(1:l),'...'
          depfile=abdir(1:l)//fnam(1:m)
          n=lenv(depfile)
          inquire (file=depfile(1:n),exist=iexi)
        endif
        if (iexi.eqv..false.) then
c look in shared abndir
          write (*,*) abnfile(1:m),' NOT FOUND.'
          abdir='/usr/local/share/mappings/abund/'
          l=lenv(abdir)
          write (*,*) ' Looking in ',abdir(1:l),'...'
          depfile=abdir(1:l)//fnam(1:m)
          n=lenv(depfile)
          inquire (file=depfile(1:n),exist=iexi)
        endif
c
        if (iexi) then
c
c     found the file...
c
c     FORMAT REQUIRED:
c
          n=lenv(depfile)
c
          write (*,*) 'FOUND: ',depfile(1:n)
c
c     otherwise, simply:
c
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max ~80 chars)
c     3) # of entries in the file (integer)  number type (integer)
c     4) two columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Depletion Factor, le 0, -ve means log depletion
c     +ve mean linear scaling factors
c
c
          luf=99
          open (unit=luf,file=depfile(1:n),status='OLD')
  150     read (unit=luf,fmt=10) (ibuf(j),j=1,19)
          ilgg=ibuf(1)
          if (ilgg(1:1).eq.'%') goto 150
          write (*,*)
          write (*,*) ' Read depletion from :'
          write (*,20) (ibuf(j),j=1,19)
          numtype=0
          read (luf,fmt='(a80)') numstring
c
c read up to two integers from the single line, using type = 0
c if numtype is *not found*
c
          read (unit=numstring,fmt='(I2,x,I2)') nentries,numtype
          if (numtype.lt.0) numtype=0
c
          if (numtype.eq.0) write (*,*) ' Reading mixed values...'
          if (numtype.eq.1) write (*,*) ' Reading log values...'
          if (numtype.eq.2) write (*,*) ' Reading linear values...'
          if (numtype.eq.3) write (*,*) ' Reading base 12 log values...'
c
          do i=1,nentries
            read (luf,*) at,zi
            if (zmap(at).ne.0) then
              atom=zmap(at)
              if (numtype.eq.0) then
                if (zi.le.0.d0) zi=10.d0**zi
              endif
              if (numtype.eq.1) then
                zi=10.d0**zi
              endif
              if (numtype.eq.3) then
                zi=10.d0**(zi-12.d0)
              endif
              dion0(atom)=zi
            endif
          enddo
          close (unit=luf)
c
        else
          write (*,160) depfile(1:m)
  160       format(' File: ',a,' NOT FOUND...'/
     &           'Try again? or cancel? (a/c) : ',$)
          read (*,170) ilgg
  170       format(a)
          ilgg=ilgg(1:1)
c
          if (ilgg.eq.'a') ilgg='A'
c
          if (ilgg.eq.'A') goto 120
        endif
      endif
c
c     loop back to display
c
      goto 30
c
c     Happy....
c
  180 continue
c
c     Finally calculate invdion
c
      do j=1,atypes
        dion(j)=dion0(j)
        zion(j)=zion0(j)*deltazion(j)*dion(j)
        invdion(j)=1.d0/dion(j)
      enddo
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to setup charge exchange reactions..
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine chacha ()
c
      include 'cblocks.inc'
c
c
c No one has fiddled with charex for 15 years, so I think it is safe to
c just use the default and continue on.
c Keeping old code in comments incase anyone goes crazy
c RSS2012
c
      return
c c
c c  going crazy...
c c
c cc
c       character ibell*4
c       character ilgg*4, il*4, ili*4, ilt*4, ils*4
c cc
c       real*8 chtes
c       integer*4 j,i1,i2
c cc
c cc  External Functions
c cc
c       real*8 densnum
c cc
c cc    ***CHOICE OF CHARGE EXCHANGE REACTIONS WITH HYDROGEN
c cc
c       ibell(1:4) = char(7)
c 340   format (a)
c       chargemode = 1
c cc
c        write(*,*) ' '
c        write(*,*) ' Standard charge exchange reactions are being used.'
c        write(*,*) ' '
c cc
c 101    format(' Disable the reactions? (y/n) : ',$)
c        write(*,101)
c cc
c       read(*,340) ilgg
c       ilgg = ilgg(1:1)
c       if (ilgg.eq.'y') ilgg = 'Y'
c       if (ilgg.eq.'n') ilgg = 'N'
c       if (ilgg.eq.'Y') chargemode = 2
c cc
c       if (chargemode.eq.1) then
c cc
c 100     format(' Do you wish to modify the reaction list? (y/n) : ',$)
c         write(*,100)
c cc
c       read(*,340) ilgg
c       ilgg = ilgg(1:1)
c       if (ilgg.eq.'y') ilgg = 'Y'
c       if (ilgg.eq.'n') ilgg = 'N'
c       if (ilgg.eq.'Y') then
c         chtes = ((densnum(1.d0)-zion(1))-zion(2))/(zion(1)+zion(2))
c         if (chtes.gt.0.5d0) then
c 503        write(*, 507)
c 507        format(//'THE USE OF CHARGE EXCHANGE REACTIONS IS NOT ',
c      &           'RECOMMENDED'/'WHEN THE ABUNDANCES OF HEAVY ELEMENTS ',
c      &           'ARE TO HIGH'/'DISCARD ALL OF THEM (y/n) : ',$)
c cc
c          read(*, 340) ilgg
c          ilgg = ilgg(1:1)
c          if (ilgg.eq.'y') ilgg = 'Y'
c          if (ilgg.eq.'n') ilgg = 'N'
c          if ((ilgg(1:1) .ne. 'Y').and.(ilgg .ne. 'N')) goto 503
c          if (ilgg(1:1) .eq. 'Y') then
c             do j = 1, 26
c               charco(1,j) = - abs(charco(1,j))
c             enddo
c          endif
c         endif
c cc
c 200        write(*, 210)
c 210        format(//' Charge exchange reactions to be included :')
c            do j = 1, 26
c               i1 = charco(1,j)
c               if (i1 .eq. 0) goto 230
c               ili = 'YES'
c               if (i1.lt.0) ili = 'NO '
c               i1 = iabs(i1)
c               i2 = charco(2,j)
c               ils = '(H )'
c               if (charco(7,j) .eq. 2.0) ils = '(HE)'
c               ilt = '<>'
c               if (charco(3,j).le.0.0) ilt = '< '
c               if (charco(4,j).le.0.0) ilt = ' >'
c               write(*, 220) ils, elem(i1),rom(i2), ilt,
c      &              elem(i1),rom(i2+1), ili
c 220           format(2x,a4,6x,a2,a3,1x,a2,5x,a2,a3,1x,':',2x,a3)
c 230        enddo
c cc
c 235        write(*, 240)
c 240        format(' Do you want to edit the list (Y/N) : ',$)
c cc
c          read(*,340) ilgg
c          ilgg = ilgg(1:1)
c          if (ilgg.eq.'y') ilgg = 'Y'
c          if (ilgg.eq.'n') ilgg = 'N'
c cc
c          if (ilgg .eq. 'N') goto 290
c          if (ilgg .ne. 'Y') goto 235
c cc
c            write(*, 250)
c 250        format(//' Include this reaction (N/CR/Y/E)? : ')
c cc
c          do 270 j = 1, 26
c            i1 = charco(1,j)
c            ili = 'YES'
c            if (i1.lt.0) ili = 'NO '
c            i1 = iabs(i1)
c            if (i1 .eq. 0) goto 270
c            i2 = charco(2,j)
c            ils = '(H )'
c            if (charco(7,j) .eq. 2.0) ils = '(HE)'
c            ilt = '<>'
c            if (charco(3,j).le.0.0) ilt = '< '
c            if (charco(4,j).le.0.0) ilt = ' >'
c            charco(1,j) = i1
c 252          write(*, 254)ils,elem(i1),rom(i2),ilt,
c      &              elem(i1),rom(i2+1),ili
c 254          format(2x,a4,6x,a2,a3,1x,a2,5x,a2,a3,1x,':',
c      &              2x,a3,5x,': ',$)
c c
c            read(*, 340) ilgg
c            ilgg = ilgg(1:1)
c            if (ilgg .eq. 'e') ilgg = 'E'
c            if (ilgg .eq. 'y') ilgg = 'Y'
c            if (ilgg .eq. 'n') ilgg = 'N'
c cc
c            il = ilgg
c            if ((ilgg.eq.' ').or.(ilgg(1:1).eq.'E')) ilgg = ili
c            if (ilgg.eq.'Y') ilgg = 'YES'
c            if (ilgg.eq.'N') ilgg = 'NO '
c            if ((ilgg.ne.'NO ').and.(ilgg.ne.'YES')) goto 269
c            if (ilgg.eq.'NO ') charco(1,j) = - i1
c            if (il(1:1).eq.'E') goto 200
c            goto 270
c 269          write(*, 148) ibell
c 148          format(a4)
c            goto 252
c 270      enddo
c          goto 200
c 290    continue
c       endif
c cc
c cc     end old mappings code
c cc
c       endif
c cc
c       return
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
c
c     General purpose routine to set the ionisation balance in
c     the population array pop (global)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine popcha (model)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 p,fr,epotmi
c
      integer*4 i1,i2,iflag,in,nentries
      integer*4 atom,at,io,luf,i,j,m,luo
c
      character blanc*80, model*32
      character*512 fnam
      character ilgg*4,ibell*4
      character*4 ibuf(19)
c
      logical iexi
c
      integer*4 lenv
c
c    ***CHOICE OF THE ELEMENT
c
      ibell(1:4)=char(7)
      luo=40
c
c
   10 format(i2)
   20 format(a)
c     Default ionisation balance, species with ionisation
c     potetials below 1 Rydberg are singly ionised, those
c     above are neutral (1 part in 10^6 to provide minimal number of
c     electrons in worst case of pure H)
c
      epotmi=iphe
c
      do i=1,atypes
        do j=3,maxion(i)
          pop(j,i)=0.0d0
        enddo
        pop(1,i)=0.990d0
        if (ipote(1,i).lt.epotmi) then
           pop(1,i) = 1.d-2
        endif
        pop(2,i)=1.0d0-pop(1,i)
      enddo
c
      ionsetup='Default ionisation'
c
   30 format(//
     &  ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Setting the ionisation state for ',a16/
     &  ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     &' Default: Elements with ionisation potentials'/
     &' above 1 Rydberg are neutral.  All others are'/
     &' singly ionised.'//)
c
      write (*,30) model
c
      call copypop (pop, pop0)
c
c     Allow an ionisation equilibrium
c     calculation to be a preionisation..
c
   40 write (*,50)
   50 format(//
     &  ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '  Choose an ionisation balance:'/
     &  ' ::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::',/
     & '    C : Calulate a preionisation balance'/
     & '    D : Default ionisation values'//
     & '    E : Enter ionisation values'/
     & '    F : read a ionisation balance from a File'/
     & '    S : Save current balance to a File'//
     & '    X : eXit with current balance'/
     & ' :: ',$)
c
      read (*,20) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'c') ilgg='C'
      if (ilgg.eq.'d') ilgg='D'
      if (ilgg.eq.'e') ilgg='E'
      if (ilgg.eq.'f') ilgg='F'
      if (ilgg.eq.'s') ilgg='S'
      if (ilgg.eq.'x') ilgg='X'
c
      if (ilgg.eq.'D') goto 400
c
      if (ilgg.eq.'E') goto 140
c
      if (ilgg.eq.'S') then
        luf=99
   51   fnam=' '
        write (*,52)
   52   format(/' File name to save: ',$)
        read (*,53,err=51) fnam
   53   format(a)
c
         m=lenv(fnam)
         open (luf,file=fnam(1:m),status='NEW')
         write (luf,'("%")')
         write (luf,'("% Ionisation Balance File")')
         write (luf,'("%")')
         write (luf,'("%")')
         write (luf,*) 'Ionisation balance by :MAPPINGS V ', theversion
         nentries=0
c
         do i=1,atypes
           do j=1,maxion(i)
             nentries=nentries+1
           enddo
         enddo
c
         write (luf,*) nentries
c
         do i=1,atypes
           do j=1,maxion(i)
   54        format(i2,1x,i2,1x,1pg14.6)
             write (luf,54) mapz(i),j,pop(j,i)
           enddo
         enddo
         close(luf)

        goto 40
      endif
c
      if (ilgg.eq.'F') then
c
c     Read balance file
c
   60   fnam=' '
        write (*,70)
   70    format(/' Enter file name : ',$)
        read (*,80,err=60) fnam
   80    format(a)
        m=lenv(fnam)
        inquire (file=fnam(1:m),exist=iexi)
        if (iexi.eqv..false.) then
          write (*,*) fnam(1:m),' NOT FOUND.'
          write (*,*) ' Looking in ',datadir(1:dtlen)
          fnam=datadir(1:dtlen)//fnam(1:m)
          m=lenv(fnam)
          inquire (file=fnam(1:m),exist=iexi)
        endif
c
c
        if (iexi) then
c
c     Found the file...
c
c     FORMAT REQUIRED:
c
c
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max 44 chars)
c     3) # of entries in the file (integer)
c     4) three columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Ion (1 is neutral etc)
c     c) Fraction (+ve numbers are fractions, -ve are logs)
c
c
c     first zero arrays
c
          do i=1,atypes
            do j=1,maxion(i)
              pop(j,i)=0.0d0
            enddo
          enddo
c
   90       format(19a4)
  100       format(' ',19a4)
c
          luf=99
          open (unit=luf,file=fnam(1:m),status='OLD')
  110     read (unit=luf,fmt=90) (ibuf(j),j=1,19)
          ilgg=ibuf(1)
          if (ilgg(1:1).eq.'%') goto 110
          write (*,*) ' Read ion balance from :'
          write (*,100) (ibuf(j),j=1,19)
          ionsetup=fnam(1:m)
          read (luf,fmt=*) nentries
          do i=1,nentries
            read (luf,*) at,io,fr
            if (zmap(at).ne.0) then
              atom=zmap(at)
              if (fr.lt.0.d0) fr=10**fr
              pop(io,atom)=fr
            endif
          enddo
          close (unit=luf)
          call copypop (pop, pop0)
c
        else
          write (*,120) fnam(1:m)
  120       format(' File: ',a,' NOT FOUND...'/
     &           'Try again? or cancel? (a/c) : ',$)
          read (*,130) ilgg
  130       format(a)
          ilgg=ilgg(1:1)
c
          if (ilgg.eq.'a') ilgg='A'
c
          if (ilgg.eq.'A') goto 60
        endif
c
        goto 40
      endif
c
      if (ilgg.eq.'C') then
        call sinsla (model)
        call copypop (pop, pop0)
        ionsetup='Calculated at start'
        goto 40
      endif
c
      if (ilgg.eq.'X') goto 410
c
      goto 40
c
c
c     Original ionisation balance entering code
c
c
  140 write (*,150) (i,elem(i),i=1,atypes)
  150 format(6(i2,1x,a2,5x))
c
      ionsetup='Entered manually'
c
  160 write (*,170)
c
  170 format(//' Give element number (0=end, 99=all) : ',$)
      read (*,10) in
      if (in.le.0) goto 410
      if (in.gt.99) goto 160
c
      i1=in
      i2=in
c
      if (in.eq.99) then
        i1=1
        i2=atypes
      endif
c
      do 390 i=i1,i2
        iflag=0
c
  180 format(//' Species versus ionisation fraction for element ',
     &' : ',a2,'  @@@@@@@')
  190   write (*,180) elem(i)
c
        do j=1,maxion(i)
  200     format(' ',a2,a6,1x,':',3x,1pg13.6)
          write (*,200) elem(i),rom(j),pop(j,i)
        enddo
c
        if (iflag.eq.0) goto 240
c
  210   write (*,220)
  220 format(' ::::::::::    Alter?(y/n) : ',$)
        read (*,230,err=210) ilgg
  230 format(a)
        ilgg=ilgg(1:1)
        if (ilgg.eq.'n') ilgg='N'
        if (ilgg.eq.'y') ilgg='Y'
c
        if (ilgg.eq.'N') goto 390
        if (ilgg.ne.'Y') goto 210
c
c
c    ***ALTERING OF INITIAL VALUES
c
  240   write (*,250)
  250 format(' Enter new values (number <= 0.0 as log, number 0-1,',
     &  ' C=conserve, X=eXit)')
c
        iflag=1
  260   do 290 j=1,maxion(i)
c
  270 format(' ',a2,a6,1x,':',3x,1pg13.6,2x,': ',$)
          write (*,270) elem(i),rom(j),pop(j,i)
c
  280 format(a)
          read (*,280) blanc
          ilgg=blanc(1:1)
          if (ilgg.eq.'x') ilgg='X'
          if (ilgg.eq.'c') ilgg='C'
c
          if (ilgg.eq.'C') goto 290
          if (ilgg.eq.'X') goto 300
c
          read (blanc,*,err=370) p
c
          if (p.lt.0.d0) p=10.d0**p
          if (p.gt.1.d0) p=1.0d0
c
          pop(j,i)=p
c
  290   continue
c
c    ***CHECK ON NORMALISATION
c
  300   fr=0.d0
        do 310 j=1,maxion(i)
  310     fr=fr+pop(j,i)
        if (dabs(fr-1.d0).gt.0.1d0) goto 330
        do 320 j=1,maxion(i)
  320     pop(j,i)=pop(j,i)/fr
        goto 350
  330   write (*,340) fr,ibell
  340  format(' Normalisation error.',
     & ' ERRNO :',1pg9.2,a1)
        goto 240
  350   write (*,360) elem(i)
  360  format(/' New normalised values for element #',' : ',a2)
c
c
        goto 190
  370   write (*,380) ibell
  380 format(' ERROR *******************************',a1)
c
        goto 260
c
  390 continue
c
      call copypop (pop, pop0)
c
      goto 40
c
c     Default ionisation balance, species with ionisation
c     potetials below 1 Rydberg are singly ionised, those
c     above are neutral
c
  400 epotmi=iphe
c
      do i=1,atypes
        do j=3,maxion(i)
          pop(j,i)=0.0d0
        enddo
        pop(1,i)=0.990d0
        if (ipote(1,i).lt.epotmi) then
           pop(1,i) = 1.d-2
        endif
        pop(2,i)=1.0d0-pop(1,i)
      enddo
c
      ionsetup='Default ionisation'
c
      call copypop (pop, pop0)
c
      goto 410
c
c
  410 call zeroemiss ! clear any residual radiation
c
      return
c
      end

      subroutine ciepops( t, dh)
c
      include 'cblocks.inc'
c
c parameter, temperature
      real*8 t, dh, de
c
c locals, dh = 1.0 de = 1.0
c
      integer*4 i,j
      real*8 tc,dhc, dec
      real*8 cab,dr,dv,epotmi,fi
      real*8 rad,trea,trec
      real*8 popz(mxion, mxelem),dift
c
      character imod*4,lmod*4
c
c           Functions
c
      real*8 frectim
      real*8 feldens
c
c local copies
c
      t=dabs(t)
      tc=t
      dhc=dh
      dec=dhc*zion(1)
c
      cab=0.0d0
      epotmi=iphe
c
      trea=epsilon
c
      call copypop (pop, pop0)
c
c setup all singly ionised
c
      do i=1,atypes
        do j=3,maxion(i)
          pop(j,i)=0.0d0
        enddo
        pop(1,i)=0.d0
        pop(2,i)=1.0d0
      enddo
c
c zeros t, restore after
c
      call zer
c
      t=tc
      de=dec
      dh=dhc
c
      jspot='NO'
      jcon='YES'
c
      lmod='SO'
      imod='ALL'
      rad=1.d38
      dr=1.0d0
      dv=0.0d0
      fi=1.d0
      wdil=0.0d0
c
      caseab(1)=cab
      caseab(2)=cab
c
      call localem (t, de, dh)
      call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff (fi, dh)
      call equion (t, de, dh)
c
        i=0
c
        trea=1.d-12
        dift=0.d0
c
   10   call copypop (pop, popz)
c
        call localem (t, de, dh)
        call totphot (t, dh, fi, rad, dr, dv, wdil, lmod)
        call zetaeff (fi, dh)
        call equion (t, de, dh)
c
        call difpop (pop, popz, trea, atypes, dift)
c
        call copypop (pop, popz)
c
        i=i+1
        de=feldens(dh,pop)

c
        if ((dift.ge.1.0d-4).or.(i.le.4)) goto 10
c
c         write (*,*) i,dift
      trec=frectim(t,de,dh)
c
      call copypop (pop, pop0)
c
      call zeroemiss ! clear any residual radiation
      return
c
      end

