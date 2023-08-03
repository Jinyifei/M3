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
c     Subject a single ion to temperature and density ranges
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ionemit
c
c
      include 'cblocks.inc'
c
      integer*4 luop
      integer*4 i,j,k,idx,kmax
      integer*4 i0,j0,i1,j1,itr0,itr1
      integer*4 nl,line,series
      integer*4 at,io
      integer*4 nt, itr, sl
      integer*4 ntsteps
      integer*4 ndhsteps
      integer*4 modeltype
      integer*4 linetype
c      character str*18
      character fn*32
      character fl*32
      character pfx*16
      character sfx*4
      character ilgg*4
      character imod*4
      character model*32,s*32

c
      real*8 t, de, dh, dr
      real*8 tmin,tdelta,logt
      real*8 dhmin,dhdelta
      real*8 cab
      real*8 br0,br1
      real*8 br(mxfetr)
      real*8 ups(mxfetr)
      real*8 emiss(200,101)
c
      real*8 fupsilontr
      real*8 ffeupsilontr

      character*32 trimid(mxtr)
c
   10 format(//' Single Multi-Level Ion Emission Model'/
     &' produced by MAPPINGS V ',a8,'  Run:',a80/)
      model='Single ion model'
      luop=23
      imod='ALL'
c
      fn=' '
      pfx='ion'
      sfx='csv'
      call newfile (pfx, 3, sfx, 3, fn)
      fl=fn(1:13)
      open (luop,file=fl,status='NEW')
c
      call zer
      t=1.0d4
      de=1.0d0
c
      itr0=1
      itr1=2
      at=zmap(1)
      io=1
c
   20 format(a)
   30 format(//' Choose a single multi-level ion model type :'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     A :  Fixed temperature, vary density'/
     &'     B :  Fixed density, vary temperature'/
     &'  :: ',$ )
   40 write (*,30)
      read (*,20) ilgg
      ilgg=ilgg(1:1)
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if ((ilgg.lt.'A').or.(ilgg.gt.'B')) goto 40
      if (ilgg.eq.'A') modeltype=0
      if (ilgg.eq.'B') modeltype=1
c
   50 format(//' Choose Collisional Excitation or Recombination :'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     A :  Collisional, Multi-level Ions '/
     &'     B :  Collisional, Massively Multi-level Iron Ions '/
     &'     C :  Recombination Spectra'/
     &'  :: ',$ )
   60 write (*,50)
      read (*,20) ilgg
      ilgg=ilgg(1:1)
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
      if ((ilgg.lt.'A').or.(ilgg.gt.'C')) goto 60
      if (ilgg.eq.'A') linetype=0
      if (ilgg.eq.'B') linetype=2
      if (ilgg.eq.'C') linetype=1
c
      if (linetype.eq.0) then
   70 format(//' Choose a multi-level species :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,70)
        do i=1,(nfmions/5)+1
          j=(i-1)*5
          kmax=min(j+5,nfmions)-j
          write (*,'(5(2x,i3,": ",a2,a6))') (j+k,elem(fmatom(j+k)),
     &     rom(fmion(j+k)),k=1,kmax)
        enddo
   80 format(' :: ',$)
        write (*,80)
c
        read (*,*) idx
        if (idx.lt.1) idx=1
        if (idx.gt.nfmions) idx=nfmions
        at=fmatom(idx)
        io=fmion(idx)
c
        write (*,*)
        write (*,*)
c
        do i=1,atypes
          do j=1,maxion(i)
            pop(j,i)=0.0d0
          enddo
        enddo
c
        pop(io,at)=1.000d0
c
        write (luop,*) 'Atomic Level Data for Ion: ',elem(at),rom(io)
        write (*,*) 'Atomic Level Data for Ion: ',elem(at),rom(io)
        nt=nfmtrans(idx)
        nl=fmnl(idx)
        write (luop,*) nl,' Levels:'
   90 format (/' Level  Term      g       E (cm^-1)'/)
  100 format ( 2x,i2,4x, a6,x,1pg9.1,x,1pg14.6)
  110 format (/' Level,  Term,      g,       E (cm^-1)'/)
  120 format ( i2,',', a6,',',1pg9.1,',',1pg14.6)
        write (luop,110)
        write (*,90)
        do j=1,nl
          write (luop,120) j,nfmterm(j,idx),wim(j,idx),eim(j,1,idx)/
     &     (plk*cls)
          write (*,100) j,nfmterm(j,idx),wim(j,idx),eim(j,1,idx)/(plk*
     &     cls)
        enddo
      endif
c
      if (linetype.eq.2) then
  130 format(//' Choose a multi-level iron species :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,130)
        do i=1,(nfeions/5)+1
          j=(i-1)*5
          kmax=min(j+5,nfeions)-j
          write (*,'(5(2x,i2,": ",a2,a6))') (j+k,elem(featom(j+k)),
     &     rom(feion(j+k)),k=1,kmax)
        enddo
  140 format(' :: ',$)
        write (*,140)
c
        read (*,*) idx
        if (idx.lt.1) idx=1
        if (idx.gt.nfeions) idx=nfeions
        at=featom(idx)
        io=feion(idx)
c
        do i=1,atypes
          do j=1,maxion(i)
            pop(j,i)=0.0d0
          enddo
        enddo
c
        pop(io,at)=1.000d0
c
        write (luop,*) 'Atomic Level Data for Ion: ',elem(at),rom(io)
        write (*,*) 'Atomic Level Data for Ion: ',elem(at),rom(io)
        nt=nfetrans(idx)
        nl=fenl(idx)
        write (luop,*) nl,' Levels:'
        write (luop,110)
        write (*,90)
        do j=1,nl
          write (luop,120) j,nfeterm(j,idx),wife(j,idx),eife(j,1,idx)/
     &     (plk*cls)
          write (*,100) j,nfeterm(j,idx),wife(j,idx),eife(j,1,idx)/(plk*
     &     cls)
        enddo
c
  150  format(//' Select Line 1 by level:',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  lower, upper  : ',$)
  160  format(i2,x,i2,x,1pg14.7)
        write (*,150)
        read (*,*) i0,j0
        itr0=nfetridx(i0,j0,idx)
        write (*,160) i0,j0,felam(itr0,idx)*1.d8
c
  170 format(//' Select Line 2 by level:',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  lower, upper  : ',$)
        write (*,170)
        read (*,*) i1,j1
        itr1=nfetridx(i1,j1,idx)
        write (*,160) i1,j1,felam(itr1,idx)*1.d8
c
        write (*,*)
        write (*,*)
      endif
c
      if (linetype.eq.1) then
  180 format(//' Choose a recombination species :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,180)
        write (*,'(5(2x,i2,": ",a2,a6))') 1,elem(zmap(1)),rom(1),2,
     &   elem(zmap(2)),rom(1),3,elem(zmap(2)),rom(2),4,elem(zmap(6)),
     &   rom(2),5,elem(zmap(7)),rom(2),6,elem(zmap(8)),rom(1),7,
     &   elem(zmap(8)),rom(2),8,elem(zmap(10)),rom(2)
  190 format(' :: ',$)
        write (*,190)
c
        read (*,*) idx
c
        at=zmap(1)
        io=1
        if (idx.lt.1) idx=1
        if (idx.gt.8) idx=8
        if (idx.eq.1) then
          at=zmap(1)
          io=1
        endif
        if (idx.eq.2) then
          at=zmap(2)
          io=1
        endif
        if (idx.eq.3) then
          at=zmap(2)
          io=2
        endif
        if (idx.eq.4) then
          at=zmap(6)
          io=2
        endif
        if (idx.eq.5) then
          at=zmap(7)
          io=2
        endif
        if (idx.eq.6) then
          at=zmap(8)
          io=1
        endif
        if (idx.eq.7) then
          at=zmap(8)
          io=2
        endif
        if (idx.eq.8) then
          at=zmap(10)
          io=2
        endif
c
      endif
c
c  optically thin:
c
      dr=1.d0
c
      if (modeltype.eq.0) then
c
c Fixed Te
c
  200 format(//' Give fixed conditions:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    T (K) ( <= 10, T as a log10 > 10 linear),'/
     &' :: ',$)
        write (*,200)
        read (*,*) t
        if (t.le.10.d0) t=10.d0**t
        cab=0.d0
  210 format(//' Give density (de = dh) variation:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    dh_min, delta, n:'/
     &'    ( dh_min, delta taken as log10)'/
     &'    ( dh = 10**(dh_min + (delta)*i), i = 0..n)'/
     &' :: ',$)
        write (*,210)
c
        read (*,*) dhmin,dhdelta,ndhsteps
c
  220 format (a)
  230 format(//' Give a name/code for this run: ',$)
        write (*,230)
        read (*,220) runname
        write (luop,10) theversion,runname
        write (*,*)
c
  240 format (/' Transition Data at :',1pg11.3,' K'//
     &,' i, j,  Eij,         lamij,       Upsij,       Aji, ',/)
  250 format (/' Transition Data at :',1pg11.3,' K'//
     &,' i j  Eij         lamij       Upsij       Aji ',/)
  260 format (i2,',',i2,4(',',1pg14.7))
  270 format (i2,x,i2,4(x,1pg14.7))
        write (luop,240) t
        write (*,250) t
        nt=nfmtrans(idx)
        if (linetype.eq.0) then
          nt=nfmtrans(idx)
          do itr=1,nt
            i=nfmlower(itr,idx)
            j=nfmupper(itr,idx)
c          if ( aim(j,i,idx).gt.0.d0)then
            write (luop,260) i,j,eim(i,j,idx),fmlam(itr,idx)*1.d8,
     &       fupsilontr(t,itr,idx),aim(j,i,idx)
            write (*,270) i,j,eim(i,j,idx),fmlam(itr,idx)*1.d8,
     &       fupsilontr(t,itr,idx),aim(j,i,idx)
c          endif
          enddo
        endif
        if (linetype.eq.2) then
          nt=nfetrans(idx)
          do itr=1,nt
            i=nfelower(itr,idx)
            j=nfeupper(itr,idx)
c          if ( aife(j,i,idx).gt.0.d0)then
            write (luop,260) i,j,eife(i,j,idx),felam(itr,idx)*1.d8,
     &       ffeupsilontr(t,itr,idx),aife(j,i,idx)
            write (*,270) i,j,eife(i,j,idx),felam(itr,idx)*1.d8,
     &       ffeupsilontr(t,itr,idx),aife(j,i,idx)
c          endif
          enddo
        endif
c
c loop over density range
c
        write (luop,*)
        write (luop,*) 'Emissivity for Ion: ',elem(at),rom(io)
        write (*,*)
        write (*,*) 'Emissivity for Ion: ',elem(at),rom(io)
        if (linetype.eq.0) then
          write (luop,'(/,a,38(x,1pg14.7,","))') 'T,  ne,',(1.0d8*
     &     fmlam(itr,idx),itr=1,nt)
          write (*,'(/,a,38(x,1pg14.7," "))') 'T,  ne,',(1.0d8*
     &     fmlam(itr,idx),itr=1,nt)
        endif
        if (linetype.eq.2) then
          write (luop,'(/,a,2(x,1pg14.7,","),a)') 'T,  ne,',(1.0d8*
     &     felam(itr0,idx)),(1.0d8*felam(itr1,idx)),'Line0/Line1'
          write (*,'(/,a,2(x,1pg14.7," "),a)') 'T    ne ',(1.0d8*
     &     felam(itr0,idx)),(1.0d8*felam(itr1,idx)),'Line0/Line1'
        endif
        do i=0,ndhsteps
          dh=10.d0**(dhdelta*i+dhmin)
          de=dh
          call zerbuf
          pop(io,at)=1.0000
          call cool (t, de, dh)
          if (linetype.eq.0) then
            do itr=1,nt
              br(itr)=fmbri(itr,idx)
            enddo
            write (luop,'(38(x,1pg14.7,","))') t,de,(br(itr),itr=1,nt)
            write (*,'(38(x,1pg14.7,","))') t,de,(br(itr),itr=1,nt)
          endif
          if (linetype.eq.2) then
            br0=febri(itr0,idx)
            br1=febri(itr1,idx)
            write (luop,'(5(x,1pg14.7,","))') t,de,br0,br1,br0/br1
            write (*,'(5(x,1pg14.7))') t,de,br0,br1,br0/br1
          endif
        enddo
        close (luop)
c
c Fixed Te
c
      endif
c
      if (modeltype.eq.1) then
c
c Fixed de
c
  280 format(//' Give fixed conditions:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    nH, (<=20 as log10. > 20.d0 as linear ne = nH)'//
     &' :: ',$)
        write (*,280)
        read (*,*) dh
        if (dh.gt.20.d0) dh=dlog10(dh)
        dh=10.d0**(dh)
        de=dh
        cab=0.d0
  290 format(//' Give Temperature variation:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    T_min, delta, n:'/
     &'    ( T_min <= 10 as log, > 10 linear)'/
     &'    ( delta always taken as log10 )'/
     &'    ( T = 10**(T_min + (delta)*i), i = 0..n)'/
     &' :: ',$)
        write (*,290)
c
        read (*,*) tmin,tdelta,ntsteps
        if (tmin.gt.10.d0) tmin=dlog10(tmin)
        if (ntsteps.gt.101) ntsteps=101
c
  300 format (a)
  310 format(//' Give a name/code for this run: ',$)
        write (*,310)
        read (*,300) runname
        write (luop,10) theversion,runname
        write (*,*)
c
        if (linetype.eq.0) then
c
c Collisional Lines
c
  320 format (/' Transition Data :'//
     &  ,' i, j,  Eij,         lamij,       Aji, ',/)
  330 format (/' Transition Data :'//
     &  ,' i j  Eij         lamij       Aji ',/)
  340 format (i2,',',i2,3(',',1pe11.4))
  350 format (i2,i2,3(x,1pe11.4))
          write (luop,320)
          write (*,330)
          nt=nfmtrans(idx)
          if (linetype.eq.0) then
            nt=nfmtrans(idx)
            do itr=1,nt
              i=nfmlower(itr,idx)
              j=nfmupper(itr,idx)
              write (luop,340) i,j,eim(i,j,idx),fmlam(itr,idx)*1.d8,
     &         aim(j,i,idx)
              write (*,350) i,j,eim(i,j,idx),fmlam(itr,idx)*1.d8,aim(j,
     &         i,idx)
            enddo
          endif
          if (linetype.eq.2) then
            nt=nfetrans(idx)
            do itr=1,nt
              i=nfelower(itr,idx)
              j=nfeupper(itr,idx)
              write (luop,340) i,j,eife(i,j,idx),felam(itr,idx)*1.d8,
     &         aife(j,i,idx)
              write (*,350) i,j,eife(i,j,idx),felam(itr,idx)*1.d8,
     &         aife(j,i,idx)
            enddo
          endif
c
c loop over temprature range for upsilons
c
          write (*,*)
          write (*,*) 'Collision Strengths, Ion: ',elem(at),rom(io)
          write (luop,*)
          write (luop,*) 'Collision Strengths, Ion: ',elem(at),rom(io)
          if (usekappa) write (*,*) 'Kappa:',nt
          if (usekappa) write (luop,*) 'Kappa:',nt
c
  360 format(/,a,36(x,a13,','))
  370 format(a,36(x,1pg13.6,','))
  380 format(39(x,1pe13.6,','))
          if (linetype.eq.0) then
            do itr=1,nt
              s=nfmid(itr,idx)
              sl=len(s)
              if (sl.gt.13) then
                trimid(itr)=s(sl-12:sl)
              else
                trimid(itr)=s(1:sl)
              endif
            enddo
            write (*,360) '              ,              ,              ,
     &',(trimid(itr),itr=1,nt)
            write (*,370) '          LogT,             T,            ne,
     &',(1.0d8*fmlam(itr,idx),itr=1,nt)
            write (luop,360) '              ,              ,
     &  ,',(trimid(itr),itr=1,nt)
            write (luop,370) '          LogT,             T,
     &ne,',(1.0d8*fmlam(itr,idx),itr=1,nt)
          endif
          if (linetype.eq.2) then
            do itr=1,nt
              s=nfeid(itr,idx)
              sl=len(s)
              if (sl.gt.13) then
                trimid(itr)=s(sl-12:sl)
              else
                trimid(itr)=s(1:sl)
              endif
            enddo
            write (*,360) '              ,              ,              ,
     &',(trimid(itr),itr=1,nt)
            write (*,370) '          LogT,             T,            ne,
     &',(1.0d8*felam(itr,idx),itr=1,nt)
            write (luop,360) '              ,              ,
     &  ,',(trimid(itr),itr=1,nt)
            write (luop,370) '          LogT,             T,
     &ne,',(1.0d8*felam(itr,idx),itr=1,nt)
          endif
          do i=0,ntsteps
            logt=(tdelta*i+tmin)
            t=10.d0**(tdelta*i+tmin)
            if (linetype.eq.0) then
              do itr=1,nt
                ups(itr)=fupsilontr(t,itr,idx)
              enddo
              write (*,380) logt,t,de,(ups(itr),itr=1,nt)
              write (luop,380) logt,t,de,(ups(itr),itr=1,nt)
            endif
            if (linetype.eq.2) then
              do itr=1,nt
                ups(itr)=ffeupsilontr(t,itr,idx)
              enddo
              write (*,380) logt,t,de,(ups(itr),itr=1,nt)
              write (luop,380) logt,t,de,(ups(itr),itr=1,nt)
            endif
          enddo
c
c loop over temprature range
c
          write (luop,*)
          write (luop,*) 'Emissivity for Ion: ',elem(at),rom(io)
          write (*,*)
          write (*,*) 'Emissivity for Ion: ',elem(at),rom(io)
          if (linetype.eq.0) then
            do itr=1,nt
              s=nfmid(itr,idx)
              sl=len(s)
              if (sl.gt.13) then
                trimid(itr)=s(sl-12:sl)
              else
                trimid(itr)=s(1:sl)
              endif
            enddo
            write (*,360) '              ,              ,              ,
     &',(trimid(itr),itr=1,nt)
            write (*,370) '          LogT,             T,            ne,
     &',(1.0d8*fmlam(itr,idx),itr=1,nt)
            write (luop,360) '              ,              ,
     &  ,',(trimid(itr),itr=1,nt)
            write (luop,370) '          LogT,             T,
     &ne,',(1.0d8*fmlam(itr,idx),itr=1,nt)
          endif
          if (linetype.eq.2) then
            do itr=1,nt
              s=nfeid(itr,idx)
              sl=len(s)
              if (sl.gt.13) then
                trimid(itr)=s(sl-12:sl)
              else
                trimid(itr)=s(1:sl)
              endif
            enddo
            write (*,360) '              ,              ,              ,
     &',(trimid(itr),itr=1,nt)
            write (*,370) '          LogT,             T,            ne,
     &',(1.0d8*felam(itr,idx),itr=1,nt)
            write (luop,360) '              ,              ,
     &  ,',(trimid(itr),itr=1,nt)
            write (luop,370) '          LogT,             T,
     &ne,',(1.0d8*felam(itr,idx),itr=1,nt)
          endif
          do i=0,ntsteps
            logt=(tdelta*i+tmin)
            t=10.d0**logt
            call zerbuf
            pop(io,at)=1.0000
            call cool (t, de, dh)
            if (linetype.eq.0) then
              do itr=1,nt
                br(itr)=fmbri(itr,idx)
              enddo
              write (*,380) logt,t,de,(br(itr),itr=1,nt)
              write (luop,380) logt,t,de,(br(itr),itr=1,nt)
            endif
            if (linetype.eq.2) then
              do itr=1,nt
                br(itr)=febri(itr,idx)
              enddo
              write (*,380) logt,t,de,(br(itr),itr=1,nt)
              write (luop,380) logt,t,de,(br(itr),itr=1,nt)
            endif
          enddo
          close (luop)
c
c Collisional Lines
c
        endif
c
        if (linetype.eq.1) then
c
c Recombination Lines
c
          write (luop,*) 'Emissivity for Ion: ',elem(at),rom(io)
          write (*,*) 'Emissivity for Ion: ',elem(at),rom(io)
c
c Recombination
c
  390 format (/' H II Transitions :'//
     &,' Line, Series, lambda, Emiss. ',/)
  400 format (i3,', ',i3,', ',1pe14.6)
  410 format (i3,',',1pe14.6',',1pe14.6)
  420 format (201(1pe14.6,','))
c
          if (at.eq.zmap(1)) then
            write (luop,390)
            write (*,390)
            do series=1,nhseries
              do line=1,nhlines
                write (luop,400) line,series,hlambda(line,series)
              enddo
            enddo
          endif
c
  430 format (/' He II Transitions :'//
     &,' Line, Series, lambda, Emiss. ',/)
c
          if (at.eq.zmap(2)) then
            if (io.eq.2) then
              write (luop,430)
              write (*,430)
              do series=1,nheseries
                do line=1,nhelines
                  write (luop,400) line,series,helambda(line,series)
                enddo
              enddo
            endif
          endif
c
  440 format (/' He I Transitions :'//
     &,' Line, lambda, Emiss. ',/)
c
          if (at.eq.zmap(2)) then
            if (io.eq.1) then
              write (luop,440)
              write (*,440)
              write (luop,*) ' HeI Legacy Lines:'
              write (*,*) ' HeI Legacy Lines:'
              do line=1,3
                write (luop,410) line,(heilam(line)*1.d8)
              enddo
              write (luop,*) ' HeI Singlets:'
              write (*,*) ' HeI Singlets:'
              do line=1,nheislines
                write (luop,410) line,heislam(line),heiseij(line)
              enddo
              write (luop,*) ' HeI Triplets:'
              write (*,*) ' HeI Triplets:'
              do line=1,nheitlines
                write (luop,410) line,heitlam(line),heiteij(line)
              enddo
            endif
          endif
c
  450 format (/' C II Transitions :'//
     &,' Line, lambda, Eij. ',/)
c
          if (at.eq.zmap(6)) then
            if (io.eq.2) then
              write (luop,450)
              write (*,450)
              do line=1,nrccii
                write (luop,410) line,rccii_lam(line),rccii_eij(line)
              enddo
            endif
          endif
c
  460 format (/' N II Transitions :'//
     &,' Line, lambda, Eij. ',/)
c
          if (at.eq.zmap(7)) then
            if (io.eq.2) then
              write (luop,460)
              write (*,460)
              do line=1,nrcnii
                write (luop,410) line,rcnii_lam(line),rcnii_eij(line)
              enddo
            endif
          endif
c
  470 format (/' O I Transitions :'//
     &,' Line, lambda, Eij. ',/)
c
          if (at.eq.zmap(8)) then
            if (io.eq.1) then
              write (luop,470)
              write (*,470)
              do line=1,nrcoi_q
                write (luop,410) line,rcoi_qlam(line),rcoi_qeij(line)
              enddo
              do line=1,nrcoi_t
                write (luop,410) line,rcoi_tlam(line),rcoi_teij(line)
              enddo
            endif
          endif
c
  480 format (/' O II Transitions :'//
     &,' Line, lambda, Eij. ',/)
c
          if (at.eq.zmap(8)) then
            if (io.eq.2) then
              write (luop,480)
              write (*,480)
              do line=1,nrcoii
                write (luop,410) line,rcoii_lam(line),rcoii_eij(line)
              enddo
            endif
          endif
c
c
          if (at.eq.zmap(10)) then
            if (io.eq.2) then
              do line=1,nrcneii
                write (luop,410) line,rcneii_lam(line),rcneii_eij(line)
              enddo
            endif
          endif
          do i=0,ntsteps
c
c loop over temperature range
c
            t=10.d0**(tdelta*i+tmin)
            call zerbuf
            pop(io+1,at)=1.0000
            call cool (t, de, dh)
c
c Recombination
c
            if (at.eq.zmap(1)) then
              write (luop,390)
              write (*,390)
              itr=0
              do series=1,nhseries
                do line=1,nhlines
                  itr=itr+1
                  emiss(itr,i)=hydrobri(line,series)
                enddo
              enddo
            endif
c
            if (at.eq.zmap(2)) then
              if (io.eq.2) then
                itr=0
                do series=1,nheseries
                  do line=1,nhelines
                    itr=itr+1
                    emiss(itr,i)=helibri(line,series)
                  enddo
                enddo
              endif
            endif
c
            if (at.eq.zmap(2)) then
              if (io.eq.1) then
                itr=0
                do line=1,3
                  itr=itr+1
                  emiss(itr,i)=heibri(line)
                enddo
                do line=1,nheislines
                  itr=itr+1
                  emiss(itr,i)=heisbri(line)
                enddo
                do line=1,nheitlines
                  itr=itr+1
                  emiss(itr,i)=heitbri(line)
                enddo
              endif
            endif
c
            if (at.eq.zmap(6)) then
              if (io.eq.2) then
                itr=0
                do line=1,nrccii
                  itr=itr+1
                  emiss(itr,i)=rccii_bbri(line)
                enddo
              endif
            endif
c
            if (at.eq.zmap(7)) then
              if (io.eq.2) then
                itr=0
                do line=1,nrcnii
                  itr=itr+1
                  emiss(itr,i)=rcnii_bbri(line)
                enddo
              endif
            endif
c
            if (at.eq.zmap(8)) then
              if (io.eq.1) then
                itr=0
                do line=1,nrcoi_q
                  itr=itr+1
                  emiss(itr,i)=rcoi_qbbri(line)
                enddo
                do line=1,nrcoi_t
                  itr=itr+1
                  emiss(itr,i)=rcoi_tbbri(line)
                enddo
              endif
            endif
c
            if (at.eq.zmap(8)) then
              if (io.eq.2) then
                itr=0
                do line=1,nrcoii
                  itr=itr+1
                  emiss(itr,i)=rcoii_bbri(line)
                enddo
              endif
            endif
c
            if (at.eq.zmap(10)) then
              if (io.eq.2) then
                itr=0
                do line=1,nrcneii
                  itr=itr+1
                  emiss(itr,i)=rcneii_bbri(line)
                enddo
              endif
            endif
c
c loop over temperature range
c
          enddo
          write (luop,*) ' Te(K), lines....'
          do i=0,ntsteps
            t=10.d0**(tdelta*i+tmin)
            write (luop,420) t,(emiss(idx,i),idx=1,itr)
          enddo
          close (luop)
c
c Recombination Lines
c
        endif
c
c Fixed de
c
      endif
c
  490 format(/
     &' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'  OUTPUT WRITTEN TO : ', a,/
     &' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
      write (*,490) fl
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subject a single ion to a temperature and density range,
c     targeting two lines, and outputting de-excitation information
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine critdens
c
c
      include 'cblocks.inc'
c
      integer*4 luop
      integer*4 i,j,k,idx,kmax,stepidx
      integer*4 i0,j0,i1,j1,itr0,itr1
      integer*4 nl
      integer*4 at,io
c      integer*4 line,series
      integer*4 nt, itr
      integer*4 ndhsteps
      integer*4 modeltype
      integer*4 linetype
c      character str*18
      character fn*32
      character fl*32
      character pfx*16
      character sfx*4
c      character ilgg*4
      character imod*4
      character model*32
c
      real*8 t,de, dh, dr
      real*8 cab
      real*8 br(3)
      real*8 f,nmin,nmax,ndelta,a0,a1
      real*8 ncrit0,ncrit1,n0,n1
      real*8 ups0,ups1,aa
      real*8 colexrate0,colexrate1
      real*8 coldexrate0,coldexrate1
      real*8 wlower,wupper
c
c      real*8 fupsilontr
c      real*8 ffeupsilontr
      real*8 fupsilonij
c
   10 format(//' Single Multi-Level Ion Critical Density Model'/
     &' produced by MAPPINGS V ',a8,'  Run:',a80/)
      model='CD ion model'
      luop=23
      imod='ALL'
c
      fn=' '
      pfx='ion'
      sfx='csv'
      call newfile (pfx, 3, sfx, 3, fn)
      fl=fn(1:13)
      open (luop,file=fl,status='NEW')
c
      call zer
      t=1.0d4
      de=1.0d0
c
      itr0=1
      itr1=2
      at=zmap(1)
      io=1
c
      modeltype=0
      linetype=0
c
      if (linetype.eq.0) then
   20 format(//' Choose a multi-level species :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,20)
        do i=1,(nfmions/5)+1
          j=(i-1)*5
          kmax=min(j+5,nfmions)-j
          write (*,'(5(2x,i3,": ",a2,a6))') (j+k,elem(fmatom(j+k)),
     &     rom(fmion(j+k)),k=1,kmax)
        enddo
   30 format(' :: ',$)
        write (*,30)
c
        read (*,*) idx
        if (idx.lt.1) idx=1
        if (idx.gt.nfmions) idx=nfmions
        at=fmatom(idx)
        io=fmion(idx)
c
        write (*,*)
        write (*,*)
c
        do i=1,atypes
          do j=1,maxion(i)
            pop(j,i)=0.0d0
          enddo
        enddo
c
        pop(io,at)=1.000d0
c
        nt=nfmtrans(idx)
        nl=fmnl(idx)
   40   format('% Atomic Level Data for Ion: ',a2,a6,/,
     &  '% ',i4,' Levels:')
        write (luop,40) elem(at),rom(io),nl
        write (*,40) elem(at),rom(io),nl
   50   format ( '%  Level  Term      g       E (cm^-1)')
   60   format ( '% ',i2,4x, a6,x,1pg9.1,x,1pg14.6)
   70   format ( '%  Level,  Term,      g,       E (cm^-1)')
   80   format ( '% ',i2,',', a6,',',1pg9.1,',',1pg14.6)
        write (luop,70)
        write (*,50)
        do j=1,nl
          write (luop,80) j,nfmterm(j,idx),wim(j,idx),eim(j,1,idx)/(plk*
     &     cls)
          write (*,60) j,nfmterm(j,idx),wim(j,idx),eim(j,1,idx)/(plk*
     &     cls)
        enddo
      endif
c
   90  format(//' Select Line 1 by level:',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  lower, upper  : ',$)
  100  format(i2,x,i2,x,1pg14.7)
      write (*,90)
      read (*,*) i0,j0
      itr0=nfmtridx(i0,j0,idx)
      write (*,100) i0,j0,fmlam(itr0,idx)*1.d8
c
  110 format(//' Select Line 2 by level:',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  lower, upper  : ',$)
      write (*,110)
      read (*,*) i1,j1
      itr1=nfmtridx(i1,j1,idx)
      write (*,100) i1,j1,fmlam(itr1,idx)*1.d8
c
      write (*,*)
c
c  optically thin:
c
      dr=1.d0
c
      if (modeltype.eq.0) then
c
c Fixed Te
c
  120 format(//' Give fixed conditions:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    T (K) ( <= 10, T as a log10 > 10 linear),'/
     &' :: ',$)
        write (*,120)
        read (*,*) t
        if (t.le.10.d0) t=10.d0**t
        f=1.d0/dsqrt(t)
c
        cab=0.d0
  130 format(//' Give density (n = de = dh) variation:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    n_min, n_max, delta:'/
     &'    ( as log10)'/
     &' :: ',$)
        write (*,130)
c
        read (*,*) nmin,nmax,ndelta
        ndhsteps=idnint((nmax-nmin)/ndelta)+1
c
  140 format (a)
  150 format(//' Give a name/code for this run: ',$)
        write (*,150)
        read (*,140) runname
        write (luop,10) theversion,runname
        write (*,*)
c
  160 format ('%  Transition Data at :',1pg11.3,' K',/,
     & '%  i, j,           Eij,         lamij,         Upsij,',
     & '           Aji,             f, gi, gj,  ne_crit')
        write (*,160) t
        write (luop,160) t
  170 format
     &('% ',i2,',',i2,5(',',1pg14.7),', '1pf5.2,', ',1pf5.2,',',1pg14.7)
        i=nfmlower(itr0,idx)
        j=nfmupper(itr0,idx)
        ups0=fupsilonij(t,i,j,idx)
        wlower=wim(i,idx)
        wupper=wim(j,idx)
        a0=aim(j,i,idx)
        ncrit0=(a0*wupper)/(rka*f*ups0)
        write (luop,170) i,j,eim(i,j,idx),fmlam(itr0,idx)*1.d8,ups0,a0,
     &   f,wlower,wupper,ncrit0
        write (*,170) i,j,eim(i,j,idx),fmlam(itr0,idx)*1.d8,ups0,a0,f,
     &   wlower,wupper,ncrit0
        i=nfmlower(itr1,idx)
        j=nfmupper(itr1,idx)
        ups1=fupsilonij(t,i,j,idx)
        wlower=wim(i,idx)
        wupper=wim(j,idx)
        a1=aim(j,i,idx)
        ncrit1=(a1*wupper)/(rka*f*ups1)
        write (luop,170) i,j,eim(i,j,idx),fmlam(itr1,idx)*1.d8,ups1,a1,
     &   f,wlower,wupper,ncrit1
        write (*,170) i,j,eim(i,j,idx),fmlam(itr1,idx)*1.d8,ups1,a1,f,
     &   wlower,wupper,ncrit1
c
c loop over density range
c
  180   format('% Fluxes and Ratios for Ion:',a2,a6)
        write (luop,180) elem(at),rom(io)
        write (*,180) elem(at),rom(io)
  190  format(a16,a16,2(', ',1pg14.7),a16,2(a48,a32))
        write (luop,190) '%      ne       ',',ne.(10^4/T)^0.5',1.0d8*
     &   fmlam(itr0,idx),1.0d8*fmlam(itr1,idx),',     ratio     ',',
     &  alpha12_0,      alpha21_0,          A21_0',',         pop1_0,
     &      pop2_0',',      alpha12_1,      alpha21_1,          A21_1','
     &,         pop1_1,         pop2_1'
c
        write (*,190) '%      ne       ',',ne.(10^4/T)^0.5',1.0d8*
     &   fmlam(itr0,idx),1.0d8*fmlam(itr1,idx),',     ratio     ',',
     &  alpha12_0,      alpha21_0,          A21_0',',         pop1_0,
     &      pop2_0',',      alpha12_1,      alpha21_1,          A21_1','
     &,         pop1_1,         pop2_1'
c
        do stepidx=0,ndhsteps
          dh=10.d0**(nmin+stepidx*ndelta)
          de=dh
c
c
          call zerbuf
          pop(io,at)=1.0000
          call cool (t, de, dh)
          br(1)=fmbri(itr0,idx)
          br(2)=fmbri(itr1,idx)
          br(3)=br(2)/br(1)
c
          itr=itr0
          i=nfmlower(itr,idx)
          j=nfmupper(itr,idx)
          wlower=wim(i,idx)
          wupper=wim(j,idx)
          aa=dexp(-eim(i,j,idx)/(rkb*t))
          colexrate0=de*rka*f*aa*(ups0/wlower)
          coldexrate0=de*((rka*f)*ups0)/wupper
          a0=aim(j,i,idx)
          n0=colexrate0/(a0+coldexrate0)
c
          itr=itr1
          i=nfmlower(itr,idx)
          j=nfmupper(itr,idx)
          wlower=wim(j,idx)
          wupper=wim(j,idx)
          colexrate1=de*rka*f*aa*(ups1/wlower)
          coldexrate1=de*((rka*f)*ups1)/wupper
          n1=colexrate1/(a1+coldexrate1)
c
  200     format('  ',1pg14.7,14(', ',1pg14.7))
          write (luop,200) de,de*((1d4/t)**0.5d0),(br(itr),itr=1,3),
     &     colexrate0,coldexrate0,a0,fmx(nfmlower(itr0,idx),idx),
     &     fmx(nfmupper(itr0,idx),idx),colexrate1,coldexrate1,a1,
     &     fmx(nfmlower(itr1,idx),idx),fmx(nfmupper(itr1,idx),idx)
          write (*,200) de,de*((1d4/t)**0.5d0),(br(itr),itr=1,3),
     &     colexrate0,coldexrate0,a0,fmx(nfmlower(itr0,idx),idx),
     &     fmx(nfmupper(itr0,idx),idx),colexrate1,coldexrate1,a1,
     &     fmx(nfmlower(itr1,idx),idx),fmx(nfmupper(itr1,idx),idx)
        enddo
        close (luop)
c
c Fixed Te
c
      endif
c
  210 format(/
     &' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'  OUTPUT WRITTEN TO : ', a,/
     &' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
      write (*,210) fl
      return
      end
