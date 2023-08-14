cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******PHOTOIONISATIONMODEL
c     DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION, RADIATION
c     PRESSURE
c
c     space steps derived from optical depth of predicted temperature
c     and ionistation, extrapolated from previous steps
c
c     CHOICE OF EQUILIBRIUM OR FINITE AGE CONDITIONS
c
c     NB. COMPUTATIONS PERFORMED IN SUBROUTINE COMPPH6
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine montph7 ()
c
      include 'cblocks.inc'
      include 'mpif.h'
c
c           Variables
c
      real*8 blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      real*8 epotmi,fin,uinit, qinit, rcin
      real*8 rstsun,wid,starlum
      real*8 scale,dnt,lam
      real*8 blum0,dflum0,dflum(8)
      real*8 mcdx,mcdy,mcdz
      real*8 dhntmp, tmp
c      real*8 fluxes(mxmonlines)
c
      integer*4 lenv      
      integer*4 l,m,n
      integer*4 ix,iy,iz
c
      integer*4 j,luin,i,kmax,k,nl
      integer*4 idx1 ,idx2, idx3, idx4, iem      
      integer*4 numstar
      integer*4 nblk(3), blockid(3)
      integer*4 mcnx0,mcny0,mcnz0
c            
c      character=jtb*4
      character carac*36,model*32
      character caract*4,ilgg*4
      character banfil*32
      character*512 fnam
      character*512 grdfile,denfile
      character jren*4,jrnu*4
c
      logical iexi      
c
c           Functions
c
      real*8 densnum,fdilu
c
      luin=10
      jcon='YES'
      jspot='NO'
      dtau0=0.025d0
      nprefix=6
      numstar=0
c
      write (*,10)
c
c     ***INITIAL IONISATION CONDITIONS
c
   10 format(///
     &     ' *********************************************************'/
     &     '  Photoionisation model Monte-Carlo selected:',/,
     &     ' *********************************************************')
c
      epotmi=iphe
c
c     set up ionisation state
c
      model='photoionisation'
      call popcha (model)
      model='Photo 6'
c
c     set random seed for the simulation
      seed = time()
c
c     artificial support for low ionisation species
c
      if (expertmode.gt.0) then
c
   20    format(a)
        carac(1:33)='   '
        j=1
        do 30 i=3,atypes
          j=1+index(carac(1:33),'   ')
          caract='   '
          if (ipote(1,i).lt.epotmi) write (caract,20) elem(i)
          carac(j:j+2)=caract//'   '
   30   continue
        i=j+2
   40   write (*,50) carac(1:i)
        ilgg='    '
c
   50    format(//' Allow the elements : ',a/
     &        ' to recombine freely? (y/n) : ',$)
c
        if (taskid.eq.0) read (*,20) ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
        ilgg=ilgg(1:1)
c
        if (ilgg.eq.'Y') goto 70
        if (ilgg.eq.'y') goto 70
        if ((ilgg.ne.'N').and.(ilgg.ne.'n')) goto 40
c
        do 60 i=3,atypes
          arad(2,i)=dabs(arad(2,i))
          if (ipote(1,i).lt.epotmi) arad(2,i)=-arad(2,i)
   60   continue
   70   continue
c
      endif
c
c
      write (*,75)
   75 format(/' Enter block in x, y & z axis : ',$)  
      if (taskid.eq.0) read (*,*) nblk(1), nblk(2), nblk(3)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(nblk,3,MPI_INTEGER4,0,MPI_COMM_WORLD,taskerr)
c
      if (nblk(1)*nblk(2)*nblk(3).lt.numtasks) then
         write (*,*) 'Wrong Block Schema'
         stop      
      endif
c

      blockid(3) = taskid/(nblk(1)*nblk(2))
      blockid(2) = (taskid-nblk(1)*nblk(2)*blockid(3))
     &             /nblk(1)
      blockid(1) = taskid-nblk(1)*nblk(2)*blockid(3)
     &            -nblk(1)*blockid(2)
c
c      blockid(1) = taskid/(nblk(3)*nblk(2))
c      blockid(2) = (taskid-nblk(3)*nblk(2)*blockid(1))
c     &             /nblk(3)
c      blockid(3) = taskid-nblk(3)*nblk(2)*blockid(1)
c     &            -nblk(3)*blockid(2)
c
cc      write (*,*) taskid, blockid
cc      call mpi_barrier(MPI_COMM_WORLD, taskerr)
cc      stop           
c
cc      do ix=1, nblk(1) do begin
cc      do iy=1, nblk(2) do begin
cc      do iz=1, nblk(3) do begin
cc         blockid = iz+(iy-1)*nblk(3)+(ix-1)*nblk(3)*nblk(2)
cc         if (blockid.eq.taskid) 
cc      enddo
cc      enddo
cc      enddo
c
c     ***SET 3D GEOMETRY
c
      write (*,80)
   80 format(/' Enter grid file name : ',$)   
c
      if (taskid.eq.0) read (*,90) fnam      
   90 format(a)
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(fnam,512,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
      write (*,*)
      m=lenv(fnam)
      grdfile=fnam(1:m)
      iexi=.false.
      inquire (file=grdfile(1:m),exist=iexi)  
c
      if (iexi.eqv..false.) then
c
         write (*,*) grdfile(1:m),' NOT FOUND.'
c
      else 
c
c     FORMAT REQUIRED:
c
c     1) number of cell wall
c     2) the cell wall value
c     3) cells are written in order as x, y and z
c
c
         write (*,*) 'FOUND: ',grdfile(1:m)
         luin=99
         open (unit=luin,file=grdfile(1:m),status='OLD') 
c
         read (luin,*) mcnx0
         mcnx = mcnx0/nblk(1)
         l = 1
         do i=1,mcnx0+1            
            read (luin,*) tmp
            if ((i.ge.blockid(1)*mcnx+1)
     &     .and.(i.le.((blockid(1)+1)*mcnx+1))) then
             cwx(l) = tmp
             l = l + 1
            endif
         enddo 
c         
         read (luin,*) mcny0
         mcny = mcny0/nblk(2)
         l = 1         
         do i=1,mcny0+1
            read (luin,*) tmp
            if ((i.ge.blockid(2)*mcny+1)
     &     .and.(i.le.((blockid(2)+1)*mcny+1))) then
             cwy(l) = tmp
             l = l + 1
            endif
         enddo
c                        
         read (luin,*) mcnz0
         mcnz = mcnz0/nblk(3)
         l = 1         
         do i=1,mcnz0+1
            read (luin,*) tmp
            if ((i.ge.blockid(3)*mcnz+1)
     &     .and.(i.le.((blockid(3)+1)*mcnz+1))) then
             cwz(l) = tmp
             l = l + 1
            endif
         enddo         
c                 
      endif    
c      
c
      write(*,95)
   95 format(//' Input Density Grid? (y/N) : ',$)
c
      if (taskid.eq.0) read (*,20) ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'Y') goto 99
      if (ilgg.eq.'y') goto 99
c
      write(*,98)
   98 format(//' Input Density Vaule : ',$)
      if (taskid.eq.0) read (*,*) dhntmp
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(dhntmp,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
      do ix=1,mcnx
      do iy=1,mcny
      do iz=1,mcnz
         dh_arr(ix,iy,iz) = dhntmp
      enddo 
      enddo 
      enddo 
c
      goto 105
c
c
   99 write(*,100)
  100 format(/' Input Density File : ',$)
      if (taskid.eq.0) read (*,20) fnam
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(fnam,512,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)

      m=lenv(fnam)
      denfile=fnam(1:m)
c
      iexi=.false.
      inquire (file=denfile(1:m),exist=iexi)  
c
      if (iexi.eqv..false.) then
c
         write (*,*) denfile(1:m),' NOT FOUND.'
         stop
c
      else 
c
         write (*,*) 'FOUND: ',denfile(1:m)
         luin=99
         open (unit=luin,file=denfile(1:m),status='OLD') 
c
         do i=1,mcnx0
         do j=1,mcny0
         do k=1,mcnz0
            read (luin,*) tmp
            if ((i.ge.(blockid(1)*mcnx+1))
     &     .and.(i.le.((blockid(1)+1)*mcnx))
     &     .and.(j.ge.(blockid(2)*mcny+1))
     &     .and.(j.le.((blockid(2)+1)*mcny))
     &     .and.(k.ge.(blockid(3)*mcnz+1))
     &     .and.(k.le.((blockid(3)+1)*mcnz))) then
c              
            ix = i-blockid(1)*mcnx
            iy = j-blockid(2)*mcny
            iz = k-blockid(3)*mcnz
c            
            dh_arr(ix,iy,iz) = tmp
            write (*,*), ix, iy, iz, i, j, k
c
            endif            
c            
         enddo 
         enddo 
         enddo 
c
      endif
c
  105 write(*,*) 'Density Reading Complete'
c
c     ***CHOOSING PHOTON SOURCE
c
      model='Photo 6'
c
      numstar=numstar+1
c
      call photsou (model)      
c
c     Possible to have no photons
c
      alphahb=2.585d-13
      alphaheb=1.533d-12
c
c     Define freqency range of spectral sampling
c     Default range is from 5-100 eV
c
      mcnuin=5.0d0
      mcnuend=100.0d0
c      
      rstar=1.d0
      astar=rstar*rstar
c
      blum=0.d0
      ilum=0.d0
      dflum0=0.0d0
c
      do i=1,infph-1
         wid=photev(i+1)-photev(i)
         blum=blum+soupho(i)*wid*evplk
         if (photev(i).ge.iph) ilum=ilum+soupho(i)*wid*evplk
         if ((photev(i).ge.mcnuin).and.(photev(i).lt.mcnuend))
     &      dflum0=dflum0+soupho(i)*wid*evplk                    
      enddo
c      
c
      blum0=blum
c      
      blum=pi*blum
      ilum=pi*ilum
      dflum0=pi*dflum0
c
      rstar=1.d0
      astar=rstar*rstar
c
      if (blum.gt.0.d0) then
c
  110    write (*,120)
  120 format(/' Define the source size or luminosity '/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    R   :  By Radius'/
     & '    L   :  By Luminosity'/
     & '    P   :  By Ionizing Photons'//
     & ' :: ',$)
         if (taskid.eq.0) read (*,90) ilgg
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
         ilgg=ilgg(1:1)
c
         if (ilgg.eq.'r') ilgg='R'
         if (ilgg.eq.'l') ilgg='L'
         if (ilgg.eq.'p') ilgg='P'
c
         if ((ilgg.ne.'R').and.(ilgg.ne.'L').and.(ilgg.ne.'P')) goto
     &   110
c
         if (ilgg.eq.'R') then
c
  130      format(/,/,' Give photoionisation source radius'
     &            ,'(Rsun=6.96e10)',/
     &            ,' (in solar units (<1.e5) or in cm (>1.e5) : ',$)
            write (*,130)
            if (taskid.eq.0) read (*,*) rstsun
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(rstsun,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
            if (rstsun.le.0.0d0) goto 110
            if (rstsun.lt.1.d5) then
              rstar=6.96d10*rstsun
            else
              rstar=rstsun
              rstsun=rstar/6.96d10
            endif
c
            astar=fpi*rstar*rstar
            blum=astar*blum
            ilum=astar*ilum
            dflum0=astar*dflum0
c            
c
  140 format(//' ****************************************************'/
     & '  Source Total Luminosity                   : ',1pg12.5/
     & '  Source Ionising (13.6eV+) Luminosity      : ',1pg12.5/
     & '  Source Ionising (13.6eV+) Photons         : ',1pg12.5/
     & '  Source Luminosity ('1pg12.5'-'1pg12.5'eV) : ',1pg12.5/
     & ' ****************************************************')
            write (*,140) blum,ilum,qht*astar
     &                   ,mcnuin,mcnuend,dflum0
c
c     end def by radius
c
          endif
c
          if (ilgg.eq.'L') then
c
  150       format(//' Total or Ionising Luminosity (T/I)',$)
            write (*,150)
            if (taskid.eq.0) read (*,90) ilgg
        call mpi_barrier(MPI_COMM_WORLD, taskerr)
        call mpi_bcast(ilgg,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
            ilgg=ilgg(1:1)
c
            if (ilgg.eq.'t') ilgg='T'
            if (ilgg.eq.'i') ilgg='I'
c
            if ((ilgg.ne.'I').and.(ilgg.ne.'T')) ilgg='I'
c
  160       format(//' Give ionising source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
  170       format(//' Give bolometric source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
            if (ilgg.eq.'T') then
              write (*,170)
            else
              write (*,160)
            endif
            if (taskid.eq.0) read (*,*) starlum
        call mpi_barrier(MPI_COMM_WORLD, taskerr)
        call mpi_bcast(starlum,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
            if (starlum.lt.1.d2) starlum=10.d0**starlum
c
c
            if (ilgg.eq.'I') then
              rstsun=dsqrt(starlum/(ilum*fpi))
            else
              rstsun=dsqrt(starlum/(blum*fpi))
            endif
c
c
  180  format(//' ****************************************************'/
     &      '  Source Radius : ',1pg12.5,' cm'/
     &      ' ****************************************************')
            write (*,180) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            dflum0=astar*dflum0
c            
            write (*,140) blum,ilum,qht*astar
     &                   ,mcnuin,mcnuend,dflum0
c     
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c     end def by luminosity
          endif
c
          if (ilgg.eq.'P') then
c
  190       format(//' Give source ionising photon rate '/
     &               ' (log (<100) or photons/s (>100)) : ',$)
            write (*,190)
            if (taskid.eq.0) read (*,*) starlum
        call mpi_barrier(MPI_COMM_WORLD, taskerr)
        call mpi_bcast(starlum,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
            if (starlum.lt.1.d2) starlum=10.d0**starlum
c
            rstsun=dsqrt(starlum/(qht*fpi))
c
            write (*,180) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            dflum0=astar*dflum0
c
            write (*,140) blum,ilum,qht*astar
     &                   ,mcnuin,mcnuend,dflum0
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c    end source radius with photons
          endif
c    end blum>0
      endif
c
c     ***SETTING MONTE CARLO RT PARAMETER
c
ccc  200 format(//' Give initial energy packet number per range' /
  200 format(//' Give initial energy packet number per bin' /
     &         ' (log (<100) or (>100)) :'/
     &         ' increasing factor & convergence plateau : ',$)
      write (*,200)
      if (taskid.eq.0) read (*,*) npck, npckinc, convrt
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(npck,1,MPI_INTEGER4,0,MPI_COMM_WORLD,taskerr)
      call mpi_bcast(npckinc,1,MPI_INTEGER4,0,MPI_COMM_WORLD,taskerr)
      call mpi_bcast(convrt,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
      if (npck.lt.1.d2) npck=10**npck
      if (convrt.gt.1) convrt=convrt*0.01d0
c  
c     ***SETTING INNER RADIUS
c
  220 write (*,225)
  225 format(//' Give the inner radius (in cm) '/
     &         ' (log (<100) or (>100)) : ',$)
      if (taskid.eq.0) read (*,*) remp
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(remp,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
      if (remp.lt.1.0d2) remp=10.d0**remp 
      if (remp.lt.rstar) remp=rstar
c
c
  230 fren=0.1d0
      diend=0.0d0
c
  240 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &          '  Boundry Conditions  '/
     &          ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,240)
c
  250 format(/'  Choose a model ending : '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Radiation bounded, HII <',2pf5.2,'%'/
     &     '    E  :   Density bounded, Distance'/
     &     ' :: ',$)
      write (*,250) fren
      if (taskid.eq.0) read (*,90) jend
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(jend,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c
      if ((jend(1:1).eq.'A').or.(jend(1:1).eq.'a')) jend='A'
      if ((jend(1:1).eq.'E').or.(jend(1:1).eq.'e')) jend='E'
c
      if ((jend.ne.'A').and.(jend.ne.'E')) goto 230
c      
c      
      if (jend.eq.'E') then
  260    write (*,270)
  270    format(/' Give the distance or radius at which the density',
     &        ' drops (in cm) : ',$)
         if (taskid.eq.0) read (*,*) diend
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(diend,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
c         diend=remp+diend
         if (diend.le.remp) goto 260
      endif
c
      if (jend.eq.'A') then
         diend=-1.0d0
  280    write (*,290) fren
  290    format(/' Change the ionised hydrogen fraction',2pf5.2,'%?'/,
     &        ' (y/N) : ',$)
         if (taskid.eq.0) read (*,90) jren 
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(jren,4,MPI_CHARACTER,0,MPI_COMM_WORLD,taskerr)
c         
         if ((jren(1:1).eq.'Y').or.(jren(1:1).eq.'y')) then
  300       write (*,310) fren
  310       format(/' New ionised hydrogen fraction (0-1)',
     &           ' (currnt',2pf5.2,'%) : ',$)
            if (taskid.eq.0) read (*,*) fren
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
      call mpi_bcast(fren,1,MPI_REAL8,0,MPI_COMM_WORLD,taskerr)
c
            if ((fren.lt.0).or.(fren.gt.1)) goto 300
         endif
      endif
c
      difma=0.01d0
      dtlma=0.01d0
      dhlma=0.010d0
c
cc      stop
c
      call mciter (difma,dtlma,dhlma,blum,blum0)
c
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   
c     Start MC Iteration
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mciter (difma,dtlma,dhlma,blum,blum0)
c
      include 'cblocks.inc'
      include 'mpif.h'
c
c     variable
c
      integer*4 ix,iy,iz
      integer*4 iel,ion,j,i,inl
      integer*4 m,incream
c      integer*4 convg_arr(mcnx,mcny,mcnz)
      real*8 convg_arr(maxnx,maxny,maxnz)
      real*8 poparr(maxnx,maxny,maxnz,2)
      real*8 cov,cov0,tmp
      real*8 totion,totpix
      real*8 difma,dtlma,dhlma
      real*8 blum,blum0
      real*8 a,b
      real*8 opac(mxinfph),recpdf(mxinfph),abio
      real*8 tei,dei,dhni
      real*8 tef,def,dhnf    
      real*8 t,de,dh    
      real*8 tetmp,detmp,dhntmp    
      real*8 dton,dedhma,dedhmi
      real*8 dhl
      real*8 energ,wid
      real*8 pres0
      real*8 dv
c      real*8 poppre(mxion, mxelem)
      character nmod*4
c
      real*8 f5007,f4363,f3727,f6300,f6584,f6720
c
      real*8 sigmt,tauso,sigmtso,inttphot,inttphot2
      real*8 dustsigmat, dusttau
      real*8 difde,difte,dhhe,difg,difdh
      real*8 q1,q2,q3,q4
      real*8 q1field(maxnx,maxny,maxnz),q2field(maxnx,maxny,maxnz)
      real*8 q3field(maxnx,maxny,maxnz),q4field(maxnx,maxny,maxnz)
      real*8 teiarr(maxnx,maxny,maxnz),deiarr(maxnx,maxny,maxnz)
      real*8 dhiarr(maxnx,maxny,maxnz)
      real*8 rec0    
      real*8 rtmp, nsum1, nsum2    
c
      real*8 feldens,fpressu,dif
c
      integer*4 lulin,luop1,luop2
c
      character*512 outfile1, outfile2, outfile3
c
      dif(a,b)=dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
c
c
      dton=1.d33
      nmod='EQUI'
      lulin=5
c
      write (*,10)  
   10 format(///
     &     ' ****************************************************'/
     &     '   MCITER SUBROUTINE IS RUN'/
     &     ' ****************************************************')
c
      outfile1 = "./output/structure.out"
      outfile2 = "./output/all_spec.out"
      outfile3 = "./output/all_spec_int.out"
c         
      if (numtasks.gt.1) then
         write (outfile1,'("./output/structure_",i0,".out")') taskid
         write (outfile2,'("./output/all_spec_",i0,".out")') taskid
         write (outfile3,'("./output/all_spec_int_",i0,".out")') taskid
      endif
c
c
      m=0
      incream=1
      cov=0.d0
c
      pop=0.d0 
      do i=1,atypes
        pop(1,i)=1.d-5
        pop(maxion(i),i)=1.d0-1.d-5                  
      enddo  
c      
      do ix = 1, mcnx
      do iy = 1, mcny
      do iz = 1, mcnz         
c
cccccc
         te_arr(ix,iy,iz)=2000.d0  ! initial guess
         de_arr(ix,iy,iz)=0.0d0
         convg_arr(ix,iy,iz) = 0.0d0
cccccc
c
         dhni=dh_arr(ix,iy,iz)          
         tei=te_arr(ix,iy,iz)
         dei=de_arr(ix,iy,iz)
c
         poparr(ix,iy,iz,1) = pop(2,zmap(1))
         poparr(ix,iy,iz,2) = pop(2,zmap(2))
c
         recpdf=0.d0
         call calpdf(tei,dei,dhni,recpdf,rec0)  
c
         rec0_arr(ix,iy,iz)=rec0
c         
         do inl=1,infph-1
            jnu_arr(ix,iy,iz,inl) = 0.0d0
            recpdf_arr(ix,iy,iz,inl)=recpdf(inl) 
            tauso=0.d0
            sigmt=0.d0
            dustsigmat=0.d0
            call crosssections (inl, tauso, sigmt, dustsigmat)
            sigmt_arr(ix,iy,iz,inl)=sigmt*dhni
         enddo  
c
      enddo
      enddo
      enddo   
c   
      goto 40
c
c     Loop starts here
   20 continue
c
      cov0=cov
      cov=0.d0
      totion=0.d0
      totpix=0.d0
c         
      write (*,*) 'start'                
c
      nsum1=0.0d0
      nsum2=0.0d0
c
c
cccccccccccccccccccccccccccc
c
      if (m.eq.1) then
         open (lulin,file=trim(outfile2),
     &    status='unknown')
c
          do inl=1,infph-1
           write (lulin,*) photev(inl)
     &     ,jnu_arr(1,1,1,inl)
     &     ,recpdf_arr(1,1,1,inl)
     &     ,sigmt_arr(1,1,1,inl)
     &     ,te_arr(1,1,1),de_arr(1,1,1)
          enddo
          close (lulin)
c
      endif
c
      if (m.eq.-1) then
         open (lulin,file=trim(outfile2),
     &    status='unknown')
c
          do inl=1,infph-1
           write (lulin,*) photev(inl)
     &     ,jnu_arr(1,2,5,inl), jnu_arr(3,2,5,inl)
     &     ,jnu_arr(2,2,3,inl), jnu_arr(2,2,7,inl)
     &     ,recpdf_arr(1,2,5,inl),recpdf_arr(3,2,5,inl)
     &     ,recpdf_arr(2,2,3,inl),recpdf_arr(2,2,7,inl)
     &     ,sigmt_arr(1,2,5,inl),sigmt_arr(3,2,5,inl)
     &     ,sigmt_arr(2,2,3,inl),sigmt_arr(2,2,7,inl)
     &     ,te_arr(1,2,5),te_arr(3,2,5)
     &     ,te_arr(2,2,3),te_arr(2,2,7)
          enddo
          close (lulin)
c
          stop
      endif
cccccccccccccccccccccccccccccc
c
c
ccc      if (cov0.ge.0.7) then 
c
         open (lulin,file=trim(outfile1),status='unknown')
         write (lulin,*) m, '  x, y, z, Te, dne, Q1, Q2, 
     &                      REC0, DLOS, ELOSS, EGAIN'     
         close(lulin)

         call mcspec (luop1,1,'ZERO')

ccc      endif
c
c
      do ix = 1, mcnx
      do iy = 1, mcny
      do 30 iz = 1, mcnz        
c
         xps = ix
         yps = iy
         zps = iz
c
         dhni=dh_arr(ix,iy,iz)          
         tei=te_arr(ix,iy,iz)
         dei=de_arr(ix,iy,iz)  
c
         inttphot=0.d0
         inttphot2=0.d0
c
         do inl=1,infph-1
            tphot(inl)=0.d0
            wid=(photev(inl+1)-photev(inl))*evplk
            tphot(inl)=jnu_arr(ix,iy,iz,inl)/fpi/wid
            inttphot=inttphot+jnu_arr(ix,iy,iz,inl)
            if (photev(inl).ge.20) 
     &      inttphot2=inttphot2+jnu_arr(ix,iy,iz,inl)
            jnu_arr(ix,iy,iz,inl)=0.d0
         enddo   
c
         intpt = inttphot
c
         if (inttphot.le.epsilon) goto 25                
c
c         write (*,*) 'call equion, taskid', taskid
         call equion (tei, dei, dhni)
c         write (*,*) 'call equion end, taskid', taskid
c
         intpt = 1
c
c         write (*,*) 'call intvec, taskid', taskid
         call intvec (tphot, q1, q2, q3, q4)         
c         write (*,*) 'call intvec end, taskid', taskid
c         
         q1field(ix,iy,iz)=q1
         q2field(ix,iy,iz)=q2
         q3field(ix,iy,iz)=q3
         q4field(ix,iy,iz)=q4
c
         ipho=1
         iphom=0 ! ionising field changes
         nmod='EQUI'         
c   
c         write (*,*) 'call teequi2, taskid', taskid, nmod
         call teequi2 (tei, tef, def, dhni, 1.d33, nmod)
c         write (*,*) 'call teequi2, taskid', taskid
c
         dhnf=dhni
c         
c        --------------------
c        convergence criteria
c
         difte=dif(tef,tei)/dtlma
         dhhe=abs(1-poparr(ix,iy,iz,1)/pop(2,zmap(1)))
         dhhe=dmax1(abs(1-poparr(ix,iy,iz,2)/pop(2,zmap(2))),dhhe)
         dhhe=dhhe/dhlma
         difg=dmax1(dhhe,difte)
         dhhe=dhhe/difma
c
         if((q1.gt.epsilon)) totpix = totpix + 1.d0
c
         if ((difg.le.1).and.(pop(2,zmap(1)) .ge.fren)
     & .and.(q1.gt.epsilon)) cov=cov+1.d0
c
         if ((pop(2,zmap(1)).ge.fren)
     & .and.(q1.gt.epsilon)) totion=totion+1.d0
         convg_arr(ix,iy,iz)=difg
c        
         te_arr(ix,iy,iz)=tef
         de_arr(ix,iy,iz)=def
c         
   25    continue 
         recpdf=0.d0
         rec0=0.d0
         call calpdf(tef,def,dhnf,recpdf,rec0) 
         rec0_arr(ix,iy,iz)=rec0          
c
c        -------------
c        Calculate Tau
         sigmtso=0.0d0
         do inl=1,infph-1
            recpdf_arr(ix,iy,iz,inl)=recpdf(inl) 
            tauso=0.d0
            sigmt=0.d0
            dustsigmat=0.d0
            call crosssections (inl, tauso, sigmt, dustsigmat)
            sigmt_arr(ix,iy,iz,inl)=sigmt*dhnf
            sigmtso=sigmtso+sigmt*tphot(inl)/inttphot
         enddo 
c
         call cool (tef,def,dhnf)                  
c
c
ccc         if (cov0.ge.0.7) then
c
c
            if (ox3.ne.0) then
               hoiii(m)=fmbri(8,ox3)
               f5007=hoiii(m)
               f4363=fmbri(10,ox3)
            endif
            if (ox2.ne.0) then
               hoii(m)=(fmbri(1,ox2)+fmbri(2,ox2))
               f3727=hoii(m)
            endif
            if (ox1.ne.0) then
               hoi(m)=fmbri(3,ox1)
               f6300=hoi(m)
            endif
            if (ni2.ne.0) then
               hnii(m)=fmbri(10,ni2)
               f6584=hnii(m)
            endif
            if (su2.ne.0) then
               hsii(m)=(fmbri(1,su2)+fmbri(2,su2))
               f6720=hsii(m)
            endif
c
            open (lulin,file=trim(outfile1),
     &        status='OLD',access='APPEND') 
c
            write (lulin,*) (cwx(ix)+cwx(ix+1))*0.5d0,
     &                   (cwy(iy)+cwy(iy+1))*0.5d0,
     &                   (cwz(iz)+cwz(iz+1))*0.5d0,
     &                   te_arr(ix,iy,iz),
     &                   de_arr(ix,iy,iz),
     &                   q1field(ix,iy,iz),
     &                   q2field(ix,iy,iz),
     &                   q3field(ix,iy,iz),
     &                   q4field(ix,iy,iz),
     &                   inttphot,
     &                   rec0_arr(ix,iy,iz),
     &                   dh_arr(ix,iy,iz)
     &        ,hbeta,f5007,f4363,f3727,f6300,f6584,f6720
     &        ,((difg.le.1).and.(pop(2,zmap(1)).ge.fren)
     &          .and.(q1field(ix,iy,iz).gt.epsilon))
     &        ,((pop(2,zmap(1)).ge.fren)
     &          .and.(q1field(ix,iy,iz).gt.epsilon))
     &        ,(pop(j,zmap(1)),j=1,maxion(zmap(1)))
     &        ,(pop(j,zmap(2)),j=1,maxion(zmap(2)))
     &        ,(pop(j,zmap(8)),j=1,9)
     &        ,(pop(j,zmap(16)),j=1,9)
     &        ,poparr(ix,iy,iz,1)
     &        ,poparr(ix,iy,iz,2)
c
            close (lulin)     
c
            dv=((cwx(ix+1)-cwx(ix))
     &         *(cwy(iy+1)-cwy(iy))
     &         *(cwz(iz+1)-cwz(iz)))
c
            call mcspec (luop1,dv,'SUMM')
c            
ccc         endif
c
c
         poparr(ix,iy,iz,1) = pop(2,zmap(1))
         poparr(ix,iy,iz,2) = pop(2,zmap(2))
c
c         
   30 continue 
      enddo
      enddo 
c
c
      if (cov0.ge.0.7) then
c
   35       format(//' ======================================'//
     &               '     Lambda(A),  E (eV) , Flux (HB=1.0),',
     &               ' Species  , Kind, Accuracy (1-5)'/
     &               ' ====================================',
     &               '=============================')
c
           open (lulin,file=trim(outfile2),
     &          status='unknown') 
           write (lulin,35)
c
           call mcspec (lulin,1,'WRIT')
c
           close (lulin) 
c
      endif
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)  
      call mpi_allreduce(cov,tmp,1,
     &     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, taskerr)
c      
      cov = tmp
c
      call mpi_allreduce(totion,tmp,1,
     &     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, taskerr) 
c
      totion = tmp
c
      call mpi_allreduce(totpix,tmp,1,
     &     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, taskerr) 
c
      totpix = tmp
c
c
      if (totion.gt.0) then
         cov = cov / totion      
      else
         if (totion.lt.0) then 
            write (*,*) 'Wrong convergence. Taskid: ', taskid
         else
            cov = 0.d0
         endif
      endif
c
      if (totion.eq.0) then 
         npck = npck * npckinc
         goto 40
      endif
c
      if (cov.ge.convrt) goto 50
c
      if ((m.gt.1).and.(abs(cov-cov0).le.0.2*cov))
     &   npck = npck * npckinc
c
   40 m=m+1
c      
      if (taskid.eq.0) then
         if (m.eq.1)
     &   open (lulin,file='./output/status.out', 
     &        status='unknown')   
c
         if (m.gt.1)
     &   open (lulin,file='./output/status.out',
     &        status='OLD',access='APPEND')
         write (lulin,*) m, cov, totion, totpix, npck
         close(lulin) 
c
c
         if (m.eq.1)
     &   open (lulin,file='./output/grid.out', 
     &        status='unknown')
         write (lulin,*) mcnx
         do i=1,mcnx+1
            write (lulin,*) i, cwx(i)
         enddo  
         write (lulin,*) mcny
         do i=1,mcny+1
            write (lulin,*) i, cwy(i)
         enddo 
         write (lulin,*) mcnz
         do i=1,mcnz+1
            write (lulin,*) i, cwz(i)
         enddo 
         close(lulin)
c
c
      endif        
c
      if (m.le.20) then
c
         call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
         call mcrt0 (m,incream,blum,blum0)   
c
         goto 20  
c
      endif
c
   50 return
c      
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to transfer energy packet
c
c            OUTPUT: Jnu(?)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mcrt0 (m,incream,blum,blum0)
c
      include 'cblocks.inc'
      include 'mpif.h'
c    
      logical   lgstr,lgdif,lgescape,lgabs
      logical   lgfin,lgrec,lgactive
      integer*4 ix,iy,iz,icont
      integer*4 absEvent
      integer*4 iphot,nphot,nphtlop,iloop
      integer*4 cwi(3),ci(3)
      integer*4 nup,nup0,nupcnt,m
      integer*4 xp, yp, zp
      integer*4 numstar,irg
      integer*4 incream,npck0
c      
      real*8    blum,blum0,wid,invnph,invblum
      real*8    rvec(3),rvectmp(3),hvec(3)
      real*8    drvec,drvectmp
      real*8    nust,nued,mcde
      real*8    mcdx,mcdy,mcdz
c
      real*8    abstau, tauloc
      real*8    passprob, random, ran2
c
      real*8    cw
      real*8    rx,ry,rz,rx0,ry0,rz0
c
      real*8    ds,dsx,dsy,dsz
      real*8    dv,vllmt,volfil
c
      real*8    a,b,c
      real*8    posorg(3)
c
      real*8    sigmt 
c
      real*8    epckarr(mxinfph,looplmt,7)
      integer*4 nuparr(mxinfph,looplmt)
      integer*4 absarr(mxinfph,looplmt)
      logical   flgarr(mxinfph,looplmt,4)
      real*8    cwmax(3),cwmin(3)
c
      real*8    epckarr0(mxinfph,looplmt,7)
      integer*4 nuparr0(mxinfph,looplmt)
      integer*4 absarr0(mxinfph,looplmt)
      logical   flgarr0(mxinfph,looplmt,4)
      real*8    cwmax0(3),cwmin0(3)
c
      integer*4 nfinevnt,nfinevnt0,ncnt,nmpi
      integer*4 nphot0
      logical   lgtag,lgdebug,lgsilence
c
c
cccccccccccccc    
      integer*4 lulin      
cccccccccccccc    
c
c
      lgdebug = .false.
      lgsilence = .true.
c
c         
      if (m.eq.1) write (*,10) 
   10 format(///
     &     ' ****************************************************'/
     &     '   MC-RT SUBROUTINE IS RUN  ')
      write (*,*) '   Loop: ', m
c
c     For radiation bounded case, there is no boundary limit      
      if (diend.eq.-1.0d0) 
     &    diend=dsqrt((cwx(mcnx+1)-cwx(1))**2.0d0
     &               +(cwy(mcny+1)-cwy(1))**2.0d0
     &               +(cwz(mcnz+1)-cwz(1))**2.0d0)
c      
cccccccccccccc    
      strpos(1) = 0.0d0
      strpos(2) = 0.0d0
      strpos(3) = 0.0d0
c      strpos(3) = -9.0000001e+18
      numstar = 1
cccccccccccccc    
c
      npck0=npck
c
      invnph=1.d0/npck 
      invblum=1.d0/blum0
c
c      nphot = npck0
c
      if (npck0.gt.looplmt) then
         nphot = looplmt
         nphtlop = npck0/nphot
         if (nphtlop*nphot.lt.npck0) nphtlop=nphtlop+1
      else
         nphot = npck0
         nphtlop = 1
      endif
c
      if (.not.lgsilence) then
         if (taskid.eq.0) then
            write (*,*) 'Nphot, Npack0, Npack'
            write (*,*) nphot, npck0, npck            
         endif
      endif
c
ccccccccccccc
c
c
      cwmin0(1) = cwx(1)
      cwmax0(1) = cwx(mcnx+1)
c
      cwmin0(2) = cwy(1)
      cwmax0(2) = cwy(mcny+1)      
c
      cwmin0(3) = cwz(1)
      cwmax0(3) = cwz(mcnz+1)
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)  
      call mpi_allreduce(cwmax0,cwmax,3,
     &     MPI_REAL8,MPI_MAX,MPI_COMM_WORLD, taskerr) 
      call mpi_allreduce(cwmin0,cwmin,3,
     &     MPI_REAL8,MPI_MIN,MPI_COMM_WORLD, taskerr)       
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      if (.not.lgsilence) then
         write (*,*) taskid,'max',
     &             cwmax(1), cwmax(2), cwmax(3),
     &             cwx(mcnx+1), cwy(mcny+1), cwz(mcnz+1)
         write (*,*) taskid,'min',  
     &             cwmin(1), cwmin(2), cwmin(3),
     &             cwx(1), cwy(1), cwz(1)
         call mpi_barrier(MPI_COMM_WORLD, taskerr)
      endif 
c
ccccccccccccc
c
c
      flgarr0 = .false.
      epckarr0 = 0                                   
      nuparr0 = 0
      absarr0 = 0
c
c
ccccccccccccc
c
      irg = 0
   15 irg = irg + 1
c
      if (irg.eq.nphtlop) then
         if (irg*looplmt.gt.npck0) then
            nphot0 = nphot
            nphot = npck0 - (nphtlop-1)*nphot0              
         endif
      endif
c
      if (taskid.eq.0)
     & write (*,*) 'round:',irg,'/',nphtlop
c
ccccccccccccc
c
c  initialize
c
      do nupcnt = 1, infph-1
         do iphot = 1, nphot    
c           lgstr = .true.          
            flgarr(nupcnt,iphot,1) = .true.  
c           lgdif = .false.
            flgarr(nupcnt,iphot,2) = .false.  
c           lgfin = .false.      
            flgarr(nupcnt,iphot,3) = .false. 
c           lgrec = .false.      
            flgarr(nupcnt,iphot,4) = .false.
c           position
            epckarr(nupcnt,iphot,1) = strpos(1)      
            epckarr(nupcnt,iphot,2) = strpos(2)      
            epckarr(nupcnt,iphot,3) = strpos(3)      
c           direction
            hvec(1)=0.0d0
            hvec(2)=0.0d0
            hvec(3)=0.0d0                  
            call photdir (hvec)
            epckarr(nupcnt,iphot,4) = hvec(1)      
            epckarr(nupcnt,iphot,5) = hvec(2)      
            epckarr(nupcnt,iphot,6) = hvec(3) 
c           Tau              
            random=0.0d0
            random=ran2(seed)            
            epckarr(nupcnt,iphot,7) = -log(1.0d0-random)                                   
c           frequency
            nuparr(nupcnt,iphot) = nupcnt
c           absevent
            absarr(nupcnt,iphot) = 0            
         enddo
      enddo
c
      nmpi = 0
c
   20 nfinevnt = 0
c
      if (.not.lgsilence) 
     &write (*,*) 'Taskid:',taskid, ' Enter RT, nmpi=', nmpi, nphot
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)   
c
      do nupcnt = 1, infph-1
         do iphot = 1, nphot     
ccc      do nupcnt = 5270, 5270
ccc         do iphot = 1, 1
c
            lgescape = .true.                 
c
            nup = nuparr(nupcnt,iphot)
c
            if (nup.eq.0) then
               lgactive = .false.
c               lgescape = .true.
               goto 50
            endif            
c
c            if ((photev(nup).lt.mcnuin).or.
c     &          (photev(nup+1).gt.mcnuend)) then
c               lgactive = .false.
c               lgescape = .true.
c               goto 50
c            endif
c      
            lgtag = .false.
            if ((nup.le.0).or.(nup.ge.mxinfph-1)) then
               lgtag = .true.
            endif
c
c      call mpi_barrier(MPI_COMM_WORLD, taskerr)          
c      call mpi_allreduce(lgtag,lgtag,1,
c     &     MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD, taskerr)
c
            if (lgtag) then
               if (.not.lgdebug)
     &         write (*,*) 'taskid',taskid,nup
     &         , nuparr(nupcnt,iphot), nuparr0(nupcnt,iphot)
     &         , nupcnt, iphot
     &         , 'lgstr', flgarr(nupcnt,iphot,1)
     &         , 'lgdif', flgarr(nupcnt,iphot,2)
     &         , 'lgfin', flgarr(nupcnt,iphot,3)
     &         , 'lgrec', flgarr(nupcnt,iphot,4)
c                call mpi_barrier(MPI_COMM_WORLD, taskerr)          
               stop
            endif
c
c
            lgstr = flgarr(nupcnt,iphot,1)  
            lgdif = flgarr(nupcnt,iphot,2)  
            lgfin = flgarr(nupcnt,iphot,3) 
            lgrec = flgarr(nupcnt,iphot,4)       
c
            absEvent = absarr(nupcnt,iphot)
            passprob = epckarr(nupcnt,iphot,7)  
c
            if (lgfin) goto 50   
c
c            
c            wid=photev(nup+1)-photev(nup)
c            mcde=soupho(nup)*wid*evplk*invnph*blum*invblum
            wid=photev(nupcnt+1)-photev(nupcnt)
            mcde=soupho(nupcnt)*wid*evplk*invnph*blum*invblum
c
            lgactive=.true.
c      
c           initialise photon position
            posorg(1)=epckarr(nupcnt,iphot,1)      
            posorg(2)=epckarr(nupcnt,iphot,2)
            posorg(3)=epckarr(nupcnt,iphot,3)
            call gtcwi(posorg,cwi)     
            xp=cwi(1)
            yp=cwi(2)
            zp=cwi(3)         
c
c           initialise photon direction
            hvec(1)=epckarr(nupcnt,iphot,4)      
            hvec(2)=epckarr(nupcnt,iphot,5)
            hvec(3)=epckarr(nupcnt,iphot,6)            
c
c
c            hvec(1)=-0.26842783623701466      
c            hvec(2)=-0.87942114561703666       
c            hvec(3)=0.39314748552513135
c
c           determine which block
            if ((xp*yp*zp).eq.0) then
               lgactive = .false.
c               lgescape = .true.
            endif
               
            if (((dabs(posorg(1)-cwx(1)).le.1d-10) 
     &      .and.(hvec(1).lt.-1d-10)) 
     &      .or.((dabs(posorg(1)-cwx(mcnx+1)).le.1d-10)
     &      .and.(hvec(1).gt.1d-10))) then
               lgactive = .false.
c               lgescape = .true.               
            endif  
c
            if (((dabs(posorg(2)-cwy(1)).le.1d-10) 
     &      .and.(hvec(2).lt.-1d-10)) 
     &      .or.((dabs(posorg(2)-cwy(mcny+1)).le.1d-10)
     &      .and.(hvec(2).gt.1d-10))) then
               lgactive = .false.
c               lgescape = .true.               
            endif  
c
            if (((dabs(posorg(3)-cwz(1)).le.1d-10) 
     &      .and.(hvec(3).lt.-1d-10)) 
     &      .or.((dabs(posorg(3)-cwz(mcnz+1)).le.1d-10)
     &      .and.(hvec(3).gt.1d-10))) then
               lgactive = .false.
c               lgescape = .true.               
            endif
c
            if (.not.lgactive) goto 50              
c
c          
c
            lgescape = .false.            
c
   30       if ((.not.lgrec).and.(.not.lgescape)
     &     .and.(absEvent.le.safelmt)) then 
c
c              update photon position
               rvec(1)=posorg(1)
               rvec(2)=posorg(2)
               rvec(3)=posorg(3)
c
c              Initialise photon Nu & Direction
               if (.not.(lgstr.or.lgdif)) then
                  write(*,*)'Neither Stellar nor Diffuse Emission'                        
                  stop
               else
c                 Initialise photon cell index                
                  ci(1)=xp
                  ci(2)=yp
                  ci(3)=zp                       
c                                        
                  lgabs=.false.
c
               endif
c               
c              distance to center star
               drvec=0.0d0
               drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                    +(rvec(2)-strpos(2))**2.0d0
     &                    +(rvec(3)-strpos(3))**2.0d0)
c
             if ((photev(nup).lt.mcnuin)
     &       .or.(photev(nup+1).gt.mcnuend)               
     &       .or.(drvec.gt.diend)) then 
                  lgescape = .true.
                  lgfin = .true.
                  goto 50
               endif
c 
c              Start path loop
c              Looking for dSx & determining Cell x-index  a
   40          if (dabs(hvec(1)).le.1d-10) then
                  dsx = cwx(mcnx+1)
               else 
                  if (hvec(1).gt.1d-10) then
                     dsx = (cwx(ci(1)+1)-rvec(1))/hvec(1)
                  else
                     if (hvec(1).lt.-1d-10) then
                        dsx = (cwx(ci(1))-rvec(1))/hvec(1)
                     else
                        write(*,*)'Wrong Dirction Vector in X'                    
                        stop
                     endif
                  endif
               endif                
c
c              Looking for dSy & determining Cell y-index             
               if (dabs(hvec(2)).le.1d-10) then
                  dsy = cwy(mcny+1)
               else 
                  if (hvec(2).gt.1d-10) then
                     dsy = (cwy(ci(2)+1)-rvec(2))/hvec(2)
                  else
                     if (hvec(2).lt.-1d-10) then
                        dsy = (cwy(ci(2))-rvec(2))/hvec(2)
                     else
                        write(*,*)'Wrong Dirction Vector in Y'                    
                        stop
                     endif
                  endif   
               endif
c
c              Looking for dSz & determining Cell z-index             
               if (dabs(hvec(3)).le.1d-10) then
                  dsz = cwz(mcnz+1)
               else 
                  if (hvec(3).gt.1d-10) then
                     dsz = (cwz(ci(3)+1)-rvec(3))/hvec(3)
                  else
                     if (hvec(3).lt.-1d-10) then
                        dsz = (cwz(ci(3))-rvec(3))/hvec(3)
                     else
                        write(*,*)'Wrong Dirction Vector in Z'                    
                        stop
                     endif
                  endif   
               endif
c
               if (lgdebug) then
               write (*,*) iphot,nup,rvec(1),rvec(2),rvec(3),
     &                     hvec(1),hvec(2),hvec(3),
     &                     lgstr,lgdif,passprob,taskid
               endif
c
c
               ds=0.0d0
               ds=min(dsx,dsy,dsz)
               dv=0.0d0
               dv=((cwx(ci(1)+1)-cwx(ci(1)))
     &            *(cwy(ci(2)+1)-cwy(ci(2)))
     &            *(cwz(ci(3)+1)-cwz(ci(3))))  
c       
               rvectmp(1)=0.0d0
               rvectmp(2)=0.0d0
               rvectmp(3)=0.0d0
               rvectmp(1)=rvec(1)+ds*hvec(1)
               rvectmp(2)=rvec(2)+ds*hvec(2)
               rvectmp(3)=rvec(3)+ds*hvec(3)
c       
               drvectmp=0.0d0
               drvectmp=dsqrt((rvectmp(1)-strpos(1))**2.0d0
     &                       +(rvectmp(2)-strpos(2))**2.0d0
     &                       +(rvectmp(3)-strpos(3))**2.0d0)
c
c
c              Whether photon reach Rin (remp) or useless cell
               if ((drvectmp.le.remp).and.(drvec.lt.remp)
     &            .and.(drvectmp.gt.rstar)) then
c
                  if ((dabs(ds-dsx).le.1d-10)
     &               .and.(hvec(1).gt.1d-10)) then
                      xp=xp+1 
                  else 
                  if ((dabs(ds-dsx).le.1d-10)
     &               .and.(hvec(1).lt.-1d-10)) then
                      xp=xp-1
                  else    
                  if ((dabs(ds-dsy).le.1d-10)
     &               .and.(hvec(2).gt.1d-10)) then
                      yp=yp+1 
                  else
                  if ((dabs(ds-dsy).le.1d-10)
     &               .and.(hvec(2).lt.-1d-10)) then
                      yp=yp-1
                  else    
                  if ((dabs(ds-dsz).le.1d-10)
     &               .and.(hvec(3).gt.1d-10)) then
                      zp=zp+1 
                  else    
                  if ((dabs(ds-dsz).le.1d-10)
     &               .and.(hvec(3).lt.-1d-10)) then
                      zp=zp-1
                  endif 
                  endif 
                  endif 
                  endif 
                  endif
                  endif                
c
                  ci(1)=xp
                  ci(2)=yp
                  ci(3)=zp 
c                  
                  rvec(1)=rvectmp(1)
                  rvec(2)=rvectmp(2)
                  rvec(3)=rvectmp(3) 
c 
                  if ((xp.le.0)
     &            .or.(yp.le.0)
     &            .or.(zp.le.0)
     &            .or.(xp.ge.(mcnx+1))
     &            .or.(yp.ge.(mcny+1))
     &            .or.(zp.ge.(mcnz+1))) then
                     lgescape = .true.
                     goto 50
                  endif             
c
                  goto 40        
               endif
c                
c              Whether photon reach the surface of star (rstar)
               if (drvectmp.le.rstar) then
                  lgescape = .true.
                  lgfin = .true.
                  goto 50               
               endif
c
c              Whether photon out of Rout (diend)
               if (drvectmp.ge.diend) then
                  a=1.0d0
                  b=2.0d0*(rvec(1)*hvec(1)
     &                    +rvec(2)*hvec(2)
     &                    +rvec(3)*hvec(3))
                  c=(drvec**2.0d0-diend**2.0d0)
                  ds=(dsqrt(b**2.0d0-4.0d0*a*c)-b)/(2.0d0*a)
               endif     
c          
c
c              Calculate local Tau
               sigmt=0.d0 
               tauloc=0.d0      
               sigmt=sigmt_arr(ci(1),ci(2),ci(3),nup)       
               if (sigmt.gt.0)
     &         tauloc=sigmt*ds
c
c               if (sigmt.gt.0)
c     &         write (*,*) iphot,nup,passprob,tauloc,taskid,
c     &         m
c
c              Whether the photon is absorbed?
               if ((tauloc.ge.passprob).and.(tauloc.gt.0.0d0))
     &         then                    
                  ds=passprob/sigmt
                  rvec(1)=rvec(1)+ds*hvec(1)
                  rvec(2)=rvec(2)+ds*hvec(2)
                  rvec(3)=rvec(3)+ds*hvec(3)
c                
                  jnu_arr(ci(1),ci(2),ci(3),nup)
     &           =jnu_arr(ci(1),ci(2),ci(3),nup)+ds*mcde/dv
c
c
c
                  if (isnan(jnu_arr(ci(1),ci(2),ci(3),nup)))
     &            then
                     write (*,*) 'abs', ds,ci(1),ci(2),ci(3)
     &                           ,ds*mcde/dv,dv,taskid
                     stop
                  endif
c
c
c
                  lgstr=.false.
                  lgdif=.true.
c
                  drvec=0.0d0
                  drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                       +(rvec(2)-strpos(2))**2.0d0
     &                       +(rvec(3)-strpos(3))**2.0d0)                  
c
                  if (drvec.ge.diend) then
                     lgescape=.true.
                     lgfin = .true.
                     goto 50
                  else                    
                     lgabs=.true.
c
                     if ((rvec(1).le.cwmin(1)).or.
     &                   (rvec(1).ge.cwmax(1)).or.
     &                   (rvec(2).le.cwmin(2)).or.
     &                   (rvec(2).ge.cwmax(2)).or.
     &                   (rvec(3).le.cwmin(3)).or.
     &                   (rvec(3).ge.cwmax(3))) then
                            lgescape=.true.
                            lgfin = .true.
                            goto 50            
                     endif
c
                     if ((rvec(1).lt.cwx(1)).or.
     &                   (rvec(1).gt.cwx(mcnx+1)).or.
     &                   (rvec(2).lt.cwy(1)).or.
     &                   (rvec(2).gt.cwy(mcny+1)).or.
     &                   (rvec(3).lt.cwz(1)).or.
     &                   (rvec(3).gt.cwz(mcnz+1))) then
                            lgescape=.true.
                            goto 50
                     endif
                  endif
c               
               else
c
                  if (dabs(ds).ge.1.0d-10) then
c
                     rvec(1)=rvec(1)+ds*hvec(1)
                     rvec(2)=rvec(2)+ds*hvec(2)
                     rvec(3)=rvec(3)+ds*hvec(3)
c
                     passprob=passprob-tauloc        
c
ccc                  if ((ci(1).eq.2).and.(ci(2).eq.13).and.(ci(3).eq.13)
ccc     &               ) then
ccc                   if (lgstr) jnuw_arr(2,nup)=jnuw_arr(2,nup)+ds*mcde/dv 
ccc                   if (lgdif) jnuw_arr(1,nup)=jnuw_arr(1,nup)+ds*mcde/dv 
ccc                  endif  
c  
c                  write (*,*) ci                
                     jnu_arr(ci(1),ci(2),ci(3),nup)
     &                =jnu_arr(ci(1),ci(2),ci(3),nup)+ds*mcde/dv        
c
                  endif
c
c
                  if (isnan(jnu_arr(ci(1),ci(2),ci(3),nup)))
     &            then
                     write (*,*) 'noabs', ds,ci(1),ci(2),ci(3)
     &                           ,ds*mcde/dv,dv,taskid
                     stop
                  endif
c
c
                  drvec=0.0d0
                  drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                       +(rvec(2)-strpos(2))**2.0d0
     &                       +(rvec(3)-strpos(3))**2.0d0)    
c
                  if ((drvec.ge.diend).or.
     &                (rvec(1).le.cwmin(1)).or.
     &                (rvec(1).ge.cwmax(1)).or.
     &                (rvec(2).le.cwmin(2)).or.
     &                (rvec(2).ge.cwmax(2)).or.
     &                (rvec(3).le.cwmin(3)).or.
     &                (rvec(3).ge.cwmax(3))) then
                     lgescape=.true.
                     lgfin = .true.
c
                     goto 50            
                  endif
c
                  if ((rvec(1).lt.cwx(1)).or.
     &                (rvec(1).gt.cwx(mcnx+1)).or.
     &                (rvec(2).lt.cwy(1)).or.
     &                (rvec(2).gt.cwy(mcny+1)).or.
     &                (rvec(3).lt.cwz(1)).or.
     &                (rvec(3).gt.cwz(mcnz+1))) then
c                    
                     lgescape=.true.
c
                     goto 50
c
                  else
c
                     if ((dabs(ds-dsx).le.1d-10)
     &                  .and.(hvec(1).gt.1d-10)) then
                        xp=xp+1 
                     else 
                     if ((dabs(ds-dsx).le.1d-10)
     &                  .and.(hvec(1).lt.-1d-10)) then
                        xp=xp-1
                     else    
                     if ((dabs(ds-dsy).le.1d-10)
     &                  .and.(hvec(2).gt.1d-10)) then
                        yp=yp+1 
                     else
                     if ((dabs(ds-dsy).le.1d-10)
     &                  .and.(hvec(2).lt.-1d-10)) then
                        yp=yp-1
                     else    
                     if ((dabs(ds-dsz).le.1d-10)
     &                  .and.(hvec(3).gt.1d-10)) then
                        zp=zp+1 
                     else    
                     if ((dabs(ds-dsz).le.1d-10)
     &                  .and.(hvec(3).lt.-1d-10)) then
                        zp=zp-1
                     endif 
                     endif 
                     endif 
                     endif 
                     endif
                     endif                   
c 
                     if ((xp.le.0)
     &               .or.(yp.le.0)
     &               .or.(zp.le.0)
     &               .or.(xp.ge.(mcnx+1))
     &               .or.(yp.ge.(mcny+1))
     &               .or.(zp.ge.(mcnz+1))) then
                        lgescape = .true.
c
                        goto 50
                     endif  
c
                     ci(1)=xp
                     ci(2)=yp
                     ci(3)=zp 
c
                     goto 40  
                  endif 
               endif
c
               if (lgabs) then
                  if(lgescape) then 
                    write(*,*)'Wrong in photon travelling'                    
                    stop  
                  endif              
                  absEvent=absEvent+1
                  posorg(1)=rvec(1)
                  posorg(2)=rvec(2)
                  posorg(3)=rvec(3)
c
c                 Initialise photon Direction                
                  hvec(1)=0.0d0
                  hvec(2)=0.0d0
                  hvec(3)=0.0d0                  
                  call photdir (hvec)
c                  
c                 Initialise photon Frequency                
c
c                  seed = time() 
                  random=ran2(seed)  
                  if (random.lt.rec0_arr(ci(1),ci(2),ci(3))) then
                      lgescape = .true.
                      lgrec = .true.
                      lgfin = .true.
                      goto 50                        
                  endif                     
                  call frqgen (nup,mcnuin,mcnuend,
     &                         ci,lgstr,lgdif,numstar)
c
c
c                 set random number for absorption event
                  random=0.0d0
                  random=ran2(seed)          
                  passprob = -log(1.0d0-random)
c
                  goto 30 
c                  
               else 
                  if(lgescape) then 
                    goto 50
                  else
                    write(*,*)'Wrong in photon travelling'                    
                    stop  
                  endif
               endif                                 
c
c           end label 30
            endif      
c
            if (absEvent.gt.safelmt) then
               lgfin = .true.
               lgescape = .true.
            endif
c
   50       continue      
c   
            if (lgescape.and.lgactive) then
c           lgstr = .true.          
               flgarr0(nupcnt,iphot,1) = lgstr  
c           lgdif = .false.
               flgarr0(nupcnt,iphot,2) = lgdif  
c           lgfin = .false.      
               flgarr0(nupcnt,iphot,3) = lgfin 
c           lgrec = .false.      
               flgarr0(nupcnt,iphot,4) = lgrec
c           position
               epckarr0(nupcnt,iphot,1) = rvec(1)      
               epckarr0(nupcnt,iphot,2) = rvec(2)      
               epckarr0(nupcnt,iphot,3) = rvec(3)      
c           direction
               epckarr0(nupcnt,iphot,4) = hvec(1)      
               epckarr0(nupcnt,iphot,5) = hvec(2)      
               epckarr0(nupcnt,iphot,6) = hvec(3) 
c           Tau        
               epckarr0(nupcnt,iphot,7) = passprob                                   
c           frequency
               nuparr0(nupcnt,iphot) = nup
c           absevent
               absarr0(nupcnt,iphot) = absEvent
c
               if (lgfin) nfinevnt = nfinevnt + 1
c               
            else
               if (lgescape.and.(.not.lgactive)) then
                  flgarr0(nupcnt,iphot,1) = .false.  
                  flgarr0(nupcnt,iphot,2) = .false.
                  flgarr0(nupcnt,iphot,3) = .false.
                  flgarr0(nupcnt,iphot,4) = .false.
                  epckarr0(nupcnt,iphot,1) = 0      
                  epckarr0(nupcnt,iphot,2) = 0      
                  epckarr0(nupcnt,iphot,3) = 0      
                  epckarr0(nupcnt,iphot,4) = 0      
                  epckarr0(nupcnt,iphot,5) = 0     
                  epckarr0(nupcnt,iphot,6) = 0 
                  epckarr0(nupcnt,iphot,7) = 0                                   
                  nuparr0(nupcnt,iphot) = 0
                  absarr0(nupcnt,iphot) = 0
               else
                  write (*,*) 'Wrong process, taskid: ', taskid
                  write (*,*) lgstr,lgdif,lgfin,lgrec
                  write (*,*) rvec(1),rvec(2),rvec(3)      
                  write (*,*) hvec(1),hvec(2),hvec(3) 
                  write (*,*) passprob,nup,absEvent  
                  write (*,*) absEvent, safelmt               
                  stop   
               endif
            endif
c
         enddo
      enddo   
c
c     MPI message interchange
c
      if (.not.lgsilence)
     & write (*,*) 'mpi start, taskid:', taskid, ' nmpi:',nmpi
      call mpi_barrier(MPI_COMM_WORLD, taskerr)  
c          
      ncnt = mxinfph*looplmt*4      
      call mpi_allreduce(flgarr0,flgarr,ncnt,
     &     MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD, taskerr)
c
      if (.not.lgsilence)
     &write (*,*) 'flg, taskid=', taskid
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      ncnt = mxinfph*looplmt*7
c      
      if (.not.lgsilence)
     &write (*,*) taskid, ncnt, mxinfph, looplmt
c
      call mpi_allreduce(epckarr0,epckarr,ncnt,
     &     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, taskerr)      
c
      if (.not.lgsilence)
     &write (*,*) 'epckarr, taskid=', taskid
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      ncnt = mxinfph*looplmt
      call mpi_allreduce(nuparr0,nuparr,ncnt,
     &     MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD, taskerr) 
c
      if (.not.lgsilence)
     &write (*,*) 'nup, taskid=', taskid
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      ncnt = mxinfph*looplmt
      call mpi_allreduce(absarr0,absarr,ncnt,
     &     MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD, taskerr)      
c
      if (.not.lgsilence)
     &write (*,*) 'abs, taskid=', taskid
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      call mpi_allreduce(nfinevnt,nfinevnt0,1,
     &     MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD, taskerr)
c
      nfinevnt = nfinevnt0
c
      if (.not.lgsilence)
     &write (*,*) 'fin, taskid=', taskid
c
      call mpi_barrier(MPI_COMM_WORLD, taskerr)
c
      if (taskid.eq.0) then
         write (*,*) 'Nfin,  Nphot @ nmpi:', nmpi
         write (*,*) nfinevnt, (infph-1)*nphot
      endif
c
      if (.not.lgsilence)
     &write (*,*) 'Taskid:',taskid, ' Leave RT, nmpi=', nmpi
c
      nmpi = nmpi + 1
c
      if ((nfinevnt.lt.(infph-1)*nphot).and.(nmpi.le.safelmt)) 
     & goto 20
c
      if (.not.lgsilence)
     &write (*,*) 'Taskid:',taskid, ' End MPI'
c
      if (irg.lt.nphtlop) goto 15
c
c      stop
c
      return
c
      end     
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get Cell Wall Index CWI(6) 
c    
c     INPUT:  rvec -- photon position vector 
c          
c     OUTPUT: cwi  -- cell wall index
c             CWI(1),CWI(2): Left & Right wall in X
c             CWI(3),CWI(4): Left & Right wall in Y
c             CWI(5),CWI(6): Left & Right wall in Z
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gtcwi(rvec,cwi)
c
      include 'cblocks.inc'
c
      real*8    rvec(3)
      integer*4 cwi(3)  
      integer*4 i
c
c     External Function
      real*8    cw
c
c
      cwi=0
c
      i=1   
   10 if (rvec(1).lt.cwx(i)) then 
         goto 20
c         write(*,*)'Wrong in photon position-X'                    
c         stop
      endif  
      if (rvec(1).lt.cwx(i+1)) then 
         cwi(1)=i                          
      else
         if (i.lt.mcnx) then
            i=i+1
            goto 10
         else 
            if (rvec(1).gt.(cwx(mcnx+1))) then
               cwi(1)=0
            else
               cwi(1)=mcnx+1
            endif
c            write(*,*)'Cannot find cell wall in X'                    
c            stop             
         endif  
      endif         
c
c
      i=1   
   20 if (rvec(2).lt.cwy(i)) then 
         goto 30
c         write(*,*)'Wrong in photon position-Y'                    
c         stop
      endif  
      if (rvec(2).lt.cwy(i+1)) then 
         cwi(2)=i                          
      else
         if (i.lt.mcny) then
            i=i+1
            goto 20
         else 
            if (rvec(2).gt.(cwy(mcny+1))) then          
               cwi(2)=0
            else
               cwi(2)=mcny+1
            endif
c            write(*,*)'Cannot find cell wall in Y'                    
c            stop             
         endif  
      endif   
c
c
      i=1   
   30 if (rvec(3).lt.cwz(i)) then 
         goto 40
c         write(*,*)'Wrong in photon position-Z'                    
c         stop
      endif  
      if (rvec(3).lt.cwz(i+1)) then 
         cwi(3)=i                          
      else
         if (i.lt.mcnz) then
            i=i+1
            goto 30
         else 
            if (rvec(3).gt.(cwz(mcnz+1))) then
               cwi(3)=0
            else
               cwi(3)=mcnz+1
            endif
c            write(*,*)'Cannot find cell wall in Z'                    
c            stop             
         endif  
      endif  
c  
   40 continue
      return
c
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Create Photon Direction
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine photdir (hvec)
c
      integer*4   i
      real*8      random, ran2
      real*8      hvec(3) 
      real*8      w,u,v,t,ang
c
c      
      random = ran2(seed)
      w = 2.*random - 1.
      t = sqrt(1.-w*w)
      random = ran2(seed)
      ang = 3.141592654*(2.*random-1.)
      u = t*cos(ang)
      v = t*sin(ang)
      hvec(1) = u
      hvec(2) = v
      hvec(3) = w   
c
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Frequency Generator
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
      subroutine frqgen (nup,nust,nued,ci,lgstr,lgdif,numstar)
c
      include 'cblocks.inc'
c
      logical   lgstr, lgdif
      integer*4 nup, i, numstar
      integer*4 ci(3)
c
      real*8    random, ran2      
      real*8    nust, nued
c
c      
   10 format(//
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Problem with Frequency Generator'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')     
c
      i = 0
   20 random = ran2(seed)
      i = i + 1
      if (i.ge.1000) then
         write (*,10) '1', i
         stop
      endif
      if ((random.ge.1.0d0).or.(random.le.0.0d0)) goto 20
c 
      nup=0     
   30 nup=nup+1
c
      if (photev(nup).lt.nust) goto 30
c      
      if (photev(nup).ge.nued) then
         write (*,*) '2', photev(nup), nued, nust        
         write (*,10) 
         stop
      endif   
c
      if ((lgstr.and.lgdif).or.((.not.lgstr).and.(.not.lgdif))) then
         write (*,10) '3', lgstr, lgdif
         stop      
      endif              
c     stellar photon      
      if (lgstr.and.(.not.lgdif)
     &   .and.(random.ge.strpdf_arr(numstar,nup))) goto 30
c     diffuse photon      
      if ((.not.lgstr).and.lgdif
     &   .and.(random.ge.recpdf_arr(ci(1),ci(2),ci(3),nup))) goto 30     
c
c
      return
      end     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Random Number Generator ran2 in Numerical Recipe
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function ran2(idum)
c
      implicit none
c
      integer*4 idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*8 ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1.0d0/im1,imm1=im1-1,
     &           ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &           ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2d-7,
     &           rnmx=1.0d0-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum=idum+im1
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif 
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if (iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
c
      return
c      
      end 
c    
