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
c
c           Variables
c
      real*8 blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      real*8 epotmi,fin,uinit, qinit, rcin
      real*8 rstsun,wid,starlum
      real*8 scale,dnt,lam
      real*8 blum0,dflum0,dflum(8)
      real*8 mcdx,mcdy,mcdz
c      real*8 fluxes(mxmonlines)
c
      integer*4 lenv      
      integer*4 l,m,n
c
      integer*4 j,luin,i,kmax,k,nl
      integer*4 idx1 ,idx2, idx3, idx4, iem      
      integer*4 numstar
c            
c      character=jtb*4
      character carac*36,model*32
      character caract*4,ilgg*4
      character banfil*32
      character*512 fnam
      character*512 grdfile
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
        read (*,20) ilgg
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
c     ***SET 3D GEOMETRY
c
      write (*,80)
   80 format(/' Enter grid file name : ',$)   
c
      read (*,90) fnam      
   90 format(a)
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
         read (luin,*) mcnx
         do l=1,mcnx+1
            read (luin,*) cwx(l)
         enddo 
         read (luin,*) mcny
         do l=1,mcny+1
            read (luin,*) cwy(l)
         enddo               
         read (luin,*) mcnz
         do l=1,mcnz+1
            read (luin,*) cwz(l)
         enddo    
c                       
      endif    

  100 format(/' Input Density File : ',$)
      read (*,20) fnam

      m=lenv(fnam)
      denfile=fnam(1:m)
c
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
c      qhlo = dlog(qht+epsilon)
c      rechy = 2.6d-13
c     reclo = dlog(rechy)
c
      alphahb=2.585d-13
      alphaheb=1.533d-12
c
c
c     Define freqency range of spectral sampling
c     Default range is from 5-100 eV
c
ccc      nrange=8
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
         read (*,90) ilgg
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
            read (*,*) rstsun
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
            read (*,90) ilgg
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
            read (*,*) starlum
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
            read (*,*) starlum
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
      read (*,*) npck, npckinc, convrt
      if (npck.lt.1.d2) npck=10.d0**npck
      if (convrt.lt.1) convrt=convrt*100.d0
c
      mcdengy0= dflum0/(npck*nrange)
c      
c  
c     ***SETTING INNER RADIUS
c
  220 write (*,225)
  225 format(//' Give the inner radius (in cm) '/
     &         ' (log (<100) or (>100)) : ',$)
      read (*,*) remp
c
      if (remp.lt.1.0d2) remp=10.d0**remp 
      if (remp.lt.rstar) remp=rstar
c
c
  230 fren=0.01d0
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
      read (*,90) jend
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
         read (*,*) diend
c         diend=remp+diend
         if (diend.le.remp) goto 260
      endif
c
      if (jend.eq.'A') then
         diend=-1.0d0
  280    write (*,290) fren
  290    format(/' Change the ionised hydrogen fraction',2pf5.2,'%?'/,
     &        ' (y/N) : ',$)
         read (*,90) jren   
c         
         if ((jren(1:1).eq.'Y').or.(jren(1:1).eq.'y')) then
  300       write (*,310) fren
  310       format(/' New ionised hydrogen fraction (0-1)',
     &           ' (currnt',2pf5.2,'%) : ',$)
            read (*,*) fren
            if ((fren.lt.0).or.(fren.gt.1)) goto 300
         endif
      endif
c
      difma=0.05d0
      dtlma=0.10d0
      dhlma=0.050d0
c
c
c      cwx(1)=mcxst-0.125d18
c      cwy(1)=mcyst
c      cwz(1)=mczst
c
c      do i=2,15
c         cwx(i) = cwx(i-1) + 0.125d18
c      enddo
c      do i=16,mcnx+1
c         cwx(i) = cwx(i-1) + 0.25d17
c      enddo  
c      do i=2,mcny+1
c         cwy(i) = cwy(i-1) + 0.25d18
c      enddo      
c      do i=2,mcnz+1
c         cwz(i) = cwz(i-1) + 0.25d18
c      enddo
c
      do i=1,mcnx+1
         write (*,*) i, cwx(i), cwy(i), cwz(i)
      enddo
cc      stop
c
c
cc      mcdy=cwy(2)-cwy(1)
cc      mcdz=cwz(2)-cwz(1)
cc      blum = .50d0*mcdz*mcdy
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
c
c     variable
c
      integer*4 ix,iy,iz
      integer*4 iel,ion,j,i,inl
      integer*4 m
c      integer*4 convg_arr(mcnx,mcny,mcnz)
      real*8 convg_arr(maxnx,maxny,maxnz)
      real*8 cov,cov0,incream,totion
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
      real*8 poppre(mxion, mxelem)
      character nmod*4
c
      real*8 f5007,f4363,f3727,f6300,f6584,f6720
c
      real*8 sigmt,tauso,sigmtso,inttphot,inttphot2
      real*8 dustsigmat, dusttau
      real*8 difde,difte,dhhe,difg
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
      dif(a,b)=dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
c
c
      dton=1.d33
      nmod='EQUI'
c
      write (*,10)  
   10 format(///
     &     ' ****************************************************'/
     &     '   MCITER SUBROUTINE IS RUN'/
     &     ' ****************************************************')
c
cccccc
      te_arr=2000.d0  ! initial guess
      dh_arr=100.d0  ! initial guess
      de_arr=100.d0  ! initial guess
      jnu_arr=0.0d0  
ccc      jnuw_arr=0.0d0  
      recpdf_arr = 0.0d0 
      rec0_arr = 0.0d0
      convg_arr = 100.0d0     
cccccc
c 
      m=0
      incream=1.d0
      cov=0.d0
c
      pop=0.d0 
      do i=1,atypes
        pop(1,i)=1.d-5
        pop(maxion(i),i)=1.d0-1.d-5                  
      enddo  
c      
c
      do ix = 1, mcnx
      do iy = 1, mcny
      do iz = 1, mcnz
         dhni=dh_arr(ix,iy,iz)          
         tei=te_arr(ix,iy,iz)
         dei=de_arr(ix,iy,iz)
c
         recpdf=0.d0
         call calpdf(tei,dei,dhni,recpdf,rec0)  
c
         rec0_arr(ix,iy,iz)=rec0
c         
         do inl=1,infph-1
            recpdf_arr(ix,iy,iz,inl)=recpdf(inl) 
            tauso=0.d0
            sigmt=0.d0
            dustsigmat=0.d0
            call crosssections (inl, tauso, sigmt, dustsigmat)
            sigmt_arr(ix,iy,iz,inl)=sigmt*dhni
         enddo  
c
c      write (*,*) ix, iy, iz,         
c     &            (cwx(ix)+cwx(ix+1))*0.5d0,
c     &            (cwy(iy)+cwy(iy+1))*0.5d0,
c     &            (cwz(iz)+cwz(iz+1))*0.5d0
c
      enddo
      enddo
      enddo      
c
      goto 40      
c      
c
c     Loop starts here
   20 continue
c
      cov0=cov
      cov=0.d0
      totion=0.d0
c         
      write (*,*) 'start'                
c
      nsum1=0.0d0
      nsum2=0.0d0
c      

      if (m.ge.15) then 
         open (lulin,file='./output/structure.out',status='unknown')
         write (lulin,*) m, '  x, y, z, Te, dne, Q1, Q2, 
     &                      REC0, DLOS, ELOSS, EGAIN'     
         close(lulin)

         call mcspec (luop1,1,'ZERO')

      endif
c
c
c
      do ix = 1, mcnx
      do iy = 1, mcny
      do 30 iz = 1, mcnz        
c
c
         dhni=dh_arr(ix,iy,iz)          
         tei=te_arr(ix,iy,iz)
         dei=de_arr(ix,iy,iz) 
c         
ccc         if ((convg_arr(ix,iy,iz).le.1).or.(dhni.le.epsilon)) goto 30       
ccc         if (dhni.le.epsilon) goto 30       
c
         call equion (tei, dei, dhni)                 
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
         if (inttphot.le.epsilon) goto 30                
c
         call intvec (tphot, q1, q2, q3, q4)         
c         
         q1field(ix,iy,iz)=q1
         q2field(ix,iy,iz)=q2
         q3field(ix,iy,iz)=q3
         q4field(ix,iy,iz)=q4
c
         call copypop (pop, poppre)
c
         ipho=1
         iphom=0 ! ionising field changes
c         write (*,*) ix,iy,iz
c         call mcteequi (tei, tef, def, dhni, nmod)     
         call teequi2 (tei, tef, def, dhni, 1.d33, nmod)
c
         dhnf=dhni
c
c        --------------------
c        convergence criteria
c
         call difhhe (pop, poppre, dhhe)
         difde=dif(def,dei)/dhlma
         difte=dif(tef,tei)/dtlma
         dhhe=dhhe/difma
         difg=dmax1(dhhe,difde,difte)   
c
c         write (*,*) dhhe, difde, difte, difg, ix, iy, iz
c
         if ((difg.le.1).and.(pop(1,1).le.fren)) cov=cov+1.d0
         if (pop(1,1).le.fren) totion=totion+1.d0
         convg_arr(ix,iy,iz)=difg
c
cc         tef=(tef+tei)/2.d0
cc         call equion (tef, def, dhnf)        
         te_arr(ix,iy,iz)=tef
         de_arr(ix,iy,iz)=def
c         
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
         call cool (tef, def, dhnf)                  
c
c
         if (m.ge.15) then
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
         open (lulin,file='./output/structure.out',
     &        status='OLD',access='APPEND') 
c
         write (lulin,*) (cwx(ix)+cwx(ix+1))*0.5d0,
     &                   (cwy(iy)+cwy(iy+1))*0.5d0,
     &                   (cwz(iz)+cwz(iz+1))*0.5d0,
     &                   tef,def,q1,q4,rec0
     &        ,hbeta,f5007,f4363,f3727,f6300,f6584,f6720     
     &        ,dlos,eloss,egain
     &        ,(pop(j,zmap(1)),j=1,maxion(zmap(1)))
     &        ,(pop(j,zmap(2)),j=1,maxion(zmap(2)))
     &        ,(pop(j,zmap(8)),j=1,9)
     &        ,(pop(j,zmap(16)),j=1,9)
     &        ,hloss,rloss,fslos,fmloss,xr3loss,xrlloss
     &        ,fflos,colos,f3loss,feloss,gcool
     &        ,pgain,cosgain,paheat,gheat,chgain,rngain  
c
         close (lulin)     
c
         dv=((cwx(ix+1)-cwx(ix))
     &      *(cwy(iy+1)-cwy(iy))
     &      *(cwz(iz+1)-cwz(iz)))
c
         call mcspec (luop1,dv,'SUMM')
c            
         endif
c
c         
   30 continue 
      enddo
      enddo 
c
c
      if (m.ge.15) then
c
   35       format(//' ======================================'//
     &               '     Lambda(A),  E (eV) , Flux (HB=1.0),',
     &               ' Species  , Kind, Accuracy (1-5)'/
     &               ' ====================================',
     &               '=============================')
c
           open (lulin,file='./output/all_spec.out',
     &          status='unknown') 
           write (lulin,35)
c
           call mcspec (lulin,1,'WRIT')
c
           close (lulin) 
c
      endif
c
c      
   40 m=m+1
c
      if (m.le.20) then
c
          npck = 1.d5
c         if (m.ge.3) npck = 1.d5        
c         if (m.ge.5) npck = 1.d7
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
c    
      logical   lgstr, lgdif, lgescape, lgabs
      integer*4 ix,iy,iz,icont
      integer*4 absEvent, seed0, seed
      integer*4 iphot,nphot,nphtlop,iloop
      integer*4 cwi(3),ci(3)
      integer*4 nup,nup0,m
      integer*4 xp, yp, zp
      integer*4 numstar,irg
c      
      real*8    blum,blum0,wid,invnph,invblum
      real*8    incream,npck0
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
      real*8    xorg,yorg,zorg
c
      real*8    sigmt 
c
c
cccccccccccccc    
      integer*4 lulin      
cccccccccccccc    
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
c
      write (*,*) diend
c      
cccccccccccccc    
      strpos(1) = 0.0d0
      strpos(2) = 0.0d0
      strpos(3) = 0.0d0
      numstar = 1
cccccccccccccc    
c
c     
ccc      irg=0
ccc   15 irg=irg+1
c
ccc      mcde=mcdengy(irg)/incream
c
      npck0=npck*incream
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
      invnph=1.d0/nphot     
      invblum=1.d0/blum0
c
cc      open (lulin,file='/Users/jyf/Desktop/code_check/engpack'
cc     &     ,status='unknown')
c
      do nup0 = 1, infph-1
c
         if ((photev(nup0).lt.mcnuin).or.(photev(nup0+1).gt.mcnuend))
     &   goto 50
  
         wid=photev(nup0+1)-photev(nup0)
         mcde=soupho(nup0)*wid*evplk*invnph*blum*invblum
c
         do iphot = 1, nphot             
c
            lgstr = .true.
            lgdif = .false.
            lgescape = .false.
c      
c           initialise photon position
            xorg=strpos(1)      
            yorg=strpos(2)
            zorg=strpos(3)
c
            call gtcwi(strpos,cwi)     
            xp=cwi(1)
            yp=cwi(2)
            zp=cwi(3)                   
c
            absEvent = 0
c
   20       if ((.not.lgescape).and.(absEvent.le.safelmt)) then 
c
c              update photon position
   30          rvec(1)=xorg
               rvec(2)=yorg
               rvec(3)=zorg
c
c              Initialise photon Nu & Direction
               if (.not.(lgstr.or.lgdif)) then
                  write(*,*)'Neither Stellar nor Diffuse Emission'                        
                  stop
               else
c                  if (lgstr) then
c                 Initialise photon Direction                
                  hvec(1)=0.0d0
                  hvec(2)=0.0d0
                  hvec(3)=0.0d0                  
                  call photdir (hvec)
ccc                  hvec(1)=1.0d0
ccc                  hvec(2)=0.0d0
ccc                  hvec(3)=0.0d0
c
c                 Initialise photon cell index                
                  ci(1)=xp
                  ci(2)=yp
                  ci(3)=zp                   
c                  endif
c
c                 Initialise photon Frequency                
                  if (lgstr.and.(.not.lgdif)) then  
c                    stellar photon   
                     nup=nup0
cc                     call frqgen (nup,nust,nued,ci,lgstr,lgdif,numstar)
                  endif
                  if ((.not.lgstr).and.lgdif) then  
c                    diffuse photon   
c                     seed=86432
                     seed = time() 
                     random=ran2(seed)  
                     if (random.lt.rec0_arr(ci(1),ci(2),ci(3))) then
                        lgescape = .true.
                        goto 20                        
                     endif                     
                     call frqgen (nup,mcnuin,mcnuend
     &                           ,ci,lgstr,lgdif,numstar)
                  endif             
                  lgabs=.false.
               endif
c              distance to center star
               drvec=0.0d0
               drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                    +(rvec(2)-strpos(2))**2.0d0
     &                    +(rvec(3)-strpos(3))**2.0d0)
c
               if ((photev(nup).lt.mcnuin).or.(drvec.gt.diend)) then 
                  lgescape = .true.
                  goto 20
               endif
c
c              set random number for absorption event
               abstau=0.0d0
               seed0=86432
               seed=seed0
               seed = time()                
               random=0.0d0
               random=ran2(seed)            
               passprob=-log(1.0d0-random)
cc               passprob=10
c 
c              Start path loop
c              Looking for dSx & determining Cell x-index             
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
c               
c              Whether photon reach Rin (remp) or useless cell
               if ((drvectmp.le.remp).and.(drvec.lt.remp)
     &            .and.(drvectmp.gt.rstar)) then
c
                  if ((ds.eq.dsx).and.(hvec(1).gt.1d-10)) then
                      xp=xp+1 
                  else 
                  if ((ds.eq.dsx).and.(hvec(1).lt.-1d-10)) then
                      xp=xp-1
                  else    
                  if ((ds.eq.dsy).and.(hvec(2).gt.1d-10)) then
                      yp=yp+1 
                  else
                  if ((ds.eq.dsy).and.(hvec(2).lt.-1d-10)) then
                      yp=yp-1
                  else    
                  if ((ds.eq.dsz).and.(hvec(3).gt.1d-10)) then
                      zp=zp+1 
                  else    
                  if ((ds.eq.dsz).and.(hvec(3).lt.-1d-10)) then
                      zp=zp-1
                  endif 
                  endif 
                  endif 
                  endif 
                  endif
                  endif                
c 
                  if ((xp.le.0).or.(yp.le.0).or.(zp.le.0)
     &            .or.(xp.ge.mcnx).or.(yp.ge.mcny).or.(zp.ge.mcnz)) then
                     lgescape = .true.
                     goto 20
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
                  goto 40        
               endif
c                
c              Whether photon reach the surface of star (rstar)
               if (drvectmp.le.rstar) then
                  lgescape = .true.
                  goto 20               
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
c              Whether the photon is absorbed?
               if (((tauloc+abstau).gt.passprob).and.(tauloc.gt.0.0d0))
     &         then                    
                  ds=(passprob-abstau)/sigmt
                  rvec(1)=rvec(1)+ds*hvec(1)
                  rvec(2)=rvec(2)+ds*hvec(2)
                  rvec(3)=rvec(3)+ds*hvec(3)
c                 
ccc                  if ((ci(1).eq.2).and.(ci(2).eq.13).and.(ci(3).eq.13)
ccc     &               ) then
ccc                   if (lgstr) jnuw_arr(2,nup)=jnuw_arr(2,nup)+ds*mcde/dv 
ccc                   if (lgdif) jnuw_arr(1,nup)=jnuw_arr(1,nup)+ds*mcde/dv 
ccc                  endif
c      
                  jnu_arr(ci(1),ci(2),ci(3),nup)
     &             =jnu_arr(ci(1),ci(2),ci(3),nup)+ds*mcde/dv
c
c
                  lgstr=.false.
                  lgdif=.true.
c
                  drvec=0.0d0
                  drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                       +(rvec(2)-strpos(2))**2.0d0
     &                       +(rvec(3)-strpos(3))**2.0d0)                  
                  if ((drvec.ge.diend).or.
     &                (rvec(1).le.cwx(1)).or.
     &                (rvec(1).ge.cwx(mcnx+1)).or.
     &                (rvec(2).le.cwy(1)).or.
     &                (rvec(2).ge.cwy(mcny+1)).or.
     &                (rvec(3).le.cwz(1)).or.
     &                (rvec(3).ge.cwz(mcnz+1))) then
                     lgescape=.true.
                  else
                     lgabs=.true.
                  endif 
c               
               else
                  rvec(1)=rvec(1)+ds*hvec(1)
                  rvec(2)=rvec(2)+ds*hvec(2)
                  rvec(3)=rvec(3)+ds*hvec(3)
c
                  abstau=abstau+tauloc        
c
ccc                  if ((ci(1).eq.2).and.(ci(2).eq.13).and.(ci(3).eq.13)
ccc     &               ) then
ccc                   if (lgstr) jnuw_arr(2,nup)=jnuw_arr(2,nup)+ds*mcde/dv 
ccc                   if (lgdif) jnuw_arr(1,nup)=jnuw_arr(1,nup)+ds*mcde/dv 
ccc                  endif  
c                  
                  jnu_arr(ci(1),ci(2),ci(3),nup)
     &             =jnu_arr(ci(1),ci(2),ci(3),nup)+ds*mcde/dv        
c
                  drvec=0.0d0
                  drvec=dsqrt((rvec(1)-strpos(1))**2.0d0
     &                       +(rvec(2)-strpos(2))**2.0d0
     &                       +(rvec(3)-strpos(3))**2.0d0)    
c
                  if ((drvec.ge.diend).or.
     &                (rvec(1).le.cwx(1)).or.
     &                (rvec(1).ge.cwx(mcnx+1)).or.
     &                (rvec(2).le.cwy(1)).or.
     &                (rvec(2).ge.cwy(mcny+1)).or.
     &                (rvec(3).le.cwz(1)).or.
     &                (rvec(3).ge.cwz(mcnz+1))) then
                     lgescape=.true.
                  else
c
                     if ((ds.eq.dsx).and.(hvec(1).gt.1d-10)) then
                        xp=xp+1 
                     else 
                     if ((ds.eq.dsx).and.(hvec(1).lt.-1d-10)) then
                        xp=xp-1
                     else    
                     if ((ds.eq.dsy).and.(hvec(2).gt.1d-10)) then
                        yp=yp+1 
                     else
                     if ((ds.eq.dsy).and.(hvec(2).lt.-1d-10)) then
                        yp=yp-1
                     else    
                     if ((ds.eq.dsz).and.(hvec(3).gt.1d-10)) then
                        zp=zp+1 
                     else    
                     if ((ds.eq.dsz).and.(hvec(3).lt.-1d-10)) then
                        zp=zp-1
                     endif 
                     endif 
                     endif 
                     endif 
                     endif
                     endif                   
c 
                     if ((xp.le.0).or.(yp.le.0).or.(zp.le.0)
     &               .or.(xp.ge.mcnx).or.(yp.ge.mcny)
     &               .or.(zp.ge.mcnz)) then
                        lgescape = .true.
                        goto 20
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
                  xorg=rvec(1)
                  yorg=rvec(2)
                  zorg=rvec(3)
                  goto 30 
               else 
                  if(lgescape) then 
                    goto 20
                  else
                    write(*,*)'Wrong in photon travelling'                    
                    stop  
                  endif
               endif                                 
c
c           end label 20
            endif         
         enddo
   50    continue      
      enddo   
c
ccc      if (irg.lt.nrange) goto 15
c
cc      close(lulin)
c
cc      stop
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
         write(*,*)'Wrong in photon position-X'                    
         stop
      endif  
      if (rvec(1).lt.cwx(i+1)) then 
         cwi(1)=i                          
      else
         if (i.lt.mcnx) then
            i=i+1
            goto 10
         else 
            write(*,*)'Cannot find cell wall in X'                    
            stop             
         endif  
      endif         
c
c
      i=1   
   20 if (rvec(2).lt.cwy(i)) then 
         write(*,*)'Wrong in photon position-Y'                    
         stop
      endif  
      if (rvec(2).lt.cwy(i+1)) then 
         cwi(2)=i                          
      else
         if (i.lt.mcny) then
            i=i+1
            goto 20
         else 
            write(*,*)'Cannot find cell wall in Y'                    
            stop             
         endif  
      endif   
c
c
      i=1   
   30 if (rvec(3).lt.cwz(i)) then 
         write(*,*)'Wrong in photon position-Z'                    
         stop
      endif  
      if (rvec(3).lt.cwz(i+1)) then 
         cwi(3)=i                          
      else
         if (i.lt.mcnz) then
            i=i+1
            goto 30
         else 
            write(*,*)'Cannot find cell wall in Z'                    
            stop             
         endif  
      endif  
c  
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
      integer*4   i, seed0
      real*8      random, ran2
      real*8      hvec(3) 
      real*8      w,u,v,t,ang
c
      seed0 = time() 
c      
      random = ran2(seed0)
      w = 2.*random - 1.
      t = sqrt(1.-w*w)
      random = ran2(seed0)
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
      integer*4 seed0, seed
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
   20 seed0 = 86432
      seed0 = time()                   
      random = ran2(seed0)
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
