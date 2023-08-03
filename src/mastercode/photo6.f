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
      subroutine photo6 ()
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      real*8 epotmi,fin,uinit,qinit,rcin
      real*8 rstsun,wid,starlum
      real*8 scale,dnt,lam
c      real*8 fluxes(mxmonlines)
c
      integer*4 j,luin,i,kmax,k,nl
      integer*4 idx1 ,idx2, idx3, idx4, iem
c      character=jtb*4
      character carac*36,model*32
      character caract*4,ilgg*4
      character banfil*32
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
c
      write (*,10)
c
c     ***INITIAL IONISATION CONDITIONS
c
   10 format(///
     &     ' *********************************************************'/
     &     '  Photoionisation model P6 selected:'//
     &     '  Diffuse Field : Robust integration: Radiation Pressure',/,
     &     ' *********************************************************')
c
      epotmi=iphe
c
c     set up ionisation state
c
      model='protoionisation'
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
c     ***CHOOSING PHOTON SOURCE AND GEOMETRICAL PARAMETERS
c
      model='Photo 6'
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
c Default outward only intgration assuming spherical symmetric nebula
c
      jtrans='OUTW'
c
   80 write (*,100)
   90 format(a)
  100 format(/' Spherical or plane II  geometry (s/p) : ',$)
      read (*,90) jgeo
      jgeo=jgeo(1:1)
c
      if ((jgeo.eq.'S').or.(jgeo.eq.'s')) jgeo='S'
      if ((jgeo.eq.'P').or.(jgeo.eq.'p')) jgeo='P'
c
      if ((jgeo.ne.'S').and.(jgeo.ne.'P')) goto 80
c
      rstar=1.d0
      astar=rstar*rstar
c
c     ***IF SYMMETRY IS plane parallel :
c
      if (jgeo.eq.'P') then
c
        jtrans='LODW'
c
        write (*,110)
  110  format(/' Radiative Transfer Mode '/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    T   :  Two Sided, outward only (symmetrical, default)'/
     & '    O   :  One Sided, 0.5 upstream attenuated (precursors)'/
     & ' :: ',$)
        read (*,90) ilgg
        ilgg=ilgg(1:1)
        if (ilgg.eq.'t') ilgg='T'
        if (ilgg.eq.'o') ilgg='O'
        if ((ilgg.ne.'O').and.(ilgg.ne.'T')) ilgg='T'
        if (ilgg.eq.'T') then
          jtrans='OUTW'
        else
          jtrans='LODW'
        endif
c     ***IF SYMMETRY IS plane parallel :
      endif
c
c     ***IF SYMMETRY IS SPHERICAL :
c
      if (jgeo.eq.'S') then
c
        jtrans='OUTW'
c
c     1st total energy in Inu
c
c
        blum=0.d0
        ilum=0.d0
c
        do i=1,infph-1
          wid=photev(i+1)-photev(i)
          blum=blum+soupho(i)*wid*evplk
          if (photev(i).ge.iph) ilum=ilum+soupho(i)*wid*evplk
        enddo
c
c     Plane Parallel x pi
c
        blum=pi*blum
        ilum=pi*ilum
c
        rstar=1.d0
        astar=rstar*rstar
c
        if (blum.gt.0.d0) then
c
  120     write (*,130)
  130 format(/' Define the source size or luminosity '/
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
     &     120
c
          if (ilgg.eq.'R') then
c
  140      format(/,/,' Give photoionisation source radius'
     &            ,'(Rsun=6.96e10)',/
     &            ,' (in solar units (<1.e5) or in cm (>1.e5) : ',$)
            write (*,140)
            read (*,*) rstsun
            if (rstsun.le.0.0d0) goto 120
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
c
  150 format(//' ****************************************************'/
     & '  Source Total Luminosity              : ',1pg12.5/
     & '  Source Ionising (13.6eV+) Luminosity : ',1pg12.5/
     & '  Source Ionising (13.6eV+) Photons    : ',1pg12.5/
     & ' ****************************************************')
            write (*,150) blum,ilum,qht*astar
c
c     end def by radius
c
          endif
c
          if (ilgg.eq.'L') then
c
  160       format(//' Total or Ionising Luminosity (T/I)',$)
            write (*,160)
            read (*,90) ilgg
            ilgg=ilgg(1:1)
c
            if (ilgg.eq.'t') ilgg='T'
            if (ilgg.eq.'i') ilgg='I'
c
            if ((ilgg.ne.'I').and.(ilgg.ne.'T')) ilgg='I'
c
  170       format(//' Give ionising source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
  180       format(//' Give bolometric source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
            if (ilgg.eq.'T') then
              write (*,180)
            else
              write (*,170)
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
  190  format(//' ****************************************************'/
     &      '  Source Radius : ',1pg12.5,' cm'/
     &      ' ****************************************************')
            write (*,190) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            write (*,150) blum,ilum,qht*astar
c
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c     end def by luminosity
c
          endif
c
          if (ilgg.eq.'P') then
c
  200       format(//' Give source ionising photon rate '/
     &               ' (log (<100) or photons/s (>100)) : ',$)
            write (*,200)
            read (*,*) starlum
            if (starlum.lt.1.d2) starlum=10.d0**starlum
c
c     P/s = qht*fpi*rstsun*rstsun
c
            rstsun=dsqrt(starlum/(qht*fpi))
c
c
            write (*,190) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            write (*,150) blum,ilum,qht*astar
c
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c    end source radius with photons
c
          endif
c
c    end blum>0
c
        endif
c
      endif
c
c
  210 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Setting the Physical Structure  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      write (*,210)
c
c     ***THERMAL STRUCTURE
c
      jthm='S'
      if (expertmode.gt.0) then
c
  220   write (*,230)
  230 format(/' Choose a Thermal Structure '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    S  : Self-consistent,  (recommended)'/
     &     '    T  : isoThermal, (non-physical, not recommended)'/
     &     ' :: ',$)
        read (*,90) jthm
c
        if ((jthm(1:1).eq.'S').or.(jthm(1:1).eq.'s')) jthm='S'
        if ((jthm(1:1).eq.'T').or.(jthm(1:1).eq.'t')) jthm='T'
        if ((jthm(1:1).eq.'F').or.(jthm(1:1).eq.'f')) jthm='F'
        if ((jthm(1:1).eq.'I').or.(jthm(1:1).eq.'i')) jthm='I'
c
        if ((jthm.ne.'S').and.(jthm.ne.'T').and.(jthm.ne.'F')
     &   .and.(jthm.ne.'I')) goto 220
c
        if (jthm.eq.'T') then
c
  240    format(/,' Give the fixed electron temperature (<=10 as log):'
     &         ,$)
          write (*,240)
          read (*,*) fixtemp
c
          if (fixtemp.le.10.d0) fixtemp=10.d0**fixtemp
c
        endif
      endif
c
c
c     ***DENSITY BEHAVIOR
c
c default partial pressure of Hydrogen, ignored if not isobaric
c
      ponk=1.d6
c
  250 write (*,260)
  260 format(/'  Choose a Density Structure '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    C  : isoChoric,  (const volume)'/
     &     '    B  : isoBaric, (const pressure)'/
     &     '    F  : Functional form, F(r)'/
     &     ' :: ',$)
      read (*,90) jden
c
      if ((jden(1:1).eq.'C').or.(jden(1:1).eq.'c')) jden='C'
      if ((jden(1:1).eq.'B').or.(jden(1:1).eq.'b')) jden='B'
      if ((jden(1:1).eq.'F').or.(jden(1:1).eq.'f')) jden='F'
c
      if ((jden.ne.'C').and.(jden.ne.'B').and.(jden.ne.'F')) goto 250
c
      if (jden.eq.'F') then
        write (*,270)
  270    format(/'  Density Function '/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &        ' n(x) = n0*[fx*exp((x-a)/r0) + fp*((x/a)**b)]+c  '/
     &        ' Give parameters n0 fx fp a b c & r0 (cgs): ',$)
        read (*,*) ofac,xfac,pfac,afac,bfac,cfac,scalen
        dhn=ofac
      endif
c
c     Estimate strom radius for nominal density
c     if we have any photons that is
c
c     if isobaric, get pressure regime
c
      if (jden.eq.'B') then
c
  280    format(/' Give TOTAL pressure ',
     &           ' (n.T = P/k, <10 as log) : ',$)
        write (*,280)
        read (*,*) ponk
        if (ponk.le.10.d0) ponk=10.d0**ponk
c
        tinner=1.d4
c
  290    format(//' Give an estimate of the initial temperature :',//
     & '   * Use 1.0e4 if not known. Initial estimates of density,',/
     & '   * and radius will be uncertain until the model starts.',/
     & '   * Ionisation uncertainty will always make nH estimates',/
     & '   * here approximate, but total pressure will be accurate.',/
     & '   * Initial Q, P/k, T and nH is solved in model,',/
     & '   * so this estimate is not crucial.',//
     & ' [T ~1e4, <10 as log] : ',$)
        write (*,290)
        read (*,*) tinner
        if (tinner.le.10.d0) tinner=10.d0**tinner
c
        dnt=ponk/tinner
        dhn=dnt/2.2d0
c
  300 format(/' ****************************************************',/
     &        '  Total number density at',1pg10.3,': ',1pg10.3,/
     &        '  Est. Hydrogen number density at',1pg10.3,': ',1pg10.3/
     &        ' ****************************************************')
        write (*,300) tinner,dnt,tinner,dhn
c
c     jden = B isobaric
c
      endif
c
  310 if (jden.eq.'C') then
  320    format(//' Give hydrogen number density : ',$)
        write (*,320)
        read (*,*) dhn
      endif
c
      dht=densnum(dhn)
c
  330 format(//' ****************************************************'/
     &        '   Mean Ion Density (N):',1pg12.5/
     &        '   Mean H  Density (dh):',1pg12.5/
     &        ' ****************************************************')
c
      write (*,330) dht,dhn
c
  340 format(/' Set filling factor (0<f<=1) : ',$)
      write (*,340)
      read (*,*) fin
      if (((fin.le.0.0d0).or.(fin.gt.1.0d0)).or.(dhn.le.0.0d0)) goto
     &310
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ((jgeo.eq.'S')) then
c
c       Spherical Geometry
c
        if (qht.gt.0.d0) then
c
c
c         vstrlo = ((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(
c     &        dht)))-dlog(fin))-reclo
c         rmax = dexp((vstrlo-1.43241d0)/3.d0)
c
c         write(*,*) "original rmax", rmax
c         dht ~ de fully ionised
          volstromhb=(qht*astar)/(dhn*dht*alphahb*fin)
          rstromhb=(3.d0*volstromhb/fpi)**0.333333333333d0
          volstromheb=(qheii*astar)/(zion(2)*dhn*dht*alphaheb*fin)
          rstromheb=(3.d0*volstromheb/fpi)**0.333333333333d0
c
          rmax=dmax1(rstromhb,rstromheb)
c
c
c         write (*,*) "rstromHB", rstromHB
c         write (*,*) "rstromHeB",rstromHeB
c
c     write (*,*) qht,dht,dhn
c
          qhdnin=qht/dht
          qhdnav=(4.d0*qht/dht)*fdilu(rstar,(rmax+rstar)*0.5d0)
c
          unin=qhdnin/cls
          unav=qhdnav/cls
c
          qhdhin=qht/dhn
          qhdhav=qhdnav*(dht/dhn)
c
          uhin=qhdhin/cls
          uhav=qhdhav/cls
c
  350     write (*,360) rstromhb,rstromheb,qht*astar,qheii*astar,qht,
     &     qheii,qhdnin,qhdhin,qhdnav,qhdhav,unin,uhin,unav,uhav
  360 format(//,
     & ' ****************************************************',/,
     & '   Filled Sphere Parameters:',/,
     & ' ****************************************************',/,
     & '  Estimated HII   Stromgren radius:',1pg10.3,' cm.',/,
     & '  Estimated HeIII Stromgren radius:',1pg10.3,' cm.',/,
     & ' ****************************************************',/,
     & '  Photon Luminosity LQH    :',1pg10.3, ' Phot/s',/,
     & '  Photon Luminosity LQHeII :',1pg10.3, ' Phot/s',/,
     & '  Photon Flux       FQH    :',1pg10.3, ' Phot/cm2/s',/,
     & '  Photon Flux       FQHeII :',1pg10.3, ' Phot/cm2/s',/,
     & ' ****************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s'/,
     & '        QHDH inner   : ',1pg12.5,' cm/s'/,
     & '          <QHDN>     : ',1pg12.5,' cm/s'/,
     & '          <QHDH>     : ',1pg12.5,' cm/s'/,
     & '        U(N) inner   : ',1pg12.5,/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <U(N)>     : ',1pg12.5,/,
     & '          <U(H)>     : ',1pg12.5,/,
     & ' ****************************************************',//
     & ' Give initial radius in terms of distance or Q(N), U(N),',
     & ' Q(H), or U(H) (d/q/n/h/u):',$)
          read (*,20) ilgg
          ilgg=ilgg(1:1)
c
          if (ilgg.eq.'D') ilgg='d'
          if (ilgg.eq.'Q') ilgg='q'
          if (ilgg.eq.'N') ilgg='n'
          if (ilgg.eq.'H') ilgg='h'
          if (ilgg.eq.'U') ilgg='u'
          if ((ilgg.ne.'q').and.(ilgg.ne.'h').and.(ilgg.ne.'d')
     &     .and.(ilgg.ne.'u').and.(ilgg.ne.'n')) goto 350
c
          if (ilgg.eq.'u') then
            write (*,370)
  370       format(//,' Give U(H) at inner radius (<=0 as log) : ',$)
            read (*,*) uinit
c
            if (uinit.le.0) uinit=10.d0**uinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(uhav/uinit)
c
          endif
c
          if (ilgg.eq.'n') then
            write (*,380)
  380       format(//,' Give U(N) at inner radius (<=0 as log) : ',$)
            read (*,*) uinit
c
            if (uinit.le.0) uinit=10.d0**uinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(unav/uinit)
c
          endif
c
          if (ilgg.eq.'q') then
            write (*,390)
  390       format(//,' Give QHDN at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(qhdnav/qinit)
c
          endif
c
          if (ilgg.eq.'h') then
            write (*,400)
  400       format(//,' Give QHDH at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(qhdhav/qinit)
c
          endif
c
          if (ilgg.eq.'d') then
c
c     give R_0 directly
c
            write (*,410)
  410       format(//' Give the initial radius ',
     &         ' (in cm (>1) or in fraction of Stromgren radius) : ',$)
            read (*,*) remp
c
            if (remp.lt.0.0d0) goto 350
            if (remp.le.1) remp=remp*rmax
            if (remp.lt.rstar) remp=rstar
c
          endif
c
c         rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax)**3)))-
c     &        1.43241d0)/3.0d0)
c
c
c         write (*,*) "rmax", rmax
c
          rstromhb=((rstromhb**3.d0)+(remp**3.d0))**0.333333333333d0
          rmax=rstromhb
c
c         write (*,*) "rstromHB", rstromHB
c
          rstromheb=((rstromheb**3.d0)+(remp**3.d0))**0.333333333333d0
c
c         write (*,*) "rstromHeB",rstromHeB
c
          blum=0.d0
          ilum=0.d0
c
          do i=1,infph-1
            wid=photev(i+1)-photev(i)
            blum=blum+soupho(i)*wid*evplk
            if (photev(i).ge.iph) ilum=ilum+soupho(i)*wid*evplk
          enddo
          blum=pi*blum
          ilum=pi*ilum
          wdpl=fdilu(rstar,remp)
          qhdnin=4.d0*qht/dht*wdpl
          unin=qhdnin/cls
          qhdhin=qhdnin*(dht/dhn)
          uhin=qhdhin/cls
          blum=blum*wdpl
          ilum=ilum*wdpl
          wdpl=fdilu(rstar,(rmax+remp)*0.5d0)
          qhdnav=4.d0*qht/dht*wdpl
          unav=qhdnav/cls
          qhdhav=qhdnav*(dht/dhn)
          uhav=qhdhav/cls
  420 format(//,
     & ' ****************************************************',/,
     & '   Partially Filled Sphere Parameters:',/,
     & ' ****************************************************',/,
     & '  Empty inner radius :',1pg12.5,' cm'/,
     & '  Outer HII   radius :',1pg12.5,' cm.',/,
     & '  Outer HeIII radius :',1pg12.5,' cm.',/,
     & ' ****************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s'/,
     & '        QHDH inner   : ',1pg12.5,' cm/s'/,
     & '          <QHDN>     : ',1pg12.5,' cm/s'/,
     & '          <QHDH>     : ',1pg12.5,' cm/s'/,
     & '        U(N) inner   : ',1pg12.5,/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <U(N)>     : ',1pg12.5,/,
     & '          <U(H)>     : ',1pg12.5,/,
     & '   Total intensity   : ',1pg12.5,' erg/s/cm2' /,
     & '   Ionizing intensity: ',1pg12.5,' erg/s/cm2'/,
     & ' ****************************************************',/)
          write (*,420) remp,rstromhb,rstromheb,qhdnin,qhdhin,qhdnav,
     &     qhdhav,unin,uhin,unav,uhav,blum,ilum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          radin=0.d0
          radout=1.d38
c
          write (*,430)
  430    format(//' Volume integration over the whole sphere?',
     &        ' (y/n) : ',$)
          read (*,90) ilgg
c
          if ((ilgg(1:1).eq.'Y').or.(ilgg(1:1).eq.'y')) ilgg='Y'
          if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg='N'
c
          if (ilgg.ne.'N') ilgg='Y'
c
          if (ilgg.eq.'N') then
  440       write (*,450)
  450       format(/' Integration through the line of sight ',
     &           'for a centered ring aperture.'/
     &           'Give inner and outer radii (cm) :')
            read (*,*) radin,radout
            if ((radin.lt.0.0d0).or.(radin.ge.(0.999d0*radout))) goto
     &       440
          endif
c
c     end with photons (qht>0)
c
        endif
c
c   end Spherical
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***IF GEOMETRY PLANE-PARALLEL :
c
      if (jgeo.eq.'P') then
        rstar=1.0d0
        astar=1.0d0
        remp=0.0d0
c
        if (qht.gt.0.d0) then
c
c     Determine number of ionising photons
c
c
c     1st total energy in Inu
c
c
          blum=0.d0
          ilum=0.d0
c
          do i=1,infph-1
            wid=photev(i+1)-photev(i)
            blum=blum+soupho(i)*wid*evplk
            if (photev(i).ge.iph) ilum=ilum+soupho(i)*wid*evplk
          enddo
c
c     Plane Parallel x pi
c
          blum=pi*blum
          ilum=pi*ilum
c
  460     write (*,470)
  470    format(/,' Give Ionizing Flux at inner edge by: ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    B  : Total Flux (erg/s/cm^2)',/,
     &     '    I  : Ionising Flux (erg/s/cm^2)',/,
     &     '    F  : Ionising Photons FQHI (phot/s/cm^2)',/,
     &     '    Q  : Ionisation parameter QHDN ',/,
     &     '    H  : Ionisation parameter QHDH ',/,
     &     '    N  : Ionisation parameter U(N) ',/,
     &     '    U  : Ionisation parameter U(H) ',/,
     &     '    X  : No change (Use current flux)',/,
     &     ' :: ',$)
          read (*,90) ilgg
c
          if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg='B'
          if ((ilgg(1:1).eq.'I').or.(ilgg(1:1).eq.'i')) ilgg='I'
          if ((ilgg(1:1).eq.'F').or.(ilgg(1:1).eq.'f')) ilgg='F'
          if ((ilgg(1:1).eq.'Q').or.(ilgg(1:1).eq.'q')) ilgg='Q'
          if ((ilgg(1:1).eq.'H').or.(ilgg(1:1).eq.'h')) ilgg='H'
          if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg='N'
          if ((ilgg(1:1).eq.'U').or.(ilgg(1:1).eq.'u')) ilgg='U'
          if ((ilgg(1:1).eq.'X').or.(ilgg(1:1).eq.'x')) ilgg='X'
          if (ilgg.eq.'B') then
            write (*,480)
  480       format(/,'Give Bolometric flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/blum
          elseif (ilgg.eq.'I') then
            write (*,490)
  490       format(/,'Give Ionising flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/ilum
          elseif (ilgg.eq.'F') then
            write (*,500)
  500       format(/,
     &      'Give Ionising Photon flux at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale/qht
          elseif (ilgg.eq.'U') then
            write (*,510)
  510     format(/,'Give U(H) at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0.d0) scale=10**scale
            scale=scale*dhn*cls/qht
          elseif (ilgg.eq.'N') then
            write (*,520)
  520     format(/,'Give U(N) at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0.d0) scale=10**scale
            scale=scale*dht*cls/qht
          elseif (ilgg.eq.'Q') then
            write (*,530)
  530       format(/,'Give QHDN at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale*dht/qht
          elseif (ilgg.eq.'H') then
            write (*,540)
  540       format(/,'Give QHDH at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale*dhn/qht
          elseif (ilgg.eq.'X') then
            scale=1.d0
          else
            goto 460
          endif
c
          do j=1,infph
            soupho(j)=soupho(j)*scale
          enddo
c
          blum=blum*scale
          ilum=ilum*scale
          qht=qht*scale
          qhdnav=qht/dht
          unav=qhdnav/cls
          qhdhav=qht/dhn
          uhav=qhdhav/cls
  550  format(//,
     & ' ****************************************************',/,
     & '  Total Flux at inner      : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Flux at inner   : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Phot Flux FQ, Ryd+  : ',1pg12.5,' (phots/s/cm^2)',/,
     & '  Ionisation parameter QHDN inner : ',1pg12.5/
     & '  Ionisation parameter QHDH inner : ',1pg12.5/
     & '  Ionisation parameter U(N) inner : ',1pg12.5/
     & '  Ionisation parameter U(H) inner : ',1pg12.5/
     & ' ****************************************************',/)
          write (*,550) blum,ilum,qht,qhdnav,qhdhav,unav,uhav
c
          volstromhb=qht/(dhn*dht*alphahb*fin)
          rstromhb=volstromhb
          rmax=rstromhb
          volstromheb=qheii/(zion(2)*dhn*dht*alphaheb*fin)
          rstromheb=volstromheb
c
c        Dilution = 0.5, qht includes geo dil of 0.5
c
  560     write (*,570)
  570    format(/' Give the geometrical dilution factor (<=0.5) : ',$)
          read (*,*) wdpl
c
          if ((wdpl.le.0.0d0).or.(wdpl.gt.0.5d0)) goto 560
c
          rmax=rmax*wdpl
c
          qhdnin=((2.0d0*qht)*wdpl)/dht
          qhdnav=qhdnin
c
          unin=((2.0d0*qht)*wdpl)/(cls*dht)
          unav=unin
c
          qhdhin=((2.0d0*qht)*wdpl)/dhn
          qhdhav=qhdhin
c
          uhin=((2.0d0*qht)*wdpl)/(cls*dhn)
          uhav=uhin
c
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  580 admach=0.0d0
      turbheatmode=0
c 440  format(//,' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     &  '  Micro-Turbulent Dissipation Heating',/,
c     &  ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     &  '     Choose Timescale Type:',/,
c     &  '     A  : Recombination Timescale.',/,
c     &  '     B  : Collisional Timescale.',/,
c     &  '     C  : Cooling Timescale.',/,
c     &  '     D  : Fixed Timescale.',/,
c     &  '     N  : No Dissipaton.',/,
c     &  '   ::', $)
c      write(*,440)
c
c      Read(*, 350), jtb
c
c      if ((jtb(1:1).eq.'A').or.(jtb(1:1).eq.'a')) jtb = 'A'
c      if ((jtb(1:1).eq.'B').or.(jtb(1:1).eq.'b')) jtb = 'B'
c      if ((jtb(1:1).eq.'C').or.(jtb(1:1).eq.'c')) jtb = 'C'
c      if ((jtb(1:1).eq.'D').or.(jtb(1:1).eq.'d')) jtb = 'D'
c      if ((jtb(1:1).eq.'N').or.(jtb(1:1).eq.'n')) jtb = 'N'
c
c      if ( jtb .eq. 'A') turbheatmode = 1
c      if ( jtb .eq. 'B') turbheatmode = 2
c      if ( jtb .eq. 'C') turbheatmode = 3
c      if ( jtb .eq. 'D') turbheatmode = 4
c
c      if ( turbheatmode .gt. 0) then
c
c 450  format(/,' Adiabatic Mach Number: ',$)
c      write(*,450)
c      read(*, *) admach
c
c      if ( admach .lt. 0.0d0) admach = 0.0d0
c
c      if ( turbheatmode .eq. 4) then
c 451  format(/,' Fixed Timescale (s): ',$)
c      write(*,451)
c      read(*, *) alphaTurbFixed
c      if (alphaTurbFixed .lt. 0.0) alphaTurbFixed = 0.d0
c      endif
c
c      endif
c
c      if ( turbheatmode .ge. 1 ) then
c 455  format(//,' ::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     &  '  Micro-Turbulent Dissipation Enabled',/,
c     &  ' ::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c      write(*,455)
c       if (turbheatmode .eq. 1) then
c       write(*,*) 'On Recombination Timescale.'
c       endif
c       if (turbheatmode .eq. 2) then
c       write(*,*) 'On Collisional Timescale.'
c       endif
c       if (turbheatmode .eq. 3) then
c       write(*,*) 'On Cooling Timescale.'
c       endif
c       if (turbheatmode .eq. 4) then
c       write(*,*) 'On Fixed Timescale:', alphaTurbFixed, 's'
c       endif
c      endif
      frlum=0.0d0
      tm00=0.0d0
c
c     ***EQUILIBRIUM OR FINITE AGE ASSUMPTION
c
  590 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Ionisation Balance Calculations'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
      write (*,590)
  600 write (*,610)
  610 format(//'  Choose type'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '     E  : Equilibrium balance.'/
     &     '     F  : Finite source lifetime.'/
     &     '     P  : Post-Equilibrium decay.'/
     &     '  :: ',$)
      read (*,90) jeq
c
      if ((jeq(1:1).eq.'E').or.(jeq(1:1).eq.'e')) jeq='E'
      if ((jeq(1:1).eq.'F').or.(jeq(1:1).eq.'f')) jeq='F'
      if ((jeq(1:1).eq.'P').or.(jeq(1:1).eq.'p')) jeq='P'
      if ((jeq(1:1).eq.'c').or.(jeq(1:1).eq.'c')) jeq='C'
c
      if ((jeq.ne.'E').and.(jeq.ne.'F').and.(jeq.ne.'P')) goto 600
c
c     forced CIE for given tprofile.
c     useful for hot x-ray bubbles.
c
      if (jeq.eq.'C') then
c
        write (*,620)
  620    format(/'  Fixed Thermal Profile:'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',//
     &     ' f(x) = t0*[fx*exp((x-ta)/r0) + fp*((x/ta)**tb)]+tc  ',/
     &     ' (with r in r0 units , t0<10K as log)',/
     &     ' Give parameters t0 ta tb tc r0 : ',$)
        read (*,*) tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
      endif
c
c     non equilibrium
c
      if ((jeq.ne.'E').and.(jeq.ne.'C')) then
        rcin=0.9
        if (jgeo.eq.'S') then
          dti=-(dlog(1.0d0-(rcin**3))/(alphahb*dht))
        else
          dti=-(dlog(1.0d0-rcin)/(alphahb*dht))
        endif
        dti=(rmax/3.d10)+dmax1(dti,rmax/3.d10)
        dthre=(1.2/dht)/3d-13
c
c     finite lifetime
c
        if (jeq.eq.'F') then
  630     write (*,640) rmax,dti
  640       format(/'  Source turn on and finite lifetime:'/
     &       ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' (Hydrogen Stromgrem radius of the order of',1pg10.3,
     &       ' cm'/' Characteristic time scale for ionisation: ',1pg7.1,
     &       ' sec)'/' Give elapsed time since source turned on (sec)'/
     &     ' and give life-time of ionising source (sec) : ',$)
          read (*,*) telap,tlife
          if ((telap.le.0.0d0).or.(tlife.le.0.0d0)) goto 630
          tm00=100.0d0
  650     write (*,660)
  660      format(//' Initial temperature of neutral gas : ',$)
          read (*,*,err=650) tm00
          if (tm00.lt.1.0d0) goto 650
        else
          tlife=0.0d0
  670     write (*,680) rmax,dthre
  680       format(/'  Source turn-off and final luminosity'/
     &       ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &      ' (Hydrogen Stromgrem radius of the order of',1pg10.3,' cm'
     &       /'Characteristic time scale for recombination: ',1pg7.1,
     &       ' sec)'/' Give elapsed time since source turned off (sec)'/
     &      ' and fractional final lumnosity of the source (>=0): ',$)
          read (*,*) telap,frlum
          if ((telap.le.0.0d0).or.(frlum.lt.0.0d0)) goto 670
        endif
      endif
c
  690 write (*,700)
  700 format(//' Give the step photon absorption fraction: ',$)
      read (*,*) dtau0
      if ((dtau0.le.0.0d0)) goto 690
c
c     standard convergence limits
c
      difma=0.05d0
      dtlma=0.10d0
      dhlma=0.050d0
c
      if (expertmode.gt.0) then
  710    format(/' Do you wish to alter the convergence criteria ?',
     &        ' (not recommended) : ',$)
        write (*,710)
        read (*,90) ilgg
        ilgg=ilgg(1:1)
        if (ilgg.eq.'n') ilgg='N'
        if (ilgg.eq.'y') ilgg='Y'
        if (ilgg.ne.'Y') ilgg='N'
c
        if (ilgg.eq.'Y') then
          write (*,720)
  720    format(/,' Set the convergence criteria, the maximum changes'
     &   ,'allowed',/
     &   ,' for three parameters: ion population, te and hydrogen'
     &   ,'density.',/
     &   ,' Convergence criteria: difma (0.05),dtlma (0.10)'
     &   ,'dhlma(0.015) : ',$)
c
          read (*,*) difma,dtlma,dhlma
        endif
c
c     end expertmode
c
      endif
c
c     ***PRINT SET UP
c
      if (jeq.eq.'E') then
        write (*,730) dtau0
  730    format(//' Model Summary :',/
     &          ,'    Mode  : Thermal and Ionic Equilibrium ',/
     &          ,'    dTau  :',0pf8.5)
      elseif (jeq.eq.'F') then
        write (*,740) telap,tlife,dtau0
  740 format(//' Model Summary :',/
     & ,'    Mode        :  Non-Equilibrium + Finite Source Life',/
     & ,'    Age         :',1pg10.3,' sec',/
     & ,'    Source Life :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      else
        write (*,750) telap,dtau0
  750 format(//' Model Summary :',/
     & ,'    Mode        :  Equilibrium + Source Switch off',/
     & ,'    Switch off  :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      endif
      if (usekappa) then
  760    format(/' Kappa Electron Distribution Enabled :'/
     &           '    Electron Kappa : ',1pg11.4)
        write (*,760) kappa
      endif
c
  770 format(/' Radiation Field and Parameters:',//,
     & t5, ' At estimated T_inner :',1pg10.3,' K'/
     & t5,' Rsou.',t20,' Remp.',t35,' Rmax',t50,' <DILU>'/
     & t5, 1pg10.3,t20,1pg10.3,t35,1pg10.3,t50,1pg10.3/
     & t5,' <Hdens>',t20,' <Ndens>',t35,' Fill Factor',/
     & t5, 1pg10.3,t20,1pg10.3,t35,0pf9.6,//
     & t5,' FQ tot',t20,' Log FQ',t35,' LQ tot',t50,' Log LQ',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' Q(N) in',t20,' Log Q(N)',t35,' <Q(N)>',t50,' Log<Q(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' Q(H) in',t20,' Log Q(H)',t35,' <Q(H)>',t50,' Log<Q(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' U(N) in',t20,' Log U(N)',t35,' <U(N)>',t50,' Log<U(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' U(H) in',t20,' Log U(H)',t35,' <U(H)>',t50,' Log<U(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3)
c
      write (*,770) tinner,rstar,remp,rmax,wdpl,dhn,dht,fin,qht,
     &dlog10(qht),astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),
     &qhdnav,dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),
     &unin,dlog10(unin),unav,dlog10(unav),uhin,dlog10(uhin),uhav,
     &dlog10(uhav)
c
c     ***TYPE OF EXIT FROM THE PROGRAM
c
  780 ielen=1
      jpoen=1
      tend=10.0d0
      fren=0.01d0
      diend=0.0d0
      tauen=0.0d0
c
c
  790 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Boundry Conditions  '/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,790)
c
  800 format(/'  Choose a model ending : '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Radiation bounded, HII <',2pf5.2,'%'/
     &     '    B  :   Ionisation bounded, XII < Y'/
     &     '    C  :   Temperature bounded, Tmin.'/
     &     '    D  :   Optical depth limited, Tau'/
     &     '    E  :   Density bounded, Distance'/
     &     '    F  :   Column Density limited, atom, ion'/
     &     '       :   '/
     &     '    R  :   (Reinitialise)'/
     &     '    G  :   (Reset geometry only)'/
     &     ' :: ',$)
      write (*,800) fren
      read (*,90) jend
c
      if ((jend(1:1).eq.'A').or.(jend(1:1).eq.'a')) jend='A'
      if ((jend(1:1).eq.'B').or.(jend(1:1).eq.'b')) jend='B'
      if ((jend(1:1).eq.'C').or.(jend(1:1).eq.'c')) jend='C'
      if ((jend(1:1).eq.'D').or.(jend(1:1).eq.'d')) jend='D'
      if ((jend(1:1).eq.'E').or.(jend(1:1).eq.'e')) jend='E'
      if ((jend(1:1).eq.'F').or.(jend(1:1).eq.'f')) jend='F'
      if ((jend(1:1).eq.'R').or.(jend(1:1).eq.'r')) jend='R'
      if ((jend(1:1).eq.'G').or.(jend(1:1).eq.'g')) jend='G'
c
      if ((jend.ne.'A').and.(jend.ne.'B').and.(jend.ne.'C')
     &.and.(jend.ne.'D').and.(jend.ne.'E').and.(jend.ne.'F')
     &.and.(jend.ne.'R').and.(jend.ne.'G')) goto 780
      if (jend.eq.'R') return
      if (jend.eq.'G') goto 580
      if ((jend.eq.'B').or.(jend.eq.'D')) then
  810   write (*,820)
  820    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 810
        ielen=zmap(ielen)
        jpoen=1
        if (arad(2,ielen).le.0.0d0) jpoen=2
      endif
c
      if (jend.eq.'B') then
  830   write (*,840) elem(ielen)
  840    format(/' Give the final ionisation fraction of ',a2,' : ',$)
        read (*,*) fren
        if ((fren.lt.0.0d0).or.(fren.gt.1.0d0)) goto 830
      endif
c
      if (jend.eq.'D') then
  850   write (*,860) elem(ielen)
  860 format(/,
     & ' Give the final optical depth at threshold of ',a2,':',$)
        read (*,*) tauen
        if (tauen.le.0.0d0) goto 850
      endif
c
      if (jend.eq.'C') then
  870   write (*,880)
  880    format(/' Give the final temperature (<10 as log): ',$)
        read (*,*) tend
        if (tend.lt.1.0d0) goto 870
      endif
c
      if (jend.eq.'E') then
  890   write (*,900)
  900    format(/' Give the distance or radius at which the density',
     &        ' drops: '/
     &        ' (in cm (>1E6) or as a fraction of the Stromgren'/
     &        ' radius (<1E6)) : ',$)
        read (*,*) diend
        if (diend.lt.1.d6) diend=diend*rmax
        diend=remp+diend
        if (diend.le.remp) goto 890
      endif
c
      if (jend.eq.'F') then
        write (*,910)
  910    format(/' Give the final column density (<100 as log): ',$)
        read (*,*) colend
        if (colend.lt.1.0d2) colend=10.d0**colend
  920   write (*,930)
  930    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 920
        ielen=zmap(ielen)
  940   write (*,950)
  950    format(/' Applies to ion stage : ',$)
        read (*,*) jpoen
        if (jpoen.le.0) goto 940
        if (jpoen.gt.maxion(ielen)) jpoen=maxion(ielen)+1
      endif
c
  960 format(//
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Output Requirements  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      write (*,960)
c
c
  970 format(//'  Choose output settings : '/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    A  :   Standard output (photnxxxx,phapnxxxx).'/
     & '    B  :   Standard + monitor 4 element ionisation.'/
     & '    C  :   Standard + all ions file.'/
     & '    D  :   Standard + final source + nu-Fnu spectrum.'/
     & '    F  :   Standard + first balance.'/
     & '    G  :   Everything'//
     & '    I  :   Standard + B + monitor 4 multi-level ion emission.'/
     & '    J  :   Standard + B + monitor up to 16 lines'/
     & ' :: ',$)
      write (*,970)
      read (*,90) ilgg
c
      if ((ilgg(1:1).eq.'A').or.(ilgg(1:1).eq.'a')) ilgg='A'
      if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg='B'
      if ((ilgg(1:1).eq.'C').or.(ilgg(1:1).eq.'c')) ilgg='C'
      if ((ilgg(1:1).eq.'D').or.(ilgg(1:1).eq.'d')) ilgg='D'
      if ((ilgg(1:1).eq.'F').or.(ilgg(1:1).eq.'f')) ilgg='F'
      if ((ilgg(1:1).eq.'G').or.(ilgg(1:1).eq.'g')) ilgg='G'
      if ((ilgg(1:1).eq.'I').or.(ilgg(1:1).eq.'i')) ilgg='I'
      if ((ilgg(1:1).eq.'J').or.(ilgg(1:1).eq.'j')) ilgg='J'
c
      jiel='NO'
      jiem='NO'
      if (ilgg.eq.'B') jiel='YES'
      if (ilgg.eq.'G') jiel='YES'
      if (ilgg.eq.'G') jiem='YES'
      if (ilgg.eq.'I') jiel='YES'
      if (ilgg.eq.'I') jiem='YES'
      if (ilgg.eq.'J') jiel='YES'
      if (ilgg.eq.'J') jlin='YES'
c
c
c     default monitor elements
c
      iel1=zmap(1)
      iel2=zmap(2)
      iel3=zmap(6)
      iel4=zmap(8)
c
      if (jiel.eq.'YES') then
  980    format(//' Give 4 elements to monitor (atomic numbers): ',$)
        write (*,980)
c
        read (*,*) iel1,iel2,iel3,iel4
        iel1=zmap(iel1)
        iel2=zmap(iel2)
        iel3=zmap(iel3)
        iel4=zmap(iel4)
c
      endif
c
      if (jiem.eq.'YES') then
c
  990 format(//' Choose 4 multi-level species #',i2,' :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,990) ,iem
        do i=1,(nfmions/5)+1
            j=(i-1)*5
          kmax=min(j+5,nfmions)-j
          write (*,'(5(2x,i2,": ",a2,a6))') (j+k,elem(fmatom(j+k)),
     &     rom(fmion(j+k)),k=1,kmax)
          enddo
 1000 format(' :: ',$)
          write (*,1000)
        read (*,*) idx1,idx2,idx3,idx4
        if (idx1.lt.1) idx1=1
        if (idx1.gt.nfmions) idx1=nfmions
        if (idx2.lt.1) idx2=1
        if (idx2.gt.nfmions) idx2=nfmions
        if (idx3.lt.1) idx3=1
        if (idx3.gt.nfmions) idx3=nfmions
        if (idx4.lt.1) idx4=1
        if (idx4.gt.nfmions) idx4=nfmions
        iemidx(1)=idx1
        iemidx(2)=idx2
        iemidx(3)=idx3
        iemidx(4)=idx4
      endif
c
      if (jlin.eq.'YES') then
c
 1100  format(//' Choose up to ',i2,' lines by wavelength (A)   :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
 1110  format(' Number of lines : ',$)
 1120  format(/' #',i2,' : ',$)
        write (*,1100) mxmonlines
        write (*,1110)
        read (*,*) nl
        nl=min(max(nl,1),mxmonlines)
        njlines=nl
        do i=1,njlines
            write (*,1120) i
            read (*,*) lam
            lam=dabs(lam)
            emlinlist(i)=lam
        enddo
        do i=1,mxmonlines
           emlindeltas(i)=0.001d0 ! will allow custom bins later
        enddo
        call speclocallineids(emlinlistatom,emlinlistion)
      endif
c
      jall='NO'
      if ((ilgg.eq.'C').or.(ilgg.eq.'G')) jall='YES'
c
      jspec='NO'
      if ((ilgg.eq.'D').or.(ilgg.eq.'G')) jspec='YES'
c
      jsou='NO'
      jpfx='psend'
c
      if ((ilgg.eq.'D').or.(ilgg.eq.'G')) then
c
        jsou='YES'
c
c     get final field file prefix
c
 1010   format (a16)
 1020   format(//' Give a prefix for final source file : ',$)
        write (*,1020)
        read (*,1010) jpfx
        nprefix=1
 1021 if ((jpfx(nprefix+1:nprefix+1).eq.' ')
     &     .or.(nprefix.ge.8)) goto 1022
          nprefix=nprefix+1
          goto 1021
 1022   continue
        jpfx=jpfx(1:nprefix)
      endif
c
      jbal='NO'
      jbfx='p6bal'
c
      if ((ilgg.eq.'F').or.(ilgg.eq.'G')) then
c
        jbal='YES'
c
c     get final field file prefix
c
 1030    format(//' Give a prefix for first ion balance file : ',$)
        write (*,1030)
        read (*,1010) jbfx
        nprefix=1
 1031 if ((jbfx(nprefix+1:nprefix+1).eq.' ')
     &     .or.(nprefix.ge.8)) goto 1032
          nprefix=nprefix+1
          goto 1031
 1032   continue
        jbfx=jbfx(1:nprefix)
      endif
c
c     get runname
c
 1040 format (a96)
 1050 format(//' Give a name/code for this run: ',$)
      write (*,1050)
      read (*,1040) runname
c
c
c     remains of old vax batch system (not used)
c
      banfil='INTERACTIVE'
c
      call compph6 (dhn, fin, banfil, difma, dtlma, dhlma)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c
c       PHOTO6: As for P5 but adds new integration scheme
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compph6 (dhn, fin, banfil, difma, dtlma, dhlma)
c
      include 'cblocks.inc'
c
      real*8 poppre(mxion, mxelem), popend(mxion, mxelem)
      real*8 ppre(mxion, mxelem),difma,dtlma,dhlma
      real*8 thist(3),dpw,dxp
      real*8 linfluxes(mxmonlines)
c
      real*8 aadn,agmax,agmu,aper,aperi,apero
      real*8 a,b,c,x,fx,fp,n0,r0,fv,e,f
      real*8 de,dedhma,dedhmi,deeq,chcksum
      real*8 dh,dhf,dhl,dhn,dhp,frp
      real*8 en,t_pk,qhratio
      real*8 f5007,f4363,f3727,f6300,f6584,f6720
c
      real*8 difg  , difmi
      real*8 difde , difte
      real*8 difpre, difend
c
      real*8 disin,dison,dlos0,dlos1,dlosav
      real*8 dr,drp,gcharge
      real*8 dtau,dtauon,dtaux,dtco,dtco0,dtco1
      real*8 dtoff,dton,dton0,dton1,durmax
      real*8 dv,dva,dvoluni,eqscal,exl,exla
      real*8 fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw,frx,fthick
      real*8 prescc,prf,protim,rad,radi,recscal,rww
      real*8 t,taux,te1a,tef,telast,tempre
      real*8 tep,tex,texm,tii0,tprop,trea,runit
      real*8 treach,treach0,treach1,tspr
      real*8 wd,wdil0,wdil1,wei
      real*8 ogas, oonh
      real*8 radpint, pres0, pres1
      real*8 popinttot
      real*8 a_v,ffuv,dnu,eng,habing,pahtest(mxinfph)
      real*8 pahsum, pahq0, pahqion, pahqfac
      real*8 q1,q2,q3,q4,q0
      real*8 def
c
      real*8 pop_initial (mxion, mxelem)
      real*8 pop_initial0(mxion, mxelem)
      real*8 pop_initial1(mxion, mxelem)
      real*8 te1_initial
      real*8 dh1_initial
      real*8 de1_initial
c
      integer*4 np,luop,lups,lusp,lusl,lupf,lunt,ir1,lulin
      integer*4 lut0,m,maxio,mma,n,nidh,niter
      integer*4 i,j,k,ndifatoms,atom
      integer*4 idx, nt, itr,init_it
c
      integer*4 luions(4),luemiss(4)
      integer*4 fuvmin,fuvmax,pahabsmax
c
      character jjd*4, pollfile*12,tab*4
      character newfil*32, filna*32, banfil*32, fnam*32, filnb*32
      character imod*4, lmod*4, nmod*4,wmod*4,ispo*4
      character linemod*4, spmod*4
      character fn*32
      character filnc*32, filnd*32, fspl*32, filnt*32, filin*32
      character filn(4)*32, filnem(4)*32
      character pfx*24,caller*4,sfx*4
      character irfile*32,chargefile*32
      logical iexi
c
c     External Functions
c
      real*8 fdilu,feldens,fradpress,fpressu,frectim,fzgas
      integer*4 lenv
c
c     internal functions
c
      real*8 dif,frad,fsight,fring
c
      dif(a,b)=dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      frad(x,n0,fx,fp,a,b,c,r0,fv,e,f)=n0*(fx*dexp((x-a)/r0)+fp*((x/a)**
     &b)+fv*(((r0-x)/e)**(-f)))+c
      fsight(aper,radi)=1.d0-(dcos(dasin(dmin1(1.d0,aper/(radi+(1.d-38*
     &aper)))))*(1.d0-(dmin1(1.d0,aper/(radi+(1.d-38*aper)))**2)))
      fring(aperi,apero,radi)=dmax1(0.d0,fsight(apero,radi)-
     &fsight(aperi,radi))
c
      jspot='NO'
      jcon='YES'
      ispo='SII'
      irfile=' '
      if (graindestructmode.eq.1) then
        grainmode=0
        grainadjust=1
      endif
      tab=char(9)
c
      ndifatoms=2
      q0=qhdhin
      te0=tinner
c
c      nt0 = time()
c      nt1 = 0
      mma=mxnsteps-1
      exla=5.d-5
      trea=1.d-5
      dton0=0.0d0
      dton1=0.0d0
      dtoff=0.0d0
      dtaux=0.d0
      drp=0.d0
      dedhmi=0.d0
      dedhma=0.d0
      treach0=0.0d0
      treach1=0.0d0
      texm=0.d0
      te0=1.d4
      te1=1.d4
      tep=1.d4
      dtco0=0.d0
      dtco1=0.d0
      dlos0=0.d0
      dlos1=0.d0
      thist(1)=0.d0
      thist(2)=0.d0
      thist(3)=0.d0
c
c     Define logical unit numbers
      luop=21
      lups=22
      lusp=23
      lupf=27
      lusl=28
      lunt=29
      lut0=0
      ir1=91
c
      ieln=4
      do i=1,ieln
        luions(i)=32+i
        luemiss(i)=42+i
      enddo
      lulin=52
c
c    ***DERIVE FILENAME : PHN**.PH6
c
      call p6fheader (newfil, banfil, fnam, filna, filnb, filnc, filnd,
     &filnt, filn, filnem, filin, luop, lupf, lusl, lusp, lups, lunt,
     &luions, luemiss, lulin, dhn, fin)
c
      chargefile=' '
c
c  Set up for PAH existence
c
      fuvmin=1
      fuvmax=1
      pahabsmax=1
      if ((grainmode.eq.1).and.(pahmode.ne.0)) then
        k=1
        do while (photev(k).le.6.d0)
          k=k+1
        enddo
c     FUVMIN = photon bin corresponding to 6 eV energy
        fuvmin=k-1
        do while (photev(k).le.iph)
          k=k+1
        enddo
c     FUVMAX = photon bin corresponding to 13.6 eV energy
        fuvmax=k-1
        do while (photev(k).le.262.012d0)
          k=k+1
        enddo
c     PAHABSMAX = photon bin corresponding to 262.0 eV energy
c       maximal PAH absorption bin
        pahabsmax=k-1
c
c Calculate photon flux in FUV range using Weingartner &
c   Draine (2001) ISRF
c
        ffuv=0.d0
        pahsum=0.d0
        do k=fuvmin,fuvmax
          eng=cphotev(k)
          dnu=(photev(k+1)-photev(k))*evplk
          pahsum=pahsum+pahnabs(k)*dnu
          if (eng.lt.9.26d0) then
            ffuv=ffuv+2.548d-18*eng**(-1.3322d0)/ev*pahnabs(k)*dnu
          elseif (eng.lt.11.2d0) then
            ffuv=ffuv+1.049d-16*eng**(-3.d0)/ev*pahnabs(k)*dnu
          else
            ffuv=ffuv+4.126d-13*eng**(-6.4172d0)/ev*pahnabs(k)*dnu
          endif
        enddo
        pahq0=ffuv/pahsum
        write (*,10) pahq0
   10 format(' ISRF photon flux :',1pg12.5)
      endif
c
c*********************************************************************
c
c   The routine does some initial setup, guesses at the parameters at
c   the inner boundary of the first spatial step, then for each spatial
c   step it guesses at the parameters at the outer boundary of the step,
c   iterates till these parameters converge, and then iterates over
c   spatial steps till the completion criteria are achieved.
c
c      write(*,*)'***COMPUTATION STARTS HERE'
c
c*********************************************************************
c
      fi=fin
      fi0=fin
      fi1=fin
c
      fidhcc=fin*dhn
c
      dv=1.0d6
      runit=(rmax-remp)/300.d0
      dr=0.d0
c
      if (dtau0.gt.0.075d0) dtau0=0.075d0
c
      dtau=dtau0*0.125d0
      difg=1.d0
c
      if (jgeo.eq.'S') then
        vunilog=3.0d0*dlog10(runit)
      else
        vunilog=dlog10(runit)
      endif
      if ((jeq.eq.'E').or.(jeq.eq.'P')) then
c     Equilibrium structure (E=equilibrium, P=post-equilibrium decay,
        nmod='EQUI'
      else
c     Equilibrium structure (F=finite source)
        nmod='TIM'
      endif
c
c  set convergence criterion parameters
c  agmax = ? ; aadn = ?
c
      if (jden.ne.'B') then
c        Density structure (B=isobaric)
        agmax=0.05d0
        aadn=2.5d0
      else
c     Density structure (C=isochroic, F=function, I=input file)
        agmax=0.80d0
        aadn=5.0d0
      endif
      agmu=agmax
c
      radpint=0.d0
      pres0=fpressu(tinner,dhn,pop0) ! uses proto-ionisation
c
c***********************************************
c
c      write(*,*)'***ITERATE, M = STEP NUMBER'
c
c***********************************************
c
      m=1
c
   20 jspot='NO'
      niter=0
      nidh=0
      difmi=1.d30
      wei=(dtau/dtau0)**0.75d0
      dtau=dtau0*wei
      if (dtau.gt.0.075d0) dtau=0.075d0
      wei=3.d0
      agmu=(agmax+(wei*agmu))/(1.d0+wei)
c
c****************************************************************
c
c      write(*,*)'IF M=1 , CHECK CONVERGENCE FOR THE DENSITY'
c      write(*,*)'AND IONIC POP. AT THE INNER BOUNDARY'
c
c****************************************************************
c
      if (m.eq.1) then
c
        te0=tinner
c
        if (graindestructmode.eq.1) te0=3.d5
c
        if (jden.eq.'B') dhn=(ponk/te0)/2.2d0 ! intial estimate
c
        pres0=fpressu(te0,dhn,pop0) ! uses proto-ionisation
c
        if (jthm.eq.'T') te0=fixtemp
c
        dr=0.d0
        drp=dr
        if (jden.ne.'F') then
c     Density structure (C=isochroic, B=isobaric, I=input file)
c     DH0 = Hydrogen density at inner spatial step boundary
          dh0=dhn
        else
c     Density structure (F=function)
          dh0=frad(remp,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)
        endif
c     Derive electron density using ionic populations in POP0.
c     DE0 = Electron density at the inner spatial step boundary
c     DEDHMA = max de/dh [total metallicity of the gas (relative to hydrogen)]
c     DEDHMI = min de/dh, allowing for ions prevented from recombining, usually
c             [total metallicity of gas for all species with
c                   -ve or 0 recombination rates from 2nd ion state??? - nope rs
        de0=feldens(dh0,pop0)
        dedhmi=0.d0
        dedhma=0.d0
        do j=1,atypes
          dedhma=dedhma+zion(j)
          if (arad(2,j).le.0.d0) dedhmi=dedhmi+zion(j)
        enddo
        call copypop (pop0, pop)
c
c     Derive geometrical dilution factor using photon source radius RSTAR
c     and radial distance from source of the inner step boundary RAD0.
        if (jgeo.eq.'S') then
          rad0=remp
          dis0=rad0
          wdil0=fdilu(rstar,rad0)
        endif
        if (jgeo.eq.'P') then
          rad0=0.d0
          dis0=remp
          wdil0=wdpl
        endif
c
c  Loop back here if hydrogen density change is too great
c
        init_it=0
c
   30   if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
        lmod='DW'
        call localem (te0, de0, dh0)
        call totphot2 (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
c
c     calculate qhdn into zone
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII,
c     q4=QATOT.
c
        call intvec (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
        if ((jgeo.eq.'P').and.(jden.eq.'B')) then
          qhratio=q0/qhdhin
          do j=1,infph
            soupho(j)=soupho(j)*qhratio
          enddo
          call localem (te0, de0, dh0)
          call totphot2 (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
          q1=0.0d0
          q2=0.0d0
          q3=0.0d0
          q4=0.0d0
          call intvec (tphot, q1, q2, q3, q4)
          qhdnin=q4/(dh0*zen)
          unin=qhdnin/cls
          qhdhin=q4/(dh0)
          uhin=qhdhin/cls
        endif
c
        if ((jgeo.eq.'S').and.(jden.eq.'B')) then
          qhratio=qhdhin/q0
          remp=remp*dsqrt(qhratio)
          rad0=remp
          dis0=rad0
          wdil0=fdilu(rstar,rad0)
          call localem (te0, de0, dh0)
          call totphot2 (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
          q1=0.0d0
          q2=0.0d0
          q3=0.0d0
          q4=0.0d0
          call intvec (tphot, q1, q2, q3, q4)
          qhdnin=q4/(dh0*zen)
          unin=qhdnin/cls
          qhdhin=q4/(dh0)
          uhin=qhdhin/cls
        endif
c
        if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 1)
c
        if ((jeq.eq.'E').or.(jeq.eq.'P').and.(jthm.eq.'S')) then
          dton0=1.d33
          call teequi2 (te0, te0, de0, dh0, dton0, nmod)
        else
          distim(1,0)=dis0
          distim(2,0)=dis0/cls
          treach0=distim(2,0)
          dton0=dmax1(0.d0,dmin1(tlife,telap-((2.d0*dis0)/cls)))
          jjd=jden
          jden='C'
          exl=-1.d0
          tex=0.95d0*tm00
          prescc=pres0+fradpress(0.0d0,dhn)
c
          call copypop (pop0, pop)
          call evoltem (tm00, te0, de0, dh0, prescc, dton0, exl, tex,
     &     lut0)
c
          jden=jjd
c
          if (dton0.le.0.d0) then
            write (*,40) dton0
   40         format(/,' Age was incorrectly set, it is shorter than',/,
     &            ' needed at the first space step , DTon0 :',1pg10.3,/)
            return
          endif
          dtco0=(1.5d0*fpressu(te0,dh0,pop))/(eloss+1.d-36)
        endif
c
        dhp=dh0
        tep=te0
        dlos0=dlos
c
        if (jden.eq.'C') then
          dh0=dhn
        elseif (jden.eq.'F') then
          dh0=frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)
        else
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c inner edge, rad pres = 0.0, iterate on nh only for
c pres0, used subsequently
c
c          dhn=ponk/te0
          pres0=fpressu(te0,dhn,pop) ! iterate to get initial press
          dhn=dhn*(ponk*rkb/pres0)
          dh0=dhn !*(prescc/(pres0))
        endif
c
        dhl=dif(dhp,dh0)
        write (*,'(1p3g14.7)') te0,dhl,remp
        init_it=init_it+1
        if ((dhl.ge.dhlma).or.
     &     ((init_it.lt.5).and.(jden.eq.'B'))) goto 30
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
        if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 1)
c
        tspr=dton0
c
        tempre=te0
        call copypop (pop, poppre)
        call copypop (pop, ppre)
c
c*****************************************************
c      write(*,*)'END OF INNER BOUNDRY FOR M = 1'
c*****************************************************
c
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII,
c     q4=QATOT.
c
        call intvec (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
      else
c
c     calculate qhdn into zone
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII,
c     q4=QATOT.
c
        call intvec (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh1*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh1)
        uhin=qhdhin/cls
c
      endif
c
c check grains, don't allow to destroy if already formed
c
      if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 0)
c
      zgas=fzgas()
      oonh=dlog10(zion(zmap(8))/zsol(zmap(8)))
      ogas=12.d0+dlog10(zion(zmap(8)))
c
c*****************************************************************
c
c      write(*,*)'***ESTIMATE OF VALUES AT THE OUTER BOUNDARY'
c
c*****************************************************************
c
      call copypop (poppre, pop)
      call copypop (poppre, popend)
c
      if (m.lt.4) then
c
c     linear extrap from previous 2 steps
c
        te1=te0+(((te0-tep)*te0)/tep)
        te1=dmax1(te0*agmu,dmin1(te0/agmu,te1))
        te1a=te1
      else
        te1a=te1
        te1=thist(3)
      endif
c
      if (jthm.eq.'T') te1=fixtemp
c
      if (te1.lt.0.d0) te1=100.d0
      if (jden.eq.'C') then
        dh1=dhn
      elseif (jden.eq.'F') then
        dh1=frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
c     error if pop is very different from equilibrium at te1
c         frp    = fradpress
c         dh1 = (dhn*prescc)/(fpressu(te1,dhn,pop)+frp)
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        pres1=fpressu(te1,dhn,pop)
        dh1=dhn*prescc/(pres1)
      endif
      de1=feldens(dh1,pop)
c
c***********************************************
c      write(*,*)'***DETERMINE SPACE STEP : DR'
c***********************************************
c
      te1_initial=te1
      dh1_initial=dh1
      de1_initial=de1
c
      call copypop (poppre, pop_initial0)
      call copypop (popend, pop_initial1)
      call copypop (pop, pop_initial)
c
   50 call copypop (pop_initial0, poppre)
      call copypop (pop_initial1, popend)
      call copypop (pop_initial, pop)
c
      t=(te0+te1)*0.5d0
      dh=(dh0+dh1)*0.5d0
      de=(de0+de1)*0.5d0
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      call equion (t, de, dh)
c
      lmod='SO'
      call totphot2 (t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
      call absdis2 (dh, fi, dtau, dr, rad0, pop)
c
c      dr = dmax1(dr,(1.0d-6*(rmax-remp)))
c
      if (jend.eq.'E') then
        if ((dis0+dr).gt.diend) dr=(diend-dis0)
      endif
      if (jgeo.eq.'S') then
        if (jden.eq.'F') then
          if ((frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)).gt.ofac) then
            dpw=dr
            if (pfac.gt.0.0d0) then
              if (bfac.gt.0.0d0) then
                dpw=dis0*((1.1d0**(1/bfac))-1.d0)
              elseif (bfac.lt.0.0d0) then
                dpw=dis0*((0.9d0**(1/bfac))-1.d0)
              endif
            endif
            dxp=dr
            if (xfac.gt.0.d0) then
              dxp=dabs(0.05d0*scalen)
            endif
            dr=1.d0/((1.d0/dpw)+(1.d0/dxp)+(1.d0/dr))
          endif
        endif
        rad1=rad0+dr
        dis1=rad1
        wdil1=fdilu(rstar,rad1)
        if (radout.lt.1.d38) then
          dvoluni=ftpi*((fring(radin,radout,rad1)*((rad1/runit)**3))-
     &     (fring(radin,radout,rad0)*((rad0/runit)**3)))
        else
          dvoluni=ftpi*(((rad1/runit)**3)-((rad0/runit)**3))
        endif
      endif
      if (jgeo.eq.'P') then
        dis1=dis0+dr
        rad1=0.d0
        wdil1=wdpl
        dvoluni=dr/runit
      endif
c
c**********************************************
c       write(*,*)'*** dr done',dr
c**********************************************
c
c******************************************************************
c
c      write(*,*)'CALCULATION OF IONIC POPULATION AT THE OUTER'
c     of current step
c
c******************************************************************
      if (jden.eq.'C') then
        dh1=dhn
      elseif (jden.eq.'F') then
        dh1=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        pres1=fpressu(te1,dhn,popend)
        dh1=dhn*prescc/(pres1)
      endif
c
      de1=feldens(dh1,popend)
c
      rww=1.d0
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
c
      t=(te0+(rww*te1))/(1.d0+rww)
      dh=(dh0+(rww*dh1))/(1.d0+rww)
      de=(de0+(rww*de1))/(1.d0+rww)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
c  Determine dust continuum emission for this region
c
      lmod='DW'
      wei=1.d0/(1.d0+rww)
      call averinto (wei, poppre, popend, pop)
      call localem (t, de, dh)
      call totphot2 (t, dh, fi, rad0, dr, dv, wdil1, lmod)
      call zetaeff (fi, dh)
c
c*****************************************************************
c      write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************
c
      if ((jeq.eq.'E').or.((jeq.eq.'P').and.(jthm.eq.'S'))) then
        dton1=1.d33
        call copypop (poppre, pop)
        call teequi2 (te1, tef, de1, dh1, dton1, nmod)
      else
c
c     Finite age, non-equilibrium ionisation. Determine time step.
c     Takes into account ionization front velocity.
c
        distim(1,m)=dis1
        distim(2,m)=distim(2,m-1)+(dr/(viofr+(1.d-35*dr)))
        tprop=distim(2,m)
        dtauon=1.0d0
c
        call absdis2 (dh, fi, dtauon, disin, 0.0d0, pop0)
        fthick=dmin1(disin,dis1-distim(1,0))
        dison=dis1-fthick
        protim=distim(2,0)
        do 60 j=m-1,0,-1
          if (distim(1,j).gt.dison) goto 60
          prf=(dison-distim(1,j))/((distim(1,j+1)-distim(1,j))+epsilon)
          protim=distim(2,j)+(prf*(distim(2,j+1)-distim(2,j)))
          goto 70
   60   continue
   70   continue
c
c
        treach1=dmax1(dis1/cls,protim+(fthick/cls))
        durmax=dmax1(0.d0,tlife+((dis1/cls)-treach1))
        dton1=dmax1(0.d0,dmin1(durmax,(telap-(dis1/cls))-treach1))
        deeq=dmax1(de1/3.d0,exla*dh1)
        eqscal=5.d0*frectim(2.d0*te1,deeq,dh1)
        if (dton1.lt.eqscal) then
c
c      *TIME STEP TOO SHORT , EQUILIBRIUM TEMPERATURE NOT REACHED
c
          jjd=jden
          jden='C'
          exl=-1.d0
          tex=0.95d0*tm00
          call copypop (pop0, pop)
          call evoltem (tm00, tef, de1, dh1, prescc, dton1, exl, tex,
     &     lut0)
          jden=jjd
        else
c
c      *TIME STEP LONG ENOUGH TO ASSUME EQUILIBRIUM TEMPERATURE
c
          call copypop (pop0, pop)
c     Compute equilibrium temperature and ionization state of the gas
          call teequi2 (te1, tef, de1, dh1, dton1, nmod)
        endif
        dtco1=(1.5d0*fpressu(tef,dh1,pop))/(eloss+epsilon)
      endif
c
c*****************************************************************
c       write(*,*)'END EQUILIBRIUM IONISATION AND TEMPERATURE .'
c       write(*,*) t,te1,tef
c*****************************************************************
c
      if (jthm.eq.'T') tef=fixtemp
c
      if (jden.eq.'C') then
        dhf=dhn
      elseif (jden.eq.'F') then
        dhf=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        pres1=fpressu(tef,dhn,pop)
        dhf=dhn*prescc/(pres1)
c       write(*,'(" P/k P_H/k:",2(1pe11.4))')
c     &          fpressu(tef,dhf,pop)/rkb, dhf*tef
      endif
      def=feldens(dhf,pop)
c
c************************************************************
c   write(*,*)'***COMPARES END VALUES WITH PREVIOUS STEP'
c************************************************************
c
      call difhhe (pop, popend, difend)
      call difhhe (pop, poppre, difpre)
c      call difpop(pop, popend, 1e-8, 5, difend)
c      call difpop(pop, poppre, 1e-8, 5, difpre)
      call copypop (pop, popend)
      difend=difend/difma
      difpre=difpre/difma
      difde=dif(def,de1)/dhlma
      difte=dif(tef,te1)/dtlma
c
      difg=dmax1(difpre,difend,difde,difte)
c
      niter=niter+1
      if (difg.le.difmi) then
        texm=te1
        dtaux=dtau
        difmi=difg
      endif
c
      if (niter.gt.1) then
        write (*,80) tef,difte,difde,difpre,difg,difend,dtau,dr
   80  format(t7,0pf8.0,5(f8.3),3(1pg9.2))
      endif
c
c***************************************************************
c  write(*,*)'AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED'
c***************************************************************
c
c
      rww=1.d0
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav=(dis0+(rww*dis1))/(1.d0+rww)
      wei=1.d0/(1.d0+rww)
      call averinto (wei, poppre, popend, pop)
c
      t_pk=fpressu(t,dh,pop)
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'P') then
        wdil=wdpl
        rad=0.d0
      endif
c
c****************************************************************
c  write(*,*)'CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP'
c****************************************************************
c
c
      if ((difg.ge.1.d0).and.(difend.gt.0.001d0)) then
        te1=tef
        if (niter.gt.6) then
c     No convergence after 6 iterations
c            write(*, 474)  tef,difte,difde,difpre,difg,difend,dtau,dr
c 474        format(' FINAL',0pf8.0,5(f8.3),3(1pg9.2))
          goto 90
        endif
c
        dtau=dtau*0.3333d0
c
        te1=te1_initial
        dh1=dh1_initial
        de1=de1_initial
c
        goto 50
c     Jump here if no spatial convergence.
   90   continue
      endif
c
      te1=tef
      dlos1=dlos
      telast=t
      frx=1.d0-pop(jpoen,ielen)
c
c*****************************************************************
c      write(*,*)'END CHECK CONVERGENCE STEP'
c*****************************************************************
c
c****************************************************************
c      write(*,*)'INTEGRATES COLUMN DENSITY AND END DIFFUSE FIELD'
c****************************************************************
c
      if (jtrans.eq.'OUTW') then
        frdw=1.d0
      endif
      if (jtrans.eq.'LOUP') then
        frdw=0.5d0
      endif
      if (jtrans.eq.'LODW') then
        frdw=0.5d0
      endif
      if (jtrans.eq.'DWUP') then
        frdw=0.5d0
      endif
c
      lmod='DW'
      call localem (t, de, dh)
      call totphot2 (t, dh, fi, rad0, dr, dv, wdil1, lmod)
      call zetaeff (fi, dh)
      if ((grainmode.eq.1).and.(irmode.ne.0)) then
        call dusttemp (t, dh, de, fi, dr, m)
      endif
      call newdif2 (t, t, dh, fi, rad0, dr, dv, 0.d0, dv, frdw, jtrans)
c
      imod='COLD'
      call sumdata (t, de, dh, fi, dvoluni, dr, disav, imod)
c
c*****************************************************************
c      write(*,*)'CALCULATION OF AVERAGE QUANTITIES '
c       EQUAL VOLUMES
c*****************************************************************
c
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
      en=dh*zen
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav=(dis0+(rww*dis1))/(1.d0+rww)
      dlosav=(dlos0+(rww*dlos1))/(1.d0+rww)
      dton=(dton0+(rww*dton1))/(1.d0+rww)
      treach=(treach0+(rww*treach1))/(1.d0+rww)
      dtco=(dtco0+(rww*dtco1))/(1.d0+rww)
      wei=1.d0/(1.d0+rww)
c
      call averinto (wei, poppre, popend, pop)
c
      t_pk=fpressu(t,dh,pop)/rkb
      frp=fradpress(dr,dh)
      radpint=radpint+frp
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'P') then
        wdil=wdpl
        rad=0.d0
      endif
c
c*************************************************************
c    ***CALCULATION OF IONIC POPULATION AND TEMPERATURE WHEN
c       OR IF SOURCE SWITCHED OFF
c*************************************************************
c
      dtoff=dmax1(0.d0,(telap-((2.d0*disav)/cls))-tlife)
      if ((jeq.ne.'E').and.(dtoff.gt.0.0)) then
c
        if (jeq.eq.'P') recscal=frectim(t,0.8d0*de,dh)
c
        wd=frlum*wdil
        lmod='SO'
        jspot='YES'
        tii0=t
c
        call totphot2 (t, dh, fi, rad, dr, dv, wd, lmod)
        call zetaeff (fi, dh)
        call evoltem (tii0, t, de, dh, prescc, dtoff, exla, tm00, lut0)
c
        if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      endif
c*********************************************************
c
c      write(*,*)'COMPUTES SPECTRUM OF THE REGION CONSIDERED'
c
c*********************************************************
c
      imod='REST'
c     Calculate total cooling rate of the plasma
      call cool (t, de, dh)
c     Integrate the line and temperature data and column density
      call sumdata (t, de, dh, fi, dvoluni, dr, disav, imod)
c
c     record line ratios
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
c    Test for PAH existence (if grainmode true)
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
c
c  Calculate habing photodissociation parameter
c
        ffuv=0.d0
        do i=fuvmin,fuvmax-1
          dnu=(photev(i+1)-photev(i))*evplk
          eng=cphote(i)
c    photons s-1 cm-2
          ffuv=ffuv+fpi*tphot(i)*dnu/eng
        enddo
        habing=ffuv/(cls*dh)
c
c Calcualte PAH ionization/Local ISRF ratio (Q factor)
c
        pahsum=0.d0
        ffuv=0.d0
        do i=fuvmin,pahabsmax-1
          dnu=(photev(i+1)-photev(i))*evplk
          eng=cphote(i)
          pahsum=pahsum+pahnabs(i)*dnu
          ffuv=ffuv+fpi*tphot(i)/eng*pahnabs(i)*dnu
        enddo
        pahqion=ffuv/pahsum
        pahqfac=pahqion/pahq0
c
c      write(*,210) Hab,pahQfac
c 210  format('        HP =',1pg12.5,', PAH ion/ISRF =',1pg12.5,/)
c
c Determine if PAHs exist
c If so, deplete gas of Carbon in PAHs
c
        if (pahactive.eq.0) then
c
          if (pahend.eq.'H') then
            if (habing.lt.pahlimit) then
              pahactive=1
            endif
          endif
          if (pahend.eq.'Q') then
            if (qhdh.lt.pahlimit) then
              pahactive=1
            endif
          endif
          if (pahend.eq.'I') then
            if (pahqfac.lt.pahlimit) then
              pahactive=1
            endif
          endif
c
        endif
c
        if (pahactive.eq.1) then
          atom=zmap(6)
          if (clinpah.eq.1) then
            zion(atom)=zion0(atom)*deltazion(atom)*dion(atom)
          else
            zion(atom)=zion0(atom)*deltazion(atom)*(dion(atom)+(1.d0-
     &       dion(atom))*pahcfrac)
          endif
        endif
c
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Special output controlled by poll files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      fhi=pop(1,1)
      fhii=pop(2,1)
c
      pollfile='balance'
      inquire (file=pollfile,exist=iexi)
      if (((jbal.eq.'YES').and.(m.eq.1)).or.(iexi)) then
        caller='P6'
        pfx=jbfx//' '
        np=lenv(pfx)
        call wbal (caller, pfx, np, pop)
      endif
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        caller='P6'
        pfx='psou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
      endif
c
      pollfile='nebphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        lmod='NEBL'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P6'
        pfx='nsou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='sophot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        lmod='SO'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P6'
        pfx='ssou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='loclphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        lmod='LOCL'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P6'
        pfx='lsou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot2 (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='speclocal'
      inquire (file=pollfile,exist=iexi)
      if (iexi) then
        caller='P6'
        pfx='local'
        np=5
        sfx='lam'
        call newfile (pfx, np, sfx, 3, fn)
        fspl=fn(1:(np+3+5))
        open (lusp,file=fspl,status='NEW')
        linemod='LAMB'
        spmod='REL'
        call speclocal (lusp, tloss, eloss, egain, dlos, t, dh, de,
     &   pop(1,1), disav, dr, linemod, spmod)
        close (lusp)
      endif
c
      if (grainmode.eq.1) then
        pollfile='IRphot'
        inquire (file=pollfile,exist=iexi)
        if ((iexi).and.(irmode.gt.0)) then
          caller='P6'
          pfx='irsou'
          np=5
          wmod='REAL'
          dva=0
          do i=1,infph-1
            pahtest(i)=0.5d0*(photev(i+1)+photev(i))*ev*irphot(i)*dr
          enddo
          call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0,
     &     qhdh, disav, dr, 0.d0, 0.d0, 1.d0, pahtest)
        endif
c
        pollfile='pahphot'
        inquire (file=pollfile,exist=iexi)
        if ((iexi).and.(pahactive.eq.1)) then
          caller='P6'
          pfx='pahs'
          np=4
          wmod='REAL'
          dva=0
          do i=1,infph
            pahtest(i)=paheng*pahflux(i)
          enddo
          call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0,
     &     qhdh, disav, dr, 0.d0, 0.d0, 1.d0, pahtest)
        endif
c
c   Infrared output (every 1 step)
c
        pollfile='IRlocal'
        inquire (file=pollfile,exist=iexi)
        if ((irmode.ne.0).and.iexi) then
          if (irfile.eq.' ') then
            caller='P6'
            np=6
            pfx='IRflux'
            sfx='sou'
            call newfile (pfx, np, sfx, 3, fn)
            irfile=fn(1:(np+3+5))
            open (ir1,file=irfile,status='NEW')
            write (ir1,100) theversion
  100      format('%  Infrared flux per region ',/,
     &            '%  MAPPINGS V ',a8,/,
     &            '%  given as:',/,
     &            '%    energy edge (eV) ,,continuumflux, ',
     &            'IRflux Fnu(erg s-1 cm-2 Hz-1Sr-1)',/,
     &            '%')
            close (ir1)
          endif
          do i=1,infph-1
            cnphot(i)=0.d0
          enddo
          call freefree (t, de, dh)
          call freebound (t, de, dh)
          call twophoton (t, de, dh)
          open (ir1,file=irfile,status='OLD',access='APPEND')
          write (ir1,110) m
          write (ir1,*) infph
  110     format(/,' Region ',i4.4)
          do i=1,infph-1
            cnphot(i)=ffph(i)+fbph(i)+p2ph(i)
            write (ir1,120) photev(i),cnphot(i)*dr,irphot(i)*dr
          enddo
  120     format(1pg14.7,' ',1pg14.7,' ',1pg14.7)
          close (ir1)
        endif
c
      endif
c
c     graincharge output
c
      if (grainmode.eq.1) then
        pollfile='graincharge'
        inquire (file=pollfile,exist=iexi)
        if (iexi) then
          if (chargefile.eq.' ') then
            caller='P6'
            np=5
            pfx='grpot'
            sfx='ph6'
            call newfile (pfx, np, sfx, 3, fn)
            chargefile=fn(1:(np+3+5))
            open (ir1,file=chargefile,status='NEW')
            write (ir1,130) theversion
  130      format('%  grain charge in each region ',/,
     &            '%  MAPPINGS V ',a8,/,
     &            '%  given as:',/,
     &            '%  m,Av. distance,Temp,dh,de,',/,
     &            '%  grain charge(allsizes) (graphite),',/,
     &            '%  grain charge(allsizes) silicate',/,
     &            '%')
            write (ir1,140) (grainrad(i),i=mindust(1),maxdust(1))
  140      format('gra radii',10(1pg11.4))
            write (ir1,*)
            write (ir1,150) (grainrad(i),i=mindust(2),maxdust(2))
  150      format('sil radii',10(1pg11.4))
            close (ir1)
          endif
          open (ir1,file=chargefile,status='OLD',access='APPEND')
          write (ir1,160) m,disav,t,dh,de
          write (ir1,170) (grainpot(i,1),i=mindust(1),maxdust(1))
          write (ir1,180) (grainpot(i,2),i=mindust(2),maxdust(2))
          close (ir1)
  160      format(/,i3,4(1pg14.5))
  170      format('gra ',30(1pg11.4))
  180      format('sil ',30(1pg11.4))
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Main output
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open (luop,file=filna,status='OLD',access='APPEND')
      open (lupf,file=filnb,status='OLD',access='APPEND')
      open (lups,file=filnc,status='OLD',access='APPEND')
      open (lunt,file=filnt,status='OLD',access='APPEND')
c
      if (jall.eq.'YES') then
        open (lusl,file=newfil,status='OLD',access='APPEND')
      endif
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
c
  190 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7),6(1x,1pg10.3))
c
        write (luop,190) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),
     &   hoiii(m),graindgr,zgas,ogas,habing,pahqfac
c
      elseif (graindestructmode.eq.1) then
c
  200 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7),4(1x,1pg10.3))
c
        write (luop,200) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),
     &   hoiii(m),graindgr,zgas,ogas
c
      else
c
  210 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7))
c
        write (luop,210) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls)
c
      endif
c
  220 format(i4,16(1pg12.5,', '))
      write (lunt,220) m,(dis0-remp),dr,disav,t,dh,de,de+en,fhi,fhii,
     &hbeta,f5007,f4363,f3727,f6300,f6584,f6720
c
      if (jeq.eq.'F') write (luop,230) fthick,treach,dtco,dton,dtoff
  230 format(t71,5(1pg9.2))
c
  240 format(t98,5(1pg9.2))
      if (jeq.eq.'P') write (luop,240) recscal,dtoff
c
      if (mod(m,20).eq.1) then
c
        if ((graindestructmode.eq.0).and.(grainmode.eq.1)
     &   .and.(pahmode.eq.1)) then
c
          if (jgeo.eq.'S') then
  250  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103,'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,'<Charge>')
            write (*,250)
          else
  260  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103,'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,'<Charge>')
            write (*,260)
          endif
c
        elseif ((graindestructmode.eq.1).and.(pahmode.eq.1)) then
c
          if (jgeo.eq.'S') then
  270  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,' Zgas',t159,'LogO+12',t167,'<Charge>')
            write (*,270)
          else
  280  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,' Zgas',t159,'LogO+12',t167,'<Charge>')
            write (*,280)
          endif
c
        elseif ((graindestructmode.eq.1).and.(pahmode.eq.0)) then
c
          if (jgeo.eq.'S') then
  290  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,' Zgas',t143,' LogO+12s'
     & ,       t150,'<Charge>')
            write (*,290)
          else
  300  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,' Zgas',t143,' LogO+12s'
     & ,       t150,'<Charge>')
            write (*,300)
          endif
c
        else
c
          if (jgeo.eq.'S') then
  310  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k', t103, 'Log<QN>', t111,'Log<UH>'
     & ,       t119,'OIII/HB')
            write (*,310)
          else
  320  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<FHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k', t103, 'Log<QN>', t111,'Log<UH>'
     & ,       t119,'OIII/HB')
            write (*,320)
          endif
c
        endif
c
      endif
c
c 465  format('    Case A-B (HI, HeII): ',2(1pg11.3))
c      write(*,465) caseab(1), caseab(2)
c
  330 format(i4,x,0pf8.0,(1pg9.2),0pf6.3,0pf6.3
     & ,6(1pg10.3),5(0pf8.3),2(0pf7.3,x),6(0pf8.3))
c
      gcharge=0.5d0*(avgpot(1)+avgpot(2))
c
      if ((graindestructmode.eq.0).and.(grainmode.eq.1)
     &.and.(pahmode.eq.1)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),f5007/hbeta,
     &   graindgr,habing,pahqfac,gcharge
c
      elseif ((graindestructmode.eq.1).and.(pahmode.eq.1)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),
     &   dlog10(qhdn/cls),f5007/hbeta,graindgr,habing,pahqfac,zgas,ogas,
     &   gcharge
c
      elseif ((graindestructmode.eq.1).and.(pahmode.eq.0)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),
     &   dlog10(qhdn/cls),f5007/hbeta,graindgr,zgas,ogas,gcharge
c
      else
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),f5007/hbeta
c
      endif
c
      if (jiel.eq.'YES') then
        do i=1,ieln
          open (luions(i),file=filn(i),status='OLD',access='APPEND')
c
  340  format(1x,4(1pg12.5,', '),31(1pg12.5,', '))
          if (i.eq.1) then
            write (luions(i),340) disav,(dis0-remp),dr,t,(pop(j,iel1),j=
     &       1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),340) disav,(dis0-remp),dr,t,(pop(j,iel2),j=
     &       1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),340) disav,(dis0-remp),dr,t,(pop(j,iel3),j=
     &       1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),340) disav,(dis0-remp),dr,t,(pop(j,iel4),j=
     &       1,maxion(iel4))
          endif
c
          close (luions(i))
        enddo
c
      endif
c
      if (jiem.eq.'YES') then
        do i=1,ieln
          open (luemiss(i),file=filnem(i),status='OLD',access='APPEND')
          idx=iemidx(i)
          nt=nfmtrans(idx)
c
  350  format(1x,i4,4(', ',1pg12.5),31(', ',1pg12.5))
c
          write (luemiss(i),350) m,(dis0-remp),dr,disav,t,hbeta,
     &     (fmbri(itr,idx),itr=1,nt)
c
          close (luions(i))
        enddo
      endif
c
      if (jlin.eq.'YES') then
          call speclocallines(linfluxes,emlinlist,njlines)
          open(lulin,file=filin,status='OLD',access='APPEND')
  355     format(1x,i4,5(', ',1pg12.5),16(',    ',1pg12.5))
          write (lulin,355) m,(dis0-remp),dr,disav,t,hbeta,
     &     (linfluxes(itr),itr=1,njlines)
          close(lulin)
      endif
c
      if (jall.eq.'YES') then
        write (lusl,360)
        write (lusl,370) m,t,de,dh,fi,dvoluni,disav,(dis0-remp),dr
  360    format(//' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  fill. Fact. ',
     &                '  dVol Unit.  ',
     &                '    dR (cm)   ',
     &                ' Dist.Ave.(cm)',
     &                ' Din-Remp.(cm)')
  370    format(i3,8(1pg14.7)/)
      endif
c
      wmod='PROP'
      call wmodel (lupf, t, de, dh, disav, wmod)
      wmod='LOSS'
      call wmodel (lups, t, de, dh, disav, wmod)
c
      if (jall.eq.'YES') then
        call wionabal (lusl, pop)
        call wionabal (lusl, popint)
      endif
c
      close (lunt)
      close (luop)
      close (lupf)
      if (jall.eq.'YES') close (lusl)
      close (lups)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***END OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***TEST ENDING CONDITIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) goto 400
      if ((jend.eq.'C').and.(telast.le.tend)) goto 400
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 400
c
      do 380 n=1,ionum
c      find first (lowest energy) cross section that matches species
        if ((atpho(n).eq.ielen).and.(ionpho(n).eq.jpoen)) goto 390
  380 continue
c sigpho contains new threshold values if verner in force.
  390 taux=sigpho(n)*popint(jpoen,ielen)
c      write(*,*) n,jpoen,ielen,sigpho(n)
      if ((jend.eq.'D').and.(taux.ge.tauen)) goto 400
      if (jend.eq.'F') then
        if (jpoen.gt.maxion(ielen)) then
          popinttot=0
          do i=1,maxion(ielen)
            popinttot=popinttot+popint(i,ielen)
          enddo
          if (popinttot.ge.colend) goto 400
        else
          write (*,*) jpoen,ielen,popint(jpoen,ielen),colend
          if (popint(jpoen,ielen).ge.colend) goto 400
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Check Visual extinction condition
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (jend.eq.'H') then
        a_v=dustav*dustint
        if (a_v.gt.a_vend) goto 400
      endif
c
      if (((((de1/dh1)-dedhmi)/dedhma).le.exla).and.(grainmode.eq.0))
     &goto 400
c
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 400
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (popend, poppre)
c
      tep=te0
      te0=tef
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Record electron temperatures in THIST for quadratic extrapolation
c     when m > 3.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      thist(1)=thist(2)
      thist(2)=thist(3)
      thist(3)=tef
c
      de0=def
      dh0=dhf
      dton0=dton1
      treach0=treach1
      dtco0=dtco1
      fi0=fi1
      rad0=rad1
      dis0=dis1
      wdil0=wdil1
      drp=dr
      dlos0=dlos1
c
      if (difg.lt.0.1d0) dtau=0.1d0*dtau/difg
      if (difg.gt.0.75d0) dtau=dtau/(1.d0+(difg-0.75d0)*5.d0)**2.d0
      dtau=dmax1(1.0d-5,dtau)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      m=m+1
c
c     stop things getting silly....
c
      if (m.gt.mma) goto 400
c
c     loop back
c
      goto 20
c
  400 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write(*,*)'MODEL ENDED ; OUTPUT RESULTS **************'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      open (luop,file=filna,status='OLD',access='APPEND')
      open (lupf,file=filnb,status='OLD',access='APPEND')
      open (lunt,file=filnt,status='OLD',access='APPEND')
      write (luop,410)
      write (lunt,410)
  410 format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &'TAUXF',t58,'Ending')
      write (luop,420) t,dis1,frx,taux,jend
  420 format(' ###########' ,t18,0pf8.0,3(1pg10.3),4x,a4,//)
      if ((jeq.ne.'E').and.(frlum.gt.0.0)) write (luop,430) frlum
  430 format(/' Final fractional luminosity of source after ',
     &'turn off :',1pg10.3)
      if ((radout.lt.huge).and.(jgeo.eq.'S')) write (luop,440) radin,
     &radout
  440 format(/' NB. :::::::::::::::Integration through the line',
     &'of sight for a ring aperture of radii :',2(1pg10.3))
      write (lupf,450) tempre,tspr
  450 format(//t4,'Preionisation conditions for step#0 ',
     &'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
      write (lupf,460) (elem(i),i=1,atypes)
  460 format(' ',t8,16(4x,a2,4x)/)
      maxio=0
      do j=1,atypes
        if (maxio.le.maxion(j)) maxio=maxion(j)
      enddo
      do j=1,maxio
        chcksum=0.0d0
        do i=1,atypes
          chcksum=chcksum+ppre(j,i)
        enddo
        if (chcksum.gt.0.0d0) then
          write (lupf,470) rom(j),(ppre(j,i),i=1,atypes)
  470    format(' ',a6,t8,16(1pg10.3))
        endif
      enddo
c
      if (jsou.eq.'YES') then
c
        caller='P6'
        pfx='spec_tot'
        np=8
        wmod='LFLM'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
        caller='P6'
        pfx=jpfx//' '
        np=lenv(pfx)
        wmod='NFNU'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c     Down stream photon field
c
        caller='P6'
        pfx=jpfx//' '
        np=lenv(pfx)
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c
c      write nebula spectrum in lam - flam (ergs/s/cm2/A)
c
        lmod='NEBL'
        call totphot2 (te1, dh1, fi, dis0, dr, 0.0d0, wdil0, lmod)
c
        caller='P6'
        pfx='spec_neb'
        np=8
        wmod='LFLM'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
      endif
c
      if (jspec.eq.'YES') then
c also do nu - nufnu
        caller='P6'
        pfx=jpfx
        np=5
        wmod='NFNU'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c     Down stream nebula only photon field
c  and sou
        caller='P6'
        pfx=jpfx
        np=5
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
      endif
c
c    ***COMPUTE AVERAGE SPECTRUM AND WRITES IT IN FILE PHN
c
      close (lupf)
c
      call avrdata
      call wrsppop (luop)
      close (luop)
      close (lunt)
c
      linemod='LAMB'
      spmod='REL'
c      call spec2 (lunt, linemod, spmod)
c
      open (lusp,file=filnd,status='OLD',access='APPEND')
      call spec2 (lusp, linemod, spmod)
      close (lusp)
c
      write (*,480) fnam
  480 format(//' Output created &&&&&& File : ',a/)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to get the file header buisness out where
c     it can be worked on, and/or modified for other headers.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine p6fheader (newfil, banfil, fnam, filna, filnb, filnc,
     &filnd, filnt, filn, filnem, filin, luop, lupf, lusl, lups,
     &lusp, lunt, luions, luemiss, lulin, dhn, fin)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 dhn,fin
      real*8 dht,zi(mxelem)
      integer*4 i, ie, j, at, io
      integer*4 idx, plen, nt, nl, itr
      integer*4 luop,lups,lupf,lusl,lusp,lunt,lulin
      integer*4 luions(4), luemiss(4)
      character fn*32
      character abundtitle*24
      character pfx*16,sfx*4
      character newfil*32, filna*32, banfil*32, fnam*32, filnb*32
      character filin*32
      character filnc*32, filnd*32, filnt*32, filn(4)*32, filnem(4)*32
c
      real*8 densnum
c
      dht=densnum(dhn)
c
      fn=' '
      pfx='photn'
      sfx='ph6'
      call newfile (pfx, 5, sfx, 3, fn)
      filna=fn(1:13)
      fnam=filna
c
c
      fn=' '
      pfx='phapn'
      sfx='ph6'
      call newfile (pfx, 5, sfx, 3, fn)
      filnb=fn(1:13)
c
      fn=' '
      pfx='phlss'
      sfx='ph6'
      call newfile (pfx, 5, sfx, 3, fn)
      filnc=fn(1:13)
c
      fn=' '
      pfx='phsem'
      sfx='ph6'
      call newfile (pfx, 5, sfx, 3, fn)
      filnt=fn(1:13)
c
      fn=' '
      pfx='emiss'
      sfx='csv'
      call newfile (pfx, 5, sfx, 3, fn)
      filin=fn(1:13)
c
      fn=' '
      pfx='spec'
      sfx='csv'
      call newfile (pfx, 4, sfx, 3, fn)
      filnd=fn(1:12)
c
      if (jall.eq.'YES') then
        fn=' '
        pfx='allion'
        sfx='ph6'
        call newfile (pfx, 6, sfx, 3, fn)
        newfil=fn(1:14)
      endif
c
      open (luop,file=filna,status='NEW')
      open (lupf,file=filnb,status='NEW')
      open (lups,file=filnc,status='NEW')
      open (lusp,file=filnd,status='NEW')
      open (lunt,file=filnt,status='NEW')
      if (jlin.eq.'YES') then
         open (lulin,file=filin,status='NEW')
      endif
      if (jall.eq.'YES') then
         open (lusl,file=newfil,status='NEW')
      endif
c
c     ***WRITES ON FILE INITIAL PARAMETERS
c
      write (luop,10) theversion,runname,banfil,fnam
      write (lupf,10) theversion,runname,banfil,fnam
      write (lups,10) theversion,runname,banfil,fnam
      write (lusp,10) theversion,runname,banfil,fnam
      write (lunt,10) theversion,runname,banfil,fnam
      if (jall.eq.'YES') then
        write (lusl,10) theversion,runname,banfil,fnam
      endif
      if (jlin.eq.'YES') then
        write (lulin,10) theversion,runname,banfil,fnam
      endif
   10 format(
     &'::::::::::::::::::::::::::::::::::::'
     &,'::::::::::::::::::::::::::::::::::::',/
     & ' Photoionisation model P6, MAPPINGS V ',a8,/
     &'::::::::::::::::::::::::::::::::::::'
     &,'::::::::::::::::::::::::::::::::::::',/
     &' (Diffuse Field , Robust Integration, Radiation Pressure)',//
     &,' Run   :, ',a96/
     &,' Input :, ',a32/
     &,' Output:, ',a32/)
      write (lupf,20)
   20 format(' This file contains the plasma properties for ',
     &'each model step,'/
     &' and the pre-ionisation conditions used.'/ )
   30 format(' This file contains the ionisation fraction for ',
     &'all atomic elements as a function of distance,'/
     &' and the abundances used. See photn for model details.'/ )
      if (jall.eq.'YES') then
        write (lusl,30)
      endif
      write (lunt,40)
   40 format(' This file contains the nebula structure summary,'/
     &' and the strong line emissivities. See photn for model details/')
   45 format(' This file line emissivities (erg/cm^3/s/sr) for up to '
     & ,i2,' selected lines as a function of distance.',
     & ' See photn file for model details/')
      if (jlin.eq.'YES') then
        write (lulin,45) mxmonlines
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     abundances file header
c
c
      abundtitle=' Initial Abundances :'
      do i=1,atypes
        zi(i)=zion0(i)*deltazion(i)
      enddo
c
      call dispabundances (luop, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
c
      call wabund (lupf)
      call wabund (lups)
      call wabund (lusp)
      call wabund (lunt)
      if (jall.eq.'YES') then
        call wabund (lusl)
      endif
c
      close (lusp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jeq.eq.'E') then
        write (luop,50) dtau0
   50 format( ' Model Summary :',//
     &,t5,' Mode  : Thermal and Ionic Equilibrium ',/
     &,t5,' dTau  :',0pf8.5,/)
      elseif (jeq.eq.'F') then
        write (luop,60) telap,tlife,dtau0
   60 format( ' Model Summary :',//
     &,t5,' Mode        :  Non-Equilibrium + Finite Source Life ',/
     &,t5,' Age         :',1pg10.3,' sec',/
     &,t5,' Source Life :',1pg10.3,' sec',/
     &,t5,' dTau        :',0pf8.5,/)
      else
        write (luop,70) telap,dtau0
   70 format( ' Model Summary :',//
     &,t5,' Mode        :  Equilibrium + Source Switch off ',/
     &,t5,' Switch off  :',1pg10.3,' sec',/
     &,t5,' dTau        :',0pf8.5,/)
      endif
      if (usekappa) then
   80 format( ' Kappa Electron Distribution Enabled :',//
     & ,t5,' Electron Kappa : ',1pg11.4,/)
        write (luop,80) kappa
      endif
c
   90 format(/' Radiation Field and Parameters:',//,
     & t5, ' At estimated T_inner :',1pg10.3,' K'/
     & t5,' Rsou.',t20,' Remp.',t35,' Rmax',t50,' <DILU>'/
     & t5, 1pg10.3,t20,1pg10.3,t35,1pg10.3,t50,1pg10.3/
     & t5,' <Hdens>',t20,' <Ndens>',t35,' Fill Factor',/
     & t5, 1pg10.3,t20,1pg10.3,t35,0pf9.6,//
     & t5,' FQ tot',t20,' Log FQ',t35,' LQ tot',t50,' Log LQ',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' Q(N) in',t20,' Log Q(N)',t35,' <Q(N)>',t50,' Log<Q(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' Q(H) in',t20,' Log Q(H)',t35,' <Q(H)>',t50,' Log<Q(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' U(N) in',t20,' Log U(N)',t35,' <U(N)>',t50,' Log<U(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     & t5,' U(H) in',t20,' Log U(H)',t35,' <U(H)>',t50,' Log<U(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3)
c
      write (luop,90) tinner,rstar,remp,rmax,wdpl,dhn,dht,fin,qht,
     &dlog10(qht),astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),
     &qhdnav,dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),
     &unin,dlog10(unin),unav,dlog10(unav),uhin,dlog10(uhin),uhav,
     &dlog10(uhav)
c
c
  100 format(t6,'MOD',t10,'Temp.',t22,'Alpha',t31,'Turn-on',
     & t39,'Cut-off'/
     &,t6,a2,x,1pg10.3,x,0pf8.3,x,0pf8.3,x,0pf8.3,//
     &,t6,'Zstar',t16,'FQHI',t27,'FQHEI',t38,'FQHEII'/
     &,t6,0pf8.4,1x,3(1pg10.3,x),/)
      write (luop,100) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (turbheatmode.eq.1) then
  110 format(/,'  Micro-Turbulent Dissipation Enabled, Mach = ',1pg10.3)
        write (*,110) admach
      endif
c
      if (grainmode.eq.1) then
        write (luop,140) galpha,amin(1),amax(1),amin(2),amax(2),
     &   graindens(1),graindens(2),bgrain,yinf
c
  120  format(  t5,' Projected dust area/H atom (graphite,silicate):',/
     &         ,t5,'  ',1pg12.4,1x,' (cm^2)',1x,1pg12.4,1x,' (cm^2)' )
        write (luop,120) siggrain(1),siggrain(2)
c
  130  format(  t5,' Composition Mass Ratios:  '/
     &         ,t5,'    Gas/H     :',1pg12.4/
     &         ,t5,'    Dust/Gas  :',1pg12.4/
     &         ,t5,'    PAH/Gas   :',1pg12.4/)
c
        write (luop,130) grainghr,graindgr,grainpgr
c
      endif
c
  140 format(   ' Dust Parameters :',//
     &         ,t5,' Alpha :',1pg10.3,/
     &         ,t5,' Graphite min:       max:',/
     &         ,t5,3x,1pg10.3,x,1pg10.3,/
     &         ,t5,' Silicate min:       max:',/
     &         ,t5,3x,1pg10.3,x,1pg10.3,/
     &         ,t5,' Graph dens  Sil dens',/
     &         ,t5,3x,1pg10.3,x,1pg10.3,/
     &         ,t5,' Bgrain      Yinf ',/
     &         ,t5,3x,1pg10.3,x,1pg10.3)
c
  150 format(   /,' Model Parameters:',//
     & ,t6,'Jden',t12,'Jgeo',t18,'Jend',t23,'Ielen',t29,
     & 'Jpoen',t38,'Fren',/
     & ,t6,3(a4,2x),2(i2,4x),0pf6.4,//
     & ,t6,'Tend',t12,'DIend',t22,'TAUen',t31,'Jeq',t38,'Teini',/
     & ,t5,0pf6.0,2(1pg10.3),x,a4,x,0pf7.1,//,
     & '::::::::::::::::::::::::::::::::::::',
     & '::::::::::::::::::::::::::::::::::::',/)
      write (luop,150) jden,jgeo,jend,ielen,jpoen,fren,tend,diend,tauen,
     &jeq,tm00
c
  160 format(/' Description of each space step:')
      write (luop,160)
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
        if (jgeo.eq.'S') then
  170 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     & t40,'dTau',t50,'<R>',t64,'Delta R',t78,'dr',t93,'nH',
     & t107,'ne',t120,'nt',t134,'F_HI',t148,'F_HII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)',t232,'Hab.P.',t246,'QPAH/ISRF')
          write (luop,170)
        else
  180 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     & t40,'dTau',t50,'<X>',t64,'Delta X',t78,'dx',t93,'nH',
     & t107,'ne',t120,'nt',t134,'F_HI',t148,'F_HII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)',t232,'Hab.P.',t246,'QPAH/ISRF')
          write (luop,180)
        endif
      else
        if (jgeo.eq.'S') then
  190 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     & t40,'dTau',t50,'<R>',t64,'Delta R',t78,'dr',t93,'nH',
     & t107,'ne',t120,'nt',t134,'F_HI',t148,'F_HII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)')
          write (luop,190)
        else
  200 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     & t40,'dTau',t50,'<X>',t64,'Delta X',t78,'dx',t93,'nH',
     & t107,'ne',t120,'nt',t134,'F_HI',t148,'F_HII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)')
          write (luop,200)
        endif
      endif
c
      if (jeq.eq.'F') write (luop,210)
  210 format(t73,'Fthick',t82,'Treach',t91,'DTCO',t100,'Dton',t109,
     &'DToff')
      if (jeq.eq.'P') write (luop,220)
  220 format(t100,'RECscal',t109,'DToff')
      close (luop)
      write (*,230) fnam
  230 format(/' Output in file : ',a/)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phsem file header
c
c     contains structure summary and emissivity of strong lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jgeo.eq.'S') then
  240  format(/,t3,'#,',t8,'[2]Delta R,',t22,'[3]dr,',t36,'[4]<R>,',
     &  t50,'[5]<Te>,',t64,'[6]nH,',t78,'[7]ne,',t92,'[8]nt,',
     & t106,'[9]F_HI,',t120,'[10]F_HII,',
     & t134,'[11]HBeta,',t148,'[12]5007,'
     & t162,'[13]4363,',t176,'[14]3727+,',t190,'[15]6300,',
     & t204,'[16]6584,',t218,'[17]6725+,')
        write (lunt,240)
      else
  250  format(/,t3,'#,',t8,'[2]Delta X,',t22,'[3]dx,',t36,'[4]<X>,',
     &  t50,'[5]<Te>,',t64,'[6]nH,',t78,'[7]ne,',t92,'[8]nt,',
     & t106,'[9]F_HI,',t120,'[10]F_HII,',
     & t134,'[11]HBeta,',t148,'[12]5007,'
     & t162,'[13]4363,',t176,'[14]3727+,',t190,'[15]6300,',
     & t204,'[16]6584,',t218,'[17]6725+,')
        write (lunt,250)
      endif
c
      close (lunt)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phapn file header
c
c     new format: contains plasma properties in line, for
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  260 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,
     &'  en',t70,'  FHI',t84,'  Density',t98,'  Pressure',t112,
     &'  Flow',t126,'  Ram Press.',t140,'  Sound Spd',t154,
     &'  Slab depth',t168,'  Alfen Spd.',t182,
     &'  Mag. Field',t196,'  Grain Pot.')
      write (lupf,260)
c
      close (lupf)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phlss file header
c
c     contains plasma cooling/heating in line, for
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  270 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,'  en',t70,
     &     '  Dloss',t84,'  Eloss',t98,'  Egain',t112,'  Coll. H',t126,
     &     '  Chrg. Ex.',t140,'  Reson',t154,
     &     '  Xreson',t168,'  Inter/fine',t182,'  He int-forb',t196,
     &     '  Forbid',t210,'  Fe II',t224,'  2Photon',t238,
     &     '  Compton',t252,'  Free-Free',t266,'  Coll. Ion',t280,
     &     '  Photo.',t294,'  Recomb.',t308,'  Cosmic',t322,
     &     '  GrainsH-C',t336,'  G.Heat',t350,'  G.Cool',t364,
     &     '  PAHs')
      write (lups,270)
c
      close (lups)
c
      if (jall.eq.'YES') close (lusl)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     monitor ions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jiel.eq.'YES') then
c
c     write ion balance files if requested
c
        do i=1,ieln
c
          if (i.eq.1) ie=iel1
          if (i.eq.2) ie=iel2
          if (i.eq.3) ie=iel3
          if (i.eq.4) ie=iel4
          fn=' '
          pfx=elem(ie)
          sfx='csv'
          call newfile (pfx, elem_len(ie), sfx, 3, fn)
          filn(i)=fn
        enddo
c
        do i=1,ieln
          open (luions(i),file=filn(i),status='NEW')
        enddo
c
c
        do i=1,ieln
c
          write (luions(i),'(" Run     : ",a96)') runname
c
  280  format(
     & ' [ 1] <R>    , [ 2] DeltaR , [ 3] dR     , [ 4] <T>    , ',
     & 31(' [',i2,'] ', a6,', '))
  290  format(
     & ' [ 1] <X>    , [ 2] DeltaX , [ 3] dX     , [ 4] <T>    , ',
     & 31(' [',i2,'] ', a6,', '))
          if (i.eq.1) then
            write (luions(i),'(" Element : ",a2)') elem(iel1)
            if (jgeo.eq.'S') then
              write (luions(i),280) (j+4,rom(j),j=1,maxion(iel1))
            else
              write (luions(i),290) (j+4,rom(j),j=1,maxion(iel1))
            endif
          endif
          if (i.eq.2) then
            write (luions(i),'(" Element: ",a2)') elem(iel2)
            if (jgeo.eq.'S') then
              write (luions(i),280) (j+4,rom(j),j=1,maxion(iel2))
            else
              write (luions(i),290) (j+4,rom(j),j=1,maxion(iel2))
            endif
          endif
          if (i.eq.3) then
            write (luions(i),'(" Element: ",a2)') elem(iel3)
            if (jgeo.eq.'S') then
              write (luions(i),280) (j+4,rom(j),j=1,maxion(iel3))
            else
              write (luions(i),290) (j+4,rom(j),j=1,maxion(iel3))
            endif
          endif
          if (i.eq.4) then
            write (luions(i),'(" Element: ",a2)') elem(iel4)
            if (jgeo.eq.'S') then
              write (luions(i),280) (j+4,rom(j),j=1,maxion(iel4))
            else
              write (luions(i),290) (j+4,rom(j),j=1,maxion(iel4))
            endif
          endif
c
        enddo
c
        do i=1,ieln
          close (luions(i))
        enddo
c
      endif
c
      if (jiem.eq.'YES') then
c
c     write ion emission files if requested
c
        do i=1,ieln
c
          idx=iemidx(i)
c
          at=fmatom(idx)
          io=fmion(idx)
          nt=nfmtrans(idx)
          nl=fmnl(idx)
          fn=' '
          pfx(1:16)='_'
          pfx(1:elem_len(at))=elem(at)
          plen=elem_len(at)+1
          pfx(plen:plen)='_'
          plen=plen+1
          pfx(plen:plen+rom_len(io))=rom(io)
          plen=plen+rom_len(io)+1
          pfx(plen:plen)='_'
          sfx='csv'
          call newfile (pfx, plen, sfx, 3, fn)
          filnem(i)=fn
c
          open (luemiss(i),file=filnem(i),status='NEW')
c
          write (luemiss(i),'("Run: ",a96)') runname
          write (luemiss(i),*) 'Emissivity (erg/cm^3/s/sr) for Ion: ',
     &     elem(at),rom(io)
  300     format(
     & ' Step [1], <X> [2], DeltaX [3], dX [4], <T> [5], <HB> [6]',
     & 40(',',1pg12.6,'[',i2,']'))
          write (luemiss(i),300) (1.0d8*fmlam(itr,idx),itr+6,itr=1,nt)
c
          close (luemiss(i))
        enddo
c
      endif
c
      if (jlin.eq.'YES') then
 440    format(/,'     ,             ,             ,',
     &   '             ,             ,             ',
     &   16(',    ',a2,a6,'    '))
       write (lulin,440)(elem(emlinlistatom(itr)),rom(emlinlistion(itr))
     &                 ,itr=1,njlines)
 450   format( ' #[1],   Delta X[2],        dx[3],',
     & '       <X>[4],       <T>[5],     HBeta[6]',
     &   16(',',f12.3,'[',i2,']'))
        write(lulin,450) (emlinlist(itr),itr+6,itr=1,njlines)
        close (lulin)
      endif
c
      return
c
      end
