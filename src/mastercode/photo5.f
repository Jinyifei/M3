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
c
c     CHOICE OF EQUILIBRIUM OR FINITE AGE CONDITIONS
c
c     CALL SUBR.: PHOTSOU,COMPPH5
c
c     NB. COMPUTATIONS PERFORMED IN SUBROUTINE COMPPH5
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine photo5 ()
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      real*8 epotmi,fin, uinit, qinit
      real*8 qhlo,rcin,rechy,reclo
      real*8 rstsun,vstrlo,wid, starlum
      real*8 presreg,dh10, scale
c
      integer*4 j,luin,i,kmax, k
      integer*4 idx1 ,idx2, idx3, idx4, iem
c
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
c
      write (*,10)
c
c     ***INITIAL IONISATION CONDITIONS
c
   10 format(///
     & ' ****************************************************'/
     & '  Photoionisation model P5 selected:'//
     & '  Diffuse Field : Outward integration: Radiation Pressure',/,
     & ' ****************************************************')
c
      epotmi=iphe
c
c     set up ionisation state
c
      model='protoionisation'
      call popcha (model)
      model='Photo 5'
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
      model='Photo 5'
c
      call photsou (model)
c
c     Possible to have no photons
c
      qhlo=dlog(qht+epsilon)
      rechy=2.6d-13
      reclo=dlog(rechy)
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
  130    format(/' Define the source size or luminosity '/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &        '    R   :  By Radius'/
     &        '    L   :  By Luminosity'/
     &        '    P   :  By Ionizing Photons'//
     &        ' :: ',$)
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
     &            ,' (in solar units (<1.e8) or in cm (>1.e8) : ',$)
            write (*,140)
            read (*,*) rstsun
            if (rstsun.le.0.0d0) goto 120
            if (rstsun.lt.1.d8) then
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
     &     '  Source Total Luminosity              : ',1pg12.5/
     &     '  Source Ionising (13.6eV+) Luminosity : ',1pg12.5/
     &     '  Source Ionising (13.6eV+) Photons    : ',1pg12.5/
     &     ' ****************************************************')
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
            rstsun=rstar/6.96d10
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
c     P/s = qht*4.d0*pi*rstsun*rstsun
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
            rstsun=rstar/6.96d10
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
  280    format(/' Give pressure regime (p/k, <10 as log) : ',$)
        write (*,280)
        read (*,*) presreg
        if (presreg.le.10.d0) presreg=10.d0**presreg
c
  290    format(//'  Give estimate of mean temperature ',
     &        ' [T ~1e4, <10 as log]'/
     &        ' (for initial size estimates only, actual temperature'/
     &        ' and structure calculated later in detail, '/
     &        ' if in doubt use 1e4 (or 4))'/
     &        ' : ',$)
        write (*,290)
        read (*,*) tinner
        if (tinner.le.10.d0) tinner=10.d0**tinner
c
        dhn=presreg/tinner
        dh10=presreg/1e4
c
  300  format(/' ****************************************************'/
     &         '  H number density at 1e4 K : ',1pg12.5/
     &         ' ****************************************************')
        write (*,300) dh10
c
c     jden = B isobaric
c
      endif
c
  310 if (jden.eq.'C') then
  320    format(//' Give nominal hydrogen density : ',$)
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
      if ((jgeo.eq.'S')) then
c
c       Spherical Geometry
c
        if (qht.gt.0.d0) then
c
          vstrlo=((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(dht)))
     &     -dlog(fin))-reclo
          rmax=dexp((vstrlo-1.43241d0)/3.d0)
c
c     write (*,*) qht,dht,dhn
c
          qhdnin=qht/dht
          qhdnav=(4.d0*qht/dht)*fdilu(rstar,(rmax+rstar)*0.5d0)
c
          uhin=qht/(cls*dhn)
          uhav=(dht/dhn)*qhdnav/cls
c
          qhdhin=qht/dhn
          qhdhav=qhdnav*(dht/dhn)
c
  350     write (*,360) rmax,qht*astar,qht,qhdnin,qhdhin,uhin,qhdnav,
     &     qhdhav,uhav
  360 format(//,
     & ' ****************************************************',/,
     & '   Filled Sphere Parameters:',/,
     & ' ****************************************************',/,
     & '  Estimated HII   Stromgren radius:',1pg10.3,' cm.',/,
     & ' ****************************************************',/,
     & '  Photon Luminosity LQH    :',1pg10.3, ' Phot/s',/,
     & '  Photon Flux       FQH    :',1pg10.3, ' Phot/cm2/s',/,
     & ' ****************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s'/,
     & '        QHDH inner   : ',1pg12.5,' cm/s'/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <QHDN>     : ',1pg12.5,' cm/s'/,
     & '          <QHDH>     : ',1pg12.5,' cm/s'/,
     & '          <U(H)>     : ',1pg12.5,/,
     & ' ****************************************************',//
     & ' Give initial radius in terms of distance or Q(N),',
     & ' Q(H), or U(H) (d/q/h/u):',$)
          read (*,20) ilgg
          ilgg=ilgg(1:1)
c
          if (ilgg.eq.'D') ilgg='d'
          if (ilgg.eq.'Q') ilgg='q'
          if (ilgg.eq.'H') ilgg='h'
          if (ilgg.eq.'U') ilgg='u'
          if ((ilgg.ne.'q').and.(ilgg.ne.'h').and.(ilgg.ne.'d')
     &     .and.(ilgg.ne.'u')) goto 350
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
            remp=((rmax+rstar)/2.d0)*dsqrt(uhav/uinit)
c
          endif
c
          if (ilgg.eq.'q') then
            write (*,380)
  380       format(//,' Give QHDN at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)/2.d0)*dsqrt(qhdnav/qinit)
c
          endif
c
          if (ilgg.eq.'h') then
            write (*,390)
  390       format(//,' Give QHDH at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)/2.d0)*dsqrt(qhdhav/qinit)
c
          endif
c
          if (ilgg.eq.'d') then
c
c     give R_0 directly
c
            write (*,400)
  400       format(//' Give the initial radius ',
     &         ' (in cm (>1) or in fraction of Stromgren radius) : ',$)
            read (*,*) remp
c
            if (remp.lt.0.0d0) goto 350
            if (remp.le.1) remp=remp*rmax
            if (remp.lt.rstar) remp=rstar
c
          endif
c
          rmax=dexp(((vstrlo+dlog(1.d0+((remp/rmax)**3)))-1.43241d0)/
     &     3.0d0)
c
c
c     end with photons (qht>0)
c
        else
c
c     no photons
c
c
c     give R_0 directly
c
  410     write (*,420)
  420        format(//' Give the initial radius:',$)
          read (*,*) remp
c
          if (remp.lt.0.0d0) goto 410
          if (remp.lt.rstar) remp=rstar
          rmax=remp*100.0d0
c
        endif
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
c     Plane Parallel (Sr=1/pi) x pi
c
        blum=pi*blum
        ilum=pi*ilum
        wdpl=fdilu(rstar,remp)
        qhdnin=4.d0*qht/dht*wdpl
        qhdhin=qhdnin*(dht/dhn)
        uhin=qhdhin/cls
        blum=blum*wdpl
        ilum=ilum*wdpl
        wdpl=fdilu(rstar,(rmax+remp)*0.5d0)
        qhdnav=4.d0*qht/dht*wdpl
        qhdhav=qhdnav*(dht/dhn)
        uhav=qhdhav/cls
  430    format(//,
     & ' ****************************************************',/,
     & '   Partially Filled Sphere Parameters:',/,
     & ' ****************************************************',/,
     & '  Empty inner radius :',1pg12.5,' cm'/,
     & '  Outer HII   radius :',1pg12.5,' cm.',/,
     & ' ****************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s'/,
     & '        QHDH inner   : ',1pg12.5,' cm/s'/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <QHDN>     : ',1pg12.5,' cm/s'/,
     & '          <QHDH>     : ',1pg12.5,' cm/s'/,
     & '          <U(H)>     : ',1pg12.5,/,
     & '   Total intensity   : ',1pg12.5,' erg/s/cm2' /,
     & '   Ionizing intensity: ',1pg12.5,' erg/s/cm2'/,
     & ' ****************************************************',/)
        write (*,430) remp,rmax,qhdnin,qhdhin,uhin,qhdnav,qhdhav,uhav,
     &   blum,ilum
c
        radin=0.d0
        radout=1.d38
c
        write (*,440)
  440    format(//' Volume integration over the whole sphere?',
     &        ' (y/n) : ',$)
        read (*,90) ilgg
c
        if ((ilgg(1:1).eq.'Y').or.(ilgg(1:1).eq.'y')) ilgg='Y'
        if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg='N'
c
        if (ilgg.ne.'N') ilgg='Y'
c
        if (ilgg.eq.'N') then
  450     write (*,460)
  460       format(/' Integration through the line of sight ',
     &           'for a centered ring aperture.'/
     &           'Give inner and outer radii (cm) :')
          read (*,*) radin,radout
          if ((radin.lt.0.0d0).or.(radin.ge.(0.999d0*radout))) goto 450
        endif
c
c   end Spherical
c
      endif
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
c
  470     write (*,480)
  480    format(/,' Give Ionizing Flux at inner edge by: ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    B  : Total Flux (erg/s/cm^2)',/,
     &     '    I  : Ionising Flux (erg/s/cm^2)',/,
     &     '    F  : Ionising Photons FQHI (phot/s/cm^2)',/,
     &     '    U  : Ionisation parameter U(H) ',/,
     &     '    Q  : Ionisation parameter QHDN ',/,
     &     '    H  : Ionisation parameter QHDH ',/,
     &     '    N  : No change (Use current flux)',/,
     &     ' :: ',$)
          read (*,90) ilgg
c
          if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg='B'
          if ((ilgg(1:1).eq.'I').or.(ilgg(1:1).eq.'i')) ilgg='I'
          if ((ilgg(1:1).eq.'F').or.(ilgg(1:1).eq.'f')) ilgg='F'
          if ((ilgg(1:1).eq.'U').or.(ilgg(1:1).eq.'u')) ilgg='U'
          if ((ilgg(1:1).eq.'Q').or.(ilgg(1:1).eq.'q')) ilgg='Q'
          if ((ilgg(1:1).eq.'H').or.(ilgg(1:1).eq.'h')) ilgg='H'
          if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg='N'
          if (ilgg.eq.'B') then
            write (*,490)
  490       format(/,'Give Bolometric flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/blum
          else if (ilgg.eq.'I') then
            write (*,500)
  500       format(/,'Give Ionising flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/ilum
          else if (ilgg.eq.'F') then
            write (*,510)
  510       format(/,
     & 'Give Ionising Photon flux at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100) scale=10**scale
            scale=scale/qht
          else if (ilgg.eq.'U') then
            write (*,520)
  520     format(/,'Give U(H) at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0) scale=10**scale
            scale=scale*dhn*cls/qht
          else if (ilgg.eq.'Q') then
            write (*,530)
  530       format(/,'Give QHDN at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0) scale=10**scale
            scale=scale*dht/qht
          else if (ilgg.eq.'H') then
            write (*,540)
  540       format(/,'Give QHDH at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0) scale=10**scale
            scale=scale*dhn/qht
          else if (ilgg.eq.'N') then
            scale=1.d0
          else
            goto 470
          endif
c
          do j=1,infph
            soupho(j)=soupho(j)*scale
          enddo
c
c
          blum=blum*scale
          ilum=ilum*scale
          qht=qht*scale
          qhdnav=qht/dht
          qhdhav=qht/dhn
          uhav=qht/(cls*dhn)
  550 format(//,
     & ' ****************************************************',/,
     & '  Total Flux at inner      : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Flux at inner   : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Phot Flux FQ, Ryd+  : ',1pg12.5,' (phots/s/cm^2)',/,
     & '  Ionisation parameter U(H) inner : ',1pg12.5/
     & '  Ionisation parameter QHDN inner : ',1pg12.5/
     & '  Ionisation parameter QHDH inner : ',1pg12.5/
     & ' ****************************************************',/)
          write (*,550) blum,ilum,qht,uhav,qhdnav,qhdhav
c
          vstrlo=((qhlo-(2.0d0*dlog(dht)))-dlog(fin))-reclo
          rmax=dexp(vstrlo)
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
          qhdhin=((2.0d0*qht)*wdpl)/dhn
          qhdhav=qhdhin
c
          uhin=((2.0d0*qht)*wdpl)/(cls*dhn)
          uhav=uhin
c
c
        endif
c
      endif
c
      if (jgeo.eq.'F') then
c
        if (qht.gt.0.d0) then
c
          vstrlo=((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(dht)))
     &     -dlog(fin))-reclo
          rmax=dexp((vstrlo-1.43241d0)/3.d0)
          qhdnin=qht/dht
          qhdnav=((2.d0*qht)/dht)*fdilu(rstar,(rmax+rstar)*0.5d0)
c
  580     write (*,590) rmax,qhdnin,qhdnav
  590   format(/,
     & ' ****************************************************'/
     & '  Estimated hydrogen Stromgren length of the'/
     & '  order of ................:',1pg10.3,' cm.'/
     & '  Impact parameters:QHDN in:',1pg10.3,/
     & '                    QHDN av:',1pg10.3,/
     & ' ****************************************************'//
     & ' Give the radius of empty zone in front of source '/
     & ' (in cm (0,>100)) or in log(cm) (>0<=100)) : ',$)
          read (*,*) remp
c
          if (remp.lt.0.0d0) goto 580
          if ((remp.le.1.d2).and.(remp.gt.0.d0)) remp=10.d0**remp
c
          rmax=dexp(((vstrlo+dlog(1.d0+((remp/rmax)**3)))-1.43241d0)/
     &     3.0d0)
c
        else
c
c     give R_0 directly
c
  600     write (*,610)
  610        format(/' Give the initial radius:',$)
          read (*,*) remp
c
          if (remp.lt.0.0d0) goto 600
          if (remp.lt.rstar) remp=rstar
          rmax=remp*100.0d0
c
        endif
c
        write (*,*) rstar,remp,rmax
        wdpl=fdilu(rstar,remp)
        qhdnin=((4.0d0*qht)*wdpl)/dht
        wdpl=fdilu(rstar,(rmax+remp)/2.0d0)
        qhdnav=((4.0d0*qht)*wdpl)/dht
c
      endif
c
c
  620 admach=0.0d0
      turbheatmode=0
c 440  format(//,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     & '  Micro-Turbulent Dissipation Heating',/,
c     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
c     & '     Choose Timescale Type:',/,
c     & '     A  : Recombination Timescale.',/,
c     & '     B  : Collisional Timescale.',/,
c     & '     C  : Cooling Timescale.',/,
c     & '     D  : Fixed Timescale.',/,
c     & '     N  : No Dissipaton.',/,
c     & '   ::', $)
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
c     &          '  Micro-Turbulent Dissipation Enabled',/,
c     &          ' ::::::::::::::::::::::::::::::::::::::::::::::::::',/)
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
  630 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Ionisation Balance Calculations'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
      write (*,630)
  640 write (*,650)
  650 format(//'  Choose type'/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '     E  : Equilibrium balance.'/
     & '     F  : Finite source lifetime.'/
     & '     P  : Post-Equilibrium decay.'/
     & '  :: ',$)
      read (*,90) jeq
c
      if ((jeq(1:1).eq.'E').or.(jeq(1:1).eq.'e')) jeq='E'
      if ((jeq(1:1).eq.'F').or.(jeq(1:1).eq.'f')) jeq='F'
      if ((jeq(1:1).eq.'P').or.(jeq(1:1).eq.'p')) jeq='P'
      if ((jeq(1:1).eq.'c').or.(jeq(1:1).eq.'c')) jeq='C'
c
      if ((jeq.ne.'E').and.(jeq.ne.'F').and.(jeq.ne.'P')) goto 640
c
c     forced CIE for given tprofile.
c     useful for hot x-ray bubbles.
c
      if (jeq.eq.'C') then
c
        write (*,660)
  660    format(/'  Fixed Thermal Profile:'/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//
     &        ' f(x) = t0*[fx*exp((x-ta)/r0) + fp*((x/ta)**tb)]+tc  '/
     &        ' (with r in r0 units , t0<10K as log)'/
     &       ' Give parameters t0 ta tb tc r0 : ',$)
        read (*,*) tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
      endif
c
c     non equilibrium
c
      if ((jeq.ne.'E').and.(jeq.ne.'C')) then
        rcin=0.9
        if (jgeo.eq.'S') then
          dti=-(dlog(1.0d0-(rcin**3))/(rechy*dht))
        else
          dti=-(dlog(1.0d0-rcin)/(rechy*dht))
        endif
        dti=(rmax/3.d10)+dmax1(dti,rmax/3.d10)
        dthre=(1.2/dht)/3d-13
c
c     finite lifetime
c
        if (jeq.eq.'F') then
  670     write (*,680) rmax,dti
  680       format(/'  Source turn on and finite lifetime:'/
     &       ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' (Hydrogen Stromgrem radius of the order of',1pg10.3,
     &       ' cm'/' Characteristic time scale for ionisation: ',1pg7.1,
     &       ' sec)'/' Give elapsed time since source turned on (sec)'/
     &     ' and give life-time of ionising source (sec) : ',$)
          read (*,*) telap,tlife
          if ((telap.le.0.0d0).or.(tlife.le.0.0d0)) goto 670
          tm00=100.0d0
  690     write (*,700)
  700      format(//' Initial temperature of neutral gas : ',$)
          read (*,*,err=690) tm00
          if (tm00.lt.1.0d0) goto 690
        else
          tlife=0.0d0
  710     write (*,720) rmax,dthre
  720       format(/'  Source turn-off and final luminosity'/
     &       ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &      ' (Hydrogen Stromgrem radius of the order of',1pg10.3,' cm'
     &       /'Characteristic time scale for recombination: ',1pg7.1,
     &       ' sec)'/' Give elapsed time since source turned off (sec)'/
     &      ' and fractional final lumnosity of the source (>=0): ',$)
          read (*,*) telap,frlum
          if ((telap.le.0.0d0).or.(frlum.lt.0.0d0)) goto 710
        endif
      endif
c
  730 write (*,740)
  740 format(//' Give the step photon absorption fraction: ',$)
      read (*,*) dtau0
      if ((dtau0.le.0.0d0)) goto 730
c
c     standard convergence limits
c
      difma=0.05d0
      dtlma=0.03d0
      dhlma=0.050d0
c
      if (expertmode.gt.0) then
  750    format(/' Do you wish to alter the convergence criteria ?',
     &        ' (not recommended) : ',$)
        write (*,750)
        read (*,90) ilgg
        ilgg=ilgg(1:1)
        if (ilgg.eq.'n') ilgg='N'
        if (ilgg.eq.'y') ilgg='Y'
        if (ilgg.ne.'Y') ilgg='N'
c
        if (ilgg.eq.'Y') then
          write (*,760)
  760      format(/,' Set the convergence criteria, the maximum changes'
     &         ,'allowed',/
     &         ,' for three parameters: ion population, te and hydrogen'
     &         ,'density.',/
     &         ,' Convergence criteria: difma (0.05),dtlma (0.03)'
     &         ,'dhlma(0.015) : ',$)
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
        write (*,770) dtau0
  770    format(//' Model Summary :',/
     &          ,'    Mode  : Thermal and Ionic Equilibrium ',/
     &          ,'    dTau  :',0pf8.5)
      else if (jeq.eq.'F') then
        write (*,780) telap,tlife,dtau0
  780  format(//' Model Summary :',/
     & ,'    Mode        :  Non-Equilibrium + Finite Source Life',/
     & ,'    Age         :',1pg10.3,' sec',/
     & ,'    Source Life :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      else
        write (*,790) telap,dtau0
  790 format(//' Model Summary :',/
     & ,'    Mode        :  Equilibrium + Source Switch off',/
     & ,'    Switch off  :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      endif
      if (usekappa) then
  800    format(/' Kappa Electron Distribution Enabled :'/
     &           '    Electron Kappa : ',1pg11.4)
        write (*,800) kappa
      endif
c
  810 format(/' Radiation Field :',//,
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
     & t5,' U(H) in',t20,' Log U',t35,' <U(H)>',t50,' Log<U>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3)
c
      write (*,810) rstar,remp,rmax,wdpl,dhn,dht,fin,qht,dlog10(qht),
     &astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),qhdnav,
     &dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),uhin,
     &dlog10(uhin),uhav,dlog10(uhav)
c
c
c     ***TYPE OF EXIT FROM THE PROGRAM
c
  820 ielen=1
      jpoen=1
      tend=10.0d0
      fren=0.01d0
      diend=0.0d0
      tauen=0.0d0
c
c
  830 format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Boundry Conditions  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,830)
c
  840 format(/'  Choose a model ending : '/
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
      write (*,840) fren
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
     &.and.(jend.ne.'R').and.(jend.ne.'G')) goto 820
      if (jend.eq.'R') return
      if (jend.eq.'G') goto 620
      if ((jend.eq.'B').or.(jend.eq.'D')) then
  850   write (*,860)
  860    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 850
        ielen=zmap(ielen)
        jpoen=1
        if (arad(2,ielen).le.0.0d0) jpoen=2
      endif
c
      if (jend.eq.'B') then
  870   write (*,880) elem(int(ielen))
  880    format(/' Give the final ionisation fraction of ',a2,' : ',$)
        read (*,*) fren
        if ((fren.lt.0.0d0).or.(fren.gt.1.0d0)) goto 870
      endif
c
      if (jend.eq.'D') then
  890   write (*,900) elem(int(ielen))
  900 format(/' Give the final optical depth at threshold of ',a2,':',$)
        read (*,*) tauen
        if (tauen.le.0.0d0) goto 890
      endif
c
      if (jend.eq.'C') then
  910   write (*,920)
  920    format(/' Give the final temperature (<10 as log): ',$)
        read (*,*) tend
        if (tend.lt.1.0d0) goto 910
      endif
c
      if (jend.eq.'E') then
  930   write (*,940)
  940    format(/' Give the distance or radius at which the density',
     &        ' drops: '/
     &        ' (in cm (>1E6) or as a fraction of the Stromgren'/
     &        ' radius (<1E6)) : ',$)
        read (*,*) diend
        if (diend.lt.1.d6) diend=diend*rmax
        diend=remp+diend
        if (diend.le.remp) goto 930
      endif
c
      if (jend.eq.'F') then
        write (*,950)
  950    format(/' Give the final column density (<100 as log): ',$)
        read (*,*) colend
        if (colend.lt.1.0d2) colend=10.d0**colend
  960   write (*,970)
  970    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 960
        ielen=zmap(ielen)
  980   write (*,990)
  990    format(/' Applies to ion stage : ',$)
        read (*,*) jpoen
        if (jpoen.le.0) goto 980
        if (jpoen.gt.maxion(ielen)) jpoen=maxion(ielen)+1
      endif
c
 1000 format(//
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Output Requirements  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      write (*,1000)
c
c
 1010 format(//'  Choose output settings : '/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '    A  :   Standard output (photnxxxx,phapnxxxx).'/
     & '    B  :   Standard + monitor 4 element ionisation.'/
     & '    C  :   Standard + all ions file.'/
     & '    D  :   Standard + final source + nu-Fnu spectrum.'/
     & '    F  :   Standard + first balance.'/
     & '    G  :   Everything'//
     & '    I  :   Standard + B + monitor 4 multi-level ion emission.'/
     & ' :: ',$)
      write (*,1010)
      read (*,90) ilgg
c
      if ((ilgg(1:1).eq.'A').or.(ilgg(1:1).eq.'a')) ilgg='A'
      if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg='B'
      if ((ilgg(1:1).eq.'C').or.(ilgg(1:1).eq.'c')) ilgg='C'
      if ((ilgg(1:1).eq.'D').or.(ilgg(1:1).eq.'d')) ilgg='D'
      if ((ilgg(1:1).eq.'F').or.(ilgg(1:1).eq.'f')) ilgg='F'
      if ((ilgg(1:1).eq.'G').or.(ilgg(1:1).eq.'g')) ilgg='G'
      if ((ilgg(1:1).eq.'I').or.(ilgg(1:1).eq.'i')) ilgg='I'
c
      jiel='NO'
      jiem='NO'
      if (ilgg.eq.'B') jiel='YES'
      if (ilgg.eq.'G') jiel='YES'
      if (ilgg.eq.'G') jiem='YES'
      if (ilgg.eq.'I') jiel='YES'
      if (ilgg.eq.'I') jiem='YES'
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
 1020    format(//' Give 4 elements to monitor (atomic numbers): ',$)
        write (*,1020)
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
 1030 format(//' Choose 4 multi-level species #',i2,' :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,1030) ,iem
        do i=1,(nfmions/5)+1
          j=(i-1)*5
          kmax=min(j+5,nfmions)-j
          write (*,'(5(2x,i2,": ",a2,a6))') (j+k,elem(fmatom(j+k)),
     &     rom(fmion(j+k)),k=1,kmax)
        enddo
 1040 format(' :: ',$)
        write (*,1040)
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
 1050    format (a8)
 1060    format(//' Give a prefix for final source file : ',$)
        write (*,1060)
        read (*,1050) jpfx
      endif
c
      jbal='NO'
      jbfx='p5bal'
c
      if ((ilgg.eq.'F').or.(ilgg.eq.'G')) then
c
        jbal='YES'
c
c     get final field file prefix
c
 1070    format(//' Give a prefix for first ion balance file : ',$)
        write (*,1070)
        read (*,1050) jbfx
      endif
c
c     get runname
c
 1080 format (a64)
 1090 format(//' Give a name/code for this run: ',$)
      write (*,1090)
      read (*,1080) runname
c
c     remains of old vax batch system (not used)
c
      banfil='INTERACTIVE'
c
      if (jden.eq.'B') dhn=dh10
c
      call compph5 (dhn, fin, banfil, difma, dtlma, dhlma)
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c
c       PHOTO5:  As for P4 but adds integrated radiation pressure
c               (altered for dust pressure)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compph5 (dhn, fin, banfil, difma, dtlma, dhlma)
c
      include 'cblocks.inc'
c
c
      real*8 poppre(mxion, mxelem), popend(mxion, mxelem)
      real*8 ppre(mxion, mxelem),difma,dtlma,dhlma
      real*8 thist(3),dpw,dxp
      real*8 popz(mxion, mxelem),difp,treap
      real*8 ctime,rtime,ptime,eqtime,stime
c
      real*8 aadn,agmax,agmu,aper,aperi,apero
      real*8 a,b,c,x,fx,fp,n0,r0,fv,e,f
      real*8 de,dedhma,dedhmi,deeq
      real*8 dh,dhf,dhl,dhn,dhp,frp
      real*8 difg,difmi,difp0,dift
      real*8 disin,dison,dlos0,dlos1,dlosav
      real*8 dr,drp,drxp,drzf,drzi
      real*8 dtau,dtauon,dtaux,dtco,dtco0,dtco1
      real*8 dtl,dtoff,dton,dton0,dton1,durmax
      real*8 dv,dva,dvoluni,eqscal,exl,exla
      real*8 fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw,frx,fthick
      real*8 prescc,prf,protim,rad,radi,recscal,rmm,rww
      real*8 t,taux,te1a,tef,telast,temp1,temp2,tempre
      real*8 tep,tex,texm,tii0,tprop,trea,runit
      real*8 treach,treach0,treach1,tspr
      real*8 wd,wdil0,wdil1,wei,weiph
      real*8 radpint, pres0, pres1
      real*8 popinttot
      real*8 a_v,ffuv,dnu,eng,habing,pahtest(mxinfph)
      real*8 pahsum, pahq0, pahqion, pahqfac
      real*8 q1,q2,q3,q4
      real*8 def
c
      integer*4 np,luop,lups,lusp,lusl,lupf,ir1
      integer*4 lut0,m,maxio,mma,n,nidh,niter
      integer*4 i,j,k,atom
      integer*4 idx, nt, itr
      integer*4 luions(4),luemiss(4)
      integer*4 fuvmin,fuvmax,pahabsmax
c
c     Old computer time calculations
c      real tarray(2)
c      real dtarr,dtime
c      integer*4 nt0,nt1, time
c
      character jjd*4, pollfile*12,tab*4
      character newfil*32, filna*32, banfil*32, fnam*32, filnb*32
      character imod*4, lmod*4, mmod*4, nmod*4,wmod*4,ispo*4
      character linemod*4, spmod*4
      character fn*32
      character filnc*32, filnd*32, fspl*32
      character filn(4)*32, filnem(4)*32
      character pfx*16,caller*4,sfx*4
      character irfile*32,chargefile*32
      logical iexi
c
c     External Functions
c
      real*8 fcolltim,fcrit,fdilu,feldens,fradpress
      real*8 fphotim,fpressu,frectim,frectim2
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
c
      irfile=' '
      tab=char(9)
c
c
c      nt0 = time()
c      nt1 = 0
      mma=mxnsteps-1
c      scalen = 1.0d16
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
      lut0=0
      ir1=91
c
      ieln=4
      do i=1,ieln
        luions(i)=30+i
        luemiss(i)=40+i
      enddo
c
c
c    ***DERIVE FILENAME : PHN**.PH5
c
      call p5fheader (newfil, banfil, fnam, filna, filnb, filnc, filnd,
     &filn, filnem, luop, lupf, lusl, lusp, lups, luions, luemiss, dhn,
     &fin)
      chargefile=' '
c
c  Set up for PAH existence
c
      fuvmin=1
      fuvmax=1
      pahabsmax=1
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
        k=1
        do while (photev(k).le.6.0)
          k=k+1
        enddo
c     FUVMIN = photon bin corresponding to 6 eV energy
        fuvmin=k-1
        do while (photev(k).le.1.359907296e+1)
          k=k+1
        enddo
c     FUVMAX = photon bin corresponding to 13.6 eV energy
        fuvmax=k-1
        do while (photev(k).le.262.012)
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
          if (eng.lt.9.26) then
            ffuv=ffuv+2.548e-18*eng**(-1.3322)/ev*pahnabs(k)*dnu
          else if (eng.lt.11.2) then
            ffuv=ffuv+1.049e-16*eng**(-3.d0)/ev*pahnabs(k)*dnu
          else
            ffuv=ffuv+4.126e-13*eng**(-6.4172)/ev*pahnabs(k)*dnu
          endif
        enddo
        pahq0=ffuv/pahsum
        write (*,10) pahq0
   10 format(/'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/
     &        'ISRF photon flux :',1pg12.5,/
     &   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
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
      fi=fin
      fi0=fin
      fi1=fin
      fidhcc=fin*dhn
      dv=1.0d6
      runit=(rmax-remp)/300.d0
      dtau=dtau0/2.d0
      if (jgeo.eq.'S') then
        vunilog=3.0d0*dlog10(runit)
      else
        vunilog=dlog10(runit)
      endif
      if ((jeq.eq.'E').or.(jeq.eq.'P').or.(jeq.eq.'C')) then
c     Equilibrium structure (E=equilibrium, P=post-equilibrium decay,
c     C=collisional ionization equilibrium)
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
      ponk=0.d0
      if (jden.eq.'B') ponk=dhn*1e4
      radpint=0.d0
c
c***********************************************
c
c      write(*,*)'***ITERATE, M = STEP NUMBER'
c
c***********************************************
c
      m=1
c
c
   20 jspot='NO'
      niter=0
      nidh=0
      difmi=1.d30
      wei=fcrit(dtau/dtau0,1.0d0)
      dtau=(dtau0+((19.d0*wei)*dtau))/(1.0d0+(19.d0*wei))
      wei=3.d0
      agmu=(agmax+(wei*agmu))/(1.d0+wei)
c
c****************************************************************
c
c      write(*,*)'IF M=1 , CHECK CONVERGENCE FOR THE DENSITY'
c      write(*,*)'AND IONIC POP. AT THE INNER BOUNDARY'
c
c****************************************************************
      if (m.eq.1) then
c
        if (jeq.eq.'C') then
c     Equilibrium structure (C=collisional ionization equilibrium)
c     TE0 = Electron temperature at inner spatial step boundary
          te0=frad(remp,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen,
     &     0.d0,1.d0,0.d0)
        else
          te0=1.d4
        endif
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
c     DEDHMA = total metallicity of the gas (relative to hydrogen)
c     DEDHMI = total metallicity of gas for all species with
c              -ve or 0 recombination rates from 2nd ion state???
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
        if (jgeo.eq.'F') then
          rad0=0.d0
          dis0=remp
          wdil0=fdilu(rstar,dis0)
        endif
        if (jgeo.eq.'P') then
          rad0=0.d0
          dis0=remp
          wdil0=wdpl
        endif
c
c  Loop back here if hydrogen density change is too great
c
   30   if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
        lmod='DW'
c     Compute local emissivity in vectors EMIDIF (diffuse field) and
c     EMILIN (line field) used by routines TOTPHOT and NEWDIF
        call localem (te0, de0, dh0)
c     Compute mean intensity of radiation in vector TPHOT using
c     source vector SOUPHO
        call totphot (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
c
        if ((jeq.eq.'E').or.(jeq.eq.'P').and.(jthm.eq.'S')) then
          dton0=1.d33
          call teequi (te0, te0, de0, dh0, dton0, nmod)
        else if ((jeq.eq.'C').or.(jthm.eq.'T')) then
c     Equilibrium structure (C=collisional ionization equilibrium)
c     Thermal structure (T=isothermal)
c     Compute ionization equilibrium at initial temperature and
c     density for all elements. Outputs are electron density DE
c     and the fractional abundances of the different species of
c     each element POP(6,11)
          call equion (te0, de0, dh0)
c
          call localem (te0, de0, dh0)
          call totphot (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
          call equion (te0, de0, dh0)
c
          call localem (te0, de0, dh0)
          call totphot (te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
          call equion (te0, de0, dh0)
c
        else
          distim(1,0)=dis0
          distim(2,0)=dis0/cls
          treach0=distim(2,0)
          dton0=dmax1(0.d0,dmin1(tlife,telap-((2.d0*dis0)/cls)))
          jjd=jden
          jden='C'
          exl=-1.d0
          tex=0.95d0*tm00
          prescc=dhn*rkb*1.d4+fradpress(0.0d0,dhn)
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
        else if (jden.eq.'F') then
          dh0=frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)
        else
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
          frp=fradpress(dr,dhn)
          prescc=dhn*rkb*1.d4+radpint+frp
          ponk=(prescc/rkb)
          pres0=fpressu(te0,dhn,pop)
          dh0=dhn*prescc/(pres0)
        endif
c
        dhl=dif(dhp,dh0)
        write (*,80) te0,dhl
        if (dhl.ge.dhlma) goto 30
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
        qhdhin=q4/(dh0)
        uhin=q4/(cls*dh0)
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
        qhdhin=q4/(dh1)
        uhin=q4/(cls*dh1)
c
      endif
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
      if (te1.lt.0.d0) te1=290.d0
      if (jden.eq.'C') then
        dh1=dhn
      else if (jden.eq.'F') then
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
        prescc=dhn*rkb*1.d4+radpint+frp
        ponk=(prescc/rkb)
        pres1=fpressu(te1,dhn,pop)
        dh1=dhn*prescc/(pres1)
      endif
      t=(te1+te0)/2.d0
      dh=(dh0+dh1)/2.d0
c
c***********************************************
c      write(*,*)'***DETERMINE SPACE STEP : DR'
c***********************************************
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      mmod='DIS'
      lmod='SO'
c
c     get new ionisation at t predicted
c     equion calculates a new de
c
c      call equion(t,de,dh)
c
      call totphot (t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
c     Compute distance DR to obtain a given total photon absorption
c     fraction DTAU
      call absdis (dh, fi, dtau, dr, rad0, pop)
   50 if (jgeo.eq.'S') then
        if (jden.eq.'F') then
          if ((frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)).gt.ofac) then
            dpw=dr
            if (pfac.gt.0.0d0) then
              if (bfac.gt.0.0d0) then
                dpw=dis0*((1.1d0**(1/bfac))-1.d0)
              else if (bfac.lt.0.0d0) then
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
        if (radout.lt.huge) then
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
      if (jgeo.eq.'F') then
        dr=1.d15
c        dtau = dtau0
        call absdis (dh, fi, dtau, dr, rad0, pop)
        dis1=dis0+dr
        rad1=0.d0
        wdil1=fdilu(rstar,dis1)
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
      else if (jden.eq.'F') then
        dh1=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
c         frp    = fradpress
c         dh1 = dhn*prescc/(fpressu(te1,dhn,popend)+frp)
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
c        write(*,'(" P Pr",2(1pg12.5))') fpressu(te1,dhn,popend),frp
        frp=fradpress(dr,dh)
        prescc=dhn*rkb*1.d4+radpint+frp
        ponk=(prescc/rkb)
        pres1=fpressu(te1,dhn,popend)
        dh1=dhn*prescc/(pres1)
      endif
      de1=feldens(dh1,popend)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi1=dmin1(1.d0,fidhcc/dh1)
      rww=1.d0
      if (jeq.eq.'C') te1=frad(dis1,tofac,txfac,tpfac,tafac,tbfac,tcfac,
     &tscalen,0.d0,1.d0,0.d0)
      t=(te0+(rww*te1))/(1.d0+rww)
      dh=(dh0+(rww*dh1))/(1.d0+rww)
      de=(de0+(rww*de1))/(1.d0+rww)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
c  Determine dust continuum emission for this region
c
c      if ((grainmode.eq.1).AND.(IRmode.ne.0)) then
c      lmod = 'ALL'
c       call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
c       call dusttemp(t,dh,de, fi, dr,m,phot0)
c      endif
c
c
      lmod='DW'
      wei=1.d0/(1.d0+rww)
c     Form weighted average of POPPRE and POPEND in POP
      call averinto (wei, poppre, popend, pop)
      call localem (t, de, dh)
      call totphot (t, dh, fi, rad0, dr, dv, wdil1, lmod)
c     Calculate effective density of photons per particle
      call zetaeff (fi, dh)
c
c*****************************************************************
c      write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************
c
      if ((jeq.eq.'E').or.((jeq.eq.'P').and.(jthm.eq.'S'))) then
        dton1=1.d33
        call copypop (popend, pop)
        call teequi (te1, tef, de1, dh1, dton1, nmod)
      else if ((jeq.eq.'C').or.(jthm.eq.'T')) then
c
        te1=frad(dis1,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen,0.d0,
     &   1.d0,0.d0)
        call copypop (popend, pop)
c
        call equion (te1, de1, dh1)
c
        i=0
c
c     PJM Aug 2008: Changed IF/GOTO loop to DO WHILE. Loop was testing
c     value of DIFT that was not set by DIFPOP. It now triggers on the
c     value of DIFP that is set by DIFPOP. It is also unnecessary to call
c     copypop(POP,POPZ) at the beginning and the end of the loop.
        i=0
        treap=1.d-12
        difp=1.d0
        do while ((difp.ge.0.001).and.(i.le.4))
c     Store current ionic populations POP in POPZ.
          call copypop (pop, popz)
c     Compute local emissivity in vectors EMIDIF and EMILIN.
          call localem (te1, de1, dh1)
c     Compute mean intensity of radiation in vector TPHOT.
          call totphot (te1, dh1, fi, rad, dr, dv, wdil, lmod)
c     Calculate effective density of photons per particle.
          call zetaeff (fi, dh1)
          call equion (te1, de1, dh1)
c     Determine average logarithmic change in ionic populations.
          call difpop (pop, popz, treap, atypes, difp)
          i=i+1
        enddo
c
        tef=te1
c
      else
c
c     Finite age, non-equilibrium ionisation. Determine time step.
c     Takes into account ionization front velocity.
c
        distim(1,m)=dis1
        distim(2,m)=distim(2,m-1)+(dr/(viofr+(1.d-35*dr)))
        tprop=distim(2,m)
        dtauon=1.0d0
        mmod='DIS'
c
        call absdis (dh, fi, dtauon, disin, 0.0d0, pop0)
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
          call teequi (te1, tef, de1, dh1, dton1, nmod)
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
      else if (jden.eq.'F') then
        dhf=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
c         frp    = fradpress
c         dhf = dhn*prescc/(fpressu(tef,dhn,pop)+frp)
c         write(*,'(" P Pr",2(1pg12.5))') fpressu(tef,dhn,pop),frp
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
        frp=fradpress(dr,dh)
        prescc=dhn*rkb*1.d4+radpint+frp
        ponk=(prescc/rkb)
        pres1=fpressu(tef,dhn,pop)
        dhf=dhn*prescc/(pres1)
c            write(*,'(" P/k Pr/k dPr/k",3(1pg12.5))')
c     &               ponk, radPint/rkb, frp/rkb
            write(*,'(" P/k P_H/k:",2(1pe11.4))') ponk, dhf*tef
      endif
      def=feldens(dhf,pop)
c
c************************************************************
c   write(*,*)'***COMPARES END VALUES WITH PREVIOUS STEP'
c************************************************************
c
      call difpop (pop, popend, trea, atypes, dift)
      call difpop (pop, poppre, trea, atypes, difp0)
      call copypop (pop, popend)
c
      drxp=dr
      dtl=dif(tef,te1)
      dhl=dif(dh1,dhf)
c
c     strict criteria
c
c     difg = dmax1(dift/difma,dtl/dtlma,dhl/dhlma)
c
c     relaxed cirteria
c
c     difg = dhl/dhlma
c
c     lenient
c
      difg=dmax1(dtl/dtlma,dhl/dhlma)
c
c     no criteria
c
c      difg = 0.d0
c
c
      niter=niter+1
      if (difg.le.difmi) then
        texm=te1
        dtaux=dtau
        difmi=difg
      endif
c
      write (*,80) tef,dhl,dtl,dift,difg,difp0,dtau,dr
   80 format(t7,0pf8.0,5(0pf8.3),3(1pg9.2))
c
c***************************************************************
c  write(*,*)'AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED'
c***************************************************************
c
      rww=1.d0
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav=(dis0+(rww*dis1))/(1.d0+rww)
      wei=1.d0/(1.d0+rww)
      call averinto (wei, poppre, popend, pop)
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'F') then
        rad=0.d0
        wdil=fdilu(rstar,disav)
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
      temp1=((de0/dh0)-dedhmi)/dedhma
      temp2=(de1/(de0+1.d-7))*((dh0/dhf)**2)
      if ((difg.ge.1.d0).and.(qhdh.gt.1d1)) then
        weiph=0.5d0
        rmm=1.d5
        if ((niter.gt.11).and.(te1.eq.texm)) then
c     No convergence of T or  H dens after 12 iterations
          write (*,90) dhl,dtl,dift,difg
   90       format(' NO SPATIAL CONVERGENCE :',4(0pf8.3))
          goto 110
        else if (niter.eq.12) then
c     Convergence of T and H dens after 12 iterations
          te1=texm
          dtau=dtaux
          difmi=0.d0
        else if (((temp1.gt.0.01d0).and.(temp2.lt.0.2d0))
     &   .or.(dhl.ge.0.07d0)) then
c     changes in ionization state or Hdens are big
          agmu=dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-((1.d0-agmu)*0.6d0))
          nidh=nidh+1
          if (nidh.lt.3) niter=niter-1
          te1=(te0+dmax1(agmu*te0,dmin1(te0/agmu,te1)))/2.d0
          weiph=0.85d0
          rmm=0.25d0
          dtau=(0.2d0*rmm)*dtau
        else if (((niter.eq.3).or.(niter.eq.6)).and.((jden.ne.'B')
     &   .or.(difg.eq.(dift/difma)))) then
c at 3 or 6 iterations and not isobaric or convergence in T slower than H dens
          agmu=dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-((1.d0-agmu)*0.65d0))
          nidh=max0(0,min0(nidh-2,idnint(1.0d0+dint(niter/3.0d0))))
          te1=(((te1+te0)+te1a)+te0)/4.0
          te1=dmax1(te0*agmu,dmin1(te0/agmu,te1))
          weiph=0.75d0
          rmm=0.4d0
          dtau=(0.2d0*rmm)*dtau
        else
c  None of the above criteria hold
          te1=tef
          weiph=0.85d0
          rmm=0.25d0
          dtau=(0.2d0*rmm)*dtau
        endif
  100   dtau=dtau/2.d0
        mmod='DIS'
        lmod='SO'
        drzi=dr
        drzf=dr
        call totphot (t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
        call absdis (dh, fi, dtau, drzi, rad0, poppre)
        call absdis (dh, fi, dtau, drzf, rad0, popend)
c
c
        dr=dexp((weiph*dlog(drzi+.1d0))+((1.d0-weiph)*dlog(drzf+.1d0)))
        if ((dr.gt.(rmm*drxp)).and.(dtau.gt.1.d-16)) goto 100
        if (dr.lt.drp*1.d-2) then
          dr=drp*1.d-2
          niter=niter+11
          te1=texm
        endif
        if (dif(dr,drxp).le.0.01d0) then
          niter=niter+11
          te1=texm
        endif
        wei=dmin1(1.d0,(1.2d0*dr)/drxp)
        call averinto (wei, popend, poppre, popend)
c
c     Loop back until space step (???) converges.
        goto 50
c
c     Jump here if no spatial convergence.
  110   continue
      else if ((difg.ge.1.d0).and.(qhdh.gt.1d4)) then
c     THIS IS THE SAME TEST AS FOR THE MAIN TEST ABOVE
        dtau=dtau/4.d0
      endif
c
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
c  Determine dust continuum emission for this region
c
      if ((grainmode.eq.1).and.(irmode.ne.0)) then
        call localem (t, de, dh)
        lmod='DW'
        call totphot (t, dh, fi, rad0, dr, dv, wdil1, lmod)
        call dusttemp (t, dh, de, fi, dr, m)
      endif
c
c
c      for diagnostic purposes
c      call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
c
      call zetaeff (fi, dh)
      call localem (t, de, dh)
      call newdif (t, t, dh, fi, rad0, dr, dv, 0.d0, dv, frdw, jtrans)
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
c
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
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
      frp=fradpress(dr,dh)
      radpint=radpint+frp
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'F') then
        rad=0.d0
        wdil=fdilu(rstar,disav)
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
      if ((jeq.ne.'E').and.(jeq.ne.'C').and.(dtoff.gt.0.0)) then
        if (jeq.eq.'P') recscal=frectim(t,0.8d0*de,dh)
        wd=frlum*wdil
        jspot='YES'
        tii0=t
c
        lmod='SO'
        call totphot (t, dh, fi, rad, dr, dv, wd, lmod)
        call evoltem (tii0, t, de, dh, prescc, dtoff, exla, tm00, lut0)
c
        if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      endif
c
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
c         hoii(m) = (fbri(1,ox2)+fbri(2,ox1))/(hbri(2)+epsilon)
c
      if (ox3.ne.0) then
        hoiii(m)=(fmbri(7,ox3)+fmbri(10,ox3))/(hbri(2)+epsilon)
      endif
      if (ox2.ne.0) then
        hoii(m)=(fmbri(3,ox1))/(hbri(2)+epsilon)
      endif
      if (ni2.ne.0) then
        hnii(m)=(fmbri(7,ni2)+fmbri(10,ni2))/(hbri(2)+epsilon)
      endif
      if (su2.ne.0) then
        hsii(m)=fmbri(1,su2)+fmbri(2,su2)/(hbri(2)+epsilon)
      endif
      if (ox1.ne.0) then
        if (ispo.eq.' OI') hsii(m)=fmbri(3,ox1)/(hbri(2)+epsilon)
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
          eng=cphotev(i)
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
        write (*,120) habing,pahqfac
  120   format(' Habing parameter =',1pg12.5,', PAH ion/ISRF =',1pg12.5)
c
c Determine if PAHs exist
c If so, deplete gas of Carbon in PAHs
c
        if (pahactive.eq.0) then
          if (pahend.eq.'H') then
            if (habing.lt.pahlimit) then
              pahactive=1
            endif
          else if (pahend.eq.'Q') then
            if (qhdh.lt.pahlimit) then
              pahactive=1
            endif
          else if (pahend.eq.'I') then
            if (pahqfac.lt.pahlimit) then
              pahactive=1
            endif
          endif
          if (pahactive.eq.1) then
            atom=zmap(6)
            if (clinpah.eq.1) then
              zion(atom)=zion0(atom)*deltazion(atom)*dion(atom)
            else
              zion(atom)=zion0(atom)*deltazion(atom)*(dion(atom)+(1.d0-
     &         dion(atom))*pahcfrac)
            endif
          endif
        endif
      endif
c
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
c
c     write balance for first step or if poll file exists
c
c
      fhi=pop(1,1)
      fhii=pop(2,1)
c
      pollfile='balance'
      inquire (file=pollfile,exist=iexi)
      if (((jbal.eq.'YES').and.(m.eq.1)).or.(iexi)) then
        caller='P5'
        pfx='p4bal'
        pfx=jbfx
        np=5
        call wbal (caller, pfx, np, pop)
      endif
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        caller='P5'
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
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P5'
        pfx='nsou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='sophot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        lmod='SO'
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P5'
        pfx='ssou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='loclphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        lmod='LOCL'
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
        caller='P5'
        pfx='lsou'
        np=4
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
        lmod='ALL'
        call totphot (t, dh, fi, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='IRphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi).and.(irmode.gt.0)) then
        caller='P5'
        pfx='irsou'
        np=5
        wmod='REAL'
        dva=0
        do i=1,infph-1
          pahtest(i)=0.5d0*(photev(i+1)+photev(i))*ev*irphot(i)*dr
        enddo
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, pahtest)
      endif
c
      pollfile='pahphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi).and.(pahactive.eq.1)) then
        caller='P5'
        pfx='pahs'
        np=4
        wmod='REAL'
        dva=0
        do i=1,infph
          pahtest(i)=paheng*pahflux(i)
        enddo
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, pahtest)
      endif
      pollfile='speclocal'
      inquire (file=pollfile,exist=iexi)
      if (iexi) then
        caller='P5'
        pfx='local'
        np=5
        sfx='ph5'
        call newfile (pfx, np, sfx, 3, fn)
        fspl=fn(1:13)
        open (lusp,file=fspl,status='NEW')
        linemod='LAMB'
        spmod='REL'
        call speclocal (lusp, tloss, eloss, egain, dlos, t, dh, de,
     &   pop(1,1), disav, dr, linemod, spmod)
        close (lusp)
      endif
c
c   Infrared output (every 1 step)
c
      pollfile='IRlocal'
      inquire (file=pollfile,exist=iexi)
      if ((irmode.ne.0).and.iexi) then
        if (irfile.eq.' ') then
          caller='P5'
          np=6
          pfx='IRflux'
          sfx='sou'
          call newfile (pfx, np, sfx, 3, fn)
          irfile=fn(1:14)
          open (ir1,file=irfile,status='NEW')
          write (ir1,130) theversion
  130      format('%  Infrared flux per region ',/,
     &            '%  MAPPINGS V ',a8,/
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
        write (ir1,140) m
        write (ir1,*) infph-1
  140    format(/,' Region ',i4.4)
        do i=1,infph-1
          cnphot(i)=ffph(i)+fbph(i)+p2ph(i)
          write (ir1,150) photev(i),cnphot(i)*dr,irphot(i)*dr
        enddo
  150    format(1pg14.7,' ',1pg14.7,' ',1pg14.7)
        close (ir1)
c         do i=1,mxinfph
c            storephot(i)=tphot(i)
c         enddo
c         call totphot(t, dh, fi, rad, dr, dv, wdil, 'LOCL')
c         call wpsoufile(caller,IRneb,
c     &        t,de,dh,
c     &        0.d0,0.d0,qhdh,
c     &        disav,dr,
c     &        0.d0,0.d0,
c     &        1.d0,tphot)
c         do i=1,mxinfph
c            tphot(i)=storephot(i)
c         enddo
      endif
c
c     graincharge output
c
      if (grainmode.eq.1) then
        pollfile='graincharge'
        inquire (file=pollfile,exist=iexi)
        if (iexi) then
          if (chargefile.eq.' ') then
            caller='P5'
            np=5
            pfx='grpot'
            sfx='ph5'
            call newfile (pfx, np, sfx, 3, fn)
            chargefile=fn(1:13)
            open (ir1,file=chargefile,status='NEW')
            write (ir1,160) theversion
  160      format('%  grain charge in each region ',/,
     &            '%  MAPPINGS V ',a8,/
     &            '%  given as:',/,
     &            '%  m,Av. distance,Temp,dh,de,',/,
     &            '%  grain charge(allsizes) (graphite),',/,
     &            '%  grain charge(allsizes) silicate',/,
     &            '%')
            write (ir1,170) (grainrad(i),i=mindust(1),maxdust(1))
  170      format('gra radii',10(1pg11.4))
            write (ir1,*)
            write (ir1,180) (grainrad(i),i=mindust(2),maxdust(2))
  180      format('sil radii',10(1pg11.4))
            close (ir1)
          endif
          open (ir1,file=chargefile,status='OLD',access='APPEND')
          write (ir1,190) m,disav,t,dh,de
          write (ir1,200) (grainpot(i,1),i=mindust(1),maxdust(1))
          write (ir1,210) (grainpot(i,2),i=mindust(2),maxdust(2))
          close (ir1)
  190      format(/,i3,4(1pg14.5))
  200      format('gra ',30(1pg11.4))
  210      format('sil ',30(1pg11.4))
        endif
      endif
c
      open (luop,file=filna,status='OLD',access='APPEND')
      open (lupf,file=filnb,status='OLD',access='APPEND')
      open (lups,file=filnc,status='OLD',access='APPEND')
      if (jall.eq.'YES') then
        open (lusl,file=newfil,status='OLD',access='APPEND')
      endif
c
  220 format(i4,1pg11.4,1x,1pg10.3,13(1pg14.7))
      write (luop,220) m,t,dlosav,disav,dis0-remp,dr,dh,de,fhi,fhii,
     &qhdnin,qhdn,caseab(1),caseab(2),difg,dtau
c
c      write(luop, 461) m, t, dlosav, dis1, disav, dr, fi, dh, de
c      write(luop, 1461) fhi, fhii, zetae, qhdh,hoiii(m),hoii(m),
c     &                  hnii(m),hsii(m),hbeta,caseab(1)
c
      if (jeq.eq.'F') write (luop,230) fthick,treach,dtco,dton,dtoff
  230 format(t71,5(1pg9.2))
c
  240 format(t98,5(1pg9.2))
      if (jeq.eq.'P') write (luop,240) recscal,dtoff
c
  250 format('    Case A-B (HI, HeII): ',2(1pg11.3))
      write (*,250) caseab(1),caseab(2)
c
c
      if (jgeo.eq.'S') then
  260 format(' #',t7,'<Te>',t15,'<DLOS>',t23,'<R>',t33,'dr',t43
     & ,'<nH>',t53,'<FHI>',t63,'<Q(N)>',t73,'OIII/HB',/
     & i4,0pf8.0,1pg9.1,x,6(1pg9.3,x))
        write (*,260) m,t,dlosav,disav,dr,dh,fhi,qhdn,hoiii(m)
      else
  270 format(' #',t7,'<Te>',t15,'<DLOS>',t23,'<R>',t33,'dr',t43
     & ,'<nH>',t53,'<FHI>',t63,'<Q(N)>',t73,'OIII/HB',/
     & i4,0pf8.0,1pg9.1,x,6(1pg9.3,x))
        write (*,270) m,t,dlosav,disav,dr,dh,fhi,qhdn,hoiii(m)
      endif
c
c
      if (jiel.eq.'YES') then
        do i=1,ieln
          open (luions(i),file=filn(i),status='OLD',access='APPEND')
c
  280  format(1x,i4,3(', ',1pg12.5),31(', ',1pg12.5))
          if (i.eq.1) then
            write (luions(i),280) m,disav,dr,t,(pop(j,iel1),j=1,
     &       maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),280) m,disav,dr,t,(pop(j,iel2),j=1,
     &       maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),280) m,disav,dr,t,(pop(j,iel3),j=1,
     &       maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),280) m,disav,dr,t,(pop(j,iel4),j=1,
     &       maxion(iel4))
          endif
c
          close (luions(i))
        enddo
c
      endif
c
c
      if (jiem.eq.'YES') then
        do i=1,ieln
          open (luemiss(i),file=filnem(i),status='OLD',access='APPEND')
          idx=iemidx(i)
          nt=nfmtrans(idx)
c
  290  format(1x,i4,4(', ',1pg12.5),31(', ',1pg12.5))
c
          write (luemiss(i),290) m,disav,dr,t,hbeta,(fmbri(itr,idx),itr=
     &     1,nt)
c
          close (luions(i))
        enddo
c
      endif
c
c
c      nt1 = time()
c
c      dtarr = dtime(tarray)
c
c 100  format('Step time (map3,sys): ',2(f8.1,1x),/)
c      write(*,100) tarray(1),tarray(2)
c
c      nt0 = nt1
c
      if (jall.eq.'YES') then
        write (lusl,300)
        write (lusl,310) m,t,de,dh,fi,dvoluni,dr,disav,(disav-remp)
  300    format(//' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  fill. Fact. ',
     &                '  dVol Unit.  ',
     &                '    dR (cm)   ',
     &                ' Dist.Ave.(cm)',
     &                ' D - Remp.(cm)')
  310    format(i3,8(1pg14.7)/)
      endif
c
      wmod='PROP'
      call wmodel (lupf, t, de, dh, disav, wmod)
      wmod='LOSS'
      call wmodel (lups, t, de, dh, disav, wmod)
c
c     experiment with timescales
c
      ctime=fcolltim(de)
      rtime=frectim2(de)
      ptime=fphotim()
c
      eqtime=(1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime=dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime=(1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime=dabs(1.d0/stime)
c
c 700  format(' Timescales: ',5(1pg11.4,1x))
c      write(*,700) ctime,rtime,ptime,eqtime,stime
c
      if (jall.eq.'YES') then
        call wionabal (lusl, pop)
        call wionabal (lusl, popint)
      endif
c
      close (luop)
      close (lupf)
      if (jall.eq.'YES') close (lusl)
      close (lups)
c
c     ***END OUTPUT
c
c**********************************************************************
c
c     ***TEST ENDING CONDITIONS
c
      do 320 n=1,ionum
c      find first (lowest energy) cross section that matches species
        if ((atpho(n).eq.ielen).and.(ionpho(n).eq.jpoen)) goto 330
  320 continue
c sigpho contains new threshold values if verner in force.
  330 taux=sigpho(n)*popint(jpoen,ielen)
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) goto 340
      if ((jend.eq.'C').and.(telast.le.tend)) goto 340
      if ((jend.eq.'D').and.(taux.ge.tauen)) goto 340
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 340
      if (jend.eq.'F') then
        if (jpoen.gt.maxion(ielen)) then
          popinttot=0
          do i=1,maxion(ielen)
            popinttot=popinttot+popint(i,ielen)
          enddo
          if (popinttot.ge.colend) goto 340
        else
          if (popint(jpoen,ielen).ge.colend) goto 340
        endif
      endif
c
c Check Visual extinction condition or end when HII<1%
c
      if (jend.eq.'H') then
        a_v=dustav*dustint
        if ((a_v.gt.a_vend).or.(de/dh.lt.1e-5)) goto 340
      endif
c
      if (((((de1/dh1)-dedhmi)/dedhma).le.exla).and.(grainmode.eq.0))
     &goto 340
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 340
c
c
c    ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c
      call copypop (popend, poppre)
c
      tep=te0
      te0=tef
c     Record electron temperatures in THIST for quadratic extrapolation
c     when m > 3.
      thist(1)=thist(2)
      thist(2)=tep
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
c**********************************************************************
c
c      write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
c
c**********************************************************************
c
      m=m+1
c      write(*,'(/)')
c
c     stop things getting silly....
c
      if (m.gt.mma) goto 340
c
c     loop back
c
      goto 20
c
c
  340 continue
c
c********************************************************************
c
c     write(*,*)'MODEL ENDED ; OUTPUT RESULTS **************'
c
c********************************************************************
c
      open (luop,file=filna,status='OLD',access='APPEND')
      open (lupf,file=filnb,status='OLD',access='APPEND')
      write (luop,350)
  350 format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &'TAUXF',t58,'Ending')
      write (luop,360) t,dis1,frx,taux,jend
  360 format(' ===========' ,t18,0pf8.0,3(1pg10.3),4x,a4)
      if ((jeq.ne.'E').and.(frlum.gt.0.0)) write (luop,370) frlum
  370 format(/' Final fractional luminosity of source after ',
     &'turn off :',1pg10.3)
      if ((radout.lt.huge).and.(jgeo.eq.'S')) write (luop,380) radin,
     &radout
  380 format(/' NB. :::::::::::::::Integration through the line',
     &'of sight for a ring aperture of radii :',2(1pg10.3))
      write (lupf,390) tempre,tspr
  390 format(//t4,'Preionisation conditions for step#0 ',
     &'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
      write (lupf,400) (elem(i),i=1,atypes)
  400 format(' ',t8,16(4x,a2,4x)/)
      maxio=0
      do j=1,atypes
        if (maxio.le.maxion(j)) maxio=maxion(j)
      enddo
      do j=1,maxio
        write (lupf,410) rom(j),(ppre(j,i),i=1,atypes)
  410    format(a7,t8,16(1pg10.3))
      enddo
c
      if (jsou.eq.'YES') then
c
        caller='P4'
        pfx=jpfx
        np=5
        wmod='NFNU'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c     Down stream photon field
c
        caller='P4'
        pfx=jpfx
        np=5
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
      endif
c
c
      if (jspec.eq.'YES') then
c
c     write nebula spectrum in nuFnu (ergs/s/cm2/sr)
c
c
        lmod='NEBL'
        call totphot (te1, dh1, fi, dis0, dr, 0.0d0, wdil0, lmod)
        caller='P4'
        pfx=jpfx
        np=5
        wmod='NFNU'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
c     Down stream nebula only photon field
c
        caller='P4'
        pfx=jpfx
        np=5
        wmod='REAL'
        dva=0
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, qhdh,
     &   disav, dr, 0.d0, 0.d0, 1.d0, tphot)
c
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
c
      open (lusp,file=filnd,status='OLD',access='APPEND')
      linemod='LAMB'
      spmod='REL'
      call spec2 (lusp, linemod, spmod)
      close (lusp)
c
      write (*,420) fnam
  420 format(//' Output created &&&&&& File : ',a/)
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
      subroutine p5fheader (newfil, banfil, fnam, filna, filnb, filnc,
     &filnd, filn, filnem, luop, lupf, lusl, lups, lusp, luions,
     &luemiss, dhn, fin)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 dhn,fin
      real*8 dht,zi(mxelem)
      integer*4 i, ie, j, at, io
      integer*4 idx, plen, nt, nl, itr
      integer*4 luop,lups,lupf,lusl,lusp
      integer*4 luions(4), luemiss(4)
      character*24 abundtitle
      character fn*32
      character pfx*16,sfx*4
      character newfil*32, filna*32, banfil*32, fnam*32, filnb*32
      character filnc*32, filnd*32, filn(4)*32, filnem(4)*32
c
      real*8 densnum
c
      dht=densnum(dhn)
c
      fn=' '
      pfx='photn'
      sfx='ph5'
      call newfile (pfx, 5, sfx, 3, fn)
      filna=fn(1:13)
      fnam=filna
c
c
      fn=' '
      pfx='phapn'
      sfx='ph5'
      call newfile (pfx, 5, sfx, 3, fn)
      filnb=fn(1:13)
c
      fn=' '
      pfx='phlss'
      sfx='ph5'
      call newfile (pfx, 5, sfx, 3, fn)
      filnc=fn(1:13)
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
        sfx='ph5'
        call newfile (pfx, 6, sfx, 3, fn)
        newfil=fn(1:14)
      endif
c
      open (luop,file=filna,status='NEW')
      open (lupf,file=filnb,status='NEW')
      open (lups,file=filnc,status='NEW')
      open (lusp,file=filnd,status='NEW')
      if (jall.eq.'YES') open (lusl,file=newfil,status='NEW')
c
c     ***WRITES ON FILE INITIAL PARAMETERS
c
      write (luop,10) theversion,runname,banfil,fnam
      write (lupf,10) theversion,runname,banfil,fnam
      write (lups,10) theversion,runname,banfil,fnam
      write (lusp,10) theversion,runname,banfil,fnam
      if (jall.eq.'YES') then
        write (lusl,10) theversion,runname,banfil,fnam
      endif
   10 format(
     &' ###################################'
     &,'####################################',/
     & '  Photoionisation model P5, MAPPINGS V ',a8,/
     &' ###################################'
     &,'####################################',/
     &' (Diffuse Field ,Outward Integration, Radiation Pressure)'//
     &' Run   : ',a/
     &' Input : ',a/
     &' Output: ',a/)
      close (lusp)
      write (lupf,20)
   20 format(' This file contains the plasma properties for ',
     &'each model step,'/
     &' and the pre-ionisation conditions used'/ )
   30 format(' This file contains the ionisation fraction for ',
     &'all atomic elements as a function of distance,'/
     &' and the abundances used'// )
      if (jall.eq.'YES') then
        write (lusl,30)
        write (lusl,40) (elem(i),zion(i),i=1,atypes)
   40    format('ABUND:',16(a2,':',1pg10.3,1x)//)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     photn file header
c
      abundtitle=' Initial Abundances :'
      do i=1,atypes
        zi(i)=zion0(i)*deltazion(i)
      enddo
      call dispabundances (luop, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jeq.eq.'E') then
        write (luop,50) dtau0
   50    format(/' Model Summary :',/
     &  ,'    Mode  : Thermal and Ionic Equilibrium :',/
     &  ,'    dTau  :',0pf8.5)
      else if (jeq.eq.'F') then
        write (luop,60) telap,tlife,dtau0
   60    format(/' Model Summary :',/
     &  ,'    Mode        :  Non-Equilibrium + Finite Source Life :',/
     &  ,'    Age         :',1pg10.3,' sec',/
     &  ,'    Source Life :',1pg10.3,' sec',/
     &  ,'    dTau        :',0pf8.5)
      else
        write (luop,70) telap,dtau0
   70    format(/' Model Summary :',/
     &  ,'    Mode        :  Equilibrium + Source Switch off :',/
     &  ,'    Switch off  :',1pg10.3,' sec',/
     &  ,'    dTau        :',0pf8.5)
      endif
      if (usekappa) then
   80    format(/' Kappa Electron Distribution Enabled :'/
     &           '    Electron Kappa : ',1pg11.4)
        write (luop,80) kappa
      endif
c
   90 format(/' Radiation Field :',//,
     &    t5,' Rsou.',t20,' Remp.',t35,' Rmax',t50,' <DILU>'/
     &    t5, 1pg10.3,t20,1pg10.3,t35,1pg10.3,t50,1pg10.3/
     &    t5,' <Hdens>',t20,' <Ndens>',t35,' Fill Factor',/
     &    t5, 1pg10.3,t20,1pg10.3,t35,0pf9.6,/
     &    t5,' FQ tot',t20,' Log FQ',t35,' LQ tot',t50,' Log LQ',/
     &    t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     &    t5,' Q(N) in',t20,' Log Q(N)',t35,' <Q(N)>',t50,' Log<Q(N)>',/
     &    t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     &    t5,' Q(H) in',t20,' Log Q(H)',t35,' <Q(H)>',t50,' Log<Q(H)>',/
     &    t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3,/
     &    t5,' U(H) in',t20,' Log U',t35,' <U(H)>',t50,' Log<U>',/
     &    t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t50,0pf8.3)
c
      write (luop,90) rstar,remp,rmax,wdpl,dhn,dht,fin,qht,dlog10(qht),
     &astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),qhdnav,
     &dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),uhin,
     &dlog10(uhin),uhav,dlog10(uhav)
c
c
c      write(luop, 430)
c 430  format(//' Description of Photon Source :'/
c     &   ' #####################################################')
c
  100 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
      write (luop,100)
  110 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
      write (luop,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (turbheatmode.eq.1) then
  120 format(/,'  Micro-Turbulent Dissipation Enabled, Mach = ',1pg10.3)
        write (*,120) admach
      endif
c
      if (grainmode.eq.1) then
        if (gdist.eq.'M') then
          write (luop,150) galpha,amin(1),amax(1),amin(2),amax(2)
        else if (gdist.eq.'P') then
          write (luop,160) galpha,amin(1),amax(1),amin(2),amax(2),
     &     graindens(1),graindens(2),bgrain,yinf
        endif
c
  130    format(//'  Projected dust area/H atom (graphite,silicate):'/
     &            ' ::::::::::::::::::::::::::::::::::::::::::::::::'/
     &            '  ',1pg12.4,1x,' (cm^2)',1x,1pg12.4,1x,' (cm^2)' )
        write (luop,130) siggrain(1),siggrain(2)
c
  140    format( /'  Composition Ratios:  '/
     &            ' ::::::::::::::::::::::::::::::::::::::::::::::::'/
     &            '   Gas/H     :',1pg12.4/
     &            '   Dust/Gas  :',1pg12.4/
     &            '   PAH/Gas   :',1pg12.4/
     &            ' ::::::::::::::::::::::::::::::::::::::::::::::::'/
     &           ' (ratios by mass) ')
c
        write (luop,140) grainghr,graindgr,grainpgr
c
      endif
  150 format(//,' Dust Parameters :',//,
     &         ' MRN size distribution :',/,
     &         ' alpha   ',' Graphite min:    max: ',
     &         '   Silicate min:    max: '/,1pg10.3,
     &         2(4x,1pg10.3,1pg10.3))
c
  160 format(//,' Dust Parameters :',//,
     &         ' Powerlaw size distribution'/,
     &         ' alpha   ',' Graphite min:    max: ',
     &         '   Silicate min:    max: '/,1pg10.3,
     &         2(4x,1pg10.3,1pg10.3),
     &         ' Graph dens  Sil dens  Bgrain    Yinf ',/,
     &         1pg10.3,2x,3(1pg10.3))
c
  170 format(//' Model Parameters:',//
     & ,' ',t3,'Jden',t9,'Jgeo',t15,'Jend',t20,'Ielen',
     & t26,'Jpoen',t35,'Fren',/
     & ,'  ',x,3(a4,2x),2(i2,4x),0pf6.4,//
     & ,' ',t3,'Tend',t9,'DIend',t19,'TAUen',t28,'Jeq',t35,'Teini',/
     & ,0pf6.0,2(1pg10.3),x,a4,x,0pf7.1)
      write (luop,170) jden,jgeo,jend,ielen,jpoen,fren,tend,diend,tauen,
     &jeq,tm00
c
      write (luop,180)
  180 format(//' Description of each space step:')
c
c      write(luop, 1440)
c 440  format(/' #',t7,'Teav',t15,'DLOSav',t24,'DISout',t33,'DISav',t42,
c     &  'dr',t51,'Fill.F.',t60,'H.Den.',t69,'El.Den.')
c 1440 format(t7,'FHI',t15,'FHII',t24,'ZETAef',t33,'QHDN',t42,'[OIII]'
c     & ,t51,'[OII]',t60,'[NII]',t69,'[SII]',t78,'Caseab H',' Caseab He')
c'QHDNin',5x,
c
c       write(luop, 461) m,t,dlosav,difg,dtau,disav,(dis0-remp),
c      &       dr,dh,de,fhi,fhii,
c      &       qhdnin, qhdn, caseab(1), caseab(2)
c
      if (jgeo.eq.'S') then
  190 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t28,'<R>',t42,'Delta R',
     & t56,'dr'
     & ,t71,'H.den',t85,'el.den',t98,'F_HI',t113,'F_HII',t126,'QHDNin',
     & t140,'<QHDN>',t155,'Caseab H',t169,'Caseab HeII'
     & ,t183,'H.P.',t197,'QPAH/ISRF')
        write (luop,190)
      else
  200 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t28,'<X>',t42,'Delta X',
     & t56,'dx'
     & ,t71,'H.den',t85,'el.den',t98,'F_HI',t113,'F_HII',t126,'Q(N)in',
     & t140,'<QHDN>',t155,'Caseab H',t169,'Caseab HeII'
     & ,t183,'H.P.',t197,'QPAH/ISRF')
        write (luop,200)
      endif
      if (jeq.eq.'F') write (luop,210)
  210 format(t73,'Fthick',t82,'Treach',t91,'DTCO',t100,'Dton',t109,
     &'DToff')
      if (jeq.eq.'P') write (luop,220)
  220 format(t100,'RECscal',t109,'DToff')
      close (luop)
      write (*,230) fnam
  230 format(//' Output in file : ',a//)
c
      if (jgeo.eq.'S') then
  240 format(' #',t7,'<Te>',t15,'<DLOS>',t24,'<R>',t33,'dr',t42
     & ,'<nH>',t51,'<FHI>',t60,'<Q(N)>',t69,'OIII/HB')
        write (*,240)
      else
  250 format(' #',t7,'<Te>',t15,'<DLOS>',t24,'<X>',t33,'dx',t42
     & ,'<nH>',t51,'<FHI>',t60,'<Q(N)>',t69,'OIII/HB')
        write (*,250)
      endif
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
          filn(i)=fn(1:(elem_len(ie)+3+5))
        enddo
c
        do i=1,ieln
          open (luions(i),file=filn(i),status='NEW')
        enddo
c
c
        do i=1,ieln
c
          write (luions(i),'(" Run: ",a80)') runname
c
  280    format('[1] <X>, [2] DeltaX, [3] dX, [4] <T>, ',
     & 31('[',i2,'] ',a6,', '))
          if (i.eq.1) then
            write (luions(i),280) (j+4,rom(j),j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),280) (j+4,rom(j),j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),280) (j+4,rom(j),j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),280) (j+4,rom(j),j=1,maxion(iel4))
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
          filnem(i)=fn(1:(plen+3+5))
c
          open (luemiss(i),file=filnem(i),status='NEW')
c
c
          write (luemiss(i),*) 'Run: ',runname
          write (luemiss(i),*) 'Emissivity (erg/cm^3/s/sr) for Ion: ',
     &     elem(at),rom(io)
  290   format(
     &  ' Step [1], <X> [2], DeltaX [3], dX [4], <T> [5], <HB> [6]',
     &  40(',',1pg12.6,'[',i2,']'))
          write (luemiss(i),290) (1.0d8*fmlam(itr,idx),itr+6,itr=1,nt)
c
          close (luemiss(i))
        enddo
c
      endif
c
      return
c
      end
