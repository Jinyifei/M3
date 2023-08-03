cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine totphot (t, dh, fi, rad, dr, dv, wd, lmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM VECTOR : TPHOT REPRESENTING THE MEAN INTENSITY
c     OF RADIATION : JNU (AT THE POINT CONSIDERED) USING
c     SOURCE VECTOR : SOUPHO(=INU AT RDIST=RSTAR).
c
c     ADDS DIFFUSE FIELD CONTAINED IN VECTORS :
c          UPDIF,DWDIF,UPLIN AND DWLIN
c     ALSO ADDS LOCAL EMISSIVITY VECTORS : EMIDIF,EMILIN
c
c     NOTE : UNITS FOR THE VECTORS EMIDIF & EMILIN ARE
c              IN NUMBERS OF PHOTONS
c              (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c
c     GENERALLY SUBR. LOCALEM  NEEDS TO BE PREVIOUSLY CALLED
c     GENERALLY SUBR. NEWDIF  NEEDS TO BE PREVIOUSLY CALLED
c
c     ABSORPTION : TAU  IS DETERMINED USING THE INTEGRATED
c     COLUMN DENSITY OF MATTER : POPINT  FOR THE SOURCE VECTOR
c     AND THE COLUMN DENS. WITHIN THE SPACE STEP ITSELF FOR ALL
c     VECTORS (FROM R UP TO R+DR  USING POPUL. IN : POP(6,11))
c
c     THE GEOMETRICAL DILUTION FACTOR : WDIL  NEEDS TO BE GIVEN.
c     GENERALLY : WDIL=(RDIST**2)/(4*RSTAR**2)  (SEE FUNCTION FDILU).
c     IN THE PLANE-PARALLEL CASE, WDIL HAS TO BE EQUAL TO 0.5.
c
c     THE DILUTION FACTOR WITHIN THE SPACE STEP ITSELF IS
c     DETERMINED FROM THE RADIUS OF CURVATURE : RAD  DEFINED
c     AT THE INNER LIMIT OF THE SPACE STEP : DR
c     FOR THE UPSTREAM VECTOR, IT CORRESPONDS (BY CONVENTION FOR
c     AN INWARD FLUX) TO A CONCENTRATION FACTOR
c     (GIVE A VALUE FOR RAD OF ZERO FOR PLANE PARALLEL GEOMETRY)
c     CALL SUBR. RESTRANS,CASAB
c
c
c     LMOD='DW'   : ADDS SOURCE VECTOR AND DOWNSTREAM VECTOR.
c     LMOD='UP'   : ADDS SOURCE VECTOR AND UPSTREAM VECTOR.
c     LMOD='ALL'  : ADDS ALL COMPONENTS.
c     LMOD='SO'   : ADDS ONLY SOURCE VECTOR (ON THE SPOT APPROX.)
c                   (IN THIS CASE SUBR. LOCALEM IS UNNECESSARY).
c     LMOD='CAB'  : THE SUBROUTINE ONLY DETERMINE CASE A,B
c                   (NO PHOTON FIELD AT ALL).
c     LMOD='LOCL' : ONLY BUILD LOCAL DIFFUSE FIELD.
c     LMOD='NEBL' : ONLY INTEGRATED NEBULA EMISSION (IE NO SOURCE)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     T    = Electron tempeature in K.
c     DH   = Total Hydrogen density
c     FI   = Filling factor
c     RAD  = Distance from photoionization source in cm
c     DR   = Width of space step in cm
c     DV   = Width of velocity step in cm s^-1
c     LMOD = Code for radiation components to include.
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t, dh, fi, rad, dr, dv, wd
      character lmod*4
      real*8 dadw,daup,den,dismul
      real*8 drh,dvh,dwf,emf,energ
      real*8 swdil
      real*8 totem,upf,wid
      real*8 tauso, tau, sigmt, curad, wlo, ulo
      real*8 tau0, dusttau0, srcf0
      real*8 dwdil,rhoion,phots,z,p,f,es,casab,fb,dfb
      real*8 pathlength
      real*8 localflux
      real*8 escape, escape0
      real*8 rf0,rf1
      real*8 tlin0,tlin1,thlin0,thlin1
      real*8 txhlin0,txhlin1,thelin0,thelin1
c
c      real*8 ti, dustloss
c
      real*8 dusttau, dustsigmat
      real*8 dusttrans, dusttotem
      real*8 downflux, upflux, retain
c
      real*8 sum1,sum2,sum3,tlocalem
      real*8 dsum1,dsum2,dsum3
c
      integer*4 inl,atom,ion,ilower
      integer*4 bincount
      integer*4 line,series
c
c           Functions
c
      real*8 densnum, localout, transferout
      real*8 fdilu,fdismul,fbowen
c
c     constants
c
      rhoion=densnum(dh)
c
c       write(*,*) 'Totphot:',t, dh, fi, rad, dr, dv, wd, lmod
c
      if ((lmod(1:2).ne.'DW').and.(lmod(1:2).ne.'UP')
     &.and.(lmod(1:3).ne.'ALL').and.(lmod(1:2).ne.'SO')
     &.and.(lmod(1:4).ne.'LOCL').and.(lmod(1:4).ne.'NEBL')
     &.and.(lmod(1:4).ne.'CONT').and.(lmod(1:3).ne.'CAB')) then
        write (*,10) lmod(1:4)
   10     format(//,'MODE IMPROPERLY DEFINED IN TOTPHOT :',a4)
        stop
      endif
      if ((fi.gt.1.0d0).or.(fi.le.0.0d0)) then
        write (*,20) fi
   20   format(//,'INCONSISTENT FILLING FACTOR IN TOTPHOT :',1pg10.3)
        stop
      endif
c
      if ((wd.gt.0.5d0).or.(wd.lt.0.0d0)) then
        write (*,30) wd
   30   format(/,'GEOMETRICAL DILUTION FACTOR IS IMPROPER IN',
     &     ' SUBR. TOTPHOT  :',1pg10.3)
        stop
      endif
c
c     Compute distance multipliers for resonance lines and
c     determine case A,B.
c
      drh=dr*0.5d0
      dvh=dv*0.5d0
c
      if (lmod.eq.'CONT') goto 100
c
      do line=1,nxr3lines
        xr3lines_emilin(2,line)=1.d0
        if (xr3lines_emilin(1,line).gt.epsilon) then
          if (xr3lines_frac(line).gt.epsilon) then
            ilower=xr3lines_imap_i(line)
            if (ilower.eq.1) then
              atom=xr3lines_at(line)
              ion=xr3lines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xr3lines_gf(line)
                es=xr3lines_egij(line)
            xr3lines_emilin(2,line) =
     &       fdismul(t,dh,drh,dvh,atom,ion,es,f)
              endif
            endif
          endif
        endif
      enddo
c
      do line=1,nxrllines
        xrllines_emilin(2,line)=1.d0
        if (xrllines_emilin(1,line).gt.epsilon) then
          if (xrllines_frac(line).gt.epsilon) then
            ilower=xrllines_imap_i(line)
            if (ilower.eq.1) then
              atom=xrllines_at(line)
              ion=xrllines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xrllines_gf(line)
                es=xrllines_egij(line)
            xrllines_emilin(2,line) =
     &       fdismul(t,dh,drh,dvh,atom,ion,es,f)
              endif
            endif
          endif
        endif
      enddo
c
C     do line=1,xlines
C       emilin(2,line)=1.d0
C       if (xbr(line).gt.epsilon) then
C         z=zion(xrat(line))
C         p=pop(xion(line),xrat(line))
C         if ((z*p.ge.pzlimit)) then
C           f=xfef(line)/xbr(line)
C           es=xejk(line)*ev
C           emilin(2,line)=fdismul(t,dh,drh,dvh,xrat(line),xion(line),
C    &       es,f)
C         endif
C       endif
C     enddo
c
      fb=0.d0
      if (zmap(2).gt.0) then
        fb=fbowen(t,dvh)
      endif
c
      atom=1
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nhseries
        do line=1,nhlines
          f=hydrogf(line,series)
          hydlin(2,line,series)=1.d0
          es=(lmev/hlambda(line,series))*ev
          if ((z*p.ge.pzlimit).and.(hbin(line,series).ne.0)
     &     .and.(series.eq.1)) then
            hydlin(2,line,series)=fdismul(t,dh,drh,dvh,atom,ion,es,f)
          endif
        enddo
      enddo
c
c     Helium.
      atom=2
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nheseries
        do line=1,nhelines
          f=hydrogf(line,series)
          hellin(2,line,series)=1.d0
          es=(lmev/helambda(line,series))*ev
c
          if ((z*p.ge.pzlimit).and.(hebin(line,series).ne.0)
     &     .and.(series.eq.1)) then
            hellin(2,line,series)=fdismul(t,dh,drh,dvh,atom,ion,es,f)
          endif
        enddo
      enddo
c
c     Heavier elements.
      do atom=3,atypes
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
        do series=1,nxhseries
          do line=1,nxhlines
            f=hydrogf(line,series)
            xhydlin(2,line,series,atom)=1.d0
            es=(lmev/xhlambda(line,series,atom))*ev
c
            if ((z*p.ge.pzlimit).and.(xhbin(line,series,atom).ne.0))
     &       then
              xhydlin(2,line,series,atom)=fdismul(t,dh,drh,dvh,atom,ion,
     &         es,f)
            endif
          enddo
        enddo
      enddo
c
c     call casab with hydrogen Lyman gamma
c
      if (zmap(1).ne.0) then
        line=3
        series=1
        caseab(1)=casab(1,line,series)
      endif
c
c     call casab with Helium Lyman gamma
c
      if (zmap(2).ne.0) then
        line=3
        series=1
        caseab(2)=casab(2,line,series)
      endif
c
c     Assume heavy series are case A
c
      do atom=3,atypes
        caseab(atom)=0.0d0
      enddo
c
      if (lmod.eq.'CAB') return
c
  100 continue
c
c     Set mode and internal dilution factors.
c
      upf=1.0d0
      dwf=1.0d0
      emf=1.0d0
      if ((lmod.ne.'UP').and.(lmod.ne.'ALL')) upf=0.0d0
      if ((lmod.ne.'DW').and.(lmod.ne.'ALL')) dwf=0.0d0
c
      if (lmod.eq.'LOCL') dwf=1.0d0
      if (lmod.eq.'NEBL') dwf=1.0d0
c
      if (lmod.eq.'SO') emf=0.0d0
C     if (lmod.eq.'NEBL') emf=0.0d0
C     if (lmod.eq.'CONT') emf=0.0d0
c
      curad=rad
      dadw=1.0d0
      daup=1.0d0
c
      if ((curad.gt.0.d0).and.(jgeo.eq.'S')) then
        wlo=(2.0d0*dlog(curad))-(2.0d0*dlog(curad+dr))
        ulo=(2.0d0*dlog(curad+dr))-(2.0d0*dlog(curad))
        dadw=dexp(wlo)
        daup=dexp(ulo)
      endif
c
      dwdil=1.d0
      if (jgeo.eq.'F') then
        dilradpho=0.d0
        dwdil=fdilu(rshock,rad-dilradpho)
      endif
c
      swdil=1.d0
      if (jgeo.eq.'P') then
        swdil=0.5d0
      endif
c
c      write(*,*) 'mode ',lmod
c      write(*,*) 'AB H He',caseab(1), caseab(2)
c      write(*,*) 'fracs',upf,dwf,emf
c      write(*,*) 'da',daup,dadw
c      write(*,*) 'dils ',wd,dwdil,swdil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Form vector TPHOT bin by bin.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      bincount=0
c       dustloss = 0.0d0
c       ti = 0.0d0
      sum1=0.d0
      dsum1=0.d0
      sum2=0.d0
      dsum2=0.d0
      sum3=0.d0
      dsum3=0.d0
c
      tlocalem=0.d0
c
      totsum=0.d0
      srcsum=0.d0
c Main loop over photon energies.
      do 50 inl=1,infph-1
        den=cphotev(inl)
c 0.5d0*(photev(inl)+photev(inl+1))
        wid=(photev(inl+1)-photev(inl))
        energ=den*ev
        tauso=0.d0
        sigmt=0.d0
        dustsigmat=0.d0
        call crosssections (inl, tauso, sigmt, dustsigmat)
c total crosssection, including dust
        sigmt=dh*sigmt
c just the dust component
        dustsigmat=dh*dustsigmat
        pathlength=dr*fi
        tau0=pathlength*sigmt
        escape0=transferout(1.d0,tau0)
        dusttau0=pathlength*dustsigmat
c     XSEC is use to determine "skip bin" for photoionization
c     calculations only, don't include dust in these.
        xsec(inl)=(sigmt-dustsigmat)*fi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Add diluted source field and integrated diffuse (source) field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Downstream, source attenuated by local opacity.
c
        srcf0=0.d0
c
        rf0=0.d0
        rf1=0.d0
c
        tlin0=0.0d0
        tlin1=0.0d0
c
        thlin0=0.0d0
        thlin1=0.0d0
c
        thelin0=0.0d0
        thelin1=0.0d0
c
        txhlin0=0.0d0
        txhlin1=0.0d0
c
        tau=tau0
        dusttau=dusttau0
        dusttrans=0.d0
c
        if ((lmod.ne.'LOCL')) then
c
c      Source field if not local and not nebula only
c
          if (lmod.ne.'NEBL') then
            srcf0=wd*swdil*soupho(inl)*dexp(-tauso)
          endif
c
c     Integrated diffuse field and source
c
          if (lmod.eq.'CONT') then
          downflux=dwf*dadw*dwdil*dwdifcont(inl)+srcf0
          upflux=upf*daup*updifcont(inl)
          else
          downflux=dwf*dadw*dwdil*dwdif(inl)+srcf0
          upflux=upf*daup*updif(inl)
          endif
c
          escape=transferout(1.d0,tau)
c
          rf0=downflux*escape
          rf1=upflux
c
          retain=1.d0-transferout(1.d0,dusttau)
          dusttrans=dusttrans+downflux*retain
c
          sum1=sum1+(downflux)*(1.d0-escape)*wid*evplk
          dsum1=dsum1+(downflux)*retain*wid*evplk
        endif
c
        if (lmod.eq.'CONT') goto 110
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now Add Integrated Resonance Lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Metal Resonance Lines
c
        if (lmod.ne.'LOCL') then
c
          do line=1,nxr3lines
            if ((xr3lines_bin(line).eq.inl)) then
              downflux=dwf*dadw*dwdil*xr3lines_dwlin(1,line)
              upflux=upf*daup*xr3lines_uplin(1,line)
c
              dismul=xr3lines_emilin(2,line)
              tau=tau0
              escape=transferout(dismul,tau)
c
              tlin0=tlin0+downflux*escape
              tlin1=tlin1+upflux
c
              dusttau=dusttau0
              retain=1.d0-transferout(dismul,dusttau)
              dusttrans=dusttrans+downflux*retain
c
              sum1=sum1+downflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+downflux*retain*wid*evplk
c
            endif
          enddo
c
          do line=1,nxrllines
            if ((xrllines_bin(line).eq.inl)) then
              downflux=dwf*dadw*dwdil*xrllines_dwlin(1,line)
              upflux=upf*daup*xrllines_uplin(1,line)
c
              dismul=xrllines_emilin(2,line)
              tau=tau0
              escape=transferout(dismul,tau)
c
              tlin0=tlin0+downflux*escape
              tlin1=tlin1+upflux
c
              dusttau=dusttau0
              retain=1.d0-transferout(dismul,dusttau)
              dusttrans=dusttrans+downflux*retain
c
              sum1=sum1+downflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+downflux*retain*wid*evplk
c
            endif
          enddo
c
C         do n=1,xlines
C           if ((xbin(n).eq.inl)) then
C             lin=n
C             downflux=dwf*dadw*dwdil*dwlin(1,lin)
C             upflux=upf*daup*uplin(1,lin)
Cc
C             dismul=emilin(2,lin)
C             tau=tau0
C             escape=transferout(dismul,tau)
Cc
C             tlin0=tlin0+downflux*escape
C             tlin1=tlin1+upflux
Cc
C             dusttau=dusttau0
C             retain=1.d0-transferout(dismul,dusttau)
C             dusttrans=dusttrans+downflux*retain
Cc
C             sum1=sum1+downflux*(1.d0-escape)*wid*evplk
C             dsum1=dsum1+downflux*retain*wid*evplk
Cc
C           endif
C         enddo
        endif
c
c   Hydrogen
c
c
        if (lmod.ne.'LOCL') then
c
          do series=1,nhseries
            do line=1,nhlines
              if (hbin(line,series).eq.inl) then
c
                downflux=dwf*dadw*dwdil*hyddwlin(1,line,series)
                upflux=upf*daup*hyduplin(1,line,series)
                dismul=hydlin(2,line,series)
                tau=tau0
                escape=transferout(dismul,tau)
                thlin0=thlin0+downflux*escape
                thlin1=thlin1+upflux
c
                dusttau=dusttau0
                retain=1.d0-transferout(dismul,dusttau)
                dusttrans=dusttrans+downflux*retain
c
                sum1=sum1+downflux*(1.d0-escape)*wid*evplk
                dsum1=dsum1+downflux*retain*wid*evplk
c
              endif
            enddo
          enddo
        endif
c
c
c   Helium
c
c
        if (lmod.ne.'LOCL') then
c
          do series=1,nheseries
            do line=1,nhelines
              if (hebin(line,series).eq.inl) then
                downflux=dwf*dadw*dwdil*heldwlin(1,line,series)
                upflux=upf*daup*heluplin(1,line,series)
                if (downflux.lt.epsilon) downflux=0.d0
                if (upflux.lt.epsilon) upflux=0.d0
                dismul=hellin(2,line,series)
                tau=tau0
                escape=transferout(dismul,tau)
                thelin0=thelin0+downflux*escape
                thelin1=thelin1+upflux
                dusttau=dusttau0
                retain=1.d0-transferout(dismul,dusttau)
                dusttrans=dusttrans+downflux*retain
c
                sum1=sum1+downflux*(1.d0-escape)*wid*evplk
                dsum1=dsum1+downflux*retain*wid*evplk
c
              endif
            enddo
          enddo
        endif
c
c   Hydrogenic heavies
c
c
        if (lmod.ne.'LOCL') then
c
          do atom=3,atypes
            if ((xhyddwlin(1,1,1,atom).gt.epsilon).or.(xhyduplin(1,1,1,
     &       atom).gt.epsilon)) then
              do series=1,nxhseries
                do line=1,nxhlines
                  if (xhbin(line,series,atom).eq.inl) then
                    downflux=dwf*dadw*dwdil*xhyddwlin(1,line,series,
     &               atom)
                    upflux=upf*daup*xhyduplin(1,line,series,atom)
                    dismul=xhydlin(2,line,series,atom)
                    tau=tau0
                    escape=transferout(dismul,tau)
                    txhlin0=txhlin0+downflux*escape
                    txhlin1=txhlin1+upflux
                    dusttau=dusttau0
                    retain=1.d0-transferout(dismul,dusttau)
                    dusttrans=dusttrans+downflux*retain
c
                    sum1=sum1+downflux*(1.d0-escape)*wid*evplk
                    dsum1=dsum1+downflux*retain*wid*evplk
c
                  endif
                enddo
              enddo
            endif
          enddo
        endif
c
  110   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now add local diffuse field if on-the-spot is not used
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***ADS LOCAL EMISSIVITY CONTRIBUTION TO THE NEW DIFFUSE FIELD
c     *NB. ONLY WHEN THE ON THE SPOT APPROX. IS NOT USED : JSPOT
c
c     Local field intensity is the fraction of the local field that
c     is  passed onto the next zone , ie (1-exp(-tau))/tau
c
c
        totem=0.0d0
        dusttotem=0.0d0
c
        if ((emf.eq.0.0d0).or.(jspot.eq.'YES')) goto 40
c
        pathlength=(dr*fi)
        energ=den*ev
        if (lmod.eq.'CONT') then
        localflux=energ*emidifcont(inl)*pathlength
        else
        localflux=energ*emidif(inl)*pathlength
        endif
c
        tau=tau0
        escape=localout(1.d0,tau)
        totem=totem+localflux*escape
        dusttau=dusttau0
        retain=1.d0-localout(1.d0,dusttau)
        dusttotem=dusttotem+localflux*retain
c
        sum1=sum1+localflux*(1.d0-escape)*wid*evplk
        dsum1=dsum1+localflux*retain*wid*evplk
c         sum2 = sum2 + localflux*escape*wid*evplk
        tlocalem=tlocalem+localflux*escape*wid*evplk
c
        if (lmod.eq.'CONT') goto 40
c
c   Resonance Lines
c
        do line=1,nxr3lines
          if ((xr3lines_bin(line).eq.inl)) then
            if ((xr3lines_emilin(1,line).gt.epsilon)) then
c
              energ=xr3lines_egij(line)
              localflux=energ*xr3lines_emilin(1,line)*pathlength
              tau=tau0
              dismul=xr3lines_emilin(2,line)
              escape=localout(dismul,tau)
              totem=totem+localflux*escape
              dusttau=dusttau0
              retain=1.d0-localout(dismul,dusttau)
              dusttotem=dusttotem+localflux*retain
c
              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+localflux*retain*wid*evplk
              tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          endif
        enddo
c
        do line=1,nxrllines
          if ((xrllines_bin(line).eq.inl)) then
            if ((xrllines_emilin(1,line).gt.epsilon)) then
c
              energ=xrllines_egij(line)
              localflux=energ*xrllines_emilin(1,line)*pathlength
              tau=tau0
              dismul=xrllines_emilin(2,line)
              escape=localout(dismul,tau)
              totem=totem+localflux*escape
              dusttau=dusttau0
              retain=1.d0-localout(dismul,dusttau)
              dusttotem=dusttotem+localflux*retain
c
              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+localflux*retain*wid*evplk
              tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          endif
        enddo
c
c
C       do m=1,xlines
C         if ((xbin(m).eq.inl).and.(emilin(1,m).gt.epsilon)) then
C           lin=m
C           energ=ev*xejk(lin)
C           localflux=energ*emilin(1,lin)*pathlength
C           tau=tau0
C           dismul=emilin(2,lin)
C           escape=localout(dismul,tau)
C           totem=totem+localflux*escape
C           dusttau=dusttau0
C           retain=1.d0-localout(dismul,dusttau)
C           dusttotem=dusttotem+localflux*retain
Cc
C           sum1=sum1+localflux*(1.d0-escape)*wid*evplk
C           dsum1=dsum1+localflux*retain*wid*evplk
C           tlocalem=tlocalem+localflux*escape*wid*evplk
Cc
C         endif
C       enddo
c
c   Hydrogen
c
        do series=1,nhseries
          do line=1,nhlines
            if ((hbin(line,series).eq.inl).and.(hydlin(1,line,series)
     &       .gt.epsilon)) then
              energ=ev*lmev/hlambda(line,series)
              localflux=energ*hydlin(1,line,series)*pathlength
              tau=tau0
              dismul=hydlin(2,line,series)
              escape=localout(dismul,tau)
              totem=totem+localflux*escape
              dusttau=dusttau0
              retain=1.d0-localout(dismul,dusttau)
              dusttotem=dusttotem+localflux*retain
c
              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+localflux*retain*wid*evplk
c               dsum2 = dsum2 + localflux*escape*wid*evplk
              tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          enddo
        enddo
c
c
c   Helium
c
        do series=1,nheseries
          do line=1,nhelines
            if ((hebin(line,series).eq.inl).and.(hellin(1,line,series)
     &       .gt.epsilon)) then
              dfb=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfb=1.d0-fb
              endif
c
c degrade local HeII 303 photons locally to OIII BF lines by an
c approximate conversion fraction....  Set to 1-0 = 1.0, no
c conversion for now
c
              energ=ev*lmev/helambda(line,series)
              localflux=dfb*energ*hellin(1,line,series)*pathlength
              tau=tau0
              dismul=hellin(2,line,series)
              escape=localout(dismul,tau)
              totem=totem+localflux*escape
              dusttau=dusttau0
              retain=1.d0-localout(dismul,dusttau)
              dusttotem=dusttotem+localflux*retain
c
              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
              dsum1=dsum1+localflux*retain*wid*evplk
c               dsum3 = dsum3 + localflux*escape*wid*evplk
              tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          enddo
        enddo
c
c   Heavy Hydrogen
c
        do atom=3,atypes
          if ((xhydlin(1,1,1,atom).gt.epsilon)) then
            do series=1,nxhseries
              do line=1,nxhlines
                if (xhbin(line,series,atom).eq.inl) then
                  energ=ev*lmev/xhlambda(line,series,atom)
                  localflux=energ*xhydlin(1,line,series,atom)*
     &             pathlength
                  tau=tau0
                  dismul=xhydlin(2,line,series,atom)
                  escape=localout(dismul,tau)
                  totem=totem+localflux*escape
                  dusttau=dusttau0
                  retain=1.d0-localout(dismul,dusttau)
                  dusttotem=dusttotem+localflux*retain
c
                  sum1=sum1+localflux*(1.d0-escape)*wid*evplk
                  dsum1=dsum1+localflux*retain*wid*evplk
c               dsum3 = dsum3 + localflux*escape*wid*evplk
                  tlocalem=tlocalem+localflux*escape*wid*evplk
c
                endif
              enddo
            enddo
          endif
        enddo
c
   40   continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        totem=totem*dwf
        dusttotem=dusttotem*dwf
c
c******************************************************************
c
c     Further attenuate the downstream field by local x-sect.
c     Upstream field is already attenuated by popint in newdif.
c
c     RF0 is the downward flux(???).
c     RF1 is the upward flux(???).
c
        tphot(inl)=rf0+rf1
c     Add metal resonance lines escaped downward flux (TLIN0) and
c     upward flux (TLIN1).
        tphot(inl)=tphot(inl)+tlin0+tlin1
c     Add hydrogen lines escaped downward flux (THLIN0) and
c     upward flux (THLIN1).
        tphot(inl)=tphot(inl)+thlin0+thlin1
c     Add helium lines escaped downward flux (THELIN0) and
c     upward flux (THELIN1).
        tphot(inl)=tphot(inl)+thelin0+thelin1
c     Add hydrogenic heavies lines escaped downward flux (TXHLIN0)
c     and upward flux (THLIN1).
        tphot(inl)=tphot(inl)+txhlin0+txhlin1
c
c     Add local emissivity contribution to the diffuse field (TOTEM).
c
        tphot(inl)=tphot(inl)+totem
c
        sphot(inl)=srcf0*escape0
c
        if (den.gt.iph) then
          totsum=totsum+tphot(inl)*wid
          srcsum=srcsum+sphot(inl)*wid
        endif
c
c     Add component of the downward flux absorbed by dust(???) and ???.
c
        dustphot(inl)=dusttrans+dusttotem
c
c******************************************************************
c
        if (lmod.ne.'LOCL') then
        if (lmod.ne.'CONT') then
        if (lmod(1:3).ne.'NEB') then
c
c     ephot is in eV
c     rhoion = densnum(dh)
c     energ = mean energy of bin in ergs
c
          energ=den*ev
          phots=(tphot(inl)/energ)*evplk*wid
c
          skipbin(inl)=.false.
          bincount=bincount+1
          if (((xsec(inl)/dh)*phots).le.epsilon) then
            skipbin(inl)=.true.
            bincount=bincount-1
          endif
c
          if (photonmode.eq.0) skipbin(inl)=.true.
c
        endif
        endif
        endif
c
   50 continue
c
c
c     ***FORCE RECALCULATION OF PHOTOIONISING RATES
c
c
      if (lmod.ne.'LOCL') then
      if (lmod.ne.'CONT') then
      if (lmod(1:3).ne.'NEB') then
        ipho=ipho+1
      endif
      endif
      endif
c
c*****************************************************************
c
c     Jspec disabled
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine totphot2 (t, dh, fi, rad, dr, dv, wd, lmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Uses rad - rempty for escape probabilites of resonance lines
c
c     T    = Electron temperature in K.
c     DH   = Total Hydrogen density
c     FI   = Filling factor
c     RAD  = Distance from photoionization source in cm
c     DR   = Width of space step in cm
c     DV   = Width of velocity step in cm s^-1
c     LMOD = Code for radiation components to include.
c
c     FORM VECTOR : TPHOT REPRESENTING THE MEAN INTENSITY
c     OF RADIATION : JNU (AT THE POINT CONSIDERED) USING
c     SOURCE VECTOR : SOUPHO(=INU AT RDIST=RSTAR).
c
c     ADDS DIFFUSE FIELD CONTAINED IN VECTORS :
c          UPDIF,DWDIF,UPLIN AND DWLIN
c     ALSO ADDS LOCAL EMISSIVITY VECTORS : EMIDIF,EMILIN
c
c     NOTE : UNITS FOR THE VECTORS EMIDIF & EMILIN ARE
c              IN NUMBERS OF PHOTONS
c              (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c
c   GENERALLY SUBR. LOCALEM  NEEDS TO BE PREVIOUSLY CALLED
c   GENERALLY SUBR. NEWDIF  NEEDS TO BE PREVIOUSLY CALLED
c
c     ABSORPTION : TAU  IS DETERMINED USING THE INTEGRATED
c     COLUMN DENSITY OF MATTER : POPINT  FOR THE SOURCE VECTOR
c     AND THE COLUMN DENS. WITHIN THE SPACE STEP ITSELF FOR ALL
c     VECTORS (FROM R UP TO R+DR  USING POPUL. IN : POP(6,11))
c
c     THE GEOMETRICAL DILUTION FACTOR : WDIL  NEEDS TO BE GIVEN.
c     GENERALLY : WDIL=(RDIST**2)/(4*RSTAR**2)  (SEE FUNCTION FDILU).
c     IN THE PLANE-PARALLEL CASE, WDIL HAS TO BE EQUAL TO 0.5.
c
c     THE DILUTION FACTOR WITHIN THE SPACE STEP ITSELF IS
c     DETERMINED FROM THE RADIUS OF CURVATURE : RAD  DEFINED
c     AT THE INNER LIMIT OF THE SPACE STEP : DR
c     FOR THE UPSTREAM VECTOR, IT CORRESPONDS (BY CONVENTION FOR
c     AN INWARD FLUX) TO A CONCENTRATION FACTOR
c     (GIVE A VALUE FOR RAD OF ZERO FOR PLANE PARALLEL GEOMETRY)
c     CALL SUBR. RESTRANS,CASAB
c
c
c     LMOD='DW'   : ADDS SOURCE VECTOR AND DOWNSTREAM VECTOR.
c     LMOD='UP'   : ADDS SOURCE VECTOR AND UPSTREAM VECTOR.
c     LMOD='ALL'  : ADDS ALL COMPONENTS.
c     LMOD='SO'   : ADDS ONLY SOURCE VECTOR (ON THE SPOT APPROX.)
c                   (IN THIS CASE SUBR. LOCALEM IS UNNECESSARY).
c     LMOD='CAB'  : THE SUBROUTINE ONLY DETERMINE CASE A,B
c                   (NO PHOTON FIELD AT ALL).
c     LMOD='LOCL' : ONLY BUILD LOCAL DIFFUSE FIELD.
c     LMOD='NEBL' : ONLY INTEGRATED NEBULA EMISSION (IE NO SOURCE)
c     LMOD='CONT' : No line emission
c     LMOD='NEBC' : ONLY INTEGRATED NEBULA EMISSION Continuum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t, dh, fi, rad, dr, dv, wd
      character lmod*4
      real*8 dadw,daup,den,dismul
      real*8 drh,dvh,dwf,emf,energ
      real*8 swdil, dr_eff
      real*8 upf,wid
      real*8 tauso, tau, sigmt, curad, wlo, ulo
      real*8 srcf0
      real*8 dwdil,rhoion,phots,casab,fb,dfb
      real*8 pathlength
      real*8 localflux
      real*8 escape,sum
      real*8 rf0,rf1
      real*8 tlin0,tlin1,thlin0,thlin1
      real*8 txhlin0,txhlin1,thelin0,thelin1
c
c      real*8 ti, dustloss
c
      real*8 dusttau, dustsigmat
      real*8 dusttrans
      real*8 downflux, upflux, retain
c
      real*8 tau0_i(mxinfph)
      real*8 tauso0_i(mxinfph)
      real*8 trans0_i(mxinfph)
      real*8 locout0_i(mxinfph)
      real*8 dusttau0_i(mxinfph)
      real*8 dtransretain0_i(mxinfph)
      real*8 dlocretain0_i(mxinfph)
c
c      real*8 sum1,sum2,sum3,tlocalem
c      real*8 dsum1,dsum2,dsum3
c
      integer*4 inl,atom,nz
      integer*4 bincount
      integer*4 line,series
      integer*4 nxlinestart, mxlinestart
c
c           Functions
c
      real*8 densnum, localout, transferout
      real*8 fdilu,fbowen
c
c     constants
c
      rhoion=densnum(dh)
c
c       write(*,*) 'Totphot:',t, dh, fi, rad, dr, dv, wd, lmod
c
      if ((lmod(1:2).ne.'DW').and.(lmod(1:2).ne.'UP')
     &.and.(lmod(1:3).ne.'ALL').and.(lmod(1:2).ne.'SO')
     &.and.(lmod(1:4).ne.'LOCL').and.(lmod(1:3).ne.'NEB')
     &.and.(lmod(1:4).ne.'CONT').and.(lmod(1:3).ne.'CAB')) then
        write (*,10) lmod(1:4)
   10     format(//,'MODE IMPROPERLY DEFINED IN TOTPHOT :',a4)
        stop
      endif
      if ((fi.gt.1.0d0).or.(fi.le.0.0d0)) then
        write (*,20) fi
   20   format(//,'INCONSISTENT FILLING FACTOR IN TOTPHOT :',1pg10.3)
        stop
      endif
c
      if ((wd.gt.0.5d0).or.(wd.lt.0.0d0)) then
        write (*,30) wd
   30   format(/,'GEOMETRICAL DILUTION FACTOR IS IMPROPER IN',
     &     ' SUBR. TOTPHOT  :',1pg10.3)
        stop
      endif
c
c     Compute distance multipliers for resonance lines and
c     determine case A,B.
c
      dr_eff=2.d0*(rad+0.5d0*dr-remp)
c       dr_eff = (rad + 0.5d0*dr - remp)
      drh=0.5d0*dr_eff
      if (jgeo.eq.'S') then
        drh=dmin1((rstromhb-remp),drh)
      endif
      dvh=dv
c
      call emilindismuls (t, dh, drh, dvh)
c
      fb=0.d0
      if (zmap(2).gt.0) then
        fb=fbowen(t,dvh)
c        write(*,*) 'tot fb',fb, 1.d0-fb
      endif
      dfb=1.d0-fb
      heiioiiibf=fb
c
c     call casab with hydrogen Lyman gamma
c
      if (zmap(1).ne.0) then
        line=3
        series=1
        caseab(1)=casab(1,line,series)
      endif
c
c     call casab with Helium Lyman gamma
c
      if (zmap(2).ne.0) then
        line=3
        series=1
        caseab(2)=casab(2,line,series)
      endif
c
c     Assume heavy series are case A
c
      do atom=3,atypes
        line=3
        series=1
        caseab(atom)=casab(mapz(atom),line,series)
      enddo
c
      if (lmod.eq.'CAB') return
c
c     Set mode and internal dilution factors.
c
      upf=1.0d0
      dwf=1.0d0
      emf=1.0d0
c
      if ((lmod.ne.'UP').and.(lmod.ne.'ALL')) upf=0.0d0
      if ((lmod.ne.'DW').and.(lmod.ne.'ALL')) dwf=0.0d0
c
      if (lmod.eq.'LOCL') dwf=1.0d0
      if (lmod(1:3).eq.'NEB') dwf=1.0d0
      if (lmod.eq.'CONT') dwf=1.0d0
c
      if (lmod.eq.'SO') emf=0.0d0
C     if (lmod(1:3).eq.'NEB') emf=0.0d0
C     if (lmod.eq.'CONT') emf=0.0d0
c
      curad=rad
      dadw=1.0d0
      daup=1.0d0
c
      if ((curad.gt.0.d0).and.(jgeo.eq.'S')) then
        wlo=(2.0d0*dlog(curad))-(2.0d0*dlog(curad+dr))
        ulo=(2.0d0*dlog(curad+dr))-(2.0d0*dlog(curad))
        dadw=dexp(wlo)
        daup=dexp(ulo)
      endif
c
      dwdil=1.d0
      if (jgeo.eq.'F') then
        dilradpho=0.d0
        dwdil=fdilu(rshock,rad-dilradpho)
      endif
c
c   Additional source dilution because sources are Inu (ie 1/pi)
c   and wdil = 0.5 in plane parallel (1/2pi) and phion integrates
c   4 pi, so we need an additional 0.5 for sources on one side
c   due to units being used (so phion effectively integrates Inu*0.5
c   sources in pp over only 2pi radians not 4). Local 1/4pi fields
c   are integrated normally and lines have their distance mulipliers
c   and are not strictly one sided.  Handled implicitly in spherical models
c   by fdilu.
c
      swdil=1.d0
      if (jgeo.eq.'P') then
        swdil=0.5d0
      endif
c
c     write(*,*) 'dils',wd,dwdil,swdil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Form vector TPHOT bin by bin.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      bincount=0
c
c      dustloss = 0.0d0
c      ti = 0.0d0
c      sum1=0.d0
c      dsum1=0.d0
c      sum2=0.d0
c      dsum2=0.d0
c      sum3=0.d0
c      dsum3=0.d0
cc
c      tlocalem=0.d0
      sum=0.d0
c
      totsum=0.d0
      srcsum=0.d0
c
      nxlinestart=1
      mxlinestart=1
      pathlength=dr*fi
c
      do inl=1,infph-1
c
        tphot(inl)=0.d0
        sphot(inl)=0.d0
c
        resflux(inl)=0.0d0
        resflux2(inl)=0.0d0        
c    
        tauso=0.d0
        sigmt=0.d0
        dustsigmat=0.d0
        call crosssections (inl, tauso, sigmt, dustsigmat)
        tauso0_i(inl)=tauso
        sigmt=dh*sigmt
        dustsigmat=dh*dustsigmat
c
        tau0_i(inl)=pathlength*sigmt
        trans0_i(inl)=transferout(1.d0,tau0_i(inl))
        locout0_i(inl)=localout(1.d0,tau0_i(inl))
c
        dusttau0_i(inl)=pathlength*dustsigmat
        dtransretain0_i(inl)=1.d0-transferout(1.d0,dusttau0_i(inl))
        dlocretain0_i(inl)=1.d0-localout(1.d0,dusttau0_i(inl))
c     XSEC is use to determine "skip bin" for photoionization
c     calculations only, don't include dust in these.
        xsec(inl)=dmax1(0.d0,sigmt-dustsigmat)*fi
c
      enddo

      do inl=1,infph-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Add diluted source field and integrated diffuse (source) field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Downstream, source attenuated by local opacity.
c
        srcf0=0.d0
        rf0=0.d0
        rf1=0.d0
c
        tlin0=0.d0
        tlin1=0.d0
c
        tauso=tauso0_i(inl)
        tau=tau0_i(inl)
        dusttau=dusttau0_i(inl)
        dusttrans=0.d0
c
        if ((lmod.ne.'LOCL')) then
c
c      Source field if not local and not nebula only
c
          if (lmod(1:3).ne.'NEB') then
            srcf0=wd*swdil*soupho(inl)*dexp(-tauso)
          endif
c
          rf0=0.d0
          rf1=0.d0
c
          if ((lmod.eq.'CONT').or.(lmod.eq.'NEBC')) then
          downflux=dwf*dadw*dwdil*dwdifcont(inl)+srcf0
          upflux=upf*daup*updifcont(inl)
          else
          downflux=dwf*dadw*dwdil*dwdif(inl)+srcf0
          upflux=upf*daup*updif(inl)
          endif
c
          if ((downflux+upflux).gt.epsilon) then
c
            escape=trans0_i(inl)
            rf0=downflux*escape
            rf1=upflux
c
            retain=dtransretain0_i(inl)
            dusttrans=dusttrans+downflux*retain
c          sum1=sum1+(downflux)*(1.d0-escape)*wid*evplk
c          dsum1=dsum1+(downflux)*retain*wid*evplk
          endif

        endif
c
        dustphot(inl)=dustphot(inl)+dusttrans
        tphot(inl)=tphot(inl)+rf0+rf1
        sphot(inl)=sphot(inl)+srcf0*trans0_i(inl)
        resflux(inl)=resflux(inl)+rf0
      enddo
c
c
      if ((lmod.eq.'CONT').or.(lmod.eq.'NEBC'))  goto 100
c
      if (lmod.ne.'LOCL') then
        do line=1,nxr3lines
          inl=xr3lines_bin(line)
          if (inl.gt.0) then
            if ((xr3lines_dwlin(1,line)+xr3lines_uplin(1,line))
     &       .gt.epsilon) then
c
              downflux=dwf*dadw*dwdil*xr3lines_dwlin(1,line)
              upflux=upf*daup*xr3lines_uplin(1,line)
              if ((downflux+upflux).gt.epsilon) then
c
                dismul=xr3lines_emilin(2,line)
                tau=tau0_i(inl)
                escape=transferout(dismul,tau)
c
                tlin0=downflux*escape
                tlin1=upflux
c
                tphot(inl)=tphot(inl)+tlin0+tlin1
c
                dusttau=dusttau0_i(inl)
                retain=1.d0-transferout(dismul,dusttau)
c              dusttrans=dusttrans+downflux*retain
                dustphot(inl)=dustphot(inl)+downflux*retain
c
c              sum1=sum1+downflux*(1.d0-escape)*wid*evplk
c              dsum1=dsum1+downflux*retain*wid*evplk
c
              endif
            endif
          endif
        enddo
C     endif
C     if (lmod.ne.'LOCL') then
        do line=1,nxrllines
          inl=xrllines_bin(line)
          if (inl.gt.0) then
            if ((xrllines_dwlin(1,line)+xrllines_uplin(1,line))
     &       .gt.epsilon) then
c
              downflux=dwf*dadw*dwdil*xrllines_dwlin(1,line)
              upflux=upf*daup*xrllines_uplin(1,line)
              if ((downflux+upflux).gt.epsilon) then
c
                dismul=xrllines_emilin(2,line)
                tau=tau0_i(inl)
                escape=transferout(dismul,tau)
c
                tlin0=downflux*escape
                tlin1=upflux
c
                tphot(inl)=tphot(inl)+tlin0+tlin1
c
                dusttau=dusttau0_i(inl)
                retain=1.d0-transferout(dismul,dusttau)
c              dusttrans=dusttrans+downflux*retain
                dustphot(inl)=dustphot(inl)+downflux*retain
c
c              sum1=sum1+downflux*(1.d0-escape)*wid*evplk
c              dsum1=dsum1+downflux*retain*wid*evplk
c
              endif
            endif
          endif
        enddo
C      endif
C     if (lmod.ne.'LOCL') then
C       do line=1,xlines
C         inl=xbin(line)
C         if (inl.gt.0) then
C           if ((dwlin(1,line)+uplin(1,line)).gt.epsilon) then
Cc
C             downflux=dwf*dadw*dwdil*dwlin(1,line)
C             upflux=upf*daup*uplin(1,line)
C             if ((downflux+upflux).gt.epsilon) then
Cc
C               dismul=emilin(2,line)
C               tau=tau0_i(inl)
C               escape=transferout(dismul,tau)
Cc
C               tlin0=downflux*escape
C               tlin1=upflux
C               tphot(inl)=tphot(inl)+tlin0+tlin1
Cc
C               dusttau=dusttau0_i(inl)
C               retain=1.d0-transferout(dismul,dusttau)
C               dustphot(inl)=dustphot(inl)+downflux*retain
Cc
Cc                sum1=sum1+downflux*(1.d0-escape)*wid*evplk
Cc                dsum1=dsum1+downflux*retain*wid*evplk
Cc
C             endif
C           endif
C         endif
C       enddo
C     endif
c
C      if (lmod.ne.'LOCL') then
c
        do series=1,nhseries
          do line=1,nhlines
            inl=hbin(line,series)
            if (inl.gt.0) then
              if ((hyddwlin(1,line,series)+hyduplin(1,line,series))
     &         .gt.epsilon) then
c
                downflux=dwf*dadw*dwdil*hyddwlin(1,line,series)
                upflux=upf*daup*hyduplin(1,line,series)
                if ((downflux+upflux).gt.epsilon) then
c                dismul   = hydlin(2,line,series)
                  dismul=1.d0
                  escape=trans0_i(inl)
                  thlin0=downflux*escape
                  thlin1=upflux
                  tphot(inl)=tphot(inl)+thlin0+thlin1
c
                  retain=dtransretain0_i(inl)
                  dustphot(inl)=dustphot(inl)+downflux*retain
c
c                sum1=sum1+downflux*(1.d0-escape)*wid*evplk
c                dsum1=dsum1+downflux*retain*wid*evplk
c
                endif
              endif
            endif
          enddo
        enddo
C     endif
Cc
C     if (lmod.ne.'LOCL') then
        do series=1,nheseries
          do line=1,nhelines
            inl=hebin(line,series)
            if (inl.gt.0) then
              if ((heldwlin(1,line,series)+heluplin(1,line,series))
     &         .gt.epsilon) then
c
                downflux=dwf*dadw*dwdil*heldwlin(1,line,series)
                upflux=upf*daup*heluplin(1,line,series)
                if ((downflux+upflux).gt.epsilon) then
c
                  dismul=1.d0
                  escape=trans0_i(inl)
                  thelin0=downflux*escape
                  thelin1=upflux
                  tphot(inl)=tphot(inl)+thelin0+thelin1
                  retain=dtransretain0_i(inl)
                  dustphot(inl)=dustphot(inl)+downflux*retain
c
c                sum1=sum1+downflux*(1.d0-escape)*wid*evplk
c                dsum1=dsum1+downflux*retain*wid*evplk
c
                endif
              endif
            endif
          enddo
        enddo
C
C     endif
C
C     if (lmod.ne.'LOCL') then
        do atom=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              inl=xhbin(line,series,atom)
              if (inl.gt.0) then
                if ((xhyddwlin(1,line,series,atom)+xhyduplin(1,line,
     &           series,atom)).gt.epsilon) then
                  downflux=dwf*dadw*dwdil*xhyddwlin(1,line,series,atom)
                  upflux=upf*daup*xhyduplin(1,line,series,atom)
                  if ((downflux+upflux).gt.epsilon) then
                    dismul=1.d0
                    escape=trans0_i(inl)
                    txhlin0=downflux*escape
                    txhlin1=upflux
                    tphot(inl)=tphot(inl)+txhlin0+txhlin1
                    retain=dtransretain0_i(inl)
                    dustphot(inl)=dustphot(inl)+downflux*retain
c
c                   sum1=sum1+downflux*(1.d0-escape)*wid*evplk
c                   dsum1=dsum1+downflux*retain*wid*evplk
c
                  endif
                endif
              endif
            enddo
          enddo
        enddo
c ! local
      endif
c
c
  100 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now add local diffuse + local lines if on-the-spot is not used
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ((emf.eq.0.0d0).or.(jspot.eq.'YES')) goto 40
c
      if ((lmod.eq.'CONT').or.(lmod.eq.'NEBC'))  goto 110
c
c        Local Resonance Lines
c
      do line=1,nxr3lines
        inl=xr3lines_bin(line)
        if (inl.gt.0) then
          if ((xr3lines_emilin(1,line).gt.epsilon)) then
            energ=xr3lines_egij(line)
            localflux=energ*xr3lines_emilin(1,line)*pathlength
            if (localflux.gt.epsilon) then
              tau=tau0_i(inl)
              dismul=xr3lines_emilin(2,line)
              escape=localout(dismul,tau)
              tphot(inl)=tphot(inl)+localflux*escape
              dusttau=dusttau0_i(inl)
              retain=1.d0-localout(dismul,dusttau)
              dustphot(inl)=dustphot(inl)+localflux*retain
c
c            sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c            dsum1=dsum1+localflux*retain*wid*evplk
c            tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          endif
        endif
      enddo
c
      do line=1,nxrllines
        inl=xrllines_bin(line)
        if (inl.gt.0) then
          if ((xrllines_emilin(1,line).gt.epsilon)) then
            energ=xrllines_egij(line)
            localflux=energ*xrllines_emilin(1,line)*pathlength
            if (localflux.gt.epsilon) then
              tau=tau0_i(inl)
              dismul=xrllines_emilin(2,line)
              escape=localout(dismul,tau)
              tphot(inl)=tphot(inl)+localflux*escape
              dusttau=dusttau0_i(inl)
              retain=1.d0-localout(dismul,dusttau)
              dustphot(inl)=dustphot(inl)+localflux*retain
c
c            sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c            dsum1=dsum1+localflux*retain*wid*evplk
c            tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
          endif
        endif
      enddo
c
C     do line=1,xlines
C       inl=xbin(line)
C       if (inl.gt.0) then
C         if (emilin(1,line).gt.epsilon) then
Cc
C           energ=ev*xejk(line)
C           localflux=energ*emilin(1,line)*pathlength
C           if (localflux.gt.epsilon) then
Cc
C             tau=tau0_i(inl)
C             dismul=emilin(2,line)
C             escape=localout(dismul,tau)
C             tphot(inl)=tphot(inl)+localflux*escape
C             dusttau=dusttau0_i(inl)
C             retain=1.d0-localout(dismul,dusttau)
C             dustphot(inl)=dustphot(inl)+localflux*retain
Cc
Cc              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
Cc              dsum1=dsum1+localflux*retain*wid*evplk
Cc              tlocalem=tlocalem+localflux*escape*wid*evplk
C           endif
C         endif
C       endif
C     enddo
c
      do series=1,nhseries
        do line=1,nhlines
          inl=hbin(line,series)
          if (inl.gt.0) then
            if (hydlin(1,line,series).gt.epsilon) then
              energ=ev*lmev/hlambda(line,series)
              localflux=energ*hydlin(1,line,series)*pathlength
              if (localflux.gt.epsilon) then
c
                dismul=1.d0
                escape=locout0_i(inl)
                tphot(inl)=tphot(inl)+localflux*escape
c
                retain=dlocretain0_i(inl)
                dustphot(inl)=dustphot(inl)+localflux*retain
c
c              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c              dsum1=dsum1+localflux*retain*wid*evplk
c              tlocalem=tlocalem+localflux*escape*wid*evplk
c
              endif
            endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          inl=hebin(line,series)
          if (inl.gt.0) then
            if (hellin(1,line,series).gt.epsilon) then
              energ=ev*lmev/helambda(line,series)
              dfb=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfb=1.d0-fb
              endif
c
c degrade local HeII 303 photons locally to OIII BF lines by an
c approximate conversion fraction....  Set to 1-0 = 1.0, no
c conversion for now
c
              localflux=dfb*energ*hellin(1,line,series)*pathlength
              if (localflux.gt.epsilon) then
                dismul=1.d0
                escape=locout0_i(inl)
                tphot(inl)=tphot(inl)+localflux*escape
                retain=dlocretain0_i(inl)
                dustphot(inl)=dustphot(inl)+localflux*retain
c
c              sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c              dsum1=dsum1+localflux*retain*wid*evplk
c              tlocalem=tlocalem+localflux*escape*wid*evplk
c
              endif
            endif
          endif
        enddo
      enddo
c
      do atom=3,atypes
        nz=mapz(atom)
        do series=1,nxhseries
          do line=1,nxhlines
            inl=xhbin(line,series,atom)
            if (inl.gt.0) then
              if (xhydlin(1,line,series,atom).gt.epsilon) then
                energ=ev*lmev/xhlambda(line,series,atom)
                localflux=energ*xhydlin(1,line,series,atom)*pathlength
                if (localflux.gt.epsilon) then
c
                  dismul=1.d0
                  escape=locout0_i(inl)
                  tphot(inl)=tphot(inl)+localflux*escape
c
                  retain=dlocretain0_i(inl)
                  dustphot(inl)=dustphot(inl)+localflux*retain
c
c                 sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c                 dsum1=dsum1+localflux*retain*wid*evplk
c                 tlocalem=tlocalem+localflux*escape*wid*evplk
c
                endif
              endif
            endif
          enddo
        enddo
      enddo
c
  110 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now add local diffuse field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do inl=1,infph-1

        energ=cphotev(inl)*ev

        if ((lmod.eq.'CONT').or.(lmod.eq.'NEBC')) then
            localflux=energ*emidifcont(inl)*pathlength
        else
            localflux=energ*emidif(inl)*pathlength
        endif

        if (localflux.gt.epsilon) then
c
          escape=locout0_i(inl)
          tphot(inl)=tphot(inl)+localflux*escape
c
          retain=dlocretain0_i(inl)
          dustphot(inl)=dustphot(inl)+localflux*retain
c
c             sum1=sum1+localflux*(1.d0-escape)*wid*evplk
c             dsum1=dsum1+localflux*retain*wid*evplk
c             tlocalem=tlocalem+localflux*escape*wid*evplk
c
        endif
      enddo
   40 continue
c
      do inl=1,infph-1
        if (cphotev(inl).gt.iph) then
          wid=(photev(inl+1)-photev(inl))
          totsum=totsum+tphot(inl)*wid
          srcsum=srcsum+sphot(inl)*wid
        endif
      enddo
c
c      write(*,*) totsum, srcsum
c
      if (lmod.ne.'LOCL') then
      if (lmod.ne.'CONT') then
      if (lmod(1:3).ne.'NEB') then
        do inl=1,infph-1
          wid=(photev(inl+1)-photev(inl))
          den=cphotev(inl)
          energ=den*ev
          phots=(tphot(inl)/energ)*evplk*wid
          skipbin(inl)=.false.
          bincount=bincount+1
          if (((xsec(inl)/dh)*phots).le.epsilon) then
            skipbin(inl)=.true.
            bincount=bincount-1
          endif
          if (photonmode.eq.0) skipbin(inl)=.true.
        enddo
        ipho=ipho+1
      endif
      endif
      endif
c
c*****************************************************************
c
c     Jspec disabled
c
      return
      end
