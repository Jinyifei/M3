cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newdif (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     &frdw, jmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO COMPUTE NEW DIFFUSE FIELD VECTORS AT THE END OF SPACE
c     STEP,USING LOCAL EMISSIVITY PRODUCED BY SUBR. LOCALEM)
c
c     DOWNSTREAM AND UPSTREAM TEMPERATURES : TDW,TUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF DISTANCE :DRDW,DRUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF VELOCITY : DVDW,DVUP
c     REM: FOR THE NON LOCAL COMPONENT,THESE QUANTITIES REFER
c     TO THE TOTAL VARIATION OF T,DR,DV THAT TAKES
c     PLACE BETWEEN THE ORIGIN AND THE POINT CONSIDERED
c     FRDW : FRACTION OF LOCAL EMISSIVITY THAT IS ADDED TO THE
c     DOWNSTREAM COMPONENT;THE COMPLEMENT IS ADDED TO
c     THE UPSTREAM VECTORS
c
c     THE DILUTION FACTOR IN CASE OF SPHERICAL SYMMETRY IS
c     DERIVED FROM THE CURVATURE RADIUS AT THE INNER LIMIT
c     OF THE SPACESTEP(=LOCAL DR)  FOR PLANE PARALLEL SYMETRY
c     THE RAD (DISTANCE) IS IGNORED, FOR CYLINDERS IT IS USED
c     TO FORM WEIGHTED SUMS FOR DIFFUSE FIELD DILUTIONS
c     BY CONVENTION,THE DOWNSTREAM VECTORS ARE RESERVED FOR
c     THE OUTWARD (DIRECTION OF INCREASING R) FLUX AND THE
c     THE UPSTREAM ONES FOR THE INWARD FLUX
c
c     NOTE : UNITS FOR THE VECTORS EMIDIF,EMILIN ARE IN NUMBERS
c     OF PHOTONS (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c
c
c     CALL SUBR.RESTRANS
c
c     JMOD='LODW' MEANS THAT ONLY THE DOWNSTREAM DIFFUSE FIELD
c     IS LOCAL.THE FRACTION ADDED TO THE UPSTREAM
c     COMPONENT IS FURTHER ATTENUATED USING
c     THE INTEGRATED COLUMN DENSITIES : POPINT
c     JMOD='LOUP' REVERSED SITUATION TO THE PREVIOUS CASE
c     JMOD='OUTW' OUTWARD ONLY APPROXIMATION WHERE ALL THE
c     LOCAL EMISSIVITY IS ADDED TO THE DOWNSTREAM
c     COMPONENT;NO UPSTREAM DIFFUSE FIELD
c     JMOD='DWUP' BOTH STREAMS ARE LOCAL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 dildw
      real*8 dilup,dismul,dr,drem,dremh,dvem,dvemh,dwem
      real*8 dwex, dwf, dwliem, energ, telc
      real*8 upem, upex, upf, upliem
      real*8 wadw, waup, wedw, weup
      real*8 tauso, tau, tau0
      real*8 sigmt, curad, rm, wlo, ulo
      real*8 tdw, tup, dh, fi, rad, drdw, dvdw, drup
      real*8 dvup,frdw,emis,z,p,f,es
      real*8 dfbu,dfbd,dfb
      real*8 fbu,fbd,fb
      real*8 drdwh,dvdwh,druph,dvuph
      real*8 dustsigmat
      real*8 pathlength
      real*8 escape, escapeso
c
      integer*4 line,series,j, atom,ion,inl
c
      character jmod*4
c
c           Functions
c
      real*8 fdismul,localout, tauline, fbowen
c
c
c     Check parameters and replace with sane values if necessary
c
c
c filling factor
c
      if ((fi.gt.1.0d0).and.(fi.le.0.0d0)) then
        write (*,10) fi
   10     format(//' INCONSISTENT FILLING FACTOR IN NEWDIF: ',1pg10.3)
        fi=1.0
      endif
      dvdw=dabs(dvdw)
      dvup=dabs(dvup)
c
      if ((tdw.lt.0.d0).or.(tup.lt.0.d0).or.(rad.lt.0.d0)
     &.or.(drdw.lt.0.d0).or.(drup.lt.0.d0)) then
        write (*,20) tdw,tup,rad,drdw,drup
   20   format(/,/,'ONE OF THE ARGUMENTS HAS A NEGATIVE VALUE IN NEWDIF'
     &         ,'****',/,5(1pg10.3))
        tdw=dabs(tdw)
        tup=dabs(tup)
        rad=dabs(rad)
        drdw=dabs(drdw)
        drup=dabs(drup)
      endif
c
      if ((jmod.eq.'LODW').or.(jmod.eq.'LOUP').or.(jmod.eq.'OUTW')
     &.or.(jmod.eq.'DWUP')) goto 40
      write (*,30) jmod
   30 format(//' MODE IMPROPERLY DEFINED IN NEWDIF : ',a4)
      stop
   40 continue
c
      if ((frdw.gt.1.d0).or.(frdw.lt.0.d0)) then
        write (*,50) frdw
   50       format(/,/,' INCONSISTENT VALUE FOR THE ARGUMENT FRDW :'
     &             ,1pg10.3)
        frdw=0.5d0
      endif
c
c     ***SET MODE AND INTERNAL DILUTION FACTORS
c
      dwf=frdw
      if (jmod.eq.'OUTW') dwf=1.0d0
      upf=1.d0-dwf
c
      dwex=1.d0
      if (jmod.eq.'LOUP') dwex=0.d0
      upex=1.d0-dwex
      if (jmod.eq.'DWUP') upex=1.d0
c
      if (jmod.eq.'LOUP') then
        dr=drup
        drem=drup
        dvem=dvup
        telc=tup
      else
        dr=drdw
        drem=drdw
        dvem=dvdw
        telc=tdw
      endif
c
c     Plane Parallel / Finite Cylinder Default
c
      curad=rad
      dildw=1.d0
      dilup=1.d0
      if (jgeo.eq.'S') then
c
c     Spherical
c
        rm=curad+drdw
        wlo=(2.d0*dlog(curad))-(2.d0*dlog(rm))
        rm=(curad+dr)-drup
        if (rm.le.0.d0) then
          write (*,60) drup,curad
   60       format(/,/,' DRUP IS LARGER THAN THE RADIUS OF CURVATURE :'
     &              ,2(1pg10.3))
          stop
        endif
        ulo=(2.d0*dlog(curad+dr))-(2.d0*dlog(rm))
        dildw=dexp(wlo)
        dilup=dexp(ulo)
c
c     Spherical
c
      endif
c
      wadw=dildw
      wedw=1.d0
      waup=1.d0
      weup=dilup
      if (jmod.ne.'LOUP') goto 70
      wadw=1.d0
      wedw=dildw
      waup=dilup
      weup=1.d0
   70 if (jmod.ne.'DWUP') goto 80
      wadw=dildw
      wedw=1.d0
      waup=dilup
      weup=1.d0
   80 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Begin Routine Proper
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***COMPUTES RESONANCE-LINES DISTANCE MULTIPLIERS
c
      dremh=drem*0.5d0
      dvemh=dvem*0.5d0
c
      do line=1,nxr3lines
        xr3lines_uplin(2,line)=1.d0
        xr3lines_dwlin(2,line)=1.d0
        xr3lines_emilin(2,line)=1.d0
        if (xr3lines_emilin(1,line).gt.epsilon) then
          if (xr3lines_frac(line).gt.epsilon) then
            atom=xr3lines_at(line)
            ion=xr3lines_ion(line)
            z=zion(atom)
            p=pop(ion,atom)
            if ((z*p.ge.pzlimit)) then
              f=xr3lines_gf(line)
              es=xr3lines_egij(line)
              xr3lines_emilin(2,line)=fdismul(telc,dh,dremh,dvemh,atom,
     &         ion,es,f)
            endif
          endif
        endif
      enddo
c
      do line=1,nxrllines
        xrllines_uplin(2,line)=1.d0
        xrllines_dwlin(2,line)=1.d0
        xrllines_emilin(2,line)=1.d0
        if (xrllines_emilin(1,line).gt.epsilon) then
          if (xrllines_frac(line).gt.epsilon) then
            atom=xrllines_at(line)
            ion=xrllines_ion(line)
            z=zion(atom)
            p=pop(ion,atom)
            if ((z*p.ge.pzlimit)) then
              f=xrllines_gf(line)
              es=xrllines_egij(line)
              xrllines_emilin(2,line)=fdismul(telc,dh,dremh,dvemh,atom,
     &         ion,es,f)
            endif
          endif
        endif
      enddo
c
C     do line=1,xlines
C       uplin(2,line)=1.d0
C       dwlin(2,line)=1.d0
C       emilin(2,line)=1.d0
C       if (xbr(line).gt.epsilon) then
C         z=zion(xrat(line))
C         p=pop(xion(line),xrat(line))
C         if ((z*p.ge.pzlimit)) then
C           f=xfef(line)/xbr(line)
C           es=xejk(line)*ev
C           emilin(2,line)=fdismul(telc,dh,dremh,dvemh,xrat(line),
C    &       xion(line),es,f)
C         endif
C       endif
C     enddo
c
c   Hydrogenic lines
c
      fb=0.0d0
      fbu=0.0d0
      fbd=0.0d0
      if (zmap(2).gt.0) then
        fb=fbowen(telc,dvemh)
        fbu=fbowen(tup,dvup)
        fbd=fbowen(tdw,dvdw)
      endif
      heiioiiibf=fb
c
      atom=1
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nhseries
        do line=1,nhlines
          f=hydrogf(line,series)
          hyduplin(2,line,series)=1.d0
          hyddwlin(2,line,series)=1.d0
          hydlin(2,line,series)=1.d0
          j=hbin(line,series)
          if (j.ne.0) then
            es=lmev/hlambda(line,series)*ev
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              hydlin(2,line,series)=fdismul(telc,dh,dremh,dvemh,atom,
     &         ion,es,f)
            endif
          endif
        enddo
      enddo
      atom=2
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nheseries
        do line=1,nhelines
          f=hydrogf(line,series)
          heluplin(2,line,series)=1.d0
          heldwlin(2,line,series)=1.d0
          hellin(2,line,series)=1.d0
          j=hebin(line,series)
          if (j.ne.0) then
            es=lmev/helambda(line,series)*ev
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              hellin(2,line,series)=fdismul(telc,dh,dremh,dvemh,atom,
     &         ion,es,f)
            endif
          endif
        enddo
      enddo
c
      do atom=3,atypes
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
        do series=1,nxhseries
          do line=1,nxhlines
            f=hydrogf(line,series)
            xhyduplin(2,line,series,atom)=1.d0
            xhyddwlin(2,line,series,atom)=1.d0
            xhydlin(2,line,series,atom)=1.d0
            j=xhbin(line,series,atom)
            if (j.ne.0) then
              es=lmev/xhlambda(line,series,atom)*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                xhydlin(2,line,series,atom)=fdismul(telc,dh,dremh,dvemh,
     &           atom,ion,es,f)
              endif
            endif
          enddo
        enddo
      enddo
c
c     do up and down if needed
c
      if (dwf.gt.0.d0) then
c
        drdwh=drdw*0.5d0
        dvdwh=dvdw*0.5d0
c
c first the general resonance lines
c
        do line=1,nxr3lines
          xr3lines_dwlin(2,line)=1.d0
          if (xr3lines_dwlin(1,line).gt.epsilon) then
            if (xr3lines_frac(line).gt.epsilon) then
              atom=xr3lines_at(line)
              ion=xr3lines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xr3lines_gf(line)
                es=xr3lines_egij(line)
                xr3lines_dwlin(2,line)=fdismul(tdw,dh,drdwh,dvdwh,atom,
     &           ion,es,f)
              endif
            endif
          endif
        enddo
c
        do line=1,nxrllines
          xrllines_dwlin(2,line)=1.d0
          if (xrllines_dwlin(1,line).gt.epsilon) then
            if (xrllines_frac(line).gt.epsilon) then
              atom=xrllines_at(line)
              ion=xrllines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xrllines_gf(line)
                es=xrllines_egij(line)
                xrllines_dwlin(2,line)=fdismul(tdw,dh,drdwh,dvdwh,atom,
     &           ion,es,f)
              endif
            endif
          endif
        enddo
c
C       do line=1,xlines
C         dwlin(2,line)=1.d0
C         if (xbr(line).gt.0.0) then
C           z=zion(xrat(line))
C           p=pop(xion(line),xrat(line))
C           if ((z*p.ge.pzlimit)) then
C             f=xfef(line)/xbr(line)
C             es=xejk(line)*ev
C             dwlin(2,line)=fdismul(tdw,dh,drdwh,dvdwh,xrat(line),
C    &         xion(line),es,f)
C           endif
C         endif
C       enddo
c
c
c   Hydrogenic lines
c
        atom=1
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
c
        do series=1,nhseries
          do line=1,nhlines
            f=hydrogf(line,series)
            hyddwlin(2,line,series)=1.d0
            j=hbin(line,series)
            if (j.ne.0) then
              es=(lmev/hlambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                hyddwlin(2,line,series)=fdismul(tdw,dh,drdwh,dvdwh,atom,
     &           ion,es,f)
              endif
            endif
          enddo
        enddo
c
        atom=2
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
c
        do series=1,nheseries
          do line=1,nhelines
            f=hydrogf(line,series)
            heldwlin(2,line,series)=1.d0
            j=hebin(line,series)
            if (j.ne.0) then
              es=(lmev/helambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                heldwlin(2,line,series)=fdismul(tdw,dh,drdwh,dvdwh,atom,
     &           ion,es,f)
              endif
            endif
          enddo
        enddo
c
        do atom=3,atypes
          ion=mapz(atom)
          z=zion(atom)
          p=pop(ion,atom)
          do series=1,nxhseries
            do line=1,nxhlines
              f=hydrogf(line,series)
              xhyddwlin(2,line,series,atom)=1.d0
              j=xhbin(line,series,atom)
              if (j.ne.0) then
                es=lmev/xhlambda(line,series,atom)*ev
                if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                  xhyddwlin(2,line,series,atom)=fdismul(tdw,dh,drdwh,
     &             dvdwh,atom,ion,es,f)
                endif
              endif
            enddo
          enddo
        enddo
c dwf
      endif
c
      if (upf.gt.0.d0) then
c
        druph=drup*0.5d0
        dvuph=dvup*0.5d0
c
c first the general resonance lines
c
        do line=1,nxr3lines
          xr3lines_uplin(2,line)=1.d0
          if (xr3lines_uplin(1,line).gt.epsilon) then
            if (xr3lines_frac(line).gt.epsilon) then
              atom=xr3lines_at(line)
              ion=xr3lines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xr3lines_gf(line)
                es=xr3lines_egij(line)
                xr3lines_uplin(2,line)=fdismul(tup,dh,druph,dvuph,atom,
     &           ion,es,f)
              endif
            endif
          endif
        enddo
c
        do line=1,nxrllines
          xrllines_uplin(2,line)=1.d0
          if (xrllines_uplin(1,line).gt.epsilon) then
            if (xrllines_frac(line).gt.epsilon) then
              atom=xrllines_at(line)
              ion=xrllines_ion(line)
              z=zion(atom)
              p=pop(ion,atom)
              if ((z*p.ge.pzlimit)) then
                f=xrllines_gf(line)
                es=xrllines_egij(line)
                xrllines_uplin(2,line)=fdismul(tup,dh,druph,dvuph,atom,
     &           ion,es,f)
              endif
            endif
          endif
        enddo
c
C       do line=1,xlines
C         uplin(2,line)=1.d0
C         if (uplin(1,line).gt.0.0) then
C           z=zion(xrat(line))
C           p=pop(xion(line),xrat(line))
C           if ((z*p.ge.pzlimit)) then
C             f=xfef(line)/xbr(line)
C             es=xejk(line)*ev
C             uplin(2,line)=fdismul(tup,dh,druph,dvuph,xrat(line),
C    &         xion(line),es,f)
C           endif
C         endif
C       enddo
c
c   Hydrogenic lines
c
        atom=1
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
c
        do series=1,nhseries
          do line=1,nhlines
            f=hydrogf(line,series)
            hyduplin(2,line,series)=1.d0
            j=hbin(line,series)
            if (j.ne.0) then
              es=(lmev/hlambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                hyduplin(2,line,series)=fdismul(tup,dh,druph,dvuph,atom,
     &           ion,es,f)
              endif
            endif
          enddo
        enddo
c
        atom=2
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
c
        do series=1,nheseries
          do line=1,nhelines
            f=hydrogf(line,series)
            heluplin(2,line,series)=1.d0
            j=hebin(line,series)
            if (j.ne.0) then
              es=(lmev/helambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                heluplin(2,line,series)=fdismul(tup,dh,druph,dvuph,atom,
     &           ion,es,f)
              endif
            endif
          enddo
        enddo
c
        do atom=3,atypes
          ion=mapz(atom)
          z=zion(atom)
          p=pop(ion,atom)
          do series=1,nxhseries
            do line=1,nxhlines
              f=hydrogf(line,series)
              xhyduplin(2,line,series,atom)=1.d0
              j=xhbin(line,series,atom)
              if (j.ne.0) then
                es=lmev/xhlambda(line,series,atom)*ev
                if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                  xhyduplin(2,line,series,atom)=fdismul(tup,dh,druph,
     &             dvuph,atom,ion,es,f)
                endif
              endif
            enddo
          enddo
        enddo
c upf
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***UPDATES DIFFUSE FIELD VECTORS UPDIF, DWDIF etc BIN BY BIN
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 90 inl=1,infph-1
c
        energ=cphote(inl)
        pathlength=dr*fi
c
        tauso=0.d0
        sigmt=0.d0
        dustsigmat=0.d0
c
        call crosssections (inl, tauso, sigmt, dustsigmat)
c
        sigmt=dh*sigmt
        dustsigmat=dh*dustsigmat
c
c tauso already includes dh, fi and distance implicitly
c
c
c     ***ATTENUATION OF THE ORIGINAL INTENSITY OF THE DIFFUSE FIELD
c
        tau0=pathlength*sigmt
        tau=tau0
c
c tauso already includes dh, fi and distance implicitly
c
c upsteam escape for total optical depth
c
        escapeso=dexp(-tauso)
c
        dwdif(inl)=dwdif(inl)*(wadw*dexp(-(dwex*tau)))
        updif(inl)=updif(inl)*(waup*dexp(-(upex*tau)))
c
        dwdifcont(inl)=dwdifcont(inl)*(wadw*dexp(-(dwex*tau)))
        dwdifcont(inl)=dwdifcont(inl)*(waup*dexp(-(upex*tau)))
c
c     distance multipliers attenuate lines for emilin etc
c
        do line=1,nxr3lines
          if (xr3lines_bin(line).eq.inl) then
            if ((xr3lines_uplin(1,line).gt.epsilon)
     &       .or.(xr3lines_dwlin(1,line).gt.epsilon)) then
              dismul=(dwex*xr3lines_dwlin(2,line))+(upex*
     &         xr3lines_uplin(2,line))
              tau=tauline(dismul,tau0)
              xr3lines_dwlin(1,line)=xr3lines_dwlin(1,line)*(wadw*dexp(-
     &         (dwex*tau)))
              xr3lines_uplin(1,line)=xr3lines_uplin(1,line)*(waup*dexp(-
     &         (upex*tau)))
            endif
          endif
        enddo
c
        do line=1,nxrllines
          if (xrllines_bin(line).eq.inl) then
            if ((xrllines_uplin(1,line).gt.epsilon)
     &       .or.(xrllines_dwlin(1,line).gt.epsilon)) then
              dismul=(dwex*xrllines_dwlin(2,line))+(upex*
     &         xrllines_uplin(2,line))
              tau=tauline(dismul,tau0)
              xrllines_dwlin(1,line)=xrllines_dwlin(1,line)*(wadw*dexp(-
     &         (dwex*tau)))
              xrllines_uplin(1,line)=xrllines_uplin(1,line)*(waup*dexp(-
     &         (upex*tau)))
            endif
          endif
        enddo
c
C       do n=1,xlines
C         if (xbin(n).eq.inl) then
C           if ((uplin(1,n).gt.epsilon).or.(dwlin(1,n).gt.epsilon))
C    &       then
C             dismul=(dwex*dwlin(2,n))+(upex*uplin(2,n))
C             tau=tauline(dismul,tau0)
C             dwlin(1,n)=dwlin(1,n)*(wadw*dexp(-(dwex*tau)))
C             uplin(1,n)=uplin(1,n)*(waup*dexp(-(upex*tau)))
C           endif
C         endif
C       enddo
c
c   Hydrogenic lines
c
        do series=1,nhseries
          do line=1,nhlines
            if (hbin(line,series).eq.inl) then
              dismul=(dwex*hyddwlin(2,line,series))+(upex*hyduplin(2,
     &         line,series))
              tau=tauline(dismul,tau0)
              hyddwlin(1,line,series)=hyddwlin(1,line,series)*(wadw*
     &         dexp(-(dwex*tau)))
              hyduplin(1,line,series)=hyduplin(1,line,series)*(waup*
     &         dexp(-(upex*tau)))
            endif
          enddo
        enddo
c
c   Helium lines
c
        do series=1,nheseries
          do line=1,nhelines
            if (hebin(line,series).eq.inl) then
              dismul=(dwex*heldwlin(2,line,series))+(upex*heluplin(2,
     &         line,series))
              tau=tauline(dismul,tau0)
c
              dfbu=1.d0
              dfbd=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfbu=1.d0-fbu
                dfbd=1.d0-fbd
c                 write(*,*)'fbu fbd',fbu, fbd
              endif
c
              heldwlin(1,line,series)=dfbd*heldwlin(1,line,series)*
     &         (wadw*dexp(-(dwex*tau)))
              heluplin(1,line,series)=dfbu*heluplin(1,line,series)*
     &         (waup*dexp(-(upex*tau)))
            endif
          enddo
        enddo
c
        do atom=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              if (xhbin(line,series,atom).eq.inl) then
                dismul=(dwex*xhyddwlin(2,line,series,atom))+(upex*
     &           xhyduplin(2,line,series,atom))
                tau=tauline(dismul,tau0)
                xhyddwlin(1,line,series,atom)=xhyddwlin(1,line,series,
     &           atom)*(wadw*dexp(-(dwex*tau)))
                xhyduplin(1,line,series,atom)=xhyduplin(1,line,series,
     &           atom)*(waup*dexp(-(upex*tau)))
              endif
            enddo
          enddo
        enddo
c
c     ***ADS LOCAL EMISSIVITY CONTRIBUTION TO THE NEW DIFFUSE FIELD
c
        tau=tau0
        escape=localout(dismul,tau)
c
        emis=energ*pathlength
c
        dwem=(dwf*emis)*escape
        upem=(upf*emis)*escapeso
c
c     Add local diffuse field through this zone
c
        dwdif(inl)=dwdif(inl)+(dwem*wedw*emidif(inl))
        updif(inl)=updif(inl)+(upem*weup*emidif(inl))
c
        dwdifcont(inl)=dwdif(inl)+(dwem*wedw*emidifcont(inl))
        updifcont(inl)=updif(inl)+(upem*weup*emidifcont(inl))
c
c     need to add dust scattering contribution,
c
c
c         if (grainmode.eq.1) then
c
c     not so good, need to assume that old tphot is up to date
c
c       dwsc=0.d0
c       upsc=0.d0
c      if ((ClinPAH.eq.1).AND.(pahactive.eq.0)) then
c        do dtype=2,numtypes
c            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
c            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
c           enddo
c         else
c        do dtype=1,numtypes
c            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
c            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
c           enddo
c         endif
c
c     Add PAH scattering
c
c          if (pahactive) then
c           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
c     &           +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
c     &           *0.5d0*(1+pahncos(inl))
c           dwsc = dwsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
c           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
c     &              +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
c     &              *0.5d0*(1-pahncos(inl))
c           upsc = upsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
c          endif
c
c
c     Add dust scattered source light to diffuse fields
c
c            dwdif(inl) = dwdif(inl)+(dwsc*(wedw*(dexp(-(upex*tauso)))))
c            updif(inl) = updif(inl)+(upsc*(weup*(dexp(-(dwex*tauso)))))
c
c         endif
c
        do line=1,nxr3lines
          if (xr3lines_bin(line).eq.inl) then
            if ((xr3lines_uplin(1,line).gt.epsilon)
     &       .or.(xr3lines_dwlin(1,line).gt.epsilon)) then
              dismul=(dwex*xr3lines_dwlin(2,line))+(upex*
     &         xr3lines_uplin(2,line))
              tau=tau0
              energ=xr3lines_egij(line)
              emis=energ*xr3lines_emilin(1,line)*pathlength
              escape=localout(dismul,tau)
              dwliem=(dwf*emis)*escape
              upliem=(upf*emis)*escapeso
              xr3lines_dwlin(1,line)=xr3lines_dwlin(1,line)+(dwliem*
     &         wedw)
              xr3lines_uplin(1,line)=xr3lines_uplin(1,line)+(upliem*
     &         weup)
            endif
          endif
        enddo
c
        do line=1,nxrllines
          if (xrllines_bin(line).eq.inl) then
            if ((xrllines_uplin(1,line).gt.epsilon)
     &       .or.(xrllines_dwlin(1,line).gt.epsilon)) then
              dismul=(dwex*xrllines_dwlin(2,line))+(upex*
     &         xrllines_uplin(2,line))
              tau=tau0
              energ=xrllines_egij(line)
              emis=energ*xrllines_emilin(1,line)*pathlength
              escape=localout(dismul,tau)
              dwliem=(dwf*emis)*escape
              upliem=(upf*emis)*escapeso
              xrllines_dwlin(1,line)=xrllines_dwlin(1,line)+(dwliem*
     &         wedw)
              xrllines_uplin(1,line)=xrllines_uplin(1,line)+(upliem*
     &         weup)
            endif
          endif
        enddo
c
c     do emilin resonance lines
c
C       do m=1,xlines
C         if (xbin(m).eq.inl) then
C           lin=m
C           dismul=(dwex*dwlin(2,lin))+(upex*uplin(2,lin))
C           tau=tau0
C           energ=ev*xejk(lin)
C           emis=energ*emilin(1,lin)*pathlength
C           escape=localout(dismul,tau)
C           dwliem=(dwf*emis)*escape
C           upliem=(upf*emis)*escapeso
C           dwlin(1,lin)=dwlin(1,lin)+(dwliem*wedw)
C           uplin(1,lin)=uplin(1,lin)+(upliem*weup)
C         endif
C       enddo
c
c Hydrogen
c
        do series=1,nhseries
          do line=1,nhlines
c
            if (hbin(line,series).eq.inl) then
c
              dismul=(dwex*hyddwlin(2,line,series))+(upex*hyduplin(2,
     &         line,series))
c
              tau=tau0
c
              energ=ev*lmev/hlambda(line,series)
              emis=hydlin(1,line,series)*pathlength*energ
c
              escape=localout(dismul,tau)
              dwliem=(dwf*emis)*escape
c
              upliem=(upf*emis)*escapeso
c
              hyddwlin(1,line,series)=hyddwlin(1,line,series)+(dwliem*
     &         wedw)
              hyduplin(1,line,series)=hyduplin(1,line,series)+(upliem*
     &         weup)
c
            endif
c
          enddo
        enddo
c
c Helium
c
        do series=1,nheseries
          do line=1,nhelines
c
            if (hebin(line,series).eq.inl) then
              dismul=(dwex*heldwlin(2,line,series))+(upex*heluplin(2,
     &         line,series))
c
              tau=tau0
              escape=localout(dismul,tau)
              energ=ev*lmev/helambda(line,series)
c
              dfb=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfb=1.d0-fb
              endif
c
c degrade local HeII 303 photons locally to OIII BF lines by an
c approximate conversion fraction....  Set to 1-0 = 1.0, no
c conversion for now
c
              emis=dfb*energ*hellin(1,line,series)*pathlength
              dwliem=(dwf*emis)*escape
              upliem=(upf*emis)*escapeso
              heldwlin(1,line,series)=heldwlin(1,line,series)+(dwliem*
     &         wedw)
              heluplin(1,line,series)=heluplin(1,line,series)+(upliem*
     &         weup)
c
            endif
c
          enddo
        enddo
c
c
c Heavy Hydrogenic
c
        do atom=3,atypes
          if (xhydlin(1,1,1,atom).gt.epsilon) then
            do series=1,nxhseries
              do line=1,nxhlines
c
                if (xhbin(line,series,atom).eq.inl) then
                  dismul=(dwex*xhyddwlin(2,line,series,atom))+(upex*
     &             xhyduplin(2,line,series,atom))
                  tau=tau0
c
                  energ=ev*lmev/xhlambda(line,series,atom)
                  emis=energ*xhydlin(1,line,series,atom)*pathlength
c
                  escape=localout(dismul,tau)
                  dwliem=(dwf*emis)*escape
c
                  upliem=(upf*emis)*escapeso
c
                  xhyddwlin(1,line,series,atom)=xhyddwlin(1,line,series,
     &             atom)+(dwliem*wedw)
                  xhyduplin(1,line,series,atom)=xhyduplin(1,line,series,
     &             atom)+(upliem*weup)
c
                endif
c
              enddo
            enddo
          endif
        enddo
c
c end vector loop inl
c
   90 continue
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
      subroutine newdif2 (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     & frdw, jmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     uses rad - r-empty instead of drdw and dvdw
c
c*******TO COMPUTE NEW DIFFUSE FIELD VECTORS AT THE END OF SPACE
c     STEP,USING LOCAL EMISSIVITY PRODUCED BY SUBR. LOCALEM)
c
c     DOWNSTREAM AND UPSTREAM TEMPERATURES : TDW,TUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF DISTANCE :DRDW,DRUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF VELOCITY : DVDW,DVUP
c     REM: FOR THE NON LOCAL COMPONENT,THESE QUANTITIES REFER
c     TO THE TOTAL VARIATION OF T,DR,DV THAT TAKES
c     PLACE BETWEEN THE ORIGIN AND THE POINT CONSIDERED
c     FRDW : FRACTION OF LOCAL EMISSIVITY THAT IS ADDED TO THE
c     DOWNSTREAM COMPONENT;THE COMPLEMENT IS ADDED TO
c     THE UPSTREAM VECTORS
c
c     THE DILUTION FACTOR IN CASE OF SPHERICAL SYMMETRY IS
c     DERIVED FROM THE CURVATURE RADIUS AT THE INNER LIMIT
c     OF THE SPACESTEP(=LOCAL DR)  FOR PLANE PARALLEL SYMETRY
c     THE RAD (DISTANCE) IS IGNORED, FOR CYLINDERS IT IS USED
c     TO FORM WEIGHTED SUMS FOR DIFFUSE FIELD DILUTIONS
c     BY CONVENTION,THE DOWNSTREAM VECTORS ARE RESERVED FOR
c     THE OUTWARD (DIRECTION OF INCREASING R) FLUX AND THE
c     THE UPSTREAM ONES FOR THE INWARD FLUX
c
c     NOTE : UNITS FOR THE VECTORS EMIDIF,EMILIN ARE IN NUMBERS
c     OF PHOTONS (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c
c
c     CALL SUBR.RESTRANS
c
c     JMOD='LODW' MEANS THAT ONLY THE DOWNSTREAM DIFFUSE FIELD
c     IS LOCAL.THE FRACTION ADDED TO THE UPSTREAM
c     COMPONENT IS FURTHER ATTENUATED USING
c     THE INTEGRATED COLUMN DENSITIES : POPINT
c     JMOD='LOUP' REVERSED SITUATION TO THE PREVIOUS CASE
c     JMOD='OUTW' OUTWARD ONLY APPROXIMATION WHERE ALL THE
c     LOCAL EMISSIVITY IS ADDED TO THE DOWNSTREAM
c     COMPONENT;NO UPSTREAM DIFFUSE FIELD
c     JMOD='DWUP' BOTH STREAMS ARE LOCAL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 dildw,dresc,f1,f2
      real*8 dilup,dismul,dr,drem,dremh,dvem,dvemh,dwem
      real*8 dwex, dwf, dwliem, energ, telc
      real*8 upem, upex, upf, upliem
      real*8 wadw, waup, wedw, weup
      real*8 tauso, tau, tau0
      real*8 sigmt, curad, rm, wlo, ulo
      real*8 tdw, tup, dh, fi, rad, drdw, dvdw, drup
      real*8 dvup,frdw,emis
      real*8 dfbu,dfbd,dfb
      real*8 fbu,fbd,fb
      real*8 drdwh,dvdwh,druph,dvuph,dr_eff
      real*8 dustsigmat
      real*8 pathlength, casab
      real*8 escape, escapeso,tau_i(mxinfph),escapeso_i(mxinfph)
c
      integer*4 line,series,atom
      integer*4 inl,nxlinestart,mxlinestart
c
      character jmod*4
c
c           Functions
c
      real*8 localout, tauline, fbowen
c
c
c     Check parameters and replace with sane values if necessary
c
c
c filling factor
c
      if ((fi.gt.1.0d0).and.(fi.le.0.0d0)) then
        write (*,10) fi
   10     format(//' INCONSISTENT FILLING FACTOR IN NEWDIF: ',1pg10.3)
        fi=1.0
      endif
      dvdw=dabs(dvdw)
      dvup=dabs(dvup)
c
      if ((tdw.lt.0.d0).or.(tup.lt.0.d0).or.(rad.lt.0.d0)
     &.or.(drdw.lt.0.d0).or.(drup.lt.0.d0)) then
        write (*,20) tdw,tup,rad,drdw,drup
   20   format(/,/,'ONE OF THE ARGUMENTS HAS A NEGATIVE VALUE IN NEWDIF'
     &         ,'****',/,5(1pg10.3))
        tdw=dabs(tdw)
        tup=dabs(tup)
        rad=dabs(rad)
        drdw=dabs(drdw)
        drup=dabs(drup)
      endif
c
      if ((jmod.eq.'LODW').or.(jmod.eq.'LOUP').or.(jmod.eq.'OUTW')
     &.or.(jmod.eq.'DWUP')) goto 40
      write (*,30) jmod
   30 format(//' MODE IMPROPERLY DEFINED IN NEWDIF : ',a4)
      stop
   40 continue
c
      if ((frdw.gt.1.d0).or.(frdw.lt.0.d0)) then
        write (*,50) frdw
   50       format(/,/,' INCONSISTENT VALUE FOR THE ARGUMENT FRDW :'
     &             ,1pg10.3)
        frdw=0.5d0
      endif
c
c     ***SET MODE AND INTERNAL DILUTION FACTORS
c
      dwf=frdw
      if (jmod.eq.'OUTW') dwf=1.0d0
      upf=1.d0-dwf
c
      dwex=1.d0
      if (jmod.eq.'LOUP') dwex=0.d0
      upex=1.d0-dwex
      if (jmod.eq.'DWUP') upex=1.d0
c
      if (jmod.eq.'LOUP') then
        dr=drup
        drem=drup
        dvem=dvup
        telc=tup
      else
        dr=drdw
        drem=drdw
        dvem=dvdw
        telc=tdw
      endif
c
      curad=rad
      dildw=1.d0
      dilup=1.d0
      dresc=(rad+0.5d0*dr-remp)
      if (jgeo.eq.'S') then
c
c     Spherical
c
        rm=curad+drdw
        wlo=(2.d0*dlog(curad))-(2.d0*dlog(rm))
        rm=(curad+dr)-drup
        if (rm.le.0.d0) then
          write (*,60) drup,curad
   60       format(/,/,' DRUP IS LARGER THAN THE RADIUS OF CURVATURE :'
     &              ,2(1pg10.3))
          stop
        endif
        ulo=(2.d0*dlog(curad+dr))-(2.d0*dlog(rm))
        dildw=dexp(wlo)
        dilup=dexp(ulo)
c
c     Spherical
c
      endif
c
      wadw=dildw
      wedw=1.d0
      waup=1.d0
      weup=dilup
      if (jmod.ne.'LOUP') goto 70
      wadw=1.d0
      wedw=dildw
      waup=dilup
      weup=1.d0
   70 if (jmod.ne.'DWUP') goto 80
      wadw=dildw
      wedw=1.d0
      waup=dilup
      weup=1.d0
   80 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Begin Routine Proper
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***COMPUTES RESONANCE-LINES DISTANCE MULTIPLIERS
      dr_eff=2.d0*(rad+0.5d0*dr-remp)
c      dr_eff = (rad + 0.5d0*dr - remp)
      dremh=0.5*dr_eff
      if (jgeo.eq.'S') then
        dremh=dmin1((rstromhb-remp),dremh)
      endif
      dvemh=dvem
      call emilindismuls (telc, dh, dremh, dvemh)
c
c   Hydrogenic lines
c
      fb=0.0d0
      fbu=0.0d0
      fbd=0.0d0
      if (zmap(2).gt.0) then
        fb=fbowen(telc,dvemh)
        fbu=fbowen(tup,dvup)
        fbd=fbowen(tdw,dvdw)
      endif
      heiioiiibf=fb
c
c     do up and down if needed
c
      if (dwf.gt.0.d0) then
        drdwh=0.5d0*dr_eff
        if (jgeo.eq.'S') then
          drdwh=0.5d0*dmin1((rstromhb-remp),drdwh)
        endif
        dvdwh=dvdw
        call dwlindismuls (tdw, dh, drdwh, dvdwh)
c dwf
      endif
c
      if (upf.gt.0.d0) then
        druph=0.5d0*dr_eff
        if (jgeo.eq.'S') then
          druph=0.5d0*dmin1((rstromhb-remp),druph)
        endif
        dvuph=dvup
        call uplindismuls (tup, dh, druph, dvuph)
c upf
      endif
c
c     call casab with hydrogenic Lyman gamma
c
      do atom=1,atypes
        line=3
        series=1
        caseab(atom)=casab(mapz(atom),line,series)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***UPDATES DIFFUSE FIELD VECTORS UPDIF, DWDIF etc BIN BY BIN
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nxlinestart=1
      mxlinestart=1
      pathlength=dr*fi
c
      do inl=1,infph-1
        tau_i(inl)=0.d0
        escapeso_i(inl)=1.d0
        energ=cphotev(inl)*ev
        tauso=0.d0
        sigmt=0.d0
        dustsigmat=0.d0
c
        call crosssections (inl, tauso, sigmt, dustsigmat)
c
        sigmt=dh*sigmt
        dustsigmat=dh*dustsigmat
c
c tauso already includes dh, fi and distance implicitly, so no simatso
c
c     ***ATTENUATION OF THE ORIGINAL INTENSITY OF THE DIFFUSE FIELD
c
        tau0=pathlength*sigmt
c
c tauso already includes dh, fi and distance implicitly
c
        escapeso=dexp(-tauso)
c

        f1 = (wadw*dexp(-(dwex*tau0)))
        f2 = (waup*dexp(-(upex*tau0)))
        dwdif(inl)=dwdif(inl)*f1
        updif(inl)=updif(inl)*f2
c
        dwdifcont(inl)=dwdifcont(inl)*f1
        updifcont(inl)=updifcont(inl)*f2
c
        tau_i(inl)=tau0
        escapeso_i(inl)=escapeso
      enddo
c
c     distance multipliers attenuate lines for emilin etc
c
      do line=1,nxr3lines
        inl=xr3lines_bin(line)
        if (inl.gt.0) then
          if ((xr3lines_dwlin(1,line)+xr3lines_uplin(1,line))
     &     .gt.epsilon) then
            dismul=(dwex*xr3lines_dwlin(2,line))+(upex*xr3lines_uplin(2,
     &       line))
            tau=tauline(dismul,tau_i(inl))
            xr3lines_dwlin(1,line)=xr3lines_dwlin(1,line)*(wadw*dexp(-
     &       (dwex*tau)))
            xr3lines_uplin(1,line)=xr3lines_uplin(1,line)*(waup*dexp(-
     &       (upex*tau)))
          endif
        endif
      enddo
c
      do line=1,nxrllines
        inl=xrllines_bin(line)
        if (inl.gt.0) then
          if ((xrllines_dwlin(1,line)+xrllines_uplin(1,line))
     &     .gt.epsilon) then
            dismul=(dwex*xrllines_dwlin(2,line))+(upex*xrllines_uplin(2,
     &       line))
            tau=tauline(dismul,tau_i(inl))
            xrllines_dwlin(1,line)=xrllines_dwlin(1,line)*(wadw*dexp(-
     &       (dwex*tau)))
            xrllines_uplin(1,line)=xrllines_uplin(1,line)*(waup*dexp(-
     &       (upex*tau)))
          endif
        endif
      enddo
c
C     do line=1,xlines
C       inl=xbin(line)
C       if (inl.gt.0) then
C         if ((dwlin(1,line)+uplin(1,line)).gt.epsilon) then
C           dismul=(dwex*dwlin(2,line))+(upex*uplin(2,line))
C           tau=tauline(dismul,tau_i(inl))
C           dwlin(1,line)=dwlin(1,line)*(wadw*dexp(-(dwex*tau)))
C           uplin(1,line)=uplin(1,line)*(waup*dexp(-(upex*tau)))
C         endif
C       endif
C     enddo
c
c   Hydrogenic lines
c
      do series=1,nhseries
        do line=1,nhlines
          inl=hbin(line,series)
          if (inl.gt.0) then
            if ((hyddwlin(1,line,series)+hyduplin(1,line,series))
     &       .gt.epsilon) then
c Hydrogenic series use dismul to lookup case A/B
c modified ratios already, so not needed for transfer here
              dismul=1.d0
              tau=tauline(dismul,tau_i(inl))
              hyddwlin(1,line,series)=hyddwlin(1,line,series)*(wadw*
     &         dexp(-(dwex*tau)))
              hyduplin(1,line,series)=hyduplin(1,line,series)*(waup*
     &         dexp(-(upex*tau)))
            endif
          endif
        enddo
      enddo
c
c   Helium lines
c
      do series=1,nheseries
        do line=1,nhelines
          inl=hebin(line,series)
          if (inl.gt.0) then
            if ((heldwlin(1,line,series)+heluplin(1,line,series))
     &       .gt.epsilon) then
              dismul=1.d0
              tau=tauline(dismul,tau_i(inl))
c
              dfbu=1.d0
              dfbd=1.d0
c             if ((line.eq.1).and.(series.eq.1)) then
c               dfbu=1.d0-fbu
c               dfbd=1.d0-fbd
c               write(*,*)'dfbu dfbd',dfbu,dfbd
c             endif
c
              heldwlin(1,line,series)=dfbd*heldwlin(1,line,series)*
     &         (wadw*dexp(-(dwex*tau)))
              heluplin(1,line,series)=dfbu*heluplin(1,line,series)*
     &         (waup*dexp(-(upex*tau)))
            endif
          endif
        enddo
      enddo
c
      do atom=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            inl=xhbin(line,series,atom)
            if (inl.gt.0) then
              if ((xhyddwlin(1,line,series,atom)+xhyduplin(1,line,
     &         series,atom)).gt.epsilon) then
                dismul=1.d0
                tau=tauline(dismul,tau_i(inl))
                xhyddwlin(1,line,series,atom)=xhyddwlin(1,line,series,
     &           atom)*(wadw*dexp(-(dwex*tau)))
                xhyduplin(1,line,series,atom)=xhyduplin(1,line,series,
     &           atom)*(waup*dexp(-(upex*tau)))
              endif
            endif
          enddo
        enddo
      enddo
c
c     ***ADS LOCAL EMISSIVITY CONTRIBUTION TO THE NEW DIFFUSE FIELD
c
      do inl=1,infph-1
        if (emidif(inl).gt.epsilon) then
          escape=localout(1.d0,tau_i(inl))
          energ=cphotev(inl)*ev
          emis=energ*pathlength
          dwem=(dwf*emis)*escape
          upem=(upf*emis)*escapeso_i(inl)
          dwdif(inl)=dwdif(inl)+(dwem*wedw*emidif(inl))
          updif(inl)=updif(inl)+(upem*weup*emidif(inl))
          dwdifcont(inl)=dwdifcont(inl)+(dwem*wedw*emidifcont(inl))
          updifcont(inl)=updifcont(inl)+(upem*weup*emidifcont(inl))
        endif
      enddo
c
      do line=1,nxr3lines
        inl=xr3lines_bin(line)
        if (inl.gt.0) then
          if (xr3lines_emilin(1,line).gt.epsilon) then
            dismul=(dwex*xr3lines_dwlin(2,line))+(upex*xr3lines_uplin(2,
     &       line))
            escape=localout(dismul,tau_i(inl))
            energ=xr3lines_egij(line)
            emis=energ*xr3lines_emilin(1,line)*pathlength
            dwliem=(dwf*emis)*escape
            upliem=(upf*emis)*escapeso_i(inl)
            xr3lines_dwlin(1,line)=xr3lines_dwlin(1,line)+(dwliem*wedw)
            xr3lines_uplin(1,line)=xr3lines_uplin(1,line)+(upliem*weup)
          endif
        endif
      enddo
c
      do line=1,nxrllines
        inl=xrllines_bin(line)
        if (inl.gt.0) then
          if (xrllines_emilin(1,line).gt.epsilon) then
            dismul=(dwex*xrllines_dwlin(2,line))+(upex*xrllines_uplin(2,
     &       line))
            escape=localout(dismul,tau_i(inl))
            energ=xrllines_egij(line)
            emis=energ*xrllines_emilin(1,line)*pathlength
            dwliem=(dwf*emis)*escape
            upliem=(upf*emis)*escapeso_i(inl)
            xrllines_dwlin(1,line)=xrllines_dwlin(1,line)+(dwliem*wedw)
            xrllines_uplin(1,line)=xrllines_uplin(1,line)+(upliem*weup)
          endif
        endif
      enddo
c
c
c     do emilin resonance lines
c
C     do line=1,xlines
C       inl=xbin(line)
C       if (inl.gt.0) then
C         if (emilin(1,line).gt.epsilon) then
C           dismul=(dwex*dwlin(2,line))+(upex*uplin(2,line))
C           escape=localout(dismul,tau_i(inl))
C           energ=ev*xejk(line)
C           emis=energ*emilin(1,line)*pathlength
C           dwliem=(dwf*emis)*escape
C           upliem=(upf*emis)*escapeso_i(inl)
C           dwlin(1,line)=dwlin(1,line)+(dwliem*wedw)
C           uplin(1,line)=uplin(1,line)+(upliem*weup)
C         endif
C       endif
C     enddo
c
c Hydrogen
c
      do series=1,nhseries
        do line=1,nhlines
          inl=hbin(line,series)
          if (inl.gt.0) then
            if (hydlin(1,line,series).gt.epsilon) then
c
c Hydrogenic series use dismul to lookup case A/B
c modified ratios already, so not needed for transfer here
c
              dismul=1.d0
              escape=localout(dismul,tau_i(inl))
              energ=ev*lmev/hlambda(line,series)
              emis=hydlin(1,line,series)*pathlength*energ
              dwliem=(dwf*emis)*escape
              upliem=(upf*emis)*escapeso_i(inl)
              hyddwlin(1,line,series)=hyddwlin(1,line,series)+(dwliem*
     &         wedw)
              hyduplin(1,line,series)=hyduplin(1,line,series)+(upliem*
     &         weup)
            endif
          endif
        enddo
      enddo
c
c Helium
c
      do series=1,nheseries
        do line=1,nhelines
          inl=hebin(line,series)
          if (inl.gt.0) then
            if (hellin(1,line,series).gt.epsilon) then
              dismul=1.d0
              escape=localout(dismul,tau_i(inl))
              energ=ev*lmev/helambda(line,series)
c
              dfb=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfb=1.d0-fb
c                write(*,*)'dfb',dfb
              endif
c
c degrade local HeII 303 photons locally to OIII BF lines by an
c approximate conversion fraction....  Set to 1-0 = 1.0, no
c conversion for now
c
              emis=dfb*energ*hellin(1,line,series)*pathlength
              dwliem=(dwf*emis)*escape
              upliem=(upf*emis)*escapeso_i(inl)
              heldwlin(1,line,series)=heldwlin(1,line,series)+(dwliem*
     &         wedw)
              heluplin(1,line,series)=heluplin(1,line,series)+(upliem*
     &         weup)
            endif
          endif
        enddo
      enddo
c
c Heavy Hydrogenic
c
      do atom=3,atypes
        do series=1,nxhseries
          do line=1,nxhlines
            inl=xhbin(line,series,atom)
            if (inl.gt.0) then
              if (xhydlin(1,line,series,atom).gt.epsilon) then
                dismul=1.d0
                escape=localout(dismul,tau_i(inl))
                energ=ev*lmev/xhlambda(line,series,atom)
                emis=energ*xhydlin(1,line,series,atom)*pathlength
                dwliem=(dwf*emis)*escape
                upliem=(upf*emis)*escapeso_i(inl)
                xhyddwlin(1,line,series,atom)=xhyddwlin(1,line,series,
     &           atom)+(dwliem*wedw)
                xhyduplin(1,line,series,atom)=xhyduplin(1,line,series,
     &           atom)+(upliem*weup)
              endif
            endif
          enddo
        enddo
      enddo
c
      return
      end
