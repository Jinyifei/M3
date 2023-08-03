cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function tauline(dismul, tau0)
c
c tau * dismul subject to an upper limit of 1.d10
c
      include 'const.inc'
c
      real*8 dismul,  tau0
      real*8 tau_eff
c
      tau_eff=tau0
c
      if (dismul.gt.1.d0) then
        tau_eff=tau0*dismul
      endif
c
      tauline=min(tau_eff,1.d10)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function transferout(dismul, tau)
c
c  Simply the attenuation of radition across a depth tau
c  with no sources along the path.
c  The fraction retained in the zone (ie absorbed) is
c  1-transferout.
c
      include 'const.inc'
c
      real*8 dismul, tau, tau_eff
      real*8 escape, tauline
c
      tau_eff=tauline(dismul,tau)
c
      escape=dexp(-tau_eff)
c
      transferout=escape
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function localout(dismul, tau)
c
c  A zone producing a uniform field.  Emission at one side
c sees a different optical depth to escape than at the other for a given directi
c The local production of energy escaping sees an average of exp(-tau)
c and the fraction that escapes is (1-exp(-tau))/tau.  The fraction
c retained in the zone is 1-escape
c
      include 'const.inc'
c
      real*8 dismul,  tau , tau_eff
      real*8 escape, tauline
c
      tau_eff=tauline(dismul,tau)
c
      escape=1.d0
      if (tau_eff.ge.1.0d-6) then
        escape=(1.d0-dexp(-tau_eff))/tau_eff
      endif
c
      localout=escape
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function meanfield(dismul, tau)
c
c A field crossing a zone depth tau, attenuated to exp(-tau)
c at the end has a mean field in the zone of (1-exp(tau))/tau
c Exactly the same as the one-side escape fraction of a locally produced
c field  Called meanfield just to distinguish the case from local
c production.
c
c The difference is that the amount escaping is exp(-tau), not
c 1 - meanfield.  This function averages a crossing field with
c no local sources.
c
      include 'const.inc'
c
      real*8 dismul,  tau , tau_eff
      real*8 mean, tauline
c
      tau_eff=tauline(dismul,tau)
c
      mean=1.d0
      if (tau_eff.ge.1.0d-6) then
        mean=(1.d0-dexp(-tau_eff))/tau_eff
      endif
c
      meanfield=mean
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
c     This is the new function to replace restrans. It will calculate
c     the radiative line transfer distance multiplier for the new
c     XLINDAT short wavelength resonance lines <2000 Angstoms.
c
c     Old resonance lines with Lambda longer than this will be simplly
c     added to the diffuse field like the inter-combination lines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******COMPUTES THE RESONANCE LINE TRANSFER USING ESCAPE
c   PROBABILITY FORMULATION
c   DV IS THE CHANGE IN VELOCITY OF THE GAS FROM R TO R+DR
c   RETURNS EFFECTIVE STEP LENGTH MULTIPLIER FOR THE LINE.
c     **REF. : CAPRIOTTI,E.R.,1965 APJ.142,1101 (EQ* 88A,B)
c          (SIGN ERRORS IN 88A HAVE BEEN CORRECTED)
c
      real*8 function fdismul(t,dh,dr,dv,atom,ion,ejk,fab)
c
c           ejk in ergs as usual
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 abu,con,conte,cte,dreff,dveff,rkte,tau
      real*8 tauther,vra,vther,x
      real*8 t, dh, dr, dv, ejk, fab , telc
      real*8 e7a,e7b,argu,erf
      real*8 a(5),f(3),e
      real*8 z, p, y
c
      integer*4 atom,ion
c
      e7a(x)=(1.0d0+((0.71d0*x)*dlog((2.0d0*x)+1.d-37)))+(((((((((a(1)*
     &x)+a(2))*x)+a(3))*x)+a(4))*x)+a(5))*x)
      e7b(x)=((0.14d0+(0.25d0/dsqrt(dlog(2.0d0*x))))+dsqrt(dlog(2.0d0*x)
     &))/(3.5449077d0*x)
      argu(z)=1.0d0/(1.0d0+(0.47047d0*z))
      erf(z)=1.0d0-(((f(1)+((f(2)+f(3)*argu(z))*argu(z)))*argu(z))*
     &dexp(-(z*z)))
c
c     constants for some functions
c
      data a/0.00048d0,-0.019876d0,0.08333d0,-0.3849d0,-0.83d0/
      data f/0.3480242d0,-0.0958798d0,0.7478556d0/
c
      e=1.d0
c
      z=zion(atom)
      p=pop(ion,atom)
c
      if (z*p.gt.1.0d-6) then
c
        dveff=dmax1(dabs(dv),1.0d6)
c      dveff = dabs(dv)
        dreff=dabs(dr)
        telc=dmax1(t,10.0d0)
        rkte=rkb*telc
        cte=rkte/mp
        con=3.0d-18
c
        y=(ejk/rkte)
        if (y.lt.maxdekt) then
          conte=con*(1.0d0-dexp(-y))
        else
          conte=con
        endif
c
        vther=dsqrt((2.0d0*cte)/atwei(atom))
        dveff=dmax1(dveff,vther)
        abu=dh*z*p*conte*fab*dreff/ejk
c
        tauther=abu/vther
        vra=dveff/vther
c
        tau=tauther
        if (vra.ge.0.01d0) tau=tauther*(erf(vra)/vra)*0.8862564664d0
        if (tau.lt.1.5d0) e=e7a(tau)
        if (tau.ge.1.5d0) e=e7b(tau)
c
        e=dmax1(1.d0,1.d0/(e+1.0d-6))
c
      endif
c
      fdismul=e
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine emilindismuls (telc, dh, dremh, dvemh)
c
      include 'cblocks.inc'
      real*8 telc,dh,dremh,dvemh
      integer line, atom, ion, series,j
      real*8 p,z,f,es,dr_eff
c
      real*8 fdismul
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
C       if (emilin(1,line).gt.epsilon) then
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
c  H
c
      atom=1
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      dr_eff=dremh
      if (jgeo.eq.'S') then
        dr_eff=dmin1((rstromhb-remp),dremh)
      endif
      do series=1,nhseries
        do line=1,nhlines
          hyduplin(2,line,series)=1.d0
          hyddwlin(2,line,series)=1.d0
          hydlin(2,line,series)=1.d0
          j=hbin(line,series)
          if (j.ne.0) then
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              es=lmev/hlambda(line,series)*ev
              f=hydrogf(line,series)
              hydlin(2,line,series)=fdismul(telc,dh,dr_eff,dvemh,atom,
     &         ion,es,f)
            endif
          endif
        enddo
      enddo
c
      atom=zmap(2)
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      dr_eff=dremh
      if (jgeo.eq.'S') then
        dr_eff=dmin1((rstromheb-remp),dremh)
      endif
      do series=1,nheseries
        do line=1,nhelines
          heluplin(2,line,series)=1.d0
          heldwlin(2,line,series)=1.d0
          hellin(2,line,series)=1.d0
          j=hebin(line,series)
          if (j.ne.0) then
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              es=lmev/helambda(line,series)*ev
              f=hydrogf(line,series)
              hellin(2,line,series)=fdismul(telc,dh,dr_eff,dvemh,atom,
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
            xhyduplin(2,line,series,atom)=1.d0
            xhyddwlin(2,line,series,atom)=1.d0
            xhydlin(2,line,series,atom)=1.d0
            j=xhbin(line,series,atom)
            if (j.ne.0) then
              if ((z*p.ge.pzlimit).and.(series.eq.1)) then
                f=hydrogf(line,series)
                es=lmev/xhlambda(line,series,atom)*ev
                xhydlin(2,line,series,atom)=fdismul(telc,dh,dremh,dvemh,
     &           atom,ion,es,f)
              endif
            endif
          enddo
        enddo
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uplindismuls (tup, dh, druph, dvuph)
c
      include 'cblocks.inc'
      real*8 tup,dh,druph,dvuph
      integer line, atom, ion, series
      real*8 p,z,f,es, dr_eff
c
      real*8 fdismul
c
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
            if (z*p.ge.pzlimit) then
              f=xr3lines_gf(line)
              es=xr3lines_egij(line)
              xr3lines_uplin(2,line)=fdismul(tup,dh,druph,dvuph,atom,
     &         ion,es,f)
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
            if (z*p.ge.pzlimit) then
              f=xrllines_gf(line)
              es=xrllines_egij(line)
              xrllines_uplin(2,line)=fdismul(tup,dh,druph,dvuph,atom,
     &         ion,es,f)
            endif
          endif
        endif
      enddo
c
C     do line=1,xlines
C       uplin(2,line)=1.d0
C       if (uplin(1,line).gt.epsilon) then
C         z=zion(xrat(line))
C         p=pop(xion(line),xrat(line))
C         if ((z*p.ge.pzlimit)) then
C           f=xfef(line)/xbr(line)
C           es=xejk(line)*ev
C           uplin(2,line)=fdismul(tup,dh,druph,dvuph,xrat(line),
C    &       xion(line),es,f)
C         endif
C       endif
C     enddo
c
c   Hydrogenic lines multipliers for case A/B interpolation
c
      dr_eff=druph
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromhb-remp),druph)
      endif
      atom=zmap(1)
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nhseries
        do line=1,nhlines
          hyduplin(2,line,series)=1.d0
          if ((z*p.ge.pzlimit).and.(series.eq.1)) then
            es=(lmev/hlambda(line,series))*ev
            f=hydrogf(line,series)
            hyduplin(2,line,series)=fdismul(tup,dh,dr_eff,dvuph,atom,
     &       ion,es,f)
          endif
        enddo
      enddo
c
      dr_eff=druph
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromheb-remp),druph)
      endif
      atom=zmap(2)
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nheseries
        do line=1,nhelines
          heluplin(2,line,series)=1.d0
          if ((z*p.ge.pzlimit).and.(series.eq.1)) then
            es=(lmev/helambda(line,series))*ev
            f=hydrogf(line,series)
            heluplin(2,line,series)=fdismul(tup,dh,dr_eff,dvuph,atom,
     &       ion,es,f)
          endif
        enddo
      enddo
c
      dr_eff=druph
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromheb-remp),druph)
      endif
      do atom=3,atypes
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
        do series=1,nxhseries
          do line=1,nxhlines
            xhyduplin(2,line,series,atom)=1.d0
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              es=lmev/xhlambda(line,series,atom)*ev
              f=hydrogf(line,series)
              xhyduplin(2,line,series,atom)=fdismul(tup,dh,dr_eff,dvuph,
     &         atom,ion,es,f)
            endif
          enddo
        enddo
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dwlindismuls (tdw, dh, drdwh, dvdwh)
c
      include 'cblocks.inc'
      real*8 tdw,dh,drdwh,dvdwh
      integer line, series, atom, ion
      real*8 p,z,f,es,dr_eff
c
      real*8 fdismul
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
     &         ion,es,f)
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
     &         ion,es,f)
            endif
          endif
        endif
      enddo
c
C     do line=1,xlines
C       dwlin(2,line)=1.d0
C       if (dwlin(1,line).gt.epsilon) then
C         z=zion(xrat(line))
C         p=pop(xion(line),xrat(line))
C         if ((z*p.ge.pzlimit)) then
C           f=xfef(line)/xbr(line)
C           es=xejk(line)*ev
C           dwlin(2,line)=fdismul(tdw,dh,drdwh,dvdwh,xrat(line),
C    &       xion(line),es,f)
C         endif
C       endif
C     enddo
c
c
c   Hydrogenic lines
c
c
      dr_eff=drdwh
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromhb-remp),drdwh)
      endif
      atom=zmap(1)
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nhseries
        do line=1,nhlines
          hyddwlin(2,line,series)=1.d0
          if ((z*p.ge.pzlimit).and.(series.eq.1)) then
            f=hydrogf(line,series)
            es=(lmev/hlambda(line,series))*ev
            hyddwlin(2,line,series)=fdismul(tdw,dh,dr_eff,dvdwh,atom,
     &       ion,es,f)
          endif
        enddo
      enddo
c
      dr_eff=drdwh
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromheb-remp),drdwh)
      endif
      atom=zmap(2)
      ion=mapz(atom)
      z=zion(atom)
      p=pop(ion,atom)
      do series=1,nheseries
        do line=1,nhelines
          heldwlin(2,line,series)=1.d0
          if ((z*p.ge.pzlimit).and.(series.eq.1)) then
c             if (heldwlin(1,line,series).gt.epsilon) then
            f=hydrogf(line,series)
            es=(lmev/helambda(line,series))*ev
            heldwlin(2,line,series)=fdismul(tdw,dh,dr_eff,dvdwh,atom,
     &       ion,es,f)
c             endif
          endif
        enddo
      enddo
c
      dr_eff=drdwh
      if (jgeo.eq.'S') then
        dr_eff=0.5d0*dmin1((rstromheb-remp),drdwh)
      endif
      do atom=3,atypes
        ion=mapz(atom)
        z=zion(atom)
        p=pop(ion,atom)
        do series=1,nxhseries
          do line=1,nxhlines
            xhyddwlin(2,line,series,atom)=1.d0
            if ((z*p.ge.pzlimit).and.(series.eq.1)) then
              es=lmev/xhlambda(line,series,atom)*ev
              f=hydrogf(line,series)
              xhyddwlin(2,line,series,atom)=fdismul(tdw,dh,dr_eff,dvdwh,
     &         atom,ion,es,f)
            endif
          enddo
        enddo
      enddo
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
c     Approximate Bowen Flourescence losses to HeII303  tau_HeLa>>>Tau_c
c   Kallman & MacCray 1980
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      real*8 function fbowen(t, dv)
c
c           ejk in ergs as usual
c
c     Gives the fraction of He II 303A that is converted to Bowen lines
c
      include 'cblocks.inc'
c
      real*8 t, dv
c
      real*8 aba, aboiii, abheii, abhei
      real*8 dveff,vra8,vra2,vther8,vther2
      real*8 telc
      real*8 con,fb,fb2,fb3
c
      real*8 taucr, taula,f(3), argu, erf, z
c
      integer*4 atom
c
      data f/0.3480242d0,-0.0958798d0,0.7478556d0/
c
      argu(z)=1.0d0/(1.0d0+(0.47047d0*z))
      erf(z)=1.0d0-(((f(1)+((f(2)+f(3)*argu(z))*argu(z)))*argu(z))*
     &dexp(-(z*z)))
c
c     constants for some functions
c
      fb=0.d0
c
      dveff=dmax1(dabs(dv),2.0d6)
      telc=dmax1(t,10.0d0)
c
      con=2.0d-3
c
      fb=0.d0
      fb2=0.d0
      fb3=0.d0
c
      atom=zmap(8)
      if (atom.gt.0) then
c
        aboiii=zion(atom)*pop(3,atom)
c
        if (aboiii>1.d-6) then
          vther8=dsqrt((2.0d0*rkb*telc)/(amu*atwei(zmap(8))))
          vther2=dsqrt((2.0d0*rkb*telc)/(amu*atwei(zmap(2))))
          vra8=dveff/vther8
          vra2=dveff/vther2
          atom=zmap(2)
          abhei=zion(2)*pop(1,atom)
          abheii=zion(atom)*pop(2,atom)
          taula=2.0d6*abheii
          if (vra2.gt.0.01d0) taula=taula*(erf(vra2)/vra2)*
     &     0.8862564664d0
c
          atom=zmap(1)
          aba=zion(atom)*pop(1,atom)+9.3d0*abhei
c
          taucr=80.d0*abheii/aboiii
c
          if (taucr.lt.1.0d6) then
c
            fb2=dmax1(0.0d0,1.0d0/(1.d0+con*(aba/aboiii)))
            fb3=dmax1(0.0d0,1.0d0/(1.d0+(taucr/taula)))
c const normalises erfx/x:
            fb=(fb2*fb3)*(erf(vra8)/vra8)*0.8862564664d0
c
          endif
        endif
c
      endif
c
      fbowen=0.d0
c fb disabled until rederived and re-implemented
c
c This can't really be solved locally, it is a global nebula
c property, and not fully compatible with the step by step
c outward only approximation.  It needs more research and probably
c a new integration method.  If turned on then it seems to give
c values, but I am not confident they mean anything I can understand.
c
      return
c
      end
