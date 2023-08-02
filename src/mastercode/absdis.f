cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES THE DISTANCE : DRTA  TO OBTAIN
c       A GIVEN TOTAL PHOTON ABSORBTION FRACTION FROM TPHOT
c       TAUAV (MMOD='DIS','LIN') takes account of dust
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine absdis ( dh, fi, absf, drta, rad, popul)
c
c
      include 'cblocks.inc'
c
c
      real*8 popul(mxion, mxelem), popt(mxion, mxelem)
      real*8 wid,wei,phots
      real*8 dh, fi, absf, drta, rad,pcros
      real*8 plos,qto, de, g1, g2, temp
      real*8 nfn(mxinfph),r2,r3,cebin,sig
      integer*4 bincount,inl,i,ie,j,k,dtype,dstmin
c
      real*8 ptime,rtime,ctime,abio,crosec
c
c       External Functions
c
      real*8 feldens,fphotim,frectim2,fcolltim,fdilu
c
c       Internal Functions
c
c      real*8 acrs,eph,at,bet,se
c
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (drta.eq.0.d0) drta=1.0d16
c
c
      call copypop (pop, popt)
      call copypop (popul, pop)
c
      de=feldens(dh,popul)
c
c      blum = 0.d0
c      do inl = 1, infph-1
c         cebin = cphotev(inl)
c         blum = blum+tphot(inl)*cebin*evplk
c      enddo
c
      ptime=fphotim()
      rtime=frectim2(de)
      ctime=fcolltim(de)
c
      write (*,*) ptime,rtime,ctime
c      r4 = rtime/(ptime+rtime)
c      write(*,*) 'Recombination Correction : ',r4
c      absf = absf*r4
c      write(*,*) 'Absorbsion Fraction : ',absf
c
      call copypop (popt, pop)
c
      plos=0.d0
      xsect=0.d0
      bincount=0
      qto=0.d0
      do inl=1,infph-1
        xsec(inl)=0.d0
        if (skipbin(inl)) goto 20
        cebin=cphotev(inl)
        nfn(inl)=tphot(inl)*cebin*evplk
        sig=0.0d0
        do 10 i=1,ionum
          if (photxsec(i,inl).gt.0.d0) then
            ie=atpho(i)
            j=ionpho(i)
            abio=zion(ie)*popul(j,ie)
            if (abio.gt.1.d-12) then
              crosec=photxsec(i,inl)
              if (crosec.gt.0.d0) sig=sig+(abio*crosec)
c          eph=cebin/ipotpho(i)
c          if (eph.lt.1.0d0) goto 10
c          crosec=acrs(eph,sigpho(i),betpho(i),spho(i))
c          sig=sig+(abio*crosec)
            endif
          endif
   10   continue
c
        if (sig.lt.0.d0) sig=0.d0
        sig=(dh*fi)*sig
        xsec(inl)=sig
        xsect=xsect+sig
   20   continue
      enddo
c
c  dust
c
      if (grainmode.eq.1) then
c
c   pahs
c
        if (pahactive.eq.1) then
          do inl=1,infph-1
            sig=0.d0
            do k=1,pahi
              if (pahion(k).gt.0) then
                pcros=pahz(k)*pahiext(inl)*pahfrac
              else
                pcros=pahz(k)*pahnext(inl)*pahfrac
              endif
              sig=sig+pcros
            enddo
            xsec(inl)=xsec(inl)+sig*dh*fi
            xsect=xsect+sig*dh*fi
          enddo
        endif
c
c  grains
c
        dstmin=1
        if ((clinpah.eq.1).and.(pahactive.eq.0)) dstmin=2
c
        do dtype=dstmin,numtypes
          do inl=1,infph-1
            sig=fi*dh*dcrosec(inl,dtype)
            xsec(inl)=xsec(inl)+sig
            xsect=xsect+sig
          enddo
        enddo
      endif
c
c      write(*,*) 'Total X-section : ',xsect
c
      do inl=1,infph-1
        wei=xsec(inl)/xsect
        wid=(photev(inl+1)-photev(inl))*evplk
        phots=(tphot(inl))*wid*wei
c         write(*,*) inl,phots,qto
c         write(*,*) tphot(inl),xsec(inl),xsect
        qto=qto+phots
      enddo
c
c      write(*,*) 'Weighted total : ',qto
c
      iter=0
   30 plos=0.d0
c
      do inl=1,infph-1
        wei=xsec(inl)/xsect
        wid=(photev(inl+1)-photev(inl))*evplk
        phots=(tphot(inl))*wid*wei
        plos=plos+(1.d0-dexp(-drta*xsec(inl)))*(phots/qto)
      enddo
c
      iter=iter+1
      if (plos.gt.0.d0) then
        r2=absf/plos
        r3=dabs(1.d0-r2)
c        write(*,*) 'Change in dr:',r2
        if (r3.gt.0.01d0) then
          drta=drta*r2
          goto 30
        endif
      else
        r2=1.d0/epsilon
        drta=r2
      endif
c
c
c     diagnosic printouts, normally commented out.
c
c      caller = 'AD'
c      pfx = 'xsect'
c      np = 5
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,xsec)
c      caller = 'AD'
c      pfx = 'nfn'
c      np = 3
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,nfn)
c      caller = 'AD'
c      pfx = 'tphot'
c      np = 5
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,tphot)
c
c
      if (plos.gt.0.d0) drta=drta*absf/plos
c
c     correct for geometric dilution changes if necessary
c
      if ((jgeo.eq.'S').or.(jgeo.eq.'F')) then
        g1=fdilu(rstar,rad)
        r2=rad+drta
        g2=fdilu(rstar,r2)
        temp=dsqrt(1.23456789d0*g2/g1)
        r2=(r2*temp)-rad
        if (r2.lt.0) r2=drta
        r3=1.d0/r2+1.d0/drta
        drta=1.d0/r3
c         write(*,*) 'Geometric Correction:',temp
      endif
      if (drta.lt.0.d0) then
        write (*,*) 'dr is negative:',drta
        stop
      endif
c
c      drta = max(1.0e13, drta)
c
c      drta = min(5.0e14,drta)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES THE DISTANCE : DRTA  TO OBTAIN
c       A GIVEN TOTAL PHOTON ABSORBTION FRACTION FROM TPHOT
c       TAUAV (MMOD='DIS','LIN') takes account of dust
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine absdis2 (dh, fi, absf, drta, rad, popul)
c
      include 'cblocks.inc'
c
      real*8 popul(mxion, mxelem)
      real*8 wid,wei,phots
      real*8 dh, fi, absf, drta, rad, pcros
      real*8 plos,qto
c      real*8 g1, g2
      real*8 r1,r2,r3,cebin,sig
      integer*4 bincount,inl,i,ie,j,k,dtype,dstmin
c
      real*8 abio,crosec
c
c       External Functions
c
c      real*8 fdilu
c
c       Internal Functions
c
c      real*8 acrs,eph,at,bet,se
c
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (drta.eq.0.d0) drta=1.0d16
c
      plos=0.d0
      xsect=0.d0
      bincount=0
      qto=0.d0
      do inl=1,infph-1
        xsec(inl)=0.d0
        if (skipbin(inl)) goto 20
        cebin=cphotev(inl)
        sig=0.0d0
        do 10 i=1,ionum
          if (photxsec(i,inl).gt.0.d0) then
            ie=atpho(i)
            j=ionpho(i)
            abio=zion(ie)*popul(j,ie)
            if (abio.gt.1.d-12) then
              crosec=photxsec(i,inl)
              if (crosec.gt.0.d0) sig=sig+(abio*crosec)
c          eph=cebin/ipotpho(i)
c          if (eph.lt.1.0d0) goto 10
c          crosec=acrs(eph,sigpho(i),betpho(i),spho(i))
c          sig=sig+(abio*crosec)
            endif
          endif
   10   continue
c
        if (sig.lt.0.d0) sig=0.d0
        sig=(dh*fi)*sig
        xsec(inl)=sig
        xsect=xsect+sig
   20   continue
      enddo
c
c  dust
c
      if (grainmode.eq.1) then
c
c   pahs
c
        if (pahactive.eq.1) then
          do inl=1,infph-1
            sig=0.d0
            do k=1,pahi
              if (pahion(k).gt.0) then
                pcros=pahz(k)*pahiext(inl)*pahfrac
              else
                pcros=pahz(k)*pahnext(inl)*pahfrac
              endif
              sig=sig+pcros
            enddo
            xsec(inl)=xsec(inl)+sig*dh*fi
            xsect=xsect+sig*dh*fi
          enddo
        endif
c
c  grains
c
        dstmin=1
        if ((clinpah.eq.1).and.(pahactive.eq.0)) dstmin=2
c
        do dtype=dstmin,numtypes
          do inl=1,infph-1
            sig=fi*dh*dcrosec(inl,dtype)
            xsec(inl)=xsec(inl)+sig
            xsect=xsect+sig
          enddo
        enddo
      endif
c
      do inl=1,infph-1
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          wid=(photev(inl+1)-photev(inl))*evplk
          phots=(tphot(inl))*wid*wei
          qto=qto+phots
        endif
      enddo
c
      iter=0
   30 plos=0.d0
c
      do inl=1,infph-1
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          wid=(photev(inl+1)-photev(inl))*evplk
          phots=(tphot(inl))*wid*wei
          plos=plos+(1.d0-dexp(-drta*xsec(inl)))*phots
        endif
      enddo
c
      plos=plos/qto
c
      iter=iter+1
      if (plos.gt.0.d0) then
        r2=absf/plos
        r3=dabs(1.d0-r2)
        if (r3.gt.0.01d0) then
          drta=drta*r2
          goto 30
        endif
      endif
c
      if (plos.gt.0.d0) drta=drta*absf/plos
c
c     correct for geometric changes if necessary in spheres
c
      if (jgeo.eq.'S') then
c
c 1/(1% of current radius
c
        r1=100.d0/(rad)
c
c harmonic mean, dominated by smallest value
c
        drta=1.d0/(r1+1.d0/drta)
c
c compute actual absorb fraction
c and return dtau which will make next step quicker
c
        plos=0.d0
c
        do inl=1,infph-1
          if (xsec(inl).gt.0.d0) then
            wei=xsec(inl)/xsect
            wid=(photev(inl+1)-photev(inl))*evplk
            phots=(tphot(inl))*wid*wei
            plos=plos+(1.d0-dexp(-drta*xsec(inl)))*phots
          endif
        enddo
c
        plos=plos/qto
        if (plos.gt.0.d0) absf=plos
c
      endif
c
c
      if (drta.lt.0.d0) then
        write (*,*) 'ABSDIS: dr is negative:',drta,plos,qto
        stop
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES THE attenuation and total absorption area for
c       a set of fraction columns (pop*r). fi , dr are not used unless
c       takes account of dust if present,
c       but only as a constant
c       over dr, no dust columns available yet.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine attsig (dh, fi, dr, popcols, attenuate, sigma)
c
      include 'cblocks.inc'
c
      real*8 dh, fi, dr
      real*8 popcols(mxion, mxelem)
      real*8 sigma(infph)
      real*8 attenuate(infph)
      real*8 sig,att,pcros
      integer*4 inl,i,ie,j,k,dtype,dstmin
c
      real*8 abio,crosec
c
      do inl=1,infph
        attenuate(inl)=1.0d0
        sigma(inl)=0.0d0
      enddo
c
      do inl=ionstartbin,infph-1
        sig=0.0d0
        do 10 i=1,ionum
          if (photxsec(i,inl).gt.0.d0) then
            ie=atpho(i)
            j=ionpho(i)
            abio=zion(ie)*popcols(j,ie)
            if (abio.gt.pzlimit) then
              crosec=photxsec(i,inl)
              sig=sig+(abio*crosec)
            endif
          endif
   10   continue
        if (sig.lt.0.d0) sig=0.d0
        sigma(inl)=dh*sig!popcolsincludesdralready
      enddo
c
c  dust
c
      if (grainmode.eq.1) then
c
c   pahs
c
        if (pahactive.eq.1) then
          do inl=1,infph-1
            sig=0.d0
            do k=1,pahi
              if (pahion(k).gt.0) then
                pcros=pahz(k)*pahiext(inl)*pahfrac
              else
                pcros=pahz(k)*pahnext(inl)*pahfrac
              endif
              sig=sig+pcros
            enddo
            sigma(inl)=sigma(inl)+(dh*fi*dr)*sig
          enddo
        endif
c
c  grains
c
        dstmin=1
        if ((clinpah.eq.1).and.(pahactive.eq.0)) dstmin=2
c
        do dtype=dstmin,numtypes
          do inl=1,infph-1
            sigma(inl)=sigma(inl)+(dh*fi*dr)*dcrosec(inl,dtype)
          enddo
        enddo
c
      endif
c
      do inl=1,infph-1
        if (sigma(inl).gt.0.d0) then
          att=dexp(-sigma(inl))
          attenuate(inl)=dmax1(dmin1(att,1.d0),0.d0)
c          if (attenuate(inl).lt.1.d0) write(*,*) inl, attenuate(inl)
        endif
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES THE weigthed fraction of the field absorbed by hydrogen
c       over drta  takes account of dust, assumes plane parallel geometry
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine habsfrac (absf, dh, fi, drta, popul)
c
      include 'cblocks.inc'
c
      real*8 popul(mxion, mxelem)
      real*8 wei, dh, fi, absf, drta, pcros
      real*8 plos,qto
c      real*8 g1, g2
      real*8 cebin,sig
      integer*4 bincount,inl,i,ie,j,k,dtype,dstmin
c
      real*8 abio,crosec
c
c       External Functions
c
c      real*8 fdilu
c
c       Internal Functions
c
c      real*8 acrs,eph,at,bet,se
c
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (drta.eq.0.d0) drta=1.0d16
c
      plos=0.d0
      xsect=0.d0
      bincount=0
      qto=0.d0
      do inl=1,infph-1
        xsec(inl)=0.d0
        if (skipbin(inl)) goto 20
        cebin=cphotev(inl)
        sig=0.0d0
        do 10 i=1,ionum
          if (photxsec(i,inl).gt.0.d0) then
            ie=atpho(i)
            j=ionpho(i)
            abio=zion(ie)*popul(j,ie)
            if (abio.gt.1.d-12) then
              crosec=photxsec(i,inl)
              if (crosec.gt.0.d0) sig=sig+(abio*crosec)
            endif
          endif
   10   continue
c
        if (sig.lt.0.d0) sig=0.d0
        sig=(dh*fi)*sig
        xsec(inl)=sig
        xsect=xsect+sig
   20   continue
      enddo
c
c  dust
c
      if (grainmode.eq.1) then
c
c   pahs
c
        if (pahactive.eq.1) then
          do inl=1,infph-1
            sig=0.d0
            do k=1,pahi
              if (pahion(k).gt.0) then
                pcros=pahz(k)*pahiext(inl)*pahfrac
              else
                pcros=pahz(k)*pahnext(inl)*pahfrac
              endif
              sig=sig+pcros
            enddo
            xsec(inl)=xsec(inl)+sig*dh*fi
            xsect=xsect+sig*dh*fi
          enddo
        endif
c
c  grains
c
        dstmin=1
        if ((clinpah.eq.1).and.(pahactive.eq.0)) dstmin=2
c
        do dtype=dstmin,numtypes
          do inl=1,infph-1
            sig=fi*dh*dcrosec(inl,dtype)
            xsec(inl)=xsec(inl)+sig
            xsect=xsect+sig
          enddo
        enddo
      endif
c
      do inl=1,infph-1
        if ( (cphotev(inl).gt.ipotev(1,1)).and.
     &       (cphotev(inl).lt.ipotev(1,2)) ) then
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          qto=qto+wei
        endif
        endif
      enddo
c
      absf=0.d0
c
      do inl=1,infph-1
        if ( (cphotev(inl).gt.ipotev(1,1)).and.
     &       (cphotev(inl).lt.ipotev(1,2)) ) then
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          absf=absf+(1.d0-dexp(-drta*xsec(inl)))*wei
        endif
        endif
      enddo
c
      absf=absf/qto
c
      return
      end
c
      subroutine absfrac (absf, dh, fi, drta, popul)
c
      include 'cblocks.inc'
c
      real*8 popul(mxion, mxelem)
      real*8 wid,wei,phots
      real*8 dh, fi, absf, drta, pcros
      real*8 plos,qto
c      real*8 g1, g2
      real*8 cebin,sig
      integer*4 bincount,inl,i,ie,j,k,dtype,dstmin
c
      real*8 abio,crosec
c
c       External Functions
c
c      real*8 fdilu
c
c       Internal Functions
c
c      real*8 acrs,eph,at,bet,se
c
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (drta.eq.0.d0) drta=1.0d16
c
      plos=0.d0
      xsect=0.d0
      bincount=0
      qto=0.d0
      do inl=1,infph-1
        xsec(inl)=0.d0
        if (skipbin(inl)) goto 20
        cebin=cphotev(inl)
        sig=0.0d0
        do 10 i=1,ionum
          if (photxsec(i,inl).gt.0.d0) then
            ie=atpho(i)
            j=ionpho(i)
            abio=zion(ie)*popul(j,ie)
            if (abio.gt.1.d-12) then
              crosec=photxsec(i,inl)
              if (crosec.gt.0.d0) sig=sig+(abio*crosec)
            endif
          endif
   10   continue
c
        if (sig.lt.0.d0) sig=0.d0
        sig=(dh*fi)*sig
        xsec(inl)=sig
        xsect=xsect+sig
   20   continue
      enddo
c
c  dust
c
      if (grainmode.eq.1) then
c
c   pahs
c
        if (pahactive.eq.1) then
          do inl=1,infph-1
            sig=0.d0
            do k=1,pahi
              if (pahion(k).gt.0) then
                pcros=pahz(k)*pahiext(inl)*pahfrac
              else
                pcros=pahz(k)*pahnext(inl)*pahfrac
              endif
              sig=sig+pcros
            enddo
            xsec(inl)=xsec(inl)+sig*dh*fi
            xsect=xsect+sig*dh*fi
          enddo
        endif
c
c  grains
c
        dstmin=1
        if ((clinpah.eq.1).and.(pahactive.eq.0)) dstmin=2
c
        do dtype=dstmin,numtypes
          do inl=1,infph-1
            sig=fi*dh*dcrosec(inl,dtype)
            xsec(inl)=xsec(inl)+sig
            xsect=xsect+sig
          enddo
        enddo
      endif
c
      do inl=1,infph-1
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          wid=(photev(inl+1)-photev(inl))*evplk
          phots=(tphot(inl))*wid*wei
          qto=qto+phots
        endif
      enddo
c
      absf=0.d0
c
      do inl=1,infph-1
        if (xsec(inl).gt.0.d0) then
          wei=xsec(inl)/xsect
          wid=(photev(inl+1)-photev(inl))*evplk
          phots=(tphot(inl))*wid*wei
          absf=absf+(1.d0-dexp(-drta*xsec(inl)))*phots
        endif
      enddo
c
      absf=absf/qto
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES THE AVERAGE DISTANCE : DRTA  TO OBTAIN
c     A GIVEN ABSORBTION : TAUAV (MMOD='DIS','LIN')
c     OR DERIVES EFFECTIVE TAU  :  TAUAV  FOR A GIVEN
c     DISTANCE : DRTA  (MMOD='TAU')
c
c     MMOD='DIS' IS HIGHLY NON LINEAR AND PRESUPPOSE EQUILIBRIUM
c     IONISATION , TAKES INTO EFFECT GEOMETRICAL DILUTION .
c
c     RAD : CURVATURE RADIUS AT THE POINT CONSIDERED
c     FI : FILLING FACTOR  ;  DH : PEAK H DENSITY
c     USES PHOTON FLUX : TPHOT  ,  IONIC POPULATIONS : POPUL(J,I)
c     AND PHOTOIONISATION CROSS SECTIONS OF ALL ELEMENTS
c     UP TO*LIMPH
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine taudist (dh, fi, tauav, drta, rad, popul, mmod)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a,abio,adto,ar,cmpos
      real*8 cof,convg,crosec,ctpos,divsi,dri,ecr
      real*8 expo,fg,gam,gammax,geo
      real*8 phma,pos,relma,relp
      real*8 rewe,rr,sc,sca,sca2,scaw,scaw2
      real*8 sig,sigm,sigmgeo,simax,steep,taunew,testab,u
      real*8 ulo,wei,weiad,weim2,weima,ratio,wlo
      real*8 radius, sigwei, sigtot
      real*8 sigen(mxinfph), popul(mxion, mxelem)
      real*8 relwei(mxinfph), rela(16)
      real*8 dh, fi, tauav, drta, rad
c
      integer*4 i,ie,inl
      integer*4 j,k,lim,n
c
      character mmod*4
c
      real*8 far
c      real*8 acrs,eph,at,bet,se
c
      far(u,a)=(a*(dexp(dmin1(24.0d0,(2.0d0*u)/a))-1.0d0))/
     &(dexp(dmin1(24.0d0,(2.0d0*u)/a))+1.0d0)
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
      if (((((mmod.eq.'DIS').or.(mmod.eq.'LIN')).or.(mmod.eq.'TAU'))
     &.or.(mmod.eq.'DPRI')).or.(mmod.eq.'LPRI')) goto 20
      write (*,10) mmod
   10 format(//'MODE IMPROPERLY DEFINED IN TAUDIST : ',a4)
      stop
   20 continue
c
      fg=0.66d0
      convg=0.005d0
      sc=1.d-34
      geo=1.0d0
      cof=convg*2.d0
c
      if ((mmod.ne.'DIS').and.(mmod.ne.'DPRI')) then
        steep=2.0d0
        expo=1.0d0
        divsi=1.0d0
        ctpos=0.0d0
        cmpos=0.4d0
        radius=0.d0
      else
        steep=150.0d0
        expo=0.7d0
        divsi=8.0d0
        ctpos=0.15d0
        cmpos=0.005d0
        radius=dabs(rad)
      endif
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (mmod.ne.'TAU') drta=0.0d0
      phma=0.0d0
      relma=0.0d0
      simax=1.d-37
      do 60 inl=1,infph-1
c        cebin=cphotev(inl)
        if (skipbin(inl)) goto 60
        sig=0.0d0
        rewe=0.0d0
        do 30 i=1,atypes
   30     rela(i)=0.0d0
        do 40 i=1,ionum
          if (photxsec(i,inl).le.0.d0) goto 40
          ie=atpho(i)
          j=ionpho(i)
          abio=zion(ie)*(popul(j,ie)**expo)
          crosec=photxsec(i,inl)
c          eph=cebin/ipotpho(i)
c          if (eph.lt.1.0d0) goto 50
c          crosec=acrs(eph,sigpho(i),betpho(i),spho(i))
          sig=sig+(abio*crosec)
          if (arad(j+1,ie).gt.0.0d0) rela(ie)=1.0d0
   40   continue
c   50   continue
c
        sig=(dh*fi)*sig
c
        do 50 i=1,atypes
   50     rewe=rewe+(rela(i)*zion(i))
c
        phma=dmax1(phma,tphot(inl))
        simax=dmax1(simax,sig*(rewe**2.d0))
        relma=dmax1(relma,rewe)
        sigen(inl)=sig
        relwei(inl)=rewe
c
   60 continue
c
c
c    ***AVERAGES SIGMA OVER PHOTON SPECTRUM
c
      simax=simax/(relma**2.d0)
      sca=1.0d0/dmax1(1.d-37,dsqrt(phma)*sc)
      sca2=1.0d0/dmax1(sc,(dsqrt(phma)*sca)*sc)
      weima=0.0d0
      weim2=0.0d0
      do 70 inl=1,infph-1
        if (skipbin(inl)) goto 70
        testab=(relwei(inl)/relma)**2.d0
        if (testab.gt.0.001d0) then
          ratio=(photev(inl+1)-photev(inl))/(photev(inl)+photev(inl+1))
          ecr=dsqrt(sigen(inl)/simax)
          wei=((dsqrt(tphot(inl))*sca)*sca2)*dsqrt(ratio)
          relp=wei*testab
          relwei(inl)=relp
          weima=dmax1(weima,ecr*relp)
          weim2=dmax1(weim2,relp)
        else
          relwei(inl)=0.0d0
        endif
   70 continue
c
      scaw=1.0d0/dmax1(1.d-35,1.d-35*weima)
      scaw2=1.0d0/dmax1(1.d-35,1.d-35*weim2)
      do 140 k=1,14
        pos=dmin1(1.d2,cmpos+(ctpos/(tauav+1.d-20)))
        sigtot=0.d0
        sigwei=0.d0
        gammax=0.d0
c
        do 80 inl=1,infph-1
          if (skipbin(inl)) goto 80
          ecr=dsqrt(sigen(inl)/simax)
          ar=drta*sigen(inl)
          gam=far(pos*ar,steep)
          gammax=dmax1(gam,gammax)
          wei=relwei(inl)
          if ((wei.gt.0.0d0).and.(ecr.gt.0.0d0)) then
            adto=dexp(((dlog(wei)+dlog(ecr))+dlog(scaw))-gam)
            weiad=dexp((dlog(wei)+dlog(scaw2))-gam)
          else
            adto=0.0d0
            weiad=0.0d0
          endif
          sigwei=sigwei+weiad
          sigtot=sigtot+adto
   80   continue
c
c    ***DETERMINE AVERAGE TAUNEW
c
        sigwei=dmax1(1.d-37,sigwei)
        if (sigwei.gt.0.d0) then
          sigm=dmax1(1.d-37,simax*(dexp(((dlog(sigtot+1.d-37)-
     &     dlog(sigwei))-dlog(scaw))+dlog(scaw2))**2))/divsi
        else
          sigm=1.d-37
        endif
        taunew=(sigm*geo)*drta
c
c    ***IMPROVES VALUE OF DRTA IF MMOD='DIS','LIN'
c    TAKES INTO ACCOUNT THE AMPLIFICATION FACTOR DUE TO GEOM.
c
        if (mmod.eq.'TAU') goto 130
        lim=5
   90   lim=lim+5
        if (radius.eq.0.d0) lim=2
        drta=tauav/sigm
        do 110 n=1,lim
          dri=drta
          geo=1.0d0
          if (radius.eq.0.d0) goto 100
          wlo=(2.d0*dlog(radius))-(2.d0*dlog(radius+(fg*drta)))
          ulo=(2.d0*dlog(radius+(fg*drta)))-(2.d0*dlog(radius))
          if (rad.lt.0.0d0) geo=dexp(wlo)
          if (rad.gt.0.0d0) geo=dexp(ulo)
  100     continue
c
          u=dble(n)*(1.0d0/(dble(lim)-1.0d0))
          if (n.eq.lim) u=1.0d0
          sigmgeo=sigm*(geo**u)
          drta=tauav/sigmgeo
          cof=dabs(drta-dri)/dmax1(drta,dri,1.d-8)
  110   continue
c
        if ((cof.gt.convg).and.(lim.lt.26)) goto 90
        if ((mmod.eq.'DPRI').or.(mmod.eq.'LPRI')) write (*,120) tauav,
     &   taunew,drta,geo,gammax,lim,mmod
c
  120 format(' ',5(1pg12.4),i4,1x,a4)
        rr=taunew/(tauav+1.d-20)
        if ((dabs(1.0d0-rr).lt.convg).and.(k.gt.2)) goto 150
        goto 140
c
c    ***REPLACES VALUE OF TAUAV USING NEW TAU IF MMOD='TAU'
c
  130   continue
        tauav=taunew
c
c     *NOW EXIT
c
        goto 150
c
c     *LOOP BACK
c
  140 continue
c
  150 continue
c
      return
      end
