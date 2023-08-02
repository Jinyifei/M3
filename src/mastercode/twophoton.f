cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine TWOPHOTON calculates the two-photon continuum
c     contribution to the local diffuse field for H and He-like ions.
c     It adds the emission to CNPHOT (in photons/cm3/s/Hz/sr).
c
c     Must be called soon after HYDRO so that COLLRATE and RECRATE are
c     up to date, i.e., from CONTCONT from LOCALEM.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine twophoton (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh,bev
      real*8 fpsiy,rkt,t12
      real*8 y,dexcoll
      real*8 r2q,ab0,ab1,tz,dp
      real*8 u,phots,omes1
c
      integer*4 atom,nz,inl
c
      real*8 x,z,z2,z12,a,ez,f
      real*8 a21,ara,cos1,frs1
      real*8 omep1,ra3,ra1,ras1
      real*8 reche,es1
c
c     Internal functions for hydrogen.
c
      real*8 qpr,qel,ahi
c
c     Functions for Helium I
c
      real*8 pol5,ahei
      real*8 qss,qsp,frs3,fusp
c
      qpr(u)=4.74d-4*(u**(-0.151d0))
      qel(u)=0.57d-4*(u**(-0.373d0))
c
      pol5(x)=dmax1(0.0d0,(((((((((-(6.5289d0*x))+41.54554d0)*x)-
     &97.135778d0)*x)+97.0517d0)*x)-32.02831d0)*x)-0.045645d0)
c
      qss(u)=1.d-8*pol5(u**0.3333333d0)
      qsp(u)=5.73d-9*pol5((0.175d0+(0.2d0*u))**0.3333333d0)
c
      frs3(u)=0.775d0*(u**0.0213d0)
      fusp(u)=0.42d0*(u**0.23d0)
c
c HI Parpia, F. A., and Johnson, W. R., 1982, Phys. Rev. A, 26, 1142.
c HI Goldman, S.P. and Drake, G.W.F., 1981, Phys Rev A, 24, 183
c
      ahi(z,a)=8.22943d0*(z**6)*(1+(a*z)*(a*z)*(3.9448d0-(a*z)*(a*z)*
     &2.040d0))/(1.d0+4.6019d0*(a*z)*(a*z))
c
c HeI Drake, G.W.F., 1986, Phys. Rev. A, 34, 2871
c
      ahei(z)=16.458762d0*((z-0.806389d0)**6)*(1.d0+(1.539d0/(z+2.5d0)**
     &2))
c
c     get scaled energies wrt Te
c
      rkt=rkb*t
      t12=1.d0/dsqrt(t)
c
c      u   = bev/(rkt)
c
      do inl=1,infph
        p2ph(inl)=0.d0
      enddo
c
      do atom=1,atypes
c
        r2q=0.d0
c
        nz=mapz(atom)
        ab1=zion(atom)*pop(nz+1,atom)
        ab0=zion(atom)*pop(nz  ,atom)
c
        if ((ab0.gt.pzlimit).and.(ab1.gt.pzlimit)) then
c
c   H like atoms
c
          z=dble(mapz(atom))
          z2=z*z
          ez=ee2p(atom)
c
          r2q=de*dh*(ab0*collrate2p(atom)+ab1*recrate2p(atom))
c
c     scale temperature
c
          tz=1.d-4*t/z2
c
c     proton & electron deexcitation collision
c
          dp=dh*zion(1)*pop(2,1) ! protons
          dexcoll=(dp*qpr(tz)+de*qel(tz))/z
c
          a21=ahi(z,alphafs)
          r2q=r2q/(1.0d0+(dexcoll/a21))
c
          if (r2q.gt.epsilon) then
          do inl=1,infph-1
            bev=cphote(inl)
            y=bev/ez
            if (y.lt.1.d0) then
              phots=2.d0*plk*y*fpsiy(y)*r2q/(fpi*bev)
              p2ph(inl)=p2ph(inl)+phots
            endif
          enddo
          endif
        endif
      enddo
c
c   Old  (wrong) He I collision strengths.
c
      omep1=2.57d0
      omes1=0.55d0*omep1
c  energy ratio for He like levels
c      xe=heien(2)/(0.75d0*ipote(1,2))
c
      do atom=2,atypes
c
c      atom=2
      nz=mapz(atom)-1
      ab1=zion(atom)*pop(nz+1,atom)
      ab0=zion(atom)*pop(nz  ,atom)
c
      if ((ab1.gt.pzlimit).and.(ab0.gt.pzlimit)) then
c
        z=dble(nz+1)
        z2=z*z
c
        ez=ee2phe(atom)
        es1=ez/rkt
c
c approx He 2p coll strengths as f(z)
c based on CHIANTI 8 mid spline values
c
        z12=dsqrt(z-1.d0)
        omes1=10.d0**(-9.9690240511d-01+z12*
     &  (-4.2553741502d-01+z12*2.1011474554d-04))
c
        r2q=0.d0
c
c        u=dabs(4.d0*t/z2) ! scaled t
c        f=dsqrt(1.0d0/u)
c
        f=t12
c
        cos1=0.d0
        if (es1.lt.loghuge) then
c
c    ***APPROX COLLISIONAL RATES FROM GROUND STATE 1S
c
          ara=rka*f
C         br0=(1.d0-1.34d0/z)
C         ara=rka*f*br0
C         cos1=ara*(omes1/(nz*nz))*dexp(-es1)*fgaunt(1,1,es1)
          cos1=ara*omes1*dexp(-es1) ! gi = 1.0, br included in omes fit
          if ( collrate2phe(atom).gt.0.d0) then
c get full t dependent rate if available
            cos1=collrate2phe(atom)
          endif
        endif
c
        if (nz.eq.1) then
c
c   Recombination into triplet and singlet systems
c
c
c   He atom only - don't yet have rec(4,2) for heavier elements.
c
          tz=1.d-4*t
          if (tz.lt.0.05d0) tz=0.05d0
          if (tz.gt.10.0d0) tz=10.d0
c
          reche=rec(4,2)
c
          frs1=(reche*fusp(tz))*(1.0d0-frs3(tz))
c
          ra3=reche*frs3(tz)
          ra1=(reche*fusp(tz))*(1.0d0-frs3(tz))
          ras1=ra1*qss(tz)
c
          r2q=de*dh*(ab0*cos1+ab1*ra1)
          dexcoll=de*(qss(tz)+qsp(tz))
          a21=ahei(z)
c
          r2q=r2q/(1.0d0+(dexcoll/a21))
        else
c
c colls only for heavy He atm, will fix when we have
c rec(4,2) equivalent rates
c
          r2q=de*dh*ab0*cos1
        endif
c
c        write(*,*) 'HeR', mapz(atom), nz, cos1*de*dh*ab0,r2q
c
        if (r2q.gt.epsilon) then
          do inl=1,infph-1
            bev=cphote(inl)
            y=bev/ez
            if (y.lt.1.d0) then
              phots=2.d0*plk*y*fpsiy(y)*r2q/(fpi*bev)
              p2ph(inl)=p2ph(inl)+phots
            endif
          enddo
        endif
c
      endif
c
c   He like atoms
c
      enddo
c
      return
c
      end
