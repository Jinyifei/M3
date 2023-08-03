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
c      rankhug solves the Rankine Hugoniot conditions between
c      any two points with full specification of v0, and
c      loss rate, returns v1.
c
c      sets te0 and te1 :pre and post temps
c      magnetic field Bmag (Gauss), populations in common pop.
c      also need dh0, H num density
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine rankhug (tpr, delpr, dhpr, vpr, bmag, tl, tstep)
c
      include 'cblocks.inc'
c
      real*8 a(5)
      real*8 tpr,delpr,dhpr,vpr,bmag
      real*8 tl, tstep
      real*8 g,lambda, en0,delta
      real*8 rv0,beta,q,pram
      real*8 cmpf,x1,x2,r1,r2,x
      real*8 cmpfNEW,rv2,cmpfHD
c
c      character rmod*4
c
c test vars
c
C      real*8 mom,ener,ent0,ent1
C      real*8 root
C      real*8 cmpfOLD
c       integer*4 i
c       character qmod*4
c       real*8 cmpfSHK,alpha,m1,y,ym,an
c
c      real*8 cmpfSHK
c test fn
c      real*8 shockcmpf
      real*8 frho
c
c     set globals
c
      vel0=vpr
      te0=tpr
      dh0=dhpr
      de0=delpr ! feldens(dhpr,pop)
      bm0=bmag
      pb0=(bm0*bm0)/epi
c
      if (vel0.lt.0.d0) then
        write (*,*) 'ERROR in RankHug: vel0 <= 0'
        write (*,*) vel0
        stop
      endif
c
c     get pressure and density terms
c
      en0=zen*dh0+de0
      pr0=en0*rkb*te0
      rho0=frho(de0,dh0)
c
      rv0=rho0*vel0
      rv2=rv0*vel0
c
c quadratic full solution for non-magnetic case, for shocks
c and flows, with or without cooling
c
      g=(gammaeos)/(gammaeos-1.d0)
      lambda=tl*tstep
      a(3)=-(g*pr0+0.5d0*rv2-lambda)/rv2
      a(2)=g*(pr0+rv2)/rv2
      a(1)=(0.5d0-g)
      delta=(a(2)*a(2))-(4.d0*(a(1)*a(3)))
      q=-0.5d0*(a(2)+dsign(1.d0,a(2))*dsqrt(delta))
      r1=a(1)/q
      r2=q/a(3)
      cmpfHD=dmax1(r1,r2)
      vel1=vel0/cmpfHD
      cmpf=cmpfHD
c
      if (bmag.gt.0.d0) then
C
C old quartic for velocity solution - bug in there somewhere, but
C no longer used so not investigated.  Momentum is conserved but not
C energy
C
C uncomment to run tests below
C     use non-magnetic velocity as initial guess for
C     quartic mag field solution
C
C       lambda=2.d0*tl*tstep/rho0
C       g=(2.d0*gammaeos)/(gammaeos-1.d0)
C       beta=pb0*4.d0
C       a(5)=(g-1.d0)
C       a(4)=-g*(pr0/rv0+vel0+(beta/(4.d0*rv0)))
C       a(3)=((vel0*vel0)+(g*(pr0/rho0))+(beta)-lambda)
C       a(2)=g*(beta*vel0/(4.d0*rho0))
C       a(1)=-1.d0*beta*vel0*vel0
CC
C       x1=vel0*0.25d0
C       x2=vel0
CC
C       call quart (a, x1, x2, root, 'NEAR')
C       vel1=root
C       vel1=dmax1(0.d0,vel1)
C       if (vel1.le.0.d0) vel1=vel0
C       cmpfOLD=vel0/(vel1+epsilon)
c        cmpf=cmpfOLD
Cc
C
C       MHD shock only quadratic solution, discards x=1 root so cant be
C       used in general flow with cooling, here for validaton purposes
C
c        cmpfSHK=shockcmpf (te0, de0, dh0, vel0, bm0)
Cc
        lambda=tl*tstep
        g=(gammaeos)/(gammaeos-1.d0)
        beta=pb0*2.d0
        pram=rv2 ! for clarity and compared to paper
c /rv2 scale to normalise numbers in cubic root finder better
        a(4)= (pb0*(2.d0-g))/rv2
        a(3)=-(g*pr0+0.5d0*pram+2.d0*pb0-lambda)/rv2
        a(2)= g*(pr0+pram+pb0)/rv2
        a(1)= (0.5d0-g)
c
        x1=0.5d0
        x2=(gammaeos+1.d0)/(gammaeos-1.d0)
c
c find the two positive roots near x1 and x2  for the cubic in a
c
        call cubic (a, x1, x2, r1, r2)
        cmpfNEW=dmax1(r1,r2) ! for S5 this is the root we always need
c
        cmpf=cmpfNEW
c
       endif
c
       x=cmpf
       vel1=vel0/x
       bm1=bm0*x
       pb1=(bm1*bm1)/epi
       rho1=rho0*x
       call invrho (rho1, de1, dh1, pop)
c
       pr1=pr0+(rho0*vel0*vel0*(1.d0-(1.d0/x)))+pb0*(1.d0-(x*x))
       te1=pr1/((zen*dh1+de1)*rkb)
C
C      ent0=pr0/(rho0**gammaeos)
C      ent1=pr1/(rho1**gammaeos)
CC
CC  test to confirm conservation laws.  cmpfOLD fails, cmpfNEW is good
CC
C      mom=(rho0*vel0*vel0+pr0+pb0)-rho1*vel1*vel1-pr1-pb1
C      mom=100.d0*mom/(rho0*vel0*vel0+pr0+pb0)
C      g=(gammaeos)/(gammaeos-1.d0)
C      lambda=tl*tstep
C      ener=((g*pr0+0.5d0*rho0*vel0*vel0+2.d0*pb0-lambda)/rho0)
C    &     -((g*pr1+0.5d0*rho1*vel1*vel1+2.d0*pb1)/rho1)
C      ener=100.d0*ener/((g*pr0+0.5d0*rho0*vel0*vel0+2.d0*pb0)/rho0)
CCC
C 100  format('RH: ',
C    & '     v:',1pg14.7,
C    & '   ∆m%:',f10.5,
C    & '   ∆E%:',f10.5,
C    & ' <L>dt:',1pg14.7,
C    & ' x_new:',1pg14.7,
C    & ' x_old:',1pg14.7)
C      write(*,100) vel0,mom,ener,lambda,cmpfNEW,cmpfOLD
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
c
c solve for just the jump compresson factor, no cooling or flow
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function shockcmpf (t, de, dh, v, bm)
c
      include 'cblocks.inc'
c
c solve for just the jump compresson factor, no cooling or flow
c
c Cubic divided out by x-1.0 -> quadratic
c
      real*8 t, de, dh, v, bm, cmp
      real*8 m1,alpha,gam
      real*8 c0,c1,c2,delta
      real*8 q,r1,r2
      real*8 pr,rho,en0,pbm,pram
      real*8 frho
c
      cmp=1.d0
      gam=gammaeos
c
      en0=zen*dh+de
      pr=en0*rkb*t
      pbm=(bm*bm)/(epi)
      rho=frho(de0,dh0)
      pram=rho*v*v
C
C pressure ratios
C     m1=v/dsqrt(gam*pr/rho)
C     an=v/dsqrt(2.d0*pbm/rho)
C     alpha=pbm/pr0
C
C     c2=2.d0*alpha*(2.d0-gam) ! a x^2
C     c1=gam*((gam-1.d0)*m1*m1 + 2.d0*(alpha+1.d0)) ! b x
C     c0=-gam*(gam+1.d0)*m1*m1 ! c
C
C     delta=(c1*c1)-(4.d0*(c2*c0))
C
C      q=1.0d0
C     r1=1.0d0
C     r2=1.0d0
C     if (delta.ge.0.d0) then
C        q=-0.5d0*(c1+dsign(1.d0,c1)*dsqrt(delta))
C       if (q.ne.0.d0) r1=c0/q
C       if (c2.ne.0.d0) r2=q/c2
C     endif
C
C     cmp=dmax1(r1,r2)
C
C     write(*,*) 'shockcmpf',m1,alpha,q,r1,r2,cmp
c
C
C  using pressures
C
      c2=2.d0*pbm*(2.d0-gam) ! a x^2
      c1= (gam-1.d0)*pram + 2.d0*gam*(pbm+pr) ! b x
      c0=-(gam+1.d0)*pram ! c
c
      delta=(c1*c1)-(4.d0*(c2*c0))
c
       q=1.0d0
      r1=1.0d0
      r2=1.0d0
      if (delta.ge.0.d0) then
         q=-0.5d0*(c1+dsign(1.d0,c1)*dsqrt(delta))
        if (q.ne.0.d0) r1=c0/q
        if (c2.ne.0.d0) r2=q/c2
      endif
c
      cmp=dmax1(r1,r2)
c
      write(*,*) 'shockcmpf2',m1,alpha,q,c0,c2,r1,r2,cmp
c
      shockcmpf=cmp
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
c     Special Purpose Newtons Method for solving a
c     well behaved quartic.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine quart (a, x1, x2, root, rmod)
c
      include 'const.inc'
c
      real*8 a(5),x1,x2,root
      real*8 d(5),x,f,df,dx,eps
      integer*4 i,itmax
      character rmod*4
c
c     find -ve curve near root wanted
c
c      write(*,*) 'Quart:', rmod
c      x1=0.d0
c
      if (rmod.eq.'FAR') then
c search for 1st -ve down from x2
        dx=(x1-x2)*0.01d0
        x=x2-dx*100
        do i=1,200
          x=x+dx
          f=a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
c          write(*,*) x,f
          if (f.lt.0.0d0) goto 10
        enddo
      else
c search for -ve up from x1
        dx=(x2-x1)*0.01d0
        x=x1
        do i=1,200
          x=x+dx
          f=a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
c          write(*,*) x,f
          if (f.lt.0.d0) goto 10
        enddo
      endif
c
c
   10 x=x-dx
c
c     get derivative
c
      d(5)=0.d0
      do i=1,4
        d(i)=a(i+1)*dble(i)
      enddo
c
      itmax=20
c
      do i=1,itmax
        f=a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
        df=d(1)+x*(d(2)+x*(d(3)+x*d(4)))
        dx=f/df
        x=x-dx
        eps=dabs(x*1.d-10)
        if (dabs(dx).lt.(eps)) goto 20
      enddo
c
   20 root=x
c
c      write(*,*) 'Quart: Roots: ', i, x1, x2,  x, x-x1
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Special Purpose Newtons Method for solving a
c     well behaved cubic.
c
c     modified rss2009
c     modified rss2010
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cubic (a, x1, x2, r1, r2)
c
c find the two positive roots near x1 and x2  for the cubic in a
c
c Solve characterisitc cubic for R-H flow with cooling
c !!!! modifies a, by normalising.
c returns both roots if possible, may the the same.
c
c
      include 'const.inc'
c
      real*8 a(5),x1,x2,r1,r2,invk
      real*8 d(5),f0,f1,f2,df,dx,eps,mxx,mnx
      integer*4 i,itmax
c
c      write(*,*) 'Cubic:', rmod
c
c find the value at x=0 for normalising the
c polynomial.
c
      invk=dabs(1.d0/a(1))  ! preserve sign of poly
      do i =1,4
         a(i)=invk*a(i)
      enddo
c
c     get derivative
c
      d(4)=0.d0
      do i=1,3
        d(i)=a(i+1)*dble(i)
      enddo
c
      mxx=dmax1(x1,x2)*1.1d0
      mnx=dmin1(x1,x2)*0.9d0
c
      dx=(mxx-mnx)*0.01d0
C         write(*,*) mxx,mnx,dx
c search for 1st sign change down from mxx
        f0=a(1)+mxx*(a(2)+mxx*(a(3)+mxx*a(4)))
        do i=0,200
          x2=mxx-dx*dble(i)
          f2=a(1)+x2*(a(2)+x2*(a(3)+x2*a(4)))
c          write(*,*) x2,f2
          if (f2.ne.dsign(f2,f0)) goto 10
          f0=f2
        enddo
   10  continue
c search for 1st sign change up from mnx
        f0=a(1)+mnx*(a(2)+mnx*(a(3)+mnx*a(4)))
        do i=0,200
          x1=mnx+dx*dble(i)
          f1=a(1)+x1*(a(2)+x1*(a(3)+x1*a(4)))
c          write(*,*) x1,f1
          if (f1.ne.dsign(f1,f0)) goto 20
          f0=f1
        enddo
   20 continue
c
      itmax=20
c
      do i=1,itmax
        f1=a(1)+x1*(a(2)+x1*(a(3)+x1*a(4)))
        df=d(1)+x1*(d(2)+x1*d(3))
        dx=f1/df
        x1=x1-dx
        eps=dabs(x1*1.d-10)
        if (dabs(dx).lt.(eps)) goto 30
      enddo
c
   30 r1=x1
c
      do i=1,itmax
        f2=a(1)+x2*(a(2)+x2*(a(3)+x2*a(4)))
        df=d(1)+x2*(d(2)+x2*d(3))
        dx=f2/df
        x2=x2-dx
        eps=dabs(x2*1.d-10)
        if (dabs(dx).lt.(eps)) goto 40
      enddo
c
   40 r2=x2
c
c      write(*,*) ' Cubic: Roots: ', i, r1, r2
c
      return
c
      end
c
      double precision function velshock2 (dh, t, bmag, tpo)
c
c estimate shock velocity to 1 in 10^7 for given shock temperature.
c using bisection.  Mindlessly solid in the face of
c terrible instability at high B!  Performance is not an issue, this
c routine is only used in setup, not in running  models.
c
c **** assumes tpo is a monotonic function of vin ****
c
      include 'cblocks.inc'
c
      real*8 t,de,dh,bmag,tpo,tol
      real*8 tinit,vin,binit,dv,v0,v1,v2,eps
      integer*4 idx,itdivs,maxits
      real*8 feldens
c save local copy in case it is a ref
      tinit=tpo
      binit=bmag
c         write(*,*) 'velshock2:',bmag*1e6,tpo
c
      maxits=100
      tol=1.0d-7
c
      de=feldens(dh,pop)
      te0=t
      te1=tpo
      dh0=dh
      de0=de
c
c   simply step in 100km/s then 10km/s then 1km/s then 0.1...
c
      dv=100.d5
c
c enough for 30,000 km/s
c find first 100km/s too fast
c
      do idx=1,300
         vin=dble(idx)*dv
         call rankhug (t, de, dh, vin, binit, 0.d0, 0.d0)
         write(*,*) 'init',idx,vin*1.d-5,te1,tinit
         if ((te1+epsilon).gt.tinit) goto 10
      enddo
c Make sure people going over 30,000km/s know what they are doing...
      write(*,*) 'ERROR: shock too fast in velshock2, tepo too large'
      stop
c
   10  continue
c
c polish with bisection
c
      itdivs=0
      v0=vin-dv
      v1=vin
   20 continue
      v2=0.5d0*(v0+v1)
      call rankhug (t, de, dh, v2, binit, 0.d0, 0.d0)
      if (te1.gt.tinit) then
          v1=v2
      else
          v0=v2
      endif
      itdivs=itdivs+1
      eps=dabs(2.d0*(v1-v0)/(v1+v0))
c        write(*,*) itdivs,v0,v1,v2,vel0,te1,tinit
      if (itdivs.gt.maxits)then
        write(*,*) 'ERROR: in velshock2, to many iterations'
        write(*,*) itdivs,v0,v1,v2,vel0,te1,tinit
        stop
      endif
      if (eps.gt.tol) goto 20
c
      velshock2=vel0
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******COMPUTES SHOCK FRONT VELOCITY USING PRESHOCK DENSITY:DHPR
c     PRESHOCK AND POSTSHOCK TEMPERATURES : TEPR,TEPO
c     MAGNETIC FIELD : HMAG  , IONIC POPULATIONS : POP
c     AND FRACTIONAL IONISATION OF HYDROGEN : XHPR
c     NB.   ALL THESE QUANTITIES NEED TO BE DEFINED
c
c     OUTPUT INVARIANTS AND VELOCITIES IN COMMON BLOCK /VELSHO/
c
c
      subroutine velshock (dhpr, xhpr, tepr, tepo, hmag)
c
      include 'cblocks.inc'
c
      real*8 dhpr, xhpr, tepr, tepo, hmag
c
c           Variables
c
      real*8 aa(10), a1(10)
      real*8 mu, muf
      real*8 al,alp,b1,b2,b3,be,bet,co1
      real*8 co2,co3,eps,ga,gam,gamf,humag,rho
      real*8 root,tmin,u1,u2,v1
c
      integer*4 ncoef
c
c           Functions
c
      real*8 feldens,fmua,frho
c
      write (*,*) 'Entering velshock'
c
      if (xhpr.le.0.d0) xhpr=pop(2,1)
      if (pop(2,1).ne.xhpr) pop(2,1)=xhpr
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      if ((pop(1,1)+xhpr).ne.1.0d0) pop(1,1)=1.0d0-xhpr
      depr=feldens(dhpr,pop)
      mu=fmua(depr,dhpr)
      rho=frho(depr,dhpr)
c  magnetic energy dens. u=B^2/8pi
      humag=(hmag*hmag)/epi
      hmago=humag/rho
c
      alp=((5.d0*rgas)*(tepo-tepr))/mu
      bet=(rgas*tepr)/mu
      gam=(rgas*tepo)/mu
c
      write (*,*) 'velshock rgas mu alp bet gam',rgas,mu,alp,bet,gam
c
      co1=alp-(2.0d0*(gam-bet))
      al=alp/co1
      be=bet/co1
      ga=gam/co1
c
      co2=-((2.0d0*al*be)+((ga-al)**2)-(be**2))
      co3=-((al*be)**2)
      co2=co2*co1
      co3=co3*co1
      co1=1.0d0
      u1=((-co2)+dsqrt((co2*co2)-((4.d0*co1)*co3)))/(2.d0*co1)
      vshoc=dsqrt(u1)
c
c     ***APPROXIMATE SHOCK VELOCITY COMPUTED : VSHOC
c     NOW TAKE MAGNETIC FIELD INTO ACCOUNT IN ITERATIVE FASHION
c
   10 u1=vshoc*vshoc
      b1=bet/u1
      b2=hmag/u1
      b3=0.0d0
      aa(1)=1.0d0
      aa(2)=-((5.0d0*((b1+b2)+1.0d0))/4.0d0)
      aa(3)=((((5.0d0*b1)/4.0d0)+b2)+0.25d0)-(b3/2.0d0)
      aa(4)=b2
      ncoef=4
      root=1.0d0
c
      call roots (aa, ncoef, root)
c
      a1(1)=1.0d0
      a1(2)=(aa(1)*root)+aa(2)
      a1(3)=(((aa(1)*root)*root)+(aa(2)*root))+aa(3)
      ncoef=3
      root=0.25d0
c
      call roots (a1, ncoef, root)
c
      vpo=root*vshoc
      u2=vpo*vpo
      v1=vshoc
      vshoc=0.5d0*dsqrt(5.d0*(gam-bet)+(4.0d0*hmag)*((vshoc/vpo)-1.0d0)+
     &u2)+(0.5d0*vshoc)
      eps=dabs((vshoc-v1)/vshoc)
      if (eps.gt.0.001d0) goto 10
c
c     ***PHYSICAL CONDITIONS , POST-SHOCK
c
      dhpo=(dhpr*vshoc)/vpo
      depo=feldens(dhpo,pop)
      mu=fmua(depo,dhpo)
      rho=frho(depo,dhpo)
      hmago=(hmago*vshoc)/vpo
      gam=(rgas*tepo)/mu
      ww=((5.d0*gam)+(vpo*vpo)+(4.d0*hmago))*0.5d0
      pres=rho*(gam+(vpo*vpo)+hmago)
      fm=rho*vpo
c
c the next expression always evals to 0.d0 for unknown reasons
c it must be a bug and needs checking
c
      qst=(2.d0*ww)-((5.d0*gam)+(4.d0*hmago)+(vpo*vpo))
c
      write(*,*) 'qst ww',qst,ww
c
c     ***COMPUTES APPROXIMATE FINAL STATE OF GAS
c
      tmin=100.0d0
      muf=fmua(0.0d0,1.0d0)
      gamf=(rgas*tmin)/muf
      vfin=((tmin*mu)*vshoc)/(tepo*muf)
c
      qfin=(2.0d0*ww)-(((5.0d0*gamf)+(vfin*vfin))+(((4.0d0*hmago)*vshoc)
     &/vfin))
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******TO EXTRACT A ROOT OF A POLYNOMIAL
c     USING BERNOUILLI'S METHOD
c
c     A IS A VECTOR OF COEFF. IN ORDER OF
c     DECREASING POWER ; TOTAL*OF COEFF.=NCOEF (<9)
c     ROOT=SOLUTION ROOT
c     FEED A GUESS FOR ROOT FIRST
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine roots (a, ncoef, root)
c
      include 'const.inc'
c
      real*8 root, unew,anew,roota,fract
      real*8 a(10), u(10)
      integer*4 ncoef,j,iter,m,mcoef
c
c      write(*,*) 'ROOTS: init guess:',root
c
      u(1)=1.d0
      u(2)=1.d0/root
      do 10 j=3,ncoef
   10   u(j)=0.d0
      iter=0
   20 unew=0.d0
      iter=iter+1
      mcoef=ncoef-1
      do 30 j=1,mcoef
   30   unew=unew-(a(j+1)*u(j))
      anew=dabs(unew)
      if (anew.lt.1.0d20) goto 50
c
c    ***RENORMALISATION OF COEFF. IF TOO LARGE
c
      do 40 j=1,ncoef
   40   u(j)=u(j)/unew
      unew=1.d0
   50 continue
c
      do 60 j=2,ncoef
        m=(ncoef-j)+1
   60   u(m+1)=u(m)
c
c    ***EXTRACTION OF ROOT ; TEST FOR CONVERGENCE
c
      u(1)=unew
      roota=root
      root=u(1)/u(2)
      fract=dabs((root-roota)/root)
      if (iter.gt.150) return
      if (fract.gt.1.0d-7) goto 20
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine isochorflow (tpr, delpr, dhpr, vpr, hmag, tl,tstep)
c
      include 'cblocks.inc'
c
      real*8 tpr,delpr,dhpr,vpr,hmag,humag
      real*8 tl, tstep
      real*8 lambda, en0
      real*8 cmpf
      real*8 u0,u1
c      real*8 x,y
c      integer*4 i
c
c     set globals
c
      vel0=vpr
      te0=tpr
      dh0=dhpr
      de0=delpr ! feldens(dhpr,pop)
      bm0=hmag
c
      en0=zen*dh0+de0
      pr0=en0*rkb*te0
      humag=(hmag*hmag)*iepi
      u0=(1.d0/(gammaeos-1.d0))*en0*rkb*te0
      lambda=tl*tstep
      u1=dmax1(epsilon,(u0-lambda))
      cmpf=u1/u0
      te1=te0*cmpf
      vel1=vel0*cmpf!constantmachnumber
      rho1=rho0
      dh1=dh0
      de1=de0
      en0=zen*dh1+de1
      pr1=en0*rkb*te1
c
      bm0=hmag
      bm1=hmag
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine isobarflow (tpr, delpr, dhpr, vpr, hmag, tl,tstep)
c
      include 'cblocks.inc'
c
      real*8 tpr,delpr,dhpr,vpr,hmag,humag
      real*8 tl, tstep
      real*8 lambda, en0
      real*8 cmpf
      real*8 u0,u1
c      real*8 x,y
c      integer*4 i
c
c     set globals
c
      vel0=vpr
      te0=tpr
      dh0=dhpr
      de0=delpr ! feldens(dhpr,pop)
      bm0=hmag
c
      en0=zen*dh0+de0
      pr0=en0*rkb*te0
      humag=(hmag*hmag)*iepi
      u0=(1.d0/(gammaeos-1.d0))*en0*rkb*te0
      lambda=tl*tstep
      u1=dmax1(epsilon,(u0-lambda))
      cmpf=u1/u0
      te1=te0*cmpf
      vel1=vel0*cmpf!constantmachnumber
      rho1=rho0/cmpf
      dh1=dh0/cmpf
      de1=de0/cmpf
      en0=zen*dh1+de1
      pr1=en0*rkb*te1
c
      bm0=hmag
      bm1=hmag
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
