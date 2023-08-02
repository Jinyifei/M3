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
c     COMPUTES COLLISIONAL IONISATION RATES FOR
c     ALL IONS ; RATES RETURNED AS VECTOR COL(6,11)
c     IN COMMON BLOCK /COLION/
c
c     This is a new routine to calculate collisional ionisation
c     rates using an algorithm based on Arnaud and Rothenflug 1985
c     and Younger 1981,1983.  It uses a five parameter fit to the
c     collision crossection and derives an expresion for the integral
c     over a maxwellian velocity distrbution in electron energies.
c
c     Numerical procedures sugested by A&R have been revised,
c     and extended to more rigorous precision.
c       The exponential integrals
c     are evaluated in the function farint.  The basis for the second
c     integral comes from Hummer 1983 and earlier references thererin.
c     The Hummer expression is replaced when x is small with a
c       more accurate convergent series.
c
c     RSS 7/90
c
c     added selector for different calculation methods
c     so the user can choose by setting the mode in ATDAT
c
c     RSS 8/90
c
c     EXPERIMENTAL: added kappa distributon rate enhancements shell by shell for
c     A&R collmode = 0 or 1 only
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine collion (telec)
c
      include 'cblocks.inc'
c
      real*8 telec, f,f0, f1, f2, a, b, c, d
      real*8 ktev, ac, tc, rateenhance
      real*8 cdi,cai,fautoi
      real*8 rate,x,y,ep,invt,ex1
      real*8 xs(20),ys(20),y2s(20)
c
      integer*4 ion ,atom, shell,is
      integer*4 ni,io,k
c
      real*8 t,sig
c
c     External Functions
c
      real*8 farint,fsplint,fcodyE1
      real*8 fkenhance
c
      rateenhance=1.d0
c
      do atom=1,atypes
        do ion=1,maxion(atom)
          col(ion,atom)=0.0d0
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New Dere code, A&R or L&M collisional ionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (collmode.eq.0) then
c
c     Dere 2007 rates, using the spline fit for y
c
        f=2.d0
        ktev=telec*rkb/ev
        do io=1,ncol
          rate=0.d0
          atom=atcol(io)
          ion=ioncol(io)
          ep=ipotcol(io)
          invt=ep/ktev
          if (invt.lt.maxdekt) then
            if (usekappa) then
              rateenhance=fkenhance(kappa,invt)
            endif
            t=1.d0/invt
            ni=nsplcol(io)
            do k=1,ni
              xs(k)=xcol(k,io)
              ys(k)=ycol(k,io)
              y2s(k)=y2col(k,io)
            enddo
            ex1=fcodyE1(invt)
            rate=ex1/(dsqrt(t*ep)*ep)
            x=1.d0-(dlog(f)/dlog(t+f))
            y=fsplint(xs,ys,y2s,ni,x)
            rate=1.0d-6*rate*y
            rate=rate*rateenhance
          endif
          col(ion,atom)=rate
C          if (mapz(atom).eq.8) then
C          write(*,*) telec,ion,mapz(atom),col(ion,atom)
C          endif
        enddo
c
      endif
c
      if (collmode.eq.1) then
c
c     Arnaud and Rothenflug rates (previous default)
c
        ktev=telec*rkb/ev
        do 30 atom=1,atypes
          do 20 ion=1,maxion(atom)-1
            sig=0.d0
            do 10 shell=1,nshells(ion,atom)
              a=aar(shell,ion,atom)
              b=bar(shell,ion,atom)
              c=car(shell,ion,atom)
              d=dar(shell,ion,atom)
              x=collpot(shell,ion,atom)/ktev
              if (x.lt.maxdekt) then
                rateenhance=1.0d0
                if (usekappa) then
                  rateenhance=fkenhance(kappa,x)
                endif
                f1=farint(1,x)
                f2=farint(2,x)
                f0=a*(1.d0-x*f1)+b*(1.d0+x-x*(2.d0+x)*f1)+c*f1+d*x*f2
                sig=sig+(dexp(-x)/x)*f0*rateenhance
              endif
   10       continue
            cdi=((6.69d-7)/(ktev**1.5d0))*sig
            cai=fautoi(telec,atom,ion)
            col(ion,atom)=dmax1(cdi+cai,0.d0)
   20     continue
          col(maxion(atom),atom)=0.0d0
   30   continue
c
      endif
c
      if (collmode.eq.2) then
c
c     Landini & Monsignori Fossi - Shull rates
c     faster but probably less accurate.....
c     will use A&R for H sequence.
c
        t=telec+1.d-1
        f=dsqrt(t)
        ktev=t*rkb/ev
        do 50 atom=1,atypes
          do 40 ion=1,maxion(atom)-1
            is=mapz(atom)-ion+1
            cdi=0.d0
            if (is.gt.1) then
              ac=acol(ion,atom)
              tc=tcol(ion,atom)
              cdi=ac*f/(1.d0+0.1d0*(t/tc))*dexp(-tc/t)
            else
              a=aar(1,ion,atom)
              b=bar(1,ion,atom)
              c=car(1,ion,atom)
              d=dar(1,ion,atom)
              x=collpot(1,ion,atom)/ktev
              if (x.lt.maxdekt) then
                rateenhance=1.0d0
                if (usekappa) then
                  rateenhance=fkenhance(kappa,x)
                endif
                f1=farint(1,x)
                f2=farint(2,x)
                f0=a*(1.d0-x*f1)+b*(1.d0+x-x*(2.d0+x)*f1)+c*f1+d*x*f2
                sig=(dexp(-x)/x)*f0*rateenhance
                cdi=((6.69d-7)/(ktev**1.5d0))*sig
              endif
            endif
            col(ion,atom)=cdi
   40     continue
   50   continue
c
      endif
c
      return
      end
