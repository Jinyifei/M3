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
c     Subroutine freefree calculates the freefree continuum contribution
c     to the local diffuse field..
c     adds the emission to ffph(inl) (in photons/cm3/s/Hz/sr)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine freefree (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh
c
      real*8 zpop,abde,ffconst,invrkt
      real*8 u, gffm, g2, zn2,t12,rkt
      real*8 phots
      integer*4 atom,ion,i,j,inl
      real*8 zsqd,xpf,meanq,energ,et
      real*8 ee(mxinfph)
      real*8 lgg2(mxion)
c
      real*8 freem
c
c      real*8 fgffspline2
c      real*8 fgfflin,fgfflogpoly,fgfflog
      real*8 fgfflog
c
c     internal functions
c
      freem(rkt,energ)=dexp(((-106.366d0-dlog(energ))-(0.5d0*dlog(rkt)))
     &-(energ/rkt))
c
c     external functions
c
c      integer*4 locate
      real*8 favcha,densnum
c
c
      t12=1.d0/dsqrt(t)
      rkt=(rkb*t)
      invrkt=1.d0/rkt
c
      if (freefreemode.eq.0) then
cc
        do inl=1,infph
          ffph(inl)=0.d0
        enddo
c
c basic free free emission for all but specialised primordial
c  or extreme photoionisation cases
c
        meanq=favcha(pop,2)
        zsqd=(de*densnum(dh))*meanq*meanq
c
        do inl=1,infph-1
c
c     get scaled energy of bin wrt te
c
          energ=cphote(inl)
          et=energ*invrkt
c
c prevent over/underflow for exp -et
c
          if (dabs(et).lt.loghuge) then
            xpf=zsqd*freem(rkt,energ)
            ffph(inl)=ffph(inl)+xpf
          endif
        enddo
c
      else
c
        do inl=1,infph
          ffph(inl)=0.d0
          energ=cphote(inl)
          et=energ*invrkt
          ee(inl)=dexp(-et)/energ
        enddo
c
c Reordered loops to maximise cache hits in fgfflog
c changed 100% cache misses to ~10% cache misses
c
        do ion=2,mxion
          zn2=dble(ion-1)
          zn2=zn2*zn2
          lgg2(ion)=dlog10(iphe*zn2*invrkt)
        enddo
c
        do inl=1,infph-1
          energ=cphote(inl)
          et=energ*invrkt
          if (dabs(et).lt.loghuge) then
            u=dlog10(et)
            j=idnint(((u+4.d0)*10.d0))+1
            j=min(max(j,1),ngffu-1)
            do ion=2,mxion
              zn2=dble(ion-1)
              zn2=zn2*zn2
              g2=lgg2(ion)
              g2=dmax1(-4.0d0,g2)
              i=idnint(((g2+4.d0)*5.d0))+1
              i=min(max(i,1),ngffg2-1)
              do atom=1,atypes
                if (ion.le.maxion(atom)) then
                  zpop=zion(atom)*pop(ion,atom)
                  if (zpop.ge.pzlimit) then
                    abde=zpop*de*dh
                    ffconst=ffk*abde*zn2*t12
                    gffm=fgfflog(1,g2,u,i,j)
                    phots=ffconst*gffm*ee(inl)
                    ffph(inl)=ffph(inl)+phots
                  endif
                endif
              enddo
            enddo
          endif
        enddo
      endif
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
      integer*4 function locate(x, np, y)
c
      implicit none
c
      real*8 x,y(np)
      integer*4 np
      integer*4 k,khi,klo
c
c     bisection search
c
      klo=1
      khi=np
   10 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (y(k).gt.x) then
          khi=k
        else
          klo=k
        endif
        goto 10
      endif
c
      locate=klo
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******TO OBTAIN FREE-FREE COOLING RATE ; FFLOS (ERG.CM-3.S-1)
c
c     **REFERENCES : ALLEN,C.W. 1973 AP.QTIES. ED.3,103
c     **REFERENCES : Sutherland 1997
c     **REFERENCES : Karzas & Latter 1961
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine frefre (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      real*8 zpop,abde,abdftz2
      real*8 zn2,t12,gffint
      real*8 lt,lr,g2,invrkt
      integer*4 ion,atom
c
      real*8 fsplint
c
      fflos=0.d0
c
      t12=dsqrt(t)
      invrkt=1.d0/(rkb*t)
      gffint=0.d0
c
      do atom=1,atypes
        do ion=2,maxion(atom)
c
          zpop=zion(atom)*pop(ion,atom)
          abde=zpop*de*dh
c
          if (zpop.ge.pzlimit) then
c
            zn2=dble(ion-1)*dble(ion-1)
            g2=ipote(ion-1,atom)*invrkt
c
            abdftz2=abde*zn2*t12
c
c     Evaluate Spline at log(g2)
c
            lt=dlog10(g2)
            lr=fsplint(gffintg2,gffinty,gffinty2,ngffint,lt)
            lr=fk*lr*abdftz2
c
            gffint=gffint+lr
            coolz(atom)=coolz(atom)+lr
            coolzion(ion,atom)=coolzion(ion,atom)+lr
c
          endif
c
        enddo
      enddo
c
      fflos=gffint
c
      if (fflos.lt.epsilon) fflos=0.d0
c
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
c   interpolates gff from table using 2D cubics in lin gff space.
c   based on Numerical Recipes II
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fgfflin(m,g2,u,il,jl)
c
      include 'cblocks.inc'
c
      integer*4 m,n,mmax,nmax
      real*8 g2,u,y,dy
      real*8 x1,x2
      parameter (nmax=20,mmax=20)
c
      integer*4 i,j,il,jl,kg2,ku
cU    USES polint
      real*8 ymtmp(mmax),yntmp(mmax)
      real*8 x1tmp(mmax),x2tmp(mmax)
c
      fgfflin=1.0
c
      kg2=min(max(il-(m-1)/2,1),ngffg2+1-m)
      ku=min(max(jl-(m-1)/2,1),ngffu+1-m)
c
      n=m-1
c
      x1=g2
      x2=u
c
      do i=kg2,kg2+n
        x1tmp(i-kg2+1)=gffg2(i)
      enddo
      do j=ku,ku+n
        x2tmp(j-ku+1)=gffu(j)
      enddo
c
      do i=kg2,kg2+n
        do j=ku,ku+n
          yntmp(j-ku+1)=10.d0**gff(i,j)
        enddo
        call polint (x2tmp, yntmp, m, x2, ymtmp(i-kg2+1), dy)
      enddo
c
      call polint (x1tmp, ymtmp, m, x1, y, dy)
c
      fgfflin=y
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
      real*8 function fgfflog(m,g2,u,il,jl)
c
      include 'cblocks.inc'
c
      integer*4 m,il,jl
      integer*4 i0,i1,j0,j1
      real*8 g2,u
      real*8 f1,f2
      real*8 gint
      real*8 y,g0,u0
c      integer*4 hitcount,misscount,totalcount
      save y,g0,u0
c      save hitcount,misscount, totalcount
c
      if (m.lt.1) then
        y=0.d0
        g0=0.d0
        u0=0.d0
        fgfflog=0.d0
c          hitcount=0
c          misscount=0
c          totalcount=0
        return
      endif
c      totalcount=totalcount+1
      if ((g0.eq.g2).and.(u0.eq.u)) then
        fgfflog=y
c          hitcount=hitcount+1
        return
      else
c          misscount=misscount+1
        g0=g2
        u0=u
c
c no bounds checking - assumes called properly
c
        i0=il
        i1=i0+1
        j0=jl
        j1=j0+1
        f1=(g2-gffg2(i0))*5.d0
        f2=(u-gffu(j0))*10.d0
        gint=(gff(i0,j0)*(1.d0-f1)*(1.d0-f2))
        gint=gint+(gff(i1,j0)*(f1)*(1.d0-f2))
        gint=gint+(gff(i0,j1)*(1.d0-f1)*(f2))
        gint=gint+(gff(i1,j1)*(f1)*(f2))
        y=10.d0**gint
      endif
      fgfflog=y
c      write(*,*) dble(hitcount)/dble(totalcount)
c
      return
      end
      real*8 function fgfflogpoly(m,g2,u,il,jl)
c
      include 'cblocks.inc'
c
      integer*4 m,n,mmax,nmax
      real*8 g2,u,y,dy
      real*8 x1,x2
      parameter (nmax=20,mmax=20)
c
      integer*4 i,j,il,jl,kg2,ku
cU    USES polint
      real*8 ymtmp(mmax),yntmp(mmax)
      real*8 x1tmp(mmax),x2tmp(mmax)
c
      fgfflogpoly=1.0
c
      kg2=min(max(il-(m-1)/2,1),ngffg2+1-m)
      ku=min(max(jl-(m-1)/2,1),ngffu+1-m)
c
      n=m-1
c
      x1=g2
      x2=u
c
      do i=kg2,kg2+n
        x1tmp(i-kg2+1)=gffg2(i)
      enddo
      do j=ku,ku+n
        x2tmp(j-ku+1)=gffu(j)
      enddo
c
      do i=kg2,kg2+n
        do j=ku,ku+n
          yntmp(j-ku+1)=gff(i,j)
        enddo
        call polint (x2tmp, yntmp, m, x2, ymtmp(i-kg2+1), dy)
      enddo
c
      call polint (x1tmp, ymtmp, m, x1, y, dy)
c
      y=10.d0**y
      fgfflogpoly=y
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NR polynomial interpolation..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine polint (xa, ya, n, x, y, dy)
      integer n,nmax
      real*8 dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 10 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
   10 continue
      y=ya(ns)
      ns=ns-1
      do 30 m=1,n-1
        do 20 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den.eq.0.d0) then
            write (*,*) 'failure in polint'
            stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
   20   continue
        if (2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
   30 continue
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
      real*8 function fgffspline2(lg2, lu)
c
      include 'cblocks.inc'
c
       real*8 lg2,lu
       real*8 y
      call splin2 (gffg2, gffu, gff, gffy2, ngffg2, ngffu, lg2, lu, y)
      fgffspline2=10.d0**y
c
      return
      end
