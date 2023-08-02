cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1979-2012 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c
c     Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes HI and HeII cooling by collisional excitation.
c     Computes HI and HeII recombination spectrum.
c     Uses intermediate case A - B for all calculations.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hydro (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh,f,omg,rate,aa, collloss
      integer*4 i,j,k,l,m,nz,clevel,line,series,atom,ion
      integer*4 nq,ni,nd,idx
c
      real*8 a21,a,tmin,tmax,telec,tz,z,z2,z3
      real*8 tlo,delo,br
      real*8 logrh0,logrh1,logth0,logth1
      real*8 cab,uomg,ez
      real*8 trecom,ratekappa
      real*8 tkappaex,skappa
      real*8 d1b,d2b,d3b
      real*8 d1a,d2a,d3a
      real*8 hb1a,hb2a,hb3a,hb4a
      real*8 egj,helos,dexcoll
      real*8 ab,ab0,ab1,abde,r2q
      real*8 fd,cfd
      real*8 ft,cft
      real*8 fa,cfa
      real*8 y(12),y2(12)
c
c      real*8 om(20)
c
      real*8 qpr,qel,fr,u,ahi
      real*8 fkenhance, fsplint
c
c     internal functions
c
      qpr(u)=4.74d-4*(u**(-0.151d0))
      qel(u)=0.57d-4*(u**(-0.373d0))
      fr(u) =0.588d0*(u**(-0.234d0)) ! fraction of collisions -> 2Phots
c
c HI Parpia, F. A., and Johnson, W. R., 1982, Phys. Rev. A, 26, 1142.
c HI Goldman, S.P. and Drake, G.W.F., 1981, Phys Rev A, 24, 183
c
      ahi(z,a)=8.22943d0*(z**6)*
     &         (1.d0+(a*z)*(a*z)*(3.9448d0-(a*z)*(a*z)*2.040d0))/
     &         (1.d0+4.6019d0*(a*z)*(a*z))
c
c     init
c
      hloss=0.0d0
c
      if (hhecollmode.eq.0) then
        do atom=1,atypes
          do j=1,16
            if ((ishheion(atom).eq.0).or.(j.gt.5)) then
              rateton(atom,j)=0.0d0
            endif
          enddo
          recrate2p(atom)=0.0d0
          if (ishheion(atom).eq.0) collrate2p(atom)=0.0d0
        enddo
      else
        do atom=1,atypes
          do j=1,16
            rateton(atom,j)=0.0d0
          enddo
          recrate2p(atom)=0.0d0
          collrate2p(atom)=0.0d0
          collrate2phe(atom)=0.0d0
        enddo
      endif
      f=dsqrt(1.0d0/(t+epsilon))
      trecom=t
      if (usekappa) then
        trecom=t*(kappa-1.5d0)/kappa
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) ' H Case A-B (0-1) in hydro2:', caseab(1)
c      write(*,*) ' He Case A-B (0-1) in hydro2:', caseab(2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get total collisional excitation up to n=16 levels H and HeII
c     etc, first levels are split by components
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      do atom=1,atypes
      do atom=1,2
c
        nz=mapz(atom)
        ab0=zion(atom)*pop(nz,atom)
c        ab1=zion(atom)*pop(nz+1,atom)
        abde=ab0*de*dh
c
        if (ab0.gt.pzlimit) then
c
        z2=dble(nz*nz)
c
        u=dabs(t)/z2
        uomg=u
c
        if (uomg.ge.5.0d5) uomg=5.0d5
c
        do j=1,20
c
c     ground level excitation
c     energy of upper level
c
          nq=colid(j,3)
          egj=ipote(nz,atom)*(1.d0-(1.d0/(nq*nq)))
          aa=egj/(rkb*t)
c
          if (aa.lt.maxdekt) then
c
            ratekappa=1.0d0
            skappa=1.d0
            tkappaex=t
c
            if (usekappa) then
              skappa=aa
              ratekappa=fkenhance(kappa,skappa)
            endif
c
            if (u.lt.7.2d4) then
              omg=ccoln(j,1)+uomg*(ccoln(j,2)+uomg*(ccoln(j,3)+uomg*
     &         ccoln(j,4)))
            else
              omg=dcoln(j,1)+uomg*(dcoln(j,2)+uomg*(dcoln(j,3)+uomg*
     &         dcoln(j,4)))
            endif
c
c ground state 2S_1/2, gi = 2, wi/w0 = 0.5
c
            rate=ratekappa*rka*f*dexp(-aa)*(omg/z2)*0.5d0
c
            collloss=abde*rate*egj
c
            if (hhecollmode.eq.0) then
c
c Anderson rates for n = 2-5, GNP for 6-16 for H He
c  GNP 2-16 for > He
c
              if (ishheion(atom).eq.0) then
c
c dont have full ion data, use GNP
c
                if (j.eq.1) then
                  collrate2p(atom)=rate
                  rate=0.d0
                  collloss=0.d0
                endif
                hloss=hloss+collloss
                coolz(atom)=coolz(atom)+collloss
                coolzion(nz,atom)=coolzion(nz,atom)+collloss
                rateton(atom,nq)=rateton(atom,nq)+rate*abde ! no eij
              endif

              if ((ishheion(atom).eq.1).and.(nq.gt.5)) then
c have full ion data up to n = 5, use GNP 6-16
                hloss=hloss+collloss
                coolz(atom)=coolz(atom)+collloss
                coolzion(nz,atom)=coolzion(nz,atom)+collloss
                rateton(atom,nq)=rateton(atom,nq)+rate*abde ! no eij
              endif

            else
c All use GNP for 2-16
              hloss=hloss+collloss
              coolz(atom)=coolz(atom)+collloss
              coolzion(nz,atom)=coolzion(nz,atom)+collloss
              if (j.eq.1) then
                collrate2p(atom)=rate
                rate=0.d0
              endif
              rateton(atom,nq)=rateton(atom,nq)+rate*abde

            endif
c
          endif
c
        enddo
        endif
      enddo
cc
c       sumr5=rateton(1,2)
c       sumr5=sumr5+rateton(1,3)
c       sumr5=sumr5+rateton(1,4)
c       sumr5=sumr5+rateton(1,5)
cc
cc       write(*,'(" Kappa T De DH ",4(1pg14.7,x))'), kappa, t, de, dh
c        write(*,'("H N=2 rate:",2(1pg14.7))') rateton(1,2), rateton(2,2)
c        write(*,'("H N=3 rate:",2(1pg14.7))') rateton(1,3), rateton(2,3)
c        write(*,'("H N=4 rate:",2(1pg14.7))') rateton(1,4), rateton(2,4)
c       write(*,'("H N=5 rate:",2(1pg14.7))') rateton(1,5), rateton(2,5)
c       write(*,'("H N=6 rate:",2(1pg14.7))') rateton(1,6), rateton(2,6)
c       write(*,'("H N=7 rate:",2(1pg14.7))') rateton(1,7), rateton(2,7)
c       write(*,'("H N=8 rate:",2(1pg14.7))') rateton(1,8), rateton(2,8)
c       write(*,'("H N=9 rate:",2(1pg14.7))') rateton(1,9), rateton(2,9)
c       write(*,'("H N=10 rate:",2(1pg14.7))')rateton(1,10),rateton(2,10)
c
        hloss=hloss+hheloss
c
c       write(*,*) 'old:',sumr5
c       write(*,*) 'old 2s',collrate2p(1),collrate2p(2)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     ***calculates HI and HeII recomb. lines
c
c     based on Hummer and Storey 1995
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     look up densities
c
      nd=9
c
      i=1
      j=nd
c
      telec=trecom+epsilon
      tlo=dlog10(telec)
      delo=dlog10(de+epsilon)
c
c      write(*,*)' hydro.f line 166: ', telec, de, tlo, delo
c
   10 k=(i+j)/2
c
      if (delo.ge.rhhe(k)) i=k
      if (delo.lt.rhhe(k)) j=k
      if ((j-i).ge.2) goto 10
c
      if (j.eq.1) j=2
      if (i.eq.nd) i=nd-1
      if (i.eq.j) i=j-1
c
      logrh0=rhhe(i)
      logrh1=rhhe(j)
      fd=(delo-logrh0)/(logrh1-logrh0)
      cfd=1.d0-fd
c
C     logrh1=rhhe(j)-logrh0
C     logrh0=(delo-logrh0)/(logrh1+epsilon)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       write(*,*)' hydro.f line 182: ', i,j,delo,rhhe(i),
c     &  rhhe(j),logrh0,logrh1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check temperatures
c
      tmin=10.d0
      tmax=1.0d6
c
      hbeta=0.d0
      halpha=0.d0
      do series=1,nhseries
        do line=1,nhlines
          hydrobri(line,series)=0.d0
        enddo
      enddo
c
      heiilr=0.d0
      do series=1,nheseries
        do line=1,nhelines
          helibri(line,series)=0.d0
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Hydrogen
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ((telec.ge.tmin).and.(telec.le.tmax)) then
c
c     ***FINDS ABSOLUTE FLUX FOR H-BETA
c
        abde=((de*dh)*zion(1))*pop(2,1)
c
        fa=caseab(1)
        cfa=1.d0-fa
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Look up H-Beta temperature
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        ni=12
        l=1
        m=ni
c
   20   if (m-l.gt.1) then
        k=(m+l)/2
        if (th42(k).gt.tlo) then
          m=k
        else
          l=k
        endif
        goto 20
        endif
c
        logth0=th42(l)
        logth1=th42(m)
c
c        logth0=(tlo-logth0)/(logth1-logth0+epsilon)
c
        ft=(tlo-logth0)/(logth1-logth0)
        cft=1.d0-ft
c c
c c     Interpolate Hbeta
c c
c         d1b=(cft*h42b(l,i))+(h42b(m,i)*ft)
c         d2b=(cft*h42b(l,j))+(h42b(m,j)*ft)
c         d3b=(cfd*d1b)+(d2b*fd)
c c
c         d1a=(cft*h42a(l,i))+(h42a(m,i)*ft)
c         d2a=(cft*h42a(l,j))+(h42a(m,j)*ft)
c         d3a=(cfd*d1a)+(d2a*fd)
c c
c c        write(*,*)  d1a, d2a, d1b, d2b
c c        write(*,*)  d3b, d3a, (fa*d3b+cfa*d3a)
c c
c c        write(*,*) logth0,tlo,logth1,ft,cft
c c
        do idx=1,ni
          y (idx)=h42b (idx,i)
          y2(idx)=h42b2(idx,i)
        enddo
        d1b=fsplint(th42,y,y2,ni,tlo)
        do idx=1,ni
          y (idx)=h42b (idx,j)
          y2(idx)=h42b2(idx,j)
        enddo
        d2b=fsplint(th42,y,y2,ni,tlo)
c
        do idx=1,ni
          y (idx)=h42a (idx,i)
          y2(idx)=h42a2(idx,i)
        enddo
        d1a=fsplint(th42,y,y2,ni,tlo)
        do idx=1,ni
          y (idx)=h42a (idx,j)
          y2(idx)=h42a2(idx,j)
        enddo
        d2a=fsplint(th42,y,y2,ni,tlo)
c
c        write(*,*)  d1a, d2a, d1b, d2b
c
        d3b=(cfd*d1b)+(d2b*fd)
        d3a=(cfd*d1a)+(d2a*fd)
c
c        write(*,*)  d3b, d3a, (fa*d3b+cfa*d3a)
c
        br=abde*ifpi*10.d0**((cfa*d3a)+(fa*d3b))
        hbeta=dmax1(0.d0,br)
c
c     write(*,*)' hydro.f line 259: HB: ', hbeta
c     Do H series
c
c      ratio and other data is only 10 deep atm
c
        ni=10
        l=1
        m=ni
c
   21   if (m-l.gt.1) then
        k=(m+l)/2
        if (th42(k).gt.tlo) then
          m=k
        else
          l=k
        endif
        goto 21
        endif
c
        logth0=th42(l)
        logth1=th42(m)
c
c        logth0=(tlo-logth0)/(logth1-logth0+epsilon)
c
        ft=(tlo-logth0)/(logth1-logth0)
        cft=1.d0-ft
c
        do series=1,nhseries
          do line=1,nhlines
C
C           hb1b=hylratsb(line,series,l,i)
C           hb2b=hylratsb(line,series,m,i)
C           hb3b=hylratsb(line,series,l,j)
C           hb4b=hylratsb(line,series,m,j)
C
C           d1b=(cft*hb1b)+(hb2b*ft)
C           d2b=(cft*hb3b)+(hb4b*ft)
C           d3b=(cfd*d1b)+(d2b*fd)
C
            do idx=1,ni
              y (idx)=hylratsb (line,series,idx,i)
              y2(idx)=hylratsb2(line,series,idx,i)
            enddo
            d1b=fsplint(th42,y,y2,ni,tlo)
            do idx=1,ni
              y (idx)=hylratsb (line,series,idx,j)
              y2(idx)=hylratsb2(line,series,idx,j)
            enddo
            d2b=fsplint(th42,y,y2,ni,tlo)
            d3b=(cfd*d1b)+(d2b*fd)
C
C           hb1a=hylratsa(line,series,l,i)
C           hb2a=hylratsa(line,series,m,i)
C           hb3a=hylratsa(line,series,l,j)
C           hb4a=hylratsa(line,series,m,j)
C
C           d1a=(cft*hb1a)+(hb2a*ft)
C           d2a=(cft*hb3a)+(hb4a*ft)
C           d3a=(cfd*d1a)+(d2a*fd)
C
            do idx=1,ni
              y (idx)=hylratsa (line,series,idx,i)
              y2(idx)=hylratsa2(line,series,idx,i)
            enddo
            d1a=fsplint(th42,y,y2,ni,tlo)
            do idx=1,ni
              y (idx)=hylratsa (line,series,idx,j)
              y2(idx)=hylratsa2(line,series,idx,j)
            enddo
            d2a=fsplint(th42,y,y2,ni,tlo)
            d3a=(cfd*d1a)+(d2a*fd)
C
            br=hbeta*10.d0**((cfa*d3a)+(fa*d3b))
            hydrobri(line,series)=dmax1(0.d0,br)
c
          enddo
        enddo
c
c     Interpolate 2S0 recombination rate
c
        d1b=(cft*r2s1b(l,i))+(r2s1b(m,i)*ft)
        d2b=(cft*r2s1b(l,j))+(r2s1b(m,j)*ft)
        d1a=(cft*r2s1a(l,i))+(r2s1a(m,i)*ft)
        d2a=(cft*r2s1a(l,j))+(r2s1a(m,j)*ft)
c
        d3b=(cfd*d1b)+(d2b*fd)
        d3a=(cfd*d1a)+(d2a*fd)
c
        br=10.d0**((cfa*d3a)+(d3b*fa))
        recrate2p(1)=dmax1(0.d0,br)
c
c        write(*,*) 'H  rec 2p',recrate2p(1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Look up He II 4686 temperature and do He II series.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        abde=((de*dh)*zion(2))*pop(3,2)
c
        fa=caseab(2)
        cfa=1.d0-fa
c
        ni=12
        l=1
        m=ni
c
   30   if (m-l.gt.1) then
        k=(m+l)/2
        if (the43(k).gt.tlo) then
          m=k
        else
          l=k
        endif
        goto 30
        endif
c
        logth0=the43(l)
        logth1=the43(m)
c
        ft=(tlo-logth0)/(logth1-logth0)
        cft=1.d0-ft
c c
c         d1b=(cft*he43b(l,i))+(he43b(m,i)*ft)
c         d2b=(cft*he43b(l,j))+(he43b(m,j)*ft)
c         d3b=(cfd*d1b)+(d2b*fd)
c c
c         d1a=(cft*he43a(l,i))+(he43a(m,i)*ft)
c         d2a=(cft*he43a(l,j))+(he43a(m,j)*ft)
c         d3a=(cfd*d1a)+(d2a*fd)
c c
        do idx=1,ni
          y (idx)=he43b (idx,i)
          y2(idx)=he43b2(idx,i)
        enddo
        d1b=fsplint(the43,y,y2,ni,tlo)
        do idx=1,ni
          y (idx)=he43b (idx,j)
          y2(idx)=he43b2(idx,j)
        enddo
        d2b=fsplint(the43,y,y2,ni,tlo)
c
        do idx=1,ni
          y (idx)=he43a (idx,i)
          y2(idx)=he43a2(idx,i)
        enddo
        d1a=fsplint(the43,y,y2,ni,tlo)
        do idx=1,ni
          y (idx)=he43a (idx,j)
          y2(idx)=he43a2(idx,j)
        enddo
        d2a=fsplint(the43,y,y2,ni,tlo)
c
c        write(*,*)  d1a, d2a, d1b, d2b
c
        d3b=(cfd*d1b)+(d2b*fd)
        d3a=(cfd*d1a)+(d2a*fd)
c
c        write(*,*)  d3b, d3a, (fa*d3b+cfa*d3a)
c
        br=abde*ifpi*10.d0**((cfa*d3a)+(fa*d3b))
c
        heiilr=dmax1(0.d0,br)
c
        do series=1,nheseries
          do line=1,nhelines
c c
c             hb1b=helratsb(line,series,l,i)
c             hb2b=helratsb(line,series,m,i)
c             hb3b=helratsb(line,series,l,j)
c             hb4b=helratsb(line,series,m,j)
c c
c             d1b=(cft*hb1b)+(hb2b*ft)
c             d2b=(cft*hb3b)+(hb4b*ft)
c             d3b=(cfd*d1b)+(d2b*fd)
c c
            do idx=1,ni
              y (idx)=helratsb (line,series,idx,i)
              y2(idx)=helratsb2(line,series,idx,i)
            enddo
            d1b=fsplint(the43,y,y2,ni,tlo)
            do idx=1,ni
              y (idx)=helratsb (line,series,idx,j)
              y2(idx)=helratsb2(line,series,idx,j)
            enddo
            d2b=fsplint(the43,y,y2,ni,tlo)
            d3b=(cfd*d1b)+(d2b*fd)
c c
c             hb1a=helratsa(line,series,l,i)
c             hb2a=helratsa(line,series,m,i)
c             hb3a=helratsa(line,series,l,j)
c             hb4a=helratsa(line,series,m,j)
c c
c             d1a=(cft*hb1a)+(hb2a*ft)
c             d2a=(cft*hb3a)+(hb4a*ft)
c             d3a=(cfd*d1a)+(d2a*fd)
c c
            do idx=1,ni
              y (idx)=helratsa (line,series,idx,i)
              y2(idx)=helratsa2(line,series,idx,i)
            enddo
            d1a=fsplint(the43,y,y2,ni,tlo)
            do idx=1,ni
              y (idx)=helratsa (line,series,idx,j)
              y2(idx)=helratsa2(line,series,idx,j)
            enddo
            d2a=fsplint(the43,y,y2,ni,tlo)
            d3a=(cfd*d1a)+(d2a*fd)
c
            br=heiilr*10.d0**((cfa*d3a)+(fa*d3b))
            helibri(line,series)=dmax1(0.d0,br)
c
          enddo
        enddo
c
c     Interpolate 2S0 recombination rate
c
        d1b=(cft*r2s2b(l,i))+(r2s2b(m,i)*ft)
        d2b=(cft*r2s2b(l,j))+(r2s2b(m,j)*ft)
        d1a=(cft*r2s2a(l,i))+(r2s2a(m,i)*ft)
        d2a=(cft*r2s2a(l,j))+(r2s2a(m,j)*ft)
c
        d3b=(cfd*d1b)+(d2b*fd)
        d3a=(cfd*d1a)+(d2a*fd)
c
        br=10.d0**((cfa*d3a)+(d3b*fa))
        recrate2p(2)=dmax1(0.d0,br)
c
      endif
c
      if (atypes.gt.2) then
c
c     Do two series for heavy atoms, based on extrapolation of case A
c      hydrogen
c
        do atom=3,atypes
          do series=1,nxhseries
            do line=1,nxhlines
              xhbeta(atom)=0.d0
              xhydrobri(line,series,atom)=0.d0
            enddo
          enddo
        enddo

        do atom=3,atypes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ***FINDS ABSOLUTE FLUX FOR Heavy Z H-BETA
c     Assume case A, this may not hold
c     in metal rich models
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          nz=mapz(atom)
          z2=dble(nz*nz)
          abde=zion(atom)*pop(nz+1,atom)
          if (abde.gt.pzlimit) then
            abde=abde*de*dh
            telec=trecom/z2
            z3=dble(nz*nz*nz)
            if ((telec.ge.tmin).and.(telec.le.tmax)) then
              tlo=dlog10(telec)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Look up equivalent H-Beta temperature
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
              ni=10
              l=1
              m=ni
c
   40         if (m-l.gt.1) then
              k=(m+l)/2
              if (th42(k).gt.tlo) then
                m=k
              else
                l=k
              endif
              goto 40
              endif
c
              logth0=th42(l)
              logth1=th42(m)
c
              ft=(tlo-logth0)/(logth1-logth0)
              cft=1.d0-ft
c
c     Interpolate high Z Hbeta
c
              d1a=(cft*h42a(l,i))+(h42a(m,i)*ft)
              d2a=(cft*h42a(l,j))+(h42a(m,j)*ft)
c
              d3a=10.d0**((cfd*d1a)+(d2a*fd))
              xhbeta(atom)=z3*d3a*abde*ifpi
c
c     Do H series by ratio
c
              do series=1,nxhseries
                do line=1,nxhlines
c
                  hb1a=hylratsa(line,series,l,i)
                  hb2a=hylratsa(line,series,m,i)
                  hb3a=hylratsa(line,series,l,j)
                  hb4a=hylratsa(line,series,m,j)

                  d1a=(cft*hb1a)+(hb2a*ft)
                  d2a=(cft*hb3a)+(hb4a*ft)

                  d3a=(cfd*d1a)+(d2a*fd)
                  br=xhbeta(atom)*10.d0**d3a
                  xhydrobri(line,series,atom)=dmax1(0.d0,br)
c
                enddo
              enddo
c  End of IF in temperature range.
            endif
c  End of abundance.
          endif
c  End of loop over atoms.
        enddo
c End of have heavies.
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Extrapolate high Z 2S rec rates by Z and T/Z^2
c     Based on H data.
c     Assume case A, this may not hold
c     in metal rich models
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (atypes.gt.2) then
        do atom=3,atypes
c
          nz=mapz(atom)
          abde=zion(atom)*pop(nz,atom)
          if (abde.gt.pzlimit) then

          z2=dble(nz*nz)
          telec=trecom+epsilon
          tz=telec/z2
          tlo=dlog10(tz)
c
          if ((tz.ge.tmin).and.(tz.le.tmax)) then
c
            ni=10
            l=1
            m=ni
c
   50       if (m-l.gt.1) then
            k=(m+l)/2
            if (th42(k).gt.tlo) then
              m=k
            else
              l=k
            endif
            goto 50
            endif
c
            logth0=th42(l)
            logth1=th42(m)
c
            ft=(tlo-logth0)/(logth1-logth0)
            cft=1.d0-ft
c
c     Interpolate 2S0 recombination rate
c
            d1a=(cft*r2s1a(l,i))+(r2s1a(m,i)*ft)
            d2a=(cft*r2s1a(l,j))+(r2s1a(m,j)*ft)
c
            d3a=10.d0**(d1a+(d2a-d1a)*fd)
            recrate2p(atom)=nz*d3a
c
          endif
          endif
c
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Add collisional excitation to recombination spectrum of H and He
c  only. Uses cascade branching rations to distribute collisions
c  to lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c H
c
      atom=1
      nz=mapz(1)
      do series=1,nhseries
        do line=1,nhlines
c
          k=series+line
          egj=iphe*((1.0/(series*series))-(1.0/(k*k)))*nz*nz
c
          do clevel=2,15
c
            if (clevel.ge.k) then
c
              cab=caseab(atom)*collhb(line,series,clevel)+(1.0-
     &         caseab(atom))*collha(line,series,clevel)
              hydrobri(line,series)=hydrobri(line,series)+cab*
     &         rateton(atom,clevel)*egj*ifpi
c
            endif
c
          enddo
c
        enddo
      enddo
c He
      atom=2
      nz=mapz(1)
      do series=1,nheseries
        do line=1,nhelines
c
          k=series+line
          egj=iphe*((1.0/(series*series))-(1.0/(k*k)))*nz*nz
c
          do clevel=2,15
c
            if (clevel.ge.k) then
c
              cab=caseab(atom)*collhb(line,series,clevel)+(1.0-
     &         caseab(atom))*collha(line,series,clevel)
              helibri(line,series)=helibri(line,series)+cab*
     &         rateton(atom,clevel)*egj*ifpi
c
            endif
c
          enddo
c
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Adds collisional excitation to recombination spectrum of heavies.
c     not used anymore as resonance cascade takes this into account,
c     the heavy hydrogenic recombination lines now have no collisions
c     components, so the lines in spec.csv can just be added.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      if (atypes.gt.3) then
c        do atom=3,atypes
cc
c          nz=mapz(atom)
cc
c          do series=1,nxhseries
c            do line=1,nxhlines
cc
c              if (xhydrobri(line,series,atom).gt.epsilon) then
cc
c                k=series+line
cc
c                egj=iphe*((1.0/(series*series))-(1.0/(k*k)))*nz*nz
cc
c                do clevel=2,15
cc
c                  if (clevel.ge.k) then
cc
c                    cab=caseab(atom)*collhb(line,series,clevel)+(1.0-
c     &               caseab(atom))*collha(line,series,clevel)
c                    xhydrobri(line,series,atom)=xhydrobri(line,series,
c     &               atom)+cab*rateton(atom,clevel)*egj*ifpi
cc
c                  endif
cc
c                enddo
cc
c              endif
c            enddo
c          enddo
cc
c        enddo
c      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set global H-beta and He II 4686 to include the collisional rates.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      halpha=hydrobri(1,2)
      lyalpha=hydrobri(1,1)
      hbeta=hydrobri(2,2)
      heiilr=helibri(1,3)
C      egj=iphe*((1.d0/(2*2))-(1.d0/(4*4)))
C      cab=caseab(1)*collhb(2,2,4)+(1.0-caseab(1))*collha(2,2,4)
C      write(*,*)' hydro.f : ', hbeta,pop(1,1),cab,rateton(1,4),egj
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Find two-photon emission, calling old hydro spec calculations in
c     the process, but not using them. They are needed for HeI though.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call hydro2p (t, de, dh)
c
      do atom=1,atypes
c
c     Test new 2 photon calcs
c
        ion=mapz(atom)
c
        z=dble(ion)
        z2=z*z
c
        ab1=zion(atom)*pop(nz+1,atom)
        ab0=zion(atom)*pop(nz,atom)

        if ((ab0.gt.pzlimit).or.(ab1.gt.pzlimit)) then
c
          ez=0.75d0*ipote(nz,atom)
c
c     scale temperature
c
          tz=1.d-4*t/z2
c
          r2q=de*dh*(ab0*fr(tz)*collrate2p(atom)+ab1*recrate2p(atom))
c
c     proton & electron deexcitation collision
c
          ab=dh*zion(1)*pop(2,1)
          dexcoll=(ab*qpr(tz)+de*qel(tz))/z
c
c     A scales ~ as nz^6
c
          a21=ahi(z,alphafs)
c          a21=8.2249d0*nz6
          r2q=r2q/(1.0d0+(dexcoll/a21))
c
          h2qbri(atom)=ez*r2q*ifpi
c
        endif
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Find intensity of important He I lines and the resultant
c     cooling rate (uses results from old hydro spectrum code).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call helioi (t, de, dh, helos)
c
      hloss=hloss+helos
      if (hloss.lt.epsilon) hloss=0.d0
c
      return
c
      end
