cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES RESONANCE LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
c     RETURNS RLOSS = RESONANCE COOLING RATE(ERG CM-3 S-1)
c     RBRI(I,J) = BRIGTHNESS OF EACH LINE (ERG CM-3 S-1 STER-1)
c     AVAILABLE IN COMMOM BLOCK /RLINE/
c     CALL SUBR. RATEC, USES FUNCTION FGAUNT
c
C     subroutine reson (t, de, dh)
Cc
C     include 'cblocks.inc'
Cc
Cc           Variables
Cc
C     real*8 t, de, dh
CC     real*8 aa,abup,e,f, y
CC     real*8 po,r12,rr,u,zi
CCc
CC     real*8 ratekappa
CCc
CC     integer*4 jel,jjj,m
Cc
Cc           Functions
Cc
CC     real*8 fgaunt
CC     real*8 fkenhance
CCc
CCc     internal functions
CCc
CC     real*8 a , x1, fr, lg, x
Cc
CC      a(x)=(((5.232d0*x)-6.652d0)*x)+2.878d0
CC     x1(x)=0.7179d0*(x**(-0.0584d0))
CC     fr(x)=0.588d0*(x**(-0.234d0))
CC     lg(x)=0.07d0*(x**(-0.085d0))
Cc
C     rloss=0.0d0
C      return
C
C     f=dsqrt(1.0d0/(t+epsilon))
C
C    ***GET R12 , HENCE COOLING RATE
C
C     do m=1,nlines
C
C       rr=0.0d0
C       jel=ielr(m)
C       jjj=ionr(m)
C       rbri(m)=0.d0
C       zi=zion(jel)
C       po=pop(jjj,jel)
C       if (zi*po.ge.pzlimit) then
C         e=e12r(m)
C         y=e/(rkb*t)
C         if ((y.gt.0.d0).and.(y.lt.maxdekt)) then
C         abup=de*dh*zi*po
C         ratekappa=1.d0
C         if (usekappa) then
C           ratekappa=fkenhance(kappa,aa)
C         endif
C         r12=((rka*f)*omres(m))*dexp(-y)*ratekappa
C         rr=((r12*abup)*e)*fgaunt(jjj,0,aa)
C         rr=rr/(0.413497d0*a(1.0d0/dble(jjj)))
C         coolz(jel)=coolz(jel)+rr
C         coolzion(jjj,jel)=coolzion(jjj,jel)+rr
C         rloss=rloss+rr
C         rbri(m)=rr*ifpi
C         endif
C
C       endif
C
C     enddo
C
C     if (rloss.lt.epsilon) rloss=0.d0
Cc
C     return
Cc
C     end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New resonance line data in XLINDAT, requires new
c     code so that the effective osc strength f' can be used
c     effectively.
c
c
c     calls the new function fresga to get gaunt factors needed.
c
c
c     refs: Landini & Monsignori Fosse 1990 A.A.Suppl.Ser. 82:229
c        Mewe 1985 A.A.Suppl.Ser. 45:11
c
c     RSS 8/90
c
c     NOTE: transition energy does not give Egj, it gives Ejk
c     and in some cases this may not be a good approximation.....
c     ie transition H6
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES RESONANCE LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     subroutine reson2 (t, de, dh)
Cc
C     include 'cblocks.inc'
Cc
Cc
Cc           Variables
Cc
C     real*8 t, de, dh
C     real*8 abup,rr,y,f,loss
C     real*8 egj,ejk,tr,fef,gbar,fresga
C     real*8 cgjb,omegab,gi,pz
Cc      real*8 e1
Cc
C     real*8 ratekappa
Cc
C     integer*4 i, ionindex, civindex, nvindex, oviindex
C     integer*4 atom,ion,isos,code
Cc
Cc           Functions
Cc
Cc      real*8 farint
C     real*8 fkenhance
Cc
Cc
C     xrloss=0.0d0
C     if (xlines.lt.1) return
Cc
C     t=dmax1(0.1d0,t)
C     f=1.d0/dsqrt(t)
Cc
C     gi=2.d0
Cc
C     civindex=0
C     do ionindex=1,nfmions
C       if (fmatom(ionindex).eq.zmap(6)) then
C         if (fmion(ionindex).eq.4) then
C           civindex=ionindex
C         endif
C       endif
C     enddo
Cc
C     nvindex=0
C     do ionindex=1,nfmions
C       if (fmatom(ionindex).eq.zmap(7)) then
C         if (fmion(ionindex).eq.5) then
C           nvindex=ionindex
C         endif
C       endif
C     enddo
Cc
C     oviindex=0
C     do ionindex=1,nfmions
C       if (fmatom(ionindex).eq.zmap(8)) then
C         if (fmion(ionindex).eq.6) then
C           oviindex=ionindex
C         endif
C       endif
C     enddo
Cc
Cc
C     do i=1,xlines
Cc
Cc     get data from arrays for line # i
Cc
Cc     Note: only line for possible species are read in
Cc     so no need to check maxion etc...
Cc
C       atom=xrat(i)
C       ion=xion(i)
C       isos=xiso(i)
Cc
C       xrbri(i)=0.d0
Cc
C       pz=zion(atom)*pop(ion,atom)
Cc
Cc     only calculate abundant species
Cc
C       if (pz.ge.pzlimit) then
Cc
Cc     note that Egj does not necesarily equal Ejk
Cc
C         egj=xegj(i)
C         ejk=xejk(i)
C         tr=xtrans(i)
C         fef=xfef(i)
Cc
Cc     get scaled energy gap to level j from ground (not k)
Cc
C         y=egj*ev/(rkb*t)
Cc
C         if (y.lt.maxdekt) then
Cc
C           ratekappa=1.d0
C           if (usekappa) then
C             ratekappa=fkenhance(kappa,y)
C           endif
Cc
Cc     get mean gaunt factor
Cc
C           gbar=fresga(atom,isos,tr,egj,t,code)
Cc
Cc     transition power rate rr
Cc
C           rr=0.d0
C           cgjb=0.d0
C           if (code.ge.0) then
C             omegab=pi8rt3*fef*gi*gbar*r_inf/egj
C           else
C             omegab=fef
C           endif
Cc
C           cgjb=rka*f*dexp(-y)
C           rr=cgjb*omegab*ratekappa
Cc
Cc     number to transition is in abup
Cc
Cc     photons cm^-3 s^-1
Cc
C           rr=de*dh*pz*rr
Cc
Cc     ergs...
Cc
C           loss=ev*ejk*rr
Cc
C           coolz(atom)=coolz(atom)+loss
C           coolzion(ion,atom)=coolzion(ion,atom)+loss
C           xrloss=xrloss+loss
C           xrbri(i)=loss*ifpi
Cc
C         endif
Cc
Cc     end population limited loop
Cc
C       endif
Cc
C     enddo
Cc
C     if (xrloss.lt.epsilon) xrloss=0.d0
Cc
C     return
Cc
C     end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fxr3omgspl(t,y,icol,idx)
c
      include 'cblocks.inc'
c
      real*8 t, y
      integer*4 idx, icol
c
      real*8 upsilon
      integer*4 l,nspl
      integer*4 omtype
      real*8 beta
      real*8 btc
      real*8 btx(mxxr3nspl)
      real*8 bty(mxxr3nspl)
      real*8 bty2(mxxr3nspl)
c
      real*8 fsplint
c
      upsilon=0.d0
      omtype=xr3col_typespl(icol,idx)
      if ((omtype.eq.3).or.(omtype.eq.13)) then
c
c should be type 3 or 13, most differences are lost on init
c as x, y, and y2 are made for all types. only some types need
c exponentiation at the end
c
        btc=xr3col_tc(icol,idx)
        beta=(t/(t+btc))
        nspl=xr3col_nspl(icol,idx)
        do l=1,nspl
          btx(l)=xr3col_x(l) ! uniform splines for all hhe data
          bty(l)=xr3col_y(l,icol,idx)
          bty2(l)=xr3col_y2(l,icol,idx)
        enddo
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
        if (omtype.eq.13) then
           upsilon=upsilon*dlog((1.d0/y)+2.71828182845905d0)
        endif
      else
        write (*,*) 'ERROR, Invalid spline type in fxr3omgspl:',omtype
        write (*,*) t,icol,idx
        stop
      endif
      fxr3omgspl=dmax1(0.d0,upsilon)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     refs: CHIANTI database 8.01
c     Values derived by RSS2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine solvexr3ion (t, de, dh, idx, bri, ni)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh, abde
      integer*4 idx, ni
      integer*4 atom, ion, is
      integer*4 icol, nc, line, nl
      integer*4 i,j,jj,kk,l
      integer*4 mspecies
      integer*4 j2phlevel,ngrnd
c
      real*8 y,omegaij,rr,cr,cgjb,pz,loss,invrkt,sum
      real*8 ratekappa,f,br,ee,invgi,emiss
      real*8 bri(mxxr3lvls,mxxr3lvls)
      real*8 eji(mxxr3lvls,mxxr3lvls)
      real*8 fl(mxxr3lvls)
c
C     integer*4 col_jj(mxxr3casc,mxxr3cols)
C     integer*4 col_kk(mxxr3casc,mxxr3cols)
C     real*8    col_br(mxxr3casc,mxxr3cols)
c
      real*8 fkenhance
      real*8 fxr3omgspl

      atom=xr3at(idx)
      if (atom.eq.0) return
      ion=xr3ion(idx)
      is=(mapz(atom)-ion+1)
      pz=zion(atom)*pop(ion,atom)
      if (pz.le.pzlimit) return
      abde=de*dh*pz

      j2phlevel=xr32S12j(idx)
c
      do i=1,ni
        do j=1,ni
          bri(i,j)=0.d0
          eji(i,j)=plk*cls*dabs(xr3ei(j,idx)-xr3ei(i,idx))
        enddo
      enddo
c
      do l=1,mxxr3lvls
        fl(l)=0.d0
      enddo
      fl(1)=1.d0
c
      mspecies=fmspindex(ion,atom)
      if (mspecies.gt.0) then
        sum=0.d0
        ngrnd=fmnl(mspecies)
        do l=1,ngrnd
          fl(l)=0.d0
          if (fmx(l,mspecies).gt.1.0d-3) then
            fl(l)=fmx(l,mspecies)
            sum=sum+fmx(l,mspecies)
          endif
        enddo
        do l=1,ngrnd
          fl(l)=fl(l)/sum
        enddo
      endif
c
      nc=nxr3ioncol(idx)
      do icol=1,nc
        nl=xr3col_nl(icol,idx)
C       do line=1,nl
C         col_jj(line,icol)=xr3col_jj(line,icol,idx)
C         col_kk(line,icol)=xr3col_kk(line,icol,idx)
C         col_br(line,icol)=xr3col_br(line,icol,idx)
C       enddo
      enddo
c
      t=dmax1(0.1d0,t)
      f=1.d0/dsqrt(t)
      invrkt=1.d0/(rkb*t)
c
      do icol=1,nc
c i = lower, j = upper
        i=xr3col_i(icol,idx)
        xr3col_fr(icol,idx)=fl(i)
        if (fl(i).gt.0.d0) then
        j=xr3col_j(icol,idx)
        ee=eji(i,j)
        y=ee*invrkt
        if ((y.gt.0.d0).and.(y.lt.maxdekt)) then
          ratekappa=1.d0
          if (usekappa) then
            ratekappa=fkenhance(kappa,y)
          endif
          invgi=xr3invgi(i,idx)
          omegaij=fxr3omgspl(t,y,icol,idx)*invgi
          cgjb=rka*f*dexp(-y)
          cr=fl(i)*cgjb*omegaij*ratekappa
          rr=abde*cr
          if (rr.gt.epsilon) then
            loss=rr*ee
            if (linecoolmode.eq.1) then
              if ((is.eq.2).and.(i.eq.1).and.
     &            (j.eq.j2phlevel)) then
                 loss=0.d0 ! dont add the 2photon loss here, testing
              endif
              if ((is.eq.1).and.(i.eq.1).and.
     &            (j.eq.j2phlevel)) then
                 loss=0.d0 ! dont add the 2photon loss here, testing
              endif
            endif
            xr3loss=xr3loss+loss
            coolz(atom)=coolz(atom)+loss
            coolzion(ion,atom)=coolzion(ion,atom)+loss
            nl=xr3col_nl(icol,idx)
            do line=1,nl
              jj=xr3col_jj(line,icol,idx) ! lower
              kk=xr3col_kk(line,icol,idx) ! higher
              br=xr3col_br(line,icol,idx)
              emiss=(rr*br*eji(jj,kk)) ! note index swap above
              bri(kk,jj)=bri(kk,jj)+emiss
            enddo
          endif
        endif
        endif
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reson3 (t, de, dh)
c
      include 'cblocks.inc'
c
c     Local  params and variables
c
      real*8 t, de, dh
      real*8 pz
c
      real*8 ion_bri(mxxr3lvls,mxxr3lvls)
c
      integer*4 i, j, j2phlevel
      integer*4 idx, atom, ion, is
      integer*4 ni,line
c
      xr3loss=0.d0
      if (nxr3ions.lt.1) return
c
      do line=1,nxr3lines
        xr3lines_bri(line)=0.d0
      enddo
      do atom=1,atypes
        if (mapz(atom).gt.2) collrate2p(atom)=0.d0
        if (mapz(atom).gt.1) collrate2phe(atom)=0.d0
      enddo
c
      do idx=1,nxr3ions
        atom=xr3at(idx)
        if (atom.ne.0) then
          ion=xr3ion(idx)
          is=mapz(atom)-ion+1
          pz=zion(atom)*pop(ion,atom)
          if (pz.ge.pzlimit) then
            j2phlevel=xr32S12j(idx)
            ni=xr3ni(idx)
            call solvexr3ion (t, de, dh, idx, ion_bri, ni)
            do i=1,ni
              do j=1,ni
                if (ion_bri(i,j).gt.epsilon) then
                  line=xr3lines_map(i,j,idx)
c
c get effective 2 photon coll rates including cascade contributions to
c j2phlevel,  j2phlevel = 0 for non H- or He-like ions and is ignored.
c Also don't let the 2photon energy appear as a monochrome line,
c the emission is spread added elsewhere in twophoton.f
c
             if ((is.eq.1).and.(i.eq.1).and.(j.eq.j2phlevel)) then
                     collrate2p(atom)=ion_bri(i,j)/xr3lines_egij(line)
                     collrate2p(atom)=collrate2p(atom)/(de*dh*pz)
                     ion_bri(i,j)=0.d0
             endif
             if ((is.eq.2).and.(i.eq.1).and.(j.eq.j2phlevel)) then
                     collrate2phe(atom)=ion_bri(i,j)/xr3lines_egij(line)
                     collrate2phe(atom)=collrate2phe(atom)/(de*dh*pz)
                     ion_bri(i,j)=0.d0
             endif
                  xr3lines_bri(line)=ion_bri(i,j)*ifpi
c
                endif
              enddo
            enddo
          endif
        endif
      enddo
      if (xr3loss.lt.epsilon) xr3loss=0.d0
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
      real*8 function fxrlomgspl(t,y,icol,idx)
c
      include 'cblocks.inc'
c
      real*8 t, y
      integer*4 idx, icol
c
      real*8 upsilon
      integer*4 l,nspl
      integer*4 omtype
      real*8 beta
      real*8 btc
      real*8 btx(mxxrlnspl)
      real*8 bty(mxxrlnspl)
      real*8 bty2(mxxrlnspl)
c
      real*8 fsplint
c
      upsilon=0.d0
      omtype=xrlcol_typespl(icol,idx)
      if ((omtype.eq.3).or.(omtype.eq.13)) then
c
c should be type 3 or 13, most differences are lost on init
c as x, y, and y2 are made for all types. only type 13 needs
c extra scaling at the end
c
        btc=xrlcol_tc(icol,idx)
        beta=(t/(t+btc))
        nspl=xrlcol_nspl(icol,idx)
        do l=1,nspl
          btx(l)=xrlcol_x(l) ! uniform splines for all hhe data
          bty(l)=xrlcol_y(l,icol,idx)
          bty2(l)=xrlcol_y2(l,icol,idx)
        enddo
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
        if (omtype.eq.13) then
           upsilon=upsilon*dlog((1.d0/y)+2.71828182845905d0)
        endif
      else
        write (*,*) 'ERROR, Invalid spline type in fxrlomgspl:',omtype
        write (*,*) t,icol,idx
        stop
      endif
      fxrlomgspl=dmax1(0.d0,upsilon)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     refs: CHIANTI database 8.01
c     Values derived by RSS2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine solvexrlion (t, de, dh, idx, bri, ni)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh, abde, sum
      integer*4 idx, ni, ngnd
      integer*4 atom, ion
      integer*4 icol, nc, line, nl
      integer*4 i,j,jj,kk,l,ngrnd
      integer*4 fespecies,mspecies
c
      real*8 y,omegaij,rr,cr,cgjb,pz,loss,invrkt
      real*8 ratekappa,f,br,ee,invgi,emiss
      real*8 bri(mxxrllvls,mxxrllvls)
      real*8 eji(mxxrllvls,mxxrllvls)
      real*8 fl(mxxrllvls)
c
C     integer*4 col_jj(mxxrlcasc,mxxrlcols)
C     integer*4 col_kk(mxxrlcasc,mxxrlcols)
C     real*8    col_br(mxxrlcasc,mxxrlcols)
c
      real*8 fkenhance
      real*8 fxrlomgspl
c
      atom=xrlat(idx)
      if (atom.eq.0) return
      ion=xrlion(idx)
c      is=(mapz(atom)-ion+1)
      pz=zion(atom)*pop(ion,atom)
      if (pz.le.pzlimit) return
c
      abde=de*dh*pz
c
      do l=1,mxxrllvls
        fl(l)=0.d0
      enddo
c
      fl(1)=1.d0
c
      mspecies=fmspindex(ion,atom)
      if (mspecies.gt.0) then
        sum=0.d0
        ngrnd=fmnl(mspecies)
        do l=1,ngrnd
          fl(l)=0.d0
          if (fmx(l,mspecies).gt.1.0d-3) then
            fl(l)=fmx(l,mspecies)
            sum=sum+fmx(l,mspecies)
          endif
        enddo
        do l=1,ngrnd
          fl(l)=fl(l)/sum
        enddo
      endif
c
      fespecies=fespindex(ion,atom)
      if (fespecies.gt.0) then
        ngnd=fenl(fespecies)
        if ((mapz(atom).eq.26).and.(ion.eq.3)) then
          ngnd=34
        endif
        sum=0.d0
        do l=1,ngnd
          fl(l)=0.d0
          if (fex(l,fespecies).gt.1.0d-3) then
             fl(l)=fex(l,fespecies)
             sum=sum+fex(l,fespecies)
          endif
        enddo
        do l=1,ngnd
          fl(l)=fl(l)/sum
        enddo
      endif
c
      do i=1,ni
        do j=1,ni
          bri(i,j)=0.d0
          eji(i,j)=plk*cls*dabs(xrlei(j,idx)-xrlei(i,idx))
        enddo
      enddo
c
      nc=nxrlioncol(idx)
      do icol=1,nc
        nl=xrlcol_nl(icol,idx)
C       do line=1,nl
C         col_jj(line,icol)=xrlcol_jj(line,icol,idx)
C         col_kk(line,icol)=xrlcol_kk(line,icol,idx)
C         col_br(line,icol)=xrlcol_br(line,icol,idx)
C       enddo
      enddo
c
      t=dmax1(0.1d0,t)
      f=1.d0/dsqrt(t)
      invrkt=1.d0/(rkb*t)
c
      do icol=1,nc
c i = lower, j = upper
        i=xrlcol_i(icol,idx)
        xrlcol_fr(icol,idx)=fl(i)
        if (fl(i).gt.0.d0) then
        j=xrlcol_j(icol,idx)
        ee=eji(i,j)
        y=ee*invrkt
        if ((y.gt.0.d0).and.(y.lt.maxdekt)) then
          ratekappa=1.d0
          if (usekappa) then
            ratekappa=fkenhance(kappa,y)
          endif
          invgi=xrlinvgi(i,idx)
          omegaij=fxrlomgspl(t,y,icol,idx)*invgi
          cgjb=rka*f*dexp(-y)
          cr=fl(i)*cgjb*omegaij*ratekappa
          rr=abde*cr
          if (rr.gt.epsilon) then
            loss=rr*ee
            xrlloss=xrlloss+loss
            coolz(atom)=coolz(atom)+loss
            coolzion(ion,atom)=coolzion(ion,atom)+loss
            nl=xrlcol_nl(icol,idx)
            do line=1,nl
              jj=xrlcol_jj(line,icol,idx) ! lower
              kk=xrlcol_kk(line,icol,idx) ! higher
              br=xrlcol_br(line,icol,idx)
              emiss=(rr*br*eji(jj,kk)) ! note index swap above
              bri(kk,jj)=bri(kk,jj)+emiss
            enddo
          endif
        endif
        endif
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine resonl (t, de, dh)
c
      include 'cblocks.inc'
c
c     Local  params and variables
c
      real*8 t, de, dh
      real*8 pz
c
      real*8 ion_bri(mxxrllvls,mxxrllvls)
c
      integer*4 i, j
      integer*4 idx, atom, ion
      integer*4 ni,line
c
      xrlloss=0.d0
      if (nxrlions.lt.1) return
c
      do line=1,nxrllines
        xrllines_bri(line)=0.d0
      enddo
c
      do idx=1,nxrlions
        atom=xrlat(idx)
        if (atom.ne.0) then
          ion=xrlion(idx)
          ni=xrlni(idx)
          pz=zion(atom)*pop(ion,atom)
          if (pz.ge.pzlimit) then
            call solvexrlion (t, de, dh, idx, ion_bri, ni)
            do i=1,ni
              do j=1,ni
                if (ion_bri(i,j).gt.epsilon) then
                  line=xrllines_map(i,j,idx)
                  xrllines_bri(line)=ion_bri(i,j)*ifpi
                endif
              enddo
            enddo
          endif
        endif
      enddo
      if (xrlloss.lt.epsilon) xrlloss=0.d0
      return
      end
