cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Apply correction factor to convert maxwell sveraged upsilons
c to kappa averaged upsilons, for the current global kappa.
c if no correction is available, return 1.0, ie no change.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fkappacol(t,i,j,ionidx)
c
      include 'cblocks.inc'
c
      real*8 t,y
      integer*4 ionidx, idx, tridx
      integer*4 i,j,l,nspl
      real*8 beta
      real*8 btc
      real*8 btx(mxnspl)
      real*8 bty(mxnspl)
      real*8 bty2(mxnspl)
c
      real*8 fsplint
c
      y=1.d0
      idx=fmkappaidx(ionidx)
      if (idx.gt.0) then
        if (i.ne.j) then
          tridx=nkaptridx(i,j,idx)
          nspl=kupsbtn(tridx,idx)
          do l=1,nspl
            btx(l)=kupssplx(l,tridx,idx)
            bty(l)=kupssply(l,kappaidx,tridx,idx)
            bty2(l)=kupssply2(l,kappaidx,tridx,idx)
          enddo
          btc=kupsbtc(tridx,idx)
          beta=(t/(t+btc))
          y=fsplint(btx,bty,bty2,nspl,beta)
          if (usekappainterp) then
            do l=1,nspl
              bty(l)=kupssply(l,kappaidx+1,tridx,idx)
              bty2(l)=kupssply2(l,kappaidx+1,tridx,idx)
            enddo
            y=(kappaa*y)+(kappab*fsplint(btx,bty,bty2,nspl,beta))
          endif
        endif
      endif
c
      fkappacol=y
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Obtain Upsilon (Maxwell Eff Coll. Strength) for a multi-level
c ion transition at a given temperature, as transition idx or i->j levels
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fupsilonij(t,i,j,ionidx)
c
      include 'cblocks.inc'
c
      real*8 t, t4, td, upsilon
      integer*4 ionidx, nt
      integer*4 i,j,l,nspl
      integer*4 omtype
      real*8 beta, invy, localt
      real*8 btc
      real*8 btx(mxnspl)
      real*8 bty(mxnspl)
      real*8 bty2(mxnspl)
c
      real*8 fsplint,fkappacol
c
      fupsilonij=0.d0
      if (i.eq.j) return
c
      localt= t
c
      nt=nfmtridx(j,i,ionidx)
c
      omtype=fmomdatatype(nt,ionidx)
      upsilon=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 0: MAPPINGS III const or powerlaw at t4 data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (omtype.eq.0) then
        t4=localt*1.0d-4
        td=tdepm(j,i,ionidx)
        if (td.ne.0.d0) then
          upsilon=omim(j,i,ionidx)*(t4**td)
        else
          upsilon=omim(j,i,ionidx)
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 1-8 BT92 spline for omega,
c Types 1, 3-5 linear, Types 2, 6-8 log10
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((omtype.ge.1).and.(omtype.le.8)) then
        nspl=fmombtn(nt,ionidx)
        do l=1,nspl
          btx(l)=fmombsplx(l,nt,ionidx)
          bty(l)=fmombsply(l,nt,ionidx)
          bty2(l)=fmombsply2(l,nt,ionidx)
        enddo
c
c  All scaling factors are converted to characteristic temperature
c for CHIANTI and other splines, so can be treated uniformly
c
        btc=fmombtc(nt,ionidx)
        beta=(localt/(localt+btc))
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
        if ((omtype.eq.2).or.((omtype.ge.6).and.(omtype.le.8))) then
c
c spline is actually fit to log(Ups), raise to get Ups..
c
          upsilon=dabs(10.d0**upsilon)
        endif
c
      endif
c
      if ((omtype.eq.13)) then
        nspl=fmombtn(nt,ionidx)
        do l=1,nspl
          btx(l)=fmombsplx(l,nt,ionidx)
          bty(l)=fmombsply(l,nt,ionidx)
          bty2(l)=fmombsply2(l,nt,ionidx)
        enddo
c
c All scaling factors are converted to characteristic temperature
c for CHIANTI and other splines, so can be treated uniformly
c
        btc=fmombtc(nt,ionidx)
        beta=(localt/(localt+btc))
        invy=rkb*localt/fmeij(nt,ionidx)
c     beta = (dx/(dx + btc))   bt92 type 2 scaling, indirect TC
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
c
c spline is actually scaled
c
        upsilon=upsilon*dlog(invy+2.71828182845905d0)
c
      endif
c
c If kappa averaged Upsilons are needed, apply correction to
c maxwell avaeraged upsilons. fkappacol returns 1.0 if no correction
c is available.
c
      if (usekappa) then
        upsilon=upsilon*fkappacol(t,i,j,ionidx)
      endif
c
      fupsilonij=dabs(upsilon)
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c glue function for upsilon, if only have transition id:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fupsilontr(t,nt,ionidx)
c
      include 'cblocks.inc'
c
      real*8 t, upsilon
      integer*4 ionidx, nt
      integer*4 i,j
c
      real*8 fupsilonij
c
      i=nfmlower(nt,ionidx)
      j=nfmupper(nt,ionidx)
      upsilon=fupsilonij(t,i,j,ionidx)
      fupsilontr=upsilon
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
c*******COMPUTES  MULTI-LEVEL ATOM COOLING FOR up to 9 level atoms
c
      subroutine multilevel (t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables, mxnl is the global multi-level n
c
      real*8 t, de, dh
c
      real*8 w(mxnl), e(mxnl, mxnl)
      real*8 a(mxnl, mxnl),om(mxnl,mxnl),re(mxnl,mxnl)
      real*8 rd(mxnl, mxnl), alph(mxnl+2, mxnl+1), x(mxnl+1)
      real*8 f,po,zi
      real*8 aa,aion,bion
      real*8 tkappa(mxnl, mxnl), ratekappa(mxnl, mxnl)
c
      integer*4 ionindex, lineindex,ju,jl,jll,maxupper
      integer*4  j, nl, nt, ion, atom,is
c
      real*8 fkenhance,fupsilonij
c
      fmloss=0.0d0
      if (nfmions.lt.1) return
c
      f=dsqrt(1.0d0/t)
c
      do ionindex=1,nfmions
c
c     only calculate abundant ions
c
c
        atom=fmatom(ionindex)
        ion=fmion(ionindex)
        is=mapz(atom)-ion+1
        nl=fmnl(ionindex)
        nt=nfmtrans(ionindex)
        do j=1,nt
          fmbri(j,ionindex)=0.d0
        enddo
        do j=1,nl
          x(j)=0.d0
          fmx(j,ionindex)=0.d0
        enddo
c
        po=pop(ion,atom)
        zi=zion(atom)
        aion=zi*po*dh
c
        if (po*zi.ge.pzlimit) then
c
          do jl=1,nl
c
c   upper | lower level depending on symmetry
c
            w(jl)=wim(jl,ionindex)
c
c prepare transtion matrices, in particular fill in E, A, Omega(T)
c
            do ju=1,nl
c  upper | lower level depending on symmtry
              e(ju,jl)=eim(ju,jl,ionindex)
              a(ju,jl)=aim(ju,jl,ionindex)
c note reverse ij:
              om(ju,jl)=fupsilonij(t,jl,ju,ionindex)
              ratekappa(ju,jl)=1.d0
              if (usekappa) then
                if (ju.ne.jl) then
                  tkappa(ju,jl)=tkexm(ju,jl,ionindex)/t
                else
                  tkappa(ju,jl)=0.0d0
                endif
                ratekappa(ju,jl)=fkenhance(kappa,tkappa(ju,jl))
              endif
            enddo
          enddo
c
c    ***SET UP COLL. EXCIT. AND DEEXCIT. MATRIX
c
          maxupper=1
          do jl=1,nl
            do ju=1,nl
              rd(ju,jl)=0.0d0
              if (ju.gt.jl) then
c                 ju = upper, jl = lower
                rd(ju,jl)=((rka*f)*om(ju,jl))/w(ju)
                rd(ju,jl)=rd(ju,jl)*ratekappa(ju,jl)
              endif
              re(ju,jl)=0.0d0
              aa=e(ju,jl)/(rkb*t)
              if (ju.lt.jl) then
                if (aa.lt.maxdekt) then
c                 ju - lower, jl = upper
                  re(ju,jl)=(((rka*f)*om(ju,jl))*dexp(-aa))/w(ju)
                  re(ju,jl)=re(ju,jl)*ratekappa(ju,jl)
                  if (ju.eq.1) maxupper=max0(jl,maxupper)
                endif
              endif
            enddo
          enddo
c
c          write(*,*) maxupper, nl
c          nl=maxupper
c
c    ***SET UP MATRIX ELEMENTS
c
          do 50 jl=1,nl
            do 40 ju=1,nl
              if (jl.eq.1) goto 10
              if (jl.eq.ju) goto 20
              alph(ju,jl)=(de*(rd(ju,jl)+re(ju,jl)))+a(ju,jl)
              goto 40
   10         alph(ju,jl)=1.0d0
              goto 40
   20         alph(ju,jl)=0.0d0
              do 30 jll=1,nl
                if (ju.eq.jll) goto 30
                alph(ju,jl)=alph(ju,jl)-(de*(rd(ju,jll)+re(ju,jll))+
     &           a(ju,jll))
   30         continue
   40       continue
   50     continue
          do jl=1,nl
            alph(nl+1,jl)=0.0d0
          enddo
          alph(nl+1,1)=1.0d0
c
c
c    ***SOLVE FOR LEVEL POPULATIONS  X
c
c          call mdiagn (alph, x, nl)
          call matsolvemulti (alph, x, nl)
c
c    ***CHECK ON NORMALISATION
c
c         da=0.0d0
c         do ju=1,nl
c           fmx(ju,ionindex)=0.d0
c           if (x(ju).lt.0.0d0) x(ju)=0.0d0
c           if (x(ju).gt.0.0d0) da=da+x(ju)
c         enddo
c
c  Save pops for resonance exictation
c
          do ju=1,nl
c          x(ju)=x(ju)/da
            fmx(ju,ionindex)=x(ju)
          enddo
c
c    ***GET ION COOLING
c
          nl=fmnl(ionindex)
          lineindex=0
          do jl=1,nl-1
            do ju=jl+1,nl
              lineindex=lineindex+1
              if (x(ju).gt.0.d0) then
                bion=x(ju)*e(ju,jl)*aion*a(ju,jl)
                fmbri(lineindex,ionindex)=bion*ifpi
                coolz(atom)=coolz(atom)+bion
                coolzion(ion,atom)=coolzion(ion,atom)+bion
                fmloss=fmloss+bion
                if ((mapz(atom).eq.8).and.(ion.eq.3)) then
                  if (lineindex.eq.8) then
                    oiii5007loss=fmbri(lineindex,ionindex)*fpi
                  endif
                endif
              endif
            enddo
          enddo
c
c     end abundant ion
c
        endif
c
c     next ion.. ionindex
c
      enddo
c
      if (fmloss.lt.epsilon) fmloss=0.d0
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
c Obtain Upsilon (Maxwell Eff Coll. Strength) for a multi-level
c iron ion transition at a given temperature, as transition
c idx or i->j levels
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function ffeupsilonij(t,i,j,ionidx)
c
      include 'cblocks.inc'
c
      real*8 t, t4, td, upsilon
      integer*4 ionidx, nt
      integer*4 i,j,l,nspl
      integer*4 omtype
      real*8 beta
      real*8 btc,invy
      real*8 btx(mxfenspl)
      real*8 bty(mxfenspl)
      real*8 bty2(mxfenspl)
c
      real*8 fsplint
c
      ffeupsilonij=0.d0
      if (i.eq.j) return
c
      nt=nfetridx(j,i,ionidx)
      omtype=feomdatatype(nt,ionidx)
      upsilon=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 0: MAPPINGS III const or powerlaw at t4 data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (omtype.eq.0) then
        t4=t*1.0d-4
        td=tdepfe(j,i,ionidx)
        if (td.ne.0.d0) then
          upsilon=omife(j,i,ionidx)*(t4**td)
        else
          upsilon=omife(j,i,ionidx)
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Type 1-8 BT92 spline for omega,
c Types 1, 3-5 linear, Types 2, 6-8 log10
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((omtype.ge.1).and.(omtype.le.8)) then
        nspl=feombtn(nt,ionidx)
        do l=1,nspl
          btx(l)=feombsplx(l,nt,ionidx)
          bty(l)=feombsply(l,nt,ionidx)
          bty2(l)=feombsply2(l,nt,ionidx)
        enddo
c
c All scaling factors are converted to characteristic temperature
c for CHIANTI and other splines, so can be treated uniformly
c
        btc=feombtc(nt,ionidx)
        beta=(t/(t+btc))
c     beta = (dx/(dx + btc))   bt92 type 2 scaling, indirect TC
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
        if ((omtype.eq.2).or.((omtype.ge.6).and.(omtype.le.8))) then
c
c spline is actually fit to log(Ups), raise to get Ups..
c
          upsilon=dabs(10.d0**upsilon)
        endif
c
      endif
      if ((omtype.eq.13)) then
        nspl=feombtn(nt,ionidx)
        do l=1,nspl
          btx(l)=feombsplx(l,nt,ionidx)
          bty(l)=feombsply(l,nt,ionidx)
          bty2(l)=feombsply2(l,nt,ionidx)
        enddo
c
c All scaling factors are converted to characteristic temperature
c for CHIANTI and other splines, so can be treated uniformly
c
        btc=feombtc(nt,ionidx)
        beta=(t/(t+btc))
        invy=rkb*t/feeij(nt,ionidx)
c     beta = (dx/(dx + btc))   bt92 type 2 scaling, indirect TC
        upsilon=fsplint(btx,bty,bty2,nspl,beta)
c
c spline is actually scaled
c
        upsilon=upsilon*dlog(invy+2.71828182845905d0)
c
      endif
c
c If kappa averaged Upsilons are needed, apply correction to
c maxwell avaeraged upsilons. fkappacol returns 1.0 if no correction
c is available.
c
c      if (usekappa) then
c        upsilon=upsilon*fkappacol(t,i,j,ionidx)
c      endif
c
      ffeupsilonij=dabs(upsilon)
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c glue function for upsilon, if only have transition id:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function ffeupsilontr(t,nt,ionidx)
c
      include 'cblocks.inc'
c
      real*8 t, upsilon
      integer*4 ionidx, nt
      integer*4 i,j
c
      real*8 ffeupsilonij
c
      i=nfelower(nt,ionidx)
      j=nfeupper(nt,ionidx)
      upsilon=ffeupsilonij(t,i,j,ionidx)
      ffeupsilontr=upsilon
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES  MULTI-LEVEL ATOM COOLING FOR iron like atoms
c
      subroutine multiiron (t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables, mxnl is the global multi-level n
c
      real*8 t, de, dh
c
      real*8 w(mxfenl), e(mxfenl,mxfenl)
      real*8 a(mxfenl, mxfenl),om(mxfenl,mxfenl),re(mxfenl,mxfenl)
      real*8 rd(mxfenl, mxfenl), alph(mxfenl+2, mxfenl+1), x(mxfenl+1)
      real*8 f,po,zi
      real*8 aa,aion,bion
      real*8 tkappa(mxfenl, mxfenl), ratekappa(mxfenl, mxfenl)
c
      integer*4 ionindex,lineindex,ju,jl,jll
      integer*4 j,nl,nt,ion,atom,maxupper
c
      real*8 fkenhance,ffeupsilonij
c
      feloss=0.0d0
c
      if (nfeions.lt.1) return
c
      f=dsqrt(1.0d0/t)
c
      do ionindex=1,nfeions
c
c     only calculate abundant ions
c
        atom=featom(ionindex)
        ion=feion(ionindex)
        nl=fenl(ionindex)
        nt=nfetrans(ionindex)
        do j=1,nt
          febri(j,ionindex)=0.d0
        enddo
        do j=1,nl
          x(j)=0.d0
          fex(j,ionindex)=0.d0
        enddo
c
        po=pop(ion,atom)
        zi=zion(atom)
        aion=zi*po*dh
c
        if (po*zi.ge.pzlimit) then
c
c
          do jl=1,nl
c
c   upper | lower level depending on symmetry
c
            w(jl)=wife(jl,ionindex)
c
c prepare transition matrices, in particular fill in E, A, Omega(T)
c
            do ju=1,nl
c  upper | lower level depending on symmetry
              e(ju,jl)=eife(ju,jl,ionindex)
              a(ju,jl)=aife(ju,jl,ionindex)
c
c note reverse ij:
c
              om(ju,jl)=ffeupsilonij(t,jl,ju,ionindex)
c
              ratekappa(ju,jl)=1.d0
              if (usekappa) then
                if (ju.ne.jl) then
                  tkappa(ju,jl)=tkexfe(ju,jl,ionindex)/t
                else
                  tkappa(ju,jl)=0.0d0
                endif
                ratekappa(ju,jl)=fkenhance(kappa,tkappa(ju,jl))
              endif
            enddo
          enddo
c
c    ***SET UP COLL. EXCIT. AND DEEXCIT. MATRIX
c
          maxupper=1
          do jl=1,nl
            do ju=1,nl
              rd(ju,jl)=0.0d0
              if (ju.gt.jl) then
c                 ju = upper, jl = lower
                rd(ju,jl)=((rka*f)*om(ju,jl))/w(ju)
                rd(ju,jl)=rd(ju,jl)*ratekappa(ju,jl)
              endif
              re(ju,jl)=0.0d0
              aa=e(ju,jl)/(rkb*t)
              if (ju.lt.jl) then
                if (aa.lt.maxdekt) then
c                   ju - lower, jl = upper
                  re(ju,jl)=(((rka*f)*om(ju,jl))*dexp(-aa))/w(ju)
                  re(ju,jl)=re(ju,jl)*ratekappa(ju,jl)
                  if (ju.eq.1) maxupper=max0(jl,maxupper)
                endif
              endif
            enddo
          enddo
c
c only compute solution up to maxdekt above ground
c
          nl=maxupper
c
c    ***SET UP MATRIX ELEMENTS
c
          do 50 jl=1,nl
            do 40 ju=1,nl
              if (jl.eq.1) goto 10
              if (jl.eq.ju) goto 20
              alph(ju,jl)=(de*(rd(ju,jl)+re(ju,jl)))+a(ju,jl)
              goto 40
   10         alph(ju,jl)=1.0d0
              goto 40
   20         alph(ju,jl)=0.0d0
              do 30 jll=1,nl
                if (ju.eq.jll) goto 30
                alph(ju,jl)=alph(ju,jl)-(de*(rd(ju,jll)+re(ju,jll))+
     &           a(ju,jll))
   30         continue
   40       continue
   50     continue
          do jl=1,nl
            alph(nl+1,jl)=0.0d0
          enddo
          alph(nl+1,1)=1.0d0
c
c    ***SOLVE FOR LEVEL POPULATIONS  X
c
          call matsolvefe (alph, x, nl)
c
c          call mdiagfe (alph, x, nl)
c
c    ***CHECK ON NORMALISATION and save pops
c
c         da=0.0d0
c         do ju=1,nl
c           fex(ju,ionindex)=0.d0
c           if (x(ju).lt.0.0d0) x(ju)=0.0d0
c           if (x(ju).gt.0.0d0) da=da+x(ju)
c         enddo
c          write(*,*) 'fe norm:',da
          do ju=1,nl
c            x(ju)=x(ju)/da
            fex(ju,ionindex)=x(ju)
          enddo
c
          nl=fenl(ionindex)
          lineindex=0
          do jl=1,nl-1
            do ju=jl+1,nl
              lineindex=lineindex+1
              if ((x(ju).gt.0.d0).and.(a(ju,jl).gt.0.d0)) then
                bion=x(ju)*e(ju,jl)*aion*a(ju,jl)
                febri(lineindex,ionindex)=bion*ifpi
                coolz(atom)=coolz(atom)+bion
                coolzion(ion,atom)=coolzion(ion,atom)+bion
                feloss=feloss+bion
              endif
            enddo
          enddo
c
c     end abundant ion
c
        endif
c
c     next ion.. ionindex
c
      enddo
c
      if (feloss.lt.epsilon) feloss=0.d0
c
      return
      end
c
