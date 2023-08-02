cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES INTERCOMBINATION LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
c     RETURNS DATA IN COMMON BLOCK /CLINE/
c     FSLOS (ERG.CM-3.S-1)
c
      subroutine inter (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
      real*8 rp,q21,q12,aa,pn2,fbr,omg
      real*8 f, pz,ba,t2,frkt
      real*8 ratekappa
      integer*4 atom,ion,m
c
      real*8 fkenhance
c
      fslos=0.0d0
      if (mlines.lt.1) return
c
c    ***COMPUTES RATES Q12,Q21 , HENCE POPULATION OF LEVEL 2
c
      f=dsqrt(1.0d0/(t+epsilon))
      t2=t*0.01d0
      frkt=1.d0/(rkb*t)
      do m=1,mlines
        atom=ielfs(m)
        ion=ionfs(m)
        pz=zion(atom)*pop(ion,atom)

        if (pz.ge.pzlimit) then
c
        fbr=0.0d0
        fsbri(m)=0.d0
c
          aa=e12fs(m)*frkt
          if (aa.lt.maxdekt) then
            ba=dexp(-aa)
            omg=omfs(m)
            ratekappa=1.d0
            if (usekappa) then
              ratekappa=fkenhance(kappa,aa)
            endif
            q12=(rka*f)*omg*ba/w1fs(m)*ratekappa
            q21=(rka*f)*omg/w2fs(m)*ratekappa
            rp=(de*q12)/(a21fs(m)+(de*q21))
            pn2=dh*pz*(rp/(1.0d0+rp))
            fbr=(pn2*a21fs(m))*e12fs(m)
            coolz(atom)=coolz(atom)+fbr
            coolzion(ion,atom)=coolzion(ion,atom)+fbr
            fslos=fslos+fbr
            fsbri(m)=fbr*ifpi
          endif
c
        endif
c
      enddo
c
      if (fslos.lt.epsilon) fslos=0.d0
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
c     Subroutine to calculate the transitions
c     from metatsable levels of highly ionised ions.
c     long lambda and low ionisations are handled in inter.
c
c     Uses the collision strength/statwei provided
c
c
c     refs: Landini & Monsignori Fosse 1990 A.A.Suppl.Ser. 82:229
c        Mewe 1985 A.A.Suppl.Ser. 45:11
c
c     RSS 9/90
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     subroutine inter2 (t, de, dh)
Cc
C     include 'cblocks.inc'
Cc
C     real*8 t, de, dh
C     real*8 egj,ejk,rr,frkt
C     real*8 cgjb,omegab
C     real*8 abde,pz,y,loss
Cc
C     real*8 ratekappa
Cc
C     integer*4 i,atom,ion
Cc
C     real*8 fkenhance
Cc
C     xiloss=0.0d0
C     if (xilines.lt.1) return
C
C     t=dmax1(0.1d0,t)
C     frkt=1.d0/(rkb*t)
Cc
C     do i=1,xilines
Cc
Cc     get data from arrays for line # i
Cc
Cc     Note: only line for possible species are read in
Cc     so no need to check maxion etc...
Cc
C       atom=xiat(i)
C       ion=xiion(i)
Cc
Cc     only calculate for abundant species
Cc
C       xrbri(i)=0.d0
Cc
C       pz=zion(atom)*pop(ion,atom)
Cc
C       if (pz.ge.pzlimit) then
Cc
C         abde=de*dh*pz
Cc
Cc     note that Egj does equal Ejk
Cc
C         egj=xiejk(i)*ev
C         y=egj*frkt
C         if (y.lt.maxdekt) then
Cc
C         ejk=egj
C         omegab=xomeg(i)
Cc
Cc     get scaled energy gap to level j from ground (not k)
Cc
C         ratekappa=1.d0
C         if (usekappa) then
C           ratekappa=fkenhance(kappa,y)
C         endif
Cc
Cc     transition power rate rr
Cc
C         rr=0.d0
C         cgjb=0.d0
C         cgjb=rka*(dexp(-y)/dsqrt(t))*omegab
Cc
C         rr=cgjb*ratekappa
Cc
Cc     number to transition abup
Cc
Cc
Cc     total power in line
Cc
C         loss=ejk*abde*rr
Cc
C         xiloss=xiloss+loss
C         coolz(atom)=coolz(atom)+loss
C         coolzion(ion,atom)=coolzion(ion,atom)+loss
Cc
C         xibri(i)=loss*ifpi
C         endif
Cc
Cc     end population limited loop
Cc
C       endif
Cc
Cc
C     enddo
Cc
C     if (xiloss.lt.epsilon) xiloss=0.d0
Cc
C     return
Cc
C     end
