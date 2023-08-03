      subroutine teequi2 (tei, tef, edens, hdens, tstep, nmod)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     more iteration on dlos
c
c     FINDS THE EQUILIBRIUM TEMPERATURE AND THE
c     CORRESPONDING IONISING STATE OF THE GAS
c     AT EQUILIBRIUM IONISATION OF AFTER TIME STEP :STEP
c     TEI : INITIAL GUESS FOR TEMPERATURE
c     TEF,EDENS: FINAL TEMPERATURE AND ELECTRONIC DENSITY
c
c     NMOD = 'EQUI'  :  EQUILIBRIUM IONISATION
c     NMOD = 'TIM'   :  INITIAL IONIC POP. EVOLVED BY TSTEP SEC.
c
c     CALL SUBR. COOL,EQUION
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 a,b,c1,c2,cc,csub,delexi,delmin,dl0
      real*8 dl1,dlodef,dlpr,dt12,dvv,fadl,ftu,ftu0
      real*8 t2,tdef,timin,tl1,tl2,tl2a,tm1,tmin,tsub
      real*8 tti,u,varf,xhyf
      real*8 tei, teinit,tef, edens, hdens, tstep
      real*8 popul(mxion, mxelem)
c
      integer*4 iconsis,isis,jok,luop,luty,n,nf
c
      character nmod*4
c
      real*8 arctanh
c
      arctanh(u)=0.5d0*dlog((1.d0+u)/(1.d0-u))
c
      luty=6
      luop=23
c
      if ((nmod.ne.'EQUI').and.(nmod.ne.'TIM')) then
        write (luty,10) nmod
   10 format(/' MODE IMPROPERLY SET FOR SUBR. TEEQUI  : ',a4)
        stop
      endif
c
c    Fudge which sets everything constant once T =30 K
c
c
      delmin=4.d-4
      nf=6
      tmin=10.d0
      timin=30.0d0
      ftu0=0.0d0
      isis=0

      teinit=dmax1(timin,tei)
      tef=teinit
      tm1=teinit
c
      tl1=dlog(tm1)
c
      call equion (tm1, edens, hdens)
c
      call cool (tm1, edens, hdens)
c
      dl0=dlos
c
      if (dabs(dabs(dlos)-1.0d0).lt.1.d-6) dlos=dsign(1.0d0-1.d-6,dlos)
      if (dabs(dlos).gt.0.06d0) nf=nf+2
      if (dabs(dlos).gt.0.20d0) nf=nf+3
      c1=arctanh(dlos)
      tl2=tl1-dlos
      t2=dmax1(dexp(tl2),4.0d0*tmin)
      tl2=dlog(t2)
      tti=t2
      dl1=dlos
      tdef=tm1
      dlodef=dl1
      varf=0.02d0+(dabs(dl1)**0.70d0)
      jok=0
c
c     ***USING 2 TEMP. (TM1,T2)  AND 2 FRACT. RESID. (DLOS1,DLOS2)
c     IT FITS THE FUNCTION : ARCTANH(DLOS)=B+A*LN(T)
c     NEXT VALUE OF T=DEXP(-B/A)
c
      do 20 n=1,nf
        dlpr=dlos
        call equion (t2, edens, hdens)
        call cool (t2, edens, hdens)
c         write(*, *) ' TEEQUI:', t2, dlos
        if (dabs(dabs(dlos)-1.0d0).lt.1.d-6) dlos=dsign(1.0d0-1.d-6,
     &   dlos)
        if (dabs(dlos).le.dabs(dlodef)) then
          dlodef=dlos
          tdef=t2
        endif
c
        if (dabs(dlos).lt.1.d-6) goto 30
c
        fadl=dmin1(100.0d0,dmax1(0.5d0,dabs(dlpr)/(dabs(dlos)+1.d-5)))
c
        delexi=5.0d-5+(1.0d-5*dmin1(5.d0,6.d-1*fadl))
        if (dabs(dlos).lt.5.d-5) jok=jok+1
        if (((dabs(dlos).lt.delexi).or.(t2.lt.tmin)).or.((jok.gt.4)
     &   .and.(dabs(dlos).lt.(dble(jok)*0.6d-5)))) goto 30
c
c didnt converge = retry...
c
        c2=arctanh(dlos)
        csub=c2-c1
        tsub=tl2-tl1
        if (dabs(tsub).lt.1.0d-16) tsub=dsign(1.0d-16,tsub)
        if (tsub.eq.0.0d0) tsub=dsign(1.0d-16,csub)
c
c     **DOES NOT ALLOW  A  TO BE NEGATIVE
c
        a=dabs(csub/tsub)
        b=c1-(a*tl1)
        cc=c1
        c1=c2
        tl1=tl2
        tm1=t2
        tl2a=-(b/(dmax1(dabs(a),dabs(b)/20.0d0)*dsign(1.d0,a)))
        dt12=tl2a-tl1
        tl2=tl1+dsign(dmin1(varf,dabs(dt12)),dt12)
        t2=dexp(tl2)
        ftu=(t2-tm1)/dmax1(t2,tm1,tmin/10.d0)
        iconsis=idint(((1.1d0*ftu0)*ftu)/(1.d-20+dabs(ftu0*ftu)))
        ftu0=ftu
c
        dvv=dabs(tl2-tl1)
        if (iconsis.lt.0) then
          isis=isis+1
          varf=varf/(2.0d0+(1.d0/(isis*isis)))
          tl2=tl1+dsign(dmin1(varf,dabs(dt12)),dt12)
          t2=dexp(tl2)
        else if (iconsis.gt.0) then
          varf=dmin1(1.1d0,dmax1(5.0d-4,(1.3d0*varf)*(dmax1(1.5d-1,
     &     dmin1(1.d0,dvv/varf))**0.2d0)))
        endif
   20 continue
c
c     ***IF CONVERGENCE POOR , USES TEMP. WITH MIN. ABS(DLOS)
c
   30 if (dabs(dlos).lt.0.01d0) goto 40
c      write(luty, 120) dl0, teinit, dlos, t2, dlodef, tdef
      t2=tdef
   40 continue
c
c     ***USES FINAL TEMPERATURE
c
      if (dabs(dlos).gt.(2.0d0*dabs(dlodef))) t2=tdef
      tef=dmax1(t2,10.d0)

      call equion (tef, edens, hdens)

      call cool (tef, edens, hdens)
c
c
      return
c
      end