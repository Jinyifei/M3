cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO DERIVE PREIONISATION CONDITIONS IN SHOCKS
c     OUTPUT IONISATION STATES IN POP(6,11) /ELABU/
c     FI : FILLING FACTOR , LTERM : OUTPUT DEVICE NUMBER
c     DH : PEAK H DENSITY , SHOCK VELOCITY : VS
c     TSMAX : MAXIMUM TIME STEP ALLOWED FOR IONISATION
c     USES PHOTON FIELD IN ARRAY : TPHOT
c     OUPUT TEMPERATURE : TEF ,  EL.DENS : DEF ,  AND
c     NUMBER OF PHOTOIONISING PHOTONS : QDIFT
c     CALL SUBR. INTVEC,TEEQUI,TAUDIST
c
c     NB.   SUBR. TOTPHOT NEEDS TO BE PREVIOUSLY CALLED
c
c
      subroutine preion (lterm, luop, tsmax, vs, fi, dh, def, tef, qtot,
     & drta)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 ct,dch,def,dh,dq,dq0,dq1,dqmi
      real*8 dqmm,dqp,dqpp,dqva,drta,dt,fhi,fi
      real*8 fma,frmm,frpp,frta,frti,frtma,frtp
      real*8 q1,q2,q3,qavail,qtot,qused,rad,t,tef,tmi
      real*8 ts0,ts1,tsm,tsmax,tstep,u,vs,wei
      real*8 xhyf,zrr,ztr
c
      integer*4 i,j,l,lterm,luop,m,nf
c
      character nmod*4, mmod*4
      character jsp*4, lpm*4
c
c           Functions
c
      real*8 ff,fr
c
      ff(u,ct)=dmin1(0.99999d0,(ct*u)+1.d-5)
      fr(u,ct)=(1.d0-((1.d0-ff(u,ct))**ct))/ff(u,ct)
c
      rad=0.0d0
      ct=0.6321d0
      frtma=1.0d1/ct
      jsp=jspot
      jspot='YES'
      tem=-1.0d0
      ts1=0.d0
      dq1=0.d0
      nf=16
      dqmi=1.0d-2
      tmi=0.075d0
      tsm=tsmax
      if (tsm.le.0.d0) tsm=1.d36
      if (tef.lt.100.d0) then
        tef=1.d4
        lpm='FIXT'
      else
        lpm='EQ'
      endif
c
c
c    ***SET POPI(6,11) TO PROTO-IONIC POULATIONS
c     AND FINDS NUMBER OF PHOTONS AVAILABLE : QTOT
c
      call copypop (propop, popi)
c
      call intvec (tphot, q1, q2, q3, qtot)
      frti=0.07d0
      mmod='DIS'
      nmod='TIM'
      l=0
      if (lterm.gt.0) write (lterm,10) dh,vs,qtot,fi
      if (luop.gt.0) write (luop,10) dh,vs,qtot,fi
   10 format(/' HDENS.,VSHOC, QTOT AND FIL.F. : ' ,4(1pg11.3))
      if (lterm.gt.0) write (lterm,20)
      if (luop.gt.0) write (luop,20)
   20 format(/' #   <TE>',t15,'QAVAIL',t25,'QUSED',t34
     &,'FHI',t43,'TSTEP',t52,'DRTA',t62,'TAU',t72,'DQ'/)
c
c    ***ITERATES ON OPTICAL DEPTH : FRTA   UNTIL THE NUMBER OF
c     ABSORBED PHOTONS : QUSED  IS EQUAL TO THE NUMBER OF
c     PHOTONS AVAILABLE : QAVAIL
c
      ztr=0.5d0
      dqpp=2.d0
      dqmm=-2.d0
      frta=frti
      frtp=frta
      frmm=1.d38
      frpp=0.d0
      dt=1.d0
      dqp=0.d0
      do 70 m=1,nf
        if (m.eq.1) then
          frta=frti
          frtp=frta
          frmm=1.d38
          frpp=0.d0
        else
          zrr=dmin1(ztr,dmax1(0.03d0,ztr*(dt/tmi)))
          if ((dq1.gt.zrr).and.(dq1.lt.dqpp)) then
            dqpp=dq1
            frpp=frta
          endif
          if ((dq1.lt.(-zrr)).and.(dq1.gt.dqmm)) then
            dqmm=dq1
            frmm=frta
          endif
c
          fma=5.d0
          if (((((dabs(dq1)-dabs(dqp))/((dabs(dq1)+dabs(dqp))+1.d-36))
     &     .lt.0.05d0).and.(dabs(dq1).gt.0.5d0)).and.(m.gt.2)) fma=
     &     12.0d0
          if (((dqp/(dq1+epsilon)).lt.0.d0).and.(m.gt.2)) then
            wei=0.5d0*(((dabs(dqp)/((dabs(dqp)+dabs(dq1))+1.d-36))/
     &       0.5d0)**0.75d0)
            dch=dexp((wei*dlog(frta))+((1.d0-wei)*dlog(frtp)))
            frtp=frta
            frta=dch
          else
            dch=(qavail/qused)*frta
            if (dch.gt.frmm) then
              wei=0.5d0*(((dabs(dqmm)/((dabs(dqmm)+dabs(dq1))+1.d-36))/
     &         0.5d0)**0.75d0)
              dch=dexp((wei*dlog(frta))+((1.0d0-wei)*dlog(frmm)))
            else if (dch.lt.frpp) then
              wei=0.5d0*(((dabs(dqpp)/((dabs(dqpp)+dabs(dq1))+1.d-36))/
     &         0.5d0)**0.75d0)
              dch=dexp((wei*dlog(frta))+((1.0d0-wei)*dlog(frpp)))
            endif
            frtp=frta
            frta=dmin1(fma*frta,dmax1(frta/fma,dch))
          endif
        endif
        dqp=dq1
c
        call copypop (popi, popf)
c
        do 60 l=1,10
          u=frta
          wei=fr(u,ct)
          call averinto (wei, popf, popi, pop)
          call taudist (dh, fi, frta, drta, rad, pop, mmod)
          call copypop (popi, pop)
          tstep=dmin1(drta/vs,tsm)
          if (l.gt.1) tstep=dmax1(ts1*0.5d0,dmin1(2.d0*ts1,tstep,tsm))
          ts0=ts1
          ts1=tstep
          dt=dabs(ts1-ts0)/(ts1+ts0)
c
          if (lpm.eq.'EQ') then
            call teequi (tef, t, def, dh, tstep, nmod)
          else
            t=tef
            call timion (t, def, dh, xhyf, tstep)
          endif
c
          call copypop (pop, popf)
c
          qused=0.d0
          do 40 i=1,atypes
            do 30 j=1,maxion(i)
              qused=qused+((zion(i)*0.5d0)*dabs(popf(j,i)-popi(j,i)))
   30       continue
   40     continue
c
          qused=(((dh*fi)*drta)*qused)*fr(u,ct)
          qavail=qtot*tstep
          dq=(qavail-qused)/((qavail+qused)+1.d-36)
          dq0=dq1
          dq1=dq
          dqva=dabs(dq0-dq1)
          fhi=popf(1,1)
          tef=t
c
          if (lterm.gt.0) write (lterm,50) m,tef,qavail,qused,fhi,tstep,
     &     drta,frta,dq
          if (luop.gt.0) write (luop,50) m,tef,qavail,qused,fhi,tstep,
     &     drta,frta,dq
c
   50 format(' ',i2,f9.0,2(1pg10.3),3(1pg9.2),1pg10.3,0pf9.4)
          if ((((dabs(dq).lt.dqmi).or.(frta.gt.frtma))
     &     .or.((tstep.ge.tsm).and.(dq.ge.dqmi))).and.(l.gt.1)) goto 90
          if (((l.gt.1).and.(dt.lt.tmi)).and.(dabs(dq).gt.(2.0*tmi)))
     &     goto 70
          if ((l.gt.1).and.(dqva.lt.(dqmi/2.0))) goto 70
          if ((l.gt.3).and.(dt.lt.(2.0*tmi))) goto 70
   60   continue
   70 continue
c
      if (dabs(dq).gt.(4.0*dqmi)) then
        if (luop.gt.0) write (luop,80)
        if (lterm.gt.0) write (lterm,80)
   80 format(/'$$$$$$$$$$$$$$$$ CONVERGENCE WAS NOT REACHED ' ,
     &'DURING DETERMINATION OF PREIONISATION' /)
      endif
c
   90 if (tstep.lt.tsm) goto 110
      if (luop.gt.0) write (luop,100) tsm
      if (lterm.gt.0) write (lterm,100) tsm
c
  100 format(/' [[[[[[[[[[[[[[[ TIME STEP LIMITED BY MAXIMUM ' ,
     &'VALUE :',1pg10.3)
  110 continue
      jspot=jsp
c
      tem=-1.d0
c
      return
      end
