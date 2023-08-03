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
c     SHOCK5: Shock model with steady Rankine-Hugoniot solution for
c     each step.  Time steps based on fastest atomic, cooling or
c     photon absorbion timescales.
c
c     Full continuum and diffuse field, auto preionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock5 ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iterations,iterationindex
   10 format(//,
     & ' *********************************************************',/,
     & '  SHOCK 5    Global Iteration: ',i2.2,' of ',i2.2)
   20 format(/,
     & '  SHOCK 5 Completed Iteration: ',i2.2,' of ',i2.2,/,
     & ' *********************************************************',/)
   30 format(/,
     & ' *********************************************************',/,
     & '  SHOCK 5 Model Completed',/,
     & ' *********************************************************',/)
c
      iterations=3
c
      call shock5setup (iterations)
      call shock5headers (iterations)
c
      iterationindex=0
      converged=0
      finalit=0
      if (iterations.le.1) finalit=1
      call compsh5 (0, iterations)
c
      if (iterations.gt.1) then
c
   40 continue
      iterationindex=iterationindex+1
c
      write (*,10) iterationindex,iterations
c
      call shock5precursor (iterationindex, iterations)
      call compsh5         (iterationindex, iterations)
      call shock5check     (iterationindex, iterations)
c
      write (*,20) iterationindex,iterations
c
      if ((converged.eq.0).and.(iterations.lt.12)) goto 40
c
c repeat a final model for outputs
c
      finalit=1
      iterationindex=iterationindex+1
      call shock5precursor (iterationindex, iterationindex)
      call compsh5         (iterationindex, iterationindex)
      call shock5check     (iterationindex, iterationindex)
c
      endif !have more than 1 iteration to try
c
      write (*,30)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock5setup (iterations)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iterations,i
      real*8 vapprox,tpr,tpo,va,eps
      real*8 feldens,fpressu,fpresse,frho,velshock2
c
   10 format(a)
c
      photonmode=1
c
      call zer ()!cleareverythingfirst
c
c     set up logical unit numbers
c
      lual=20
      luop=21
      lurt=22
      ludy=23
      lusp=24
      lupf=25
      lupb=26
      lucl=27
      lupc=28
      lupt=29
c
      ieln=4
      do i=1,ieln
        luions(i)=29+i
      enddo
c
c     zero arrays and define modes: all ions, no on-thespot,
c     full continuum
c     calculations..
c
      call zer
c
      jspot='NO'
      jcon='YES'
      jspec='NO'
      ispo='SII'
      allmod='NO'
      fclmod='NO'
      tsrmod='NO'
      ratmod='NO'
      dynmod='NO'
      bandsmod='NO'
      vmod='NONE'
      mypfx='v100sh'
      nprefix=6
c
      specmode='DW'
      jgeo='P'
      jtrans='LODW'
      jden='B'
c
      subname='Shock 5'
      jnorm=0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTRO - Define Input Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   20 format(///,
     & ' SHOCK 5 model selected:',/,
     & ' Shock code with auto-preionisation.',/,
     & ' *********************************************************')
      write (*,20)
c
c     get ionisation balance to work with...
c
      subname='proto-ionisation'
      call popcha (subname)
      call copypop (pop0, pop_neu)
      call copypop (pop0, pop_pre)
      call copypop (pop0, pop_pre0)
      call zeroemiss
      subname='Shock 5'
c
c     set up current field conditions
c
      call photsou (subname)
      do i=1,infph
          prefield(i)=soupho(i)
      enddo
c
c     Shock model preferences
c
   30 format(///,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Setting the shock conditions',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,30)
c
      stype='v'
   40 format(//,
     & ' Choose Shock Jump Parameter:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    T : Shock jump in terms of Temperature jump',/,
     & '    V : Shock jump in terms of flow Velocity.',//,
     & ' :: ',$)
c
   50 write (*,40)
      read (*,10) ilgg
      stype=ilgg(1:1)
      if (stype.eq.'T') stype='t'
      if (stype.eq.'V') stype='v'
c
      if ((stype.ne.'t').and.(stype.ne.'v')) goto 50
c
      magparam=1.0d0
      malpha=1.0d0
      mageta=1.0d0
      mmach=0.d0
c
      mtype='b'
   59 format(//,
     & ' Choose Shock Magnetic Parameter Type:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    B : Magnetic Field in microGauss at shockfront',/,
     & '    A : Magnetic Alpha_0 Pmag/Pgas in protostate',/,
     & '    M : Magnetic Alfven Mach Number, Ma',/,
     & '    C : Magnetic Alpha Pmag/Pgas at shockfront',/,
     & '    R : Magnetic Eta 2Pmag/Pram =1/Ma^2 at shockfront',/,
     & ' :: ',$)
   60 write (*,59)
      read (*,10) ilgg
      mtype=ilgg(1:1)
      if (mtype.eq.'A') mtype='a'
      if (mtype.eq.'B') mtype='b'
      if (mtype.eq.'C') mtype='c'
      if (mtype.eq.'R') mtype='r'
      if (mtype.eq.'M') mtype='m'
c
      if ((mtype.ne.'a').and.
     &    (mtype.ne.'b').and.
     &    (mtype.ne.'c').and.
     &    (mtype.ne.'m').and.
     &    (mtype.ne.'r')) goto 60
c
       if ( mtype.eq.'b') then
   61   format(//,' Specify the magnetic field at shockfront, B: ',/,
     &  ' (microgauss): ', $)
          write (*,61)
          read (*,*) bmag
          bmag=dabs(bmag)
          magparam=bmag
          bmag=bmag*1e-6
          pmag=(bmag*bmag)/epi
       endif
c
       if ( mtype.eq.'a') then
c     Preshock magnetic fields expressed as Alpha = Pmag/Pgas
c     Alpha = 1.0 is equipartition, Alpha = 10.0 is strong magnetic
c     Alpha = 0.1 is weak magnetic field. Alpha = 0.0 is no magnetic
c     Alpha = 1.0/Beta where Beta is the usual magnetic parameter, but
c     Alpha allows Alpha = 0.0 for the non-magnetic limit
   62   format(//,' Specify the magnetic Alpha_0 in proto-state',/,
     &  ' ( >= 0.0: Pmag/Pgas): ', $)
        write (*,62)
        read (*,*) magparam
        malpha=magparam
        if (magparam.lt.0.0d0) then
          bmag=-magparam
          bmag=bmag*1e-6
          pmag=(bmag*bmag)/epi
          mtype='b'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
         endif
       endif
c
       if ( mtype.eq.'c') then
c     Preshock magnetic fields expressed as Alpha = Pmag/Pgas
c     Alpha = 1.0 is equipartition, Alpha = 10.0 is strong magnetic
c     Alpha = 0.1 is weak magnetic field. Alpha = 0.0 is no magnetic
c     Alpha = 1.0/Beta where Beta is the usual magnetic parameter, but
c     Alpha allows Alpha = 0.0 for the non-magnetic limit
   63   format(//,' Specify the magnetic Alpha at shockfront',/,
     &  ' ( >= 0.0: Pmag/Pgas): ', $)
        write (*,63)
        read (*,*) magparam
        malpha=magparam
        if (magparam.lt.0.0d0) then
          bmag=-magparam
          bmag=bmag*1e-6
          pmag=(bmag*bmag)/epi
          mtype='b'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
       endif
c
       if ( mtype.eq.'m') then
c     Preshock magnetic fields expressed as Alpha = Pmag/Pgas
c     Alpha = 1.0 is equipartition, Alpha = 10.0 is strong magnetic
c     Alpha = 0.1 is weak magnetic field. Alpha = 0.0 is no magnetic
c     Alpha = 1.0/Beta where Beta is the usual magnetic parameter, but
c     Alpha allows Alpha = 0.0 for the non-magnetic limit
   55   format(//,' Specify the magnetic Alfven Mach Number, Ma',/,
     &  ' ( v/sqrt(Pmag/rho), Ma>= 0.0 ) ', $)
        write (*,55)
        read (*,*) magparam
        mmach=magparam
        if (magparam.lt.0.0d0) then
          bmag=-magparam
          bmag=bmag*1e-6
          pmag=(bmag*bmag)/epi
          mtype='b'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
       endif
c
       if ( mtype.eq.'r') then
c     Preshock magnetic fields expressed as Alpha_r = Pmag/Pram
c     Pram = rho*vs*vs, Pmag = (B^2)/8pi Alpha_r = B^2/8pi rho*vs*vs
c     B = sqrt(8pi rho*vs*vs)
c     Alpha allows Alpha = 0.0 for the non-magnetic limit
   64   format(//,' Specify the magnetic Eta at shockfront',/,
     &  ' ( >= 0.0: 2Pmag/Pram): ', $)
        write (*,64)
        read (*,*) magparam
        mageta=magparam
        if (magparam.lt.0.0d0) then
          bmag=dabs(magparam)*1e-6
          pmag=(bmag*bmag)/epi
          mtype='b'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
       endif
c
c     Shock in term of flow velocity
c
      wdil=0.5d0
      if (stype.eq.'v') then
   70   format(//,' Specify the proto-shock conditions:',/,
     &  '  T (K), dh (N), v (< 1e5 km/s, >= 1e5 cm/s)',/,
     &  ' (T<10 taken as a log, dh <= 0 taken as a log),',/,
     &  ' : ',$)
        write (*,70)
        read (*,*) t,dh,ve
c
        dr=1.d0
c
        if (ve.lt.1.d5) ve=ve*1.d5
        if (t.le.10.d0) t=10.d0**t
        if (dh.le.0.d0) dh=10.d0**dh
        if (wdil.gt.0.5d0) wdil=0.5d0
        wdilt0=0.d0
c
        pgas=fpressu(t,dh,pop_neu)
        de=feldens(dh,pop_neu)
        rh_neu=frho(de,dh)
        pram=rh_neu*ve*ve
c
        if (mtype.eq.'b') then
          bmag=magparam*1.0d-6
          pmag=(bmag*bmag)/epi
        endif
        if (mtype.eq.'a') then
          pmag=malpha*pgas
          bmag=dsqrt(epi*pmag)
        endif
        if (mtype.eq.'c') then
          pmag=malpha*pgas
          bmag=dsqrt(epi*pmag)
        endif
        if (mtype.eq.'m') then
          va=ve/mmach
          Pmag=(0.5*rh_neu*va*va)
          bmag=dsqrt(epi*pmag)
        endif
        if (mtype.eq.'r') then
          pmag=0.5d0*mageta*pram
          bmag=dsqrt(epi*pmag)
        endif
c
        vshoc=ve
c
        bm_neu=bmag
        te_neu=t
        de_neu=de
        dh_neu=dh
        vs_neu=ve
        pr_neu=fpresse(te_neu,de_neu,dh_neu)
        rh_neu=frho(de_neu,dh_neu)
c
c     fill in other default parameters
c
        dt=0.d0
        tloss=0.d0
        rmod='NEAR'
        call rankhug (t, de, dh, ve, bmag, tloss, dt)
        tm00=te0
      endif
c
c     Shock in terms of temp jump
c
      if (stype.eq.'t') then
c
   80   format(//,' Specify the proto-shock conditions:',/,
     &  ' (T<10 taken as a log, dh <= 0 taken as a log)',/,
     &  ' T (K), dh (N) : ',$)
        write (*,80)
        read (*,*) tpr,dh
        wdil=0.5d0
c
        if (tpr.le.10.d0) tpr=10.d0**tpr
        if (dh.le.0.0d0) dh=10.d0**dh
c
   90   format(//,' Specify the post-shock temperature:',/,
     &  ' (T<10 taken as a log) T (K) : ',$)
        write (*,90)
        read (*,*) tpo
c
        if (tpo.le.10.d0) tpo=10.d0**tpo
c
        te0=tpr
        te1=tpo
        te_neu=tpr
        dh_neu=dh
        pgas=fpressu(te_neu,dh_neu,pop_neu)
        de_neu=feldens(dh_neu,pop_neu)
        rh_neu=frho(de_neu,dh_neu)
c
        if (mtype.eq.'b') then
          bmag=magparam*1.0d-6
          pmag=(bmag*bmag)/epi
        endif
        if (mtype.eq.'a') then
          pmag=malpha*pgas
          bmag=dsqrt(epi*pmag)
        endif
        if (mtype.eq.'c') then
          pmag=malpha*pgas
          bmag=dsqrt(epi*pmag)
        endif
c
        vapprox= velshock2 (dh_neu, tpr, 0.d0, tpo)
        rh_neu=frho(de_neu,dh_neu)
        pram=rh_neu*vapprox*vapprox
        if (mtype.eq.'r') then
          pmag=0.5d0*mageta*pram
          bmag=dsqrt(epi*pmag)
        endif
c
   91    continue
         bm0=bmag
c
          vshoc=velshock2 (dh_neu, tpr, bmag, tpo)
          pram=rh_neu*vshoc*vshoc
c
           if (mtype.eq.'m') then
            va=vshoc/mmach
            Pmag=(0.5*rh_neu*va*va)
            bmag=dsqrt(epi*pmag)
            eps=dabs(2.0*(bmag-bm0)/(bmag+bm0))
            if (eps.gt.1.0d-6 ) goto 91
          endif
c
          if (mtype.eq.'r') then
c
c iterate for rampressure/alpha_r
c
            pmag=0.5d0*mageta*pram
            bmag=dsqrt(epi*pmag)
C           write(*,*) vshoc,pram,pmag,bm0,bmag,te1
            eps=dabs(2.0*(bmag-bm0)/(bmag+bm0))
            write(*,*) eps,bmag,bm0
            if (eps.gt.1.0d-6 ) goto 91
          endif
c
        dt=0.d0
        tloss=0.d0
        rmod='NEAR'
c        call isochorflow (t, de, dh, ve, bmag, rmod, tloss, dt)
c        call isobarflow (t, de, dh, ve, bmag, rmod, tloss, dt)
        call rankhug (te_neu, de_neu, dh_neu, vshoc,
     &         bmag, tloss, dt)
c
        tm00=te0
        bm_neu=bmag
        vs_neu=vshoc
        pr_neu=fpresse(te_neu,de_neu,dh_neu)
        rh_neu=frho(de_neu,dh_neu)
c
      endif
c
c Shock intial state kept for iterations, above numbers are
c reused from step to step
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0!shockvelocity
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
c
      pmag=bp0
      pgas=pr0
      pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaeos*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.0d0*pmag/rho0)
      malpha=pmag/pgas
      gaseta=gammaeos*pgas/pram
      mageta=2.d0*pmag/pram
c
      call shocksummary (6)
c
c set intial convergence vars
c
      cmpf0=cmpf
      psi0=0.0d0
      psi=0.0d0
      te_pre0=te_pre
      te_pst0=te_pst
      de_pre0=de_pre
      de_pst0=de_pst
      dh_pre0=dh_pre
      dh_pst0=dh_pst
      vs_pre0=vs_pre!shockvelocity
      vs_pst0=vs_pst
c
      rh_pre0=rh_pre
      rh_pst0=rh_pst
      pr_pre0=pr_pre
      pr_pst0=pr_pst
      bm_pre0=bm_pre
      bm_pst0=bm_pst
c
      if (alphacoolmode.eq.1) then
  100  format(//,
     &  ' Powerlaw Cooling enabled  Lambda T6 ^ alpha :',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Give index alpha, Lambda at 1e6K',/,
     &  '    (Lambda<0 as log)',/,
     &  ' :: ',$)
        write (*,100)
        read (*,*) alphaclaw,alphac0
        if (alphac0.lt.0.d0) alphac0=10.d0**alphac0
      endif
c
c     Diffuse field interaction
c
c
  110 format(//,
     & ' Choose Diffuse Field Option :',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    F  : Full diffuse field interaction (Default).',/,
     & '    Z  : Zero diffuse field interaction.',/,
     & ' :: ',$)
      write (*,110)
  120 read (*,10) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'z') ilgg='Z'
      if (ilgg.eq.'f') ilgg='F'
c
      if ((ilgg.ne.'Z').and.(ilgg.ne.'F')) goto 120
c
      if (ilgg.eq.'Z') photonmode=0
      if (ilgg.eq.'F') photonmode=1
c
  130 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Calculation Settings ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,130)
c
c     Set up time step limits
c
      tmod='A'
c
c no longer used directly, my return in the future
c
      atimefrac=0.050d0
      cltimefrac=0.050d0
      abtimefrac=1.0d0
      utime=0.d0
c
  140 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Boundary Conditions ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,140)
c
c     Choose ending
c
  150 format(/,
     & '  Choose Ending :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard ending, 1% weighted ionisation',/,
     & '    B  : Species ionisation limit',/,
     & '    C  : Temperature limit',/,
     & '    S  : Temperature limit, and >95% neutral',/,
     & '    D  : Distance limit',/,
     & '    E  : Time limit',/,
     & '    F  : Thermal balance limit',/,
     & '    G  : Heating limit',/,
     & ' :: ',$)
c
  160 write (*,150)
      read (*,10) jend
      jend=jend(1:1)
c
      if (jend.eq.'a') jend='A'
      if (jend.eq.'b') jend='B'
      if (jend.eq.'c') jend='C'
      if (jend.eq.'s') jend='S'
      if (jend.eq.'d') jend='D'
      if (jend.eq.'e') jend='E'
      if (jend.eq.'f') jend='F'
      if (jend.eq.'g') jend='G'
c
      if ((jend.ne.'A').and.(jend.ne.'B').and.(jend.ne.'C')
     &.and.(jend.ne.'D').and.(jend.ne.'E').and.(jend.ne.'F')
     &.and.(jend.ne.'G').and.(jend.ne.'S')) goto 160
c
c     Secondary info:
c
      if (jend.eq.'B') then
  170   format(//,' Give atom, ion and limit fraction: ',/,
     &   ' (eg 1 2 0.05 = stop when HII < 0.05): ',$)
        write (*,170)
        read (*,*) ielen,jpoen,fren
      endif
      if ((jend.eq.'C').or.(jend.eq.'S')) then
  180   format(//,' Give final temperature (K > 10, log <= 10): ',$)
        write (*,180)
        read (*,*) tend
        if (tend.le.10.d0) tend=10.d0**tend
      endif
      if (jend.eq.'D') then
  190   format(//,' Give final distance (cm > 100, log<=100): ',$)
        write (*,190)
        read (*,*) diend
        if (diend.le.100.d0) diend=10.d0**diend
      endif
      if (jend.eq.'E') then
  200   format(//,' Give time limit (s > 100, log<=100): ',$)
        write (*,200)
        read (*,*) timend
        if (timend.le.100.d0) timend=10.d0**timend
      endif
c
  210   format(//,
     & '  Choose the minimum number of shock iterations:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  ' (<= 1, no iterations) : ',$)
      write (*,210)
      read (*,*) iterations
      if (iterations.le.0) iterations=1
c
  220 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Output Requirements ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,220)
c
c     get final field file prefix
c
      mypfx='.               '
  230 format (a16)
  240 format(//,
     & ' Give a prefix for output files (max 8 chars): ')
      write (*,240)
      read (*,230) mypfx
      nprefix=1
  250 if ((mypfx(nprefix+1:nprefix+1).eq.' ').or.(nprefix.ge.8)) goto
     &260
      nprefix=nprefix+1
      goto 250
  260 continue
      mypfx=mypfx(1:nprefix)
c
c     default monitor elements
c
      iel1=zmap(1)
      iel2=zmap(2)
      iel3=zmap(6)
      iel4=zmap(8)
c
  270 format(//,' Choose output settings : ',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard output. ',/,
     & '    B  : Standard plus ion balance files.',/,
     & '    C  : Standard plus dynamics file.',/,
     & '    D  : Standard plus rates file.',/,
     & '    E  : Standard plus final down stream field.',/,
     & '    F  : Standard plus upstream field at each step.',/,
     & '    J  : Cooling Components by Elements file.',/,
     & '       :',/,
     & '    G  : B + E',/,
     & '    H  : B + F',/,
     & '       :',/,
     & '    I  : Everything',//,
     & ' :: ',$)
  280 write (*,270)
      read (*,10) ilgg
      ilgg=ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
      if (ilgg.eq.'d') ilgg='D'
      if (ilgg.eq.'e') ilgg='E'
      if (ilgg.eq.'f') ilgg='F'
      if (ilgg.eq.'g') ilgg='G'
      if (ilgg.eq.'h') ilgg='H'
      if (ilgg.eq.'i') ilgg='I'
      if (ilgg.eq.'j') ilgg='J'
c
      if ((ilgg.ne.'A').and.(ilgg.ne.'B').and.(ilgg.ne.'C')
     &.and.(ilgg.ne.'D').and.(ilgg.ne.'E').and.(ilgg.ne.'F')
     &.and.(ilgg.ne.'G').and.(ilgg.ne.'H').and.(ilgg.ne.'I')
     &.and.(ilgg.ne.'J')) goto 280
c
c     Upstream fields
c
      lmod='N'
      if ((ilgg.eq.'F').or.(ilgg.eq.'H').or.(ilgg.eq.'I')) then
        lmod='Y'
      endif
c
c     F-Lambda plots, ionising and optical
c
      jspec='NO'
      if ((ilgg.eq.'E').or.(ilgg.eq.'G').or.(ilgg.eq.'I')) then
        jspec='YES'
      endif
c
c     timescales and rates file
c
      ratmod='N'
      if ((ilgg.eq.'D').or.(ilgg.eq.'I')) then
        ratmod='Y'
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+5)='rates'
        np=nprefix+5
        sfx='sh5'
        call newfile (pfx, np, sfx, 3, fn)
        fr=fn(1:np+8)
        open (lurt,file=fr,status='NEW')
      endif
c
c     Cooling Components File
c
      fclmod='N'
      if ((ilgg.eq.'J').or.(ilgg.eq.'I')) then
c      lmod='Y'
c      jspec='YES'
        jnorm=0
  290  format (' Cooling File Normalisation,',/
     &         ' (0=ne.nH, 1=nH^2, 2=ne.ni, 3=n^2, 4=ne^2): ',$)
        write (*,290)
        read (*,*) jnorm
        if (jnorm.lt.0) jnorm=0
        if (jnorm.gt.4) jnorm=0
        fclmod='Y'
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+4)='cool'
        np=nprefix+4
        sfx='csv'
        fcl=' '
        call newfile (pfx, np, sfx, 3, fn)
        fcl=fn(1:np+8)
        open (lucl,file=fcl,status='NEW')
      endif
c
c     flow dynamics file
c
      dynmod='N'
      if ((ilgg.eq.'C').or.(ilgg.eq.'I')) then
        dynmod='Y'
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+3)='dyn'
        np=nprefix+3
        sfx='sh5'
        fd=' '
        call newfile (pfx, np, sfx, 3, fn)
        fd=fn(1:np+8)
        open (ludy,file=fd,status='NEW')
      endif
c
      pfx=mypfx(1:nprefix)
      pfx(nprefix+1:nprefix+4)='spec'
      np=nprefix+4
      sfx='csv'
      fsp=' '
      call newfile (pfx, np, sfx, 3, fn)
      fsp=fn(1:np+8)
      open (lusp,file=fsp,status='NEW')
c
      pfx=mypfx(1:nprefix)
      pfx(nprefix+1:nprefix+6)='pcspec'
      np=nprefix+6
      sfx='csv'
      fpc=' '
      call newfile (pfx, np, sfx, 3, fn)
      fpc=fn(1:np+8)
      open (lupc,file=fpc,status='NEW')
c
c     ionbalance files
c
      allmod='n'
      tsrmod='n'
c
      if (  (ilgg.eq.'B').or.(ilgg.eq.'G')
     &  .or.(ilgg.eq.'H').or.(ilgg.eq.'I') ) then
c
        if ((ilgg.eq.'B').or.(ilgg.eq.'G').or.(ilgg.eq.'H')) then
c
  300    format(//,' Record all ions file (Y/N)?: ',$)
          write (*,300)
          read (*,10) allmod
          allmod=allmod(1:1)
          if (allmod.eq.'y') allmod='Y'
        endif
c
        if (ilgg.eq.'I') allmod='Y'
c
        if (allmod.eq.'Y') then
          pfx=mypfx(1:nprefix)
          pfx(nprefix+1:nprefix+7)='_allion'
          np=nprefix+7
          sfx='sh5'
          call newfile (pfx, np, sfx, 3, fn)
          fa=fn(1:np+8)
          open (lual,file=fa,status='NEW')
        endif
c
        if (ilgg.eq.'B') then
  310    format(//,' Record particular ion balances (Y/N)?: ',$)
          write (*,310)
          read (*,10) tsrmod
          tsrmod=tsrmod(1:1)
          if (tsrmod.eq.'y') tsrmod='Y'
        endif
c
        if (ilgg.eq.'I') tsrmod='Y'
c
        if (tsrmod.eq.'Y') then
  320       format(//,' Give 4 elements by atomic number: ',$)
          write (*,320)
          read (*,*) iel1,iel2,iel3,iel4
          iel1=zmap(iel1)
          iel2=zmap(iel2)
          iel3=zmap(iel3)
          iel4=zmap(iel4)
        endif
      endif
c
c     get screen display mode
c
  330 format(//,' Choose runtime screen display:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard display ',/,
     & '    B  : Detailed Slab display.',/,
     & '    C  : Full Display',/,
     & '       : (Full slab display + timescales).',/,
     & ' :: ',$)
      write (*,330)
      read (*,10) ilgg
      ilgg=ilgg(1:1)
      vmod='MINI'
c
      if (ilgg.eq.'a') ilgg='A'
      if (ilgg.eq.'b') ilgg='B'
      if (ilgg.eq.'c') ilgg='C'
c
      if (ilgg.eq.'A') vmod='MINI'
      if (ilgg.eq.'B') vmod='SLAB'
      if (ilgg.eq.'C') vmod='FULL'
c
c     get runname
c
  340 format (a80)
  350 format(//,' Give a name/code for this run: ',$)
      write (*,350)
      read (*,340) runname
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shock5headers (iterations)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i,j,iterations
      character tab*1
      real*8 feldens
c
      tab=','
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Headers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion balance files if requested
c
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (i.eq.1) ie=iel1
          if (i.eq.2) ie=iel2
          if (i.eq.3) ie=iel3
          if (i.eq.4) ie=iel4
          fn=' '
          pfx=elem(ie)
          pfx=mypfx(1:nprefix)
          pfx(nprefix+1:nprefix+1)='_'
          pfx(nprefix+2:nprefix+2+elem_len(ie))=elem(ie)
          np=nprefix+1+elem_len(ie)
          sfx='csv'
          call newfile (pfx, np, sfx, 3, fn)
          filn(i)=fn(1:np+8)
        enddo
c
        do i=1,ieln
          open (luions(i),file=filn(i),status='NEW')
        enddo
c
        do i=1,ieln
          write (luions(i),*) fl,', ',runname
c
   10 format('Step [1], <X> [2], dX [3], <T> [4]',31(',',a6,'[',i2,']'))
          if (i.eq.1) then
            write (luions(i),10) (rom(j),j+4,j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),10) (rom(j),j+4,j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),10) (rom(j),j+4,j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),10) (rom(j),j+4,j=1,maxion(iel4))
          endif
c
        enddo
c
        do i=1,ieln
          close (luions(i))
        enddo
c
c     end trsmod = 'Y'
      endif
c
      de=feldens(dh,pop)
c
      fpb=' '
      if (bandsmod.eq.'Y') then
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+6)='_bands'
        np=nprefix+6
        sfx='sh5'
        call newfile (pfx, np, sfx, 3, fn)
        fpb=fn(1:np+8)
        open (lupb,file=fpb,status='NEW')
      endif
c
      fl=' '
      pfx='shck_'
      pfx(6:5+nprefix)=mypfx(1:nprefix)
      np=5+nprefix
      sfx='sh5'
      call newfile (pfx, np, sfx, 3, fn)
      fl=fn(1:np+8)
c
      open (luop,file=fl,status='NEW')
c
      fpt=' '
      pfx='prec_'
      pfx(6:5+nprefix)=mypfx(1:nprefix)
      np=5+nprefix
      sfx='sh5'
      call newfile (pfx, np, sfx, 3, fn)
      fpt=fn(1:np+8)
c
      open (lupt,file=fpt,status='NEW')
c
c     Title
c
   20 format(//,
     &       ' SHOCK 5: Steady Rankine-Hugoniot Shock Code: ',/,
     &       ' ============================================',/,
     &       ' Diffuse Field, Full Continuum Calculations.',/,
     &       ' Global Shock-Precursor Iterations: ',i2,/,
     &       ' Calculated by MAPPINGS V ',a8)
      write (*,20) iterations,theversion
      write (luop,20) iterations,theversion
      write (lupt,20) iterations,theversion
      write (lusp,20) iterations,theversion
      write (lupc,20) iterations,theversion
      if (bandsmod.eq.'Y') write (lupb,20) iterations,theversion
      if (ratmod.eq.'Y') write (lurt,20) iterations,theversion
      if (dynmod.eq.'Y') write (ludy,20) iterations,theversion
      if (allmod.eq.'Y') write (lual,20) iterations,theversion
      if (fclmod.eq.'Y') write (lucl,20) iterations,theversion
   30 format(//' Run  :',a,/
     &         ' File :',a)
      write (luop,30) runname,fl
      write (lupt,30) runname,fl
      write (lusp,30) runname,fl
      write (lupc,30) runname,fl
      if (bandsmod.eq.'Y') write (lupb,30) runname,fl
      if (ratmod.eq.'Y') write (lurt,30) runname,fl
      if (dynmod.eq.'Y') write (ludy,30) runname,fl
      if (allmod.eq.'Y') write (lual,30) runname,fl
      if (fclmod.eq.'Y') write (lucl,30) runname,fl
c
c     Write Model Parameters
c
   40 format(//,' Model Parameters:',/,
     &        ' =================',//,
     &        ' Abundances     : ',a128,/,
     &        ' Pre-ionisation : ',a128,/,
     &        ' Photon Source  : ',a128)
   50 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12)
   60 format(/,' Charge Exchange: ',a12,/,
     &        ' Photon Mode    : ',a12,/,
     &        ' Collision calcs: ',a12,/,
     &        ' Electron Kappa : ',1pg11.4)
   70 format(//' ',t2,'Jden',t9,'Jgeo',t16,'Jtrans',
     &            t23,'Jend',t29,'Ielen',t35,
     &         'Jpoen',t43,'Fren',t51,'Tend',t57,'DIend',t67,
     &         'TAUen',t77,'Jeq',t84,'Teini')
   80 format('  ',4(a4,3x),2(i2,4x),0pf6.4,0pf6.0,
     &           2(1pg10.3),3x,a4,0pf7.1)
   90 format(//' Photon Source'/
     &         ' ============='/)
  100 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
  110 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
C
      call protostate (6)
      call protostate (luop)
      call protostate (lupt)
      if (bandsmod.eq.'Y') call protostate (lupb)

      if (ratmod.eq.'Y') call protostate (lurt)
      if (dynmod.eq.'Y') call protostate (ludy)
      if (allmod.eq.'Y') call protostate (lual)
      if (fclmod.eq.'Y') call protostate (lucl)
c
      write (luop,40) abnfile,ionsetup,srcfile
      write (lupt,40) abnfile,ionsetup,srcfile
      write (lusp,40) abnfile,ionsetup,srcfile
      write (lupc,40) abnfile,ionsetup,srcfile
c
      if (bandsmod.eq.'Y') write (lupb,40) abnfile,ionsetup,srcfile
      if (ratmod.eq.'Y') write (lurt,40) abnfile,ionsetup,srcfile
      if (dynmod.eq.'Y') write (ludy,40) abnfile,ionsetup,srcfile
      if (allmod.eq.'Y') write (lual,40) abnfile,ionsetup,srcfile
      if (fclmod.eq.'Y') write (lucl,40) abnfile,ionsetup,srcfile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     abundances file header
c
c
      abundtitle=' Initial Abundances :'
      do i=1,atypes
        zi(i)=zion0(i)*deltazion(i)
      enddo
      call dispabundances (luop, zi, abundtitle)
      call dispabundances (lupt, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
      call dispabundances (lupt, zion, abundtitle)
c
      call wabund (lusp)
      call wabund (lupc)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      cht='New Set'
      if (chargemode.eq.2) cht='Disabled'
      if (chargemode.eq.1) cht='Old Set'
      pht='Normal'
      if (photonmode.eq.0) pht='Zero Field, Decoupled'
      clt='Dere 2007 Collisions'
      if (collmode.eq.1) clt='Old A&R Collision Methods'
      if (collmode.eq.1) clt='LMS Collision Methods'
c
      if (usekappa) then
        write (luop,60) cht,pht,clt,kappa
        write (lupt,60) cht,pht,clt,kappa
      else
        write (luop,50) cht,pht,clt
        write (lupt,50) cht,pht,clt
      endif
      write (luop,70)
      write (luop,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,diend,
     &tauen,tmod,tm00
      write (luop,90)
      write (luop,100)
      write (luop,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      write (lupt,70)
      write (lupt,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,diend,
     &tauen,tmod,tm00
      write (lupt,90)
      write (lupt,100)
      write (lupt,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (bandsmod.eq.'Y') then
        write (lupb,50) cht,pht,clt
        write (lupb,70)
        write (lupb,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lupb,90)
        write (lupb,100)
        write (lupb,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (ratmod.eq.'Y') then
        write (lurt,50) cht,pht,clt
        write (lurt,70)
        write (lurt,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lurt,90)
        write (lurt,100)
        write (lurt,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (dynmod.eq.'Y') then
        write (ludy,50) cht,pht,clt
        write (ludy,70)
        write (ludy,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (ludy,90)
        write (ludy,100)
        write (ludy,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
      if (allmod.eq.'Y') then
c
        write (lual,50) cht,pht,clt
        write (lual,70)
        write (lual,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lual,90)
        write (lual,100)
        write (lual,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
      if (fclmod.eq.'Y') then
c
        write (lucl,50) cht,pht,clt
        write (lucl,70)
        write (lucl,80) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lucl,90)
        write (lucl,100)
        write (lucl,110) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      endif
c
c     Shock parameters
c
  120 format(//,' Initial Jump Conditions:',/,
     &         ' =========================')
      write (luop,120)
      write (lupt,120)
      if (bandsmod.eq.'Y') write (lupb,120)
      if (ratmod.eq.'Y') write (lurt,120)
      if (dynmod.eq.'Y') write (ludy,120)
      if (allmod.eq.'Y') write (lual,120)
      if (fclmod.eq.'Y') write (lucl,120)
c
  130 format(//'  T0',1pg14.7,' V0',1pg14.7,/,
     &         ' RH0',1pg14.7,' P0',1pg14.7,' B0',1pg14.7,//,
     &         '  T1',1pg14.7,' V1',1pg14.7,/,
     &         ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7)
c
      dr=0.d0
      dv=vel1-vel0
c
      write (*,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &rho1,pr1,bm1*1.d6
      write (luop,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &rho1,pr1,bm1*1.d6
      write (lupt,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &rho1,pr1,bm1*1.d6
      if (bandsmod.eq.'Y') then
        write (lupb,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-
     &   5,rho1,pr1,bm1*1.d6
      endif
      if (ratmod.eq.'Y') then
        write (lurt,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-
     &   5,rho1,pr1,bm1*1.d6
      endif
      if (dynmod.eq.'Y') then
        write (ludy,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-
     &   5,rho1,pr1,bm1*1.d6
      endif
      if (fclmod.eq.'Y') then
        write (lucl,130) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-
     &   5,rho1,pr1,bm1*1.d6
      endif
c
      close (lusp)
      close (lupc)
c
      wdilt0=te0
      t=te0
      ve=vel0
      dr=1.0
      dv=ve/100.d0
      fi=1.d0
      rad=1.d38
      if (wdil.eq.0.5) rad=0.d0
c
c     get the electrons...
c
      de=feldens(dh,pop)
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
c
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
      endif
      call cool (t, de, dh)
c
c
  140 format(//,' Precursor Conditions and Ionisation State',/
     &          ' =========================================',/)
      write (luop,140)
      write (lupt,140)
      if (bandsmod.eq.'Y') write (lupb,140)
      if (fclmod.eq.'Y') write (lucl,140)
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
      wmod='FILE'
      call wmodel (luop, t, de, dh, dr, wmod)
      call wmodel (lupt, t, de, dh, dr, wmod)
      if (bandsmod.eq.'Y') call wmodel (lupb, t, de, dh, dr, wmod)
      if (ratmod.eq.'Y') then
        call wmodel (lurt, t, de, dh, dr, wmod)
  150 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,'  en',t70,
     &     '  Dloss',t84,'  Eloss',t98,'  Egain',t112,'  Coll. H',t126,
     &     '  Chrg. Ex.',t140,'  Reson',t154,
     &     '  Xreson',t168,'  Inter/fine',t182,'  He int-forb',t196,
     &     '  Forbid',t210,'  Fe II',t224,'  2Photon',t238,
     &     '  Compton',t252,'  Free-Free',t266,'  Coll. Ion',t280,
     &     '  Photo.',t294,'  Recomb.',t308,'  Cosmic',t322,
     &     '  GrainsH-C',t336,'  G.Heat',t350,'  G.Cool',t364,
     &     '  PAHs')
        write (lurt,150)
      endif
c
      call wionabal (luop, pop)
      call wionabal (lupt, pop)
      if (bandsmod.eq.'Y') call wionabal (lupb, pop)
      if (allmod.eq.'Y') call wionabal (lual, pop)
      if (fclmod.eq.'Y') call wionabal (lucl, pop)
c
      close (luop)
      close (lupt)
  160 format(17(a12,a1))
      if (bandsmod.eq.'Y') then
        write (lupb,160) 'Te',tab,'de',tab,'dh',tab,'en',tab,'FHI',tab,'
     &FHII',tab,'mu',tab,'tloss',tab,'Lambda',tab,'ff/total',tab,'B0.0-0
     &.1keV',tab,'B0.1-0.5keV',tab,'B0.5-1.0keV',tab,'B1.0-2.0eV',tab,'B
     &2.0-10.0keV',tab,'Ball'
        close (lupb)
      endif
  170 format(//'Mean Zone Values'/
     &         '================'/)
  180 format(10(a12,a2),30(a12,a2),33(a12,a2))
  190 format(10(a12,a2),30(a12,a2),33(1pg12.5,a2))
c
  200 format( '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '=======================================================',
     &        '======================================')
c
      if (fclmod.eq.'Y') then
        write (lucl,170)
c      write (lucl,560)
c      write (lucl,570)
        if (jnorm.eq.0) then
          write (lucl,180) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &.nH)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLy
     &alpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.1) then
          write (lucl,180) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(nH
     &^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLya
     &lpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.2) then
          write (lucl,180) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &.ni)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLy
     &alpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.3) then
          write (lucl,180) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(n^
     &2)',tab,(elem(j),tab,j=1,atypes),'Eloss/n2',tab,'Egain/n2'
     &   ,tab,'Netloss/n2'
c     &2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLyal
c     &pha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        if (jnorm.eq.4) then
          write (lucl,180) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'FHI   ',tab,'FHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLya
     &lpha'
c     &    (elem(j),tab,j=1,atypes)
        endif
        write (lucl,190) '(K)',tab,'(/cm^3)',tab,'(/cm^3)',tab,'(/cm^3)'
     &   ,tab,'(g/cm^3)',tab,' ',tab,' ',tab,'(amu)',tab,'(erg/cm^3/s)',
     &   tab,'(erg cm^3/s)',tab,('(erg cm^3/s)',tab,j=1,atypes),'(erg cm
     &^3/s)',tab,'(erg cm^3/s)',tab,'(erg cm^3/s)'
c     &    (elem(j),tab,j=1,atypes)
        write (lucl,200)
      endif
      if (ratmod.eq.'Y') close (lurt)
      if (dynmod.eq.'Y') close (ludy)
      if (allmod.eq.'Y') close (lual)
      if (fclmod.eq.'Y') close (lucl)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fieldsummary ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i
      real*8 blum,ilum,qoii,xlum,qall
      real*8 q1,q2,q3,q4
c
      blum=0.d0
      ilum=0.d0
      qoii=0.d0
      xlum=0.d0
      qall=0.d0
c
      do i=1,infph-1
        wid=photev(i+1)-photev(i)
        blum=blum+tphot(i)*wid*evplk
        qall=qall+tphot(i)*wid*evplk/(cphotev(i)*ev)
        if (photev(i).ge.iph) ilum=ilum+tphot(i)*wid*evplk
        if (photev(i).ge.100.d0) xlum=xlum+tphot(i)*wid*evplk
        if (photev(i).ge.3.51211d+01) then
          qoii=qoii+tphot(i)*wid*evplk/(cphotev(i)*ev)
        endif
      enddo
c
      blum=fpi*blum
      ilum=fpi*ilum
      xlum=fpi*xlum
      qoii=fpi*qoii
      qall=fpi*qall
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    Assuming plane parallel geometry and dilution of 0.5 and intvec
c    assumes 1/4pi units and soupho is 1/pi plane parallel units.
c    intvec is usually called on tphot which is 1/4pi
c
      q1=0.0d0
      q2=0.0d0
      q3=0.0d0
      q4=0.0d0
c
      call intvec (tphot, q1, q2, q3, q4)
c
      qhi=q1!*0.25d0
      qhei=q2!*0.25d0
      qheii=q3!*0.25d0
      qht=q4!*0.25d0
c
      write (*,10) blum,qall,ilum,qht,xlum,qhi,qhei,qheii,qoii
      open (luop,file=fl,status='OLD',access='APPEND')
      write (lupt,10) blum,qall,ilum,qht,xlum,qhi,qhei,qheii,qoii
      close (luop)
c      write (luop,10) blum,qall,ilum,qht,xlum,qhi,qhei,qheii,qoii
      write (lupt,10) blum,qall,ilum,qht,xlum,qhi,qhei,qheii,qoii
   10 format(/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Shock  photoionisation source: '/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & ' : Total intensity   : ',1pg12.5,'   (ergs/s/cm^2)      :'/
     & ' : Total Photons     : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' : Ion.  intensity   : ',1pg12.5,'   (ergs/s/cm^2)      :'/
     & ' : FQ tot   (>1Ryd)  : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' : X-Ray >0.1keV int.: ',1pg12.5,'   (ergs/s/cm^2)      :'/
     & ' : FQHI  (1-1.8Ryd)  : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' : FQHeI (1.8-4Ryd)  : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' : FQHeII   (>4Ryd)  : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' : FQOII  (>35.12eV) : ',1pg12.5,'   (phots/cm^2/s)     :'/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine protostate (lunit)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 fmua
      integer*4 lunit
c
   10   format(//
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '  Shock Proto-State Properties',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Velocity       :',1pg12.5,' km/s',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Proto Mach Number          :',1pg12.5,/,
     &  '    Proto Alfven Mach Number   :',1pg12.5,/,
     &  '    Proto Mag Alpha (Pmag/Pgas):',1pg12.5,/,
     &  '    Proto Gas Eta  (gPgas/Pram):',1pg12.5,/,
     &  '    Proto Mag Eta  (2Pmag/Pram):',1pg12.5,//,
     &  '    Proto        T :',1pg12.5,' K',/,
     &  '    Proto        ne:',1pg12.5,' cm^-3',/,
     &  '    Proto        nH:',1pg12.5,' cm^-3',/,
     &  '    Proto        d :',1pg12.5,' g/cm^-3',/,
     &  '    Proto      Pgas:',1pg12.5,' dyne/cm^2',/,
     &  '    Proto        mu:',1pg12.5,' a.m.u.',/,
     &  '    Proto       FHI:',1pg12.5,/,
     &  '    Proto      FHII:',1pg12.5,/,
     &  '    Proto      FHeI:',1pg12.5,/,
     &  '    Proto     FHeII:',1pg12.5,/,
     &  '    Proto    FHeIII:',1pg12.5,/,
     &  '    Proto        B :',1pg12.5,' microGauss',/,
     &  '    Proto      Pmag:',1pg12.5,' dyne/cm^2',/,
     &  '    Proto      Pram:',1pg12.5,' dyne/cm^2',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c
      en=zen*dh_neu
      cspd=dsqrt(gammaeos*pr_neu/rh_neu)
      wmol=rh_neu/(en+de_neu)
      mu=fmua(de_neu,dh_neu)
c
      pram=rh_neu*vel0*vel0
      pgas=pr_neu
      pmag=(bm_neu*bm_neu)/epi
c
      machnumber=vel0/dsqrt(gammaeos*pr_neu/rh_neu)
      alfvennumber=vel0/dsqrt(2.0d0*pmag/rh_neu)
      malpha=pmag/pgas
      gaseta=gammaeos*pgas/pram
      mageta=2.d0*pmag/pram
c
      write (lunit,10) vel0*1.0d-5,machnumber,alfvennumber,
     & malpha,gaseta,mageta,
     & te_neu,de_neu,dh_neu,rh_neu,
     & pr_neu,mu,pop_neu(1,zmap(1)),pop_neu(2,zmap(1)),
     & pop_neu(1,zmap(2)),pop_neu(2,zmap(2)),pop_neu(3,zmap(2)),
     & bm_neu*1.d6,pmag,pram
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shocksummary (lunit)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 fmua,fpresse,frho
      integer*4 lunit
c
c uses global pop
c
c
c pre = preshock , pst = postshock
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0 !shock velocity
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      cmpf=rho1/rho0
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
c
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
      pgas=pr0
      pmag=bp0
      pram=rh_pre*vs_pre*vs_pre
c
      machnumber=vel0/dsqrt(gammaeos*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.d0*pmag/rho0)
      malpha=pmag/pgas
      gaseta=gammaeos*pgas/pram
      mageta=2.d0*pmag/pram
c
   11   format(
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Velocity       :',1pg12.5,' km/s',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Preshock Mach Number          :',1pg12.5,/,
     &  '    Preshock Alfven Mach Number   :',1pg12.5,/,
     &  '    Preshock Mag Alpha (Pmag/Pgas):',1pg12.5,/,
     &  '    Preshock Gas Eta  (gPgas/Pram):',1pg12.5,/,
     &  '    Preshock Mag Eta  (2Pmag/Pram):',1pg12.5,//,
     &  '    Preshock     T :',1pg12.5,' K',/,
     &  '    Preshock     ne:',1pg12.5,' cm^-3',/,
     &  '    Preshock     nH:',1pg12.5,' cm^-3',/,
     &  '    Preshock     d :',1pg12.5,' g/cm^-3',/,
     &  '    Preshock   Pgas:',1pg12.5,' dyne/cm^2',/,
     &  '    Preshock     mu:',1pg12.5,' a.m.u.',/,
     &  '    Preshock    FHI:',1pg12.5,/,
     &  '    Preshock   FHII:',1pg12.5,/,
     &  '    Preshock   FHeI:',1pg12.5,/,
     &  '    Preshock  FHeII:',1pg12.5,/,
     &  '    Preshock FHeIII:',1pg12.5,/,
     &  '    Preshock     B :',1pg12.5,' microGauss',/,
     &  '    Preshock   Pmag:',1pg12.5,' dyne/cm^2',/,
     &  '    Preshock   Pram:',1pg12.5,' dyne/cm^2',/)
   12   format(
     &  '    Postshock Compression Factor   :',1pg12.5,/,
     &  '    Postshock    T :',1pg12.5,' K',/,
     &  '    Postshock    ne:',1pg12.5,' cm^-3',/,
     &  '    Postshock    nH:',1pg12.5,' cm^-3',/,
     &  '    Postshock    v :',1pg12.5,' km/s',/,
     &  '    Postshock    d :',1pg12.5,' g/cm^-3',/,
     &  '    Postshock  Pgas:',1pg12.5,' dyne/cm^2',/,
     &  '    Postshock    B :',1pg12.5,' microGauss',/,
     &  '    Postshock  Pmag:',1pg12.5,' dyne/cm^2',/,
     &  '    Postshock  Pram:',1pg12.5,' dyne/cm^2',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
      en=zen*dh0
      pr0=fpresse(te0,de0,dh0)
      rho0=frho(de0,dh0)
      cspd=dsqrt(gammaeos*pr0/rho0)
      wmol=rho0/(en+de0)
      mu=fmua(de0,dh0)
C
      pgas=pr0
      pmag=(bm0*bm0)/epi
      pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaeos*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.0d0*pmag/rho0)
      malpha=pmag/pgas
      gaseta=gammaeos*pgas/pram
      mageta=2.d0*pmag/pram
C
      write (lunit,11) vel0*1.0d-5,machnumber,alfvennumber,
     & malpha,gaseta,mageta,
     & te0,de0,dh0,rho0,
     & pr0,mu,pop(1,zmap(1)),pop(2,zmap(1)),
     & pop(1,zmap(2)),pop(2,zmap(2)),pop(3,zmap(2)),
     & bm0*1.d6,bp0,pram
C
      pram=rho1*vel1*vel1
      write (lunit,12) cmpf,te1,de1,dh1,vel1*1d-5,
     & rho1,pr1,bm1*1.d6,bp1,pram
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock5jump ()
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 tj,dej,dhj,tpo,va,eps
      real*8 feldens,fpressu,frho,velshock2
c
c     Now get revised shock solution for changed pre shock values:
c     Shock in term of flow velocity
c
      if (stype.eq.'v') then
c
        dr=1.d0
c
        if (mtype.eq.'b') then
          bmag=magparam*1.0d-6
          pmag=(bmag*bmag)/epi
        endif
c
        if (mtype.eq.'a') then
          pgas=fpressu(te_neu,dh_neu,pop_neu)
          pmag=magparam*pgas
          bmag=dsqrt(epi*pmag)
        endif
c
        if (mtype.eq.'c') then
          pgas=fpressu(te_pre,dh_pre,pop_pre)
          pmag=magparam*pgas
          bmag=dsqrt(epi*pmag)
        endif
c
        rh_pre=frho(de_pre,dh_pre)
        rh_neu=frho(de_neu,dh_neu)
c
        pram=rh_pre*vs_pre*vs_pre
c
        if (mtype.eq.'m') then
          va=vs_pre/mmach
          Pmag=(0.5*rh_pre*va*va)
          bmag=dsqrt(epi*pmag)
        endif
c
        if (mtype.eq.'r') then
          pmag=0.5d0*mageta*pram
          bmag=dsqrt(epi*pmag)
        endif
c
        bm_pre=bmag
c
        pgas=fpressu(te_pre,dh_pre,pop_pre)
        vshoc=vs_pre
        te0=te_pre
        tj=te_pre
        dhj=dh_pre
        dej=de_pre
        tloss=0.0d0
        dt=0.0d0
c
        rmod='NEAR'
c        call isochorflow (t, de, dh, vshoc, bmag, rmod, tloss, dt)
c        call isobarflow (t, de, dh, vshoc, bmag, rmod, tloss, dt)
        call rankhug (tj, dej, dhj, vshoc, bmag, tloss, dt)
        tm00=te0
      endif
c
c     Shock in term of temperature
c
      if (stype.eq.'t') then
c
        dr=1.d0
c
        if (mtype.eq.'b') then
          bmag=magparam*1.0d-6
          pmag=(bmag*bmag)/epi
        endif
c
        if (mtype.eq.'a') then
          pgas=fpressu(te_neu,dh_neu,pop_neu)
          pmag=magparam*pgas
          bmag=dsqrt(epi*pmag)
        endif
c
        if (mtype.eq.'c') then
          pgas=fpressu(te_pre,dh_pre,pop_pre)
          pmag=magparam*pgas
          bmag=dsqrt(epi*pmag)
        endif
c
        rh_neu=frho(de_neu,dh_neu)
        pram=rh_neu*vs_pre*vs_pre
c
        if (mtype.eq.'m') then
          va=vs_pre/mmach
          Pmag=(0.5*rh_pre*va*va)
          bmag=dsqrt(epi*pmag)
        endif
c
        if (mtype.eq.'r') then
          pmag=0.5d0*mageta*pram
          bmag=dsqrt(epi*pmag)
        endif
c
        bm_pre=bmag
c
        pgas=fpressu(te_pre,dh_pre,pop_pre)
        vshoc=vs_pre
        t=te_pre
        tpo=te_pst
        dh=dh_pre
        de=feldens(dh,pop)
        tloss=0.0d0
        dt=0.0d0
c
  100   bm0=bmag
c
        vshoc= velshock2 (dh, t, bmag, tpo)
c
        pram=rh_neu*vshoc*vshoc
c
        if (mtype.eq.'m') then
          va=vshoc/mmach
          Pmag=(0.5*rh_pre*va*va)
          bmag=dsqrt(epi*pmag)
          eps=dabs(2.0*(bmag-bm0)/(bmag+bm0))
          if (eps.gt.1.0d-6 ) goto 100
        endif
c
        if (mtype.eq.'r') then
c
c iterate for rampressure/alpha_r
c
          pmag=0.5d0*mageta*pram
          bmag=dsqrt(epi*pmag)
          eps=dabs(2.0*(bmag-bm0)/(bmag+bm0))
          if (eps.gt.1.0d-6 ) goto 100
        endif
c
        vs_neu=vel0
        bm_neu=bmag
        tm00=te0
      endif
c
c Shock intial state kept for iterations, above numbers are
c reused from step to step
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
c
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
      pmag=bp0
      pgas=pr0
      pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaeos*pr0/rho0)
      malpha=pmag/pgas
      gaseta=gammaeos*pgas/pram
      mageta=2.d0*pmag/pram
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock5check (its, maxits)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 its,maxits
   10 format(/,
     & ' *********************************************************',/,
     & '  SHOCK 5  Convergence Test, It.: ',i2.2,' of ',i2.2)
   20 format(
     & '  Result: CONVERGED',/,
     & ' *********************************************************',/)
   30 format(
     & '  Result: NOT CONVERGED',/,
     & ' *********************************************************',/)
   40 format(
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Psi (Q/v): ',  1pg11.4,'  `Psi      : ',  1pg11.4,/,
     & '    Compress : ',  1pg11.4,'  `Compress : ',  1pg11.4,/,
     & '    T_pre    : ',  1pg11.4,'  `T_pre    : ',  1pg11.4,/,
     & '    T_shock  : ',  1pg11.4,'  `T_shock  : ',  1pg11.4,/,
     & '    ne_pre   : ',  1pg11.4,'  `ne_pre   : ',  1pg11.4,/,
     & '    DelH/He  : ',  1pg11.4,'%',/,
     & '    RMS      : ',  1pg11.4,'%',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
      real*8 delhhe,rmserr,term
c
      converged=0
      if (its.le.1) then
c
c save current shock for next test
c
      call copypop (pop_pre, pop_pre0)
      psi0=psi
      te_pre0=te_pre
      te_pst0=te_pst
      de_pre0=de_pre
      de_pst0=de_pst
      dh_pre0=dh_pre
      dh_pst0=dh_pst
      vs_pre0=vs_pre
      vs_pst0=vs_pst
c
      rh_pre0=rh_pre
      rh_pst0=rh_pst
      pr_pre0=pr_pre
      pr_pst0=pr_pst
      bm_pre0=bm_pre
      bm_pst0=bm_pst
c
      psi0=psi
      cmpf0=cmpf
      return
c
      else
c
c uses global pop
c
      call difhhe (pop_pre, pop_pre0, delhhe)
      write (*,10) its,maxits
      term=2.d0*(cmpf-cmpf0)/(cmpf+cmpf0)
      rmserr=term*term
      term=2.d0*(te_pre-te_pre0)/(te_pre+te_pre0)
      rmserr=rmserr+(term*term)
      term=2.d0*(te_pst-te_pst0)/(te_pst+te_pst0)
      rmserr=rmserr+(term*term)
      term=2.d0*(de_pre-de_pre0)/(de_pre+de_pre0)
      rmserr=rmserr+(term*term)
      rmserr=rmserr+(delhhe*delhhe)
      rmserr=dsqrt(rmserr/6.d0)
      write (*,40) psi,psi0,cmpf,cmpf0,te_pre,te_pre0,te_pst,te_pst0,
     &de_pre,de_pre0,delhhe*100.d0,rmserr*100.d0
      if (rmserr.lt.1e-3) then
        converged=1
        write (*,20)
      else
        converged=0
        if (its.ge.maxits) maxits=maxits+1
        write (*,30)
      endif
c
c save current shock for next test
c
      call copypop (pop_pre, pop_pre0)
      psi0=psi
      te_pre0=te_pre
      te_pst0=te_pst
      de_pre0=de_pre
      de_pst0=de_pst
      dh_pre0=dh_pre
      dh_pst0=dh_pst
      vs_pre0=vs_pre
      vs_pst0=vs_pst
c
      rh_pre0=rh_pre
      rh_pst0=rh_pst
      pr_pre0=pr_pre
      pr_pst0=pr_pst
      bm_pre0=bm_pre
      bm_pst0=bm_pst
c
      psi0=psi
      cmpf0=cmpf
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine shock5precursor (iteration, maxits)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      logical iexi
      integer*4 iteration,maxits
      integer*4 itcount
      real*8 tf
      character nmod*4,ratemode*4
      character tab*1
c
      integer*4 i,idx,nfs,n100,nrms,t0lim
      real*8 dtime,dtimer,pre_par,qh
      real*8 drr,tstep,absf,rms,invnfs
      real*8 diff,xhfinal
      real*8 rmserr,rmslimit,te_0,de_0,dh_0
c front structure arrays
      real*8 frdr
      real*8 f(mxifsteps),att(mxifsteps)
      real*8 x(mxifsteps),nte(mxifsteps),ne(mxifsteps)
      real*8 nh(mxifsteps),fh(mxifsteps)
      real*8 foi(mxifsteps),foii(mxifsteps),foiii(mxifsteps)
      real*8 p1(mxion, mxelem)
      real*8 p2(mxion, mxelem)
      real*8 popfr(mxifsteps, mxion, mxelem)
      real*8 popintfr(mxifsteps, mxion, mxelem)
      real*8 src(mxinfph), attcol(mxinfph), sigcol(mxinfph)
c
      real*8 feldens,frectim3
      real*8 densnum
c
c      open (luop,file=fl,status='OLD',access='APPEND')
      open (lupt,file=fpt,status='OLD',access='APPEND')
   10 format(
     & ' *********************************************************',/,
     & '  SHOCK 5 Precursor Iteration: ',i2.2,' of ',i2.2)
      write (*,10) ,iteration,maxits
      if (finalit.gt.0) then
        write (lupt,10) ,iteration,maxits
        open (lupc,file=fpc,status='OLD',access='APPEND')
        write (lupc,10) ,iteration,maxits
        close (lupc)
      endif
c
      tab=','
      jspot='NO'
      jcon='YES'
      ratemode='ALL'
c
      rmslimit=5.0d-2
      if (iteration.gt.1) rmslimit=2.0d-2
      if (finalit.gt.0)   rmslimit=5.0d-3
c
c make final upstream field and save
c (compare with one created in main loop)
c
      fi=1.0d0
      wdil=0.5d0
c
      t=te(step)
      de=deel(step)
      dh=dhy(step)
      bmag=bmg(step)
c
      rad=dist(step)
      dr=dist(step)-dist(step-1)
      dv=veloc(step)-veloc(step-1)
c
      call localem (t, de, dh)
      specmode='UP'
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
c debugging purposes
c      caller='S5'
c      pfx='prec_up'
c      np=7
cc
c      wmod='REAL'
c      call wpsou (caller, pfx, np, wmod, t, de, dh, veloc(step), dv
c     &           , vshoc, rad, dr, timlps(step), dt, wdil, tphot)
c
      do i=1,infph
        src(i)=tphot(i)
        tphot(i)=src(i)*wdil
      enddo
c
c  lets see what we got:
c
      call fieldsummary
c
c   compare ion front velocity to shock velocity
c
      ve=vshoc
      viofr=qht/densnum(dh_neu) ! all ions
      qh=qht/dh_neu ! just nH
      pre_par=viofr/ve
c
   20 format(/,'  SHOCK 5 Initial Precursor Parameter: ',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '     v_shock       : ',  1pg12.5,' km/s',/,
     & '     Q_ions        : ',  1pg12.5,' km/s',/,
     & '     Psi (Q/v)     : ',  1pg12.5,/,
     & '     Q_H           : ',  1pg12.5,' km/s',/,
     & '     Psi_H (Q_H/v) : ',  1pg12.5,/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,20)
     &  vshoc*1.d-5,viofr*1.d-5,viofr/vshoc,
     &  qh*1.d-5,qh/vshoc
      if (finalit.gt.0) then
      write (lupt,20)
     &  vshoc*1.d-5,viofr*1.d-5,viofr/vshoc,
     &  qh*1.d-5,qh/vshoc
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c setup for precursor iteration, fill arrays with protoionisation
c and opacities
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c to allow outer edge of ion front to evolve in fast ion front case,
c but reset to neutral for each global interation
c
      te_0=te_neu
      de_0=de_neu
      dh_0=dh_neu
c
      call copypop (pop_pre, pop)
c
      t=te_pre
      de=de_pre
      dh=dh_pre
c
      call allrates (t, ratemode)
c
c  stromgren length if front detatches, recomb rate if ionised gas
c
      dtimer=frectim3(dh)*pre_par
      drr=dtimer*ve
c
      call copypop (pop_neu, pop)
c
      t=te_0
      de=de_0
      dh=dh_0
      ve=vshoc
c
      call allrates (t, ratemode)
c
c  go nfs* distance that absorbs 10% of the field in the neutral medium
c  tau = ~ 10.38 for 100 steps ~6.5 for 50 steps
c
c set up absorbsion fractions and attenuations arrays
c setup ion balance arrays too
c
      nfs=mxifsteps
      invnfs=1.d0/dble(nfs)
      absf=0.10d0
      call absdis2 (dh, fi, absf, dr, 0.0d0, pop_neu)
      frdr=nfs*dr
c
c
c index 1 is zone at the source, index nfs is outer edge of front zone
c
      call copypop (pop_neu, p1)
      call copypopstep (p1, 1, popfr)
      call clearpop (p1)
      call copypopstep (p1, 1, popintfr)
      f(1)=absf
      fh(1)=pop_neu(2,1)
      foi(1)=pop_neu(1,zmap(8))
      foii(1)=pop_neu(2,zmap(8))
      foiii(1)=pop_neu(3,zmap(8))
      att(1)=1.0d0
      x(1)=0.d0
      nh(1)=dh_0
      ne(1)=de_0
      nte(1)=te_0
      do i=2,nfs
        call copypop (pop_neu, p1)
        call copypopstep (p1, i, popfr)
        call scalepop (p1, dr)
        call copysteppop (i-1, popintfr, p2)
        call addpop (p1, p2)
        call copypopstep (p2, i, popintfr)
        fh(i)=pop_neu(2,1)
        foi(i)=pop_neu(1,zmap(8))
        foii(i)=pop_neu(2,zmap(8))
        foiii(i)=pop_neu(3,zmap(8))
        f(i)=absf
        att(i)=att(i-1)*(1.d0-f(i-1))
        x(i)=(i-1)*dr
        nh(i)=dh_0
        ne(i)=de_0
        nte(i)=te_0
      enddo
c
      dtime=frdr/vshoc
      nfs=nint(dlog(1.d-3)/dlog(1.d0-absf))
      nfs=min(nfs,mxifsteps)
c
   30 format(/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Precursor Timescales: '/
     & '    : ',i4,' zone neutral optical depth   : ',1pg11.4,/
     & '    : ',i4,' zone neutral crossing time   : ',1pg11.4,/
     & '    : ',i4,' zone neutral length          : ',1pg11.4,//
     & '    Ionised recombination time x Q/v : ',1pg11.4,/
     & '    Ionised recombination distance   : ',1pg11.4,/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,30) nfs,-dlog(att(nfs)),
     &             nfs,dtime,
     &             nfs,frdr,dtimer,drr
      if (finalit.gt.0) then
C       write (luop,30) nfs,-dlog(att(nfs)),
C    &                  nfs,dtime,
C    &                  nfs,frdr,dtimer,drr
        write (lupt,30) nfs,-dlog(att(nfs)),
     &                  nfs,dtime,
     &                  nfs,frdr,dtimer,drr
      endif
c
      dr=frdr*invnfs
      tstep=(frdr/vshoc)*invnfs
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Iterate, full solution for Psi < 1, inner edge only for Psi >1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   40 format(/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & '  Precursor Flow: ',/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write (*,40)
      if (finalit.gt.0) then
c        write (luop,40)
        write (lupt,40)
      endif
c
      xhfinal=pop(1,1)
      itcount=0
      t0lim=0
c
   50 call copypop (pop_neu, pop)
c
      t=te_0
      de=de_0
      dh=dh_0
c
      write (*,'(" Initial:     ",3(1pg11.4))') t,de,pop(1,1)
      if (finalit.gt.0) then
c        write (luop,'(" Initial:     ",3(1pg11.4))') t,de,pop(1,1)
        write (lupt,'(" Initial:     ",3(1pg11.4))') t,de,pop(1,1)
      endif
      nmod='TIM'
      rms=0.d0
c
c integrate from outer edge to source
c
      rmserr=0.d0
      nrms=min(max(3,nfs/3),30)
      do i=nfs,1,-1
c
c full column absorption applied to spectrum
c
        call copysteppop (i, popintfr, p1)!getinnerfractioncolumns
        call attsig (dh, fi, 0.d0, p1, attcol, sigcol)
        do idx=1,infph
          tphot(idx)=src(idx)*attcol(idx)*wdil ! +prefield(i)*wdil
          skipbin(idx)=.false.
          ipho=1
          iphom=0
          tem=0.d0
        enddo
c redo rates
        call allrates (t, ratemode)
c evolve in attenuated field to inner boundary
        call copypop (pop, p1)
        call teequi2 (t, tf, de, dh, tstep, nmod)
        call copypop (pop, p2)
        call averinto (0.5d0, p1, p2, p1)
c save step state at point, set nte after comparing with tf below
        x(i)=(i-1)*dr
        nh(i)=dh
        ne(i)=feldens(dh,p1)
c find new abs fraction in current zone for next iteration
        call habsfrac(absf, dh, fi, dr, p1)
        f(i)=absf
C
        if (i.le.nrms) then
          call copysteppop (i, popfr, p1)
          call difhhe (pop, p1, diff)
          diff=dabs(1.d0-(10.d0**diff))
          rmserr=rmserr+(diff*diff)
          diff=0.5d0*(tf-nte(i))/(tf+nte(i))
          rmserr=rmserr+(diff*diff)
        endif
c save step temp at point
        nte(i)=tf
c save pop after step = closer to src, ith pop is inner edge of ith zone
        call copypopstep (pop, i, popfr)
        fh(i)=pop(2,1)
        foi(i)=pop(1,zmap(8))
        foii(i)=pop(2,zmap(8))
        foiii(i)=pop(3,zmap(8))
        de=feldens(dh,pop)
        t=tf
C       write (*,'(x,5(1pg11.5,x))') (nfs-i+1)*tstep*invnfs,(i-1)*dr,tf,
C    &   de,pop(1,1)
C       if (finalit.gt.0) then
C         write (lupt,'(x,5(1pg11.5,x))') (nfs-i+1)*tstep*invnfs,(i-1)*
C    &     dr,tf,de,pop(1,1)
C        endif
      enddo
c
      att(1)=1.0d0
      do i=2,nfs
        att(i)=att(i-1)*(1.d0-f(i-1))
      enddo
c
C  65  format(/
C    & ' step  Time (s)   Dist (cm)   Temp. (K)  n_e (cm^-3)',
C    & '   FHII         <abs>       <att>',/,
C    & ' :::::::::::::::::::::::::::::::::::::::::::::::::::',
C    & ':::::::::::::::::::::::::::::::::::::')
   65 format(//,a4,12(a1,2x,a10,1x),/
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ':::::::::::::::::::::::::::::::::::::')
   66 format(i5.3,12(a1,1pg13.6))
        write (*,65) ' Step',tab,'Dist.(cm)',tab,'Time (s)',
     &   tab, 'Te (K)',tab,'n_e(cm^-3)',tab,'n_H(cm^-3)',
     &   tab,'FHI',tab,'FHII',tab,'<abs>',tab,'<att>',
     &   tab,'FOI',tab,'FOII',tab,'FOIII'
c
      if (finalit.gt.0) then
        write (lupt,65) ' Step',tab,'Dist.(cm)',tab,'Time (s)',
     &   tab, 'Te (K)',tab,'n_e(cm^-3)',tab,'n_H(cm^-3)',
     &   tab,'FHI',tab,'FHII',tab,'<abs>',tab,'<att>',
     &   tab,'FOI',tab,'FOII',tab,'FOIII'
      endif
      do i=nfs,1,-1
        write (*,66) i,tab,x(i),tab,(nfs-i+1)*tstep*invnfs,
     &    tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &    tab,fh(i),tab,f(i),tab,att(i),
     &    tab,foi(i),tab,foii(i),tab,foiii(i)
        if (finalit.gt.0) then
        write (lupt,66) i,tab,x(i),tab,(nfs-i+1)*tstep*invnfs,
     &    tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &    tab,fh(i),tab,f(i),tab,att(i),
     &    tab,foi(i),tab,foii(i),tab,foiii(i)
        endif
      enddo
c
      rmserr=dsqrt(rmserr/dble(nrms))
c
   67  format(
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ':::::::::::::::::::::::::::::::::::::',/,
     & '   RMS change %: ',1pg11.4, ' Iteration: ',i2,/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ':::::::::::::::::::::::::::::::::::::')
      write (*,67) 100.d0*rmserr,itcount+1
      if (finalit.gt.0) then
      write (lupt,67) 100.d0*rmserr,itcount+1
      endif
c
        t0lim=0
c
      if (popfr(nfs,1,1).lt.0.95d0) then
        nfs=nfs+10
        nfs=min0(nfs,mxifsteps)
      else if (att(nfs).gt.0.2d0) then
        nfs=nfs+50
        nfs=min0(nfs,mxifsteps)
      else if (nte(nfs).gt.1000.d0) then
        t0lim=1
        nfs=nfs+nint(nte(nfs)*0.015d0)+100
        nfs=min0(nfs,mxifsteps)
      else if ((nte(nfs).le.10.d0).and.(att(nfs).lt.0.01d0)) then
        n100=nfs
        do i=nfs,1,-1
          if ((nte(i).le.10.d0).and.(att(i).lt.0.01d0)) n100=i
        enddo
        nfs=max0(n100,3)
      endif
c
      att(1)=1.0d0
      do i=2,nfs
        att(i)=att(i-1)*(1.d0-f(i-1))
      enddo
c
c this will all hopefully all inline with gfortran -flto
c popintfr is columns of previous steps so starts at 0
c
      call clearpop (p1)
      call copypopstep (p1, 1, popintfr)!columns up to inneredge ofzone
      do i=2,nfs
c pops at end time of previous zone, nearer source
        call copysteppop (i-1, popfr, p1)
c pops at end time of this zone, further from source
        call copysteppop (i, popfr, p2)
c average of pops from start to end in prev zone
        call averinto (0.5d0, p1, p2, p1)
        call scalepop (p1, dr)
        call copysteppop (i-1, popintfr, p2)
        call addpop (p1, p2)
        call copypopstep (p2, i, popintfr)
      enddo
c
      xhfinal=pop(1,1)
      itcount=itcount+1
c
c     poll for terminate file
c
      inquire (file='terminate',exist=iexi)
      if (iexi) goto 100
c
      if (
     &   ((rmserr.gt.rmslimit).and.(itcount.lt.mxpcits).and.(nfs.gt.2))
     &  .or.(itcount.lt.2)
     &  .or.((nfs.lt.mxifsteps).and.(t0lim.gt.0))
     &   ) goto 50
c
 100  continue
c
      te_pre=t
      de_pre=de
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      put final/inner balance into preionisation array
c
      call copysteppop (1, popfr, pop_pre)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Now get revised shock solution for changed pre shock values:
c     Shock in terms of flow velocity
c
      call shock5jump
c
c and show the result
c
      call shocksummary (6)
c      call shocksummary (luop)
      if (finalit.gt.0) then
        call shocksummary (lupt)
        open (lupc,file=fpc,status='OLD',access='APPEND')
        call shocksummary (lupc)
        close (lupc)
c
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70 format('  SHOCK 5 Final Precursor Parameter: ',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '     v_shock       : ',  1pg12.5,' km/s',/,
     & '     Q_ions        : ',  1pg12.5,' km/s',/,
     & '     Psi (Q/v)     : ',  1pg12.5,/,
     & '     Q_H           : ',  1pg12.5,' km/s',/,
     & '     Psi_H (Q_H/v) : ',  1pg12.5,/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
c redo Psi as vel0 has changes slightly
c
      psi=viofr/vel0
      write (*,70)
     &  vel0*1.d-5,viofr*1.d-5,psi,
     &  qh*1.d-5,qh/vel0
c
      if (finalit.gt.0) then
c
        write (lupt,70)
     &  vel0*1.d-5,viofr*1.d-5,psi,
     &  qh*1.d-5,qh/vel0
c
        open (lupc,file=fpc,status='OLD',access='APPEND')
        write (lupc,70)
     &  vel0*1.d-5,viofr*1.d-5,psi,
     &  qh*1.d-5,qh/vel0
        close (lupc)
c
        call copysteppop (nfs, popintfr, p1)!get inner fraction columns
        call attsig (dh0, 1.d0, 0.d0, p1, attcol, sigcol)
        do idx=1,infph
          tphot(idx)=src(idx)*attcol(idx)*wdil
        enddo
c
c        pfx=mypfx(1:nprefix)
c        pfx(nprefix+1:nprefix+5)='pcups'
c        np=nprefix+5
c        wmod='REAL'
c        caller='S5'
c        call wpsou (caller, pfx, np, wmod, te0, de0, dh0, vel0, 0.d0
c     &    , vshoc, frdr, dr, dtime, tstep, 1.d0, tphot)
c
c        call multizone (nfs, x, nte, ne, nh, popfr, popintfr, src)
        call multizone (nfs, x, nte, ne, nh, popfr, popintfr)
c
c     Final Field in tphot, away from shock front source
c
c
        call copysteppop (nfs, popintfr, p1)!get inner fraction columns
        call attsig (dh0, fi, 0.d0, p1, attcol, sigcol)
        do idx=1,infph
c  dont include wdil, added below, total src into tphot, less absorption
          tphot(idx)=src(idx)*attcol(idx)
        enddo
c
        caller='S5'
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+5)='pc_up'
        np=nprefix+5
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, vshoc,
     &    x(nfs), dr, 0.d0, 0.d0, wdil, tphot) ! note wdil
c
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, vshoc,
     &    x(nfs), dr, 0.d0, 0.d0, wdil, tphot)
c
        wmod='NFNU'
        call wpsou (caller, pfx, np, wmod, t, de, dh, 0.d0, 0.d0, vshoc,
     &    x(nfs), dr, 0.d0, 0.d0, wdil, tphot)
c
        call wrsppop (lupt)
c
        open (lupc,file=fpc,status='OLD',access='APPEND')
        spmod='REL'
        linemod='LAMB'
        call spec2 (lupc, linemod, spmod)
        close (lupc)
c
      endif
c
c  Now clear tphot so it doesn't accumulate in each iteration
c
      call zerbuf
c
c set *both* the initial population arrays and go back to computing
c post shock flow.
c
      call copypop (pop_pre, pop)
      call copypop (pop_pre, pop0)
      t=te_pre
      de=dh_pre
      dh=dh_pre
      close (lupt)
c
c make sure final jump is consistent
c
      call shock5jump
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Mulitzone emission
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine multizone (n, x, nte, ne, nh, popfr, popintfr)
c      subroutine multizone (n, x, nte, ne, nh, popfr, popintfr, src)
c
c Compute the emission from a multi-zone structure saved in the
c x, nte, ne, nh, popfr, popintfr arrays. A general low level routine.
c src not used so far, but it will in the future when we also add
c diffuse field
c
      include 'cblocks.inc'
c
      integer*4 n,idx
      real*8 x(mxifsteps),nte(mxifsteps),nh(mxifsteps),ne(mxifsteps)
      real*8 popfr(mxifsteps, mxion, mxelem)
      real*8 popintfr(mxifsteps, mxion, mxelem)
c      real*8 src(mxinfph)
c
      real*8 t,tdw,tup,de,dh,dvdw,dvup,fi
      real*8 gdil,x0,dx,dxdw,dxup
      real*8 p1(mxion, mxelem), p2(mxion, mxelem)
c
      character imod*4,specmode*4,subname*32
c
      real*8 feldens
c
      gdil=0.5d0
      fi=1.d0
c
      call zerbuf ()
c global modes
      jspot='NO'
      jcon='YES'
      jgeo='P'
      jtrans='DWUP'
      jden='C'
c local only to mulitzone
      specmode='UP'
      subname='MultiZone'
c
      do idx=1,n-1
c
        dx=x(idx+1)-x(idx)
        x0=x(idx)
c
        call copysteppop (idx, popintfr, popint)
c
        call copysteppop (idx, popfr, p1)
        call copysteppop (idx+1, popfr, p2)
        call averinto (0.5d0, p1, p2, pop)!average ions between points
        t=0.5d0*(nte(idx)+nte(idx+1))
        dh=nh(idx)
        de=feldens(dh,pop)
        ne(idx)=de
c
        call localem (t, de, dh)
        call totphot2 (t, dh, fi, x0, dx, 0.d0, gdil, specmode)
        call zetaeff (fi, dh)
c
c     accumulate spectrum
c
        tdw=nte(idx+1)
        tup=nte(idx)
        dxdw=dx
        dvdw=0.d0
        dxup=x0
        dvup=0.d0
c
        call newdif2 (tdw, tup, dh, fi, x0, dxdw, dvdw, dxup, dvup,
     &   gdil, jtrans)
c
        imod='ALL'
        call sumdata (t, de, dh, fi, dx, dx, x0, imod)
      enddo
c
      call avrdata
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compsh5 (iteration, maxits)
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iteration, maxits
c
      integer*4 i,j,idx
      logical iexi
      real*8 n,invn,f,nexp
      real*8 b0,b1,b2,b3,b4,blum,binlum,pe
      real*8 hdt,tscale,dt0
c
c structural markers
c
      real*8 timMark(7),disMark(7),bmMark(7)
      real*8 prMark(7),rhMark(7),dhMark(7),deMark(7)
      real*8 templim,ft,cft
      integer*4 tempidx
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           Functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 feldens,frho,fmua,fdynamictimestep
c
      character tab*1
      tab=','
c
        do idx=1,7
          timMark(idx)=0.0d0
          disMark(idx)=0.0d0
          bmMark(idx)=0.0d0
          prMark(idx)=0.0d0
          rhMark(idx)=0.0d0
          dhMark(idx)=0.0d0
          deMark(idx)=0.0d0
        enddo
c
c global modes
c
      jspot='NO'
      jcon='YES'
      jgeo='P'
      jtrans='LODW'
      jden='B'
c local to S5
      specmode='DW'
      subname='Shock 5'
      write(*,*) iteration,maxits
c
   10 format(
     & ' *********************************************************',/,
     & '  SHOCK 5     Shock Iteration: ',i2.2,' of ',i2.2,/,
     & ' *********************************************************')
      write (*,10) iteration,maxits
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,10) iteration,maxits
        call shocksummary (luop)
        close (luop)
        open (lusp,file=fsp,status='OLD',access='APPEND')
        write (lusp,10) iteration,maxits
        call shocksummary(lusp)
        close (lusp)
        if (bandsmod.eq.'Y') then
        open (lupb,file=fpb,status='OLD',access='APPEND')
        call shocksummary (lupb)
        close (lupb)
        endif
        if (ratmod.eq.'Y') then
        open (lurt,file=fr,status='OLD',access='APPEND')
        call shocksummary (lurt)
        close (lurt)
        endif
        if (dynmod.eq.'Y') then
        open (ludy,file=fd,status='OLD',access='APPEND')
        call shocksummary (ludy)
        close (ludy)
        endif
        if (allmod.eq.'Y') then
        open (lual,file=fa,status='OLD',access='APPEND')
        call shocksummary (lual)
        close (lual)
        endif
        if (fclmod.eq.'Y') then
        open (lucl,file=fcl,status='OLD',access='APPEND')
        call shocksummary (lucl)
        close (lucl)
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Begin Main Calculation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 format(//'  T0',1pg14.7,' V0',1pg14.7,/,
     &         ' RH0',1pg14.7,' P0',1pg14.7,' B0',1pg14.7,//,
     &         '  T1',1pg14.7,' V1',1pg14.7,/,
     &         ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7)
c
      write (*,20) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,rho1,
     &pr1,bm1*1.d6
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,20) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &   rho1,pr1,bm1*1.d6
        close (luop)
        open (lusp,file=fsp,status='OLD',access='APPEND')
        write (lusp,20) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &   rho1,pr1,bm1*1.d6
        close (lusp)
      endif
c
      wdilt0=te1
      t=te1
      dh=dh1
      de=de1
      dr=0.d0
      ve=vel1
      vpo=vel1
      bmag=bm1
      bm0=bm1
      dv=ve*0.01d0
      fi=1.d0
      rad=1.d38
      if (wdil.eq.0.5d0) rad=0.d0
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
      call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
      endif
c
      call cool (t, de, dh)
      netloss=(eloss-egain)
c
      step=1
c
   30 format(//,a4,a1,14(3x,a6,4x,a1))
   40 format(i4.3,a1,14(1pg13.6,a1))
c
   50 format(/,
     &     ' #    Te(K)       ne(/cm3)    nH(/cm3)    ni(/cm3)    B',
     &'(G)        FHI         Time(s)     dt(s)       Dist',
     &   '(cm)    dr(cm)      v(cm/s)     Tloss(erg/cm3/s)')
      if (vmod.eq.'MINI') then
        write (*,50)
      endif
   60 format (i4,6(1x,1pg11.4),6(1x,1pg11.4))
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,30) 'Step',tab,'Dist.',tab,'Time',tab,'Te',tab,'de',
     &   tab,'dh',tab,'en',tab,'mu',tab,'Lambda',tab,'Tloss',tab,'dlos',
     &   tab,'FHI',tab,'FHII',tab,'Bbar',tab,'dr',tab
        close(luop)
      endif
c
      en=zen*dh_neu
      cspd=dsqrt(gammaeos*pr_neu/rh_neu)
      wmol=rh_neu/(en+de_neu)
      mu=fmua(de_neu,dh_neu)
c
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,40) -1,tab,0.d0,tab,0.d0,tab,te_neu,tab,de_neu,tab,
     &   dh_neu,tab,en,tab,mu,tab,0.0d0,tab,0.0d0,tab,0.d0,tab,
     &   pop_neu(1,1),tab,pop_neu(2,1),tab,bm_neu,tab,0.d0,tab
        close(luop)
      endif
      if (vmod.eq.'MINI') then
        write (*,60) -1,te_neu,de_neu,dh_neu,en,bm_neu
     &     ,pop_neu(1,1),0.d0,0.d0,0.0d0,0.d0,vel0,0.0d0
      endif
c
      en=zen*dh_pre
      cspd=dsqrt(gammaeos*pr_pre/rh_pre)
      wmol=rh_pre/(en+de_pre)
      mu=fmua(de_pre,dh_pre)
c
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,40) 0,tab,0.d0,tab,0.d0,tab,te_pre,tab,de_pre,tab,
     &   dh_pre,tab,en,tab,mu,tab,0.0d0,tab,0.0d0,tab,0.d0,tab,
     &   pop_pre(1,1),tab,pop_pre(2,1),tab,bm_pre,tab,0.d0,tab
        close(luop)
      endif
c
      if (vmod.eq.'MINI') then
        write (*,60) 0,te_pre,de_pre,dh_pre,en,bm_pre,
     &   pop(1,1),0.d0,0.d0,0.0d0,0.d0,vel0,0.0d0
      endif
c
      en=zen*dh1
      ue=(1.d0/(gammaeos-1.d0))*(en+de)*rkb*te1
      tnloss=tloss/(en*de1)
c
      cspd=dsqrt(gammaeos*pr1/rho1)
      wmol=rho1/(en+de1)
      mu=fmua(de1,dh1)
c
C     if (finalit.gt.0) then
C       write (luop,40) step,tab,dist(step),tab,timlps(step),tab,t,tab,
C    &   de1,tab,dh1,tab,en,tab,mu,tab,tloss/(dh*dh),tab,tloss,tab,dlos,
C    &   tab,pop(1,1),tab,pop(2,1),tab,bm0,tab,dr,tab
c
C       close (luop)
C     endif
c
C     if (vmod.eq.'MINI') then
C       write (*,50)
C       write (*,60) step,t,de,dh,en,bm0,pop(1,1),timlps(step),dt,
C    &   dist(step),dr,veloc(step),tloss
C     endif
c
c
      tscale=0.001d0
      dt=fdynamictimestep(t,dh,dist(step+1),ve,pop,netloss)
      dt=tscale*dt
      hdt=0.5d0*dt
c
c     if ((tmod.eq.'M').and.(dt.gt.utime)) dt=utime
c     if ((tmod.eq.'N').and.(dt.lt.utime)) dt=utime
c     if (tmod.eq.'S') dt=utime
c
      rhotot=frho(de,dh)
      wmol=rhotot/(en+de)
      mu=fmua(de,dh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Setup Initial first step in post shock-front gas(#1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      rmod='NEAR'
c      call isochorflow (t, de, dh, ve, bmag, rmod, netloss, dt)
c      call isobarflow (t, de, dh, ve, bmag, rmod, netloss, dt)
      call rankhug (t, de, dh, ve, bmag, netloss, dt)
c
      dr=dt*(vel0+vel1)*0.5d0
      dv=vel1-vel0
c
c      write(*,200) tloss,dr,dv
c
c     init arrays first step
c
      xh(1)=pop(2,1)
      xh(2)=pop(2,1)
      veloc(1)=vel0
      veloc(2)=(vel0+vel1)*0.5d0
c
      te(1)=te0
      te(2)=t
      deel(1)=de0
      deel(2)=de
      dhy(1)=dh0
      dhy(2)=dh
      bmg(1)=bm0
      bmg(2)=(bm0+bm1)*0.5d0
c
      dist(1)=0.0d0
      dist(2)=dr
      timlps(1)=0.0d0
      timlps(2)=dt
c
      step=1
c
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,40) step,tab,dist(step),tab,
     &   timlps(step),tab,te1,tab,
     &   de1,tab,dh1,tab,en,tab,mu,tab,tloss/(dh*dh),
     &   tab,tloss,tab,dlos,
     &   tab,pop(1,1),tab,pop(2,1),tab,bm1,tab,dr,tab
        close (luop)
      endif
c
      if (vmod.eq.'MINI') then
        write (*,50)
        write (*,60) step,te1,de1,dh1,dh1*zen,bm1,pop(1,1),
     &   timlps(step),0.d0,
     &   dist(step),0.d0,veloc(step),tloss
      endif
c
      count=0
c
      call copypop (pop0, pop)
      if (finalit.gt.0) then
        caller='S5'
        pfx=mypfx(1:nprefix)
        np=nprefix
        call wbal (caller, pfx, np, pop)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SHOCK Step/Iteration reentry point:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70 count=0
      rad=dist(step+1)
c
c save init ionisation and fields and get intial cooling/heating rates
c
      call copypop (pop0, pop)
      de0=feldens(dh0,pop0)
      if (photonmode.ne.0) then
        call localem (te0, de0, dh0)!callscool
        call totphot2 (te0, dh0, fi, rad, 0.0d0, 0.d0, wdil, specmode)
        call zetaeff (fi, dh0)
        call cool (te0, de0, dh0)!includesheatingifany
      else
        call cool (te0, de0, dh0)!includesheatingifany
      endif
c
c get intial time step dt, and hdt = 0.5d0*dt, init t = te0
c
      tscale=dmin1(1.d0,0.01d0*dble(step*step)) ! gentle leadin
      dt0=fdynamictimestep(te0,dh0,rad,vel0,pop,netloss)
      dt=tscale*dt0
      hdt=0.5d0*dt
c
c      evolve ionisation over hdt to get pop at midpoint
c
   80 call copypop (pop0, pop)
      call timion (te0, de0, dh0, xhii, hdt)
c
c     get new approx netloss rate at midpoint
c
      if (photonmode.ne.0) then
        call localem (te0, de0, dh0)!callscool
        call totphot2 (te0, dh0, fi, rad, hdt*vel0, 0.5d0*dv, wdil,
     &   specmode)
        call zetaeff (fi, dh0)
        call cool (te0, de0, dh0)!includesheatingifany
      else
        call cool (te0, de0, dh0)!includesheatingifany
      endif
c
c common interation in partially ionised x-ray heated tails
c
      if (dabs(netloss).gt.0.d0) then
        f=dabs((eloss-egain)/netloss)
        if ((f.ge.10.0d0).and.(count.lt.6).and.(dabs(dlos).ge.1.d-3))
     &   then
          tscale=tscale*0.333d0
          dt=tscale*dt0
          hdt=0.5d0*dt
c        write(*,*) 'Pre Reduce step:',f,dt,tscale
          count=count+1
          goto 80
        endif
      endif
      count=0
c
      netloss=(eloss-egain)
c
c loss dampening if close to heating/cooling balance
c prevents wild instabilty in x-ray neq tails
c
      nexp=4.d0
      f=(dlos**nexp)/((5.0d-4)**nexp+dlos**nexp)
c      write(*,*) dlos,f
c
c go back to initial ionisation state and go forward on midpoint
c loss rate, reduced if close to heating balance
c
      call copypop (pop0, pop)
      de0=feldens(dh0,pop0)
c
      rmod='NEAR'
c      call isochorflow (te0, de0, dh0, vel0, bmag, rmod, netloss*f, dt)
c      call isobarflow (te0, de0, dh0, vel0, bmag, rmod, netloss*f, dt)
      call rankhug (te0, de0, dh0, vel0, bmag, netloss*f, dt)
      dr=dt*(vel0+vel1)*0.5d0
      dv=vel1-vel0
      f=dabs(0.5d0*((te1-te0)/(te1+te0)))
c
c hopefully rare fail:  print warning
c
      if (f.gt.0.05d0) then
        tscale=tscale*0.333d0
        dt=tscale*dt0
        hdt=0.5d0*dt
        write (*,*) 'WARNING ! Post Te Reduce step:',dt,tscale
        count=0
        goto 80
      endif
c
c      write(*,*) 'Full Step:',te0,te1,vel0,vel1
c
c mean values T,ne,nh, and move ionisation forward the whole step
c
      t=(te0+te1)*0.5d0
      dh=(dh0+dh1)*0.5d0
      de=feldens(dh,pop)
c
c calculate the radiation field and atomic rates over whole step
c
      if (photonmode.ne.0) then
        call localem (t, de, dh)!callscool
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
        call zetaeff (fi, dh)
        call cool (t, de, dh)
      else
        call cool (t, de, dh)
      endif
c
c advance ionisation a whole timestep at mean t
c
      call timion (t, de, dh, xhii, dt)
c
c save end ionisation
c
      call copypop (pop, pop1)
c
c mean losses/gains at mean ionisation state/temp
c
      call averinto (0.5d0, pop0, pop1, pop)
      t=(te0+te1)*0.5d0
      dh=(dh0+dh1)*0.5d0
      de=feldens(dh,pop)
c
      if (photonmode.ne.0) then
        call localem (t, de, dh)!callscool
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
        call zetaeff (fi, dh)
        call cool (t, de, dh)
      else
        call cool (t, de, dh)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c End step iterations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      step=step+1
c
c     record step
c
      te(step)=t
      dhy(step)=dh
      deel(step)=de
      bmg(step)=(bm0+bm1)*0.5d0
c
      xh(step)=pop(2,1)
      veloc(step)=(vel1+vel0)*0.5d0
c
      dist(step)=dist(step-1)+dr!fromnowon,distisendofstep
      timlps(step)=timlps(step-1)+dt
c
   90 format(i4,1x,i4)
  100 format(1pg14.7,'K',1x,1pg14.7,'cm/s',1x,1pg14.7,'g/cm3',1x,
     &1pg14.7,'dyne/cm2',1x,1pg14.7,'Gauss',1x,1pg14.7,'ergs/cm3/s')
  110 format(2(1pg14.7,' cm',1x),2(1pg14.7,' s',1x))
  120 format(1pg14.7,' ergs/cm^3',1x,1pg14.7,' /cm^3',1x,
     &1pg14.7,' cm/s',1x,1pg14.7,' g/particle ')
      en=zen*dh
c
      ue=(1.d0/(gammaeos-1.d0))*(en+de)*rkb*t
      tnloss=tloss/(en*de)
c
      cspd=dsqrt(gammaeos*(pr1+pr0)/(rho1+rho0))
      wmol=(rho1+rho0)/(2.d0*(en+de))
      mu=fmua(de,dh)
c
      mb=(bm0+bm1)*0.5
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
c structural markers
c
      templim=1.0d7
  125 continue
      tempidx=idnint(dlog10(templim))
      if ((te1.le.templim).and.(te0.gt.templim)) then
          ft=(templim-te1)/(te0-te1)
          cft=1.d0-ft
          timMark(tempidx)=cft*timlps(step-1)+ft*timlps(step-1)
          disMark(tempidx)=cft*dist(step-1)+ft*dist(step-1)
          bmMark(tempidx)=cft*bm0+ft*bm1
          prMark(tempidx)=cft*pr0+ft*pr1
          rhMark(tempidx)=cft*rho0+ft*rho1
          dhMark(tempidx)=cft*dh0+ft*dh1
          deMark(tempidx)=cft*de0+ft*de1
      endif
      templim=templim*0.1d0
      if (templim.gt.10.d0) goto 125
cc
      if (vmod.eq.'MINI') then
        if (mod(step,20).eq.1) write (*,50)
        write (*,60) step,t,de,dh,en,mb,pop(1,1),timlps(step),dt,
     &   dist(step),dr,veloc(step),tloss
      endif
c
      if (ratmod.eq.'Y') then
        open (lurt,file=fr,status='OLD',access='APPEND')
        wmod='LOSS'
        call wmodel (lurt, t, de, dh, dr, wmod)
        close (lurt)
      endif
c
      if (dynmod.eq.'Y') then
        open (ludy,file=fd,status='OLD',access='APPEND')
        write (ludy,90) step,count
        write (ludy,110) dist(step),dr,timlps(step),dt
        write (ludy,100) te0,vel0,rho0,pr0,bm0,l0
        write (ludy,100) te1,vel1,rho1,pr1,bm1,l11
        write (ludy,120) ue,en,cspd,wmol
        close (ludy)
      endif
c
      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
        write (luop,40) step,tab,dist(step),tab,timlps(step),tab,t,tab,
     &   de,tab,dh,tab,en,tab,mu,tab,tloss/(dh*dh),tab,tloss,tab,dlos,
     &   tab,pop(1,1),tab,pop(2,1),tab,(0.5*(bm0+bm1)),tab,dr,tab
        close (luop)
      endif
c
  130 format(43(1pg12.5,a2))
      if (fclmod.eq.'Y') then
        do i=1,atypes
          meanion(i)=0.d0
          ionbar=0.d0
          do j=1,maxion(i)
            ionbar=ionbar+(pop(j,i)*j)
          enddo
          meanion(i)=dmax1(ionbar-1.0d0,0.0d0)
        enddo
c
        open (lucl,file=fcl,status='OLD',access='APPEND')
        en=zen*dh
        n=(de*dh)
        if (jnorm.eq.0) n=(de*dh)
        if (jnorm.eq.1) n=(dh*dh)
        if (jnorm.eq.2) n=(de*en)
        if (jnorm.eq.3) n=((de+en)*(de+en))
        if (jnorm.eq.4) n=(de*de)
        invn=1.d0/n
        write (lucl,130) t,tab,de,tab,dh,tab,en,tab,(0.5*(rho0+rho1)),
     &   tab,pop(1,1),tab,pop(2,1),tab,mu,tab,tloss,tab,tloss*invn,(tab,
     &   coolz(j)*invn,j=1,atypes),tab,eloss*invn,tab,egain*invn,tab,
     &   netloss*invn
        close (lucl)
c
      endif
c
      if (allmod.eq.'Y') then
        open (lual,file=fa,status='OLD',access='APPEND')
        write (lual,140)
        write (lual,150) step,t,de,dh,vel1,rho1,pr1,dist(step),
     &   timlps(step)
  140 format(//,' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  V1(cm/s)    ',
     &                ' Rho1(g/cm^3) ',
     &                ' Pr1(erg/cm^3)',
     &                '   Dist.(cm)  ',
     &                ' Elps. Time(s)')
  150 format(1x,i4,8(1pg14.7)/)
c
        call wionabal (lual, pop)
        call wionabal (lual, popint)
        close (lual)
      endif
c
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          open (luions(i),file=filn(i),status='OLD',access='APPEND')
c
  160 format(1x,i4,3(', ',1pg12.5),31(', ',1pg12.5))
          if (i.eq.1) then
            write (luions(i),160) step,dist(step)-(dr*0.5),dr,t,(pop(j,
     &       iel1),j=1,maxion(iel1))
          endif
          if (i.eq.2) then
            write (luions(i),160) step,dist(step)-(dr*0.5),dr,t,(pop(j,
     &       iel2),j=1,maxion(iel2))
          endif
          if (i.eq.3) then
            write (luions(i),160) step,dist(step)-(dr*0.5),dr,t,(pop(j,
     &       iel3),j=1,maxion(iel3))
          endif
          if (i.eq.4) then
            write (luions(i),160) step,dist(step)-(dr*0.5),dr,t,(pop(j,
     &       iel4),j=1,maxion(iel4))
          endif
c
          close (luions(i))
c
        enddo
c
c     end .tsr files
c
      endif
c
c
c     poll for photons file
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
c
      if ((lmod.eq.'Y').or.(iexi)) then
c
c     Local photon field
c
        specmode='LOCL'
        call localem (t, de, dh)
        rad=dist(step)
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
        caller='S5'
        pfx=mypfx(1:nprefix)
        pfx(nprefix+1:nprefix+5)='plocl'
        np=nprefix+5
        wmod='NORM'
        dva=vel1-vel0
        en=zen*dh
        n=(de*dh)
        if (jnorm.eq.0) n=(de*dh)
        if (jnorm.eq.1) n=(dh*dh)
        if (jnorm.eq.2) n=(de*en)
        if (jnorm.eq.3) n=((de+en)*(de+en))
        if (jnorm.eq.4) n=(de*de)
        invn=1.d0/n
        call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva, vshoc,
     &   dist(step), dr, timlps(step), dt, invn, tphot)
c
c       pfx=mypfx(1:nprefix)
c       pfx(nprefix+1:nprefix+5)='_slec'
c       np=nprefix+5
c       sfx='sh5'
c       call newfile (pfx, np, sfx, 3, fn)
c       fpf=fn(1:np+8)
c       write (*,*) 'fpf''',fa,''''
c       open (lupf,file=fpf,status='NEW')
c       spmod='ABS'
c       linemod='LAMB'
c       call speclocal (lupf, tloss, eloss, egain, dlos, t, de, dh,
c    &   pop(1,1), dist(step+1), dr, linemod, spmod)
c       close (lupf)
c
c     reset mode and tphot
c
        specmode='DW'
        rad=dist(step+1)
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
      endif
c
      if (bandsmod.eq.'Y') then
        open (lupb,file=fpb,status='OLD',access='APPEND')
c
c
        b0=0.d0
        b1=0.d0
        b2=0.d0
        b3=0.d0
        b4=0.d0
        blum=0.d0
        do idx=1,infph-1
          wid=photev(idx+1)-photev(idx)
          pe=photev(idx)
          binlum=tphot(idx)*wid*evplk*fpi/dr
          blum=blum+binlum
          if ((pe.gt.0.0d0).and.(pe.le.100.0d0)) b0=b0+binlum
          if ((pe.gt.100.0d0).and.(pe.le.500.0d0)) b1=b1+binlum
          if ((pe.gt.500.0d0).and.(pe.le.1000.0d0)) b2=b2+binlum
          if ((pe.gt.1000.0d0).and.(pe.le.2000.0d0)) b3=b3+binlum
          if ((pe.gt.2000.0d0).and.(pe.le.10000.0d0)) b4=b4+binlum
        enddo
c
        b0=b0/tloss
        b1=b1/tloss
        b2=b2/tloss
        b3=b3/tloss
        b4=b4/tloss
        blum=blum/tloss
c
        if (b0.lt.epsilon) b0=0.d0
        if (b1.lt.epsilon) b1=0.d0
        if (b2.lt.epsilon) b2=0.d0
        if (b3.lt.epsilon) b3=0.d0
        if (b4.lt.epsilon) b4=0.d0
        if (blum.lt.epsilon) blum=0.d0
c
  170 format(17(1pg12.5,a1))
        write (lupb,170) t,tab,de,tab,dh,tab,en,tab,pop(1,1),tab,pop(2,
     &   1),tab,mu,tab,tloss,tab,tloss*invn,tab,fflos/tloss,tab,b0,tab,
     &   b1,tab,b2,tab,b3,tab,b4
c
        close (lupb)
      endif
c
c     get mean ionisation state for step
c
      rdis=dist(step)-dr*0.5d0!middleofstep,dist(step)isnowend
c      cmpf=rho1/rho0
      fi=1.0d0
c
c     accumulate spectrum
c
      tdw=te0
      tup=te1
      drdw=dr
      dvdw=vel0-vel1
      drup=dist(step)
      dvup=dsqrt(vpo*vel0)
      frdw=0.5d0
c
      call localem (t, de, dh)
      if (photonmode.ne.0) then
        call zetaeff (fi, dh)
        call newdif2 (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup,
     &   frdw, jtrans)
      endif
c
      imod='ALL'
      call sumdata (t, de, dh, fi, dr, dr, rdis, imod)
c
c     record line ratios
c
      if (ox3.ne.0) then
        hoiii(step)=(fluxm(6,ox3)+fluxm(8,ox3))/(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then
        hoii(step)=(fluxm(1,ox2)+fluxm(2,ox2))/(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then
        hnii(step)=(fluxm(7,ni2)+fluxm(10,ni2))/(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then
        hsii(step)=(fluxm(1,su2)+fluxm(1,su2))/(fluxh(2)+epsilon)
      endif
      if (ox1.ne.0) then
        if (ispo.eq.' OI') hsii(step)=fluxm(3,ox1)/(fluxh(1)+epsilon)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set up initial conditions for the next zone
c     (= end conditions in previous zone)
c     Calculate new cooling times and next step
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop1, pop0)
      call copypop (pop1, pop)
c
      t=te1
      te0=te1
      rho0=rho1
      pr0=pr1
      vel0=vel1
      dh0=dh1
      de0=de1
      bm0=bm1
      bmag=bm0
      count=0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Program endings
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ionisation, finish when mean ionisation < 1%
c
c      if (jend.eq.'A') then
c
        totn=0.d0
c
        totn=zen
c
        ionbar=0.d0
c
        do i=1,atypes
c
          do j=2,maxion(i)
            ionbar=ionbar+(pop(j,i)*zion(i))
          enddo
c
        enddo
c
        ionbar=ionbar/totn
c
        if ((jend.eq.'A').and.(ionbar.lt.0.01d0)) goto 180
c
c      endif
c
c     Specific species ionisation limit
c
c
      if (jend.eq.'B') then
        if (pop(jpoen,ielen).lt.fren) goto 180
      endif
c
c stop first iteration at heating==cooling if maxits>1
c
      if ((iteration.le.1)
     &  .and.(maxits.gt.1)
     &  .and.(dlos.lt.0.5d0)) goto 180
c
c     Temperature limit
c
      if (jend.eq.'C') then
        if (t.lt.tend) goto 180
      endif
c
c     Temperature limit, plus 95% neutral
c
      if ((jend.eq.'S').and.(t.lt.tend)) then
        if (ionbar.lt.0.05d0) goto 180
      endif
c
c     Distance Limit
c
      if ((jend.eq.'D').and.(dist(step+1).ge.diend)) goto 180
c
c
c     Time Limit
c
      if ((jend.eq.'E').and.(timlps(step).ge.timend)) goto 180
c
c
c     thermal balance dlos<1e-2
c
      if ((jend.eq.'F').and.(dlos.lt.1.d-2)) goto 180
c
c
c     cooling function test, finish when tloss goes -ve
c
      if ((jend.eq.'G').and.(tloss.lt.0.d0)) goto 180
c
c     poll for terminate file
c
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 180
c
c     otherwise go to normal iteration loop, if interlocks
c     permit
c
      if ((t.gt.100.d0).and.(step.lt.(mxnsteps-1))) then
        goto 70
      endif
c
  180 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     End model
c     write out spectrum etc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Final Downstream Field
c
      if (finalit.gt.0) then
        if (jspec.eq.'YES') then
c
          caller='S5'
          pfx=mypfx(1:nprefix)
          pfx(nprefix+1:nprefix+7)='spec_dw'
          np=nprefix+7
          dva=0
c
          wmod='LFLM'
          call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva,
     &     vshoc, dist(step), dr, timlps(step), dt, 0.5d0, tphot)
c
          pfx=mypfx(1:nprefix)
          pfx(nprefix+1:nprefix+4)='spdw'
          np=nprefix+4
          wmod='REAL'
          call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva,
     &     vshoc, dist(step), dr, timlps(step), dt, 0.5d0, tphot)
c
          wmod='NFNU'
          call wpsou (caller, pfx, np, wmod, t, de, dh, vel1, dva,
     &     vshoc, dist(step), dr, timlps(step), dt, 0.5d0, tphot)
c
        endif
c
c     Upstream photon field
c
        specmode='UP'
        rad=dist(step+1)
        call totphot2 (t, dh, fi, rad, dr, dv, wdil, specmode)
c
        caller='S5'
        pfx=mypfx(1:nprefix)
        np=nprefix
        dva=vel1-vel0
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &   vshoc, dist(step), dr, timlps(step), dt, wdil, tphot)
c
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &   vshoc, dist(step), dr, timlps(step), dt, wdil, tphot)
c
        if (lmod.eq.'Y') then
c
          wmod='NFNU'
          call wpsou (caller, pfx, np, wmod, te1, de1, dh1, vel1, dva,
     &     vshoc, dist(step), dr, timlps(step), dt, wdil, tphot)
c
        endif
      endif
c
c     dynamics
c
      call avrdata

      if (finalit.gt.0) then
        open (luop,file=fl,status='OLD',access='APPEND')
c
  190 format(//,' Model ended: , Distance:, ',1pg14.7,
     &', Time:, ',1pg14.7,', Temp:, ',1pg12.5,/)
        write (luop,190) dist(step),timlps(step),t
        write (*,190) dist(step),timlps(step),t
c
c structural markers
c
  191 format(' Marker',i1,':, ',1pe12.5,' (K)',/,
     &       '      t',i1,':, ',1pg12.5,/,
     &       '      d',i1,':, ',1pg12.5,/,
     &       '      B',i1,':, ',1pg12.5,/,
     &       '      P',i1,':, ',1pg12.5,/,
     &       '    rho',i1,':, ',1pg12.5,/,
     &       '     nH',i1,':, ',1pg12.5,/,
     &       '     ne',i1,':, ',1pg12.5,/)
        do idx=7,1,-1
          if ( timMark(idx).gt.0.d0) then
          templim=10.0d0**dble(idx)
          write(luop,191) idx,templim
     &              ,idx,timMark(idx)
     &              ,idx,disMark(idx)
     &              ,idx,bmMark(idx)
     &              ,idx,prMark(idx)
     &              ,idx,rhMark(idx)
     &              ,idx,dhMark(idx)
     &              ,idx,deMark(idx)
          write(*,191) idx,templim
     &              ,idx,timMark(idx)
     &              ,idx,disMark(idx)
     &              ,idx,bmMark(idx)
     &              ,idx,prMark(idx)
     &              ,idx,rhMark(idx)
     &              ,idx,dhMark(idx)
     &              ,idx,deMark(idx)
          endif
        enddo
c
        call wrsppop (luop)
        close (luop)
c
        open (lusp,file=fsp,status='OLD',access='APPEND')
        spmod='REL'
        linemod='LAMB'
        call spec2 (lusp, linemod, spmod)
        close (lusp)

      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fdynamictimestep(t,dh,x,vel,p,netloss)
      include 'cblocks.inc'
c
      real*8 t,dh,de,vel,netloss
      real*8 p(mxion, mxelem)
      real*8 en, pr, ue, fi, x
      real*8 cltime, absdr, abstime
      real*8 ctime, rtime, dt
      real*8 feldens, fcolltim, frectim2
      character specmode*4
c
      de=feldens(dh,p)
      en=zen*dh+de
      pr=en*rkb*t
      ue=(1.d0/(gammaeos-1.d0))*pr
c
c initial cooling rate for initial timestep guess
c
      if (photonmode.ne.0) then
        specmode='DW'
        fi=1.d0
        call localem (t, de, dh)!callscool
        call totphot2 (t, dh, fi, x, 0.d0, 0.0d0, wdil, specmode)
        call zetaeff (fi, dh)
c
c call cool again to get egain rates
c
        call cool (t, de, dh)
      else
        call cool (t, de, dh)
      endif
      netloss=(eloss-egain)
c
      cltime=(ue/(dabs(netloss)+epsilon))
      absdr=1.0d38
      abstime=1.d0/epsilon
      if (photonmode.ne.0) then
        absdr=1.0d18/dh0
        call absdis2 (dh, 1.0d0, 0.5d0, absdr, 0.0d0, p)
        abstime=absdr/(dabs(vel)+epsilon)
      endif
      ctime=fcolltim(de)
      rtime=frectim2(de0)
c
c weighted harmonic sum, coeffs determined empirically
c
      dt=1.d0/((0.5d0/abstime)+(1.d0/(cltime*0.025d0))+(1.d0/(rtime*
     &0.05d0))+(1.d0/(ctime*0.01d0)))
c
c      write (*,100) t,ue,netloss,cltime,ctime,rtime,abstime,dt
c
      fdynamictimestep=dt
c
      return
      end
