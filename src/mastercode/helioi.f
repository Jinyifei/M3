cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DERIVES INTENSITY OF IMPORTANT HEI LINES AND THE RESULTANT
c     USES FUNCTION FGAUNT
c     COOLING RATE : HELOS
c
c
      subroutine helioi (t, de, dh, helos)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh, helos, u
      real*8 abup,ara,cop1,cos1,f,frp1,frs1
      real*8 omep1,omes1,ra21,ra3,rap1,rap3,ras1
      real*8 reche,tns3,xp1,xs1
      real*8 a2s31s1
      real*8 a2p31s1
      real*8 at1s
      real*8 a2s11s1
c
C     real*8 ome1s12p3
C     real*8 ome1s12s3
C     real*8 ome1s12s1
      integer*4 atom
c new code:
      integer*4 idx,upid,n,i
      real*8 t4,lnte,invt,abde,emiss
      real*8 a,b,c,d
      real*8 ai,bi,ci
      real*8 cr(20),sum,coolf
c DATA
c new BT splines to replace polynomial fit, but these
c reproduce the poly fits as near as possible.
c TODO. breaks rule of no code data, will move to
c data file soon - RSS
c
      real*8 beta
      real*8 xs(21)
      real*8 qss_y(21)
      real*8 qss_y2(21)
      real*8 qsp_y(21)
      real*8 qsp_y2(21)
      real*8 qsp3_y(21)
      real*8 qsp3_y2(21)
      real*8 s_x(11)
C     real*8 s12_y(11)
C     real*8 s12_y2(11)
      real*8 s13_y(11)
      real*8 s13_y2(11)
C     real*8 s14_y(11)
C     real*8 s14_y2(11)
C     real*8 s15_y(11)
C     real*8 s15_y2(11)
C     real*8 s16_y(11)
C     real*8 s16_y2(11)
      real*8 s17_y(11)
      real*8 s17_y2(11)
c
c
c           Functions
c
      real*8 fsplint
      real*8 qss,qsp,qsp3,frs3,fusp
c
c old poly fit
c      pol5(x)=dmax1(0.0d0,(((((((((-(6.5289d0*x))+41.54554d0)*x)-
c     &97.135778d0)*x)+97.0517d0)*x)-32.02831d0)*x)-0.045645d0)
c
      qss(u)=10.d0**(fsplint(xs,qss_y,qss_y2,21,    (u/(u+ 1.0d0)) ))
      qsp(u)=10.d0**(fsplint(xs,qsp_y,qsp_y2,21,    (u/(u+ 5.0d0)) ))
      qsp3(u)=10.d0**(fsplint(xs,qsp3_y,qsp3_y2,21, (u/(u+10.0d0)) ))
c old fns
c      qss(u)=1.d-8*pol5(u**0.3333333d0)
c      qsp(u)=5.73d-9*pol5((0.175d0+(0.2d0*u))**0.3333333d0)
c      qsp3(u)=3.70d-7*pol5((0.218d0+(0.1d0*u))**0.3333333d0)
c
      frs3(u)=0.775d0*(u**0.0213d0)
      fusp(u)=0.42d0*(u**0.23d0)

      data xs/0.0d0,0.05d0,0.10d0,0.15d0,0.20d0,0.25d0,0.30d0,
     & 0.35d0,0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,
     & 0.75d0,0.80d0,0.85d0,0.90d0,0.95d0,1.00d0/
c QSS  TC 1e4
      data qss_y/-1.3836238000e+01,-1.2632563296e+01,-1.1428888592e+01,
     &-1.0225213887e+01,-9.0494866438e+00,-8.2444050668e+00,
     &-7.9280646843e+00,-7.7692104500e+00,-7.6662966339e+00,
     &-7.5947805889e+00,-7.5438466106e+00,-7.5084542528e+00,
     &-7.4861301161e+00,-7.4758225038e+00,-7.4770530403e+00,
     &-7.4887852989e+00,-7.5065419490e+00,-7.5242392468e+00,
     &-7.5419365446e+00,-7.5596338424e+00,-7.5773311402e+00/
      data qss_y2/-1.0410826790e+00, 2.0821653580e+00,-7.2875787531e+00,
     & 2.7068152054e+01,-1.6805893778e+02,-2.4438199980e+02,
     &-2.7391929832e+01,-2.4017036557e+01,-1.0796927622e+01,
     &-8.1499035957e+00,-6.0004180751e+00,-5.1483133039e+00,
     &-4.7700593494e+00,-4.6111078585e+00,-4.4770663367e+00,
     &-2.6847598348e+00, 7.5756607572e-01,-2.0305894813e-01,
     & 5.4669716806e-02,-1.5619919087e-02, 7.8099595437e-03/
c
c QSP  TC 5e4
      data qsp_y/-1.1418618188e+01,-1.0165844261e+01,-8.9205948696e+00,
     &-8.3915228136e+00,-8.1778269085e+00,-8.0451923253e+00,
     &-7.9518068838e+00,-7.8819630966e+00,-7.8284067675e+00,
     &-7.7871681140e+00,-7.7561437672e+00,-7.7342460187e+00,
     &-7.7211092751e+00,-7.7167990942e+00,-7.7214115512e+00,
     &-7.7341213404e+00,-7.7507037845e+00,-7.7665641466e+00,
     &-7.7824245087e+00,-7.7982848709e+00,-7.8141452330e+00/
      data qsp_y2/-6.0562005486e+01, 1.2112401097e+02,-4.4199292384e+02,
     &-7.1977920554e+01,-2.6998156099e+01,-1.4576627609e+01,
     &-8.8932735462e+00,-6.3502485266e+00,-4.7956317875e+00,
     &-4.0296457635e+00,-3.6001212385e+00,-3.4737052026e+00,
     &-3.5314697113e+00,-3.5841664322e+00,-3.5461955198e+00,
     &-1.6646487686e+00, 9.1041883424e-01,-2.4402976835e-01,
     & 6.5700239171e-02,-1.8771428333e-02, 9.3857141667e-03/
c
c QSP3 TC 1e5
      data qsp3_y/-9.7584274136e+00,-7.5072564741e+00,-6.6902003997e+00,
     &-6.4393262879e+00,-6.2907979028e+00,-6.1877344571e+00,
     &-6.1108566376e+00,-6.0513899788e+00,-6.0047277068e+00,
     &-5.9683378997e+00,-5.9407395959e+00,-5.9213341539e+00,
     &-5.9099536775e+00,-5.9068056593e+00,-5.9120754548e+00,
     &-5.9249277373e+00,-5.9411945133e+00,-5.9560807460e+00,
     &-5.9709669786e+00,-5.9858532113e+00,-6.0007394440e+00/
      data qsp3_y2/4.7857978428e+02,-9.5715956857e+02,-9.1817186257e+01,
     &-3.4408396646e+01,-1.6178971240e+01,-9.9915729551e+00,
     &-6.7002398198e+00,-4.9942534457e+00,-4.0532747173e+00,
     &-3.4465634449e+00,-3.2600794229e+00,-3.1759871833e+00,
     &-3.2958892839e+00,-3.3983553610e+00,-3.3134421521e+00,
     &-1.5458448305e+00, 1.3020370743e+00,-3.4899954650e-01,
     & 9.3961351751e-02,-2.6846100500e-02, 1.3423050250e-02/
c
c CHIANTI 8 HeI collisions, set TC=6e4
c
c    01        0.000 01 1S0
c    02   159856.078 03 3S1
c    03   166277.542 01 1S0
c    04   169086.870 05 3P2
c    05   169086.946 03 3P1
c    06   169087.934 01 3P0
c    07   171135.000 03 1P1
c
      data s_x/0.00d0,0.10d0,0.20d0,0.30d0,0.40d0,
     &0.50d0,0.60d0,0.70d0,0.80d0,0.90d0,1.00d0/
c     i,j,aji,TC: 01 02 1.73000d-04 5.00000d+04
C     data s12_y/0.248200d-01,0.661519d-01,0.685075d-01,0.668438d-01,
C    &0.655323d-01,0.637782d-01,0.606376d-01,0.550652d-01,
C    &0.456252d-01,0.295315d-01,0.829342d-03/
C     data s12_y2/0.000000d+00,-.608841d+01,0.967842d+00,-.194561d+00,
C    &0.217146d-01,-.157816d+00,-.222353d+00,-.411923d+00,
C    &-.450419d+00,-.177867d+01,0.000000d+00/
c     i,j,aji,TC: 01 03 5.09400d+01 5.00000d+04
      data s13_y/0.992700d-02,0.320098d-01,0.373748d-01,0.404687d-01,
     &0.438688d-01,0.475921d-01,0.511368d-01,0.542704d-01,
     &0.569714d-01,0.592824d-01,0.612600d-01/
      data s13_y2/0.000000d+00,-.258740d+01,0.318984d+00,-.512395d-01,
     &0.697093d-01,-.337198d-01,-.419377d-01,-.451578d-01,
     &-.370373d-01,-.407352d-01,0.000000d+00/
c     i,j,aji,TC: 01 04 3.93000d-01 5.00000d+04
C     data s14_y/0.148700d-02,0.946138d-02,0.138628d-01,0.176216d-01,
C    &0.207036d-01,0.229644d-01,0.246010d-01,0.258522d-01,
C    &0.268539d-01,0.276840d-01,0.283900d-01/
C     data s14_y2/0.000000d+00,-.552463d+00,0.660593d-01,-.973232d-01,
C    &-.828570d-01,-.639771d-01,-.357035d-01,-.245118d-01,
C    &-.159499d-01,-.146247d-01,0.000000d+00/
c     i,j,aji,TC: 01 05 2.33000d+02 5.00000d+04
C     data s15_y/0.107500d-02,0.565209d-02,0.836571d-02,0.105047d-01,
C    &0.126091d-01,0.146954d-01,0.163562d-01,0.168962d-01,
C    &0.155370d-01,0.109733d-01,0.400687d-03/
C     data s15_y2/0.000000d+00,-.275371d+00,-.165985d-01,-.298198d-02,
C    &0.774003d-02,-.388715d-01,-.107522d+00,-.203482d+00,
C    &-.218109d+00,-.846803d+00,0.000000d+00/
c     i,j,aji,TC: 01 06 0.00000d+00 5.00000d+04
C     data s16_y/0.346200d-03,0.188106d-02,0.279491d-02,0.349303d-02,
C    &0.422034d-02,0.496029d-02,0.551883d-02,0.566644d-02,
C    &0.517237d-02,0.362741d-02,0.129966d-03/
C     data s16_y2/0.000000d+00,-.903388d-01,-.112561d-01,0.593262d-02,
C    &0.503345d-02,-.184831d-01,-.399424d-01,-.683086d-01,
C    &-.718276d-01,-.274916d+00,0.000000d+00/
c     i,j,aji,TC: 01 07 1.80000d+09 5.00000d+04
      data s17_y/0.274300d-02,0.100499d-01,0.175928d-01,0.255418d-01,
     &0.343190d-01,0.448145d-01,0.592707d-01,0.982710d-01,
     &0.207894d+00,0.488473d+00,0.191756d+01/
      data s17_y2/0.000000d+00,0.285931d-01,0.273014d-01,0.105796d+00,
     &0.464563d-01,0.739346d+00,-.627440d+00,0.164969d+02,
     &-.229865d+02,0.178022d+03,0.000000d+00/

c
      u=t/1.d4
c relax!
c      if (u.lt.0.1d0) u=0.1d0
c      if (u.gt.10.d0) u=10.d0
c
      f=dsqrt(1.0d0/(t+epsilon))
      helos=0.0d0
      atom=zmap(2)
c
c old values:
C     a2s31s1=1.13d-4 !
C     a2p31s1=1.76d+2
C     at1s=1.76d+2+1.13d-4
C     a2s11s1=5.13d+1 !
Cc
C     ome1s12s3=6.87d-2 !
C     ome1s12s1=3.61d-2 !
C     ome1s12p3=2.27d-2 !
C     omep1=2.57d0 ! wrong!!
C     omes1=0.55d0*omep1 ! ditto
cc
c
      a2s31s1=1.73d-04  ! CHIANTI 8
      a2p31s1=2.33d+02+3.93000d-01  ! CHIANTI 8
      at1s=a2p31s1+a2s31s1
      a2s11s1=5.094d+01 ! 2 photon CHIANTI 8
c
      beta=(t/(t+5e4))
      omep1=fsplint(s_x,s17_y,s17_y2,11,beta)
      omes1=fsplint(s_x,s13_y,s13_y2,11,beta)
c
      reche=de*dh*zion(atom)*pop(2,atom)*rec(4,atom)
      frs1=reche*fusp(u)*(1.0d0-frs3(u))
      frp1=reche*(1.0d0-fusp(u))*(1.0d0-frs3(u))
c
c    ***COLLISIONAL RATES FROM GROUND STATE 1S
c
      abup=de*dh*zion(atom)*pop(1,atom)
      ara=abup*rka*f
c
      xs1=heien(2)/(rkb*t)
      xp1=heien(3)/(rkb*t)
c
      cos1=(ara*omes1*dexp(-xs1))! *fgaunt(1,1,xs1) omes1 is complete
      cop1=(ara*omep1*dexp(-xp1))! *fgaunt(1,1,xp1) omep1 is complete
c
      collrate2phe(atom)=cos1/abup
      ee2phe(atom)=heien(2)
c
c    ***RATES VIA TRIPLET STATES
c
      ra3=reche*frs3(u)*0.333d0
      tns3=ra3/(a2s31s1+(de*(qss(u)+qsp(u))))
      ra21=a2s31s1*tns3
c
      ras1=de*tns3*qss(u)
      rap1=de*tns3*qsp(u)
      rap3=de*tns3*qsp3(u)
c
c    ***TOTAL RATES AND APPROPRIATE COOLING
c
c     2S(T) > 1S(S)
c
      hein(1)=(heien(1)*ra21)*ifpi
      heibri(1)=hein(1)
c
c     2S(T) > 2S(S)
c
      hein(2)=(heien(2)*((ras1+frs1)+cos1))*ifpi

c
c     2S(T) > 2P(S)
c
      hein(3)=(heien(3)*((rap1+frp1)+cop1))*ifpi
      heibri(2)=hein(3)
c
c     2S(T) > 2P(T)
c
      hein(4)=((heien(4)*rap3)*ifpi)+hel0(6,j108m)
      heibri(3)=hein(4)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c New helium I effective recombination coeffs + collisional enhancement
c
c
      atom=zmap(2)
c
      coolf=0.d0
c
      do idx=1,nheislines
        heisbri(idx)=0.d0
      enddo
c
      do idx=1,nheitlines
        heitbri(idx)=0.d0
      enddo
c
      if (zion(atom)*pop(2,atom).gt.pzlimit) then
c
c First, collision contributions
c
        t4=1.0d-4*t
        if (t4.lt.0.01d0) t4=0.01d0
        if (t4.gt.5.0d0) t4=5.0d0
c
        do idx=1,nheis
          n=nheiscr(idx)
          cr(idx)=0.d0
          if (n.gt.0) then
            sum=0.d0
            do i=1,n
              ai=heiscra(i,idx)
              bi=heiscrb(i,idx)
              ci=heiscrc(i,idx)
              sum=sum+(ai*(t4**bi)*dexp(ci/t4))
            enddo
            cr(idx)=1.d0/(1.d0+((3552.d0*(t4**(-0.55d0)))/de))
            cr(idx)=cr(idx)*sum
          endif
c
c collisions of He I are now all in XR3DATA collisional cascades
c          cr(idx)=cr(idx)+1.d0
          cr(idx)=1.d0
c      write(*,*) t,idx,cr(idx)
        enddo
        abde=de*dh*zion(atom)*pop(2,atom)
        lnte=dlog(t)
        invt=1.d0/(epsilon+t)
c
c then singlet recomb
c
        do idx=1,nheislines
          a=heisreccoef(1,idx)
          b=heisreccoef(2,idx)
          c=heisreccoef(3,idx)
          d=heisreccoef(4,idx)
          emiss=(a+(b*lnte*lnte)+(c*lnte)+(d/lnte))*invt*1.0d-25
          emiss=dmax1(0.d0,emiss)
          heisbri(idx)=abde*emiss*ifpi
          upid=heisupid(idx)
          if (cr(upid).gt.1.0d0) then
            heisbri(idx)=heisbri(idx)*cr(upid)
            coolf=coolf+heisbri(idx)*(cr(upid)-1.d0)
          endif
        enddo
c
c then triplet colls
c
        do idx=1,nheit
          n=nheitcr(idx)
          cr(idx)=0.d0
          if (n.gt.0) then
            sum=0.d0
            do i=1,n
              ai=heitcra(i,idx)
              bi=heitcrb(i,idx)
              ci=heitcrc(i,idx)
              sum=sum+(ai*(t4**bi)*dexp(ci/t4))
            enddo
            cr(idx)=1.d0/(1.d0+((3552.d0*(t4**(-0.55d0)))/de))
            cr(idx)=cr(idx)*sum
          endif
c
c collisions of He I are now all in XR3DATA collisional cascades
c          cr(idx)=cr(idx)+1.d0
          cr(idx)=1.d0
c      write(*,*) t,idx,cr(idx)
        enddo
c
c then triplet recomb
c
        do idx=1,nheitlines
          a=heitreccoef(1,idx)
          b=heitreccoef(2,idx)
          c=heitreccoef(3,idx)
          d=heitreccoef(4,idx)
          emiss=(a+(b*lnte*lnte)+(c*lnte)+(d/lnte))*invt*1.0d-25
          emiss=dmax1(0.d0,emiss)
          heitbri(idx)=abde*emiss*ifpi
          upid=heitupid(idx)
          if (cr(upid).gt.1.0d0) then
            heitbri(idx)=heitbri(idx)*cr(upid)
            coolf=coolf+heitbri(idx)*(cr(upid)-1.d0)
          endif
        enddo
c
c      write(*,*) 'HeI coll losses: ',helos, coolf
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
