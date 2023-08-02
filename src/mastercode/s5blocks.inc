cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Shock 5 shared Variables
c     Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 dr,drdw,drup,dvdw,dvup,frdw,dt,absdr
      real*8 tdw,tlim,tup,xhii,bp0,bp1,psi
      real*8 t,de,dh,dv,netloss,rad,fi,wdilt0
      real*8 tt0,l0,l11,l12,l1,cmpf,ue0,ue1
      real*8 r0,p0,v0,tnloss,cspd,wmol,mb,malpha,pram,mageta,gaseta
      real*8 pgas,pmag,magparam,bmag,ve,en,ue,rdis,dva,mmach
      real*8 pop1(mxion,mxelem),totn,ionbar,abstime
      real*8 press,mu,rhotot,cltime,utime
      real*8 ptime,ctime,rtime,stime,eqtime
      real*8 atimefrac,cltimefrac,abtimefrac
      real*8 timend,teq,dheq,tl0,tl1,tt1
      real*8 band1,band2,band3,band4,band5,band6,bandall,wid
      real*8 meanion(mxelem),zi(mxelem)
      real*8 machnumber,alfvennumber
      real*8 prefield(mxinfph)
c
c Shock intial state kept for iterations
      real* 8  te_pre,te_pst,de_pre,de_pst,dh_pre,dh_pst
      real* 8  vs_pre,vs_pst,rh_pre,rh_pst,pr_pre,pr_pst
      real* 8  bm_pre,bm_pst,pop_pre(mxion,mxelem)
c Shock state from last iteration for convergence
      real* 8  te_pre0,te_pst0,de_pre0,de_pst0,dh_pre0,dh_pst0
      real* 8  vs_pre0,vs_pst0,rh_pre0,rh_pst0,pr_pre0,pr_pst0
      real* 8  psi0,cmpf0
      real* 8  bm_pre0,bm_pst0,pop_pre0(mxion,mxelem)
c Shock proto state
      real* 8  bm_neu,te_neu,de_neu,dh_neu,vs_neu,
     &         rh_neu,pr_neu,pop_neu(mxion,mxelem)
c
      common /shock5init/
     & te_pre,te_pst,de_pre,de_pst,dh_pre,dh_pst,
     & vs_pre,vs_pst,rh_pre,rh_pst,pr_pre,pr_pst,
     & bm_pre,bm_pst,pop_pre,
     & te_pre0,te_pst0,de_pre0,de_pst0,dh_pre0,dh_pst0,
     & vs_pre0,vs_pst0,rh_pre0,rh_pst0,pr_pre0,pr_pst0,
     & bm_pre0,bm_pst0,pop_pre0,psi0,cmpf0,
     & te_neu,de_neu,dh_neu,vs_neu,pop_neu,bm_neu,rh_neu,pr_neu
c
      common /shock5R/
     & dr,drdw,drup,dvdw,dvup,frdw,dt,absdr,
     & tdw,tlim,tup,xhii,bp0,bp1,psi,
     & t,de,dh,dv,netloss,rad,fi,wdilt0,
     & tt0,l0,l11,l12,l1,cmpf,ue0,ue1,mmach,
     & r0,p0,v0,tnloss,cspd,wmol,mb,malpha,pram,mageta,gaseta,
     & pgas,pmag,magparam,bmag,ve,en,ue,rdis,dva,machnumber,
     & meanion,zi,pop1,totn,ionbar,abstime,alfvennumber,
     & press,mu,rhotot,cltime,utime,
     & ptime,ctime,rtime,stime,eqtime,
     & atimefrac,cltimefrac,abtimefrac,
     & timend,teq,dheq,tl0,tl1,tt1,
     & band1,band2,band3,band4,band5,band6,bandall,wid,
     & prefield

c
c      real*8 feii_flux,feii_bri,feii_ratio,feii_bratio
c
      integer*4 ie,luop,lusp,lupb,lucl,lupc,jnorm
      integer*4 count,step,np,lurt,ludy,lual,lupf,lupt
      integer*4 luions(4)
      integer*4 converged, finalit
      common /shock5I/ ie,jnorm,luop,lusp,lupb,lucl,lupc,converged,
     & count,step,np,lurt,ludy,lual,lupf,luions,finalit,lupt
c
      character abundtitle*24
      character fn*32,filn(4)*32,fpc*32,fpt*32
      character specmode*4,rmod*4,imod*4,tsrmod*4,spmod*4
      character fl*32,lmod*4,ilgg*4,tmod*4,linemod*4,bandsmod*4
      character fd*32,fr*32,fa*32,fsp*32,allmod*4,ratmod*4,ispo*4
      character pfx*16,sfx*4,caller*4,wmod*4,pollfile*12,vmod*4
      character cht*32,pht*32,clt*32,dynmod*4,mypfx*16,subname*32
      character fpf*64,fpb*64,fcl*64,fclmod*4,stype*4,mtype*4
      common /shock5C/ abundtitle,fn,filn,fpc,fpt,
     & specmode,rmod,imod,tsrmod,spmod,
     & fl,lmod,ilgg,tmod,linemod,bandsmod,
     & fd,fr,fa,fsp,allmod,ratmod,ispo,
     & pfx,sfx,caller,wmod,pollfile,vmod,
     & cht,pht,clt,dynmod,mypfx,subname,
     & fpf,fpb,fcl,fclmod,stype,mtype
c
