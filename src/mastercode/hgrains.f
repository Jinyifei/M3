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
c     Heating due to grains: grain photo-electron ejection.
c                          : +PAH contribution
c                          : +collisional grain cooling
c
c     Assumes tphot (J_nu) is up to date.
c
c     Output in gheat, gcool, paheat
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hgrains (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh
      real*8 den,wid,phots(mxinfph),b,je,acc_h,yi
      real*8 eu,jec,jpr,elmax,prmax,phi
      real*8 eng,engmax,englow
      real*8 je1,jec1,jpr1,rydmax,qabs
      real*8 pmax,f,es
      real*8 gam,del,delmax
      real*8 sigma,mred,dgrad
      real*8 se,sp,spf,rkt,invrkt
      real*8 u, u1, u2, umid, y1, y2, ymid
c
      integer minbin,maxbin
      integer*4 i,j,nbins,check,usolv,k,dtype
c
   10 format(5(1pg11.4,1x))
c
c     Grain initialization
c
      do dtype=1,numtypes
        avgpot(dtype)=0.d0
      enddo
      gheat=0.d0
      gcool=0.d0
c
      if (linecoolmode.eq.1) return
c
      if (grainmode.eq.1) then
c
c     Init neutral grains.
c
        if (yinf.ne.0.d0) then
c     Some global variables are copied to local variables for readibilty.
c     Standard parameters: B (cutoff) and Yinf (yield factor),
c     (B in eV), Qabs as fraction of geometric size,
c     acc_H : accomodation factor.
          acc_h=haccom
          b=bgrain
          yi=yinf
c     Sticking factors for electrons (se) and protons (sp)
          se=segrain*dexp(-t/2.0d5)
          sp=spgrain
c     Electron and proton currents, less phi factor.
c     These need to be multiplied by phi factor for different
c     potentials.
          elmax=-de*se*dsqrt((8*rkb*t)/(pi*me))
          prmax=pop(2,1)*dh*sp*dsqrt((8*rkb*t)/(pi*mp))
c
c     Cooling collision rates. Separate from grains.
c
          delmax=de*dh*((2.d0*rkb*t)**(1.5d0))/dsqrt(pi*me)
          mred=dsqrt(me/mh)
          rkt=rkb*t
          invrkt=1.d0/rkt
c
c     Do for all types and grainsize distribution.
c
          do 70 dtype=1,numtypes
            do 60 k=mindust(dtype),maxdust(dtype)
c
              dgrad=gradedge(k+1)-gradedge(k)
              sigma=dustsig(k,dtype)*dgrad
              rydmax=iph
              nbins=0
              pmax=0.d0
c
c      total = 0.d0
c
c     Calculate phots(inl) as same thruout (only diff bins accessed)
c     and get totalphot and maxphot
c
              minbin=0
              maxbin=0
              do i=1,infph-1
                if ((cphotev(i).gt.b)) then
c
                  qabs=absorp(i,k,dtype)
c
                  if (qabs.gt.0.d0) then
                    if (minbin.eq.0) minbin=i
                    den=cphotev(i)
                    wid=(photev(i+1)-photev(i))*evplk
c
c     only look where absorption will occur
c
                    phots(i)=qabs*fpi*tphot(i)*wid/(den*ev)
                    phots(i)=yi*(1.d0-b/den)*phots(i)
c
c      sets max bin to last bin where qabs > 0
                    if (phots(i).gt.pmax) pmax=phots(i)
                    maxbin=i
c
c           total = total+ phots(i)
c
c           end Qabs>0
c
                  endif
c
                endif
              enddo
c
c     write(*,*) photev(minbin), photev(maxbin)
c
c     Now solve for charge balance to get grain potential
c     photo-electron rate:
c
c     Use binary search method to get grain charge such that
c     f(U)=e-sticking + p-sticking + photoelectric=0
c
              check=0
c              nr=rydmax
c              nr=10*((nr/10)+1)
c              u1=-1.d0*dble(nr)
c              u2=1.d0*dble(nr)
              u1=-2.d0
              u2=5.d0
c     Get f(U1), check that its +ve, 0, or -ve.
c     If +ve continue with binary search. If -ve, double U1 and try again.
c     If 0 (ie U=U1), exit search.
c     (Note: f(U) only equals zero at one point)
   20         u=u1
              eu=ev*u1
              phi=eu*invrkt
c     U is -ve
              jec1=elmax*dexp(phi)+epsilon
              jpr1=prmax*(1.d0-phi)
c     Get photo currents.
c     Loop through whole vector in case vector changes one day.
c     Only want bins below 8 IPH and above cutoff.
              je1=0.d0
              do j=minbin,maxbin-1
c     As U -ve, all photoelectrons escape (spf=1).
                je1=je1+phots(j)
              enddo
              y1=je1+jpr1+jec1
c               write (69,*) 'U1'
c               write (69,100) U, je1, jpr1, jec1, y1
              if ((y1.lt.0.d0).and.(check.lt.11)) then
                u2=u1
                u1=1.5d0*u1
                check=check+1
                goto 20
              else if (dabs(y1).lt.1.0d-4) then
                u2=u1
                goto 50
              else if (check.ge.11) then
                write (*,*) 'Grain Charge Error1:lower limit'
                write (*,10) u2,u1,je1,jpr1,jec1
                write (*,10) elmax,prmax,phi
                write (*,10) de,t,se
                stop
              endif
c
c    Get f(U2), check that its +ve,0, or -ve.
c    If -ve, continue with binary search. If +ve, double U2 and try again.
c    If 0 (ie U=U2), exit search
              if (u2.lt.0d0) goto 40
c
              check=0
   30         u=u2
              eu=ev*u2
              phi=eu*invrkt
c  U +ve.
              jec1=elmax*(1.d0+phi)
              jpr1=prmax*dexp(-phi)+epsilon
c
              je1=0.d0
c
              do j=minbin,maxbin-1
                if (cphotev(j).gt.(b+u)) then
                  eng=cphotev(j)
                  engmax=(eng-b)
                  englow=u
c   U +ve, so some fraction of photoelectrons can't escape potential.
                  spf=1.d0
                  if (englow.gt.(0.5d0*engmax)) then
                    es=((engmax-englow)/engmax)
                    spf=4.d0*es*es
                    spf=dmax1(0.d0,dmin1(spf,1.d0))
                  else if (englow.gt.0.d0) then
                    es=(englow/engmax)
                    spf=1.d0-4.d0*es*es
                    spf=dmax1(0.d0,dmin1(spf,1.d0))
                  endif
                  je1=je1+phots(j)*spf
                endif
              enddo
c
              y2=je1+jpr1+jec1
              if ((y2.gt.0.d0).and.(check.lt.11)) then
                u1=u2
                u2=1.75d0*u2
                check=check+1
                goto 30
              else if (dabs(y2).lt.1.0d-4) then
                u1=u2
                goto 50
              else if (check.ge.11) then
                print *,'Grain Charge Error2: upper limit'
                write (*,10) u2,u1,je1,jpr1,jec1
                write (*,10) elmax,prmax,phi
                write (*,10) de,t,se
                stop
              endif
c
c     Now that F(U1) > 0 > F(U2) (ie U1 < U < U2) use binary search to find U
c     to several decimal places.
   40         continue
c              if (u1.lt.-2.d0) write(*,*) 'bracketed',u1,u2
c
              do usolv=1,12
                umid=(u1+u2)*0.5d0
                u=umid
                eu=ev*u
                phi=eu*invrkt
                if (phi.ge.0.d0) then
                  jec1=elmax*(1.d0+phi)
                  jpr1=prmax*dexp(-phi)+epsilon
                else
                  jec1=elmax*dexp(phi)+epsilon
                  jpr1=prmax*(1.d0-phi)
                endif
                je1=0.d0
                do j=minbin,maxbin-1
c     ephot > B due to Minbin.
                  if ((cphotev(j).gt.(b+u))) then
                    eng=cphotev(j)
                    engmax=(eng-b)
                    englow=u
                    spf=1.d0
                    if (englow.gt.(0.5d0*engmax)) then
                      es=((engmax-englow)/engmax)
                      spf=4.d0*es*es
                      spf=dmax1(0.d0,dmin1(spf,1.d0))
                    else if (englow.gt.0.d0) then
                      es=(englow/engmax)
                      spf=1.d0-(4.d0*es*es)
                      spf=dmax1(0.d0,dmin1(spf,1.d0))
                    endif
                    je1=je1+phots(j)*spf
                  endif
                enddo
c
                ymid=je1+jpr1+jec1
c
c                  write (69,100) U, je1, jpr1, jec1, ymid
c
                f=dabs(ymid)/(dabs(je1)+dabs(jpr1)+dabs(jec1))
                if (f.lt.1.0d-3) goto 50
c
                if (ymid.gt.0.d0) then
c  Too much PE heating
                  u1=umid
                else if (ymid.lt.0.d0) then
c  Electron current too large
                  u2=umid
                else
                  goto 50
                endif
c
              enddo
c
   50         continue
c
              u=(u1+u2)*0.5d0
c
c            if (Usolv.lt.10) write(*,*) Usolv,f,u,u1,u2
c
c     Then grainpot(k,dtype) = U.
c
c
c     Have balance potential in U, get heating.
c
              eu=ev*u
              phi=eu*invrkt
c
c               write (*,*) 'grains:',U,eu,phi
c
              if (phi.ge.0.d0) then
                jec=elmax*(1.d0+phi)
                jpr=prmax*dexp(-phi)+epsilon
              else
                jec=elmax*dexp(phi)+epsilon
                jpr=prmax*(1.d0-phi)
              endif
c
c     Get heating term, (ergs not photons).
c
              je=0.d0
              gam=0.d0
              do j=minbin,maxbin-1
                if ((photev(j).gt.(b+u))) then
                  den=cphotev(j)
                  engmax=(den-(b+u))
                  englow=-u
c
c     Calculate average e- escape energy.
c     (Note: eng=eng/spf, phots=phots*spf but cancelled out here)
c
                  if (englow.ge.0) then
c U -ve
                    spf=1.d0
                    eng=0.5d0*(engmax+englow)
                  else if ((engmax+englow).gt.0) then
c  engmax > abs(englow)
                    spf=1.d0-2.d0*(englow/(engmax-englow))**2
                    eng=0.5d0*(-4.d0*englow**3/(3.d0*(engmax-englow)**2)
     &               +engmax+englow)
                  else
c  engmax < abs(englow)
                    spf=2.d0*(engmax/(engmax-englow))**2
                    eng=2.d0*engmax**3/(3.d0*(engmax-englow)**2)
                  endif
                  gam=gam+eng*ev*phots(j)
                  je=je+phots(j)*spf
                endif
              enddo
c
              gam=gam*dh*sigma
c
c     Collisional cooling.
c
c     Have phi so...
c
              del=delmax*sigma
              if (phi.ge.0.d0) then
                del=del*(se*(2.d0+phi)+mred*dexp(-phi)*(acc_h+acc_h+phi)
     &           )
              endif
              if (phi.lt.0.d0) then
                del=del*((se*(2.d0-phi)*dexp(phi))+(mred*acc_h*(2.d0-
     &           phi)))
              endif
c
c     Grain size/type.
c
              grainpot(k,dtype)=u
              avgpot(dtype)=avgpot(dtype)+u
              gheat=gheat+gam
              gcool=gcool+del
c
c               write(*,'(3(1pg12.5,x))') grainrad(k), U, U/grainrad(k)
c  End loop over grain sizes
c
   60       continue
c
            avgpot(dtype)=avgpot(dtype)/dble(maxdust(dtype)-
     &       mindust(dtype)+1)
c
c            write (*,*) dtype, avgpot(dtype)
c End loop over grain types
c
   70     continue
        endif
      endif
c
c      write (*,1069) gcool, gheat
c 1069 format ('graincool: ',1pg14.7,' grainheat: ',1pg14.7)
c
      if (dabs(gcool).lt.epsilon) gcool=0.d0
      if (dabs(gheat).lt.epsilon) gheat=0.d0
cc
      return
c
      end
