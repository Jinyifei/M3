cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012 Ralph S. Sutherland, Michael A. Dopita
c                       Luc Binette and Brent Groves
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     heating due to PAHS  : photo-electron ejection.
c
c     assumes tphot (J_nu) is up to date.
c
c     output in paheat (erg s^-1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hpahs (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh
      real*8 den,wid,pahrad,sigma
      real*8 phi,ryd4,b
      real*8 phots,photsn(mxinfph),photsi(mxinfph)
      real*8 jpe(pahi),jec(pahi)
      real*8 elmax,se,ec,qtmfac
      real*8 englow,engmax,engesc
      real*8 pahstate(pahi+2,pahi+1),pahsolve(pahi+1)
      real*8 pahheat,pahcool
c      real*8 sum,beta
      logical pahbin(mxinfph)
c
      integer*4 i,j
c
      paheat=0.d0
      paheng=0.d0
      pahheat=0.d0
      pahcool=0.d0
c
      if ((grainmode.eq.1).and.(pahmode.eq.1).and.(pahactive.eq.1))
     &then
c
c
c     assuming coronene (Nc=24)
c
        pahnc=468
c
c  Use Coronene Nc=24, radius ~ 5.85 A, S.A = (2/pi)* pi*a**2 (2/pi)
c   due to disc orientation sigma= 6.85E-15
c    from http://ois.nist.gov/pah/sp922_Detail.cfm?ID=362
c     new uses 10A Draine & Li Grain
c
c *4.d0*pi??
c
        sigma=(2/pi)*(10.d0*10.d0)*1.d-16
c
        pahrad=10.0d-8
c
        ryd4=4.d0*iph
c
        b=pahip(1)
c
        se=0.5d0
c
c
c  Maximum electon sticking rate
c
        elmax=de*se*dsqrt((8*rkb*t)/(pi*me))*sigma
c
c     Note: e^2/C=3.42 eV where C=2a/pi = capacitance of disk
c
        ec=(eesu**2*pi)/(2*pahrad)
        qtmfac=1.d0/(1.d0+(27d-08/pahrad)**0.75)
        phi=ec/(rkb*t)
        ec=ec/ev
c
        do j=1,infph-1
c
          pahbin(j)=.false.
          if ((photev(j).gt.b).and.(photev(j).lt.ryd4)) then
c
            den=cphote(j)
            wid=(photev(j+1)-photev(j))*evplk
            phots=fpi*tphot(j)*wid/den
            photsn(j)=pahnabs(j)*phots
            photsi(j)=pahiabs(j)*phots
            if ((photsn(j).gt.0).or.(photsi(j).gt.0)) then
              pahbin(j)=.true.
            endif
          endif
        enddo
c
c  Calculate electon collision rate and photoionisation rates
c
c     Negative PAHs
c
        do i=1,2
          jpe(i)=0.d0
c
c  Electron capture to that level
          jec(i+1)=elmax*dexp(pahion(i+1)*phi)
c
c  Photoionisation from that level
          do j=1,infph-1
            if (pahbin(j)) then
              if (photev(j).gt.pahip(i)) then
                jpe(i)=jpe(i)+photsn(j)*pahyield(j,i)
              endif
            endif
          enddo
        enddo
c
c     Neutral PAHs
c
        i=3
        jpe(i)=0.d0
c
c  Electron capture to that level
        jec(i+1)=elmax*(1.d0+(pahion(i+1)*phi))
c
c  Photoionisation from that level
        do j=1,infph-1
          if (pahbin(j)) then
            if (photev(j).gt.pahip(i)) then
              jpe(i)=jpe(i)+photsn(j)*pahyield(j,i)
            endif
          endif
        enddo
c
c     Positive PAHs
c
        do i=4,4
          jpe(i)=0.d0
c
c  Electron capture to that level
          jec(i+1)=elmax*(1.d0+(pahion(i+1)*phi))
c
c  Photoionisation from that level
          do j=1,infph-1
            if (pahbin(j)) then
              if (photev(j).gt.pahip(i)) then
                jpe(i)=jpe(i)+photsi(j)*pahyield(j,i)
              endif
            endif
          enddo
        enddo
c
c     Solve ionisation state of PAHm, np = 5 levels
c
        do i=2,6
          do j=1,7
            pahstate(j,i)=0.d0
          enddo
        enddo
c     sum f(Z)=1.d0
        do i=1,6
          pahstate(i,1)=1.d0
        enddo
c     Jpe-Jec=0
        do i=2,5
          pahstate(i-1,i)=jpe(i-1)
          pahstate(i,i)=-jec(i)
        enddo
c solve
c        call mdiag (5, pahstate, pahsolve)
        call  matsolve (pahstate, pahsolve, 5, 5)
c Calculate (+ve) ionisation fraction
        do i=1,pahi
          pahz(i)=pahsolve(i)
        enddo
c
c Calculate heating & cooling & energy absorbed
c
c Uses both Dopita Sutherland (2000) ApJ 539 742
c and Weingartner & Draine(2001) ApJS 134, 263
c engesc =average escape energy
c        = ef(E)/f(E) from W&D
c
        pahheat=0.d0
        pahcool=0.d0
c Negative PAHs
        do i=1,2
          englow=-(pahion(i)+1)*ec*qtmfac
          do j=1,infph-1
            den=cphotev(j)
            wid=(photev(j+1)-photev(j))*evplk
            paheng=paheng+pahz(i)*pahnabs(j)*fpi*tphot(j)*(1.d0-
     &       pahyield(j,i))*wid
c                 sum=sum+pahZ(i)*pahnabs(j)*tphot(j)*wid*fpi
            if ((pahbin(j)).and.(photev(j).gt.pahip(i))) then
              engmax=den-pahip(i)+englow
              engesc=0.5d0*(engmax+englow)
              pahheat=pahheat+pahz(i)*photsn(j)*pahyield(j,i)*engesc*ev
            endif
          enddo
          pahcool=pahcool+pahz(i+1)*elmax*rkb*t*(2-pahion(i+1)*phi)*
     &     dexp(pahion(i+1)*phi)
        enddo
c Neutral PAHs
        i=3
        englow=-(pahion(i)+1)*ec
        do j=1,infph-1
          den=cphotev(j)
          wid=(photev(j+1)-photev(j))*evplk
          paheng=paheng+pahz(i)*pahnabs(j)*fpi*tphot(j)*(1.d0-
     &     pahyield(j,i))*wid
c                 sum=sum+pahZ(i)*pahnabs(j)*tphot(j)*wid*fpi
          if ((pahbin(j)).and.(photev(j).gt.pahip(i))) then
            engmax=den-pahip(i)
c                beta=1.d0/(engmax-englow)**2
c                engesc=(engmax+englow-(4.d0*beta*englow**3)/3.d0)
c     &                  /(2.d0*(1.d0-2.d0*beta*englow**2))
            engesc=0.5d0*engmax*(engmax-2*englow)/(engmax-3*englow)
            pahheat=pahheat+pahz(i)*photsn(j)*pahyield(j,i)*engesc*ev
          endif
        enddo
        pahcool=pahcool+pahz(i+1)*elmax*rkb*t*(2+pahion(i+1)*phi)
c
c Positive PAHs
c
        i=4
        englow=-(pahion(i)+1)*ec
        do j=1,infph-1
          den=cphotev(j)
          wid=(photev(j+1)-photev(j))*evplk
          paheng=paheng+pahz(i)*pahiabs(j)*fpi*tphot(j)*(1.d0-
     &     pahyield(j,i))*wid
c                sum=sum+pahZ(i)*pahiabs(j)*tphot(j)*wid*fpi
          if ((pahbin(j)).and.(photev(j).gt.pahip(i))) then
            engmax=den-pahip(i)
c                beta=1.d0/(engmax-englow)**2
c                engesc=(engmax+englow-(4.d0*beta*englow**3)/3.d0)
c     &                /(2.d0*(1.d0-2.d0*beta*englow**2))
            engesc=0.5d0*engmax*(engmax-2*englow)/(engmax-3*englow)
            pahheat=pahheat+pahz(i)*photsi(j)*pahyield(j,i)*engesc*ev
c
          endif
        enddo
        pahcool=pahcool+pahz(i+1)*elmax*rkb*t*(2+pahion(i+1)*phi)
c
        i=5
        do j=1,infph-1
          wid=(photev(j+1)-photev(j))*evplk
          paheng=paheng+pahz(i)*pahiabs(j)*fpi*tphot(j)*(1.d0-
     &     pahyield(j,i))*wid
c                sum=sum+pahZ(i)*pahiabs(j)*tphot(j)*wid*4.d0*pi
        enddo
c
        paheat=(pahheat-pahcool)*pahfrac*dh
c
c  Testing grain charge - print output
c
c
c            write(*,*) 'PAHheat', pahfrac
c            write(*,100) (pahZ(i),i=1,5)
c            write(*,100) (Jpe(i),i=1,5)
c            write(*,100) (Jec(i),i=1,5)
c            write(*,100) paheng*dh*pahfrac,elmax,pahheat,pahcool,paheat
c           write(*,100) sum, paheng+pahheat, sum-(paheng+pahheat)
c
      endif
c
      if (dabs(pahheat).lt.epsilon) pahheat=0.d0
      if (dabs(pahcool).lt.pahcool) pahcool=0.d0
c
      return
c
      end
