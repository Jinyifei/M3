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
c     This routine returns the non-relativistic net compton
c     heating.
c
c     ref: Krolik, McKee and Tarter Ap. J. 249:422-442
c
c     Thanks to Wang Chi Lin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compton (t, de)
c
      include 'cblocks.inc'
c
c
      real*8 t, de
      real*8 sigmat,f1,eav
      real*8 sumfe, sumf
      real*8 energ,wid
      integer*4 i
c
c
      cmplos=0.d0
      cmpcool=0.d0
      cmpheat=0.d0
      if (linecoolmode.eq.1) return
c
c     relativistic (2hv/mc^2) term neglected
c
      sigmat=6.6524d-25
c
c     get totalflux sumf, and mean photon energy eav
c
      sumf=0.0d0
      sumfe=0.0d0
      do i=1,infph-1
        energ=ev*cphotev(i)
        wid=(photev(i+1)-photev(i))
        f1=tphot(i)*wid*evplk
        sumfe=sumfe+f1*energ
        sumf=sumf+f1
      enddo
c
      if (sumf.gt.0.d0) then
c
        eav=sumfe/sumf
c
c     get the loss rate...
c
        cmplos=((sigmat*sumf)/(me*cls*cls))*de*((4*rkb*t)-eav)
        cmpcool=((sigmat*sumf)/(me*cls*cls))*de*((4*rkb*t))
        cmpheat=((sigmat*sumf)/(me*cls*cls))*de*(-eav)
      endif
c
      if (expertmode.gt.0) then
        write (*,10) 'Te,Ctot,ht,cl:',t,cmplos,cmpcool,cmpheat
   10    format (a15,1p4e15.4)
      endif
c
      return
      end
