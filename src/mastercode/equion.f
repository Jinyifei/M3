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
c*******TO FIND IONISATION EQUILIBRIUM AT A GIVEN TEMP.& DENS.
c     FOR ALL ELEMENTS  ;  OUTPUT ELECTRONIC DENSITY  :  DE
c     AND THE FRACTIONAL ABUNDANCE OF THE DIFFERENT SPECIES
c     OF EACH ELEMENT  :  POP(6,11)  IN COMMON BLOCK /ELABU/
c     CALL SUBROUTINES IOHYD,IOBAL
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine equion (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 popzero(mxion, mxelem)
      real*8 t, de, dh, tstep, xhy, difma
      real*8 difm1, difm2
      real*8 treh, trea, zm, xhyf, dif1
      real*8 dia, dift, dfh, dft
      integer*4 nf, n, m, iel,inttemp
c
      character mod*4, nff*4, nel*4
c
      tstep=1.d37
      mod='EQUI'
      xhy=-1.0d0
      nf=7
      difma=1.d-3
      difm1=1.d-4
      difm2=0.08d0
      treh=1.d-8
c
c
      trea=1.d-5
      nff='ALL'
      zm=0.0d0
      do 10 iel=3,atypes
   10   zm=zm+zion(iel)
c
c
c    ***ITERATES TO FIND IONIC POPULATIONS AT EQUILIBRIUM
c
      do 40 m=1,7
        call copypop (pop, popzero)
        do 20 n=1,nf
c
          call iohyd (dh, xhy, t, tstep, de, xhyf, mod)
c
          call difpop (pop, popzero, treh, 1, dif1)
c
          dia=dmax1(0.0d0,dlog(5.0d0/(zm+1.d-5))/20.0d0)*pop(2,1)
          if (dif1.ge.dia) then
            nel='ALL'
          else if ((m.eq.1).and.(n.lt.nf)) then
            nel='HE'
          else if (m.gt.2) then
            nel='ALL'
          endif
c
          if (de.lt.pzlimit) de=pzlimit
          call iobal (mod, nel, de, dh, xhyf, t, tstep)
          call difpop (pop, popzero, trea, atypes, dift)
c
          call copypop (pop, popzero)
c
          if ((dift.lt.difm1).and.(nel.eq.'HE')) goto 30
          if ((dift.lt.difm2).and.(nel.eq.'ALL')) goto 30
c
   20   continue
c
   30   if (((nel.eq.'ALL').and.(dif1.lt.difma)).and.(dift.lt.difm2))
     &   goto 60
        call iohyd (dh, xhy, t, tstep, de, xhyf, mod)
        call iobal (mod, nff, de, dh, xhyf, t, tstep)
        inttemp=1
        call difpop (pop, popzero, treh, inttemp, dif1)
c
        call difpop (pop, popzero, trea, atypes, dift)
        if ((dif1.le.difma).and.(dift.le.difm2)) goto 60
c
   40 continue
c
      if ((dif1.le.difma).and.(dift.le.difm2)) goto 60
      dfh=dif1/difma
      dft=dift/difm2
C      if (dmax1(dfh,dft).gt.5.d0) then
C      write (*,50) dfh,dft
C   50 format(' Slow convergence for equil. ionisation:','DFH:'
C     &     ,1pg9.2,'   DFT:',1pg9.2)
C      endif
c
   60 continue
      return
      end
