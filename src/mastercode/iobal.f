cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     version 2.2, extended ionisation stages
c     Vector subsampling optimisation
c     rationalised loops
c     new charge exchange
c
c     RSS 4/91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DERIVES IONIC ABUNDANCES AFTER TIME : TSTEP
c     FOR THE ATOMIC ELEMENTS HEAVIER THAN HYDROGEN .
c     RETURNS RESULTS IN POP(6,11) AND DERIVATIVES IN DNDT(6,11)
c
c     WHEN ARGUMENT MOD SET TO 'EQUI' ,IT MAKES SURE THAT THE
c     SYSTEM REACHES EQUILIBRIUM (TSTEP IRRELEVANT)
c     THE ARGUMENT NEL DETERMINE IF IT COMPUTES FOR ALL ELEMENTS
c     OR ONLY FOR HELIUM .
c
c     THE AVERAGE ELECTRONIC DENSITY : DE  MUST BE GIVEN AS WELL
c     AS THE FRACTIONAL IONISATION OF HYDROGEN : XHY (OR POP2,1))
c     DEFAULT VALUE (USING POP(2,1)) ASSIGNED TO XHY IF SET TO -1
c     CALL SUBROUTINES SDIFEQ,ALLRATES,SPOTAP,IONAB,IONSEC
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine iobal (mod, nel, de, dh, xhy, t, tstep)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 abio,abne,de2,dep,dey,dpio,drec,dyn
      real*8 fhi,fhii
      real*8 ph,ph2,sigab,timt,yh
      real*8 de, dh, xhy, t, tstep
c
      real*8 reco(mxion), pion(mxion), ab(mxion), adndt(mxion)
      real*8 pionau(mxion),rech(mxion),pich(mxion)
      real*8 rc(mxion), pn(mxion), a(mxion), adnt(mxion), pa(mxion)
c
      integer*4 i,iel,ies,j,jdo,ji,jjj
      integer*4 nom,mnde,nde,mxde,mdelta,ion,at
c
      character jjmod*4, mod*4, nel*4
c
c           Functions
c
      real*8 feldens
c
      if (mod.eq.'EQUI') then
        timt=1.d37
      else
        timt=tstep
        if (timt.lt.0.0d0) then
          write (*,10)
   10         format(/' STOP : NEGATIVE TIME STEP'//)
          stop
        endif
      endif
c
      fhii=pop(2,1)
      fhi=pop(1,1)
      if ((xhy.lt.0.0d0).or.(xhy.gt.1.0d0)) goto 20
      fhii=xhy
      fhi=1.0d0-xhy
   20 continue
c
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod='ALL'
      call allrates (t, jjmod)
c
c
c    ***CALCULATE RATES DUE TO SECONDARY IONISATION
c
      call ionsec
c
      dyn=1.d12
c      dyn = 1.d10
      jjj=atypes
      if (nel.eq.'HE') jjj=2
c
      do 90 iel=2,jjj
        if (zion(iel).le.0.0d0) goto 90
        nde=maxion(iel)
        nom=0
        jdo=1
c
c    ***SET INITIAL ABUNDANCES IN AB(6) FOR ELEMENT : IEL
c
        do i=1,nde
          ab(i)=pop(i,iel)*zion(iel)*dh
        enddo
c
c     reentry point.
c
c
  100   continue
        do i=1,mxion
          rech(i)=0.d0
          pich(i)=0.d0
          reco(i)=0.d0
          pion(i)=0.d0
          pionau(i)=0.d0
        enddo
c
c    ***SET POPULATION RATES FOR EACH IONIC SPECIES OF ELEM : IEL
c
        do i=2,nde
          reco(i-1)=rec(i,iel)*de
        enddo
c
c    ***SET DEPOPULATION RATES FOR EACH IONIC SPECIES OF ELEM : IEL
c
        do i=1,nde-1
          pion(i)=rphot(i,iel)+(col(i,iel)*de)
        enddo
c
c     Auger rates, off diagnal terms.
c
        do i=1,nde-2
          pionau(i)=auphot(i,iel)
        enddo
c
c    ***ADD CORRECTION DUE TO SECONDARY IONISATION
c
        do i=1,min0(3,nde-1)
          pion(i)=pion(i)+rasec(i,iel)
        enddo
c
c    ***ADD CHARGE EXCHANGE RATES WITH HYDROGEN
c
        if (chargemode.eq.2) then
c
c     ancient mappings charge rates
c
          do 40 j=1,nchxold
            if (idint(charte(1,j)).ne.iel) goto 40
            ji=idint(charte(2,j))
            if (reco(ji).le.0.0d0) goto 40
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c     He III not considered
c
c
            ies=idint(charte(5,j))
            abne=(dh*zion(ies))*pop(1,ies)
            abio=(dh*zion(ies))*pop(2,ies)
            drec=abne*charte(4,j)
            dpio=abio*charte(3,j)
c
            rech(1)=reco(ji)+drec
            pich(1)=pion(ji)+dpio
c
c      *COMPRESS RECH TO FIT INTO DYNAMIC RANGE (IF NECESSARY)
c
            if (reco(ji).ge.(rech(1)/dyn)) goto 30
            if ((mod.eq.'EQUI').and.(nom.lt.1)) then
              jdo=2
              goto 40
            endif
            reco(ji)=reco(ji)*dyn
            pion(ji)=pich(1)*(reco(ji)/rech(1))
c
            goto 40
   30       pion(ji)=pich(1)
            reco(ji)=rech(1)
   40     continue
c
        endif
c
c     end old reactions
c
        if (chargemode.eq.1) then
c
c     Legacy AR1985 Rates
c     First, recombination reactions
c     with neutral H and He
c
          do at=1,nlegacychxr
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c
            if (chxrlegacyat(at).eq.iel) then
              if (chxrlegacy(at).gt.0.d0) then
                ies=chxrlegacyx(at)
                ion=chxrlegacyio(at)-1
                abne=dh*zion(ies)*pop(1,ies)
                drec=abne*chxrlegacy(at)
                rech(ion)=rech(ion)+drec
              endif
            endif
c
c     end at
c
          enddo
c
c     Ionising reactions
c     with ionised H and He
c
c
          do at=1,nlegacychxi
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c
            if (chxilegacyat(at).eq.iel) then
              if (chxilegacy(at).gt.0.d0) then
                ies=chxilegacyx(at)
                ion=chxilegacyio(at)
                abio=(dh*zion(ies))*pop(2,ies)
                dpio=abio*chxilegacy(at)
                pich(ion)=pich(ion)+dpio
              endif
            endif
c
c     end at
c
          enddo
c
c     renormalise to fit dynamic range
c
          do ion=1,maxion(iel)
c
            rech(ion)=reco(ion)+rech(ion)
            pich(ion)=pich(ion)+pion(ion)
c
            if (reco(ion).ge.(rech(ion)/dyn)) goto 50
            if ((mod.eq.'EQUI').and.(nom.lt.1)) then
              jdo=2
              goto 60
            endif
            reco(ion)=reco(ion)*dyn
            pion(ion)=pich(ion)*(reco(ion)/rech(ion))
c
            goto 60
   50       pion(ion)=pich(ion)
            reco(ion)=rech(ion)
   60       continue
          enddo
c
        endif
c
        if (chargemode.eq.0) then
c
c     New Kingdon Rates
c     First, recombination reactions
c     with neutral H and He
c
          do at=1,nchxr
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c
            if (chxrat(at).eq.iel) then
              if (chxr(at).gt.0.d0) then
                ies=chxrx(at)
                ion=chxrio(at)-1
                abne=dh*zion(ies)*pop(1,ies)
                drec=abne*chxr(at)
                rech(ion)=rech(ion)+drec
              endif
            endif
c
c     end at
c
          enddo
c
c     Ionising reactions
c     with ionised H and He
c
c
          do at=1,nchxi
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c
            if (chxiat(at).eq.iel) then
              if (chxi(at).gt.0.d0) then
                ies=chxix(at)
                ion=chxiio(at)
                abio=(dh*zion(ies))*pop(2,ies)
                dpio=abio*chxi(at)
                pich(ion)=pich(ion)+dpio
              endif
            endif
c
c     end at
c
          enddo
c
c     renormalise to fit dynamic range
c
          do ion=1,maxion(iel)
c
            rech(ion)=reco(ion)+rech(ion)
            pich(ion)=pich(ion)+pion(ion)
c
            if (reco(ion).ge.(rech(ion)/dyn)) goto 70
            if ((mod.eq.'EQUI').and.(nom.lt.1)) then
              jdo=2
              goto 80
            endif
            reco(ion)=reco(ion)*dyn
            pion(ion)=pich(ion)*(reco(ion)/rech(ion))
c
            goto 80
   70       pion(ion)=pich(ion)
            reco(ion)=rech(ion)
   80       continue
          enddo
c
c     end char exchange rates
c
        endif
c
c
c    ***ALTER RECOMB. RATES FOR HE WHEN ON THE SPOT APPROX. USED
c
        if ((jspot.eq.'YES').and.(iel.eq.zmap(2))) then
c
          call spotap (de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
          reco(1)=reco(1)-((de*(1.0d0-yh))*(rec(2,2)-rec(4,2)))
          reco(2)=reco(2)-(de*(rec(3,2)-rec(5,2)))
c
        endif
c
c    ***COMPUTES NEW IONIC ABUNDANCES FOR ELEMENT : IEL
c
c     first calculation done in full
c
c      if (nom.eq.0) then
c
c
c      call ionab(reco, pion, pionau, ab, adndt, nde, timt)
c
c      else
c
c
c      if (mapz(iel).eq.2) then
c         write(*,*) 'Helium rates (t,de,dh):',t,de,dh
c         do i = 1,nde
c            write(*,*) rom(i),pion(i),reco(i)
c         enddo
c      endif
c
        mnde=1
        mxde=2
c
c     find max and min ionisation stages that matter
c
c
        do i=1,nde
          if (pop(i,iel).ge.1.d-6) mxde=i
          if (pop(nde-i+1,iel).ge.1.d-6) mnde=nde-i+1
          if (pion(i).gt.0.d0) then
            if (reco(i)/pion(i).lt.1.d3) mxde=i
          endif
          if (pionau(i).gt.0.d0) then
            if (reco(i)/pionau(i).lt.1.d3) mxde=i
          endif
c
c     pull up the minimun, esp for Nickel and Iron with their
c     fairly flat rates
c
c            if (reco(nde-i+1).gt.0.d0) then
c               if (pion(nde-i+1)/reco(nde-i+1).lt.1.d3) mnde = nde-i+1
c            endif
        enddo
c
        mxde=min0(mxde+1,nde)
        mnde=max0(mnde-1,1)
c
c         if  (iel.eq.2) write(*,*) elem(iel),mnde,mxde
c         write(*,*) 'Min rates:',pion(mnde),reco(mnde)
c         write(*,*) 'Max rates:',pion(mxde),reco(mxde)
c         write(*,*) 'min/max pops:',pop(mnde,iel),pop(mxde,iel)
c
        mdelta=mxde-mnde+1
c
c     copy section into ab
c
        sigab=0.d0
        do i=1,mdelta
          a(i)=ab(i+mnde-1)
          sigab=sigab+pop(i+mnde-1,iel)
          pn(i)=pion(i+mnde-1)
          rc(i)=reco(i+mnde-1)
          pa(i)=pionau(i+mnde-1)
        enddo
        if (sigab.le.0.99d0) then
          write (*,*) 'IOBAL VECTOR SAMPLING OPT. FAILED'
          write (*,*) 'Significant population missed:',1.d0-sigab
          write (*,*) elem(iel),'min:',mnde,'Max:',mxde
          do i=1,maxion(iel)
          write(*,*) i,pop(i,iel),ab(i),pion(i),reco(i),pionau(i)
          enddo
          stop
        endif
c
c     do ionab with shortened vector
c
        call ionab (rc, pn, pa, a, adnt, mdelta, timt)
c
c     copy results back into their correct positions
c
        do i=1,mdelta
          adndt(i+mnde-1)=adnt(i)
          ab(i+mnde-1)=a(i)
        enddo
c
        if (mnde.gt.1) then
          do i=1,mnde-1
            adndt(i)=0.0d0
            ab(i)=0.0d0
          enddo
        endif
c
        if (mxde.lt.maxion(iel)) then
          do i=mxde+1,maxion(iel)
            adndt(i)=0.d0
            ab(i)=0.0d0
          enddo
        endif
c
c
c     endif
c
c    ***CHECK IF IT HAS TO RE-DO COMPUTATION TO REACH EQUILIBRIUM
c
        nom=nom+1
        if ((mod.eq.'EQUI').and.(nom.lt.jdo)) goto 100
c
c    ***COPY FINAL IONIC POPULATIONS AND TIME DERIVATIVE FOR EL. IEL
c
        do i=1,nde
          dndt(i,iel)=adndt(i)
          pop(i,iel)=(ab(i)/(dh*zion(iel)))
          if (pop(i,iel).lt.pzlimit) then
            pop(i,iel)=0.d0
            dndt(i,iel)=0.d0
          endif
        enddo
c
   90 continue
c
c
c    ***UPDATES ELECTRONIC DENSITY
c
      de=feldens(dh,pop)
c
      return
      end
