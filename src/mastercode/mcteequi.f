      subroutine mcteequi (tei, tef, edens, hdens, nmod)
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
      real*8 tei, tef, edens, hdens
      real*8 tetmp, teHi, teLow
      real*8 reseng, resengHi, resengLow
c
      integer*4 n,nf
c
      character nmod*4
c
c
ccc      write (*,*) 'mcteequi'
c
      if ((nmod.ne.'EQUI').and.(nmod.ne.'TIM')) then
        write (*,10) nmod
   10 format(/' MODE IMPROPERLY SET FOR SUBR. TEEQUI  : ',a4)
        stop
      endif
c
c    Fudge which sets everything constant once T =30 K
c
c
      nf=25
      tetmp = tei
c
      reseng = 0.d0
      resengLow = 0.d0
      resengHi = 0.d0
      teLow = 0.d0
      teHi = 0.d0
c
      do 20 n=1,nf
c
        call equion (tetmp, edens, hdens)
        call cool (tetmp, edens, hdens)
c
        reseng = egain - eloss
c
        if (dabs(reseng).ge.0.01d0*egain) then
c           
           if (reseng.lt.0) then
              resengLow = reseng
              teLow = tetmp
              if (teHi.ne.0) then
                 tetmp = teLow - (teHi - teLow) * resengLow
     &                  /(resengHi - resengLow) 
              else
                 tetmp = tetmp * 0.83333333d0
              endif
           else
              resengHi = reseng
              teHi = tetmp
              if (teLow.ne.0) then
                 tetmp = teLow - (teHi - teLow) * resengLow
     &                  /(resengHi - resengLow) 
              else
                 tetmp = tetmp * 1.2d0
              endif
           endif
c
           if (n.ge.6) then 
              if (dabs(teHi - teLow).le.(0.002d0*(teHi+teLow))) then
                 if (reseng.lt.0) then
                    tetmp = teHi - 100.d0
                    teHi = 0.d0
                 else
                    tetmp = teLow + 100.d0
                    teLow = 0.d0
                 endif
              endif
           endif   
c
           if (tetmp.le.0) tetmp = 10.d0                   
c
        else  
           goto 30
        endif
c
   20 continue
c
      if (teHi.eq.0.d0) then 
         tetmp = teLow
      else
         if (teLow.eq.0.d0) then
            tetmp = teHi
         else
            tetmp = (teHi + teLow)*0.5d0
         endif
      endif
c 
   30 continue
c
c     ***USES FINAL TEMPERATURE
c 
      tef = tetmp
c
      call equion (tef, edens, hdens)

      call cool (tef, edens, hdens)
c
      return
c
      end