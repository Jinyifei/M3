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
c*******CALCULATES NET GAIN ARISING FROM HEATING DUE TO THE
c     ON THE SPOT APPROX. MINUS THE COOLING DUE TO RECOMB.
c     CALL SUBROUTINE SPOTAP
c
c     **REFERENCE : SEATON,M.J. (1959) MNRAS V.119,P.81
c                 HUMMER,D.G. AND SEATON,M.J. (1964)
c                    MNRAS,V.127,P.230
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine netgain (t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t, de, dh, rcfa,rcfb,z,uz2
      real*8 aa,ab,betr,de2,dep,dey,en,fhi
      real*8 ph,ph2,rat,rc,rg,spo,u,yh,rheat
c
      integer*4 i,j
c
      rcfa(uz2)=0.7714d0*(uz2**(-0.0390d0))
      rcfb(uz2)=0.6412d0*(uz2**(-0.1146d0))
c
      u=dmax1(t,10.d0)*1.0d-4
c
      en=rkb*t
      betr=0.0d0
      rheat=0.d0
      rg=0.0d0
c
c    ***IF NECESSARY,ADS HEATING DUE TO THE ON THE SPOT APPROX.
c
      if (jspot.eq.'NO') goto 20
      fhi=pop(1,1)
      call spotap (de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
      aa=(dh*zion(1))*pop(2,1)
      rc=aa*rec(3,1)
      aa=(dh*zion(2))*pop(2,2)
      rg=rg+((aa*de)*((dep*rec(4,2))+(dey*(rec(2,2)-rec(4,2)))))
      spo=aa*((yh*(rec(2,2)-rec(4,2)))+(ph*rec(4,2)))
      aa=(dh*zion(2))*pop(3,2)
      spo=spo+((aa*ph2)*rec(5,2))
      rat=spo/(rc+1.d-30)
      if (rat.le.1.0d0) goto 10
      spo=spo/rat
      de2=de2/rat
   10 rg=rg+(((aa*de)*de2)*rec(5,2))
   20 continue
c
c    ***FINDS COOLING/HEATING DUE TO RECOMBINATION
c
      do 40 i=1,atypes
        do 30 j=1,maxion(i)-1
          z=dble(j)
          uz2=u/(z*z)
          ab=zion(i)*pop(j+1,i)
          if (ab.gt.pzlimit) then
            if ((i.le.2).and.(jspot.eq.'YES')) then
              rheat=(((de*dh*ab)*rcfb(uz2))*rec(j+maxion(i),i))
            else
              rheat=(((de*dh*ab)*rcfa(uz2))*rec(j+1,i))
            endif
            heatz(i)=heatz(i)-(en*rheat)
            heatzion(j+1,i)=heatzion(j+1,i)-(en*rheat)
            betr=betr+rheat
          endif
   30   continue
   40 continue
c
      rngain=(-(en*betr))+rg
c
      if (dabs(rngain).lt.epsilon) rngain=0.d0
c
      return
      end
