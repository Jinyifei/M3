cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE HEATING DUE TO PHOTOIONISATION
c     USING IONISATION FRACTION OF HYDROGEN, IT INTERPOLATES
c     THE RATES HEAPH(5,ion,elem) CALCULATED IN SUBR. PHION
c
c
      subroutine pheat (de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 aa,bb,cc,del,xe,xmi,pz
      real*8 qval,temp,val
      real*8 x1,xc,xlast,xlo,y1,y2,y3,ztot
      real*8 de, dh , quad
c
      integer*4 i,i1,i2,i3,ic,iel,ion
c
      quad(del,y1,y2,y3)=y1+((del*((((4.0d0*y2)-(3.0d0*y1))-y3)+(((y3-
     &(2.0d0*y2))+y1)*del)))*0.5d0)
c
      pgain=0.0d0
      xmi=2.d-5
      ztot=1.0d0
      do i=2,atypes
        ztot=ztot+zion(i)
      enddo
c
c    ***USES INTERPOLATION , xe < 1.0
c
      xe=dmin1(1.0d0,dmax1((de/dh)/ztot,xmi))
      if (xe.lt.0.999d0) then
        xlo=dlog10(xe+epsilon)
        xlast=-4.0d0
        x1=dint(xlo)
        xc=x1
        if (xc.lt.xlast) xc=xlast
        if (x1.lt.-2.0d0) x1=-2.0d0
c
        del=x1-xlo
c
        i1=idnint(1.0d0+dabs(x1))
        i1=min(3,max(1,i1))
        i2=i1+1
        i3=i1+2
        ic=idnint(1.0d0+dabs(xc))
        do iel=1,atypes
          do ion=1,maxion(iel)-1
            pz=zion(iel)*pop(ion,iel)
            if (pz.gt.pzlimit) then
              aa=heaph(i1,ion,iel)
              bb=heaph(i2,ion,iel)
              cc=heaph(i3,ion,iel)
              temp=0.d0
              if (cc.gt.0.0d0) then
                qval=dexp(quad(del,dlog(aa),dlog(bb),dlog(cc)))
                val=dmax1(0.0d0,dmin1(qval,heaph(ic,ion,iel)))
                temp=(dh*pz*val)
              else if (aa.gt.0.0d0) then
                qval=quad(del,aa,bb,cc)
                val=dmax1(0.0d0,dmin1(qval,heaph(ic,ion,iel)))
                temp=(dh*pz*val)
              endif
              pgain=pgain+temp
              heatz(iel)=heatz(iel)+temp
              heatzion(ion,iel)=heatzion(ion,iel)+temp
            endif
          enddo
        enddo
c
      else
c
c    ***NO INTERPOLATION BECAUSE FHII ~1.0
c
        do iel=1,atypes
          do ion=1,maxion(iel)-1
            pz=zion(iel)*pop(ion,iel)
            if (pz.gt.pzlimit) then
              temp=dh*pz*heaph(1,ion,iel)
              pgain=pgain+temp
              heatz(iel)=heatz(iel)+temp
              heatzion(ion,iel)=heatzion(ion,iel)+temp
            endif
          enddo
        enddo
      endif
c
      if (dabs(pgain).lt.epsilon) pgain=0.d0
c
      return
      end
