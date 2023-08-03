cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******TO COMPUTE RADIATIVE AND DIELECTRONIC
c     RECOMBINATION RATES ; UNITS : CM**3 SEC**-1
c     RETURNS VECTOR REC(ion,atom) IN COMMON BLOCK /RECOMB/
c
c     *FOR H AND HE , RATES IN REC(3,1),REC(4,2) AND REC(5,2)
c     FOR THE ON THE SPOT APPROXIMATION (WHEN USED)
c     REC(4,2) is always used for HeI spectrum and heilos
c     and is essential even when OTS is not used.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recom (t)
c
      include 'cblocks.inc'
c
      real*8 t,u,u4,ar,xr
      real*8 ad,bd,ta,tb
      real*8 tp,tp4
      real*8 a,b,c,d,e,f,r
      real*8 ch,ct,et
      real*8 fcor,adlt
      real*8 rate,drate
c      real*8 arfc
c      real*8 adfc
c
      real*8 en(8),co(8)
      real*8 a1,a2,a3,a4,a5,a6
      real*8 tt3,tt4,tt5,tt6,t32
c
      integer*4 ion,atom,i,ir,maxio,chh,isos
      integer*4 rrtype,drtype
c
c      Old code functions
c
      real*8 arf, adf, adltf, hifa, flowt, fkramer, fheirec,fnoradrec
c
      arf(u4,r,e)=r*(u4**(-e))
c
      adf(u,ad,bd,ta,tb)=(((1.0d0+(bd*dexp(-(tb/u))))*dexp(-(ta/u)))*ad)
     &*(u**(-1.5d0))
c
      adltf(u4,a,b,c,d,f)=((1.0d-12*((((a/u4)+b)+(c*u4))+((d*u4)*u4)))*
     &(u4**(-1.5d0)))*dexp(-(f/u4))
c
      hifa(u4)=((0.4288d0+(0.5d0*dlog(1.5789d1/u4)))+(0.469d0*
     &((1.5789d1/u4)**(-(0.33333333d0)))))/1.9954d0
c
      flowt(u4,ct,et)=dmax1(1.d0,((u4/dsqrt(ct))**(-dmax1(0.0d0,et-
     &0.5d0)))/dmax1(1.d0,hifa(u4/dsqrt(ct))))
c
c     keep t positive
c
      tp=dmax1(10.d0,t)
c
      if (usekappa) then
c
c  correct energy temp if kappa is in use
c
        tp=tp*(kappa-1.5d0)/kappa
c
      endif
c
      tp4=tp*1.d-4
      t32=tp**(-1.5d0)
c
      if ((recommode.eq.0).or.(recommode.eq.2)) then
c
c new std recombination calculations
c
c      write(*,*) 'recommode=0, new RR code'
c
        do ir=1,nrec
          rate=0.d0
          atom=atrec(ir)
          ion=ionrec(ir)
          isos=mapz(atom)-ion+1
          if (recommode.eq.2) then
            if (noradid(ion,atom).ne.0) then
              rate=fnoradrec(ion,atom,tp)
              goto 10
            endif
          endif
C           if (mapz(atom).eq.6) then
C             rate=fnoradrec(ion,atom,tp)
C             goto 10
C           endif
C           if (mapz(atom).eq.7) then
C             rate=fnoradrec(ion,atom,tp)
C             goto 10
C           endif
C           if (mapz(atom).eq.8) then
C             rate=fnoradrec(ion,atom,tp)
C             goto 10
C           endif
C           endif
C         endif
c
c fe II-VI are always NORAD totals if avaiable, even recommode!=2
c
          if ((mapz(atom).eq.26).and.(ion.le.6)) then
            if (noradid(ion,atom).ne.0) then
              rate=fnoradrec(ion,atom,tp)
              goto 10
            endif
          endif
c
          if (isos.eq.0) then
c
c Nahar-NORAD total hydrogenic recombnination
c includes 3,1 and 5,2 rates
c
            rate=fkramer(ion,atom,tp)
            goto 10
          endif
          if ((mapz(atom).eq.2).and.(ion.eq.2)) then
c
c Nahar-NORAD total HeI recombnination
c includes 4,2 rate
c
            rate=fheirec(ion,atom,tp)
            goto 10
          endif
          rrtype=typerrec(ir)
          if (rrtype.gt.0) then
            a1=arrec(1,ir)
            a2=arrec(2,ir)
            if (rrtype.eq.2) then
              rate=a1*tp4**(-a2)
              goto 10
            endif
            a3=arrec(3,ir)
            if (rrtype.eq.3) then
              rate=a1*tp4**(-(a2+(a3*dlog10(tp4))))
              goto 10
            endif
            a4=arrec(4,ir)
            if (rrtype.eq.4) then
              tt3=dsqrt(tp/a3)
              tt4=dsqrt(tp/a4)
              rate=a1/(tt3*((tt3+1.d0)**(1.d0-a2))*((tt4+1.d0)**(1.d0+
     &         a2)))
              goto 10
            endif
            if (rrtype.eq.6) then
              a5=arrec(5,ir)
              a6=arrec(6,ir)
              tt4=sqrt(tp/a4)
              tt5=sqrt(tp/a5)
              tt6=a2+(a3*dexp(-a6/tp))
              rate=a1/(tt4*((tt4+1.d0)**(1.d0-tt6))*((tt5+1.d0)**(1.d0+
     &         tt6)))
            endif
          endif
c
   10     continue
          rate=dmax1(0.d0,rate)
          rec(ion,atom)=rate
        enddo
c (recommode.eq.0) new std RR calc
      endif
c
      if (recommode.eq.1) then
c
c Original code, atom ion loops
c
c      write(*,*) 'recommode=1, old RR code'
c
        do atom=1,atypes
          rec(1,atom)=0.0d0
          maxio=maxion(atom)
          if (mapz(atom).eq.2) maxio=4
          do ion=2,maxio
            isos=mapz(atom)-ion+1
            rate=0.0d0
            ar=arad(ion,atom)
            xr=xrad(ion,atom)
            if (isos.eq.0) then
c
c also sets OTS rates for H I and He II
c
              rate=fkramer(ion,atom,tp)
              goto 20
            endif
            if ((mapz(atom).eq.2).and.(ion.eq.4)) then
c    He I total - n=1 rate for OTS
              ar=2.7d-13
              xr=0.787d0
            endif
            if (ar.lt.0.0d0) goto 20
            rate=arf(tp4,ar,xr)
   20       continue
            rate=dmax1(0.d0,rate)
            rec(ion,atom)=rate
          enddo
        enddo
c
c end if (recommode.eq.1) old calcs
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Add Dielectronic contributions independently
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (drrecmode.eq.0) then
c
c new std recombination calculations
c includes radiative and dielectronic recombination
c
c      write(*,*) 'drecmode=0, new DR code'
c
        do ir=1,nrec
          drate=0.d0
          atom=atrec(ir)
          ion=ionrec(ir)
          isos=mapz(atom)-ion+1
c
c these total rates already include DR
c
          if (recommode.eq.2) then
            if (noradid(ion,atom).ne.0) then
              goto 40
            endif
          endif
C           if (mapz(atom).eq.6) then
C             goto 40
C           endif
C           if (mapz(atom).eq.7) then
C             goto 40
C           endif
C           if (mapz(atom).eq.8) then
C             goto 40
C           endif
C         endif
c
c fe II-VI are always NORAD totals if avaiable,even recommode!=2
c
          if ((mapz(atom).eq.26).and.(ion.le.6)) then
            if (noradid(ion,atom).ne.0) then
              goto 40
            endif
          endif
c
          if (isos.eq.0) goto 40
          if ((mapz(atom).eq.2).and.(ion.eq.2).and.(recommode.eq.0))
     &     goto 40
c
          drtype=typedrec(ir)
          do i=1,8
            co(i)=adrec(i,ir)
            en(i)=bdrec(i,ir)
          enddo
          drate=0.d0
          if (drtype.gt.0) then
            if (drtype.eq.1) then
              do i=1,8
                if (en(i).le.0.d0) goto 30
                drate=drate+co(i)*dexp(-en(i)/tp)
              enddo
   30         continue
              drate=drate*t32
              goto 40
            endif
            if (drtype.eq.2) then
            drate=t32*(co(1)*dexp(-co(3)/tp))*(1.d0+co(2)*dexp(-co(4)/
     &         tp))
            endif
            if (drtype.eq.3) then
c
c compatibility case for special data file that reproduces old code
c
            drate=t32*(co(1)*dexp(-co(3)/tp))*(1.d0+co(2)*dexp(-co(4)/
     &         tp))
              if (tp.lt.1.d5) then
                a2=arrec(2,ir)
                ch=dble(max(1,ion-1))
                fcor=flowt(tp4,ch,a2)
                fcor=dmax1(0.d0,fcor)
                adlt=adltf(tp4,en(1),en(2),en(3),en(4),en(5))
                drate=drate+dmax1(0.0d0,adlt)
                rec(ion,atom)=rec(ion,atom)*fcor
              endif
            endif
          endif
c
   40     continue
c
          drate=dmax1(0.d0,drate)
          rec(ion,atom)=rec(ion,atom)+drate
c
        enddo
c
c end if (drrecmode.eq.0) new DR calcs
c
      endif
c
      if (drrecmode.eq.1) then
c
c Original DR code
c
c      write(*,*) 'drecmode=1, old DR code'
c
        do atom=1,atypes
          rec(1,atom)=0.0d0
          maxio=maxion(atom)
          if (mapz(atom).eq.2) maxio=4
          do ion=2,maxio
            isos=mapz(atom)-ion+1
            drate=0.0d0
            fcor=1.d0
            if (isos.eq.0) goto 50
            ar=arad(ion,atom)
            xr=xrad(ion,atom)
            if (ar.lt.0.0d0) goto 50
            ad=adi(ion,atom)
            bd=bdi(ion,atom)
            ta=t0(ion,atom)
            tb=t1(ion,atom)
            if ((tp.lt.1.d5).and.(cdilt(ion,atom).gt.0.d0)) then
              adlt=0.d0
              a=adilt(ion,atom)
              b=bdilt(ion,atom)
              c=cdilt(ion,atom)
              d=ddilt(ion,atom)
              f=fdilt(ion,atom)
              chh=ion-1
              if (chh.gt.maxion(atom)-1) chh=chh-(maxion(atom)-1)
              ch=dble(chh)
              fcor=flowt(tp4,ch,xr)
              fcor=dmax1(0.d0,fcor)
              rec(ion,atom)=rec(ion,atom)*fcor
              adlt=adltf(tp4,a,b,c,d,f)
              adlt=dmax1(0.d0,adlt)
              drate=adf(tp,ad,bd,ta,tb)+adlt
            else
              drate=adf(tp,ad,bd,ta,tb)
            endif
   50       continue
            drate=dmax1(0.d0,drate)
            rec(ion,atom)=rec(ion,atom)+drate
          enddo
        enddo
c
c end if (drrecmode.eq.1) old DR calcs
c
      endif
c
      return
c
      end
