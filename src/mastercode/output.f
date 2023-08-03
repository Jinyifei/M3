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
c     writes out a machine readable ionisation balance.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine wbal (caller, pfx, np, p)
c
      include 'cblocks.inc'
c
      real*8 p(mxion,mxelem)
c
      character caller*4,tab*4
      character fn*32
      character pfx*16,sfx*4,fl*32
      integer*4 lups,i,np,nentries,j
c
      tab=char(9)
      lups=99
c
c     write out balance file
c
      fn=' '
      sfx='bln'
      call newfile (pfx, np, sfx, 3, fn)
      fl=fn(1:(np+3+5))
c
      open (lups,file=fl,status='NEW')
c
c
      write (lups,'("%")')
      write (lups,'("% Ionisation Balance File")')
      write (lups,'("%")')
      write (lups,'("%")')
c
      write (lups,10) runname
   10 format('% RUN:',a80)
      write (lups,'("%")')
c
      write (lups,*) 'Produced by ',caller,' :MAPPINGS V ',theversion
c
      nentries=0
c
      do i=1,atypes
        do j=1,maxion(i)
          nentries=nentries+1
        enddo
      enddo
c
      write (lups,*) nentries
c
      do i=1,atypes
        do j=1,maxion(i)
   30       format(i2,1x,i2,1x,1pg14.6)
          write (lups,30) mapz(i),j,p(j,i)
        enddo
      enddo
      close (lups)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c   Writes just the abundances
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wabund (luop)
c
      include 'cblocks.inc'
c
      integer*4 luop,i
c
      real*8 totmass, xi(mxelem)
c
  100 format(' Abundances:',/,
     &       ' El.,log[X/H],12+[X/H],log[Del],',
     &       'log[Dep],log[Gas],12+[Gas], Mass Fr.')
  110 format(  2x,a2,6(',',0pf8.3)',',1pe10.3)
c
  120 format(/' Base Mass XYZ:',/,
     &        ' Mass  X, ',0pf10.5,', Y,',0pf10.5,', Z,',0pf10.5)
c
  130 format( ' Gas Phase Mass XYZ:',/,
     &        ' Mass  X, ',0pf10.5,', Y,',0pf10.5,', Z,',0pf10.5)
c
      write(luop,100)
c
      totmass=0.d0
      do i=1,atypes
        totmass=totmass+zion0(i)*atwei(i)
      enddo
      do i=1,atypes
        xi(i)=zion0(i)*atwei(i)/totmass
      enddo

      do i=1,atypes
      write (luop,110) elem(i),dlog10(zion0(i)),12.d0+dlog10(zion0(i)),
     &          dlog10(deltazion(i)), dlog10(dion(i)),
     &          dlog10(zion(i)),12.d0+dlog10(zion(i)),xi(i)
      enddo
c
      write (luop,120) xi(1),xi(2),(1.d0-(xi(1)+xi(2)))
c
      totmass=0.d0
      do i=1,atypes
        totmass=totmass+zion(i)*atwei(i)
      enddo
      do i=1,atypes
        xi(i)=zion(i)*atwei(i)/totmass
      enddo

      write (luop,130) xi(1),xi(2),(1.d0-(xi(1)+xi(2)))
c
      write (luop,*) ' '
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    These are the standarised methods for writing out
c    element abundances and population data.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wionabal (luop, po)
c
c
      include 'cblocks.inc'
c
      integer*4 luop,ionmax,i,j
      real*8 po(mxion, mxelem)
      real*8 chcksum
c
      call wabund (luop)
c
      write (luop,10) (elem(i),i=1,atypes)
   10 format(/t8,30(4x,a2,4x))
c
c    find highest ionisation
c
      ionmax=0
      do i=1,atypes
        if (maxion(i).ge.ionmax) ionmax=maxion(i)
      enddo
      if (ionmax.gt.31) ionmax=31
c
c
  110 format(' ',a6,2x,30(1pg10.3))
      do j=1,ionmax
        chcksum=0.0d0
        do i=1,atypes
          chcksum=chcksum+po(j,i)
        enddo
        if (chcksum.gt.epsilon) then
          write (luop,110) rom(j),(po(j,i),i=1,atypes)
        endif
      enddo
c
c
      return
c
      end
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
c
c    This will become the standarised method for writing out
c    just population data.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wionpop (luop, po)
c
c
      include 'cblocks.inc'
c
      integer*4 luop,ionmax,i,j
      real*8 po(mxion, mxelem),chcksum
c
c    if luop > 0 then writing to a file
c
      write (luop,10) (elem(i),i=1,atypes)
   10 format(' ',t8,30(4x,a2,4x)/)
c
c    find highest ionisation
c
      ionmax=0
c
      do i=1,atypes
        if (maxion(i).ge.ionmax) ionmax=maxion(i)
      enddo
c
c
   20 format(' ',a6,2x,30(1pg10.3))
      do j=1,ionmax
        chcksum=0.0d0
        do i=1,atypes
          chcksum=chcksum+po(j,i)
        enddo
        if (chcksum.gt.0.0d0) then
          write (luop,20) rom(j),(po(j,i),i=1,atypes)
        endif
      enddo
c
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     wmodel:write out current plasma model conditions.
c
c     wmod = 'SCRN'   for vertical format, luop is not used.
c     wmod = 'FILE'   for side by side format, luop is
c                     assumed to be open.
c     wmod = 'PROP'   properties only, inline,luop
c                     assumed to be open.
c     wmod = 'LOSS'   losses and gains only in line, luop
c                     assumed to be open.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wmodel (luop, t, de, dh, dstep, wmod)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh,dstep,rl,fsl,sint,aspd
      real*8 q2los,cspd,ram,mag
      real*8 trec,tcoll,press,en,ue,cts,rhotot,wmol
      real*8 a(50)
      integer*4 luop,i,atom
      character wmod*4
c
c           Functions
c
      real*8 fcietim,fpressu,frectim,frho
c
      press=fpressu(t,dh,pop)
      rhotot=frho(de,dh)
c
      trec=frectim(t,de,dh)
      tcoll=fcietim(t,de,dh)
c
      ram=rhotot*vel0*vel0
c
      q2los=0.d0
      do atom=1,atypes
      q2los=q2los+h2qbri(atom)
      enddo
c
      en=zen*dh
c
      ue=1.5d0*(en+de)*rkb*t
      cts=ue/tloss
c
      cspd=dsqrt((5.d0/3.d0)*press/rhotot)
c
      mag=bm0
      aspd=dsqrt(mag*mag/(fpi*rhotot))
c
      wmol=rhotot/(en+de)
c
      rl=rloss+xr3loss+xrlloss
      fsl=fslos+fmloss
c
      if (wmod.eq.'LOSS') then
c
   10 format(27(1pe13.6,x))
c
        do i=1,27
          a(i)=0.d0
        enddo
        if (dabs(dstep).gt.ioepsilon) a(1)=dstep
        if (dabs(t).gt.ioepsilon) a(2)=t
        if (dabs(de).gt.ioepsilon) a(3)=de
        if (dabs(dh).gt.ioepsilon) a(4)=dh
        if (dabs(en).gt.ioepsilon) a(5)=en
        if (dabs(dlos).gt.ioepsilon) a(6)=dlos
        if (dabs(eloss).gt.ioepsilon) a(7)=eloss
        if (dabs(egain).gt.ioepsilon) a(8)=egain
        if (dabs(hloss).gt.ioepsilon) a(9)=hloss
        if (dabs(chgain).gt.ioepsilon) a(10)=chgain
        if (dabs(xr3loss+xrlloss).gt.ioepsilon) a(11)=xr3loss+xrlloss
        if (dabs(rloss).gt.ioepsilon) a(12)=rloss
        if (dabs(fsl).gt.ioepsilon) a(13)=fsl
C       if (dabs(xrloss).gt.ioepsilon) a(14)=xrloss
        if (dabs(fmloss).gt.ioepsilon) a(15)=fmloss
        if (dabs(feloss).gt.ioepsilon) a(16)=feloss
        if (dabs(q2los).gt.ioepsilon) a(17)=q2los
        if (dabs(cmplos).gt.ioepsilon) a(18)=cmplos
        if (dabs(fflos).gt.ioepsilon) a(19)=fflos
        if (dabs(colos).gt.ioepsilon) a(20)=colos
        if (dabs(pgain).gt.ioepsilon) a(21)=pgain
        if (dabs(rngain).gt.ioepsilon) a(22)=rngain
        if (dabs(cosgain).gt.ioepsilon) a(23)=cosgain
        if (dabs(gheat-gcool).gt.ioepsilon) a(24)=gheat-gcool
        if (dabs(gheat).gt.ioepsilon) a(25)=gheat
        if (dabs(gcool).gt.ioepsilon) a(26)=gcool
        if (dabs(paheat).gt.ioepsilon) a(27)=paheat
c
        write (luop,10) (a(i),i=1,27)
c
      endif
      if (wmod.eq.'PROP') then
c
   20 format(15(1pe13.6,x))
c
        do i=1,15
          a(i)=0.d0
        enddo
        if (dabs(dstep).gt.ioepsilon) a(1)=dstep
        if (dabs(t).gt.ioepsilon) a(2)=t
        if (dabs(de).gt.ioepsilon) a(3)=de
        if (dabs(dh).gt.ioepsilon) a(4)=dh
        if (dabs(en).gt.ioepsilon) a(5)=en
        if (dabs(pop(1,1)).gt.ioepsilon) a(6)=pop(1,1)
        if (dabs(rhotot).gt.ioepsilon) a(7)=rhotot
        if (dabs(press).gt.ioepsilon) a(8)=press
        if (dabs(vel0).gt.ioepsilon) a(9)=vel0
        if (dabs(ram).gt.ioepsilon) a(10)=ram
        if (dabs(cspd).gt.ioepsilon) a(11)=cspd
        if (dabs(dstep).gt.ioepsilon) a(12)=dstep
        if (dabs(aspd).gt.ioepsilon) a(13)=aspd
        if (dabs(bm0).gt.ioepsilon) a(14)=bm0
        if (dabs((avgpot(1)+avgpot(2))*0.5d0).gt.ioepsilon) a(15)=
     &   (avgpot(1)+avgpot(2))*0.5d0
        write (luop,20) (a(i),i=1,15)
      endif
c
      if (wmod.eq.'FILE') then
        do i=1,9
          a(i)=0.d0
        enddo
        if (dabs(tloss).gt.ioepsilon) a(1)=tloss
        if (dabs(eloss).gt.ioepsilon) a(2)=eloss
        if (dabs(tloss).gt.ioepsilon) a(3)=tloss
        if (dabs(egain).gt.ioepsilon) a(4)=egain
        if (dabs(dlos).gt.ioepsilon) a(5)=dlos
        if (dabs(hloss).gt.ioepsilon) a(6)=hloss
        if (dabs(chgain).gt.ioepsilon) a(7)=chgain
        if (dabs(xrlloss+xr3loss).gt.ioepsilon) a(8)=xrlloss+xr3loss
        if (dabs(rloss).gt.ioepsilon) a(9)=rloss
        write (luop,30) (a(i),i=1,9)
c
   30 format(/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS  :',1pe11.3,' : EFF. LOSS   :',1pe11.3,' : ',
     &': TOTAL LOSS:', 1pe12.5,' :::::::::::::::::::::::::::::'/
     &': EFF. GAIN   :',1pe11.3,' : FRAC. RESID.:',1pe11.3,' : ',
     &': COLEXC H. :', 1pe12.5,' : CHARGE EX.  :', 1pe12.5,' :'/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &': CASC. RES.:', 1pe12.5,' : OLD RESON.  :', 1pe12.5,' :')
c
        do i=1,27
          a(i)=0.d0
        enddo
        if (dabs(t).gt.ioepsilon) a(1)=t
        if (dabs(rhotot).gt.ioepsilon) a(2)=rhotot
        if (dabs(fsl).gt.ioepsilon) a(3)=fsl
        if (dabs(f3loss).gt.ioepsilon) a(4)=f3loss
        if (dabs(en).gt.ioepsilon) a(5)=en
        if (dabs(de).gt.ioepsilon) a(6)=de
        if (dabs(fmloss).gt.ioepsilon) a(7)=fmloss
        if (dabs(feloss).gt.ioepsilon) a(8)=feloss
        if (dabs(pop(1,1)).gt.ioepsilon) a(9)=pop(1,1)
        if (dabs(press).gt.ioepsilon) a(10)=press
        if (dabs(q2los).gt.ioepsilon) a(11)=q2los
        if (dabs(cmplos).gt.ioepsilon) a(12)=cmplos
        if (dabs(vel0).gt.ioepsilon) a(13)=vel0
        if (dabs(ram).gt.ioepsilon) a(14)=ram
        if (dabs(fflos).gt.ioepsilon) a(15)=fflos
        if (dabs(colos).gt.ioepsilon) a(16)=colos
        if (dabs(cspd).gt.ioepsilon) a(17)=cspd
        if (dabs(dstep).gt.ioepsilon) a(18)=dstep
        if (dabs(pgain).gt.ioepsilon) a(19)=pgain
        if (dabs(rngain).gt.ioepsilon) a(20)=rngain
        if (dabs(aspd).gt.ioepsilon) a(21)=aspd
        if (dabs(mag).gt.ioepsilon) a(22)=mag
        if (dabs(cosgain).gt.ioepsilon) a(23)=cosgain
        if (dabs(paheat).gt.ioepsilon) a(24)=paheat
        if (dabs(avgpot(1)).gt.ioepsilon) a(25)=avgpot(1)
        if (dabs(avgpot(2)).gt.ioepsilon) a(26)=avgpot(2)
        if (dabs(gheat-gcool).gt.ioepsilon) a(27)=(gheat-gcool)
        write (luop,40) (a(i),i=1,27)
c
   40 format(
     &': TEMP.       :',1pe11.4,' : DENSITY     :',1pe11.3,' : ',
     &': INTER.    :', 1pe12.5,' : 3 LEV. FINE :', 1pe12.5,' :'/
     &': N IONS      :',1pe11.3,' : N ELECTRONS :',1pe11.3,' : ',
     &': MULTI-LVL.:', 1pe12.5,' : MULTI_IRON  :', 1pe12.5,' :'/
     &': FRAC. N. H  :',1pe11.3,' : PRESSURE    :',1pe11.3,' : ',
     &': 2PHOTON   :', 1pe12.5,' : COMPTON     :', 1pe12.5,' :'/
     &': FLOW VELOC. :',1pe11.3,' : RAM PRESS.  :',1pe11.3,' : ',
     &': FREFRE    :', 1pe12.5,' : COLION      :', 1pe12.5,' :'/
     &': SOUND SPEED :',1pe11.3,' : SLAB DEPTH  :',1pe11.3,' : ',
     &': PGAIN     :', 1pe12.5,' : RNGAIN      :', 1pe12.5,' :'/
     &': ALFEN SPEED :',1pe11.3,' : MAG. FIELD  :',1pe11.3,' : ',
     &': COSMIC    :', 1pe12.5,' : PHOTO.  PAH :', 1pe12.5,' :'/
     &': GRA GRN POT :',1pe11.3,' : SIL GRN POT :',1pe11.3,' : ',
     &': PHOT. GRN :', 1pe12.5,' :'/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
c
      endif
c
      if (wmod.eq.'SCRN') then
c
        write (*,50) tloss,eloss,egain,dlos,t,rhotot,en,de,pop(1,1),
     &   press,vel0,ram,cspd,dstep
        write (*,60) aspd,mag,avgpot(1),avgpot(2),sint,wdil,teff,alnth,
     &   turn,cut
   50 format(/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS  :',1pe11.3,' : EFF. LOSS   :',1pe11.3,' :'/
     &': EFF. GAIN   :',1pe11.3,' : FRAC. RESID.:',1pe11.3,' :'/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TEMP.       :',1pe11.4,' : DENSITY     :',1pe11.3,' :'/
     &': N IONS      :',1pe11.3,' : N ELECTRONS :',1pe11.3,' :'/
     &': FRAC. N. H  :',1pe11.3,' : PRESSURE    :',1pe11.3,' :'/
     &': FLOW VELOC. :',1pe11.3,' : RAM PRESS.  :',1pe11.3,' :'/
     &': SOUND SPEED :',1pe11.3,' : SLAB DEPTH  :',1pe11.3,' :')
   60 format(
     &': ALFEN SPEED :',1pe11.3,' : MAG. FIELD  :',1pe11.3,' :'/
     &': GRA GRN POT.:',1pe11.3,' : SIL GRN POT.:',1pe11.3,' :'//
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': SOURCE INT. :',1pe11.3,' : DILUT. F.   :',1pe11.3,' :'/
     &': STAR TEMP.  :',1pe11.3,' : ALPHA       :',1pe11.3,' :'/
     &': TURN-ON     :',1pe11.3,' : CUT-OFF     :',1pe11.3,' :')
c
        write (*,70) tloss,hloss,chgain,xr3loss,xrlloss,rloss,
     &   f3loss,fmloss,feloss,q2los,cmplos,fflos,colos,
     &   pgain,rngain,cosgain,paheat,(gheat-gcool)
   70 format(
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS:', 1pe12.5,' :::::::::::::::::::::::::::::'/
     &': COLEXC H. :', 1pe12.5,' : CHARGE EX.  :', 1pe12.5,' :'/
     &': CASC. RES.:', 1pe12.5,' : X. RESON.   :', 1pe12.5,' :'/
     &': OLD RES. .:', 1pe12.5,' : 3 LEV. FINE :', 1pe12.5,' :'/
     &': MULIT-LVL.:', 1pe12.5,' : MULTI-IRON  :', 1pe12.5,' :'/
     &': 2PHOTON   :', 1pe12.5,' : COMPTON     :', 1pe12.5,' :'/
     &': FREFRE    :', 1pe12.5,' : COLION      :', 1pe12.5,' :'/
     &': PGAIN     :', 1pe12.5,' : RNGAIN      :', 1pe12.5,' :'/
     &': COSMIC    :', 1pe12.5,' : PHOTO. PAH  :', 1pe12.5,' :'/
     &': DUST  H-C :', 1pe12.5,/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
c
      endif
c
      return
c
      end
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
c    Output Vector Units: Lambda(E)  = ergs cm3 /s/eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wemiss2 (caller, pfx, np, t, de, dh, dr, scale, tp)
c
      include 'cblocks.inc'
c
      real*8 tp(mxinfph),tl(mxinfph),t,de,dh,dr
      real*8 wid,scale
c
      integer*4 lups,i,j,np
c
      character caller*4,tab*4
      character fn*32
      character pfx*16,sfx*4,fps*32
c
      tab=char(9)
      lups=44
c
      fps=' '
c
c     output upstream field photon source file
c
      fn=' '
      sfx='emi'
      call newfile (pfx, np, sfx, 3, fn)
      fps=fn(1:(np+3+5))
c
      open (lups,file=fps,status='NEW')
      i=infph-1
c
      write (lups,10) scale
   10 format('% EMISSIVITY SPECTRUM ',/,
     &       '% UNITS: eV (bin centre) vs',/,
     &       '%      : 1e23/ne^2 * L(bin)',/,
     &       '%      : ',1pe11.4,' * L(bin) (ergs cm3/s)')
c
   20 format('%:::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '% TOTAL LOSS:',1pe12.5,/,'% EFF. LOSS  :',1pe12.5,/,
     &     '% EFF.  GAIN:',1pe12.5,/,'% FRAC. RESID:',1pe12.5,/,
     &     '% TEMP.     :',1pe12.5,/,'% EL. DENSITY:',1pe12.5,/,
     &     '% H. DENSITY:',1pe12.5)
c
      write (lups,20) tloss,eloss,egain,dlos,t,de,dh
c
      write (lups,*) 'Produced by ',caller,' :MAPPINGS V ',theversion
   30 format(' Temperature(K) =',1pe12.5,' Log(Ne)=',1pe12.5)
      write (lups,30) t,dlog10(de)
c
      do j=1,infph-1
c       wid=photev(j+1)-photev(j)
        if (scale.gt.0.d0) then
        wid=1.d0*evplk
        tl(j)=(scale*fpi*tp(j)*wid/dr)
        if (tl(j).lt.ioepsilon) then
          tl(j)=0.d0
        endif
        else
          tl(j)=0.d0
        endif
c   40   format(1pe14.7,x,1pe14.7)
c        write (lups,40) cphotev(j),tl(j)
      enddo
c
c 10 col output for genx
c
  50  format(10(x,1pe11.4))
      write (lups,50) (tl(j), j=1,(infph-1))
c
      close (lups)
      return
c
      end
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
c     standard routine to write out a source vector, given the
c     field in tphot.
c
c     Vector Units: Fnu  = ergs/s/cm2/Hz/sr
c
c     External files are all 3D flux units 1/(4pi) because they are
c     generally written from diffuse field vectors.  2D (1/pi) units
c     are only used internally for source vectors and are converted
c     when files are read in in photsou.f
c
c     wmod = 'REAL' use total emission from slab, don't use dr
c     wmod = 'NORM' then normalise to 1cm slab, use dr
c     wmod = 'NFNU' then write nu v nuFnu (ergs/s/cm2/sr)
c     wmod = 'LFLM' then write lam (A) v Flam (ergs/s/cm2/A)
c
c     Now puts out two columns in the files: energy in eV and flux.
c
c     same as wpsou but uses given file name and appends to file if old
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wpsoufile (caller, fname, wmod, t, de, dh, vl, dv, vl0,
     & di, dr, eltim, tst, scale, tp)
c
      include 'cblocks.inc'
c
      real*8 tp(mxinfph),tl(mxinfph),t,de,dh,dr
      real*8 bv,di,eltim,tst,vl,dv,vl0,scale
      real*8 blum, ilum, wid, clam, lambda
c
      integer*4 lups,i,j
c      integer*4 np
c
      character caller*4,wmod*4,tab*4
      character fname*32
c
      logical iexi
c functions
      real*8 fnair
c
      tab=char(9)
      lups=44
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     total energy in Inu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      blum=0.d0
      ilum=0.d0
c
c     cvt to 1/4pi units in tl and sum
c
      do i=1,infph-1
        tl(i)=0.d0
        if (tp(i).ge.ioepsilon) then
          tl(i)=scale*tp(i)
        endif
c
        wid=photev(i+1)-photev(i)
        blum=blum+tl(i)*wid*evplk
        if (photev(i).ge.iph) ilum=ilum+tl(i)*wid*evplk
      enddo
c
      blum=fpi*blum
      ilum=fpi*ilum
c
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    Assuming plane parallel geometry and dilution of 0.5 and intvec
c    assumes 1/4pi units in vector.
c
      qhi=0.d0
      qhei=0.d0
      qheii=0.d0
      qht=0.d0
c
      call intvec (tl, qhi, qhei, qheii, qht)
c
      i=infph-1
c
      inquire (file=fname,exist=iexi)
      if (.not.iexi) then
        open (lups,file=fname,status='NEW')
        if (wmod.eq.'LFLM') then
          write (lups,5) runname,wmod,t,de,dh
    5    format('% LFLAM SPECTRUM '/
     &        '% RUN:',a80/
     &        '% MODE:',a4/
     &        '% UNITS: Lambda (Air) vs Flambda (ergs/s/cm2/A)'/
     &        '% UNITS: bin centre wavelengths'/
     &        '% TEMP.    :',1pe11.4/
     &        '% El. DENS.:',1pe11.4/
     &        '% H.  DENS.:',1pe11.4)
        else if (wmod.eq.'NFNU') then
          write (lups,10) runname,wmod,t,de,dh
   10    format('% NUFNU SPECTRUM '/
     &        '% RUN:',a80/
     &        '% MODE:',a4/
     &        '% UNITS: nu (Hz) v nuFnu (Hz v ergs/s/cm2/sr)'/
     &        '% UNITS: bin lower/left edge frequencies'/
     &        '% TEMP.    :',1pe11.4/
     &        '% El. DENS.:',1pe11.4/
     &        '% H.  DENS.:',1pe11.4)
        else if (wmod.eq.'NORM') then
          write (lups,15) runname,wmod,t,de,dh
   15    format('% COOLING SPECTRUM '/
     &        '% RUN:',a80/
     &        '% MODE:',a4/
     &        '% UNITS: E (eV) vs L(E) (ergs/s cm3)'/
     &        '% (4pi nu Fnu/dr)/(n^2)'/
     &        '% UNITS: bin lower/left edge energies'/
     &        '% TEMP.    :',1pe11.4/
     &        '% El. DENS.:',1pe11.4/
     &        '% H.  DENS.:',1pe11.4)
        else if (wmod.eq.'XRAY') then
          write (lups,20) runname,wmod,t,de,dh
   20    format('% CMFGEN XRAY SPECTRUM '/
     &        '% RUN:',a80/
     &        '% MODE:',a4/
     &        '% UNITS: E (eV) vs F_E (ergs/s/cm2/ev)'/
     &        '% UNITS: bin centre energies'/
     &        '% TEMP.    :',1pe11.4/
     &        '% El. DENS.:',1pe11.4/
     &        '% H.  DENS.:',1pe11.4)
        else
          write (lups,25)runname,wmod
   25    format('%PHOTON SOURCE FILE '/
     &        '% DIFFUSE FIELD PLUS SOURCE'/
     &        '% RUN:',a80/
     &        '% MODE:',a4/
     &        '% UNITS:eV vs Fnu (ergs/s/cm2/Hz/sr)'/
     &        '% Fnu in 3D units (1/(4pi)) sr not '/
     &        '% 2D source (1/pi) units. '/
     &        '% bin energies are lower/left edge of bins.'/
     &        '% fluxes are average over bin.')
        endif
      else
        open (lups,file=fname,status='OLD',access='APPEND')
      endif
      if (wmod.eq.'REAL') then
   30 format('%:::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '% TOTAL LOSS  :',1pe11.3,'% EFF. LOSS   :',1pe11.3/
     &     '% EFF.  GAIN  :',1pe11.3,'% FRAC. RESID.:',1pe11.3/
     &     '% TEMP.       :',1pe11.4,'% H DENSITY   :',1pe11.3/
     &     '% ELECTR. DENS:',1pe11.3,'% FR. NEUT. H :',1pe11.3/
     &     '% SHOCK VELOC.:',1pe11.3,'% FLOW VELOC. :',1pe11.3/
     &     '% DELTA VELOC :',1pe11.3,'% DISTANCE.   :',1pe11.3/
     &     '% SLAB DEPTH  :',1pe11.3,'% ELAPSED TIME:',1pe11.3/
     &     '% TIME STEP   :',1pe11.3,'% SOURCE MOD. :',a11    /
     &     '% DILUT. F.   :',1pe11.3,'% SRC. TEMP.  :',1pe11.3/
     &     '% ALPHA       :',1pe11.3,'% CUT-OFF     :',1pe11.3)
c
      write (lups,30) tloss,eloss,egain,dlos,t,dh,de,pop(1,1),vl0,vl,dv,
     &di,dr,eltim,tst,iso,wdil,teff,alnth,cut
c
   40 format('%:::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '%: Total intensity   : ',1pe12.5,' (ergs/s/cm^2)      :'/
     &     '%: Ion.  intensity   : ',1pe12.5,' (ergs/s/cm^2)      :'/
     &     '%: FQHI  (1-1.8Ryd)  : ',1pe12.5,' (phots/cm^2/s)     :'/
     &     '%: FQHeI (1.8-4Ryd)  : ',1pe12.5,' (phots/cm^2/s)     :'/
     &     '%: FQHeII   (>4Ryd)  : ',1pe12.5,' (phots/cm^2/s)     :'/
     &     '%: FQ tot   (>1Ryd)  : ',1pe12.5,' (phots/cm^2/s)     :'/
     &     '%:::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
      write (lups,40) blum,ilum,qhi,qhei,qheii,qht
      endif
c
      write (lups,*) 'Produced by ',caller,' :MAPPINGS V ',theversion
      if (wmod.eq.'REAL') then
      write (lups,*) fieldversion
      write (lups,*) infph-1
      endif
c
      do j=1,infph-1
c
        bv=cphotev(j)*evplk  ! nu
        clam=1.0d8*cls/(bv*evplk)
        tl(j)=0.d0
c
        if ((tp(j).ge.ioepsilon).and.(dr.gt.0.d0)) then
          if (wmod.eq.'NORM') tl(j)=fpi*(scale*bv*tp(j)/dr)
          if (wmod.eq.'XRAY') tl(j)=fpi*(scale*evplk*tp(j)/dr)
          if (wmod.eq.'REAL') tl(j)=(scale*tp(j))
          if (wmod.eq.'NFNU') tl(j)=(scale*bv*tp(j))
          if (wmod.eq.'LFLM') tl(j)=fpi*(scale*bv*tp(j)/clam)
        endif
c
        if (tl(j).lt.ioepsilon) then
          tl(j)=0.d0
        endif
c
   50    format(1pe14.7,' ',1pe14.7)
c
        if (wmod.eq.'LFLM') then
          lambda = 1.0d8*cls/(cphotev(j)*evplk)
          if ( lambda.ge. 3300.d0) then
          if ( lambda.le.10000.d0) then
          if (tl(j).gt.0.d0) then
            lambda = lambda/fnair(lambda) ! cvt to std air wavelengths
            write (lups,50) lambda,tl(j)
          endif
          endif
          endif
        else if (wmod.eq.'NFNU') then
          if (tl(j).gt.ioepsilon) then
            write (lups,50) photev(j)*evplk,tl(j)
          endif
        else if (wmod.eq.'NORM') then
            write (lups,50) photev(j),tl(j)
        else if (wmod.eq.'XRAY') then
          if (tl(j).gt.0.d0) then
            write (lups,50) cphotev(j),tl(j)
          endif
        else
            write (lups,50) photev(j),tl(j)
        endif
c
      enddo
c
      if (wmod.eq.'REAL') write (lups,50) photev(infph),0.d0
      if (wmod.eq.'NORM') write (lups,50) photev(infph),0.d0
      if (wmod.eq.'NFNU') write (lups,50) photev(infph),0.d0
c
      close (lups)
c
      return
c
      end
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
c     standard routine to write out a source vector, given the
c     field in tphot.
c
c     Vector Units: Fnu  = ergs/s/cm2/Hz/sr
c
c     External files are all 3D flux units 1/(4pi) because they are
c     generally written from diffuse field vectors.  2D (1/pi) units
c     are only used internally for source vectors and are converted
c     when files are read in in photsou.f
c
c     wmod = 'REAL' use total emission from slab, don't use dr, .sou
c     wmod = 'NORM' then normalise to 1cm slab, use dr, .fnu
c     wmod = 'NFNU' then write nu v nuFnu (ergs/s/cm2/sr), .nfn
c     wmod = 'LFLM' then write lam (A) v Flam (ergs/s/cm2/A) .lam
c
c     Now puts out two columns in the files: energy in eV and flux.
c
c     old form with just the prefix, calls wpsoufile above
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wpsou (caller, pfx, np, wmod, t, de, dh, vl, dv, vl0,
     &di, dr, eltim, tst, scale, tp)
c
      include 'cblocks.inc'
c
      real*8 tp(mxinfph),t,de,dh,dr
      real*8 di,eltim,tst,vl,dv,vl0,scale
c
      integer*4 np
c
      character caller*4,wmod*4
      character fn*32
      character pfx*16,sfx*4,fps*32
c
      fps=' '
      if (wmod.eq.'LFLM') then
c
        fn=' '
        sfx='lam'
        call newfile (pfx, np, sfx, 3, fn)
        fps=fn(1:(np+3+5))
c
      else if (wmod.eq.'NFNU') then
c
        fn=' '
        sfx='nfn'
        call newfile (pfx, np, sfx, 3, fn)
        fps=fn(1:(np+3+5))
c
      else if (wmod.eq.'NORM') then
c
        fn=' '
        sfx='emi'
        call newfile (pfx, np, sfx, 3, fn)
        fps=fn(1:(np+3+5))
c
      else if (wmod.eq.'XRAY') then
c
        fn=' '
        sfx='dat'
        call newfile (pfx, np, sfx, 3, fn)
        fps=fn(1:(np+3+5))
c
      else
c
c     output upstream field photon source file
c
        fn=' '
        sfx='sou'
        call newfile (pfx, np, sfx, 3, fn)
        fps=fn(1:(np+3+5))
c
      endif

      call wpsoufile (caller, fps, wmod, t, de, dh, vl, dv, vl0,
     & di, dr, eltim, tst, scale, tp)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrspec(caller, pfx, np, t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t,de,dh
      integer*4 np
      character caller*4
      character pfx*16

c
      integer*4 lu
c
      character fsp*32
      character fn*32, sfx*4
      character linemod*4, spmod*4
c
      lu=99
c
      fn=' '
      sfx='csv'
      call newfile (pfx, np, sfx, 3, fn)
      fsp=fn(1:(np+3+5))
      open (lu,file=fsp,status='NEW')
   10 format(
     &'::::::::::::::::::::::::::::::::::::'
     &,'::::::::::::::::::::::::::::::::::::',/
     & ' Spectrum File, MAPPINGS V ',a8,/
     &'::::::::::::::::::::::::::::::::::::'
     &,'::::::::::::::::::::::::::::::::::::',/
     &,t5,' Run   : ,',a)
      write (lu,10) theversion,runname
   20  format('%::::::::::::::::::::::::::'/
     &'% TEMP.       :,' ,1pe11.4/'% H DENSITY   :,' ,1pe11.4/
     &'% ELECTR. DENS:,' ,1pe11.4/'% ION DENSITY :,' ,1pe11.4/
     &'%::::::::::::::::::::::::::')
      write (lu,20) t,dh,de,zen
      write (lu,*) 'Produced by ',caller
      spmod='REL'
      linemod='LAMB'
      call spec2 (lu, linemod, spmod)
      close (lu)
      return
      end
c
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
c*******TO OUTPUT ON MEDIUM*LUOP  THE EMISSION SPECTRUM AND
c    THE AVERAGE TEMPERATURES AND IONIC POPULATIONS
c    CALL SUBROUTINE SPECTRUM
c
c
c     modified for varable numbers of elemnts and ionisation stages
c
c     RSS 8/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wrsppop (luop)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 achh,dhas,dhtt
      integer*4 luop
      character mode*4
c
c     FOR population total calculation
c      real*8 popintot(mxion, mxelem)
c      integer*4 i,j
c
c           Functions
c
      real*8 densnum,favcha
c
      achh=favcha(pam,1)
      dhtt=densnum(1.d0)
      dhas=deav/((dhtt*achh)+1.d-38)
      write (luop,10) zgas
c
c
   10 format(/' ',t9,'Abundances of the elements relative',
     &' to Hydrogen and their average state of ionisation :','  (Zgas='
     &,f7.4,' Zsun)'/)
      call wionabal (luop, pam)
c
c
      write (luop,20)
   20 format(/' ',t9,'Average distances for a given state ',
     &'of ionisation :'/)
c
      call wionpop (luop, rdisa)
c
c
      write (luop,30)
   30 format(/' ',t9,'Average ionic temperatures :'/)
c
      call wionpop (luop, teav)
c
      write (luop,40)
   40 format(/' ',t9,'Integrated column densities:'/)
c
      call wionpop (luop, popint)
c
c
      write (luop,50) deav,dhas,achh,dhtt
   50 format(/,' ',t9,'Average ionic electron dens. :',2x,'<ne>='
     & ,1pe10.3,3x,'<nH>=',1pe10.3,3x,'<CHARGE>=',1pe10.3,3x,'num/H :'
     & ,1pe11.4,/)
c
      call wionpop (luop, deam)
c
      write (luop,60)
   60 format(/' ',t9,
     &        'Legacy Fit Line Flux Weighted averaged quantities :')
      write (luop,70)
   70 format(' ',t11,
     &       '<T>NII',t21,'<ne>NII',t31,'<T>OIII',t41,'<ne>OIII')
      write (luop,80) tnii,denii,toiii,deoiii
   80 format(' ',t8,2(0pf10.0,1pe10.3)/)
c
      write (luop,90)
   90 format(' ',t9,'Photon field averaged quantities :'/)
      write (luop,100)
  100 format(' ',t11,'<ZETAE>',t21,'<QHDH>'/)
      write (luop,110) zetaeav,qhdha
  110 format(' ',t8,2(1pe10.3)/)
c
c    ***OUTPUT EMISSION SPECTRUM RELATIVE TO H-BETA
c
      mode='REL'
c
      call spectrum (luop, mode)
c
      return
      end
