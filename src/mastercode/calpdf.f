c
c     emidif, emidifcont is flux not photon flux, different 
c     from other subroutine, BE CAREFUL!
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine calpdf (t, de, dh, recpdf, rec0)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t, de, dh, telc
      real*8 xp(mxion, mxelem)
      real*8 energ,rkt,phots
      real*8 wid, ell
      real*8 recpdf(mxinfph),rectot
      real*8 normff(mxinfph),embri(mxinfph)
      real*8 maxrec,invmaxrec,rec0
      real*8 fbowen,dfb,fb,casab     
c
c      integer*4 ie,ii2,jn
      integer i,j,inl,nz
      integer line,series,atom,ion,trans
c
      integer*4 luso
c
      integer*4 ix,iy,iz
c
c
      maxrec=0.0d0
      rec0=0.0d0
c
      fb=0.d0
      if (zmap(2).gt.0) then
        fb=fbowen(t,1.0d6)
c        write(*,*) 'tot fb',fb, 1.d0-fb
      endif
      dfb=1.d0-fb
c
      if (zmap(1).ne.0) then
        line=3
        series=1
        caseab(1)=1.0d0        
c        caseab(1)=casab(1,line,series)
      endif
c
c     call casab with Helium Lyman gamma
c
      if (zmap(2).ne.0) then
        line=3
        series=1
        caseab(2)=1.0d0
c        caseab(2)=casab(2,line,series)
      endif
c
c     Assume heavy series are case A
c
      do atom=3,atypes
        line=3
        series=1
        caseab(atom)=casab(mapz(atom),line,series)
      enddo
c
c      
      call localem(t,de,dh)
c
      embri = 0.d0      
c
c
cc      open (luso,file='./calpdf.txt',status='unknown') 
ccc      open (luso,file='./embri.txt',status='unknown') 
c     
c        Local Resonance Lines
c
      do line=1,nxr3lines
        inl=xr3lines_bin(line)
        if (inl.gt.0) then
          if ((xr3lines_emilin(1,line).gt.epsilon)) then
            energ=xr3lines_egij(line)
            embri(inl)=embri(inl)+energ*xr3lines_emilin(1,line)
cc            write (luso,*) photev(inl), 0, 0, line, 'xr3lines'
          endif
        endif
      enddo
c
      do line=1,nxrllines
        inl=xrllines_bin(line)
        if (inl.gt.0) then
          if ((xrllines_emilin(1,line).gt.epsilon)) then
            energ=xrllines_egij(line)
            embri(inl)=embri(inl)+energ*xrllines_emilin(1,line)
cc            write (luso,*) photev(inl), 0, 0, line, 'xrllines'
          endif
        endif
      enddo
c
      do series=1,nhseries
        do line=1,nhlines
          inl=hbin(line,series)
          if (inl.gt.0) then
            if (hydlin(1,line,series).gt.epsilon) then
              energ=ev*lmev/hlambda(line,series)
              embri(inl)=embri(inl)+energ*hydlin(1,line,series)
cc            write (luso,*) photev(inl), 0, series, line, 'hseries'              
            endif
          endif
        enddo
      enddo
c
      do series=1,nheseries
        do line=1,nhelines
          inl=hebin(line,series)
          if (inl.gt.0) then
            if (hellin(1,line,series).gt.epsilon) then
              energ=ev*lmev/helambda(line,series)
              dfb=1.d0
              if ((line.eq.1).and.(series.eq.1)) then
                dfb=1.d0-fb
              endif
c
c degrade local HeII 303 photons locally to OIII BF lines by an
c approximate conversion fraction....  Set to 1-0 = 1.0, no
c conversion for now
c
              embri(inl)=embri(inl)+dfb*energ*hellin(1,line,series)
c
cc            write (luso,*) photev(inl), 0, series, line, 'heseries'              
c              
            endif
          endif
        enddo
      enddo
c
      do atom=3,atypes
        nz=mapz(atom)
        do series=1,nxhseries
          do line=1,nxhlines
            inl=xhbin(line,series,atom)
            if (inl.gt.0) then
              if (xhydlin(1,line,series,atom).gt.epsilon) then
c
                energ=ev*lmev/xhlambda(line,series,atom)
                embri(inl)=embri(inl)+energ*xhydlin(1,line,series,atom)
c
cc            write (luso,*) photev(inl), atom, series, line, 'xhseries'              
c
              endif
            endif
          enddo
        enddo
      enddo       
c
c
cc      do inl=1,infph-1
cc         energ=cphotev(inl)*ev
cc         write (luso,*) photev(inl), emidif(inl)*energ+embri(inl)         
cc      enddo
cc      close (luso)
cc      stop
c
c
      if (photev(1).ge.mcnuin) then
c        
         energ=cphotev(1)*ev
         wid=evplk*(photev(2)-photev(1))
c
         recpdf(1)=emidif(1)*energ*wid
         maxrec=recpdf(1)
c         
         if (embri(1).gt.epsilon) then
            recpdf(1)=recpdf(1)+embri(1)*wid
            maxrec=recpdf(1)
         endif         
c         
      else
c        
         energ=cphotev(1)*ev
         wid=evplk*(photev(2)-photev(1))
c                 
         rec0=rec0+emidif(1)*energ*wid
c
         if (embri(1).gt.epsilon) rec0=rec0+embri(1)*wid                     
c
      endif
c
      do inl=2,infph-1
c        
         if (photev(inl).lt.mcnuin) then
c
            energ=cphotev(inl)*ev
            wid=evplk*(photev(inl+1)-photev(inl))   
c
            rec0=rec0+emidif(inl)*energ*wid
c
            if (embri(inl).gt.epsilon) rec0=rec0+embri(inl)*wid                         
c
            goto 10
c
         endif
c
         if (photev(inl+1).ge.mcnuend) goto 20
c
         energ=cphotev(inl)*ev
         wid=evplk*(photev(inl+1)-photev(inl))         
c
         recpdf(inl)=recpdf(inl-1)+emidif(inl)*energ*wid
         if (embri(inl).gt.epsilon) 
     &      recpdf(inl)=recpdf(inl-1)+(emidif(inl)*energ+embri(inl))*wid       
c
         maxrec=recpdf(inl)
c
   10 continue
c
      enddo
c
   20 continue
c
c
      invmaxrec=1.d0/maxrec
c
c      ell=hloss+rloss+fslos+fmloss+xr3loss+xrlloss
c      ell=ell+fflos+colos
c      ell=ell+f3loss+feloss+gcool
c      ell=ell*ifpi
c
c      write (*,*) rec0, ell, maxrec      
c
c      rec0=rec0+ell

      rec0=rec0/(rec0+maxrec)
c
      do inl=1,infph
         recpdf(inl)=recpdf(inl)*invmaxrec
      enddo
c
      return
c
      end
