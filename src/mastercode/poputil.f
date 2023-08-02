cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c 
c   copyright 1979-2012
c 
c       Version v5.1.03
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      subroutine showpop (popin)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem)
      integer*4 at, ion
      do at=1,atypes
        do ion=1,maxion(at)
          if (popin(ion,at).gt.0.d0) then
            write (*,'(x,a2,x,i2.2,x,i2.2,x,1pg11.4)') elem(at),mapz(at)
     &       ,ion,popin(ion,at)
          endif
        enddo
      enddo
      return
      end
      subroutine copypop (popin, popout)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem), popout(mxion, mxelem)
      integer*4 at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popout(ion,at)=popin(ion,at)
      return
      end
      subroutine addpop (popin, popout)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem), popout(mxion, mxelem)
      integer*4 at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popout(ion,at)=popin(ion,at)+popout(ion,at)
      return
      end
      subroutine clearpop (popin)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem)
      integer*4 at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popin(ion,at)=0.0d0
      return
      end
      subroutine scalepop (popin, x)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem),x
      integer*4 at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popin(ion,at)=popin(ion,at)*x
      return
      end
      subroutine copysteppop (step, popin, popout)
      include 'cblocks.inc'
      real*8 popin(mxifsteps, mxion, mxelem), popout(mxion, mxelem)
      integer*4 step,at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popout(ion,at)=popin(step,ion,at)
      return
      end
      subroutine copypopstep (popin, step, popout)
      include 'cblocks.inc'
      real*8 popin(mxion, mxelem), popout(mxifsteps, mxion, mxelem)
      integer*4 step,at, ion
      do 10 at=1,atypes
        do 10 ion=1,mxion
   10   popout(step,ion,at)=popin(ion,at)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c 
c   copyright 1979-2012
c 
c       Version v5.1.03
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO OUTPUT AVERAGE OF POPW & POPCO IN POPOUT
c     THE WEIGHT OFF POPW IS : WEI (0<=WEI<=1)
c 
      subroutine averinto (wei, popw, popco, popout)
c 
      include 'cblocks.inc'
c 
      real*8 popw(mxion, mxelem), popco(mxion, mxelem)
     &                ,popout(mxion, mxelem)
      real*8 wei, weiw, weico
c 
      integer*4 i, j
c 
      weiw=dmax1(0.0d0,dmin1(1.0d0,wei))
      weico=dmax1(0.d0,1.d0-weiw)
c 
      do 10 i=1,atypes
        do 10 j=1,maxion(i)
        popout(j,i)=(weiw*popw(j,i))+(weico*popco(j,i))
   10 continue
c 
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c 
c   copyright 1979-2012
c 
c       Version v5.1.03
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c*******TO DERIVE RELATIVE CHANGE IN IONIC POPULATIONS : DIF
c     LIM : NUMBER OF LAST ELEMENTS TO BE INCLUDED
c     TRE : TRESHOLD FOR MEASURING CHANGE
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      subroutine difpop (popin, popfi, tre, lim, dif)
c 
      include 'cblocks.inc'
c 
c 
      real*8 popin(mxion, mxelem), popfi(mxion, mxelem)
      real*8 tre, dif, thresh, maxdif, tmdif,sum,p1,p2
      integer*4 num, iel, ion, lim, matom, mion
c 
      dif=0.0d0
      sum=0.d0
      num=0
      thresh=tre
c 
      matom=0
      mion=0
      maxdif=0.d0
c 
      if (thresh.le.epsilon) thresh=epsilon
      do iel=1,lim
        do ion=1,maxion(iel)
          if ((popin(ion,iel).ge.thresh).and.(popfi(ion,iel).ge.thresh))
     &      then
            p1=dlog10(popin(ion,iel))
            p2=dlog10(popfi(ion,iel))
            tmdif=dabs(p1-p2)
            if (tmdif.ge.maxdif) then
              maxdif=tmdif
              matom=iel
              mion=ion
            endif
            sum=sum+tmdif
            num=num+1
          endif
        enddo
      enddo
c 
      dif=sum/(num+epsilon)
c 
c     if (expertmode.gt.0) then
c       write (*,*) 'Max difference: ',elem(matom),rom(mion),'=',maxdif
c     endif
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
c       Version v5.1.03
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c*******TO DERIVE RELATIVE CHANGE IN IONIC POPULATIONS : DIF
c     LIM : NUMBER OF LAST ELEMENTS TO BE INCLUDED
c     TRE : TRESHOLD FOR MEASURING CHANGE
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      subroutine difhhe (popin, popfi, dif)
c 
      include 'cblocks.inc'
c 
c 
      real*8 popin(mxion, mxelem), popfi(mxion, mxelem)
      real*8 dif, thresh, maxdif, tmdif, sum, p1, p2
      integer*4 num, iel, ion, lim, matom, mion
c 
      dif=0.0d0
      sum=0.d0
      num=0
      lim=2
      thresh=1.0d-9
c 
      matom=0
      mion=0
      maxdif=0.d0
c 
      do iel=1,lim
        do ion=1,maxion(iel)
          if ((popin(ion,iel).ge.thresh).and.(popfi(ion,iel).ge.thresh))
     &      then
            p1=dlog10(popin(ion,iel))
            p2=dlog10(popfi(ion,iel))
            tmdif=dabs(p1-p2)
            if (tmdif.ge.maxdif) then
              maxdif=tmdif
              matom=iel
              mion=ion
            endif
            sum=sum+tmdif
            num=num+1
          endif
        enddo
      enddo
c 
      if (num.gt.0) dif=sum/(num+epsilon)
c 
      if (expertmode.gt.0) then
        write (*,*) 'Max HHe difference: ',elem(matom),rom(mion),'=',
     &   maxdif
      endif
c 
      return
c 
      end
