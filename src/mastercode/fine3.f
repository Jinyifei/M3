cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES  3-LEVEL ATOMS.
c     currently used for fine structure
c     for some neutral atoms.
c
c       uses proton and H atom collisions for neutralatom
c       fine structure lines, neglecting electons.
c       Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fine3 (t, de, dh)
c
      include 'cblocks.inc'
c
      real*8 t, de, dh
c
      real*8 w(3), e(3, 3), a(3, 3),om(3,3),tk3(3,3),re(3,3)
      real*8 rd(3, 3), alph(5, 4), x(4),t2,t4
      real*8 dhi,dhii
      real*8 aa,aion,bion,f,po,zi
      real*8 tkappa(3,3), ratekappa(3,3)
c
      integer*4 il,jc,jl,jll,m
      integer*4 i,j,k,l,nl, atom, ion
c
      real*8 fkenhance
c
      nl=3
      f3loss=0.0d0
c
c turn off proton coll cooling for line only test, even though these
c are lines....
c
      if (linecoolmode.eq.1) return
c
      f=dsqrt(1.0d0/t)
      t2=t*1.0d-2
      t4=t*1.0d-4
c
      ion=nf3ions
      do i=1,ion
c
c     only calculate abundant ions
c
        do j=1,nf3trans
          f3bri(j,i)=0.d0
        enddo
c
        atom=f3atom(i)
        ion=f3ion(i)
        po=pop(ion,atom)
        zi=zion(atom)
c
        if (po*zi.ge.pzlimit) then
c
c
          dhi=pop(1,zmap(1))*dh
          dhii=pop(2,zmap(1))*dh
c
          do j=1,nl
            w(j)=wi3(j,i)
          enddo
c
c     clear arrays
c
          do k=1,nl
            do l=1,nl
              a(k,l)=0.d0
              om(k,l)=0.d0
              tk3(k,l)=0.d0
              e(k,l)=0.d0
            enddo
          enddo
c
c     gaps
c
          tk3(1,2)=ei3(1,i)
          tk3(1,3)=ei3(2,i)
          tk3(2,3)=ei3(3,i)
          e(1,2)=ei3(1,i)*rkb
          e(1,3)=ei3(2,i)*rkb
          e(2,3)=ei3(3,i)*rkb
c
c     and symmetrically for the algorithm used here
c
          tk3(2,1)=tk3(1,2)
          tk3(3,1)=tk3(1,3)
          tk3(3,2)=tk3(2,3)
          e(2,1)=e(1,2)
          e(3,1)=e(1,3)
          e(3,2)=e(2,3)
c
c     trans probs upper triangle
c
          a(2,1)=ai3(1,i)
          a(3,1)=ai3(2,i)
          a(3,2)=ai3(3,i)
c
c     calculate oms from de-exitaion rates
c     again symmetrical for algorithm
          om(1,2)=0.d0
          om(1,3)=0.d0
          om(2,3)=0.d0
c
c     type 1 uses 100s of K and proton (electron for ions)
c          + hydrogen rates
c
c     type 2 uses 10000s of K in two powers for electron rates only
c
c     add rates proportionally, so multiply by dhi later will scale..
c
c
          if (f3type(i).eq.1) then
c
            if (f3ion(i).eq.1) then
c
c     proton rates for neutrals
c
              if (dhi.gt.0.d0) then
                om(1,2)=dhii/dhi*gam3(1,1,i)*(t2**pow3(1,1,i))
                om(1,3)=dhii/dhi*gam3(1,2,i)*(t2**pow3(1,2,i))
                om(2,3)=dhii/dhi*gam3(1,3,i)*(t2**pow3(1,3,i))
              endif
c
            else
c
c     electron rate for ions
c
              if (dhi.gt.0.d0) then
                om(1,2)=de/dhi*gam3(1,1,i)*(t2**pow3(1,1,i))
                om(1,3)=de/dhi*gam3(1,2,i)*(t2**pow3(1,2,i))
                om(2,3)=de/dhi*gam3(1,3,i)*(t2**pow3(1,3,i))
              endif
c
            endif
c
c     add hydrogen proportionally, so multiply by dhi later will scale..
c
            om(1,2)=om(1,2)+gam3(2,1,i)*(t2**pow3(2,1,i))
            om(1,3)=om(1,3)+gam3(2,2,i)*(t2**pow3(2,2,i))
            om(2,3)=om(2,3)+gam3(2,3,i)*(t2**pow3(2,3,i))
c
c     and symmetrically for the algorithm used here
c
            om(2,1)=om(1,2)
            om(3,1)=om(1,3)
            om(3,2)=om(2,3)
c
          endif
c
          if (f3type(i).eq.2) then
c
            if (dhi.gt.0.d0) then
              if (t4.le.1.d0) then
                om(1,2)=de/dhi*gam3(1,1,i)*(t4**pow3(1,1,i))
                om(1,3)=de/dhi*gam3(1,2,i)*(t4**pow3(1,2,i))
                om(2,3)=de/dhi*gam3(1,3,i)*(t4**pow3(1,3,i))
              endif
              if (t4.gt.1.d0) then
                om(1,2)=de/dhi*gam3(2,1,i)*(t4**pow3(2,1,i))
                om(1,3)=de/dhi*gam3(2,2,i)*(t4**pow3(2,2,i))
                om(2,3)=de/dhi*gam3(2,3,i)*(t4**pow3(2,3,i))
              endif
            endif
c
          endif
c
c     convert to effective omegas, so the old routine it can convert
c     back...
c
          do j=1,nl
            do m=1,nl
              om(m,j)=om(m,j)/(rka*f)
              ratekappa(m,j)=1.d0
              if (usekappa) then
                if (m.ne.j) then
                  tkappa(m,j)=tk3(m,j)/t
                else
                  tkappa(m,j)=0.0d0
                endif
                ratekappa(m,j)=fkenhance(kappa,tkappa(m,j))
              endif
            enddo
          enddo
c
c     and symmetrically for the algorithm used here
c
          om(2,1)=om(1,2)
          om(3,1)=om(1,3)
          om(3,2)=om(2,3)
c
c     if ((mapz(f3atom(i)).eq.8).and.(f3ion(i).eq.1)) then
c     c
c
c     write(*,'(/,5hfine3,/3(x,1pg11.3)/)') w(1),w(2),w(3)
c     c
c     write(*,'(2(x,1pg11.3))') gam3(1,1,i),pow3(1,1,i)
c     write(*,'(2(x,1pg11.3))') gam3(1,2,i),pow3(1,2,i)
c     write(*,'(2(x,1pg11.3)/)') gam3(1,3,i),pow3(1,3,i)
c
c      write(*,'(3(x,1pg11.3))') e(1,1),e(1,2),e(1,3)
c      write(*,'(3(x,1pg11.3))') e(2,1),e(2,2),e(2,3)
c      write(*,'(3(x,1pg11.3)/)') e(3,1),e(3,2),e(3,3)
c
c      write(*,'(3(x,1pg11.3))') tk3(1,1),tk3(1,2),tk3(1,3)
c      write(*,'(3(x,1pg11.3))') tk3(2,1),tk3(2,2),tk3(2,3)
c      write(*,'(3(x,1pg11.3)/)') tk3(3,1),tk3(3,2),tk3(3,3)
c
c     write(*,'(3(x,1pg11.3))') a(1,1),a(1,2),a(1,3)
c     write(*,'(3(x,1pg11.3))') a(2,1),a(2,2),a(2,3)
c     write(*,'(3(x,1pg11.3)/)') a(3,1),a(3,2),a(3,3)
c
c      write(*,'(3(x,1pg11.3))') om(1,1),om(1,2),om(1,3)
c      write(*,'(3(x,1pg11.3))') om(2,1),om(2,2),om(2,3)
c      write(*,'(3(x,1pg11.3)/)') om(3,1),om(3,2),om(3,3)
c
c     endif
c
c
c     continue with old forbid style routine, modified for
c     three rather than 5 levels.  No approximations here,
c     transitions considered between upper levels as well.
c
c
c     ***SET UP COLL. EXCIT. AND DEEXCIT. MATRIX
c
c
          do jl=1,nl
            do jc=1,nl
              rd(jc,jl)=0.0d0
              if (jc.gt.jl) then
                rd(jc,jl)=((rka*f)*om(jc,jl))/w(jc)
                rd(jc,jl)=rd(jc,jl)*ratekappa(jc,jl)
              endif
              re(jc,jl)=0.0d0
              aa=e(jc,jl)/(rkb*t)
              if (aa.lt.maxdekt) then
                if (jc.lt.jl) then
                  re(jc,jl)=(((rka*f)*om(jc,jl))*dexp(-aa))/w(jc)
                  re(jc,jl)=re(jc,jl)*ratekappa(jc,jl)
                endif
              endif
            enddo
          enddo
c
c     ***SET UP MATRIX ELEMENTS
c
          do 50 jl=1,nl
            do 40 jc=1,nl
              if (jl.eq.1) goto 10
              if (jl.eq.jc) goto 20
              alph(jc,jl)=(dhi*(rd(jc,jl)+re(jc,jl)))+a(jc,jl)
              goto 40
   10         alph(jc,jl)=1.0d0
              goto 40
   20         alph(jc,jl)=0.0d0
              do 30 jll=1,nl
                if (jc.eq.jll) goto 30
                alph(jc,jl)=alph(jc,jl)-(dhi*(rd(jc,jll)+re(jc,jll))+
     &           a(jc,jll))
   30         continue
   40       continue
   50     continue
c
          do jl=1,nl
            alph(4,jl)=0.0d0
          enddo
          alph(4,1)=1.0d0
c
c     ***SOLVE FOR LEVEL POPULATIONS  X
c
          jl=3
c
c          call mdiag3 (jl, alph, x)
          call  matsolve (alph, x, 3, 3)
C
C     ***CHECK ON NORMALISATION
C
C         da=0.0d0
C         do jc=1,nl
C           if (x(jc).lt.0.0d0) x(jc)=0.0d0
C           if (x(jc).gt.0.0d0) da=da+x(jc)
C         enddo
C         do jc=1,nl
C           x(jc)=x(jc)/da
C         enddo
Cc
c     ***GET ION COOLING
c
          aion=zi*po*dh
          il=1
          do jl=1,nl-1
            do jc=jl+1,nl
              bion=x(jc)*e(jc,jl)*aion*a(jc,jl)
              f3bri(il,i)=bion*ifpi
              il=il+1
              coolz(atom)=coolz(atom)+bion
              coolzion(ion,atom)=coolzion(ion,atom)+bion
              f3loss=f3loss+bion
            enddo
          enddo
c
c     end abundant ion
c
        endif
c
c     next ion..
c
      enddo
c
      if (f3loss.lt.epsilon) f3loss=0.d0
c
      return
      end
