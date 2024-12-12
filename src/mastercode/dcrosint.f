      subroutine dcrosint (inl,kabp,ksct)
c
c     Integrate Dust kappa absoption and kappa scatter 
c     over grain distribution
c
c
      include 'cblocks.inc'
c
      real*8 kabp,ksct,dgrad
      integer*4 k,inl,dtype
c
      kabp = 0.0d0
      ksct = 0.0d0
c      
      do dtype=1,numtypes
       do k=mindust(dtype),maxdust(dtype)
c       
         dgrad=gradedge(k+1)-gradedge(k)
c
         kabp = kabp + absorp(inl,k,dtype)
     &        * dustsig(k,dtype) * dgrad
         ksct = ksct + scatter(inl,k,dtype)
     &        * dustsig(k,dtype) * dgrad
       
       enddo
      enddo
c
      return 
c      
      end
