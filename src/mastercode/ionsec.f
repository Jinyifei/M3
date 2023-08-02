cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS V.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1979-2012
c
c       Version v5.1.13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******CALCULATES IONISATION RATE DUE TO SECONDARY ELECTRONS
c     FOR THE ELEMENTS INCLUDED IN MATRIX SECEL
c     OUTPUT RATES IN RASEC(3,11)
c
c     REF.: SHULL,M.J. APJ.234,P761 (1979)
c           BERGERON,J.,SOUFFRIN,S. A&A 14,P167
c
c
      subroutine ionsec ()
c
c
      include 'cblocks.inc'
c
c
c           Variables
c
      real*8 a1,a2,abr,abwe,aij,ath,atx,b1
      real*8 b2,deno,eah,eau1,eau2,eax,eij,emin
      real*8 eps,epsh,epsij,fhiieff,frio
      real*8 poef,poij,ps,valel,wei,wij,wth,wtx
      real*8 x,zshull,ztot
      real*8 sig,beff, rate
c
      integer*4 i,j,k
      integer*4 atom,ion
c
      sig(eps)=dmin1(1.0d0-(1.0d0/eps),1.0d0/eps)+epsilon
      beff(x)=(0.0352d0*(1.d0-((x)**(-(0.1d0*dlog10(x))))))/dmax1(1.d-
     &20,1.d0-x)
c
      eau1=iph
      eau2=30.0d0
      emin=iph
c
      zshull=1.1d0
      do 20 i=1,atypes
        do 10 j=1,3
          rasec(j,i)=0.0d0
   10   continue
   20 continue
c
      if ((anr(1,1,1).le.0.0d0).and.(wnr(1,1,1).le.0.0d0)) goto 70
      wth=0.0d0
      wtx=0.0d0
      ath=0.0d0
      atx=0.0d0
c
      do 40 i=1,atypes
        do 30 j=1,maxion(i)-1
          abr=zion(i)*pop(j,i)
          ath=ath+(abr*anr(1,j,i))
          atx=atx+(abr*anr(2,j,i))
          wth=wth+(abr*wnr(1,j,i))
          wtx=wtx+(abr*wnr(2,j,i))
   30   continue
   40 continue
c
      if (ath.le.0.0d0) goto 70
      eah=(ath/(wth+1.d-38))/ev
      eax=(atx/(wtx+1.d-38))/ev
      deno=dlog(eau2)-dlog(eau1)
      b1=(dlog(atx+1.d-38)-dlog(ath+1.d-38))/deno
      a1=dlog(ath+1.d-38)-(b1*dlog(eau1))
      b2=(dlog(eax+1.d-38)-dlog(eah+1.d-38))/deno
      a2=dlog(eah+1.d-38)-(b2*dlog(eau1))
      ztot=0.0d0
      ps=0.0d0
      wei=0.0d0
c
c
      do 50 k=1,25
        atom=secat(k)
        if (atom.le.0) goto 50
        ion=secio(k)
        poij=ipote(ion,atom)/ev
        poef=dmax1(poij,emin)
        eij=dexp(a2+(b2*dlog(poef)))+poef
        aij=dexp(a1+(b1*dlog(poef)))
        epsij=eij/poij
        epsh=eij/eau1
        valel=dble(secel(k))
        wij=((aij*(valel*((eau1/poij)**2)))*sig(epsij))/sig(epsh)
        secra(k)=wij
        ztot=ztot+zion(atom)
        abwe=zion(atom)*wij
        ps=ps+(abwe*pop(ion,atom))
        wei=wei+abwe
   50 continue
c
      fhiieff=dmax1(1.d-20,1.0d0-(ps/(wei+1.d-36)))
      frio=((zshull*beff(fhiieff))/ztot)/ev
c
      do 60 k=1,25
        atom=secat(k)
        if (atom.le.0) goto 60
        ion=secio(k)
        rate=dmax1(frio*secra(k),0.d0)
        if (rate.lt.1.d-28) rate=0.d0
        rasec(ion,atom)=rate
   60 continue
c
   70 return
      end
