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
c     Subroutine freebound calculates the free-bound continuum contribution
c     to the local diffuse field..  H like atoms are treated with lookuptables
c     with muliple levels.  Other atoms are derived from photo-crossections
c     with an additions excited level contribution summed into the first
c     edge above ground.
c
c     adds the emission to cnphot (in photons/cm3/s/Hz/sr)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine freebound (t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      real*8 t,de,dh,en2,abde
      real*8 gfbi,zpop,energ
      real*8 rydu,rz3,t12
      real*8 z2,z4,fbg,lgip,ip
      real*8 tz2,lte,ft0,ft1,ltelo,ltehi
      real*8 fj0,fj1,elo,ehi,gsum
      real*8 f1,df,ri,fbconst
      real*8 lgcebin,eekt,t32
      integer*4 i,atom,ion,isos
c
c      real*8 at,bet,se,eph
      real*8 abio,crosec
      real*8 ebinc,enu
      real*8 excerg,phots,een2
      real*8 phoi,rkt,statwei,et
c
c      real*8 fbloss,lum,wid
c
      integer*4 io,ie,j,inl,n0
      integer*4 k,kthi,ktlo
      integer*4 jelo,jehi,fb
c
c     Functions
c
      real*8 hir
c
c
c  -138.446097d0 = ln((h^3)/(me^3 c^2)  (me/2)^(1.5)  4/sqrt(¹) /(4¹))
c
c     Osterbrok Appendix 1
c
c Free-Bound Milne relation for photo-recombinations, in log space.
c
      hir(rkt,enu,excerg)=dexp((-138.446097d0+(2.0d0*dlog(enu)))-(1.5d0*
     &dlog(rkt))-(excerg/rkt))
c      acrs(eph,at,bet,se)=(at*(bet+((1.0d0-bet)/eph)))/(eph**se)
c
      rkt=rkb*t
      rydu=(iphe)/(rkt)
      t12=1.d0/dsqrt(t)
      t32=t**(-1.5d0)
c
      do inl=1,infph
        fbph(inl)=0.d0
      enddo
c
c      fbloss=0.d0
c
c     The Reimann zeta function for argument 3
c
      rz3=1.202056903159594d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ground level fb continuum from photo X-sections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ground level free-bound continuum. Transplant from localem.
c     This is better because it is based on the actual photo
c     cross-sections in use, rather than an H approximation.
c     **ADS RECAPTURE TO GROUND STATE FOR ELEMENTS UP TO*LIMPH
c     STATWEI = RATIO OF STATISTICAL WEIGHT OF RECOMB. TO ION
c     EX. HYDR.: STATWEI 2/1=2
c
      fbg=0.0d0
      fb=0
c
      do io=1,ionum
c
        ie=atpho(io)
        j=ionpho(io)
c
        abio=zion(ie)*pop(j+1,ie)
        statwei=stwtpho(io)
        if (statwei.gt.0.d0) fb=1
        if (phionmode.eq.0) fb=isfbpho(io)
c
        if ((abio.gt.pzlimit).and.(fb.eq.1)) then
          phoi=ipotpho(io)
          abio=abio*de*dh*statwei
          do inl=photbinstart(io),infph-1
            ebinc=cphotev(inl)
            if (ebinc.ge.phoi) then
              energ=ebinc*ev
              et=(energ)/(rkt)
              if (dabs(et).lt.loghuge) then
                enu=ebinc*evplk
                crosec=photxsec(io,inl)
                excerg=ev*(ebinc-phoi)
                fbg=crosec*abio*hir(rkt,enu,excerg)
                fbph(inl)=fbph(inl)+fbg
              endif
            endif
          enddo
c
        endif
c
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Excited level fb continuum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Add recombination t higher levels
c     loop through all ions present
c
      do atom=1,atypes
        do ion=1,maxion(atom)-1
          zpop=zion(atom)*pop(ion+1,atom)
          if (zpop.ge.pzlimit) then
c
            abde=zpop*de*dh
            fbconst=ffk*abde
c
            isos=(mapz(atom)-ion)+1
c
            if (isos.eq.1) then
c
c free bound continuum to n = 2 .. 10 summing 350 levels for hydro species
c
              z2=dble(ion)
              z2=z2*z2
              z4=z2*z2
c
              ip=ipotev(ion,atom)
              lgip=dlog10(ip)
c
              tz2=t/z2
              lte=dlog10(tz2)
c
c find temperature columns in logTeff for interpolation
c by bisection search, different for each atom
c
              ktlo=1
              kthi=81
   10         if (kthi-ktlo.gt.1) then
                k=(kthi+ktlo)/2
                if (fblte(k).gt.lte) then
                  kthi=k
                else
                  ktlo=k
                endif
                goto 10
              endif
              ltelo=fblte(ktlo)
              ltehi=fblte(kthi)
c
c get table temperature axis interpolation fractions
c
              ft1=(lte-ltelo)/(ltehi-ltelo)
              ft1=dmin1(1.d0,dmax1(0.0d0,ft1))
              ft0=1.d0-ft1
c
c       write(*,*) atom,ltelo, ltehi, ft0, ft1
c
c  now loop through inl, getting j and interpolation fractions for gamma
c  table from recalculated arrays
c
              do inl=1,infph-1
                phots=0.d0
c
                jelo=jfb(atom,inl)
c
                if (jelo.gt.0) then
c
                  eekt=fbee(atom,inl)/rkt
c
                  if (eekt.lt.maxdekt) then
c
                    lgcebin=lgcphotev(inl)
                    jehi=jelo+1
c
                    elo=lgip+fbenu(jelo)
                    ehi=lgip+fbenu(jehi)
c
                    fj1=(lgcebin-elo)/(ehi-elo)
                    fj1=dmin1(1.d0,dmax1(0.0d0,fj1))
                    fj0=1.d0-fj1
c
c Bi-linear log-log interpolation from 4 nearest points
c
                    gsum=ft0*fj0*fbloggamma(ktlo,jelo)
                    gsum=gsum+ft1*fj0*fbloggamma(kthi,jelo)
                    gsum=gsum+ft0*fj1*fbloggamma(ktlo,jehi)
                    gsum=gsum+ft1*fj1*fbloggamma(kthi,jehi)
c                    if (mapz(atom).eq.1) write (*,*) lgcebin,gsum
c
                    phots=1.0d-34*z4*dexp(-eekt)*t32*(10.d0**gsum)
                    phots=abde*phots/cphote(inl)
c
                  endif
c
                endif
c
                fbph(inl)=fbph(inl)+phots
c
              enddo
c
            else if (isos.eq.2) then
c
c free bound continuum to edge = 1.. 50  helio species
c
              z2=dble(ion-1)
              z2=z2*z2
              z4=z2*z2
c
              ip=ipotev(ion,atom)
              lgip=dlog10(ip)
c
              tz2=t/z2
              lte=dlog10(tz2)
c
c find temperature columns in logTeff for interpolation
c by bisection search, different for each atom
c
              ktlo=1
              kthi=76
   20         if (kthi-ktlo.gt.1) then
                k=(kthi+ktlo)/2
                if (hefblte(k).gt.lte) then
                  kthi=k
                else
                  ktlo=k
                endif
                goto 20
              endif
              ltelo=hefblte(ktlo)
              ltehi=hefblte(kthi)
c
c get table temperature axis interpolation fractions
c
              ft1=(lte-ltelo)/(ltehi-ltelo)
              ft1=dmin1(1.d0,dmax1(0.0d0,ft1))
              ft0=1.d0-ft1
c
c       write(*,*) atom,ltelo, ltehi, ft0, ft1
c
c  now loop through inl, getting j and interpolation fractions for gamma
c  table from recalculated arrays
c
              do inl=1,infph-1
                phots=0.d0
c
                jelo=jhefb(atom,inl)
c
                if (jelo.gt.0) then
c
                  eekt=hefbee(atom,inl)/rkt
c
                  if (eekt.lt.maxdekt) then
c
                    lgcebin=lgcphotev(inl)
                    jehi=jelo+1
c
                    elo=lgip+hefbenu(jelo)
                    ehi=lgip+hefbenu(jehi)
c
                    fj1=(lgcebin-elo)/(ehi-elo)
                    fj1=dmin1(1.d0,dmax1(0.0d0,fj1))
                    fj0=1.d0-fj1
c
c Bi-linear log-log interpolation from 4 nearest points
c
                    gsum=ft0*fj0*hefbloggamma(ktlo,jelo)
                    gsum=gsum+ft1*fj0*hefbloggamma(kthi,jelo)
                    gsum=gsum+ft0*fj1*hefbloggamma(ktlo,jehi)
                    gsum=gsum+ft1*fj1*hefbloggamma(kthi,jehi)
c                    if (mapz(atom).eq.2) write (*,*) lgcebin,gsum
c
                    phots=1.0d-34*z4*dexp(-eekt)*t32*(10.d0**gsum)
                    phots=abde*phots/cphote(inl)
c
                  endif
c
                endif
c
                fbph(inl)=fbph(inl)+phots
c
              enddo
c
            else
c
c     free bound continuum for n = 2 non-hydro species
c     using old Mewe 1986 methods, will update later
c
c     get secondary ionisation edge data
c
              en2=ipot2e(ion,atom)
              een2=dexp(en2/rkt)
c
              z2=dble(ion)
              z2=z2*z2
              z4=z2*z2
c
c     calculate f1 factor (a sum over n)
c     term f1 (reimann zeta(3) - sum n^-3 up to no) eqn 14
c
              f1=1.d0
              n0=fn0(ion,atom)
              do i=2,n0
                ri=dble(i)
                f1=f1+1.0d0/(ri*ri*ri)
              enddo
c
              f1=(rz3-f1)
c
              df=2.d0*z4*f1*een2
              gfbi=1.1d0*rydu*df
c
              do inl=1,infph-1
c
                energ=cphote(inl)
c
                if (energ.ge.en2) then
c
                  et=(energ)/(rkt)
c
                  if (dabs(et).lt.maxdekt) then
c
c    (t12? yes, gfbi contains another 1/T factor in rydu)
c    const units include T_6 ->T and 1/A -> Hz
c
                    phots=fbconst*gfbi*t12*dexp(-et)/energ
                    fbph(inl)=fbph(inl)+phots
c
c              end et limit
c
                  endif
c
c               end en2 limit
c
                endif
c
c              end radiation loop
c
              enddo
c
            endif
c
c     end if loop for zpop > pzlimit
c
          endif
c
c     end ion and atom loops
c
        enddo
c
      enddo
c
c test
c
C      lum=0.d0
C      do inl=1,infph-1
C        energ=cphote(inl)
C        wid=(photev(inl+1)-photev(inl))*evplk ! save a div, use *eV/plk
C        lum=lum+fbph(inl)*cphote(inl)*wid
C      enddo
C      fbloss=lum*fpi
c
      return
c
      end
