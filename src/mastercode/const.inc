cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   include file for global constants.
c   (2014 CODATA recommended values, IAU2015/2012 Astron Const)
c
c   Version v5.1.13 MAPPINGS V
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
c 3D Monte-Carlo Module parameters
c
      integer*4 safelmt
      integer*4 looplmt
      integer*4 mcnstar
      integer*4 maxnx,maxny,maxnz
c      real*8    mcdx,mcdy,mcdz
      real*8    mcxst,mcyst,mczst
      real*8    mcxctr,mcyctr,mczctr
c
c Math-Floating point limits and constants
c
      real*8 pi,ftpi
      real*8 rt2,euler
      real*8 twopi,itwopi
      real*8 fpi,ifpi
      real*8 epi,iepi
      real*8 pi8rt3
      real*8 epsilon, logepsilon
      real*8 huge, loghuge, maxdekt
      real*8 ioepsilon, lg10ioepsilon
c
c Physical-Chemical constants
c
      real*8 plk,cls,ev,eesu
      real*8 me,mp,amu,mh,mh1,alphafs
      real*8 bohr,stefan,nav,rgas
      real*8 r_inf,r_infe,t_ryd
      real*8 iph,iphe,iphlam
      real*8 lmev,evplk,pzlimit,log10pz
      real*8 rka,rkb,ffk,fk,kbb
      real*8 kev,kkev,kevhz
      real*8 gammaEOS
c
c Astro constants
c
      real*8 aunit, pars
      real*8 rsun,msun,lsun,tsun
      real*8 grav,rgm,ledd
c
c     maximum array sizes
c
      integer*4 mxinfph,mxelem,mxion
      integer*4 mxdust,mxdtype,pahi
      integer*4 mxtempbin,mxnnorad
      integer*4 mxnl,mxtr,mxmspecs,mxnspl
      integer*4 mxkappas,mxspeclines
      integer*4 mxnsteps,dinfph
      integer*4 mxnegdes,mxnrec,mxncol
      integer*4 mxfespecs,mxfenl,mxfetr
      integer*4 mxifsteps,mxpcits,mxmonlines

c
c Hydrogenic Spectra, n lines and n series
c
      integer*4 mxnhlines,mxnhs   ! 16/6
      integer*4 mxnhelines,mxnhes ! 16/6
      integer*4 mxnxhlines,mxnxhs ! 10/2
c
c V3 HHe Resonance Collisions data
c
      integer*4 mxhheions
      integer*4 mxhhelvls
      integer*4 mxhhecols
      integer*4 mxhhenspl
c
      integer*4 mxfenspl
c
c V3 Xray Resonance lines data
c
      integer*4 mxxr3ions
      integer*4 mxxr3lvls
      integer*4 mxxr3cols
      integer*4 mxxr3casc
      integer*4 mxxr3lines
      integer*4 mxxr3nspl
c
c V4 Xray Resonance lines data
c
      integer*4 mxxrlions
      integer*4 mxxrllvls
      integer*4 mxxrlcols
      integer*4 mxxrlcasc
      integer*4 mxxrllines
      integer*4 mxxrlnspl
c
c legacy line arrays
c
c      integer*4 mxxlines
c      integer*4 mxxilines
c      integer*4 mxxhelines
c
c stellar libraries
c
      integer*4 mxflxbins
      integer*4 mxoptxbins
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Universal constants: CODATA 2014 Values
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   speed of light in vacuum (cm/s)
c
      parameter  (cls = 2.99792458d10)
c
c   planck's constant (erg seconds)
c
      parameter  (plk = 6.626070040d-27)
c
c     1 electron Volt in ergs
c
      parameter (ev = 1.6021766208d-12)
c other const for optxagnf
      parameter (kev   = 1.6021766208d-9)
      parameter( kkev  = 1.16045220604d+7)
      parameter( kevhz = 2.417989262d17)
c
c     1 electron Volt/h
c
      parameter (evplk = ev/plk)
c
c   electron volt in esu units ev*c*1e-8 statVolt
c
      parameter (eesu = 4.803204672997659d-10)
c
c     rest mass of an electron me (g)
c
      parameter (me = 9.10938356d-28)
c
c     rest mass of a proton mp (g)
c
      parameter (mp = 1.672621898d-24)
c
c     Nav Avogadro's number
c
      parameter (nav = 6.022140857d23)
c
c   boltzman's constant k, ergs per kelvin
c
      parameter (rkb = 1.38064852d-16)
c
c     atomic mass unit assuming C12 = 12
c
      parameter (amu = 1.660539040d-24)
c
c     mean mass of a H atom : 1.00794 amu NIST
c
      parameter (mh = 1.6737237199776d-24)
c
c     mass of a 1H atom : 1.007825 amu
c
      parameter (mh1 = 1.673532757988d-24)
c
c   Rydberg Constant R_inf
c
      parameter (r_inf  = 13.605693009d0)
      parameter (r_infe = r_inf*ev)
      parameter (t_ryd  = r_infe/rkb)
c
c   Hydrogen Rydberg Constant, R_H
c
c   Ionisation potential of hydrogen (ev) Primary Def.
c   109678.7718246123 cm-1 NIST = 1.3598434005136E+01 eV
c
      parameter (iph = 1.3598434005136d1)
c
c
c   Ionisation potential of hydrogen (ergs)
c
      parameter (iphe = iph*ev)
c
c
c   Wavelength of Hydrogen Limit (Angstroms) Secondary
c   911.753495656776 using 1/NIST value
c
      parameter (iphlam = 9.11753495656776d2)
c
c   Wavelength (Angstoms) to Energy (eV) Primary Def.
c   12398.41973964072 using CODATA 2014
c
      parameter (lmev = 12398.41973964072)
c
c alphafs = (2πe^2 /(hc))
c
      parameter (alphafs = 7.2973525664d-3)
c
c   First Bohr radius  (cm)
c
c   a0 = h^2/(4 pi^2 me e^2)
c
c
      parameter (bohr = 5.2917721067d-9)
c
c   Stefan-Boltzmann constant
c   5.670367e-8
c   (ergs cm^-2 s^-1 K^-4)
c   also:
c      2 pi^5 k^4
c   = ------------
c      15 c^2 h^3
c
      parameter (stefan = 5.670367d-5)
c
c   collection of constants for coll.
c         rate = sqrt((2*pi*hbar^4)/(k me^3))
c
c   in ergs and cm units (SI use d-12)
c
      parameter (rka = 8.629133595419446d-6)
c
c     continuum emission constant
c
c    16/3 (pi/6)^0.5 (eesu^6 / (c^3 me^2)) (me/k)^0.5
c  = 16/3 (pi/(6 k me^3))^0.5  (e^2 / c )^3
c
      parameter(ffk = 5.444367585780054d-39)
c
c
c     total FreeFree emission constant
c
c    (16/(3 hbar me)) ((2 pi  k )/(3 me))^(1/2)  ((eesu^2)/c)^3
c
c
      parameter(fk = 1.425556429377619d-27)
c
c constant for planck function 1/((cls**2)*(plk**3))
c
      parameter ( kbb = 3.824639688992588d57)
c
c   universal gas constant
c
      parameter ( rgas = 8.3144598d+7 )
c
c ideal atomic gas poly index gamma, gammaEOS
c
      parameter ( gammaEOS = 5.d0/3.d0 )
c
c astro quantities, IAU2012, IAU2015
c
      parameter ( aunit = 1.49597870700d13)  !cm IAU2012
      parameter ( pars  = 3.085677581491368d18) !cm IAU2012
      parameter ( rsun  = 6.957d10)       !cm  IAU2015
      parameter ( msun  = 1.988475d33 )   !g derived from IAU2015 GMsun
      parameter ( lsun  = 3.828d33 )      !ergs/s  IAU2015
c
c L = (4*pi*rsun^2stefan*tsun^4):
c tsun = (lsun/(4*pi*rsun^2*stefan))^0.25
c 5772 K IAU2015 + CODATA2014
c
      parameter ( tsun = 5772.00d0 ) ! K
c
c
c universal gravity constant
c
      parameter( grav = 6.67408d-8)
c gravitational radius per solar mass , GM/c^2 cm
      parameter( rgm  = 147662.4729623529d0)
c Eddington Luminosity per solar mass
c
c 4*pi*G*Msun*mp*c/Thomson
c
      parameter( ledd = 1.257064901d+38)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Mathematical constants
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   pi and four thirds pi
c   pi with datan to match best to compiler math library
c
      parameter ( pi = 4.d0*datan(1.0d0), ftpi = (4.d0/3.d0)*pi )
      parameter (twopi = 2.d0*pi        , itwopi = 1.d0/(2.d0*pi))
      parameter (fpi = 4.d0*pi          , ifpi =  1.d0/(4.d0*pi))
      parameter (epi = 8.d0*pi          , iepi =  1.d0/(8.d0*pi))
      parameter (pi8rt3=epi*0.5773502691896258d0) ! 8pi/sqrt(3)

c
c   Euler' constant to 20 places
c
      parameter ( euler = 0.57721566490153286061d0 )
c
c   square root of 2
c
      parameter ( rt2 = 1.4142135623730950488d0 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Program constants
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   epsilon is a small number to prevent division by zero
c   and nasty logs. it is also used to test quantities to see
c   if they are worth calculating (see intvec), 1.0d-38 was chosen
c   to match 1E-38 in single precision on a VAX, note:
c   smallest single precision number on sun is ~1.17e-38
c
      parameter (epsilon =1.0d-200, logepsilon=dlog(epsilon))
      parameter (huge    =1.0d+200, loghuge   =dlog(huge)   )
      parameter (ioepsilon=1.0d-99, lg10ioepsilon=dlog10(ioepsilon))
c
c     limit for abundant ions: ab * pop
c
      parameter (pzlimit = 1.0d-12)
      parameter (log10pz = dlog10(pzlimit))
      parameter (maxdekt = 50.d0) ! max DeltaE/kt for collisions
c
c     maximum number of atoms (elements)
c
      parameter (mxelem = 30)
c
c     maximum number of ionization stages
c
      parameter (mxion = 31)
c
c     maximum number of kappa values in spline tables
c
      parameter (mxkappas = 9)
c
c     maximum number of dust sizes (36=MRN, 81=FULL)
c     need to change DUSTDAT files as well (.mrn <==>.full)
c
      parameter (mxdust = 81)
c
c     maximum number of dust types (ie graphite,silicate etc)
c
      parameter (mxdtype = 2)
c
c     Number of PAH ions
c
      parameter(pahi=5)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     radiative and dielectronic recombination rates vector size
c     must matchup in data files, 209 in older 16
c     elements, 465 for full 30 elements
c     mxnnorad is for the subset of NORAD ions, = 26 for C,N,O, Fe 2-6
c
      parameter (mxnrec  = 465)
      parameter (mxnnorad = 26)
c
c     Dere 2007 collisional ionisation rates vector size
c     must matchup with mxnrec data files
c
      parameter (mxncol = 465)
c
c     photoionsation edges vector size, not same as nrec ncol etc
c     696 for std 16 elements, 1696 for full 30 elements
c
      parameter (mxnegdes = 1696)
c
c     radiation vector size, up to 8192 for max photdat
c     1802 for Brents old std and dust must include at least the
c     bin edges for the edges in use.
c
      parameter (mxinfph = 10240)
c
c     model stucture step array limits
c
      parameter (mxnsteps = 5120)
c
c number of steps in S5 ionisation front flow models
c
      parameter (mxifsteps = 1024)
c
c maxmimum number of precursor iterations
c
      parameter (mxpcits = 32)
c
c maximum lines to monitor radially in an P6/7 model
c
      parameter (mxmonlines = 16)
c
      parameter (mxtempbin = 401)
c
      parameter (dinfph = 401)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  multi-level and other collisional line arrays
c
c longest spec file list, including all possible lines and some
c padding for safety
c
      parameter (mxspeclines = 100000)
c
c     maximum number of species in multi-level atoms
c     126 for std 16. 173 for full set
c
      parameter (mxmspecs = 174)
c
c     maximum number of levels in multi-level atoms
c
      parameter (mxnl = 9)
c     maximum number of spline nodes in coll str splines
c     minimum 17 for NORAD rec rate splines
c
      parameter (mxnspl    = 17)
c
c     maximum number of transitions in multi-level atoms
c     mxnl = 5 =>10, mxnl=6 =>15, mxnl=9, => 36, (N*(N-1))/2
c
      parameter (mxtr     = 36)
c
c     equiv numbers for Fe like multi-level atoms
c     new solver can handle up to 322 level fe 3 tested.
c     old mdiag gave up around 50-60 levels
c
      parameter (mxfespecs = 6    ) ! 6
      parameter (mxfenl    = 217  ) ! 108 std density 217 high for fe3
      parameter (mxfetr    = 23436) ! 5778 std, 23436 high. per ion
      parameter (mxfenspl  = 11   ) ! 11 point BT type 3 splines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter (mxnhlines = 15)
      parameter (mxnhs     =  6)
c
      parameter (mxnhelines = 15)
      parameter (mxnhes     =  6)
c
      parameter (mxnxhlines = 10)
      parameter (mxnxhs     =  2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Resonance line cascades - huge memory.
c
      parameter (mxxr3ions=189)     ! 126 std, 189 Full, max ions
      parameter (mxxr3lvls= 86)     !  82 max levels per ion
      parameter (mxxr3cols=504)     ! 504 max collisions per ion
      parameter (mxxr3casc= 72)     !  72 max lines in single cascade
      parameter (mxxr3lines=18522)  ! 18522 full 13990 std prefs
      parameter (mxxr3nspl= 11)     !  11 point BT type 3 splines
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Resonance line cascades - huge memory, was split into
c atom size groups to make multi-d array more efficient.
c
c set for full prefs std density by default
c
c mxxrlcols:
c   5532 std density
c   4389 high density
c
c mxxrllines:
c   56718 std prefs std density
c   58747 full prefs std density
c   57728 std prefs high density
c   59757 full prefs high density
c
      parameter (mxxrlions=45  )   ! 40 std 45 Full ions
      parameter (mxxrllvls=668 )   ! 668 std max lvls per ion
      parameter (mxxrlcols=5532)   ! 5532/4389 max coll. per ion
      parameter (mxxrlcasc=86  )   ! 86 std max lines in single cascade
      parameter (mxxrllines=59757) ! 56718 ss 58747 fs prefs, 59757 fh
      parameter (mxxrlnspl=11)     ! 11 point BT type 3 splines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c HI and HeII collisional cooling (not cascades)
c
      parameter (mxhheions=2 ) ! 2
      parameter (mxhhelvls=25) ! 25
      parameter (mxhhecols=24) ! 24
      parameter (mxhhenspl=11) ! h like colls
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Legacy line data array sizes
c
c      parameter (mxxlines=4) ! 920
c      parameter (mxxilines=9) ! 32
c      parameter (mxxhelines=26) ! 26
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c buffer size to read in all atmospheres, can be reduced to
c 25000 if not using CMFGEN inputs.  A major memory hog!
c
      parameter( mxflxbins = 560000) ! up to 560000 bins - CMFGEN!
c
      parameter( mxoptxbins = 5000)  ! number of energy bins in OPTXAGNF model
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Monte-Carlo Module Parameter
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter (safelmt=10)
      parameter (looplmt=1e3)
c
      parameter (mcnstar = 3)
ccc      parameter (maxnx = 13)
ccc      parameter (maxny = 13)
ccc      parameter (maxnz = 13)
      parameter (maxnx = 41)
      parameter (maxny = 41)
      parameter (maxnz = 41)
c
ccc      parameter (mcdx = 3.0d18)
ccc      parameter (mcdy = 3.0d18)
ccc      parameter (mcdz = 3.0d18) 
cc      parameter (mcdx = .20d18)
cc      parameter (mcdy = .20d18)
cc      parameter (mcdz = .20d18)      
c
      parameter (mcxctr = 0)
      parameter (mcyctr = 0)
      parameter (mczctr = 0)  
c
ccc      parameter (mcxst = -1.950d19-mcxctr)
ccc      parameter (mcyst = -1.950d19-mcyctr)
ccc      parameter (mczst = -1.950d19-mczctr) 
cc      parameter (mcxst = -0.270d19-mcxctr)
cc      parameter (mcyst = -0.270d19-mcyctr)
cc      parameter (mczst = -0.270d19-mczctr)  
      parameter (mcxst = 0.0d0)
      parameter (mcyst = -0.375d18)
      parameter (mczst = -0.375d18)   
c
