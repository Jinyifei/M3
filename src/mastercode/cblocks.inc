cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MAPPINGS V v5.1.13 COMMON BLOCK GLOBAL VARIABLES
c
c     common blocks are all in this one include file
c     this ensures identical declarations for all files.
c
c     Need documentation and more organisation....
c     will document each variable as they are edited.
c
c     (c) 1979-2016
c
c     Version v5.1.13
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'const.inc'
c
c     Need some global enviroment variables
c
c     theversion - the version
c
c     dtlen: actual length of datadir string, set once in mapinit
c     datadir: location of data files,
c              './' or '/usr/local/share/mappings/'
c
c     filename: scratch array for files in general
c
      character*8 theversion
      integer*4 dtlen
      character*512 datadir
      character*512 filename
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   declare ATDAT block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   zion: gas phase abundance fractions
c   deltazion: abundance adjustments (mulipliers)
c   dion: depletion (dust) factors (mulipliers)
c   zion = zion0*deltazion*dion
c   zion0: base total (undepleted) abundance fractions
c   zsol: solar values
c   invdion: inverse of depletion factors
c   zen: sum of all species abundance fractions
c           ie zen*dh = total number of ions
c   atwei: atomic weights
c   maxion: maximum ionisation stage for a given atom type
c   zmap&mapz : convert atomic numbers to array entries and v v
c   elem:   Atomic symbols
c   rom :   Ion roman numerals
c
c   atypes: number of different atoms used (up to mxelem)
c
      integer*4 atypes
      integer*4 maxion(mxelem)
      integer*4 zmap(mxion),mapz(mxion)
      integer*4 elem_len(mxelem),rom_len(mxion)
c
      character elem(mxelem)*4,rom(mxion)*8
c
      real*8 atwei(mxelem),zen
      real*8 zion(mxelem),zion0(mxelem)
      real*8 deltazion(mxelem),zsol(mxelem)
      real*8 dion(mxelem),dion0(mxelem)
      real*8 invdion(mxelem)
c
c   pop: ion fractional populations
c   dndt: rate of change of pop with time step
c   popint: integrated column densities
c
      real*8 pop(mxion,mxelem)
      real*8 dndt(mxion,mxelem)
      real*8 popint(mxion,mxelem)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   col: collisional ionisation rates by ion and atom
c
      real*8 col(mxion,mxelem)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Dere 2007 collisional ionisation coeffs
c   see collion.f and COLLDAT2
c
      integer*4 ncol,atcol(mxncol),ioncol(mxncol),nsplcol(mxncol)
      real*8 xcol(21,mxncol),ycol(21,mxncol),y2col(21,mxncol)
      real*8 ipotcol(mxncol)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Arnaud and Rothenflug collisional ionisation coeffs
c   see collion.f and COLLDAT
c      declare ARCOL block variables
c   arrays by shells,ion,elem + nshells to control loops
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   col: collisional ionisation rates by ion and atom
c
      real*8 collpot(3,mxion,mxelem)
      real*8 aar(3,mxion,mxelem),bar(3,mxion,mxelem)
      real*8 car(3,mxion,mxelem),dar(3,mxion,mxelem)
      integer*4 nshells(mxion,mxelem)
c
c   Landini& Monsignori Fosse & Shull coll. ion. fits
c
      real*8 acol(mxion,mxelem),tcol(mxion,mxelem)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Original mappings ionisation coeffs
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   (also uses ipote from photio)
c
c      real*8 beta(6,11),effn(6,11)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare PHOTIO block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   ionum : number of photoioisation edges in actual use
c   ipotev: ion valence ionisation potential (ev)
c   ipote: ion valence ionisation potential (ergs)
c   epot2ev: ion ionisation potetial from no+1 level (ev)
c   epot2e: ion ionisation potetial from no+1 level (ev)
c       (used in ffbgau.f)
c       atpho,ionpho: atom and ion entry id for photoionisation
c                     data arrays
c       ipotpho: ionisation edge potential
c       augpho: Auger electron potential
c       sigpho: threshold crossection
c       betpho,spho: beta and s for crossection as f(e) fit
c       stwtpho: stat weight for ion ground state (2J+1) RS coupling
c   secel: secondary ionisation data array
c   xsec : mean photo crossection for each bin
c   xsect: total crossection over ionisising band
c   .
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Common photoionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer*4 ionum
c
      real*8 ipote(mxion,mxelem)
      real*8 ipot2e(mxion,mxelem)
      real*8 ipotev(mxion,mxelem)
      real*8 ipot2ev(mxion,mxelem)
      real*8 secra(26),heaph(5,mxion,mxelem)
      integer*4 phionmode,ph2mode
      integer*4 secat(26),secio(26),secel(28)
      integer*4 nhipx
      real*8 hipxlge (1024)
      real*8 hipxlgs (1024)
      real*8 hipxlgs2(1024)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Monte-Carlo variable
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Ionization/Emission parameter 
c
c     sigmt_arr:  opacity array
c     strPDF_arr:  continuum energy distribution
c     recPDF_arr:  continuum energy distribution
c     totLine_arr: energy fraction in non-ionizing lines: 
c     jnu_arr:   intensity field
c
c     hd_arr:  H density array
c     de_arr:  electron density array
c     te_arr:  temperature array
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c
      real*8 convrt,mcdengy0
      real*8 mcdengy(8)
      real*8 mcnuin,mcnuend   
      real*8 mcnu(7)        
      real*8 sigmt_arr(maxnx,maxny,maxnz,mxinfph)
      real*8 rec0_arr(maxnx,maxny,maxnz)
      real*8 recpdf_arr(maxnx,maxny,maxnz,mxinfph)
      real*8 strpdf_arr(mcnstar,mxinfph)
c      real*8 totline_arr(maxnx,maxny,maxnz)
      real*8 jnu_arr(maxnx,maxny,maxnz,mxinfph)
      real*8 jnuw_arr(2,mxinfph)
      real*8 resflux(mxinfph),resflux2(mxinfph)
c      real*8 xarr(maxnx,maxny,maxnz)
c      real*8 yarr(maxnx,maxny,maxnz)
c      real*8 zarr(maxnx,maxny,maxnz)
      real*8 dh_arr(maxnx,maxny,maxnz)
      real*8 de_arr(maxnx,maxny,maxnz)
      real*8 te_arr(maxnx,maxny,maxnz) 
      real*8 volfilar(maxnx,maxny,maxnz)      
      real*8 cwx(maxnx+1)      
      real*8 cwy(maxny+1)      
      real*8 cwz(maxnz+1)      
c
      real*8 strpos(3)    
c
      integer*4 nrange,npck,npckinc     
      integer*4 taskerr,taskid,numtasks
      integer*4 mcnx,mcny,mcnz
      integer*4 xps, yps, zps
      real*8    intpt
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Verner Yakovlev Ferland Photionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer*4 atpho(mxnegdes)
      integer*4 ionpho(mxnegdes)
      integer*4 shellpho(mxnegdes)
      integer*4 lquantpho(mxnegdes)
      integer*4 nelecpho(mxnegdes)
      integer*4 hasph2(mxnegdes)
      integer*4 isfbpho(mxnegdes)
      integer*4 wrpho(mxnegdes)
      integer*4 wipho(mxnegdes)
c
      real*8 ipotpho(mxnegdes)
      real*8 ph1pho(5,mxnegdes)
      real*8 ph2pho(7,mxnegdes)
      real*8 ph2limitpho(mxnegdes)
c
      integer*4 photbinstart(mxnegdes)
      real*8 photxsec(mxnegdes,mxinfph)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Original sigma-Beta-s MAPPINGS I photoionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 augpho(mxnegdes)
      real*8 sigpho(mxnegdes)
      real*8 betpho(mxnegdes)
      real*8 spho(mxnegdes)
      real*8 stwtpho(mxnegdes)
      real*8 anr(2,mxion,mxelem)
      real*8 wnr(2,mxion,mxelem)
      real*8 rasec(mxion,mxelem)
      real*8 rphot(mxion,mxelem)
      real*8 auphot(mxion,mxelem)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Karzas and Latter,1961,ApJSS,6,167 free-bound continuum data
c up to n = 6.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     nklgfb:  number of principle quantum nnumber (6)
c     nklfbe:  number cross section entries,E/E0 (41)
c
c     klpxe   : cross section energy ratio values,ln values
c     klgfb   : freebound gaunt factors,ln values,(i_e,l+1,n)
c     klgfby2 : freebound gaunt factor row spline coefficients
c
      integer*4 nklgfb,nklfbe
      real*8 klpxe(41)
      real*8 klgfb(41,6,6)
      real*8 klgfby2(41,6,6)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare RECOMB block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   arad,xrad: radiative recombination rate coeffs by ion and atom
c   adi,bdi,t0,t1: dielectronic recombination coeffs
c   a-f dilt: low T dielectronic rate coeffs
c   rec: total recombination rates by ion and atom
c
      integer*4 recommode,drrecmode
      real*8 arad(mxion,mxelem),xrad(mxion,mxelem)
      real*8 adi(mxion,mxelem),t0(mxion,mxelem)
      real*8 bdi(mxion,mxelem),t1(mxion,mxelem)
      real*8 adilt(mxion,mxelem),bdilt(mxion,mxelem)
      real*8 cdilt(mxion,mxelem),ddilt(mxion,mxelem)
      real*8 fdilt(mxion,mxelem)
      real*8 rec(mxion,mxelem)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New recombiation variables, for DRRECDAT and RRRECDAT
c
c     can use CHIANTI, VF or compatibility data files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 nrec,atrec(mxnrec),ionrec(mxnrec)
      integer*4 typerrec(mxnrec),typedrec(mxnrec)
      real*8 arrec(6,mxnrec),adrec(8,mxnrec),bdrec(8,mxnrec)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NORAD total recombination splines for C, N, O, recmode = 2
c mxnspl must = 17 minimum
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer*4 nnorad, noradid(mxion,mxelem),nodesndrec(mxnnorad)
      integer*4 atndrec(mxnnorad),ionndrec(mxnnorad),typendrec(mxnnorad)
      real*8 tendrec(mxnspl,mxnnorad)
      real*8 andrec (mxnspl,mxnnorad)
      real*8 a2ndrec(mxnspl,mxnnorad)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare PHOTDAT block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   photev: energies of photon bin edges in eV: from file PHOTDAT
c   cphotev: energies of photon bin centers in eV: from file PHOTDAT
c   phote : energies of photon bin edges in egs: from file PHOTDAT
c   cphote : energies of photon bin centers in egs: from file PHOTDAT
c   lg prefix,precompute dlog10 versions
c   soupho: source photon vector
c   tphot: total local photon vector (ergs)
c   emidif,updif,dwdif: diffuse field photon vectors
c       (local,upstream and downstream)
c   emlin,uplin,dwnlin: resonance emission line arrays with distance
c       multipliers (local,upstream and downstream)
c   widbinev: widths of photon bin edges in eV: from file PHOTDAT
c   skipbin: logical array to optimise calculation of photionisation
c       to skip bin if it has too few photons to bother with.
c   dilradpho effective source radius for finite extent plane parallel
c   models as a funtion of energy,goes with rshock in s4
c   infph: number of bins in vectors
c   linpos: line bin position information.
c
c       ffphot,fbphot,p2phot: continuum spectra (in ergs)
c                            (free-free,fre-bnd,2photon)
c   dustphot: photon field absorbed by dust locally
c
      real*8 photev(mxinfph)
      real*8 lgphotev(mxinfph)
      real*8 cphotev(mxinfph)
      real*8 lgcphotev(mxinfph)
      real*8 phote(mxinfph)
      real*8 lgphote(mxinfph)
      real*8 cphote(mxinfph)
      real*8 lgcphote(mxinfph)
      real*8 widbinev(mxinfph)
c
      real*8 tphot(mxinfph)
      real*8 sphot(mxinfph)
      real*8 soupho(mxinfph),wdil
      real*8 difftot(mxinfph)
      real*8 dilradpho,rshock
      real*8 emidif(mxinfph)
      real*8 updif(mxinfph)
      real*8 dwdif(mxinfph)
      real*8 emidifcont(mxinfph)
      real*8 updifcont(mxinfph)
      real*8 dwdifcont(mxinfph)
      real*8 xsec(mxinfph),xsect
      real*8 dustphot(mxinfph)
      real*8 srcsum,totsum
      integer*4 infph,linpos(2,60)
c
      integer*4 fieldversion,ionstartbin
      logical skipbin(mxinfph)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   declare gaunt factor calculation globals
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   arg,brg,crg,drg,erg: resonace gaunt factor parameters
c
c
c        integer*4 clkup(19,19)
c        real*8 arg(38),brg(38),crg(38),drg(38),erg(38)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare CHARXC block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   charge exchange reaction rate variables.....
c   charco and charte:Ancient mappings: chargemode = 2
      integer*4 nchxold
      real*8 charco(7,26),charte(5,26)
c
c   ref Arnaud and Rothenflug 1985,A.A.Suppl.Ser. 60:425
c   Legacy Old rates: chargemode = 1
c
c   recombination reactions
c
c   nlegacychxr: number of reactions used
c   chxrlegacyat: element #
c   chxrlegacyion: ion stage
c   chxrlegacyx: reactant (H or He)
c   chxrlegacytemp: suggested temp range
c   chxrlegacycos: coefficients
c       1:a
c       2:b
c       3:c
c       4:d
c
c   chxilegacy series equivalent
c   chxilegacycos: coefficients
c       1:a
c       2:b
c       3:c
c       4:delta E (eV)
c
      integer*4 nlegacychxr,chxrlegacyat(55)
      integer*4 chxrlegacyio(55),chxrlegacyx(55)
      real*8 chxrlegacytemp(2,55),chxrlegacycos(4,55)
      real*8 chxrlegacy(55)
c
c   ionisation reactions
c
      integer*4 nlegacychxi,chxilegacyat(15)
      integer*4 chxilegacyio(15),chxilegacyx(15)
      real*8 chxilegacytemp(2,15),chxilegacycos(4,15)
      real*8 chxilegacy(15)
c
c   ref Kingdon and Ferland, Stancil et al ORNL database 2001
c   New rates: chargemode = 0
c
      integer*4 nchxr,chxrat(200),chxrio(200),chxrx(200)
      integer*4 chxrnfit(200),chxrfitid(200)
      real*8 chxrtemp(2,200),chxrcos(4,200),chxr(200)
      real*8 chxrde(200)
c
c   ionisation reactions
c
      integer*4 nchxi,chxiat(40),chxiio(40),chxix(40)
      integer*4 chxinfit(40),chxifitid(40)
      real*8 chxitemp(2,40),chxicos(5,40),chxi(40)
      real*8 chxide(40)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       COOLING block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   cooling rate quantities (ergs)
c   alphaCloss,alphaClaw,alphaC0: artificial powerlaw losses ,
c      index and normalisation used if alphacoolmode = 1
c   colos: collisional ionisation losses
c   fflos: free-free radiation losses
c   cploss: compton losses := cmpcool,cmpheat separate heating
c      and cooling processes
c   chgain: charge exchange heating
c   pgain: photionisation heating
c   rngain: recombination heating
c   gheat: photoelectric heating by large grains
c   tloss :total absolute loss
c   eloss: effective loss
c   egain: effective heating
c   wgain: heating due to pdv work
c   dlos: measure of equilibration
c   tem: current electron teperature
c   fixtemp: fixed electron temp for jthm = 'T' models
c   cosgain:cosmic ray hydrogen heating
c   crate: cosmic ray primary event rate (~10^-17 s^-1)
c   cosphi: cosmic ray secondary ionisations
c
      real*8 tloss,eloss,egain,dlos,fixtemp,wgain
      real*8 colos,fflos,chgain,pgain,rngain,tem
      real*8 cosgain,crate,cosphi,cmplos,cmpcool,cmpheat
      real*8 alphacloss,alphaclaw,alphac0
      real*8 oiii5007loss
c
      real*8 coolzion(mxion,mxelem)
      real*8 heatzion(mxion,mxelem)
      real*8 coolz( mxelem),heatz(mxelem)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Dust Grains
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Bgrain: threshold parameter(large}
c   Yinf: yield parameter (larges)
c   haccom: accomodation parameter
c   segrain,spgrain: sticking for electrons and protons
c   siggrain: area of large grains per H atom
c   Avgrpot: Average grain potential in volts
c   grainpot: grain potential for each size
c   amin,amax: grain size limits (microns)
c   galpha: distribution index,default to n^-3.5
c   gheat,gcool:heating and cooling due to grains
c   paheat:  heating due to PAH molecules
c   pahNc: mean number of C atoms in PAH molecules
c   pahfrac: # of PAHs per H
c   pahCfrac: % of C in PAHs
c   pahyield: the PE yield curve for all PAH ion states
c   pahion: ionisation state of PAH
c   pahIP: Ionisation potential of that PAH state
c   pahZ: fraction of PAHs in that ion state
c   paheng: energy absorbed per sec by PAH (=energy emitted)
c   pahn/iabs: PAH absorbtion curve for positive and neutral/negative
c        states
c   pahn/isca: PAH scattering as above
c   pahn/icos: PAH scattering angle (as above)
c   pahint:integrated column of PAH
c    ionising vector,using same bins.
c   extinct(i,k,j) dust extinction curve for grains,size k and typej,
c    (Q ext values)
c   scatter(i,k,j) dust scattering curves for grains,size k and typej,
c    (Qsca values)
c   absorp(i,k,j) dust absorption curve for grains,size k and typej,
c    (Q abs values)
c   Gcos(i,k,j) dust scattering angle curve for grains,size k and type
c    (<cos(theta)> value)
c
c   dustint:integrated column of Hydrogen for dust extinction
c   dustbinmax: number of dust radius bins,
c   mindust,maxdust: bins corresponding to min and max grainsizes
c
c   npahx,npahy : number of points in the curves
c
c   grainrad,gradedge: grain radius bin centres and edges
c   grainghr,graindgr,grainpgr composition ratios
c   dustnum,dustsig: number and crossection distribution of dust
c   dcrosec(inl): extinction cross-sections by dust at wvl inl
c   dfscacros,dbscacros:forward and back scattering crossections
c
c   IRphot: infrared flux
c
c   ClinPAH: Condition for Carbon to be linked to PAH existence
c
      real*8 grainrad(mxdust),gradedge(mxdust+1)
      real*8 dustmass(mxdtype),siggrain(mxdtype)
      real*8 dustnum(mxdust,mxdtype),dustsig(mxdust,mxdtype)
      real*8 amin(mxdtype),amax(mxdtype),graindens(mxdtype)
      real*8 bgrain,yinf,haccom
      real*8 gheat,gcool
      real*8 segrain,spgrain
      real*8 grainmaxtemp,grainscaletemp
      real*8 grainmaxq,grainscaleq
      real*8 grainghr,graindgr,grainpgr
      real*8 pahint(pahi),dustint
      real*8 grainpot(mxdust,mxdtype),avgpot(mxdtype),galpha
      real*8 extinct(mxinfph,mxdust,mxdtype),
     &           absorp(mxinfph,mxdust,mxdtype)
      real*8 scatter(mxinfph,mxdust,mxdtype),
     &   gcos(mxinfph,mxdust,mxdtype)
      real*8 dcrosec(mxinfph,mxdtype),
     &           dfscacros(mxinfph,mxdtype),dbscacros(mxinfph,mxdtype)
      real*8 irphot(mxinfph),dustav,a_vend
      integer*4 dustbinmax,mindust(mxdtype),maxdust(mxdtype),numtypes
      integer*4 irtemp,clinpah
c
c PAHs
c
      real*8 paheat,pahfrac,pahcfrac,pahnc
      real*8 pahip(pahi),pahz(pahi),pahion(pahi),pahlimit
      real*8 pahyield(mxinfph,pahi),paheng,pahflux(mxinfph)
      real*8 pahnabs(mxinfph),pahnsca(mxinfph),
     &           pahncos(mxinfph),pahnext(mxinfph)
      real*8 pahiabs(mxinfph),pahisca(mxinfph),
     &           pahicos(mxinfph),pahiext(mxinfph)
      character pahend
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 engy(mxinfph)
      real*8 phot(mxinfph)
      real*8 irflux(mxinfph)
      real*8 cabs(mxinfph)
      real*8 photabs(mxinfph)
      real*8 dustsigt(mxinfph)
      real*8 collheat(mxinfph)
      real*8 edphot(dinfph)
      real*8 dengy(dinfph)
      real*8 dde(dinfph)
      real*8 v_e(dinfph)
      real*8 dabsphot(dinfph),dele
      real*8 bbemiss(mxinfph,mxtempbin)
      real*8 t_grid(mxtempbin)
      real*8 t_edge(mxtempbin)
      real*8 h_grid(mxtempbin)
      real*8 hmin(mxtempbin)
      real*8 hmax(mxtempbin)
      real*8 deltah(mxtempbin)
      real*8 invh(mxtempbin,mxtempbin)
      real*8 tr_matrix(mxtempbin,mxtempbin)
      real*8 t_prob(mxtempbin)
      real*8 grainvol(mxdust)
      real*8 atom_no(mxdtype,mxdust)
      real*8 flxdist(mxinfph,mxdust)
      real*8 d_rad(mxdust)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      ** LINE DATA
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Kappa distributions
c      stored in /kappa/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 kappa
      real*8 kappaa,kappab
      integer*4 kappaidx
      logical usekappa,usekappainterp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare XLINE block vaiables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 nialines(mxion,mxelem)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   new resonance (permitted) line data
c   xlam: wavelengths (units?)
c   xrbri: line brightnesses (ergs)
c   rloss: resonance line heat losses (in rline)
c   920 ytransitions included from M&L 1990
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      integer*4 xlines,minxbin,xrat(mxxlines)
c      integer*4 xion(mxxlines),xiso(mxxlines),xbin(mxxlines)
c      real*8 xrlam(mxxlines),xrbri(mxxlines),xbr(mxxlines),xrloss
c      real*8 xejk(mxxlines),xfef(mxxlines),xtrans(mxxlines),brr(19,19)
c      real*8 xegj(mxxlines),fluxx(mxxlines)
c      real*8 emilin(2,mxxlines),uplin(2,mxxlines),dwlin(2,mxxlines)
c
c      integer*4 xilines,xiat(mxxilines)
c      integer*4 xiion(mxxilines),xiiso(mxxilines),xibin(mxxilines)
c      real*8 xilam(mxxilines),xibri(mxxilines),xiloss
c      real*8 xiejk(mxxilines),xiegj(mxxilines)
c      real*8 xomeg(mxxilines),xitr(mxxilines)
c      real*8 xilin(2,mxxilines),xiuplin(2,mxxilines)
c      real*8 xidwlin(2,mxxilines),fluxxi(mxxilines)
c
c     integer*4 xhelines,xheat(mxxhelines)
c     integer*4 xheion(mxxhelines),xhebin(mxxhelines)
c     real*8 xhelam(mxxhelines),xhebri(mxxhelines),xheloss
c     real*8 xhejk(mxxhelines),xhef(mxxhelines),xhebra(16)
c     real*8 xhelin(2,mxxhelines),xheuplin(2,mxxhelines)
c     real*8 xhedwlin(2,mxxhelines)
c     real*8 fheif(mxxhelines)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c New resonance Collision + Cascade data
c data computed from CHIANTI 8.0.2  and incorperated
c as multi-level ions, from which a linear line list is
c derived.  Each ion has multiple cascades from each collision
c with mulitple contributions to each line, which are summed
c and added to the line list for transfer just like previous resonance
c lines.
c
c Line list is derived from multi-level atomic data at startup.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c overall mode, for std (<1e8)  and high (>1e8) densites - 1
c
c
c prefix: xr3 globals
c
c   xr3loss  - version 3 xray resonance cooling, sum of collisions
c   nxr3ions - int number of ions with branching cascade collision data
c   xr3at(idx)    - int mappings atom id, idx = 1.. nxr3ions
c   xr3ion(idx)   - int ion stage 1...  , idx = 1.. nxr3ions
c   xr3ni(idx)    - int number of levels of ion idx,idx = 1.. nxr3ions
c   xr32S12j(idx)    - int first 2s 2S1/2 level of ion idx, 0 for non H or He
c   xr3ei(ni,idx)      - real energy of level ni in cm^-1
c   xr3gi(ni,idx)      - real stat weight of level ni
c   xr3invgi(ni,idx)      - real  1/stat weight of level ni
c   xr3id(ni,idx)      - char*16 name of level ni
c   xr3eij(ni,ni,idx)  - real energy diff from j->i in cm^-1
c                        computed at startup, maps to line list
c                        below, only computed locally at startup
c   xr3egij(ni,ni,idx)  - real energy diff from j->i in ergs,
c                       only computed locally at startup, use line
c                       map xr3lines_egij to look up eij from lineid
c
c   xr3gfij(ni,ni,idx) - real abs osc str for lines in cascades, inc gi
c                        computed at startup, maps to line list
c                        below, also computed at startup
c
c   nxr3ioncol(idx) - number of collision splines, excitations
c      xr3col_i   (icol,idx)  - int starting level, usually 1
c      xr3col_j   (icol,idx)  - int target (upper) level, <= ni
c      xr3col_fr  (icol,idx)  - fraction of ground pop in target
c      xr3col_tc  (icol,idx) - real spline scale temperature
c      xr3col_nspl(icol,idx)    - int number of spline points, default 9
c      xr3col_typespl(icol,idx) - int spline type, default 2
c      xr3col_x (nspl)          - real coll upsilon spline x coord
c      xr3col_y (nspl,icol,idx)   - real coll upsilon spline y coord
c      xr3col_y2(nspl,icol,idx)   - real coll upsilon spline y2 coord
c                                 up to 9 type 3 BT92 spline points
c      xr3col_nl(icol,idx)      - int number of branching lines
c      xr3col_jj(line,icol,idx) - int upper level of branch
c      xr3col_kk(line,icol,idx) - int lower level of branch
c      xr3col_br(line,icol,idx) - real fraction of coll *rate* in line
c                                 multiply by Ejj,kk to get line bright
c
c multiple lines are summed in a local ni x ni array, and then mapped to
c the main line list, and stored there like other line lists
c
c   nxr3lines                   - int, total line count for 1D line list
c                                  computed at startup
c   xr3lines_map(i,j,idx)       - int, line list index (lidx)
c                                 for transition i -> j in ion idx.
c                                 computed at startup, only from lines
c                                 with branching ratios
c   xr3lines_imapidx(nt)        - inverse map from line to idx
c   xr3lines_imap_i(nt)         - inverse map from line to lower level
c   xr3lines_imap_j(nt)         - inverse map from line to upper level
c
c   xr3lines_at  (lidx)         - int mappings atom id
c   xr3lines_ion (lidx)         - int ion stage 1...
c   xr3lines_tran(lidx)         - char*24 transition id
c
c   xr3lines_eij (lidx)         - real energy of line, cm^-1
c   xr3lines_egij(lidx)         - real energy of line, ergs
c   xr3lines_evij(lidx)         - real energy of line, eV , not used
c   xr3minxbin                  - int, first bin to look in
c   xr3lines_bin (lidx)         - int photdat spectrum bin index that
c                                 contains the line
c
c   xr3lines_lam (lidx)         - real wavelength A (vac Ritz)
c   xr3lines_aji (lidx)         - real Trans Probability
c   xr3lines_frac(lidx)         - >0 = dismul population fraction
c   xr3lines_gf  (lidx)         - real abs osc str, incl stat weight
c                                  computed at startup
c   xr3lines_bri (lidx)         - real brightness - computed as we go
c   xr3lines_flux(lidx)         - real integrated flux
c
c   xr3lines_emilin(2,lidx)     - real line emission in phots, for trans
c   xr3lines_uplin (2,lidx)     - real line upstream in phots, for trans
c   xr3lines_dwlin (2,lidx)     - real line downstr in phots, for trans
c
c Li sequence, idx = 20, ni = 15, icol = 14 max coll strengths
c , line = 12 max cascade lines, constants in const.inc
c other sequences to add later RSS2015
c
c
c V3 Xray Resonance lines data
c
      real*8 xr3loss
c
c multi level ion data
c
      integer*4 nxr3ions
      integer*4 xr3at(mxxr3ions),xr3ion(mxxr3ions),xr3ni(mxxr3ions)
      integer*4 nxr3ioncol(mxxr3ions),xr32s12j(mxxr3ions)
      integer*4 xr3col_i(mxxr3cols,mxxr3ions)
      integer*4 xr3col_j(mxxr3cols,mxxr3ions)
      integer*4 xr3col_nspl(mxxr3cols,mxxr3ions)
      integer*4 xr3col_typespl(mxxr3cols,mxxr3ions)
      integer*4 xr3col_nl(mxxr3cols,mxxr3ions)
      integer*4 xr3col_jj(mxxr3casc,mxxr3cols,mxxr3ions)
      integer*4 xr3col_kk(mxxr3casc,mxxr3cols,mxxr3ions)
      real*8 xr3col_fr(mxxr3cols,mxxr3ions)
      real*8 xr3col_tc(mxxr3cols,mxxr3ions)
      real*8 xr3col_x (mxxr3nspl) ! uniform 11 point splines for all
      real*8 xr3col_y (mxxr3nspl,mxxr3cols,mxxr3ions)
      real*8 xr3col_y2(mxxr3nspl,mxxr3cols,mxxr3ions)
      real*8 xr3col_br(mxxr3casc,mxxr3cols,mxxr3ions)
      real*8 xr3ei  (mxxr3lvls,mxxr3ions)
      real*8 xr3gi  (mxxr3lvls,mxxr3ions)
      real*8 xr3invgi(mxxr3lvls,mxxr3ions)
c
c      real*8 xr3eij (mxxr3lvls,mxxr3lvls,mxxr3ions) !now local to init
c      real*8 xr3egij(mxxr3lvls,mxxr3lvls,mxxr3ions) !now local to init
c      real*8 xr3aji (mxxr3lvls,mxxr3lvls,mxxr3ions) !now local to init
c      real*8 xr3gfij(mxxr3lvls,mxxr3lvls,mxxr3ions) !now local to init
c      real*8 xr3frij(mxxr3lvls,mxxr3lvls,mxxr3ions) !now local to init
c
c line list data
c
      integer*4 nxr3lines
c
      integer*4 xr3lines_at  (mxxr3lines)
      integer*4 xr3lines_ion (mxxr3lines)
c
      integer*4 xr3lines_map(mxxr3lvls,mxxr3lvls,mxxr3ions)
      integer*4 xr3lines_imapidx(mxxr3lines)
      integer*4 xr3lines_imap_i(mxxr3lines)
      integer*4 xr3lines_imap_j(mxxr3lines)
c
c      real*8 xr3lines_eij (mxxr3lines) not used
      real*8 xr3lines_egij(mxxr3lines) ! ergs
c      real*8 xr3lines_evij(mxxr3lines) not used
      integer*4 xr3minxbin
      integer*4 xr3lines_bin (mxxr3lines)
c
      real*8 xr3lines_lam (mxxr3lines)
      real*8 xr3lines_aji (mxxr3lines)
      real*8 xr3lines_frac(mxxr3lines)
c
      real*8 xr3lines_gf  (mxxr3lines)
c
      real*8 xr3lines_bri (mxxr3lines)
      real*8 xr3lines_flux(mxxr3lines)
c
      real*8 xr3lines_emilin(2,mxxr3lines)
      real*8 xr3lines_uplin (2,mxxr3lines)
      real*8 xr3lines_dwlin (2,mxxr3lines)
c
      character*16 xr3id(mxxr3lvls,mxxr3ions)
      character*32 xr3lines_tran(mxxr3lines)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c V3 Xray Resonance lines data for large ions, same as xr3 data but for
c fewer, larger atom models  xr3 -> xrl
c
      real*8 xrlloss
c
c multi level ion data
c
      integer*4 nxrlions
      integer*4 xrlat(mxxrlions),xrlion(mxxrlions),xrlni(mxxrlions)
      integer*4 nxrlioncol(mxxrlions)
      integer*4 xrlcol_i(mxxrlcols,mxxrlions)
      integer*4 xrlcol_j(mxxrlcols,mxxrlions)
      integer*4 xrlcol_nspl(mxxrlcols,mxxrlions)
      integer*4 xrlcol_typespl(mxxrlcols,mxxrlions)
c
      integer*4 xrlcol_nl(mxxrlcols,mxxrlions)
      integer*4 xrlcol_jj(mxxrlcasc,mxxrlcols,mxxrlions)
      integer*4 xrlcol_kk(mxxrlcasc,mxxrlcols,mxxrlions)
c
      real*8 xrlcol_fr(mxxrlcols,mxxrlions)
      real*8 xrlcol_tc(mxxrlcols,mxxrlions)
c
      real*8 xrlcol_x (mxxrlnspl) ! uniform 11 point splines for all
      real*8 xrlcol_y (mxxrlnspl,mxxrlcols,mxxrlions)
      real*8 xrlcol_y2(mxxrlnspl,mxxrlcols,mxxrlions)
c
      real*8 xrlcol_br(mxxrlcasc,mxxrlcols,mxxrlions)
c
      real*8 xrlei  (mxxrllvls,mxxrlions)
      real*8 xrlgi  (mxxrllvls,mxxrlions)
      real*8 xrlinvgi(mxxrllvls,mxxrlions)
c
c      real*8 xrleij (mxxrllvls,mxxrllvls,mxxrlions) !now local to init
c      real*8 xrlegij(mxxrllvls,mxxrllvls,mxxrlions) !now local to init
c      real*8 xrlaji (mxxrllvls,mxxrllvls,mxxrlions) !now local to init
c      real*8 xrlgfij(mxxrllvls,mxxrllvls,mxxrlions) !now local to init
c      real*8 xrlfrij(mxxrllvls,mxxrllvls,mxxrlions) !now local to init
c
c line list data
c
      integer*4 nxrllines
c
      integer*4 xrllines_at  (mxxrllines)
      integer*4 xrllines_ion (mxxrllines)
c
      integer*4 xrllines_map(mxxrllvls,mxxrllvls,mxxrlions)
      integer*4 xrllines_imapidx(mxxrllines)
      integer*4 xrllines_imap_i(mxxrllines)
      integer*4 xrllines_imap_j(mxxrllines)
c
c      real*8 xrllines_eij (mxxrllines) not used
      real*8 xrllines_egij(mxxrllines) ! ergs
c      real*8 xrllines_evij(mxxrllines) not used
      integer*4 xrlminxbin
      integer*4 xrllines_bin (mxxrllines)
c
      real*8 xrllines_lam (mxxrllines)
      real*8 xrllines_aji (mxxrllines)
      real*8 xrllines_frac(mxxrllines)
c
      real*8 xrllines_gf  (mxxrllines)
c
      real*8 xrllines_bri (mxxrllines)
      real*8 xrllines_flux(mxxrllines)
c
      real*8 xrllines_emilin(2,mxxrllines)
      real*8 xrllines_uplin (2,mxxrllines)
      real*8 xrllines_dwlin (2,mxxrllines)
c
      character*16 xrlid(mxxrllvls,mxxrlions)
      character*32 xrllines_tran(mxxrllines)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c V3 HI sequence Resonance collisional data
c
      real*8 hheloss
c
c multi level ion data
c
      integer*4 nhheions
      integer*4 ishheion(mxelem) ! 0: not a full coll ion, 1: full data
      integer*4 hheat(mxhheions),hheion(mxhheions)
      integer*4 hheni(mxhheions) ! # of levels index
      integer*4 hhen (mxhhelvls,mxhheions) ! index -> quantum n for h
      integer*4 nhheioncol(mxhheions)
      integer*4 hhecol_i(mxhhecols,mxhheions)
      integer*4 hhecol_j(mxhhecols,mxhheions)
      integer*4 hhecol_nspl(mxhhecols,mxhheions)
      integer*4 hhecol_typespl(mxhhecols,mxhheions)
c      integer*4 hhecol_nl(mxhcols,mxhheions)
c      integer*4 hhecol_jj(mxxr3casc,mxhcols,mxhheions)
c      integer*4 hhecol_kk(mxxr3casc,mxhcols,mxhheions)
      real*8 hhecol_tc(mxhhecols,mxhheions)
      real*8 hhecol_x (mxhhenspl)
      real*8 hhecol_y (mxhhenspl,mxhhecols,mxhheions)
      real*8 hhecol_y2(mxhhenspl,mxhhecols,mxhheions)
c      real*8 hhecol_br(mxxr3casc,mxxr3cols,mxhheions)
      real*8 hheei   (mxhhelvls,mxhheions)
      real*8 hhegi   (mxhhelvls,mxhheions)
      real*8 hheinvgi(mxhhelvls,mxhheions)
c
c      real*8 hheeij (mxhhelvls,mxhhelvls,mxhheions) ! cm^-1
c      real*8 hheegij(mxhhelvls,mxhhelvls,mxhheions) ! ergs
c      real*8 hheaji (mxhhelvls,mxhhelvls,mxhheions)
c      real*8 hhegfij(mxhhelvls,mxhhelvls,mxhheions)
c
      character*16 hheid(mxhhelvls,mxhheions)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare RLINE block vaiables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   resonance (permitted) line data
c   rlam: wavelengths (units?)
c   rbri: line brightnesses (ergs)
c   rloss: resonance line heat losses
c   lrbin : location in diffuse field bins for each line
c    rlin(2,60),ruplin(2,60),rdwlin(2,60) lines with distance
c    multipliers like emilin etc
c    fluxr: integrated fluxes
c
      real*8 rlam(60),rbri(60),rloss
      integer*4 lrbin(60)
      integer*4 nlines,ielr(60),ionr(60)
      real*8 rlin(2,60),ruplin(2,60),rdwlin(2,60)
      real*8 fluxr(60)
c
c   resonance line strengths data
c
      real*8 e12r(60),fabsr(60)
      real*8 omres(60),tkexres(60)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare CLINE block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   semi-forbidden intercombination line information
c
c   fslam: wavelengths
c   fsbri: brightnesses
c   fslos :semi forbidden cooling rate (ergs)
c   intam: ion/atom identification
c       (will be changed to numerical arrays
c       and use the global ion & roman numeral character arrays)
c   lcbin : location in diffuse field for each line
c
      real*8 fslam(80),fsbri(80),fslos
      integer*4 lcbin(80)
      integer*4 mlines,ielfs(80),ionfs(80)
      character intam(2,80)*4
      real*8 fslin(2,80),fsuplin(2,80),fsdwlin(2,80)
      real*8 fluxi(80)
c
c      semi-forbidden line strengths data
c
      real*8 e12fs(80),w1fs(80)
      real*8 w2fs(80),a21fs(80)
      real*8 omfs(80),tkexfs(80)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare F3LINE block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   three level atom line information
c   f3lam: wavelengths
c   f3bri: brightnesses
c   f3loss: forbidden cooling rate (ergs)
c   f3ion,f3atom: forbiden line ion identification,
c   read from threedat
c   f3type: type of three level lines,1 = finestructure
c   2 = forbidden + others
c
      real*8 f3lam(3,60),f3bri(3,60),f3loss
      integer*4 f3ion(60),f3atom(60),nf3ions,nf3trans
      integer*4 f3type(60),f3bin(3,60)
      real*8 f3lin(2,3,60),f3uplin(2,3,60),f3dwlin(2,3,60)
      real*8 fluxf3(3,60)
c
c   three level line strengths data
c
      real*8 wi3(3,60),ei3(3,60),tkex3(3,60)
      real*8 ai3(3,60),gam3(2,3,60)
      real*8 pow3(2,3,60)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare FMLINE multi-level atom block variables
c      for 5 - 9 level atoms,10 - 36 transitions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Multi Level line information
c   fmlam: wavelengths
c   fmeij: Delta Eij in transition list form
c   fmbri: brightnesses
c   fmx  : level populations
c   fmloss: forbidden cooling rate (ergs)
c   fmion,fmatom,fmnl: fmspindex line ion identifications and num levels
c   read from MULTIDAT
c   fluxm: integrated lines for sline block
c
      integer*4 nfmions
      integer*4 nfmtrans(mxmspecs)
      integer*4 fmion(mxmspecs)
      integer*4 fmatom(mxmspecs)
      integer*4 fmspindex(mxion,mxelem)
      integer*4 fmnl(mxmspecs)
      integer*4 fmkappaidx(mxmspecs)
      integer*4 fmbin(mxtr,mxmspecs)
      integer*4 nfmtridx(mxnl,mxnl,mxmspecs)
      integer*4 nfmlower(mxtr,mxmspecs)
      integer*4 nfmupper(mxtr,mxmspecs)
      integer*4 fmomdatatype(mxtr,mxmspecs)
      integer*4 fmombtn     (mxtr,mxmspecs)
      character*18 nfmid(mxtr,mxmspecs)
      character*8  nfmterm(mxnl,mxmspecs)
c
      real*8 fmloss
      real*8 fmlam   (mxtr,mxmspecs)
      real*8 fmeij   (mxtr,mxmspecs)
      real*8 fmbri   (mxtr,mxmspecs)
      real*8 fluxm   (mxtr,mxmspecs)
c
      real*8 fmombtc   (mxtr,mxmspecs)
      real*8 fmombsplx (mxnspl,mxtr,mxmspecs)
      real*8 fmombsply (mxnspl,mxtr,mxmspecs)
      real*8 fmombsply2(mxnspl,mxtr,mxmspecs)
c
      integer*4 nkappaions
      integer*4 nkappas
      integer*4 kapatom(mxmspecs)
      integer*4 kapion(mxmspecs)
      integer*4 kapnl(mxmspecs)
      integer*4 nkaptrans(mxmspecs)
      integer*4 nkaptridx(mxnl,mxnl,mxmspecs)
      integer*4 kupstype (mxtr,mxmspecs)
      integer*4 kupsbtn  (mxtr,mxmspecs)
      integer*4 nkaplower(mxtr,mxmspecs)
      integer*4 nkapupper(mxtr,mxmspecs)
c
      real*8 kappas    (mxkappas)
      real*8 kupsbtc   (mxtr,mxmspecs)
      real*8 kupssplx  (mxnspl,mxtr,mxmspecs)
      real*8 kupssply  (mxnspl,mxkappas,mxtr,mxmspecs)
      real*8 kupssply2 (mxnspl,mxkappas,mxtr,mxmspecs)
c
      real*8 fmx(mxnl,mxmspecs)
c
c line distance mulitpliers for strong lines
c
      real*8 fmlin  (2,mxtr,mxmspecs)
      real*8 fmuplin(2,mxtr,mxmspecs)
      real*8 fmdwlin(2,mxtr,mxmspecs)
c
c   Multi-Level line strengths data
c
      real*8 wim(mxnl,mxmspecs)
      real*8 eim(mxnl,mxnl,mxmspecs)
      real*8 aim(mxnl,mxnl,mxmspecs)
      real*8 omim(mxnl,mxnl,mxmspecs)
      real*8 tdepm(mxnl,mxnl,mxmspecs)
      real*8 tkexm(mxnl,mxnl,mxmspecs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc      declare FELINE block variables
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   Iron  multiplet line information, initially for FeII, could extend
cc   for transition species in general, hence featom used
cc
cc
c
c   Multi Level line information
c   felam: wavelengths
c   feeij: Delta Eij in transition list form
c   febri: brightnesses
c   fex  : level populations
c   feloss: forbidden cooling rate (ergs)
c   feion,featom,fenlevels:  line ion identification and num levels
c   read from FEIIDAT
c   fluxfe: integrated lines for sline block
c
      integer*4 nfeions
      integer*4 nfetrans(mxfespecs)
      integer*4 feion(mxfespecs)
      integer*4 featom(mxfespecs)
      integer*4 fespindex(mxion,mxelem)
      integer*4 fenl(mxfespecs)
      integer*4 febin(mxfetr,mxfespecs)
      integer*4 nfetridx(mxfenl,mxfenl,mxfespecs)
      integer*4 nfelower(mxfetr,mxfespecs)
      integer*4 nfeupper(mxfetr,mxfespecs)
      integer*4 feomdatatype(mxfetr,mxfespecs)
      integer*4 feombtn   (mxfetr,mxfespecs)
      character*24 nfeid(mxfetr,mxfespecs)
      character*8  nfeterm(mxfenl,mxfespecs)
c
      real*8 feloss
      real*8 felam   (mxfetr,mxfespecs)
      real*8 feeij   (mxfetr,mxfespecs)
      real*8 febri   (mxfetr,mxfespecs)
      real*8 fluxfe   (mxfetr,mxfespecs)
c
      real*8 feombtc   (mxfetr,mxfespecs)
      real*8 feombsplx (mxfenspl,mxfetr,mxfespecs)
      real*8 feombsply (mxfenspl,mxfetr,mxfespecs)
      real*8 feombsply2(mxfenspl,mxfetr,mxfespecs)
c
c line distance multipliers for strong lines
c
      real*8 felin  (2,mxfetr,mxfespecs)
      real*8 feuplin(2,mxfetr,mxfespecs)
      real*8 fedwlin(2,mxfetr,mxfespecs)
c
c   Multi-Level line strengths data
c
      real*8 fex(mxfenl,mxfespecs)
c
      real*8 wife(mxfenl,mxfespecs)
      real*8 eife(mxfenl,mxfenl,mxfespecs)
      real*8 aife(mxfenl,mxfenl,mxfespecs)
      real*8 omife(mxfenl,mxfenl,mxfespecs)
      real*8 tdepfe(mxfenl,mxfenl,mxfespecs)
      real*8 tkexfe(mxfenl,mxfenl,mxfespecs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Hydrogen data.....
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   New for hydrec2:
c
c   hspline: cubic spline fits to accurate hydrogenic recombination
c    rate for total: hspline(1,x,coeffs),sum of 210 levels:
c    hspline(2,x,coeffs),and n=1 ,l = 0: hspline(3,x,coeffs)
c    x = log(t),result = log(rate)
c
c    Nahar total recombination new NR style spline computed at startup
c    from the original data points.
c    Hydrogenic
c    nhrec
c    hreclgt,           log10 T/z^2
c    hreclga,hreclga2    log10 total recom coeff
c    hreclga1,hreclga12  log10 n=1 recom coeff
c
c    HeliumI
c    nehrec
c    hereclgt,           log10 T/z^2
c    hereclga,hereclga2    log10 total recom coeff
c    hereclga1,hereclga12  log10 n=1 recom coeff
c
c   hcoln: n -> n' collisional rate fits
c   dnl: level term for nl->nl' collisions
c
c   Line series:
c
c   hlambda(line,series) ,ie Hbeta = hlambda(2,2) (Balmer line 2)
c
c
c    th42(12),the43(12)   temperatures for hydrogen & Helium spec calc
c      data points
c    rhhe(9)  densities (n) for  hydrogen & Helium spec calc data
c      points (same for h and he)
c    h42a(12,9),h42b(12,9) lookup tables for hydrogen hbeta emission
c    he43a(12,9),he43b(12,9) lookup tables for helium 4686 emission
c    Use powerlaw interpolations and extrapolations
c    hylratsa(line,series,temps,rho) id = same order as lambdas:
c      ie lyman series then balmer etc
c    hylratsb(line,series,temps,rho)      ten lines each for six series
c    helratsa(line,series,temps,rho)
c    helratsb(line,series,temps,rho)
c
c    spline coeffs in temperature direction:
c
c    hylratsa2(line,series,temps,rho)
c    hylratsb2(line,series,temps,rho)
c    helratsa2(line,series,temps,rho)
c    helratsb2(line,series,temps,rho)
c
c    collrate2p(mxelem): collision rate to 2S level of H-like atoms
c    collrate2phr(mxelem): collision rate to 2S level of He-like atoms
c       from ground only
c    rateton(mxelem,n): total collision rate to level n of H-like atoms
c
c
      real*8 balcoe(5,16),hbema(5,6),hbemb(5,6),hherat(14)
      real*8 hheen(14),hlam(10),hbri(10),h2ql
      real*8 h2qav,hei2qa,heii2qa,fluxh(10),fhbeta
c
      real*8 collrate2p(mxelem),ee2p(mxelem)
      real*8 collrate2phe(mxelem),ee2phe(mxelem)
      real*8 recrate2p(mxelem),h2qflux(mxelem),h2qbri(mxelem)
      real*8 r2s1a(10,9),r2s1b(10,9),rateton(mxelem,16)
      real*8 r2s2a(12,9),r2s2b(12,9)
c
c      real*8 hspline(3,18,5)
c
      real*8 hloss,caseab(mxelem)
c
      integer*4 nhrec
      real*8 hreclgt(121)
      real*8 hreclga(121),hreclga2(121)
      real*8 hreclga1(121),hreclga12(121)
c
      integer*4 nherec
      real*8 hereclgt(81)
      real*8 hereclga(81),hereclga2(81)
      real*8 hereclga1(81),hereclga12(81)
c
      real*8 hydrogf(mxnhlines,mxnhs),rhhe(9)
c
      real*8 hbeta,halpha,lyalpha
      integer*4 nhlines,nhseries,hbin(mxnhlines,mxnhs)
      real*8 hlambda(mxnhlines,mxnhs)
      real*8 hydrobri(mxnhlines,mxnhs)
      real*8 hydlin(2,mxnhlines,mxnhs)
      real*8 hyduplin(2,mxnhlines,mxnhs)
      real*8 hyddwlin(2,mxnhlines,mxnhs)
c
      real*8 th42(12)
      real*8 h42a (12,9),h42b (12,9)
      real*8 h42a2(12,9),h42b2(12,9)
      real*8 hylratsa (mxnhlines,mxnhs,10,9)
      real*8 hylratsb (mxnhlines,mxnhs,10,9)
      real*8 hylratsa2(mxnhlines,mxnhs,10,9)
      real*8 hylratsb2(mxnhlines,mxnhs,10,9)
      real*8 hydroflux(mxnhlines,mxnhs)
c
      integer*4 nhelines,nheseries,hebin(mxnhelines,mxnhes)
      real*8 helambda(mxnhelines,mxnhes)
      real*8 helibri(mxnhelines,mxnhes)
      real*8 hellin(2,mxnhelines,mxnhes)
      real*8 heluplin(2,mxnhelines,mxnhes)
      real*8 heldwlin(2,mxnhelines,mxnhes)
c
      real*8 the43(12)
      real*8 he43a(12,9),he43b(12,9)
      real*8 he43a2(12,9),he43b2(12,9)
      real*8 helratsa(mxnhelines,mxnhes,12,9)
      real*8 helratsb(mxnhelines,mxnhes,12,9)
      real*8 helratsa2(mxnhelines,mxnhes,12,9)
      real*8 helratsb2(mxnhelines,mxnhes,12,9)
      real*8 heliflux(mxnhelines,mxnhes)
c
      real*8 heiioiiibf,heiioiiibfsum
c
      integer*4 nxhlines,nxhseries,xhbin(mxnxhlines,mxnxhs,mxelem)
      real*8 xhbeta(mxelem)
      real*8 xhlambda(mxnxhlines,mxnxhs,mxelem)
      real*8 xhydrobri(mxnxhlines,mxnxhs,mxelem)
c
      real*8 xhydroflux(mxnxhlines,2,mxelem)
      real*8 xhydlin(2,mxnxhlines,2,mxelem)
      real*8 xhyduplin(2,mxnxhlines,2,mxelem)
      real*8 xhyddwlin(2,mxnxhlines,2,mxelem)
c
c
c   Collisional Excitation
c
c   ccoln(56,4),dcoln(56,4): collision strength fits for hydrogenic
c     lines
c   colid(56,4): collisional level ids for ccoln and dcoln
c   collha(line,series,coll.level),collhb(line,series,coll.level):
c     Contribution to line in series from collisions to coll.level
c     for case a and b,assuming 'n method'
c
c
      real*8 ccoln(56,4),dcoln(56,4)
      integer*4 colid(56,4)
      real*8 collha(mxnhlines,mxnhs,15),collhb(mxnhlines,mxnhs,15)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Helium data......
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 heii2ql,heien(4),qren(4)
      real*8 hein(4),heiilr,hel0(6,12)
      real*8 heilam(16),heibri(16),fluxhei(16)
      integer*4  heibin(16),nheilines
c
c new He I effective recomb coeffs plus coll enhancement
c
      integer*4 nheis, nheit  ! number of level in singlet and triplets
      integer*4 nheiscr(14)  ! number of coeffs for collisions
      integer*4 nheitcr(14)  ! number of coeffs for collisions
      integer*4 nheislines, nheitlines ! number of lines computed
      integer*4 heisupid(20), heisloid(20)
      integer*4 heitupid(20), heitloid(20)
      integer*4 heisbin(20), heitbin(20)
c
      real*8 heisej(14),heitej(14)    ! levels
      real*8 heiseij(20),heiteij(20) ! lines
      real*8 heislam(20),heitlam(20) ! lines
      real*8 heisbri(20),fluxheis(20)
      real*8 heitbri(20),fluxheit(20)
c
      real*8 heisreccoef(4,20)    ! singlet line coefs per line
      real*8 heiscra(6,14)        ! singlet coll coeff per level
      real*8 heiscrb(6,14)        ! singlet coll coeff per level
      real*8 heiscrc(6,14)        ! singlet coll coeff per level
c
      real*8 heitreccoef(4,20)    ! triplet line coefs per line
      real*8 heitcra(6,14)        ! triplet coll coeff per level
      real*8 heitcrb(6,14)        ! triplet coll coeff per level
      real*8 heitcrc(6,14)        ! triplet coll coeff per level
c
      character*8 heislvl(20),heitlvl(20) ! level names
      character*24 heislamid(20),heitlamid(20) ! transition names
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Heavy Element Recombination
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c CIII -> CII recombination, CII spectrum of Doublets
c The 16 best multiplets with selected individual red lines
c
      integer*4 nrccii
      integer*4 rccii_bin(20)
      real*8 rccii_eij(20)
      real*8 rccii_lam(20)
      real*8 rccii_logte(15)        ! log_10 te for rates
      !case1=a,2=b,coef0-13,lineid
      real*8 rccii_coeffy (2,15,20) ! log_10 coeffs from tables
      real*8 rccii_coeffy2(2,15,20) ! 2nd drv of log_10 cs from tables
      character*24 rccii_tr(20)
c
c lines and fluxes
c case A - B
      real*8 rccii_abri(20),fluxrccii_a(20)
      real*8 rccii_bbri(20),fluxrccii_b(20)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c OII -> OI recombination, OI spectrum of Triplets and Quintets
c
c Data
c
      integer*4 nrcoi_q
      integer*4 rcoi_qbin(36)
      real*8 rcoi_qeij(36)
      real*8 rcoi_qlam(36)
      real*8 rcoi_qcoef(2,3,36)  ! case 1=a, 2=b, coef, id
      character*16 rcoi_qtr(36)
      integer*4 nrcoi_t
      integer*4 rcoi_tbin(30)
      real*8 rcoi_teij(30)
      real*8 rcoi_tlam(30)
      real*8 rcoi_tcoef(2,3,30)  ! case, coef, id
      character*24 rcoi_ttr(30)
c
c lines and fluxes
c case A
      real*8 rcoi_qabri(36),fluxrcoi_qa(36)
      real*8 rcoi_tabri(30),fluxrcoi_ta(30)
c case B
      real*8 rcoi_qbbri(36),fluxrcoi_qb(36)
      real*8 rcoi_tbbri(30),fluxrcoi_tb(30)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c OIII -> OII recombination, OII spectrum of Doublets Quartets
c and Intermediates.
c
c The 35 multiplets including the case insensitive 4f intermediate
c coupling levels, and the 3d D and P quartets and doublets, for
c a total of 200 lines. May edit back to smaller list later
c
      integer*4 nrcoii
      integer*4 rcoii_bin(200)
      integer*4 rcoii_mltid(200)
      real*8 rcoii_eij(200)
      real*8 rcoii_lam(200)
      real*8 rcoiin20_y (3,8,200) ! ne=1e2 case 1=a, 2=b, data te, id
      real*8 rcoiin40_y (3,8,200) ! ne=1e4 case 1=a, 2=b, data te, id
      real*8 rcoiin50_y (3,8,200) ! ne=1e5 case 1=a, 2=b, data te, id
      real*8 rcoiin60_y (3,8,200) ! ne=1e6 case 1=a, 2=b, data te, id
      real*8 rcoiin20_y2(3,8,200) ! ne=1e2 case 1=a, 2=b, spline y2, id
      real*8 rcoiin40_y2(3,8,200) ! ne=1e4 case 1=a, 2=b, spline y2, id
      real*8 rcoiin50_y2(3,8,200) ! ne=1e5 case 1=a, 2=b, spline y2, id
      real*8 rcoiin60_y2(3,8,200) ! ne=1e6 case 1=a, 2=b, spline y2, id
      character*24 rcoii_tr(200)
c
c lines and fluxes
c case A - B - C
c
      real*8 rcoii_abri(200),fluxrcoii_a(200)
      real*8 rcoii_bbri(200),fluxrcoii_b(200)
      real*8 rcoii_cbri(200),fluxrcoii_c(200)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NIII -> NII recombination, NII spectrum of Triplets and Intermediates
c The 8 best multiplets with individual lines, with density dependence
c Data fits as splines
c
      integer*4 nrcnii
      integer*4 rcnii_bin(55)
      real*8 rcnii_eij(55)
      real*8 rcnii_lam(55)
      real*8 rcniin20_y (2,8,55) ! ne=1e2 case 1=a, 2=b, data te, id
      real*8 rcniin30_y (2,8,55) ! ne=1e3 case 1=a, 2=b, data te, id
      real*8 rcniin40_y (2,8,55) ! ne=1e4 case 1=a, 2=b, data te, id
      real*8 rcniin50_y (2,8,55) ! ne=1e5 case 1=a, 2=b, data te, id
      real*8 rcniin20_y2(2,8,55) ! ne=1e2 case 1=a, 2=b, spline y2, id
      real*8 rcniin30_y2(2,8,55) ! ne=1e3 case 1=a, 2=b, spline y2, id
      real*8 rcniin40_y2(2,8,55) ! ne=1e4 case 1=a, 2=b, spline y2, id
      real*8 rcniin50_y2(2,8,55) ! ne=1e5 case 1=a, 2=b, spline y2, id
      character*24 rcnii_tr(55)
c
c lines and fluxes
c case A - B
      real*8 rcnii_abri(55),fluxrcnii_a(55)
      real*8 rcnii_bbri(55),fluxrcnii_b(55)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NeIII -> NeII recombination, NeII spectrum of Doublets and Quartets
c The 6 best multiplets with individual lines via branching ratios
c Data
c
      integer*4 nrcneii
      integer*4 rcneii_bin(38)
      real*8 rcneii_eij(38)
      real*8 rcneii_lam(38)
      real*8 rcneii_y (2,8,38)  ! case 1=a, 2=b, spline, id
      real*8 rcneii_y2(2,8,38)  ! case 1=a, 2=b, spline, id
      character*24 rcneii_tr(38)
c
c lines and fluxes
c case A - B
      real*8 rcneii_abri(38),fluxrcneii_a(38)
      real*8 rcneii_bbri(38),fluxrcneii_b(38)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Continuum (free-free,f-b and two photon) emmision data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   zn0 = effective charge of ion
c   ezz = 2 photon 1s-2s energy gaps (from grotian diags) (eV)
c   fn0 = quantum nunber of the ground state
c   zt0 = # of vacancies available
c   ffph,fbph,p2ph: arrays of continuum component photons
c   cnphot: total continuum array (photons),uses photev bins
c
c   gff(41,81),gffy2(41,81),gffg2(41),gffu(81) -free free gaunt
c                                               factor lookup table
c   gfb(20,50),gfbe(50),gfbtz(20) -free bound gaunt factor lookup table
c
c   gffintspl: (41,5) cubic spline fits to accurate integrated
c     freefree emission gfactor
c   x = log(gamma^2),result = <gff(g^2)>
c   coeffs: 1 = log g2,2 = a0,3 = a1,4=a2,5 = a3
c
      real*8 zn0(mxion,mxelem),ezz(mxion,mxelem)
      integer*4 fn0(mxion,mxelem),zt0(mxion,mxelem)
      real*8 ffph(mxinfph),fbph(mxinfph),p2ph(mxinfph)
      real*8 cnphot(mxinfph)
c
      real*8 gff(41,81),gffy2(41,81),gffg2(41),gffu(81)
      real*8 gffintspl(41,5)
      integer*4 ngffg2,ngffu,ngffint
c      integer*4 ngfbe,ngfbtz
      real*8 gffintg2(41),gffinty(41),gffinty2(41)
c
c      real*8 gfb(20,50),gfbe(50),gfbtz(20)
c
c New Hydrogenic FB data
c
      integer*4 nfbe
      real*8 fbloggamma(81,98)
      real*8 fblte(81),fbenu(98),fben(98)
      real*8 fbee(mxelem,mxinfph)
      integer*4 nfbenu(98)
      integer*4 jfb(mxelem,mxinfph)
c
c
c New Heliogenic FB data
c
      integer*4 nhefbe
      real*8 hefbloggamma(76,123)
      real*8 hefblte(76),hefbenu(123),hefben(123)
      real*8 hefbee(mxelem,mxinfph)
      integer*4 nhefbenu(123)
      integer*4 jhefb(mxelem,mxinfph)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Summed (integerated) model quantities
c
c      declare SLINE block variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 deav,teav(mxion,mxelem),rdisa(mxion,mxelem),
     &           pam(mxion,mxelem)
c
c      declare TONLIN block variables
c
      real*8 vunilog
c
c      declare DENLIN block variables
c
      real*8 roiii,rnii,weoiii
      real*8 wenii,deoiii,denii,toiii,tnii
      real*8 deam(mxion,mxelem),dhav,roii
      real*8 rsii,weoii,wesii,toii,tsii
      real*8 deoii,desii,qtosoh
      real*8 qhdh,qhdn,uhdh,zetae,qhdha,zetaeav
c
c   declare PHOBAT block variables
c
c   ofac,xfac,pfac,vfac,afac,bfac,cfac,efac,ffacm & scalen are
c   density law parameters:
c   f(x) = ofac*[xfac*exp((x-afac)/scalen) + pfac*((x/afac)**bfac)
c                vfac*(((r0-x)/efac)**-ffac)]+cfac
c   tofac,txfac,tpfac,tafac,tbfac,tcfac are the same temperature
c      profile
c   aef: accretion mass-luminosity efficiency approx 0.1
c   mdot: accretion rate
c   frlum: final luminosity fraction after time dependent source turns
c          off
c   radin,radout: inner and outer radii for ring aperture integration
c   ponk,tinner: pressure, intial t for P6/7
c   ponk0,ponk1 : inner and outer pressures in P7
c   rp0,palpha : scale radius and slope for P(r) in P7
c
      real*8 dhavv,viofr,aef
      real*8 pop0(mxion,mxelem),wdpl,telap,tlife,remp,rmax
      real*8 alphahb,alphaheb
      real*8 volstromhb,volstromheb
      real*8 rstromhb,rstromheb
      real*8 fren,tauen,colend,tend,diend,dtau0
      real*8 qhdnin,qhdnav,qhdhin,qhdhav,uhin,uhav,unin,unav
      real*8 tm00,admach,alphaturbfixed
      real*8 ofac,xfac,pfac,afac,bfac,cfac,scalen
      real*8 vfac,efac,ffac
      real*8 tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
      real*8 frlum,radout,radin
      real*8 ponk,tinner
      real*8 ponk0,ponk1,rp0,palpha
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      **SHOCK FRONT DATA :
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare VELSHO block variables
c
      real*8 dhpo,vshoc,vpo,vfin,hmago
      real*8 fm,dens,pres,ww,qst,qfin,depr
c
c
c   Declare flow and stucture variables
c   pr0,1: pressures at point a and b.
c   rho0,1: densities (gcm^-3)
c   te0,1: electron temps
c   vel0,1: flow velocities (cm s-1)
c   dh0,1: Hydrogen number densities
c   de0,1: electron number densities
c   bm0,1: magnetic field (Gauss)
c   pb0,1: magnetic pressure (B^2/8pi)
c
      real*8 pr0,rho0,te0,vel0,dh0,de0,bm0,pb0
      real*8 pr1,rho1,te1,vel1,dh1,de1,bm1,pb1
      real*8 dis0,dis1,disav,rad0,rad1
      real*8 gm,flowconst
c
c      declare SHOARR block variables
c
c   arrays to store shock models,mxnsteps
c
      real*8 depo
      real*8 dist(mxnsteps),timlps(mxnsteps)
      real*8 veloc(mxnsteps),xh(mxnsteps)
      real*8 te(mxnsteps),deel(mxnsteps),dhy(mxnsteps)
      real*8 bmg(mxnsteps)
      real*8 hoiii(mxnsteps),hoii(mxnsteps),hoi(mxnsteps)
      real*8 hnii(mxnsteps),hsii(mxnsteps)
      real*8 distim(2,0:mxnsteps)
c      real*8 abn(mxelem,mxnsteps)

      real*8 popi(mxion,mxelem)
      real*8 popf(mxion,mxelem)
      real*8 popini(mxion,mxelem)
      real*8 propop(mxion,mxelem)
      real*8 cav(5),diluf,xhmin,fval,texi
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Global operating mode flags
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   freefreemode: 0 or 1 : simple or advanced free free
c                          continuum (def 0)
c   collmode: 0, 1 (2): collisional ionisation calculation mode(def 0)
c   dielcmode: 0 or 1 :dielectr. recomb. calculation mode(def 0)
c   photonmode: 0 or 1: no coupling=0,normal coupling = 1
c        (should be 1)
c   chargemode: 0 or 1: new charex = 0,oldcharex = 1 ,off =2 (default 0)
c   chargeheatmode 1 = heating on, 0 = heating off (default 1)
c   grainmode: 0 or 1: nograins = 0,include grains = 1 (default 0)
c   graindestructmode: 0 or 1,allow grain to form or destruct
c                      at run time,default 0
c   grainadjust : 0 or 1,allow grain mixture to be recalculated
c                 in adjust grains
c   pahmode: 1 or 0: if PAH molecule physics is being uesd,default 0
c   pahactive: 1 or 0: if PAH molecules have formed 1,
c                      if not,0 default (0)
c   IRmode: 0,1 or 2: nograins = 0,if grainsand not want IR or
c        dusttemp=0,want IR =1,want IR and dusttempdist=2
c        (default = 0)
c   hhecollmode: 0 New H/He coll or 1: Palla 1989 Hcoll excite for n=1-5
c        (default 0)
c   linecoolmode: 0 normal total cooling or 1: only collisional line
c                 cooling, for testing with CHIANTI (default 0)
c   alphacoolmode: 0 or 1: physical =0,or powerlaw=1 cooling
c        (default 0)
c   expertmode: 0 or greater,expert diagnostic mode (default 0)
c
      integer*4 expertmode,collmode,dielcmode,photonmode
      integer*4 chargemode,chargeheatmode
      integer*4 alphacoolmode,hhecollmode,linecoolmode,radpressmode
      integer*4 irmode,turbheatmode,freefreemode
      integer*4 grainmode,pahmode,pahactive,graindestructmode
      integer*4 grainadjust
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      ** General CHARACTER  and INTEGER type common blocks
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare JBLOCK block variables
c
c   iphom,ipho,jhe2p,j108m,ielen,jpoen
c   istf,loop,iter,iel1,iel2
c   iel1,iel2,iel3,iel4 elements to follow
c   iemidx(4) mm ions to follow
c   nsteps: total number of steps in structure arrays (ie deel(mxnsteps))
c   jspot,jpre,jfin,jend
c   jden,jpres,jgeo,jeq,iso: structure settings
c     jden - density; C=isoChoric,B=isoBaric,F=function,I=input file
c     jgeo - geometry; S=spherical,P=plane,F=??
c     jtrans - transfer mode,'OUTW','LOUP','LODW','DWUP'
c     jeq  - equil.; E=equilibrium,F=finite source,P=post equil decay,
c            (C=CIE)
c     iso  - ??
c   jthm: thermal structure - should be 'S' for self consistent
c         (or T for isothermal)
c
c   jcon : flag to enable full continuum calculation
c         (used in photo4,sh2,sh3)
c   jspec: flag to enable optical flambda spectrum integration
c           (in newdif,totphot,photo4)
c   jall: flag to enable production of allion file
c   jiel: flag to enable monitoring of 4 elements
c   jiem: flag to enable monitoring of 4 ion emission
c   jlin: flag to follow 12 lines by radius
c   jsou: flag to enable production of final source file
c   jpfx: final source file prefix
c
c   jbal: flag to enable production of balance files
c   jbfx: balance file prefix
c   gdist:  grain distribution
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   general run info:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   runname: code for a given run,included in all output files
c   ionsetup: Ionisation balance setup
c   abnfile: name of abundance file or 'default'
c   deltafile: name of abundance differences file or 'default'
c   srcfile: name of photon source file or 'none'
c   depfile: depletion pattern file or 'default'
c   structfile: radius - dh file or 'none'
c   nprefix: chars in filname prefixes
c   njlines: number of emission lines monitored
c
      integer*4 iphom,ipho,jhe2p,j108m
      integer*4 ielen,jpoen,ieln,nprefix
      integer*4 istf,loop,iter,nsteps,njlines
      integer*4 emlinlistatom(mxmonlines),emlinlistion(mxmonlines)
      integer*4 iel1,iel2,iel3,iel4,iemidx(4)
      integer*4 ox1,ox2,ox3,ni2,su2,lygline
      real*8 emlinlist(mxmonlines),emlindeltas(mxmonlines)
      character jspot*4,jpre*4,jfin*4,jend*4,jsou*4,jall*4
      character jbal*4,jbfx*8,jiel*4,jiem*4,jthm*4,gdist*4
      character jden*4,jpres*4,jgeo*4,jtrans*4,jlin*4
      character jeq*4,iso*4,jcon*4,jspec*4,jpfx*16
      character runname*512,ionsetup*512,abnfile*512,srcfile*512
      character depfile*512,structfile*512,deltafile*512
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      stellar atmosphere data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   star blocks 16 is nothing to do with atoms
c
      real*8 tno4(16,11),tno5(16,11),tmet(8,11)
      real*8 tedge(16,3)
c
c
      real*8 turn,alnth,teff,cut,zgas,zstar,rstar,astar
c
c   source total photon fluxes (not just stars)
c
      real*8 qht,qhi,qhei,qheii
c
c   ** Globals for avrdata,summdata,and z* routines
c
c   declare ACCUM variables
c
      real*8 tlosac,tforbi,fheiilr
      real*8 stlwave(mxflxbins)
      real*8 stlwev(mxflxbins)
      real*8 stlhnu(mxflxbins)
      real*8 stlphnu(mxflxbins)
      real*8 stlalpha(mxflxbins)
      real*8 stlsflux(mxflxbins)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      OPTXAGNF agn Done 2012 spectrum model
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer*4 neoptx      ! actual number of energy bins
      real*8 optxparam(12)  ! optxagnf 11 paramters + 1 padding
      real*8 optxear(0:mxoptxbins)
      real*8 optxtotal(mxoptxbins),optxdisk(mxoptxbins)
      real*8 optxcoro(mxoptxbins),optxnont(mxoptxbins)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    BEGIN common block declarations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Monte-Carlo variables
c
      common /mcvar1/  convrt,mcdengy0,mcdengy,npck,npckinc
     &                ,mcnuin,mcnuend,mcnu,mcnx,mcny,mcnz
     &                ,xps,yps,zps,intpt
c
      common /mcvar2/  sigmt_arr,rec0_arr,recpdf_arr,strpdf_arr
     &                ,jnu_arr,jnuw_arr,volfilar,resflux,resflux2 
     &                ,cwx,cwy,cwz
c
      common /mcvar3/  dh_arr,de_arr,te_arr 
c
      common /mcvar4/  nrange,taskerr,taskid,numtasks
c
c
c   Below are original MAPPINGS variables
c
c
      common /environ/ theversion,dtlen,datadir,filename
c
      common /modes1/  expertmode,photonmode
     &                ,collmode,dielcmode,radpressmode
     &                ,phionmode,recommode,drrecmode
     &                ,chargemode,chargeheatmode
     &                ,freefreemode,turbheatmode
     &                ,alphacoolmode,hhecollmode,linecoolmode
     &                ,grainmode,irmode,pahmode,pahactive
     &                ,graindestructmode,grainadjust
c
      common /atdat1/ atypes,maxion,zmap,mapz
     &                ,elem,rom,elem_len,rom_len
c
      common /atdat2/ atwei,zen,zion0,zion,zsol
     &                ,deltazion,dion0,dion,invdion
     &                ,pop,dndt,popint
c
c   loss rates
c
      common /coolr1/ tloss,eloss,egain,dlos
     &                ,colos,fflos,chgain,pgain,rngain,wgain
     &                ,gheat,gcool,paheat,cosgain
     &                ,cmplos,cmpcool,cmpheat,oiii5007loss
c
      common /coolr2/ fixtemp,tem,crate,cosphi
     &                ,alphacloss,alphaclaw,alphac0
c
      common /coolr3/ coolzion,heatzion,coolz,heatz
c
c   processes
c
c     common /phovec1/ tphot,sphot,wdil,emidif,updif,dwdif
c    &                ,emilin,uplin,dwlin,soupho,dilradpho,difftot
c    &                ,dustphot,srcsum,totsum
      common /phovec1/ tphot,sphot,wdil,emidif,updif,dwdif
     &                ,emidifcont,updifcont,dwdifcont
     &                ,soupho,dilradpho,difftot
     &                ,dustphot,srcsum,totsum
c
      common /phovec2/ infph,linpos,fieldversion
     &                ,skipbin,ionstartbin
      common /phovec3/ photev,cphotev,phote,cphote
     &                ,lgphotev,lgcphotev,lgphote,lgcphote
     &                ,widbinev
c
      common /photint/ ionum,atpho,ionpho,photbinstart,shellpho
     &                ,lquantpho,hasph2,isfbpho,wrpho,wipho,ph2mode
     &                ,nelecpho,nhipx
      common /photn1/ photxsec,stwtpho,ipotpho,augpho
     &                ,ipote,ipot2e,ipotev,ipot2ev,ph2limitpho
      common /photn2/ ph1pho,ph2pho,sigpho,betpho,spho
      common /photn3/ secat,secio,secel,xsec,heaph,anr,xsect
     &                ,wnr,rasec ,rphot,auphot,secra
      common /photh1/ hipxlge,hipxlgs,hipxlgs2
c
      common /colio1/ ncol,atcol,ioncol,nsplcol
      common /colio2/ xcol,ycol,y2col,ipotcol
c
      common /recom2/ nrec,atrec,ionrec,typerrec,typedrec
      common /recom3/ arrec,adrec,bdrec
      common /recom4/ nnorad,atndrec,ionndrec,
     &                typendrec,noradid,nodesndrec
      common /recom5/ tendrec,andrec,a2ndrec
c
c Ancient rates
c
      common /charex1/ charco,charte
c
c A&R rates
c
      common /charex2/ nlegacychxr,chxrlegacyat,chxrlegacyio
     &               ,chxrlegacyx,nlegacychxi,chxilegacyat,chxilegacyio
     &               ,chxilegacyx
      common /charex3/ chxrlegacytemp,chxrlegacycos,chxrlegacy
     &                ,chxilegacytemp,chxilegacycos,chxilegacy
c
c Kingdon rates
c
      common /charex4/ nchxr,chxrat,chxrio,chxrx
     &                ,chxrnfit,chxrfitid
     &                ,nchxi,chxiat,chxiio,chxix
     &                ,chxinfit,chxifitid
c
      common /charex5/ chxrtemp,chxrcos,chxr,chxrde
     &                ,chxitemp,chxicos,chxi,chxide
c
c Old data
c
      common /arcol1/ collpot,aar,bar,car,dar,nshells,acol,tcol,col
c
      common /recom1/ arad,xrad,adi,t0,bdi,t1
     &                ,adilt,bdilt,cdilt,ddilt,fdilt,rec
c
c      common /resga1/  arg,brg,crg,drg,erg,clkup
c
c
      common /klfbgaunt/ nklgfb,nklfbe,klpxe,klgfb,klgfby2
c
c   dust grains
c
      common /grain1/ grainrad,gradedge,dustmass,siggrain
     &                ,amin,amax,dustnum
     &                ,bgrain,yinf,haccom,galpha
     &                ,segrain,spgrain,grainpot,avgpot,dustsig
     &                ,grainmaxtemp,grainscaletemp
     &                ,grainmaxq,grainscaleq,graindens
     &                ,grainghr,graindgr,grainpgr
     &                ,dustint,extinct,absorp,scatter,gcos
     &                ,dustbinmax,mindust,maxdust,numtypes
     &                ,irphot,dcrosec,dfscacros,dbscacros
     &                ,dustav,a_vend,irtemp,clinpah
c
      common /pahs/   pahip,pahz,pahion,pahlimit
     &                ,pahfrac,pahcfrac,pahnc
     &                ,pahyield,pahflux,paheng,pahint
     &      ,pahnabs,pahnsca,pahncos,pahnext
     &      ,pahiabs,pahisca,pahicos,pahiext,pahend
c
c NewIR dust parameters
c
      common /dustpar1/  engy,phot,cabs,photabs,dustsigt,collheat
      common /dustpar2/  edphot,dengy,dde,v_e,dabsphot,dele
      common /dustpar3/  bbemiss,t_grid,t_edge,tr_matrix,t_prob
      common /dustpar4/  h_grid,hmin,hmax,deltah,invh
      common /dustpar5/  grainvol,atom_no
      common /dustpar6/  flxdist,d_rad,irflux
c
c   line calculation globals
c
c      common /xline1/ xrlam,xrbri,xbr,xrloss
c     &                ,xejk,xfef,xtrans,xegj
c     &                ,xlines,minxbin,xrat,xion,xiso,xbin
c
c V3 Xray Resonance lines data
c
      common /xr3r/ xr3loss,xr3col_tc,xr3col_fr,xr3col_br
     &             ,xr3col_x,xr3col_y,xr3col_y2
     &             ,xr3ei,xr3gi,xr3invgi
c
      common /xr3rl/ xr3lines_egij
     &             ,xr3lines_lam,xr3lines_gf,xr3lines_bri,xr3lines_aji
     &             ,xr3lines_flux,xr3lines_emilin,xr3lines_frac
     &             ,xr3lines_uplin,xr3lines_dwlin
c
      common /xr3i/ nxr3ions,xr3at,xr3ion,xr3ni,nialines
     &             ,nxr3ioncol,xr3col_i,xr3col_j,xr32s12j
     &             ,xr3col_nspl,xr3col_typespl,xr3col_nl
     &             ,xr3col_jj,xr3col_kk,nxr3lines,xr3lines_map
     &             ,xr3lines_imapidx,xr3lines_imap_i,xr3lines_imap_j
     &             ,xr3lines_at,xr3lines_ion,xr3minxbin,xr3lines_bin
c
      common /xr3c/ xr3id,xr3lines_tran
c
c
c V3 Xray Resonance lines data large ions
c
      common /xrlr/ xrlloss,xrlcol_tc,xrlcol_fr,xrlcol_br
     &             ,xrlcol_x,xrlcol_y,xrlcol_y2
     &             ,xrlei,xrlgi,xrlinvgi
c
      common /xrlrl/ xrllines_egij
     &             ,xrllines_lam,xrllines_gf,xrllines_bri,xrllines_aji
     &             ,xrllines_flux,xrllines_emilin,xrllines_frac
     &             ,xrllines_uplin,xrllines_dwlin
c
      common /xrli/ nxrlions,xrlat,xrlion,xrlni
     &             ,nxrlioncol,xrlcol_i,xrlcol_j
     &             ,xrlcol_nspl,xrlcol_typespl,xrlcol_nl
     &             ,xrlcol_jj,xrlcol_kk,nxrllines,xrllines_map
     &             ,xrllines_imapidx,xrllines_imap_i,xrllines_imap_j
     &             ,xrllines_at,xrllines_ion,xrlminxbin,xrllines_bin
c
      common /xrlc/ xrlid,xrllines_tran
c
c      common /xilin1/ xilam,xibri,xiloss,brr
c     &                ,xiejk,xomeg,xitr,xiegj
c     &                ,xilin,xiuplin,xidwlin
c     &                ,xilines,xiat,xiion,xiiso,xibin
c
c     common /xhlin1/ xhelam,xhebri,xheloss
c    &                ,xhejk,xhef,xhebra
c    &                ,xhelin,xheuplin,xhedwlin,fheif
c    &                ,xhelines,xheat,xheion,xhebin
c
      common /rline1/ rlam,rbri,rloss
     &                ,rlin,ruplin,rdwlin
     &                ,e12r,omres,tkexres,fabsr
     &                ,nlines,ielr,ionr,lrbin
c
      common /cline1/ fslam,fsbri,fslos
     &                ,fslin,fsuplin,fsdwlin
     &                ,e12fs,w1fs,w2fs,a21fs,omfs,tkexfs
     &                ,mlines,ielfs,ionfs,lcbin,intam
c
      common /kappa/  kappa,kappaa,kappab
     &                ,kappaidx,usekappa,usekappainterp
c
c      common /fline1/  flam,fbri,floss,splidx5
c     &                ,wi,ei,ai,omi,tdep,tkex
c     &                ,fion,fatom,nfions,nftrans
c     &                ,flin,fuplin,fdwlin,fbin
c
      common /fmline1/  fmlam,fmeij,fmbri,fmloss
     &                ,fmuplin,fmdwlin,fmlin,fmx
     &                ,wim,eim,aim,omim,tdepm,tkexm
c
      common /fmline2/ nfmions,nfmtrans
     &                ,fmion,fmatom,fmnl,fmspindex
     &                ,fmkappaidx,fmbin
     &                ,nfmterm,nfmid,nfmtridx
     &                ,nfmlower,nfmupper
     &                ,fmomdatatype,fmombtn
c
      common /fmline3/ fmombtc,fmombsplx
     &                ,fmombsply,fmombsply2
c
      common /fmline4/ nkappaions,kapion,kapatom
     &                ,kapnl,nkaptrans,nkaptridx
     &                ,nkaplower,nkapupper
     &                ,nkappas,kupstype,kupsbtn
      common /fmline5/  kappas,kupsbtc,kupssplx
     &                ,kupssply,kupssply2
c
      common /f3line1/ f3ion,f3atom,f3type,nf3ions,nf3trans
     &                ,wi3,ei3,ai3,gam3,pow3,tkex3
     &                ,f3lam,f3bri,f3loss,f3bin
     &                ,f3lin,f3uplin,f3dwlin
c
      common /feline1/  felam,feeij,febri,feloss
     &                ,feuplin,fedwlin,felin,fex
     &                ,wife,eife,aife,omife,tdepfe,tkexfe
c
      common /feline2/ nfeions,nfetrans
     &                ,feion,featom,fenl,febin,fespindex
     &                ,nfetridx,nfelower,nfeupper,feombtn
     &                ,feomdatatype,nfeterm,nfeid
c
      common /feline3/ feombtc,feombsplx
     &                ,feombsply,feombsply2
c
      common /hline1/  balcoe,hbema,hbemb,hherat,hheen,hlam,hbri
     &                ,h2ql,hbeta,caseab,hloss,halpha,lyalpha,rateton
      common /hline2/  hydlin,hyduplin,hyddwlin
     &                ,hellin,heluplin,heldwlin
     &                ,heiioiiibf,heiioiiibfsum
      common /hline3/  collrate2p,recrate2p,ee2p,h2qflux,h2qbri
     &                ,r2s1a,r2s1b,r2s2a,r2s2b,collrate2phe,ee2phe
      common /hline4/  xhlambda,xhydroflux,xhydrobri,xhbeta
     &                ,xhydlin,xhyduplin,xhyddwlin
     &                ,ccoln,dcoln,collha,collhb,colid
      common /hlinei/  nhlines,nhseries
     &                ,nhelines,nheseries
     &                ,nxhlines,nxhseries
     &                ,hbin,hebin,xhbin
      common /hline5/  hlambda,helambda
     &                ,hylratsa,hylratsb,helratsa,helratsb
     &                ,hylratsa2,hylratsb2,helratsa2,helratsb2
     &                ,hydroflux,heliflux,hydrobri,helibri,hydrogf
c
c New H Collisions Data up to n=5 lines data
c
      common /hher/ hheloss,hhecol_tc
     &             ,hhecol_x,hhecol_y,hhecol_y2
     &             ,hheei,hhegi,hheinvgi
c
      common /hhei/ nhheions,hheat,hheion,hheni,hhen,ishheion
     &             ,nhheioncol,hhecol_i,hhecol_j
     &             ,hhecol_nspl,hhecol_typespl
c
      common /hhec/ hheid
c
      common /hrecom1/  nhrec,nherec
      common /hrecom2/  hreclgt,hreclga,hreclga2,hreclga1,hreclga12
     &                ,hereclgt,hereclga,hereclga2,hereclga1,hereclga12
      common /hrecom3/  rhhe,th42,the43,
     &                  h42a,h42b,he43a,he43b,
     &                  h42a2,h42b2,he43a2,he43b2
c
      common /helio1/  heii2ql,heien,qren,hein,heiilr,hel0
     &                ,heilam,heibri,fluxhei
     &                ,heibin,nheilines
c
      common /helio2/  nheis,nheit,nheiscr,nheitcr,nheislines,
     &nheitlines                ,heisupid,heisloid,heitupid,heitloid
     &                          ,heisbin,heitbin
c
      common /helio3/  heisej,heitej,heiseij,heiteij,heislam,heitlam
     &                ,heisreccoef,heiscra,heiscrb,heiscrc
     &                ,heitreccoef,heitcra,heitcrb,heitcrc
     &                ,heisbri,fluxheis,heitbri,fluxheit
c
      common /helioc/  heislvl,heitlvl,heislamid,heitlamid
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Heavy Element Recombination
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /reccii/ rccii_eij,rccii_lam,rccii_logte
     &                ,rccii_coeffy,rccii_coeffy2
     &                ,rccii_abri,fluxrccii_a
     &                ,rccii_bbri,fluxrccii_b
     &                ,rccii_tr,nrccii
     &                ,rccii_bin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /recnii/ rcnii_eij,rcnii_lam
     &                ,rcniin20_y ,rcniin30_y
     &                ,rcniin40_y ,rcniin50_y
     &                ,rcniin20_y2,rcniin30_y2
     &                ,rcniin40_y2,rcniin50_y2
     &                ,rcnii_abri,fluxrcnii_a
     &                ,rcnii_bbri,fluxrcnii_b
     &                ,rcnii_tr,nrcnii
     &                ,rcnii_bin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /recoi/   nrcoi_q,nrcoi_t
     &                ,rcoi_qeij,rcoi_qlam,rcoi_qcoef,rcoi_qtr
     &                ,rcoi_teij,rcoi_tlam,rcoi_tcoef,rcoi_ttr
     &                ,rcoi_qabri,fluxrcoi_qa
     &                ,rcoi_tabri,fluxrcoi_ta
     &                ,rcoi_qbbri,fluxrcoi_qb
     &                ,rcoi_tbbri,fluxrcoi_tb
     &                ,rcoi_qbin,rcoi_tbin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /recoii/ rcoii_eij,rcoii_lam
     &                ,rcoiin20_y ,rcoiin40_y
     &                ,rcoiin50_y ,rcoiin60_y
     &                ,rcoiin20_y2,rcoiin40_y2
     &                ,rcoiin50_y2,rcoiin60_y2
     &                ,rcoii_abri,fluxrcoii_a
     &                ,rcoii_bbri,fluxrcoii_b
     &                ,rcoii_cbri,fluxrcoii_c
     &                ,rcoii_tr,rcoii_mltid,nrcoii
     &                ,rcoii_bin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /recneii/ rcneii_eij,rcneii_lam
     &                ,rcneii_y, rcneii_y2
     &                ,rcneii_abri,fluxrcneii_a
     &                ,rcneii_bbri,fluxrcneii_b
     &                ,rcneii_tr,nrcneii
     &                ,rcneii_bin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   continuum calculations
c
c
      common /cndat1/ zn0,ezz,cnphot
     &               ,fn0,zt0,ffph,fbph,p2ph
     &               ,gff,gffy2,gffu,gffg2
     &               ,gffintg2,gffinty,gffinty2,gffintspl
     &               ,ngffg2,ngffu,ngffint
c
c     &               ,gfb,gfbe,gfbtz
c     &               ,ngfbe,ngfbtz
c
      common /cndat2/ fbloggamma,fblte
     &               ,fbenu,fben,fbee
     &               ,nfbenu,jfb,nfbe
      common /cndat3/ hefbloggamma,hefblte
     &               ,hefbenu,hefben,hefbee
     &               ,nhefbenu,jhefb,nhefbe
c
c     SUMMED LINES
c
      common /sline1/ fluxr,fluxi,fluxm,h2qav,hei2qa,heii2qa
     &                ,fluxh,fhbeta,pam,deav,teav
     &                ,rdisa,vunilog,fluxf3,fluxfe
c     &                ,fluxx,fluxxi,fluxfe
      common /tonlin/ roiii,rnii,weoiii,wenii,deoiii,denii,toiii,tnii
      common /denlin/ deam,dhav,roii,rsii,weoii,wesii,toii,tsii
     &                ,deoii,desii,qtosoh,qhdh,qhdn,zetae
     &                ,qhdha,zetaeav,dhavv,viofr,uhdh
c
      common /accum1/ tlosac,tforbi,fheiilr
c
c
c   PHOBAT
c
      common /phobat/ pop0,wdpl,telap,tlife,remp,rmax,distim
     &                ,volstromhb,volstromheb
     &                ,rstromhb,rstromheb
     &                ,alphahb,alphaheb
     &                ,fren,tauen,colend,tend,diend,dtau0
     &                ,tm00,aef,admach,alphaturbfixed
     &                ,qhdnin,qhdnav,qhdhin,qhdhav,uhin,uhav,unin,unav
     &                ,frlum,radout,radin
     &                ,ponk,tinner
     &                ,ponk0,ponk1,rp0,palpha
     &                ,ofac,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac
     &           ,ffac,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
c
c   shock wave blocks
c
      common /velsho/ dhpo,vshoc,vpo,vfin,hmago,fm,pres,dens
     &                ,ww,qst,qfin,depr,depo,rshock
      common /flow1/ pr0,rho0,te0,vel0,dh0,de0,bm0,pb0
     &                ,pr1,rho1,te1,vel1,dh1,de1,bm1,pb1
     &                ,gm,flowconst
c
      common /struc1/ dis0,dis1,disav,rad0,rad1

      common /shoarr/ dist,timlps
     &                ,veloc,xh,te,deel,dhy,bmg
     &                ,hoiii,hoii,hoi,hnii,hsii
c      real*8 abn(mxelem,mxnsteps)

      common /shoarr2/ popi,popf,popini,propop
     &                ,cav,diluf,xhmin,fval,texi
c
c
c   Character  and integer global variables
c
c
      common /jints1/  iphom,ipho,jhe2p,j108m,nchxold
     &                ,ielen,jpoen,nprefix,njlines
     &                ,istf,loop,iter,ieln,nsteps
     &                ,iel1,iel2,iel3,iel4,iemidx
     &                ,ox1,ox2,ox3,ni2,su2,lygline
     &                ,emlinlistatom,emlinlistion
c
      common /jreals/ emlinlist,emlindeltas
c
      common /jchars/ jspot,jpre,jfin,jend,jsou,jall,jbal,jbfx,jlin
     &                ,jden,jpres,jgeo,jeq,iso,jcon,jspec,jtrans
     &                ,jpfx,jiel,jiem,jthm ,gdist,runname,ionsetup
     &                ,abnfile,srcfile,deltafile
     &                ,depfile,structfile
c
c   star blocks
c
      common /star1/ turn,alnth,teff,cut,zgas,zstar,rstar,astar
     &                ,qht,qhi,qhei,qheii
      common /star2/ tno4,tno5,tmet,tedge
c
      common /stellib/ stlwave,stlwev,stlhnu,stlphnu,stlalpha,stlsflux
c
      common /optxr/optxparam,optxear
     &             ,optxtotal,optxdisk
     &             ,optxcoro,optxnont
      common /optxi/ neoptx
c
