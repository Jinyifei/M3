#
# Make MAPPINGS V
#
#    v5.1.13
#-------------------------------
#--- Output (executable) name --
#-------------------------------
#
OUTNAME = m3
#
#-------------------------------
#---------- Directories --------
#-------------------------------
#
EXEDIR = ../lab/
INCDIR = includes/
CODDIR = workcode/
BACKDIR = mastercode/
#
#-------------------------------
#---------- Compilers ----------
#-------------------------------
#
#  comment out and in as needed
#   for different machines
#
# GCC FORTRAN - gfortran generic
# tested with v4.8 but should work with v4.5 or newer
# on OSX, use MacPorts.org to get gcc4.8 or newer.
# In Linux, use normal repositories
#
#  GCC FORTRAN v4.8
#
#FC     = gfortran -std=legacy
#LDR    = -Ofast
#OPTS   = -c ${LDR} -I${INCDIR}
#LIB    =
#
#  GCC FORTRAN v4.9/5.x - added fpe suppression switch, otherwise as v4.8
#
#
#WARN   = -Wall
#WARN   = -Wunused-variable -Wconversion -Wmaybe-uninitialized
WARN   =
#
##FC     = gfortran -mcmodel=large -std=legacy
#LDR    =  ${WARN} -g -O1 -fbounds-check -ffpe-summary='none'
##LDR    = ${WARN} -Ofast -flto -ffpe-summary='none'
##OPTS   = -c ${LDR} -I${INCDIR}
##LIB    =
#
FC     = mpif77  -mcmodel=large -g -Wall -Wno-unused-dummy-argument
##FC     = /usr/lib64/mpich/bin/mpif77  -mcmodel=large 
LDR    = -O3
OPTS   = -c ${LDR} -I${INCDIR}
LIB    =
#
#     Generic FORTRAN
#
#FC     = f77
#LDR    = -O3
#OPTS   = -c ${LDR} -I${INCDIR}
#LIB    =
#
#     intel FORTRAN
#
#FC     = ifort
#LDR    = -O3
#OPTS   = -c ${LDR} -I${INCDIR}
#LIB    =
#
#
.SUFFIXES:
.SUFFIXES: .f .o
#
#-------------------------------
#-- include/header files -------
#-------------------------------
#
INCS =  ${INCDIR}cblocks.inc \
	${INCDIR}const.inc \
	${INCDIR}s5blocks.inc
#
#-------------------------------
#-- object files ---------------
#-------------------------------
#
OBJ = ${CODDIR}mappings.o \
	${CODDIR}allrates.o \
	${CODDIR}absdis.o \
	${CODDIR}avrdata.o \
	${CODDIR}calpdf.o \
	${CODDIR}changes.o \
	${CODDIR}charex.o \
	${CODDIR}cheat.o \
	${CODDIR}collion.o \
	${CODDIR}coloss.o \
	${CODDIR}compton.o \
	${CODDIR}coolc.o \
	${CODDIR}cool.o \
	${CODDIR}cosmic.o \
	${CODDIR}crosssections.o \
	${CODDIR}dusttemp.o \
	${CODDIR}equion.o \
	${CODDIR}evoltem.o \
	${CODDIR}fine3.o \
	${CODDIR}findtde.o \
	${CODDIR}freebound.o \
	${CODDIR}freefree.o \
	${CODDIR}functions.o \
	${CODDIR}grainpar.o \
	${CODDIR}heavyrec.o \
	${CODDIR}helioi.o \
	${CODDIR}hgrains.o \
	${CODDIR}hpahs.o \
	${CODDIR}hhecoll.o \
	${CODDIR}hydrec.o \
	${CODDIR}hydro.o \
	${CODDIR}hydro2p.o \
	${CODDIR}inter.o \
	${CODDIR}interpol.o \
	${CODDIR}intvec.o \
	${CODDIR}iobal.o \
	${CODDIR}iohyd.o \
	${CODDIR}ionemit.o \
	${CODDIR}ionab.o \
	${CODDIR}ionsec.o \
	${CODDIR}kappainit.o \
	${CODDIR}localem.o \
	${CODDIR}mapinit.o \
	${CODDIR}mdiag.o \
        ${CODDIR}mcteequi.o \
	${CODDIR}montph7.o \
        ${CODDIR}mcspec.o \
	${CODDIR}multilevel.o \
	${CODDIR}neqc.o \
	${CODDIR}netgain.o \
	${CODDIR}newdif.o \
	${CODDIR}newfil.o \
	${CODDIR}optxagnf.o \
	${CODDIR}output.o \
	${CODDIR}pheat.o \
	${CODDIR}phion.o \
	${CODDIR}phocrv.o \
	${CODDIR}photo4.o \
	${CODDIR}photo5.o \
	${CODDIR}photo6.o \
	${CODDIR}photo7.o \
	${CODDIR}photsou.o \
	${CODDIR}poputil.o \
	${CODDIR}preion.o \
	${CODDIR}qua.o \
	${CODDIR}rankhug.o \
	${CODDIR}ratec.o \
	${CODDIR}rebin.o \
	${CODDIR}recom.o \
	${CODDIR}reson.o \
	${CODDIR}sdifeq.o \
	${CODDIR}shock2.o \
	${CODDIR}shock3.o \
	${CODDIR}shock4.o \
	${CODDIR}shock5.o \
	${CODDIR}sinsla.o \
	${CODDIR}slab.o \
	${CODDIR}spectrum.o \
	${CODDIR}spotap.o \
	${CODDIR}sumdata.o \
	${CODDIR}teequi.o \
	${CODDIR}timion.o \
	${CODDIR}timtqui.o \
	${CODDIR}totphot.o \
	${CODDIR}transferline.o \
	${CODDIR}twophoton.o \
	${CODDIR}zer.o \
	${CODDIR}zetaeff.o
#
#-------------------------------
#-------  implicit rules  -----
#-------------------------------
#
.f.o:
	${FC} ${OPTS} -o $*.o $<
#
#-------------------------------
#------  explicit rules  ------
#-------------------------------
#
#
#-------------------------------
#---------  Targets  ----------
#-------------------------------
#
help:
	@echo MAPPINGS V v5.1.13 make options:
	@echo ' '
	@echo Type "'make help'    to see this menu"
	@echo Type "'make build'   to create executable from scratch and clean"
	@echo Type "'sudo make install'   to install the built executable into /usr/local"
	@echo ' '
	@echo Type "'make prepare' to copy over '*.f' and '*.inc' files"
	@echo Type "'make compile' to create executable"
	@echo Type "'make backup'  to backup '*.f' and '*.inc' files"
	@echo Type "'make clean'   to backup '*.f' and '*.inc', remove '*.o'"
	@echo Type "'make listing' to create listing of code"
#
#-----------------------------------------------------------
#
prepare:
	@echo ' Copying over .f files'
	@cp ${BACKDIR}*.f ${CODDIR}
	@echo ' Copying over include files'
	@cp ${BACKDIR}*.inc ${INCDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
clean:
	@echo ' Removing object files'
	@rm -f ${OBJ}
	@echo ' Backing up and cleaning up .f files'
	@mv $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up and cleaning up include files'
	@mv $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
backup:
	@echo ' Backing up .f files'
	@cp $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up include files'
	@cp $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
listing:
	@cat Makefile map.prefs ${INCS} $(OBJ:.o=.f) > MV.lis
#
#------------------------------------------------------------
#
compile:${EXEDIR}${OUTNAME}
#
${EXEDIR}${OUTNAME}: ${INCS} ${OBJ}
	${FC} ${LDR} -o ${EXEDIR}${OUTNAME} ${OBJ} ${LIB}

#
#------------------------------------------------------------
#
build:
	@echo ' Building ${OUTNAME}:'
	@echo ' Copying over .f files'
	@cp ${BACKDIR}*.f ${CODDIR}
	@echo ' Copying over include files'
	@cp ${BACKDIR}*.inc ${INCDIR}
	@echo ' Cleaning up src area'
	@echo ' Compiling...'
	@make compile
	@echo ' Removing object files'
	@rm -f ${OBJ}
	@echo ' Backing up and cleaning up .f files'
	@mv $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up and cleaning up include files'
	@mv $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#------------------------------------------------------------
#
install:
	[ -d /usr/local/bin ] || mkdir -p /usr/local/bin
	[ -d /usr/local/share/mappings ] || mkdir -p /usr/local/share/mappings
	cp ${EXEDIR}${OUTNAME} /usr/local/bin/
	cp ${EXEDIR}mapStd.prefs /usr/local/share/mappings/
	cp -r ${EXEDIR}data /usr/local/share/mappings/
	cp -r ${EXEDIR}abund /usr/local/share/mappings/
	cp -r ${EXEDIR}atmos /usr/local/share/mappings/
	@echo ' Installed ${OUTNAME} into /usr/local/bin and /usr/local/share/mappings'
#
installcode:
	[ -d /usr/local/bin ] || mkdir -p /usr/local/bin
	cp ${EXEDIR}${OUTNAME} /usr/local/bin/
	@echo ' Installed ${OUTNAME} into /usr/local/bin'
#
#------------------------------------------------------------
#
uninstall:
	rm -f /usr/local/bin/${OUTNAME}
	rm -rf /usr/local/share/mappings
	@echo ' Uninstalled /usr/local/bin/${OUTNAME} and /usr/local/share/mappings'

