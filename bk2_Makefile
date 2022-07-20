# --- COMPILER ----------------------------------------
#CC = mpicc
#CC = /usr/lib64/mpi/gcc/openmpi2/bin/mpicc
CC = /usr/lib/hpc/gnu7/mpi/openmpi/3.1.4/bin/mpicc

# --- CFLAGS -----------------------------------------
CFLAGS_gnu = -std=gnu99 -Wall -pedantic -O3 -ffast-math -fopenmp -lblas -llapack 
#CFLAGS_intel = -std=gnu99 -Wall -pedantic -O3  -xHOST -qopenmp 
#CFLAGS = $(CFLAGS_intel)
CFLAGS = $(CFLAGS_gnu)

# --- DO NOT CHANGE -----------------------------------
CPP = cpp
MAKEDEP = $(CPP) -MM
SRCDIR = src
BUILDDIR = build
BINDIR=bin
LIBDIR=lib
INCDIR=include
TSTDIR=tests
DOCDIR=doc
GSRCDIR = $(BUILDDIR)/gsrc
SRC = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.c,$(wildcard $(SRCDIR)/*.c)))
TSTS = $(patsubst %.c,%,$(wildcard $(TSTDIR)/*.c))
LIB =  $(LIBDIR)/libDDalphaAMG.a $(LIBDIR)/libDDalphaAMG_devel.a $(INCDIR)/DDalphaAMG.h
SRCGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.c))
GSRCFLT = $(patsubst %_generic.c,$(GSRCDIR)/%_float.c,$(SRCGEN))
GSRCDBL = $(patsubst %_generic.c,$(GSRCDIR)/%_double.c,$(SRCGEN))
GSRC = $(patsubst %,$(GSRCDIR)/%,$(SRC)) $(GSRCFLT) $(GSRCDBL)
HEA = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.h,$(wildcard $(SRCDIR)/*.h)))
HEAGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.h))
GHEAFLT = $(patsubst %_generic.h,$(GSRCDIR)/%_float.h,$(HEAGEN))
GHEADBL = $(patsubst %_generic.h,$(GSRCDIR)/%_double.h,$(HEAGEN))
GHEA = $(patsubst %,$(GSRCDIR)/%,$(HEA)) $(GHEAFLT) $(GHEADBL)
OBJ = $(patsubst $(GSRCDIR)/%.c,$(BUILDDIR)/%.o,$(GSRC))
OBJDB = $(patsubst %.o,%_devel.o,$(OBJ))
DEP = $(patsubst %.c,%.dep,$(GSRC))

# --- FLAGS FOR HDF5 ---------------------------------
# H5FLAGS=-DHAVE_HDF5 /usr/include
# H5LIB=-lhdf5 -lz

# --- FLAGS FOR LIME ---------------------------------
LIMEDIR = /home/ramirez/installs/qio/bin
LIMEFLAGS = -DHAVE_LIME -I$(LIMEDIR)/include
LIMELIB = $(LIMEDIR)/lib64/liblime.a

# Available flags:
# -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING
# -DSINGLE_ALLREDUCE_ARNOLDI
# -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS -DDEBUG

OPT_VERSION_FLAGS = $(CFLAGS) $(LIMEFLAGS) $(H5FLAGS) -DPARAMOUTPUT -DTRACK_RES -DOPENMP -DPROFILING
#OPT_VERSION_FLAGS += -DGCRODR
#OPT_VERSION_FLAGS += -DSINGLE_ALLREDUCE_ARNOLDI -DPIPELINED_ARNOLDI
#OPT_VERSION_FLAGS += -DBLOCK_JACOBI -DLOC_POLYPREC
#OPT_VERSION_FLAGS += -DPOLYPREC
#OPT_VERSION_FLAGS += $(LIBSMUMPS)

DEVEL_VERSION_FLAGS = $(CFLAGS) $(LIMEFLAGS) -DDEBUG -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS -DOPENMP
#DEVEL_VERSION_FLAGS += -DGCRODR
#DEVEL_VERSION_FLAGS += -DSINGLE_ALLREDUCE_ARNOLDI -DPIPELINED_ARNOLDI
#DEVEL_VERSION_FLAGS += -DBLOCK_JACOBI -DLOC_POLYPREC
#DEVEL_VERSION_FLAGS += -DPOLYPREC
#DEVEL_VERSION_FLAGS += $(LIBSMUMPS)

#---------------------------------------------------

# clean integration of MUMPS within DDalphaAMG

MUMPSDIR = /usr/local/sw/mumps-5.4.1/
MUMPS_LIBS = $(MUMPSDIR)lib/

LMETISDIR = /usr/local/sw/metis-5.1.0/install/lib/
LMETIS=-L$(LMETISDIR) -lmetis
LPMETISDIR = /usr/local/sw/parmetis-4.0.3/install/lib/
LPMETIS=-L$(LPMETISDIR) -lparmetis

LSCOTCHDIR = /usr/local/sw/scotch-6.1.1/lib/
LSCOTCH = -L$(LSCOTCHDIR) -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr

LPORD=-L$(MUMPS_LIBS) -lpord
LIBMUMPS_COMMON = -L$(MUMPS_LIBS)/ -lmumps_common

LORDERINGS = $(LPMETIS) $(LMETIS) $(LSCOTCH) $(LPORD) -L/usr/lib/hpc/gnu7/mpi/openmpi/3.1.4/lib64/ -lmpi -lmpi_mpifh -lmpi_usempif08 -lmpi_usempi_ignore_tkr

# Enable the following three to activate MUMPS within DDalphaAMG
LIBSMUMPS = -L$(MUMPS_LIBS) -lcmumps -ldmumps -lmumps_common -lpord -lsmumps -lzmumps $(LIBMUMPS_COMMON) $(LORDERINGS) -lpthread -lz
OPT_VERSION_FLAGS += -DMUMPS_ADDS
DEVEL_VERSION_FLAGS += -DMUMPS_ADDS

#---------------------------------------------------

# some old (unused) linkages

#LIBBLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#SCALAP=-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lblacs
#LIBPAR = $(SCALAP) $(LIBBLAS)
##Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
#CDEFS   = -DAdd_

##Begin Optimized options
#OPTF    = -O -nofor_main -qopenmp -Dintel_ -DALLOW_NON_INIT
#OPTL    = -O -nofor_main -qopenmp
#OPTC    = -O -qopenmp
##End Optimized options

#---------------------------------------------------

# lapack and scalapack linkages

LAPACK_DIR = dependencies/lapack-3.9.0
LAPACKE_DIR = $(LAPACK_DIR)/LAPACKE
LAPACKE_INCLUDE = $(LAPACKE_DIR)/include
BLASLIB      = $(LAPACK_DIR)/librefblas.a
LAPACKLIB    = $(LAPACK_DIR)/liblapack.a
LAPACKELIB   = $(LAPACK_DIR)/liblapacke.a
LAPACK_LIBRARIES = $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB)

#SPBLAS_DIR = dependencies/spblas
#SPBLASLIB = $(SPBLAS_DIR)/libsparseblas.a
#SPBLAS_LIBRARIES = $(SPBLASLIB)
SPBLAS_LIBRARIES = 

SCALAPACK_DIR = /usr/lib/hpc/gnu7/openmpi3/scalapack/2.0.2/
SCALAPACK_INCLUDE = -I$(SCALAPACK_DIR)/include/
SCALAPACK_LIBRARIES = -L$(SCALAPACK_DIR)/lib64/ -lscalapack -lblacs

#---------------------------------------------------

all: execs library exec-tests
execs: $(BINDIR)/DDalphaAMG #$(BINDIR)/DDalphaAMG_devel
library: $(LIB)
exec-tests: $(TSTS)
documentation: $(DOCDIR)/user_doc.pdf
install: copy

.PHONY: all wilson library
.SUFFIXES:
.SECONDARY:

$(BINDIR)/DDalphaAMG : $(OBJ) 
	mpif90 $(OPT_VERSION_FLAGS) -o $@ $(OBJ) $(H5LIB) $(LIMELIB) $(SCALAPACK_LIBRARIES) $(LAPACK_LIBRARIES) $(SPBLAS_LIBRARIES) -lm -lgfortran $(LIBSMUMPS)
#		$(MUMPS_LIBRARIES)

DDalphaAMG : $(BINDIR)/DDalphaAMG
	ln -sf $(BINDIR)/$@ $@

#DDalphaAMG_devel: $(BINDIR)/DDalphaAMG_devel
#	ln -sf $(BINDIR)/$@ $@

#$(BINDIR)/DDalphaAMG_devel : $(OBJDB)
#	$(CC) -g $(DEVEL_VERSION_FLAGS) -o $@ $(OBJDB) $(H5LIB) $(LIMELIB) $(SCALAPACK_LIBRARIES) $(LAPACK_LIBRARIES)  -lm -lgfortran
#		$(MUMPS_LIBRARIES)

$(LIBDIR)/libDDalphaAMG.a: $(OBJ)
	ar rc $@ $(OBJ)
	ar d $@ main.o
	ranlib $@

$(LIBDIR)/libDDalphaAMG_devel.a: $(OBJDB)
	ar rc $@ $(OBJDB)
	ar d $@ main.o
	ranlib $@

$(TSTDIR)/%: $(LIB) $(TSTDIR)/%.c
	$(CC) -DAdd_ $(CFLAGS) -o $@ $@.c -I$(INCDIR) -I$(LAPACKE_INCLUDE) -L$(LIBDIR) -lDDalphaAMG $(LIMELIB) $(SCALAPACK_LIBRARIES) $(LAPACK_LIBRARIES) -lm -lgfortran $(LIBSMUMPS)

$(DOCDIR)/user_doc.pdf: $(DOCDIR)/user_doc.tex $(DOCDIR)/user_doc.bib
	( cd $(DOCDIR); pdflatex user_doc; bibtex user_doc; pdflatex user_doc; pdflatex user_doc; )

$(INCDIR)/%: $(SRCDIR)/%
	cp $(SRCDIR)/`basename $@` $@

$(BUILDDIR)/%.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) -DAdd_ $(OPT_VERSION_FLAGS) -I$(LAPACKE_INCLUDE) -c $< -o $@ $(LIBSMUMPS)

$(BUILDDIR)/%_devel.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) -g -DAdd_ $(DEVEL_VERSION_FLAGS) -I$(LAPACKE_INCLUDE) -c $< -o $@ $(LIBSMUMPS)

$(GSRCDIR)/%.h: $(SRCDIR)/%.h $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

$(GSRCDIR)/%.c: $(SRCDIR)/%.c $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

%.dep: %.c $(GHEA)
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' > $@
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1_devel.o $@ : ,g' >> $@

copy: $(BINDIR) $(LIBDIR) $(INCDIR) 
	cp -r $(BINDIR)/ $(PREFIX)
	cp -r $(LIBDIR)/ $(PREFIX)
	cp -r $(INCDIR)/ $(PREFIX)

clean:
	rm -f $(BUILDDIR)/*.o
	rm -f $(GSRCDIR)/*

cleanall: clean
	rm -f $(BINDIR)/*
	rm -f $(LIBDIR)/*
	rm -f $(INCDIR)/*
	rm -f $(SRCDIR)/*~
	rm -f $(TSTDIR)/*~
	rm -f $(TSTS)
	rm -f *~

-include $(DEP)
