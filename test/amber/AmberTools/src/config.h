#  Amber configuration file.
#  Created at Wed May 17 12:19:33 PDT 2017 via ./configure -mpi --with-python /home/robin/miniconda/bin/python --python-install local gnu

###############################################################################

# (1)  Location of the installation

BASEDIR=$(HOME)/amber
AMBER_PREFIX=$(HOME)/amber
BINDIR=$(HOME)/amber/bin
LIBDIR=$(HOME)/amber/lib
INCDIR=$(HOME)/amber/include
DATDIR=$(HOME)/amber/dat
LOGDIR=$(HOME)/amber/logs
AMBER_SOURCE=$(HOME)/amber/

###############################################################################


#  (2) If you want NAB to search additional libraries by default, add them
#      to the FLIBS variable here.  (External libraries can also be linked into
#      NAB programs simply by including them on the command line; libraries
#      included in FLIBS are always searched.)

FLIBS= -lpbsa -lfftw3 -larpack -llapack -lblas  -lnetcdf  -lgfortran -w 
FLIBS_PTRAJ= -larpack -llapack -lblas   -lgfortran -w
FLIBSF= -larpack -llapack -lblas  
FLIBS_FFTW3=  -lfftw3
###############################################################################

#  (3)  Modify any of the following if you need to change, e.g. to use gcc
#        rather than cc, etc.

SHELL=/bin/sh
INSTALLTYPE=parallel
BUILDAMBER=amber

#  Set the C compiler, etc.

#  The configure script should be fine, but if you need to hand-edit,
#  here is some info:

#   Example:  CC-->gcc; LEX-->flex; YACC-->yacc (built in byacc)
#     Note: If your lexer is "really" flex, you need to set
#     LEX=flex below.  For example, on some distributions,
#     /usr/bin/lex is really just a pointer to /usr/bin/flex,
#     so LEX=flex is necessary.  In general, gcc seems to need flex.

#   The compiler flags CFLAGS and CXXFLAGS should always be used.
#   By contrast, *OPTFLAGS and *NOOPTFLAGS will only be used with
#   certain files, and usually at compile-time but not link-time.
#   Where *OPTFLAGS and *NOOPTFLAGS are requested (in Makefiles,
#   makedepend and depend), they should come before CFLAGS or
#   CXXFLAGS; this allows the user to override *OPTFLAGS and
#   *NOOPTFLAGS using the BUILDFLAGS variable.

#   AMBERBUILDFLAGS provides a hook into all stages of the build process.
#   It can be used to build debug versions, invoke special features, etc.
#   Example:  make AMBERBUILDFLAGS='-O0 -g' sander
#
CC=gcc
CFLAGS=-fPIC  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DHASGZ -DHASBZ2 -D__PLUMED_HAS_DLOPEN $(CUSTOMBUILDFLAGS)  $(AMBERBUILDFLAGS)
CNOOPTFLAGS=
COPTFLAGS=-O3 
AMBERCFLAGS= $(AMBERBUILDFLAGS)
WARNFLAGS=-Wall -Wno-unused-function

CXX=g++
CPLUSPLUS=g++
CXXFLAGS=-fPIC  $(CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
CXXNOOPTFLAGS=
CXXOPTFLAGS=-fPIC -O3 
AMBERCXXFLAGS= $(AMBERBUILDFLAGS)

NABFLAGS= $(AMBERBUILDFLAGS)
PBSAFLAG=$(AMBERBUILDFLAGS)

FP_FLAGS=

LDFLAGS= $(CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
AMBERLDFLAGS=$(AMBERBUILDFLAGS)

LEX=   flex
YACC=  yacc
AR=    ar rv
M4=    m4
RANLIB=ranlib

#  Set the C-preprocessor.  Code for a small preprocessor is in
#    ucpp-1.3;  it gets installed as $(BINDIR)/ucpp;

CPP=ucpp -l

#  These variables control whether we will use compiled versions of BLAS
#  and LAPACK (which are generally slower), or whether those libraries are
#  already available (presumably in an optimized form).

LAPACK=install
BLAS=install
F2C=skip

#  These variables determine whether builtin versions of certain components
#  can be used, or whether we need to compile our own versions.

UCPP=install
C9XCOMPLEX=skip

#  For Windows/cygwin, set SFX to ".exe"; for Unix/Linux leave it empty:
#  Set OBJSFX to ".obj" instead of ".o" on Windows:

SFX=
OSFX=.o
MV=mv
RM=rm
CP=cp
WINE=

#  Information about Fortran compilation:

FC=gfortran
FFLAGS= -fPIC  $(LOCALFLAGS) $(CUSTOMBUILDFLAGS) -I$(INCDIR) $(NETCDFINC)  $(AMBERBUILDFLAGS)
FNOOPTFLAGS= -O0
FOPTFLAGS= -O3
AMBERFFLAGS=$(AMBERBUILDFLAGS)
FREEFORMAT_FLAG= -ffree-form
LM=-lm
FPP=cpp -traditional -P
FPPFLAGS= -DBINTRAJ -DEMIL (CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
AMBERFPPFLAGS=$(AMBERBUILDFLAGS)
FCREAL8=-fdefault-real-8
NOFORTRANMAIN=-lgfortran -w
FWARNFLAGS=-Wall -Wno-unused-function

XHOME= /usr
XLIBS= 
MAKE_XLEAP=install_xleap

NETCDF=$(INCDIR)/netcdf.mod
NETCDFLIB=$(HOME)/miniconda/envs/test-environment/lib/libnetcdf.a -L$(HOME)/miniconda/envs/test-environment/lib -lhdf5 -lhdf5_hl -lmfhdf
NETCDFINC=-I$(HOME)/miniconda/envs/test-environment/include
PNETCDFLIB=
PNETCDFINC=
PNETCDFDEF=
FFTWLIB=
SANDERAPI_LIB=
SANDERAPI_DEF=
SANDERAPI_DEP=
BUILD_SANDERAPI=skip_sanderapi

EMIL=EMIL
EMILLIB=$(LIBDIR)/libemil.a -lstdc++ 

ZLIB=-lz
BZLIB=-lbz2

HASFC=yes
MTKPP=
XBLAS=
FFTW3=$(LIBDIR)/libfftw3.a
MDGX=parallel

COMPILER=gnu
MKL=
MKL_PROCESSOR=
READLINE=readline/libreadline.a

#CUDA Specific build flags
NVCC=
PMEMD_CU_INCLUDES=
PMEMD_CU_LIBS=
PMEMD_CU_DEFINES=
PBSA_CU_LIBS=

#PMEMD Specific build flags
PMEMD_F90=mpif90 -DMPI   -DBINTRAJ -DEMIL -DPUBFFT
PMEMD_FOPTFLAGS=-O3 $(AMBERBUILDFLAGS)
PMEMD_CC=mpicc
PMEMD_COPTFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DMPI  $(AMBERBUILDFLAGS)
PMEMD_FLIBSF=  $(LIBDIR)/libemil.a -lstdc++ 
PMEMD_LD=mpif90 $(AMBERBUILDFLAGS)
LDOUT= -o 
PMEMD_GNU_BUG303=-fno-tree-vectorize

#for NAB:
MPI=mpi

#1D-RISM
RISM=no

#3D-RISM NAB
RISMSFF=
SFF_RISM_INTERFACE=
TESTRISMSFF=

#3D-RISM SANDER
RISMSANDER=
SANDER_RISM_INTERFACE=
FLIBS_RISMSANDER=

#for EMIL:
EMIL_MPIFLAGS=-DUSE_MPI

#PUPIL
PUPILLIBS=-lrt -lm -lc -L${PUPIL_PATH}/lib -lPUPIL -lPUPILBlind

#Python interpreter we are using and install options
PYTHON=/home/robin/miniconda/bin/python
PYTHON_INSTALL=--prefix=$(AMBER_PREFIX)
SKIP_PYTHON=no

PYSANDER=skip
PYTRAJ=pytraj
MAKE_SAXS=install

#For LIO QM GPU Library
LIOLIBS=

# OS-specific rules for making shared objects
SHARED_SUFFIX=.so
MAKE_SHARED=-shared

# PLUMED related variables:
PLUMED_INCLUDE_FILE=
PLUMED_LOAD=
PLUMED_DEPENDENCIES=
