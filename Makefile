#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.c source:tools:main:test
vpath %.o build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
#CC       = gcc
CC       = icc
#CC       = pgcc

# your Fortran compiler:
#FC       = gfortran
FC       = ifort #-nofor_main

# your tool for creating static libraries:
AR        = ar rv

# your tool for creating shared libraries:
LD        = ld -shared

# Your python interpreter.
# In order to use Python 3, you can manually
# substitute python3 to python in the line below, or you can simply
# add a compilation option on the terminal command line:
# "PYTHON=python3 make all" (THanks to Marius Millea for pyhton3
# compatibility)
PYTHON ?= python

# your optimization flag
#OPTFLAG = -O4 -ffast-math #-march=native
#OPTFLAG = -Ofast -ffast-math #-march=native
#OPTFLAG = -fast
OPTFLAG = -O3 -march=native

# your openmp flag (comment for compiling without openmp)
OMPFLAG   = -qopenmp
#OMPFLAG   = -mp -mp=nonuma -mp=allcores -g
#OMPFLAG   = -openmp

# all other compilation flags
CCFLAG = -fPIC -g
FCFLAG = -fPIC -g
LDFLAG = -fPIC -g

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. hyrec or ../hyrec)
HYREC = hyrec

########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'# -DTHERMO_DBUG

# where to find include files *.h
INCLUDES = -I../include

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL = bispev.o fpbisp.o fpbspl.o

# any dependencies
array.o: bispev.o
bispev.o: fpbisp.o
fpbisp.o: fpbspl.o

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../$(HYREC)/include 
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o energy_injection.o hyrec_21cm.o
endif

vpath %.f dierckx

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.f .base
	cd $(WRKDIR);$(FC) $(OPTFLAG) $(OMPFLAG) $(FCFLAG) $(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o common_pbh.o

SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

NONLINEAR = nonlinear.o

LENSING = lensing.o

OUTPUT = output.o

CLASS = class.o

TEST_LOOPS = test_loops.o

TEST_LOOPS_OMP = test_loops_omp.o

TEST_DEGENERACY = test_degeneracy.o

TEST_TRANSFER = test_transfer.o

TEST_NONLINEAR = test_nonlinear.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_SIGMA = test_sigma.o

TEST_HYPERSPHERICAL = test_hyperspherical.o

TEST_BICUBIC_SPLINE = test_bicubic_spline.o

TEST_STEPHANE = test_stephane.o

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(SOURCE) $(OUTPUT))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_DEGENERACY) $(TEST_LOOPS) $(TEST_TRANSFER) $(TEST_NONLINEAR) $(TEST_PERTURBATIONS) $(TEST_THERMODYNAMICS))))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(CLASS))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE)
H_ALL = $(addprefix include/, common.h common_pbh.h svnversion.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
PRE_ALL = cl_ref.pre clt_permille.pre
INI_ALL = explanatory.ini lcdm.ini
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py

all: class libclass.so classy

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT))

libclass.so: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT)
	$(CC) -shared -o $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT))

class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) -lm

test_sigma: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_SIGMA)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_sigma $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_loops_omp: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS_OMP)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_stephane: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_STEPHANE)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_degeneracy: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_DEGENERACY)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_nonlinear: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_NONLINEAR)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) -lm

test_bicubic_spline: $(TOOLS) $(TEST_BICUBIC_SPLINE) $(EXTERNAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_bicubic_spline $(addprefix build/,$(notdir $^)) -lm -lgfortran


tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

classy: libclass.a python/classy.pyx python/cclassy.pxd
ifdef OMPFLAG
	cp python/setup.py python/autosetup.py
else
	grep -v "lgomp" python/setup.py > python/autosetup.py
endif
	cd python; export CC=$(CC); $(PYTHON) autosetup.py install || $(PYTHON) autosetup.py install --user
	rm python/autosetup.py

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.*
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build
