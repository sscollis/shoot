#==============================================================================
#  Makefile for shoot (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 8/3/2022
#
#==============================================================================
NAME   = shoot
#
# Turn on NR for now (needs license to use)
#
# USE_NR = 1
#
# Turn on ODEINT if desired
#
ifdef USE_ODEINT
  USE_NR = 1
  DEFINES += -DUSE_ODEINT
endif
#
OPT      = -O2
DEBUG    = -g -fbounds-check
DEFINES  += $(ADDONS)
TRAP     = -ffpe-trap=invalid
#TRAP    = -ffpe-trap=invalid,zero,overflow 
FFLAGS   = -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
           -cpp -std=legacy $(DEFINES) $(OPT) $(TRAP) $(DEBUG) -c
F90FLAGS = -fdefault-real-8 -fdefault-double-8 -cpp $(DEFINES) $(OPT) \
           $(TRAP) $(DEBUG) -c
OFLAGS   = $(OPT) $(DEBUG) 
#
# Tailor setup for OpenBLAS from HomeBrew depending on platform (Apple Intel,
# Apple Silicon, or Linux)
#
ifdef USE_LOCAL_OPENBLAS
LIB      = -L$(HOME)/local/OpenBLAS/lib -lopenblas
else ifdef USE_HOMEBREW_OPENBLAS
LIB      = -L/usr/local/opt/openblas/lib -lopenblas
else ifdef USE_APPLEBREW_OPENBLAS
LIB      = -L/opt/homebrew/opt/openblas/lib -lopenblas
else ifdef USE_LINUXBREW_OPENBLAS
LIB      = -L/home/linuxbrew/.linuxbrew/opt/openblas/lib -lopenblas
else
LIB      = -L$(HOME)/local/OpenBLAS/lib -lopenblas
endif
#LIB     = -L/usr/local/opt/openblas/lib -lopenblas
FC       = gfortran
#
# Three different ways to read and work with mean flow
#
MEAN   = mean.o
MEAN2D = mean-2d.o
MEANBL = mean-bl.o

.SUFFIXES: .f90

MODS = global.o stencils.o

OBJECTS = shoot.o matrix.o input.o getmat.o makename.o spline.o \
initial.o calch.o solve.o adjsolv.o adjini.o output.o parallel.o

OBJS1 = conte.o parder.o adjder.o rk4.o bslib1.o bslib2.o
#
# Use Numerical-Recepies only if you have a valid license
#
ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr
else
  $(warning SHOOT currently requires that you build with USE_NR defined)
  $(error See README.md for details))
  #$(info build will fail at link stage.)
endif

OBJS2 = grad.o grad2.o g1.o

OBJS3 = nonpar.o

all: $(NAME) $(NAME)-2d $(NAME)-bl 

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(MEAN) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEAN) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)
	\cp $(NAME) $(NAME).exe

$(NAME)-2d: $(MODS) $(OBJECTS) $(OBJS1) $(MEAN2D) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEAN2D) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)-2d
	\cp $(NAME)-2d $(NAME)-2d.exe

$(NAME)-bl: $(MODS) $(OBJECTS) $(OBJS1) $(MEANBL) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEANBL) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)-bl 
	\cp $(NAME)-bl $(NAME)-bl.exe

$(OBJECTS): global.o
	$(FC) $(F90FLAGS) $*.f90

$(MODS):
	$(FC) $(F90FLAGS) $*.f90

$(OBJS2): stencils.o
	$(FC) $(F90FLAGS) $*.f90

$(OBJS3):
	$(FC) -DIMSL $(F90FLAGS) $*.f90

.f90.o:
	$(FC) $(F90FLAGS) $*.f90

.f.o:
	$(F77) $(FFLAGS) $*.f

clean:
	$(RM) *.o *.mod $(NAME) $(NAME)-2d $(NAME)-bl *.exe
