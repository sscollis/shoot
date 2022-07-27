#==============================================================================
#  Makefile for shoot (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 5/9/2022
#
#==============================================================================
NAME   = shoot
#
# Turn on NR for now (needs license to use)
#
USE_NR = 1
#
# Turn on ODEINT if desired
#
ifdef USE_ODEINT
  USE_NR = 1
  DEFINES += -DUSE_ODEINT
endif
#
DEBUG    = -g -fbounds-check
#TRAP    = -ffpe-trap=invalid,zero,overflow 
FFLAGS   = -O2 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
           -cpp -std=legacy $(DEFINES) $(TRAP) $(DEBUG) -c
F90FLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -cpp $(DEFINES) \
           $(TRAP) $(DEBUG) -c
OFLAGS   = -O2 $(DEBUG) 
LIB      = -L/usr/local/opt/openblas/lib -lopenblas
FC       = gfortran
#
# Three different ways to read and work with mean flow
#
MEAN = mean.o
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
  $(info build will fail at link stage.)
endif

OBJS2 = grad.o grad2.o g1.o

OBJS3 = nonpar.o

all: $(NAME) $(NAME)-2d $(NAME)-bl 

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(MEAN) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEAN) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)

$(NAME)-2d: $(MODS) $(OBJECTS) $(OBJS1) $(MEAN2D) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEAN2D) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)-2d

$(NAME)-bl: $(MODS) $(OBJECTS) $(OBJS1) $(MEANBL) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(MEANBL) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB) -o $(NAME)-bl 

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
	$(RM) *.o *.mod $(NAME)
