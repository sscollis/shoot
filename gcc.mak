#==============================================================================
#  Makefile for shoot (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 5/9/2022
#
#==============================================================================
NAME     = shoot
DEBUG    =
FFLAGS   = -O2 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
           -cpp -c -std=legacy $(DEBUG)
F90FLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -cpp -c $(DEBUG)
OFLAGS   = -O2 $(DEBUG) -o $(NAME)
#LIB     = -lcomplib.sgimath /usr/people/collis/lib/bslib/bslib.a
LIB      = -L/usr/local/opt/openblas/lib -lopenblas
FC       = gfortran

#
# Currently the type of meanflow is set at build time
#
MEAN = mean.o
ifdef MEAN_2D
  MEAN = mean2d.o
else
  ifdef MEAN_BL
    MEAN = mean_bl.o
  endif
endif

.SUFFIXES: .f90

MODS = global.o stencils.o

OBJECTS = shoot.o matrix.o input.o getmat.o makename.o spline.o \
initial.o calch.o solve.o adjsolv.o adjini.o output.o parallel.o

OBJS1 = conte.o parder.o adjder.o $(MEAN) rk4.o bslib1.o bslib2.o
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

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3)
	$(FC) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB)

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
