#==============================================================================
#  Makefile for shoot (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 1-29-97
#
#  Currently needs IMSL library
#
#==============================================================================
NAME     = shoot
DEBUG    =
FFLAGS   = -O2 -r8 -cpp -c $(DEBUG)
F90FLAGS = -O2 -r8 -cpp -c $(DEBUG)
OFLAGS   = -O2 -r8 $(DEBUG) -o $(NAME)
LIB      = -lcomplib.sgimath /usr/people/collis/lib/bslib/bslib.a
COMP     = f90

.SUFFIXES: .f90 

MODS = global.o stencils.o

OBJECTS = shoot.o matrix.o input.o getmat.o makename.o spline.o \
initial.o calch.o solve.o adjsolv.o adjini.o output.o parallel.o

OBJS1 = odeint.o conte.o parder.o adjder.o mean2d.o rtsafe.o

OBJS2 = grad.o grad2.o g1.o

OBJS3 = nonpar.o

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3)
	$(COMP) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB)

$(OBJECTS): global.o
	$(COMP) $(F90FLAGS) $*.f90

$(MODS):
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS2): stencils.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS3):
	$(COMP) -DIMSL $(F90FLAGS) $*.f90

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	f77 $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod
