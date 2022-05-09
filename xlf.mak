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
FFLAGS   = -qrealsize=8 -qsuppress=cmpmsg -qsuffix=f=f:cpp=f -qfixed=120 \
	   -c $(DEBUG)
F90FLAGS = -qrealsize=8 -qsuppress=cmpmsg -qsuffix=f=f90:cpp=f90 -WF,-DXLF\
	   -c $(DEBUG)
OFLAGS   = 
LIB      = -L/Users/sscoll/dist/atlas/lib -llapack -latlas -lg2c \
	   -Wl,-framework,accelerate
F90COMP  = xlf90
FCOMP    = xlf

.SUFFIXES: .f90 

MODS = global.o stencils.o

OBJECTS = shoot.o matrix.o input.o getmat.o makename.o spline.o \
initial.o calch.o solve.o adjsolv.o adjini.o output.o parallel.o

OBJS1 = odeint.o conte.o parder.o adjder.o mean2d.o rtsafe.o

OBJS2 = grad.o grad2.o g1.o

OBJS3 = nonpar.o

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3)
	$(F90COMP) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3) \
	$(LIB) -o $(NAME)

$(OBJECTS): global.o

$(MODS):

$(OBJS2): stencils.o

$(OBJS3):

.f90.o:
	$(F90COMP) $(F90FLAGS) $*.f90 

.f.o:
	$(FCOMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod
