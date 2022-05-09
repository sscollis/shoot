#==============================================================================
#  Makefile for shoot (Cray)
#
#  Author:  Scott Collis
#
#  Revised: 9-26-96 
#
#  For perftrace, compile f77 with -F flag, 
#  and f90 with -ef flag, and link with -lperf
#==============================================================================
NPROC  = 4
NAME   = shoot 
DEBUG  = 
FFLAGS = -c $(DEBUG)
F90FLAGS = -eZ -DCRAY -DIMSL -c $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -L/usr/local/lib -limsl
COMP   = f90

.SUFFIXES: .f90 

MODS = global.o stencils.o

OBJECTS = shoot.o matrix.o input.o getmat.o makename.o spline.o initial.o \
calch.o solve.o adjsolv.o adjini.o output.o parallel.o

OBJS1 = odeint.o conte.o parder.o adjder.o mean2d.o rtsafe.o

OBJS2 = grad.o grad2.o g1.o
  
OBJS3 = nonpar.o
  
$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) $(OBJS3) $(LIB)

$(OBJECTS): global.o
	$(COMP) $(F90FLAGS) -p global.o $*.f90

($MODS):
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS2): stencils.o
	$(COMP) $(F90FLAGS) -p stencils.o $*.f90

$(OBJS3):
	$(COMP) $(F90FLAGS) $*.f90

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	cf77 $(FFLAGS) $*.f

clean:
	/bin/rm *.o
