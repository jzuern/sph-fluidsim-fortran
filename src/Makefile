# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# path to src files
SOURCEPATH = src

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = program

# "make" builds all
all: $(PROGRAMS)

#$part4
# Using Fortran MODULES: prog4.f90 USEs a Fortran module defined
# inside mod.f90, this is similar to the include case (prog2), but,
# since there is no standard naming convention for compiled module
# files in f90, the dependency is more easily built on the object
# files, because when mod.o is generated, one is sure that
# any_module_inside_mod.mod has been newly generated as well; mod.o
# must also be linked in when building the executable, so the
# dependency on mod.o is added also for prog4, as in the external
# procedure case (prog3):
program.o: module1.o gnufor2.o
program: module1.o gnufor2.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

run:
	./program

veryclean: clean
	rm -f *~ $(PROGRAMS)
