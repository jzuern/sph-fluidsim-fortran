# The compiler
FC = f95

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check

# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS +=

# path to src files
SOURCEPATH = src


# libraries needed for linking: openMP, fortranposix, gnuplotfortran
LDFLAGS =-fopenmp -lfortranposix -lgnuplotfortran

# List of executables to be built within the package
PROGRAMS = program

# "make" builds all
all: $(PROGRAMS)
program.o:  gnufor2.o \
						util.o \
						linkedlists.o \
						sphfunctions.o  \
						integrate.o

program:    gnufor2.o \
						util.o \
					  linkedlists.o \
						sphfunctions.o \
						integrate.o


# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	# $(FC) $(FCFLAGS) -o $@ $^
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: $(SOURCEPATH)/%.f90
	# $(FC) $(FCFLAGS) -c $<
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
run:
	./program

veryclean: clean
	rm -f *~ $(PROGRAMS)
