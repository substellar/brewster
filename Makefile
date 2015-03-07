#$preamble
# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check 
#FCFLAGS = -frecord-marker=4 -ffixed-line-length-132 -fdefault-double-8 -fdefault-real-8 
FCFLAGS = -O2
#FCFLAGS = -fno-leading-underscore
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS = -L/usr/local/lib/

# List of executables to be built within the package
PROGRAMS = main

# "make" builds all
all: $(PROGRAMS)

#  build the main to an executable program
#shiner: main.o

# some dependencies
common_arrays_mod.o: sizes_mod.o define_types_mod.o
define_types_mod.o: sizes_mod.o 
atmos_ops_mod.o: sizes_mod.o phys_const_mod.o
gas_mixing_mod.o sizes_mod.o phys_const_mod.o define_types_mod.o
main.o: sizes_mod.o  define_types_mod.o common_arrays_mod.o phys_const_mod.o atmos_ops_mod.o gas_mixing_mod.o

main: sizes_mod.o common_arrays_mod.o define_types_mod.o phys_const_mod.o atmos_ops_mod.o gas_mixing_mod.o main.o





# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

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
.PHONY: clean 

clean:
	rm -f *.o *.mod *.MOD
