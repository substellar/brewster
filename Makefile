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
# FC = ifort
# flags for debugging or for maximum performance, comment as necessary

# Debug version
#FCFLAGS =  -fPIC -fbounds-check -fbacktrace -Og -frecord-marker=4
#F77FLAGS = -frecord-marker=4 -fPIC -Og -fbounds-check -fbacktrace -fdefault-double-8 -fdefault-real-8 -std=legacy


# Run version
FCFLAGS = -O2 -fPIC -ffree-line-length-none -frecord-marker=4 -fbounds-check 

F77FLAGS = -O2 -fPIC -fbounds-check -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -std=legacy

# Intel compiler version
# FCFLAGS = -O2 -fPIC
# F77FLAGS = -O2 -fPIC -autodouble


# F77FLAGS =  -ffixed-line-length-132 -fdefault-double-8 -fdefault-real-8 -g -Og -fbounds-check -fbacktrace

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include
F77FLAGS += -I/usr/include
#F7FLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS = -L/usr/local/lib/

# List of executables to be built within the package
#PROGRAMS = forward

# "make" builds all
all: $(PROGRAMS)

#  build the main to an executable program
#shiner: main.o

# some dependencies
common_arrays_mod.o: sizes_mod.o define_types_mod.o
define_types_mod.o: sizes_mod.o 
atmos_ops_mod.o: sizes_mod.o common_arrays_mod.o phys_const_mod.o define_types_mod.o
gas_opacity_mod.o: sizes_mod.o  common_arrays_mod.o phys_const_mod.o define_types_mod.o atmos_ops_mod.o tristan.o
clouds_mod.o: sizes_mod.o common_arrays_mod.o define_types_mod.o phys_const_mod.o  
setup_RT_mod.o:sizes_mod.o common_arrays_mod.o define_types_mod.o phys_const_mod.o atmos_ops_mod.o bits_for_disort_f77.o DISORT.o dsolver.o gfluxi.o 

main_mod.o: sizes_mod.o common_arrays_mod.o define_types_mod.o phys_const_mod.o atmos_ops_mod.o tristan.o gas_opacity_mod.o clouds_mod.o RDI1MACH.o LINPAK.o ErrPack.o BDREF.o bits_for_disort_f77.o DISORT.o dtridgl.o dsolver.o gfluxi.o setup_RT_mod.o 

marv.o: sizes_mod.o  define_types_mod.o common_arrays_mod.o phys_const_mod.o atmos_ops_mod.o tristan.o gas_opacity_mod.o clouds_mod.o RDI1MACH.o LINPAK.o ErrPack.o BDREF.o bits_for_disort_f77.o DISORT.o dtridgl.o dsolver.o gfluxi.o  setup_RT_mod.o  main_mod.o





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

# Some rules for building f77
%.o: %.f
	$(FC) $(F77FLAGS) -c $<

#f90mods:
#	$(FC) $(FCFLAGS) -c f90wrap_*.f90

f77mods:
	$(FC) $(F77FLAGS) -c *.f

# now for python wrap
#f90wrap:
#	f90wrap -m forwardmodel *.f90


libfile:
	$(FC) -fPIC -shared -O2 *.o -o libmarvin.so
#	$(FC) -fPIC -shared -O2 *.o -o libmarvin.so

pysig:
	f2py -m forwardmodel -h forwardmodel.pyf sizes_mod.f90 marv.f90

pymod:
	f2py --fcompiler=gfortran --f90flags="-O2" -I/usr/include -L/usr/local/lib -c libmarvin.so forwardmodel.pyf marv.f90
#	f2py --fcompiler=intelem --f90flags="-O2" -I/usr/include -L/usr/local/lib -c libmarvin.so forwardmodel.pyf marv.f90

ciasig:
	f2py -m ciamod -h ciamod.pyf sizes_mod.f90 read_cia.f90

ciamod:
	f2py --fcompiler=gfortran --f90flags="-O2" -I/usr/include -L/usr/local/lib -c ciamod.pyf sizes_mod.f90 read_cia.f90
#	f2py --fcompiler=intelem --f90flags="-O2" -I/usr/include -L/usr/local/lib -c ciamod.pyf read_cia.f90

bbconvsig:
	f2py -m bbconv -h bbconv.pyf bbconv.f90

bbconv:
	f2py --fcompiler=gfortran --f90flags="-O2" -I/usr/include -L/usr/local/lib -c bbconv.pyf bbconv.f90

cloudpostsig:
	f2py -m cloudpost -h cloudpost.pyf cloudpost.f90

cloudpost:
	f2py --fcompiler=gfortran --f90flags="-O2 -fPIC -frecord-marker=4 -fbounds-check " -I/usr/include -L/usr/local/lib -c cloudpost.pyf cloudpost.f90

# Utility targets

.PHONY: clean 

clean:
	rm -f *.o *.mod *.MOD f90wrap*  *.pyc *.pyf *.so
