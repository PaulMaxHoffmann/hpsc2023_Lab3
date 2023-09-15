
# -----------------------------------------------------
# On local Ubuntu
# -----------------------------------------------------
#
# Be sure to use MPI compiler wrappers.
# Whenever compiling an MPI program, you should use the MPI wrappers:
#
#    C - mpicc
#    C++ - mpiCC, mpicxx, mpic++
#    FORTRAN - mpifort, mpif77, mpif90
#
# These wrappers do all of the dirty work for you of making sure that
# all of the appropriate compiler flags, libraries, include directories,
# library directories, etc. are included when you compile your program.


# -----------------------------------------------------
# On Summit:
# -----------------------------------------------------
#
# module load intel
# module load impi

# -----------------------------------------------------
# Set compiler based on platform
# -----------------------------------------------------

# Default compiler is Intel MPI compiler

COMPILER := mpicc
UBUNTU   := ubuntu

# Change if on ubuntu

ifeq ($(strip $(LOGNAME)),scott)
COMPILER := mpicxx
endif

# -----------------------------------------------------
# Make lookup
# -----------------------------------------------------


lookup: lookup.cpp 
	$(COMPILER) lookup.cpp -o lookup

lookup_mpi: lookup_mpi.cpp 
	$(COMPILER) lookup_mpi.cpp -o lookup_mpi

