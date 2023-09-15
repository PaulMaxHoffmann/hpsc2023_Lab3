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
#
# In  ubuntu, to install OpenMPI:
#
# sudo  apt-get update               <-- In Windows Ubuntu, needed this first
# sudo  apt-get install  openmpi-bin
# 
# Be sure you have g++ installed, too!
#
# Once  installed, the build command below will work.  Note  that the
# std=c++11 is requied to use the timers.
#

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

main: main.cpp
	$(COMPILER) -fopenmp main.cpp -o main -std=c++11
#	g++ -fopenmp main.cpp -o main



