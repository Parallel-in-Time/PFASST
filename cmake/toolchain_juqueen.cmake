# Toolchain File for JUQUEEN @ Juelich Supercomputing Center
# BG/Q system
#
# Remark: IBM XL/C++ compiler is not supported by PFASST++ due to lack of C++11 support
#
# IMPORTANT!!
# Make sure you have loaded only the following modules:
#  - gcc/4.8.1
#  - fftw3/3.3.3
#  - python3/3.4.2
#
# As well make sure you have built boost::program_options according to the instructions from the
# PFASST++ documentation on JUQUEEN.
#

message(STATUS "Please make sure you have loaded at least gcc/4.8.1")

set(CMAKE_C_COMPILER /bgsys/local/gcc/4.8.1/bin/mpigcc)
set(CMAKE_CXX_COMPILER /bgsys/local/gcc/4.8.1/bin/mpig++)
set(CMAKE_Fortran_COMPILER /bgsys/local/gcc/4.8.1/bin/mpigfortran)

set(MPI_ROOT /bgsys/local/gcc/4.8.1)
set(MPI_C_COMPILER ${MPI_ROOT}/bin/mpigcc)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpig++)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpigfortran)

set(PYTHON_EXECUTABLE /bgsys/local/python3/3.4.2/bin/python3)
