# - Find FFTW
# Find the native FFTW3 includes and library
#
# FFTW_INCLUDE_PATH - where to find fftw3.h
# FFTW_LIBRARIES - List of libraries when using FFTW.
# FFTW_FOUND - True if FFTW found.

if(FFTW_INCLUDE_PATH)
    # Already in cache, be silent
    set(FFTW_FIND_QUIETLY TRUE)
endif(FFTW_INCLUDE_PATH)

find_path(FFTW_INCLUDE_PATH fftw3.h
    HINTS ${FFTW3_INCLUDE}
        ${FFTW3_INC}
        ${FFTW3_ROOT}/include
        $ENV{FFTW3_INCLUDE}
        $ENV{FFTW3_INC}
        $ENV{FFTW_INC} # Edison
        $ENV{FFTW3_ROOT}/include
)

find_library(FFTW_LIBRARIES NAMES fftw3 fftw3l fftw3_mpi fftw3l_mpi fftw3_omp fftw3l_omp
    HINTS ${FFTW3_LIB}
        ${FFTW3_ROOT}/lib64
        ${FFTW3_ROOT}/lib
        $ENV{FFTW3_LIB}
        $ENV{FFTW_DIR} # Edison
        $ENV{FFTW3_ROOT}/lib64
        $ENV{FFTW3_ROOT}/lib
)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_PATH)

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDE_PATH)
