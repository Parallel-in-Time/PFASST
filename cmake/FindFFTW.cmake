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

find_path(FFTW3_INCLUDE_PATH fftw3.h
    HINTS ${FFTW3_ROOT}
        ${FFTW3_INCLUDE}
        $ENV{FFTW3_ROOT}
        $ENV{FFTW3_INCLUDE}
)

find_library(FFTW3_LIBRARIES NAMES fftw3
    HINTS ${FFTW3_ROOT}
        ${FFTW3_LIB}
        ${FFTW3_LIBS}
        $ENV{FFTW3_ROOT}
        $ENV{FFTW3_LIB}
        $ENV{FFTW3_LIBS}
)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_PATH)

mark_as_advanced(FFTW3_LIBRARIES FFTW3_INCLUDE_PATH)
