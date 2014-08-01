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

find_path(FFTW_INCLUDE_PATH fftw3.h)

find_library(FFTW_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_PATH)

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDE_PATH)
