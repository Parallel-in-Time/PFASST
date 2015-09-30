# - Find Eigen3
# Find the native Eigen3 includes
#
# Eigen3_INCLUDE_PATH - where to find signature_of_eigen3_matrix_library
# Eigen3_FOUND - True if Eigen3 found.

if(Eigen3_INCLUDE_PATH)
    # Already in cache, be silent
    set(Eigen3_FIND_QUIETLY TRUE)
endif(Eigen3_INCLUDE_PATH)

find_path(Eigen3_INCLUDE_PATH
    signature_of_eigen3_matrix_library
    PATH_SUFFIXES eigen3 Eigen3 Eigen
    HINTS $ENV{EIGEN3_DIR}/include/eigen3 # Edison
)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3 DEFAULT_MSG Eigen3_INCLUDE_PATH)

mark_as_advanced(Eigen3_INCLUDE_PATH)
