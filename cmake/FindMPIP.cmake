# - Find MPIP
# Find the mpiP library
#
# MPIP_LIBRARIES - List of libraries when using MPIP.
# MPIP_FOUND - True if mpiP was found.

if(MPIP_LIBRARIES)
    # Already in cache, be silent
    set(MPIP_FIND_QUIETLY TRUE)
endif(MPIP_LIBRARIES)

find_library(MPIP_LIBRARIES NAMES mpiP mpip
    HINTS ${MPIP_ROOT}/lib
        ${MPIP_ROOT}/lib64
        $ENV{MPIP_LIB}
        $ENV{MPIP_ROOT}/lib
        $ENV{MPIP_ROOT}/lib64
        $ENV{MPIP_DIR}/lib
        $ENV{MPIP_DIR}/lib64
)
list(APPEND MPIP_LIBRARIES m bfd iberty unwind)

# handle the QUIETLY and REQUIRED arguments and set MPIP_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPIP DEFAULT_MSG MPIP_LIBRARIES)

mark_as_advanced(MPIP_LIBRARIES)
