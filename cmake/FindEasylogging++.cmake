# - Find Easylogging++
# Find the native Easylogging++ includes
#
# Easylogging++_INCLUDE_PATH - where to find the single header file
# EASYLOGGING++_FOUND - True if Easyloggingpp found.

if(Easylogging++_INCLUDE_PATH)
    # Already in cache, be silent
    set(Easylogging++_FIND_QUIETLY TRUE)
endif(Easylogging++_INCLUDE_PATH)

find_path(Easylogging++_INCLUDE_PATH
    easylogging++.h
    PATH_SUFFIXES easylogging++
)

# handle the QUIETLY and REQUIRED arguments and set EASYLOGGING++_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Easylogging++
    DEFAULT_MSG
    Easylogging++_INCLUDE_PATH
)

mark_as_advanced(Easylogging++_INCLUDE_PATH)
