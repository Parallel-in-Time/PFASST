# - Find Easyloggingpp
# Find the native Easyloggingpp includes
#
# Easyloggingpp_INCLUDE_PATH - where to find the single header file
# Easyloggingpp_FOUND - True if Easyloggingpp found.

if(Easyloggingpp_INCLUDE_PATH)
    # Already in cache, be silent
    set(Easyloggingpp_FIND_QUIETLY TRUE)
endif(Easyloggingpp_INCLUDE_PATH)

find_path(Easyloggingpp_INCLUDE_PATH
    easylogging++.h
    PATH_SUFFIXES easyloggingpp
)

# handle the QUIETLY and REQUIRED arguments and set Easyloggingpp_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Easyloggingpp
    DEFAULT_MSG
    Easyloggingpp_INCLUDE_PATH
)

mark_as_advanced(Easyloggingpp_INCLUDE_PATH)
