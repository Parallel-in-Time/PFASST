if(SCOREP_FIND_QUIETLY)
    set(SCOREP_FIND_QUIETLY TRUE)
endif()

find_program(
    SCOREP_EXECUTABLE
    NAMES scorep
    DOC "Score-P scalable performance measurement for parallel codes"
)

if(SCOREP_EXECUTABLE)
    set(SCOREP_FOUND "YES")
else()
    set(SCOREP_FOUND "NO")
endif()

find_package_handle_standard_args(SCOREP FOUND_VAR SCOREP_FOUND REQUIRED_VARS SCOREP_EXECUTABLE)

mark_as_advanced(SCOREP_EXECUTABLE)
