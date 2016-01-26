function(add_to_string_list in_string out_string)
#   message(STATUS "Adding '${ARGN}' to '${in_string}'")
    foreach(element ${ARGN})
#       message(STATUS "  processing: '${element}'")
        if(NOT "${element}" STREQUAL "")
            # escape regex commands in new element
            string(REPLACE "\\" "\\\\" elem_escaped "${element}")
            string(REGEX REPLACE "([][.?*+|()$^-])" "\\\\\\1" elem_escaped "${elem_escaped}")

            # only append if not already in string
            if("${in_string}" MATCHES "${elem_escaped}")
#               message(STATUS "    '${element}' found as '${elem_escaped}'")
            else()
                set(in_string "${in_string} ${element}")
#               message(STATUS "    '${element}' appended")
            endif()
        endif()
    endforeach(element)
#   message(STATUS "Full string is now: '${in_string}'")
    set(${out_string} ${in_string} PARENT_SCOPE)
endfunction(add_to_string_list)

macro(msg_not_installed dependency_name)
    message("!!!!")
    message("!  ${dependency_name} will not be installed with PFASST++.")
    message("!  When using PFASST++ with your own code, please make sure to also add")
    message("!  ${dependency_name} manually to your project.")
    message("!!!!")
endmacro(msg_not_installed)

macro(update_site_config)
    set(WARNING_COMMENT "/*\n * DO NOT ALTER THIS FILE\n *\n * It will get rewritten on CMake's next run\n */")
    execute_process(COMMAND git describe --dirty
                    OUTPUT_VARIABLE pfasst_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    configure_file(
        "${pfasst_SOURCE_DIR}/cmake/site_config.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/pfasst/site_config.hpp"
    )
    add_definitions(-DPFASST_USE_LOCAL_CONFIG)
endmacro(update_site_config)
