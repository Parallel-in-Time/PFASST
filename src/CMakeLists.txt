# add this source directory to PFASST's includes so headers are found here
list(APPEND pfasst_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# enable debug build if required
if(${CMAKE_BUILD_TYPE} MATCHES Debug OR pfasst_DEBUG)
    add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Wall")
    add_definitions(-DPFASST_DEBUG)
endif()

## set the sources for PFASST
## (uncomment the following lines if PFASST becomes not a header-only library one day)
# set(pfasst_SOURCES
# )

# include_directories(
#     ${3rdparty_INCLUDES}
#     ${pfasst_INCLUDES}
# )
# message("3rdparty_INCLUDES: ${3rdparty_INCLUDES}")
# 
# add_library(pfasst ${pfasst_SOURCES})
# target_link_libraries(pfasst
#   ${3rdparty_DEPENDEND_LIBS}
#   ${pfasst_DEPENDEND_LIBS}
# )

set(pfasst_INCLUDES ${pfasst_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_LIBS ${pfasst_DEPENDEND_LIBS} PARENT_SCOPE)
