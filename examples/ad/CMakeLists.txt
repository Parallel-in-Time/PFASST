message(STATUS "  advection")
include_directories(
    ${fftw3_INCLUDES}
    ${pfasst_INCLUDES}
)

set(advection_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

add_executable(advection ${advection_SOURCES})
add_dependencies(advection fftw3)
target_link_libraries(advection
    ${fftw3_LIBS}
)
set_target_properties(advection
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples
)
