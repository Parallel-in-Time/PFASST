message(STATUS "  advection_diffusion")
include_directories(
    ${fftw3_INCLUDES}
    ${pfasst_INCLUDES}
)

set(advection_diffusion_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/advection_diffusion.cpp
)

add_executable(advdiff ${advection_diffusion_SOURCES})
add_dependencies(advdiff fftw3)
target_link_libraries(advdiff
    ${fftw3_LIBS}
)
set_target_properties(advdiff
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples
)
