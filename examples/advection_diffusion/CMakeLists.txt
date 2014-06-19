message(STATUS "  advection_diffusion")
include_directories(
    ${3rdparty_INCLUDES}
    ${fftw3_INCLUDES}
    ${pfasst_INCLUDES}
)

set(ex1_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ex1.cpp
)

add_executable(ex1 ${ex1_SOURCES})
add_dependencies(ex1 fftw3)
target_link_libraries(ex1
    ${fftw3_LIBS}
)
set_target_properties(ex1
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples
)

set(ex2_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ex2.cpp
)
add_executable(ex2 ${ex2_SOURCES})
add_dependencies(ex2 fftw3)
target_link_libraries(ex2
    ${fftw3_LIBS}
)
set_target_properties(ex2
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples
)

set(ex3_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ex3.cpp
)
add_executable(ex3 ${ex3_SOURCES})
add_dependencies(ex3 fftw3)
target_link_libraries(ex3
    ${fftw3_LIBS}
)
set_target_properties(ex3
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples
)
