message(STATUS "  advection_diffusion")
include_directories(
    ${3rdparty_INCLUDES}
    ${fftw3_INCLUDES}
    ${pfasst_INCLUDES}
)

set(advec_examples
    vanilla_sdc
    serial_mlsdc
    serial_mlsdc_autobuild
)

foreach(example ${advec_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    add_dependencies(${example} fftw3)
    target_link_libraries(${example}
        ${fftw3_LIBS}
    )
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/examples/advection_diffusion
    )
endforeach(example)
