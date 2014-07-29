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

set(advec_mpi_examples)
if(${pfasst_WITH_MPI})
    set(advec_mpi_examples
        mpi_pfasst
    )
endif()

set(all_advec_examples
    ${advec_examples}
    ${advec_mpi_examples}
)

foreach(example ${all_advec_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    add_dependencies(${example} fftw3)
    target_link_libraries(${example}
        ${fftw3_LIBS}
    )
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/advection_diffusion
    )
endforeach(example)
