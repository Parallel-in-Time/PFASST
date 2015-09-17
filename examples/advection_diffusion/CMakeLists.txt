set(examples_to_install ${examples_to_install})

message(STATUS "  advection_diffusion")
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(advec_examples
    vanilla_sdc
    serial_mlsdc
)

set(advec_mpi_examples)

if(${pfasst_WITH_MPI})
    set(advec_mpi_examples
        mpi_pfasst
    )
    set(all_advec_examples ${advec_examples} ${advec_mpi_examples})
else()
    set(all_advec_examples ${advec_examples})
endif()

foreach(example ${all_advec_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    if(NOT FFTW_FOUND)
        add_dependencies(${example} fftw3)
    endif()
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
    endif()
    target_link_libraries(${example}
        ${pfasst_DEPENDEND_LIBS}
        ${3rdparty_DEPENDEND_LIBS}
        ${FFTW_LIBRARIES}
    )
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${all_advec_examples} PARENT_SCOPE)
