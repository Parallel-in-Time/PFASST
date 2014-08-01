set(all_example_runs ${all_example_runs})
set(all_advec_example_runs)

message(STATUS "  advection_diffusion")
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
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
    if(NOT FFTW_FOUND)
      add_dependencies(${example} fftw3)
    endif()
    target_link_libraries(${example}
        ${FFTW_LIBRARIES}
    )
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/advection_diffusion
    )
endforeach(example)

message(STATUS "    creating make targets to run examples")
foreach(example ${advec_examples})
    add_custom_target(run_example_advec_${example}
        COMMAND ${example}
        DEPENDS ${example}
        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/advection_diffusion
        COMMENT "Running Advection-Diffusion example: ${example}" VERBATIM
    )
    list(APPEND all_advec_example_runs run_example_advec_${example})
    message(STATUS "      - run_example_advec_${example}")
endforeach(example)

if(${pfasst_WITH_MPI})
    message(STATUS "    creating make targets to run MPI examples")
    foreach(example ${advec_mpi_examples})
        add_custom_target(run_example_advec_${example}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${example} ${MPIEXEC_POSTFLAGS}
            DEPENDS ${example}
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/advection_diffusion
            COMMENT "Running Advection-Diffusion example: ${example}" VERBATIM
        )
        list(APPEND all_advec_example_runs run_example_advec_${example})
        message(STATUS "      - run_example_advec_${example}")
    endforeach(example)
endif()

add_custom_target(run_example_advec_all
  DEPENDS ${all_advec_example_runs}
  COMMENT "Running all Advection-Diffusion examples" VERBATIM
)
message(STATUS "    - run_example_advec_all")

set(all_example_runs ${all_advec_example_runs} PARENT_SCOPE)
