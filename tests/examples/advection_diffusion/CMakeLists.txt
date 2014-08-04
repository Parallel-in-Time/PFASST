# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(TESTS
    test_advection_diffusion
)

if(${pfasst_WITH_MPI})
    set(MPI_TESTS
        test_mpi_advection_diffusion
    )
endif()

foreach(test ${TESTS})
    message(STATUS "  ${test}")
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} googlemock)
    if(NOT FFTW_FOUND)
      add_dependencies(${test} googlemock fftw3)
    endif()
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${FFTW_LIBRARIES}
        ${pfasst_DEPENDED_LIBS}
    )
    set_target_properties(${test}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/advection_diffusion
    )
    add_test(NAME ${test}
        COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/advection_diffusion/${test} --gtest_output=xml:${test}_out.xml
    )
endforeach(test)

if(${pfasst_WITH_MPI})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(test ${MPI_TESTS})
        message(STATUS "  ${test}")
        add_executable(${test} ${test}.cpp)
        add_dependencies(${test} googlemock)
        if(NOT FFTW_FOUND)
          add_dependencies(${test} googlemock fftw3)
        endif()
	if(MPI_COMPILE_FLAGS)
            set_target_properties(${test} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
        if(MPI_LINK_FLAGS)
            set_target_properties(${test} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
        target_link_libraries(${test}
            ${3rdparty_DEPENDEND_LIBS}
            ${FFTW_LIBRARIES}
            ${pfasst_DEPENDED_LIBS}
            ${MPI_CXX_LIBRARIES}
        )
        set_target_properties(${test}
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/advection_diffusion
        )
        add_test(NAME ${test}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/advection_diffusion/${test} ${MPIEXEC_POSTFLAGS}
        )
    endforeach(test)
endif()
