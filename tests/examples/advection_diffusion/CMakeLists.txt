# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${TESTS_3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

if(${pfasst_WITH_MPI})
    set(MPI_TESTS
        test_mpi_advection_diffusion
        test_mpi_advection_diffusion_conv
    )

    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(test ${MPI_TESTS})
        message(STATUS "  ${test}")
        add_executable(${test} ${test}.cpp)
        if(NOT FFTW_FOUND)
            add_dependencies(${test} fftw3)
        endif()
        if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
            add_dependencies(${test} ${pfasst_DEPENDEND_TARGETS})
        endif()
        if(${pfasst_TESTS_NUM_DEPENDEND_TARGETS} GREATER 0)
            add_dependencies(${test} ${pfasst_TESTS_DEPENDEND_TARGETS})
        endif()
        if(MPI_COMPILE_FLAGS)
            if(pfasst_WITH_GCC_PROF AND ${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
                set_target_properties(${test}
                    PROPERTIES COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
                               LINK_FLAGS "-fprofile-arcs"
                )
            endif()
        endif()
        target_link_libraries(${test}
            ${3rdparty_DEPENDEND_LIBS}
            ${TESTS_3rdparty_DEPENDEND_LIBS}
            ${FFTW_LIBRARIES}
            ${pfasst_DEPENDEND_LIBS}
        )
        add_test(NAME ${test}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_BINARY_DIR}/tests/examples/advection_diffusion/${test} --gtest_output=xml:${test}_out.xml
        )
    endforeach(test)
endif()

set(TESTS
    test_advection_diffusion
    test_advection_diffusion_conv
)

foreach(test ${TESTS})
    message(STATUS "  ${test}")
    add_executable(${test} ${test}.cpp)
    if(NOT FFTW_FOUND)
        add_dependencies(${test} fftw3)
    endif()
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${test} ${pfasst_DEPENDEND_TARGETS})
    endif()
    if(${pfasst_TESTS_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${test} ${pfasst_TESTS_DEPENDEND_TARGETS})
    endif()
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${TESTS_3rdparty_DEPENDEND_LIBS}
        ${FFTW_LIBRARIES}
        ${pfasst_DEPENDEND_LIBS}
    )
    if(pfasst_WITH_GCC_PROF AND ${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
        set_target_properties(${test}
            PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs"
                       LINK_FLAGS "-fprofile-arcs"
        )
    endif()
    add_test(NAME ${test}
        COMMAND ${CMAKE_BINARY_DIR}/tests/examples/advection_diffusion/${test} --gtest_output=xml:${test}_out.xml
    )
endforeach(test)
