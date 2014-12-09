# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)

set(TESTS
    test_boris_particle_util
#     test_boris_physics
    test_boris_sweeper
)

foreach(test ${TESTS})
    message(STATUS "  ${test}")
    add_executable(${test} ${test}.cpp)
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${test} ${pfasst_DEPENDEND_TARGETS})
    endif()
    if(${pfasst_TESTS_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${test} ${pfasst_TESTS_DEPENDEND_TARGETS})
    endif()
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${TESTS_3rdparty_DEPENDEND_LIBS}
        ${pfasst_DEPENDEND_LIBS}
    )
    if(pfasst_WITH_GCC_PROF AND ${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
        set_target_properties(${test}
            PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs"
                       LINK_FLAGS "-fprofile-arcs"
        )
    endif()
    add_test(NAME ${test}
        COMMAND ${CMAKE_BINARY_DIR}/tests/examples/boris/${test} --gtest_output=xml:${test}_out.xml
    )
endforeach(test)

target_link_libraries(test_boris_sweeper simple_physics_solver)
