# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)

set(TESTS
    test_boris_particle
    test_boris_physics
    test_boris_sweeper
)

foreach(test ${TESTS})
    message(STATUS "  ${test}")
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} googlemock)
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${pfasst_DEPENDED_LIBS}
    )
    set_target_properties(${test}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/boris
    )
    add_test(NAME ${test}
        COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/examples/boris/${test} --gtest_output=xml:${test}_out.xml
    )
endforeach(test)
