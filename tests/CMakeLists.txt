# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)

set(TESTS
    test-quadrature
)

foreach(test ${TESTS})
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} googlemock)
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${pfasst_DEPENDED_LIBS}
    )
    set_target_properties(${test} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${pfasst_BINARY_DIR}/tests)
    add_test(NAME ${test}
        COMMAND ${CMAKE_BINARY_DIR}/tests/${test}
    )
endforeach(test)
