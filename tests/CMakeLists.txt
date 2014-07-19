# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${fftw3_INCLUDES}
    ${pfasst_INCLUDES}
)

set(TESTS
    test-quadrature
    test-advection-diffusion
)

foreach(test ${TESTS})
    message(STATUS "  ${test}")
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} googlemock fftw3)
    target_link_libraries(${test}
        ${3rdparty_DEPENDEND_LIBS}
        ${fftw3_LIBS}
        ${pfasst_DEPENDED_LIBS}
    )
    set_target_properties(${test}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests
    )
    add_test(NAME ${test}
        COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/${test} --gtest_output=xml:${test}_out.xml
    )
endforeach(test)
