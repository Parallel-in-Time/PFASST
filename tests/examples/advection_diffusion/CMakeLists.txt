# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(TESTS
    test_advection_diffusion
)

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
