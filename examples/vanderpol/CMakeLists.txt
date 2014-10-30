set(all_example_runs ${all_example_runs})
set(all_vanderpol_example_runs)

message(STATUS "  vanderpol")
include_directories(
    ${pfasst_INCLUDES}
    ${3rdparty_INCLUDES}
)

set(vanderpol_examples
    vdp_sdc
)

foreach(example ${vanderpol_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    target_link_libraries(${example}
        ${3rdparty_DEPENDEND_LIBS}
    )
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/vanderpol
    )
endforeach(example)

message(STATUS "    creating make targets to run examples")
foreach(example ${vanderpol_examples})
    add_custom_target(run_example_vanderpol_${example}
        COMMAND ${example}
        DEPENDS ${example}
        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/vanderpol
        COMMENT "Running van der Pol example: ${example}" VERBATIM
    )
    list(APPEND all_scalar_example_runs run_example_vanderpol_${example})
    message(STATUS "      - run_example_vanderpol_${example}")
endforeach(example)

add_custom_target(run_example_vanderpol_all
  DEPENDS ${all_vanderpol_example_runs}
  COMMENT "Running all van der Pol examples" VERBATIM
)
message(STATUS "    - run_example_vanderpol_all")

set(all_example_runs ${all_example_runs} ${all_vanderpol_example_runs} PARENT_SCOPE)
