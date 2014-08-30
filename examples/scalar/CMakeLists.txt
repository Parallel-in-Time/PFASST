set(all_example_runs ${all_example_runs})
set(all_scalar_example_runs)

message(STATUS "  scalar")
include_directories(
    ${pfasst_INCLUDES}
    ${3rdparty_INCLUDES}
)

set(scalar_examples
    scalar_sdc
)

foreach(example ${scalar_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/scalar
    )
endforeach(example)

message(STATUS "    creating make targets to run examples")
foreach(example ${scalar_examples})
    add_custom_target(run_example_scalar_${example}
        COMMAND ${example}
        DEPENDS ${example}
        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/scalar
        COMMENT "Running Scalar example: ${example}" VERBATIM
    )
    list(APPEND all_scalar_example_runs run_example_scalar_${example})
    message(STATUS "      - run_example_scalar_${example}")
endforeach(example)

add_custom_target(run_example_scalar_all
  DEPENDS ${all_scalar_example_runs}
  COMMENT "Running all Scalar examples" VERBATIM
)
message(STATUS "    - run_example_scalar_all")

set(all_example_runs ${all_example_runs} ${all_scalar_example_runs} PARENT_SCOPE)
