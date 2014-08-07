set(all_example_runs ${all_example_runs})
set(all_boris_example_runs)

message(STATUS "  boris")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)

set(boris_examples
    boris_sdc
)

foreach(example ${boris_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/boris
    )
endforeach(example)


message(STATUS "    creating make targets to run examples")
foreach(example ${boris_examples})
    add_custom_target(run_example_boris_${example}
        COMMAND ${example}
        DEPENDS ${example}
        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/boris
        COMMENT "Running Boris example: ${example}" VERBATIM
    )
    list(APPEND all_boris_example_runs run_example_boris_${example})
    message(STATUS "      - run_example_boris_${example}")
endforeach(example)


add_custom_target(run_example_boris_all
  DEPENDS ${all_boris_example_runs}
  COMMENT "Running all Boris examples" VERBATIM
)
message(STATUS "    - run_example_boris_all")

set(all_example_runs ${all_boris_example_runs} PARENT_SCOPE)
