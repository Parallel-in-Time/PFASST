message(STATUS "  scalar")
include_directories(
    ${pfasst_INCLUDES}
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
