set(all_example_runs)

add_subdirectory(advection_diffusion)
add_subdirectory(scalar)
add_subdirectory(vanderpol)

add_custom_target(run_example_all
    DEPENDS ${all_example_runs}
    COMMENT "Running all examples" VERBATIM
)
message(STATUS "  - run_example_all")
