set(examples_to_install ${examples_to_install})

message(STATUS "  boris")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)
add_subdirectory(bindings)

set(boris_examples
    boris_sdc
)

foreach(example ${boris_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    target_link_libraries(${example}
        ${3rdparty_DEPENDEND_LIBS}
        simple_physics_solver
    )
    add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS} simple_physics_solver)
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${boris_examples} PARENT_SCOPE)
