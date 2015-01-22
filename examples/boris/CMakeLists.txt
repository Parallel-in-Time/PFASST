set(examples_to_install ${examples_to_install})

message(STATUS "  boris")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)
add_subdirectory(bindings)

set(boris_examples
    boris_sdc
    boris_mlsdc
)

foreach(example ${boris_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    target_link_libraries(${example}
        ${3rdparty_DEPENDEND_LIBS}
        simple_physics_solver
    )
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
    endif()
    add_dependencies(${example} simple_physics_solver)
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${boris_examples} PARENT_SCOPE)
