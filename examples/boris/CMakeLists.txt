set(examples_to_install ${examples_to_install})

message(STATUS "  boris")

option(pfasst_BORIS_SAME_LEVELS "Use identical levels with MLSDC" OFF)

include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
)
add_subdirectory(bindings)

if(${pfasst_WITH_MPI})
    set(boris_examples
        boris_pfasst
    )
else()
    set(boris_examples
        boris_sdc
        boris_mlsdc
    )
endif()

foreach(example ${boris_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    target_link_libraries(${example}
        ${pfasst_DEPENDEND_LIBS}
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

if(pfasst_BORIS_SAME_LEVELS)
    foreach(example "boris_mlsdc" "boris_pfasst")
        if(TARGET ${example})
            get_target_property(definitions ${example} COMPILE_DEFINITIONS)
            if(${definitions})
                list(APPEND definitions "BORIS_SAME_LEVELS")
            else()
                set(definitions "BORIS_SAME_LEVELS")
            endif()
            set_target_properties(${example} PROPERTIES COMPILE_DEFINITIONS ${definitions})
        endif()
    endforeach(example)
endif()


add_feature_info(Dummy-Boris-MLSDC
    pfasst_BORIS_SAME_LEVELS
    "use identical levels for Boris-MLSDC"
)

set(examples_to_install ${examples_to_install} ${boris_examples} PARENT_SCOPE)
