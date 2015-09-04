set(examples_to_install ${examples_to_install})

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
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
    endif()
    target_link_libraries(${example}
        ${pfasst_DEPENDEND_LIBS}
        ${3rdparty_DEPENDEND_LIBS}
    )
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${scalar_examples} PARENT_SCOPE)
