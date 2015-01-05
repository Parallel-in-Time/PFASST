add_library(simple_physics_solver
    simple_physics_solver.cpp
)

if(pfasst_INSTALL_EXAMPLES)
    install(TARGETS simple_physics_solver
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif(pfasst_INSTALL_EXAMPLES)
