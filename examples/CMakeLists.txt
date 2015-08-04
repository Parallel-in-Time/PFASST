set(examples_to_install ${examples_to_install})

add_subdirectory(advection_diffusion)

set(examples_to_install ${examples_to_install} PARENT_SCOPE)
