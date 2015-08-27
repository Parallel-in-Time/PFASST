set(examples_to_install ${examples_to_install})

add_subdirectory(advec_diff)
add_subdirectory(heat1d)

set(examples_to_install ${examples_to_install} PARENT_SCOPE)
