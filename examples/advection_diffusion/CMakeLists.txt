message(STATUS "  advection_diffusion")
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(advec_examples
    vanilla_sdc
    serial_mlsdc
    serial_mlsdc_autobuild
)

foreach(example ${advec_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    if(NOT FFTW_FOUND)
      add_dependencies(${example} fftw3)
    endif()
    target_link_libraries(${example}
        ${FFTW_LIBRARIES}
    )
    set_target_properties(${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/advection_diffusion
    )
endforeach(example)
