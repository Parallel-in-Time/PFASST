# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
message(STATUS "--------------------------------------------------------------------------------")

if(pfasst_BUILD_EXAMPLES)
    message(STATUS "FFTW3")
    set(fftw3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3")
    set(fftw3_INSTALL_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3-install")
    
    ExternalProject_Add(
        fftw3
        LIST_SEPARATOR " "
        URL http://fftw.org/fftw-3.3.4.tar.gz
        URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
        TIMEOUT 10
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ${fftw3_SOURCE_DIR}/configure CC=${CMAKE_C_COMPILER} --prefix=${fftw3_INSTALL_DIR} CXX=${CMAKE_CXX_COMPILER}
        BUILD_COMMAND make
        TEST_COMMAND ""
        INSTALL_DIR ${fftw3_SOURCE_DIR}-install
        INSTALL_COMMAND make install
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
    )
    # Specify include dir
    set(fftw3_INCLUDES ${fftw3_INSTALL_DIR}/include PARENT_SCOPE)
    set(fftw3_LIBS ${fftw3_INSTALL_DIR}/lib64/libfftw3.a PARENT_SCOPE)

    if(pfasst_BUILD_TESTS)
        message(STATUS "--------------------------------------------------------------------------------")
    endif()
endif()

if(pfasst_BUILD_TESTS)
    message(STATUS "Google Testing Framework (gtest & gmock)")

    # Add gmock
    ExternalProject_Add(
      googlemock
      SVN_REPOSITORY http://googlemock.googlecode.com/svn/trunk/
      SVN_REVISION -r 449  # release 1.7.0
      TIMEOUT 10
      # Disable SVN update
      UPDATE_COMMAND ""
      PATCH_COMMAND ""
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
      -Dgtest_force_shared_crt=ON
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -Dgmock_build_tests=${gtest_BUILD_TESTS}
      # Disable install step
      INSTALL_COMMAND ""
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
    )

    # Specify include dir
    ExternalProject_Get_Property(googlemock source_dir)
    list(APPEND 3rdparty_INCLUDES ${source_dir}/include ${source_dir}/gtest/include)

    ExternalProject_Get_Property(googlemock binary_dir)
    set(Suffix ".a")
    set(Pthread "-pthread")

    list(APPEND 3rdparty_DEPENDEND_LIBS ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${Suffix})
    list(APPEND 3rdparty_DEPENDEND_LIBS ${Pthread})
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
