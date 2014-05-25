# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})

if(pfasst_BUILD_TESTS)
    message(STATUS "3rd party: Google Testing Framework (gtest & gmock)")
    message(STATUS "--------------------------------------------------------------------------------")

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
      -Dgmock_build_tests=ON
      # Disable install step
      INSTALL_COMMAND ""
      # Run tests
      TEST_BEFORE_INSTALL ON
      TEST_AFTER_INSTALL ON
      TEST_COMMAND make test
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

    message(STATUS "--------------------------------------------------------------------------------")
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
