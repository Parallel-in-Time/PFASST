# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
message(STATUS "--------------------------------------------------------------------------------")

message(STATUS "Boost")

if(pfasst_BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS OFF)
    set(pfasst_BUILD_SHARED_LIBS ON)
    set(BUILD_SHARED_LIBS ON)
else()
    set(Boost_USE_STATIC_LIBS ON)
    set(pfasst_BUILD_SHARED_LIBS OFF)
    set(BUILD_SHARED_LIBS OFF)
endif()
set(Boost_USE_MULTITHREADED ON)

find_package(Boost REQUIRED)
list(APPEND 3rdparty_INCLUDES ${Boost_INCLUDE_DIRS})


if(pfasst_BUILD_EXAMPLES)
    find_package(FFTW)

    if(NOT FFTW_FOUND)
      message(STATUS "FFTW3 not found on your system.")
      message(STATUS "Going to download and compile it automatically.")
      set(fftw3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3")
      set(fftw3_INSTALL_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3-install")

      ExternalProject_Add(
          fftw3
          LIST_SEPARATOR " "
          URL http://fftw.org/fftw-3.3.4.tar.gz
          URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
          TIMEOUT 30
          UPDATE_COMMAND ""
          PATCH_COMMAND ""
          BUILD_IN_SOURCE ON
          CONFIGURE_COMMAND ${fftw3_SOURCE_DIR}/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS} --prefix=${fftw3_INSTALL_DIR} --libdir=${fftw3_INSTALL_DIR}/lib
          BUILD_COMMAND make
          TEST_COMMAND ""
          INSTALL_DIR ${fftw3_SOURCE_DIR}-install
          INSTALL_COMMAND make install
          LOG_DOWNLOAD ON
          LOG_CONFIGURE ON
          LOG_BUILD ON
      )
      # Specify include dir
      set(FFTW_INCLUDE_PATH ${fftw3_INSTALL_DIR}/include PARENT_SCOPE)
      set(FFTW_LIBRARIES ${fftw3_INSTALL_DIR}/lib/libfftw3.a PARENT_SCOPE)
    else()
      set(FFTW_INCLUDE_PATH ${FFTW_INCLUDE_PATH} PARENT_SCOPE)
      set(FFTW_LIBRARIES ${FFTW_LIBRARIES} PARENT_SCOPE)
    endif()
    set(FFTW_FOUND ${FFTW_FOUND} PARENT_SCOPE)

    if(pfasst_BUILD_TESTS)
        message(STATUS "--------------------------------------------------------------------------------")
    endif()
endif()

if(pfasst_BUILD_TESTS)
    message(STATUS "Google Testing Framework (gtest & gmock)")

    # Add gmock
    ExternalProject_Add(
      googlemock
      URL http://googlemock.googlecode.com/files/gmock-1.7.0.zip
      URL_MD5 073b984d8798ea1594f5e44d85b20d66
      TIMEOUT 30
      UPDATE_COMMAND ""
      PATCH_COMMAND ""
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_C_COMPILE=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
      -DGTEST_USE_OWN_TR1_TUPLE=ON
      -Dgtest_force_shared_crt=ON
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
