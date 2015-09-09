# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS})
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS})

set(3rdparty_DOWNLOAD_TIMEOUT "60")

message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Boost")

set(Boost_USE_MULTITHREADED ON)
if(${pfasst_BUILD_SHARED_LIBS})
    set(Boost_USE_STATIC_LIBS OFF)
else()
    set(Boost_USE_STATIC_LIBS ON)
endif()
set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0" ${Boost_ADDITIONAL_VERSIONS})

if(${compiler_version_available} AND ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    set(Boost_COMPILER "-clang${CMAKE_CXX_COMPILER_VERSION_MAJOR}${CMAKE_CXX_COMPILER_VERSION_MINOR}")
endif()

find_package(Boost 1.53.0 COMPONENTS program_options REQUIRED)
set_package_properties(Boost
    PROPERTIES
      DESCRIPTION "C++ libraries"
      URL "http://boost.org/"
      TYPE REQUIRED
      PURPOSE "command line parameter parsing"
)
message(STATUS "Boost include path: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library path: ${Boost_LIBRARY_DIRS}")
list(APPEND 3rdparty_INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND 3rdparty_DEPENDEND_LIBS ${Boost_LIBRARIES})


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Eigen3")

find_package(Eigen3)
set_package_properties(Eigen3
    PROPERTIES
      DESCRIPTION "linear algebra library"
      URL "http://eigen.tuxfamily.org/"
      TYPE OPTIONAL
      PURPOSE "matrix and vector types and arithmetics"
)

if(NOT EIGEN3_FOUND)
    set(EIGEN3_SOURCE_URL "http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2")
    set(EIGEN3_SOURCE_MD5 "fc2e814ae449d16b331f7e1f4e272bd3")

    message(STATUS "Eigen3 not found on your system")
    message(STATUS " going to download sources from:")
    message(STATUS "  ${EIGEN3_SOURCE_URL}")
    set(Eigen3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/Eigen3")

    ExternalProject_Add(
        Eigen3
        LIST_SEPARATOR " "
        URL ${EIGEN3_SOURCE_URL}
        URL_MD5 ${EIGEN3_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        TEST_COMMAND ""
        INSTALL_DIR ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE OFF
        LOG_BUILD OFF
    )

    # Specify include dir
    set(Eigen3_INCLUDE_PATH "${Eigen3_SOURCE_DIR}")
    set(Eigen3_FOUND TRUE)

    list(APPEND pfasst_DEPENDEND_TARGETS Eigen3)

    msg_not_installed("Eigen3")
endif()

list(APPEND 3rdparty_INCLUDES "${Eigen3_INCLUDE_PATH}")


if(pfasst_BUILD_EXAMPLES)
    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "FFTW3")

    find_package(FFTW)
    set_package_properties(FFTW
        PROPERTIES
          DESCRIPTION "FFTW3"
          URL "http://fftw.org/"
          TYPE OPTIONAL
          PURPOSE "fast Fourier transformation"
    )

    if(NOT FFTW_FOUND OR NOT FFTW_LIBRARIES)
        set(FFTW3_SOURCE_URL "http://fftw.org/fftw-3.3.4.tar.gz")
        set(FFTW3_SOURCE_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

        message(STATUS "FFTW3 not found on your system")
        message(STATUS " going to download and compile it from sources retrieved from:")
        message(STATUS "  ${FFTW3_SOURCE_URL}")
        set(fftw3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3")
        set(fftw3_INSTALL_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3-install")

        ExternalProject_Add(
            fftw3
            LIST_SEPARATOR " "
            URL ${FFTW3_SOURCE_URL}
            URL_MD5 ${FFTW3_SOURCE_MD5}
            TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}
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
        set(FFTW_INCLUDE_PATH ${fftw3_INSTALL_DIR}/include)
        set(FFTW_LIBRARIES ${fftw3_INSTALL_DIR}/lib/libfftw3.a)

        set(FFTW_FOUND TRUE)

        list(APPEND pfasst_DEPENDEND_LIBS FFTW)

        msg_not_installed("FFTW3")
    else()
        set(FFTW_INCLUDE_PATH ${FFTW_INCLUDE_PATH} PARENT_SCOPE)
        set(FFTW_LIBRARIES ${FFTW_LIBRARIES} PARENT_SCOPE)
    endif()
    set(FFTW_FOUND ${FFTW_FOUND} PARENT_SCOPE)
endif()

if(pfasst_BUILD_TESTS)
    set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES})
    set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS})

    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "Google Testing Framework (gtest & gmock)")

    find_package(GMock)
    set_package_properties(GMock
        PROPERTIES
          DESCRIPTION "GMock"
          URL "https://code.google.com/p/googlemock/"
          TYPE OPTIONAL
          PURPOSE "Google testing and mocking framework"
    )

    if(NOT GMOCK_FOUND)
        set(GMOCK_SOURCE_URL "http://googlemock.googlecode.com/files/gmock-1.7.0.zip")
        set(GMOCK_SOURCE_MD5 "073b984d8798ea1594f5e44d85b20d66")

        message(STATUS " going to download and compile it from sources retrieved from:")
        message(STATUS "   ${GMOCK_SOURCE_URL}")

        # Add gmock
        ExternalProject_Add(
          googlemock
          URL ${GMOCK_SOURCE_URL}
          URL_MD5 ${GMOCK_SOURCE_MD5}
          TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}
          UPDATE_COMMAND ""
          PATCH_COMMAND ""
          CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILE=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DGTEST_USE_OWN_TR1_TUPLE=ON
            -Dgtest_force_shared_crt=ON
            -Dgmock_build_tests=OFF
          # Disable install step
          INSTALL_COMMAND ""
          # Wrap download, configure and build steps in a script to log output
          LOG_DOWNLOAD ON
          LOG_CONFIGURE ON
          LOG_BUILD ON
        )

        # Specify include dir
        ExternalProject_Get_Property(googlemock source_dir)
        list(APPEND TESTS_3rdparty_INCLUDES ${source_dir}/include ${source_dir}/gtest/include)

        ExternalProject_Get_Property(googlemock binary_dir)
        set(Suffix ".a")
        set(Pthread "-pthread")
        list(APPEND TESTS_3rdparty_DEPENDEND_LIBS ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${Suffix})
        list(APPEND TESTS_3rdparty_DEPENDEND_LIBS ${Pthread})
        list(APPEND pfasst_TESTS_DEPENDEND_TARGETS googlemock)
    else()
        list(APPEND TESTS_3rdparty_DEPENDEND_LIBS ${GMOCK_LIBRARIES} "-pthread")
        list(APPEND TESTS_3rdparty_INCLUDES ${GMOCK_INCLUDE_DIRS})
    endif(NOT GMOCK_FOUND)
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS} PARENT_SCOPE)
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS} PARENT_SCOPE)
set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES} PARENT_SCOPE)
