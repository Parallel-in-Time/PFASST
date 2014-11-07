# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS})

message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Boost")

set(Boost_USE_MULTITHREADED ON)
set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0" ${Boost_ADDITIONAL_VERSIONS})

if(${compiler_version_available} AND ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    set(Boost_COMPILER "-clang${CMAKE_CXX_COMPILER_VERSION_MAJOR}${CMAKE_CXX_COMPILER_VERSION_MINOR}")
endif()

if(${CMAKE_VERBOSE_MAKEFILE})
    set(Boost_DEBUG ON)
    set(Boost_DETAILED_FAILURE_MSG ON)
endif()

find_package(Boost COMPONENTS program_options REQUIRED)
message(STATUS "Boost include path: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library path: ${Boost_LIBRARY_DIRS}")
list(APPEND 3rdparty_INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND 3rdparty_DEPENDEND_LIBS ${Boost_LIBRARIES})


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Eigen3")

find_package(Eigen3)

if(NOT EIGEN3_FOUND)
    set(EIGEN3_SOURCE_URL "http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2")
    message(STATUS "Eigen3 not found on your system")
    message(STATUS "  going to download sources")
    message(STATUS "    download URL: ${EIGEN3_SOURCE_URL}")
    set(Eigen3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/Eigen3")

    ExternalProject_Add(
        Eigen3
        LIST_SEPARATOR " "
        URL ${EIGEN3_SOURCE_URL}
        URL_MD5 fc2e814ae449d16b331f7e1f4e272bd3
        TIMEOUT 30
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
    list(APPEND pfasst_DEPENDEND_TARGETS Eigen3)
endif()

list(APPEND 3rdparty_INCLUDES "${Eigen3_INCLUDE_PATH}")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "easyloggingpp")

set(EASYLOGGINGPP_SOURCE_URL "https://github.com/easylogging/easyloggingpp/releases/download/v9.75/easyloggingpp_v9.75.tar.gz")
# message(STATUS "easyloggingpp not found on your system")
message(STATUS "  going to download sources")
message(STATUS "    download URL: ${EASYLOGGINGPP_SOURCE_URL}")
set(easyloggingpp_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/easyloggingpp")

ExternalProject_Add(
    easyloggingpp
    LIST_SEPARATOR " "
    URL ${EASYLOGGINGPP_SOURCE_URL}
    URL_MD5 3e37462b30cb6ee49faeed6a94d1bc1e
    TIMEOUT 30
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
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
set(easyloggingpp_INCLUDE_PATH "${easyloggingpp_SOURCE_DIR}")
list(APPEND 3rdparty_INCLUDES "${easyloggingpp_INCLUDE_PATH}")


if(pfasst_BUILD_EXAMPLES)
    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "FFTW3")

    find_package(FFTW)

    if(NOT FFTW_FOUND)
        set(FFTW3_SOURCE_URL "http://fftw.org/fftw-3.3.4.tar.gz")
        message(STATUS "FFTW3 not found on your system")
        message(STATUS "  going to download and compile it from sources")
        message(STATUS "    download URL: ${FFTW3_SOURCE_URL}")
        set(fftw3_SOURCE_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3")
        set(fftw3_INSTALL_DIR "${pfasst_BINARY_DIR}/3rdparty/src/fftw3-install")

        ExternalProject_Add(
            fftw3
            LIST_SEPARATOR " "
            URL ${FFTW3_SOURCE_URL}
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
endif()

if(pfasst_BUILD_TESTS)
    set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES})
    set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS})

    message(STATUS "--------------------------------------------------------------------------------")
    set(GMOCK_SOURCE_URL "http://googlemock.googlecode.com/files/gmock-1.7.0.zip")
    message(STATUS "Google Testing Framework (gtest & gmock)")
    message(STATUS "  going to download and compile it from sources")
    message(STATUS "    download URL: ${GMOCK_SOURCE_URL}")

    # Add gmock
    ExternalProject_Add(
      googlemock
      URL ${GMOCK_SOURCE_URL}
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
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS} PARENT_SCOPE)
set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES} PARENT_SCOPE)
