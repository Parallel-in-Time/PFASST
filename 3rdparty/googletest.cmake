set(gtest_SOURCE_NAME "release-1.7.0.tar.gz")
message(STATUS "  using version: 1.7.0  (file: ${gtest_SOURCE_NAME})")

set(gtest_SOURCE_MD5 "4ff6353b2560df0afecfbda3b2763847")

if(NOT gtest_SOURCES)
    set(gtest_SOURCE_LOCATION "https://github.com/google/googletest/archive")
    message(STATUS "  going to download from: ${gtest_SOURCE_LOCATION}")
    message(STATUS "    HINT: specify -Dgtest_SOURCES=/local/directory/with/${gtest_SOURCE_NAME} to avoid redownloading")
else()
    set(gtest_SOURCE_LOCATION ${gtest_SOURCES})
    message(STATUS "  using pre-downloaded sources at: ${gtest_SOURCE_LOCATION}")
endif()


set(gmock_SOURCE_NAME "release-1.7.0.tar.gz")
set(gmock_SOURCE_MD5 "13c3b4a57ad575763deb73fc0ad96e07")

if(NOT gmock_SOURCES)
    set(gmock_SOURCE_LOCATION "https://github.com/google/googlemock/archive")
    message(STATUS "  going to download from: ${gmock_SOURCE_LOCATION}")
    message(STATUS "    HINT: specify -Dgmock_SOURCES=/local/directory/with/${gmock_SOURCE_NAME} to avoid redownloading")
else()
    set(gmock_SOURCE_LOCATION ${gmock_SOURCES})
    message(STATUS "  using pre-downloaded sources at: ${gmock_SOURCE_LOCATION}")
endif()

set(gmock_PREFIX "${3rdparty_PREFIX}/gmock")
set(gmock_TMP_PREFIX "${3rdparty_TMP_PREFIX}/gmock")

set(gmock_TMP_DIR "${gmock_TMP_PREFIX}/tmp")
set(gmock_STAMP_DIR "${gmock_TMP_PREFIX}/stamp")
set(gmock_DOWNLOAD_DIR "${gmock_PREFIX}/download")
set(gmock_SOURCE_DIR "${gmock_PREFIX}/src")
set(gmock_BINARY_DIR "${gmock_TMP_PREFIX}/build")

set(gtest_PREFIX "${3rdparty_PREFIX}/gtest")
set(gtest_TMP_PREFIX "${3rdparty_TMP_PREFIX}/gtest")

set(gtest_TMP_DIR "${gtest_TMP_PREFIX}/tmp")
set(gtest_STAMP_DIR "${gtest_TMP_PREFIX}/stamp")
set(gtest_DOWNLOAD_DIR "${gtest_PREFIX}/download")
# gtest sources should be in a subdir of gmock sources
set(gtest_SOURCE_DIR "${gmock_PREFIX}/gtest")
# and gtest is build by gmock
set(gtest_BINARY_DIR "${gmock_BINARY_DIR}/gtest")

# Add gtest -- only download
ExternalProject_Add(gtest
    TMP_DIR "${gtest_TMP_DIR}"
    STAMP_DIR "${gtest_STAMP_DIR}"
    DOWNLOAD_DIR "${gtest_DOWNLOAD_DIR}"
    SOURCE_DIR "${gtest_SOURCE_DIR}"
    BINARY_DIR "${gtest_BINARY_DIR}"
    INSTALL_DIR ""

    EXCLUDE_FROM_ALL ON

    URL ${gtest_SOURCE_LOCATION}/${gtest_SOURCE_NAME}
    URL_MD5 ${gtest_SOURCE_MD5}
    TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""

    LOG_DOWNLOAD ON
)

# Add gmock
ExternalProject_Add(gmock
    TMP_DIR "${gmock_TMP_DIR}"
    STAMP_DIR "${gmock_STAMP_DIR}"
    DOWNLOAD_DIR "${gmock_DOWNLOAD_DIR}"
    SOURCE_DIR "${gmock_SOURCE_DIR}"
    BINARY_DIR "${gmock_BINARY_DIR}"
    INSTALL_DIR ""

    DEPENDS gtest
    EXCLUDE_FROM_ALL ON

    URL ${gmock_SOURCE_LOCATION}/${gmock_SOURCE_NAME}
    URL_MD5 ${gmock_SOURCE_MD5}
    TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILE=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -Dgtest_USE_OWN_TR1_TUPLE=ON
        -Dgtest_force_shared_crt=ON
        -Dgtest_build_tests=OFF
    BUILD_COMMAND make -j4
    INSTALL_COMMAND ""

    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

set(gtest_INCLUDE_DIR ${gtest_SOURCE_DIR}/include)
set(gtest_LIBRARY_DIR ${gtest_BINARY_DIR})
set(gmock_INCLUDE_DIR ${gmock_SOURCE_DIR}/include)
set(gmock_LIBRARY_DIR ${gmock_BINARY_DIR})

add_library(gtest_lib STATIC IMPORTED GLOBAL)
set_property(TARGET gtest_lib
    PROPERTY IMPORTED_LOCATION ${gtest_LIBRARY_DIR}/libgtest.a
)
add_dependencies(gtest_lib gmock)

add_library(gmock_lib STATIC IMPORTED GLOBAL)
set_property(TARGET gmock_lib
    PROPERTY IMPORTED_LOCATION ${gmock_LIBRARY_DIR}/libgmock.a
)
add_dependencies(gmock_lib gmock)

message(STATUS "  include path: ${gtest_INCLUDE_DIR}")
message(STATUS "  include path: ${gmock_INCLUDE_DIR}")
message(STATUS "  library path: ${gtest_LIBRARY_DIR}")
message(STATUS "  library path: ${gmock_LIBRARY_DIR}")

set(GOOGLETEST_INCLUDE_DIRS ${gtest_INCLUDE_DIR} ${gmock_INCLUDE_DIR})
set(GOOGLETEST_LIBRARY_DIRS ${gtest_LIBRARY_DIR} ${gmock_LIBRARY_DIR})

set(CMAKE_THREAD_PREFER_PTHREAD)
find_package(Threads QUIET)
