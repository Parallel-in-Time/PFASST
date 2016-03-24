if(pfasst_SYSTEM_BOOST)
    message(STATUS "trying to use system's Eigen3")

    find_package(Eigen3 QUIET REQUIRED)

    message(STATUS "  include path: ${Eigen3_INCLUDE_PATH}")
else()
    set(Eigen3_SOURCE_NAME "3.2.5.tar.bz2")
    set(Eigen3_SOURCE_MD5 "21a928f6e0f1c7f24b6f63ff823593f5")
    message(STATUS "  using version: 3.2.5  (file: ${Eigen3_SOURCE_NAME})")

    if(NOT Eigen3_SOURCES)
        set(Eigen3_SOURCE_LOCATION "https://bitbucket.org/eigen/eigen/get")
        message(STATUS "  going to download from: ${Eigen3_SOURCE_LOCATION}")
        message(STATUS "    HINT: specify -DEigen3_SOURCES=/local/directory/with/${Eigen3_SOURCE_NAME} to avoid redownloading")
    else()
        set(Eigen3_SOURCE_LOCATION ${Eigen3_SOURCES})
        message(STATUS "  using pre-downloaded sources at: ${Eigen3_SOURCE_LOCATION}")
    endif()

    set(Eigen3_PREFIX "${3rdparty_PREFIX}/Eigen3")
    set(Eigen3_TMP_PREFIX "${3rdparty_TMP_PREFIX}/Eigen3")

    set(Eigen3_TMP_DIR "${Eigen3_TMP_PREFIX}/tmp")
    set(Eigen3_STAMP_DIR "${Eigen3_TMP_PREFIX}/stamp")
    set(Eigen3_DOWNLOAD_DIR "${Eigen3_PREFIX}/download")
    set(Eigen3_SOURCE_DIR "${Eigen3_PREFIX}/src")
    set(Eigen3_BINARY_DIR "${Eigen3_TMP_PREFIX}/build")
    set(Eigen3_INSTALL_DIR "${Eigen3_PREFIX}/install")

    ExternalProject_Add(Eigen3
        LIST_SEPARATOR " "
        TMP_DIR "${Eigen3_TMP_DIR}"
        STAMP_DIR "${Eigen3_STAMP_DIR}"
        DOWNLOAD_DIR "${Eigen3_DOWNLOAD_DIR}"
        SOURCE_DIR "${Eigen3_SOURCE_DIR}"
        # BINARY_DIR "${Eigen3_BINARY_DIR}"  # not used as in-source build is specified
        INSTALL_DIR "${Eigen3_INSTALL_DIR}"

        EXCLUDE_FROM_ALL ON

        URL ${Eigen3_SOURCE_LOCATION}/${Eigen3_SOURCE_NAME}
        URL_MD5 ${Eigen3_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        TEST_COMMAND ""
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )

    # Specify include dir
    set(Eigen3_INCLUDE_PATH "${Eigen3_SOURCE_DIR}")
    message(STATUS "  include path: ${Eigen3_INCLUDE_PATH}")

    list(APPEND pfasst_DEPENDEND_TARGETS Eigen3)
endif()

set(Eigen3_INCLUDE_PATH ${Eigen3_INCLUDE_PATH} PARENT_SCOPE)
