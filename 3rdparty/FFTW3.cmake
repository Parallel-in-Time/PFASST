if(pfasst_SYSTEM_FFTW3)
    message(STATUS "trying to use system's FFTW3")
    find_package(FFTW QUIET REQUIRED)

    message(STATUS "  include path: ${FFTW3_INCLUDE_PATH}")
    message(STATUS "  libraries: ${FFTW3_LIBRARIES}")
else()
    set(fftw3_SOURCE_NAME "fftw-3.3.4.tar.gz")
    set(fftw3_SOURCE_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")
    message(STATUS "  using version: 3.3.4  (file: ${fftw3_SOURCE_NAME})")

    if(NOT fftw3_SOURCES)
        set(fftw3_SOURCE_LOCATION "http://www.fftw.org")
        message(STATUS "  going to download from: ${fftw3_SOURCE_LOCATION}")
        message(STATUS "    HINT: specify -Dfftw3_SOURCES=/local/directory/with/${fftw3_SOURCE_NAME} to avoid redownloading")
    else()
        set(fftw3_SOURCE_LOCATION ${fftw3_SOURCES})
        message(STATUS "  using pre-downloaded sources at: ${fftw3_SOURCE_LOCATION}")
    endif()

    set(fftw3_PREFIX "${3rdparty_PREFIX}/fftw3")
    set(fftw3_TMP_PREFIX "${3rdparty_TMP_PREFIX}/fftw3")

    set(fftw3_TMP_DIR "${fftw3_TMP_PREFIX}/tmp")
    set(fftw3_STAMP_DIR "${fftw3_TMP_PREFIX}/stamp")
    set(fftw3_DOWNLOAD_DIR "${fftw3_PREFIX}/download")
    set(fftw3_SOURCE_DIR "${fftw3_PREFIX}/src")
    set(fftw3_BINARY_DIR "${fftw3_TMP_PREFIX}/build")
    set(fftw3_INSTALL_DIR "${fftw3_PREFIX}/install")

    ExternalProject_Add(fftw3
        LIST_SEPARATOR " "
        TMP_DIR "${fftw3_TMP_DIR}"
        STAMP_DIR "${fftw3_STAMP_DIR}"
        DOWNLOAD_DIR "${fftw3_DOWNLOAD_DIR}"
        SOURCE_DIR "${fftw3_SOURCE_DIR}"
        # BINARY_DIR "${fftw3_BINARY_DIR}"  # not used as in-source build is specified
        INSTALL_DIR "${fftw3_INSTALL_DIR}"

        EXCLUDE_FROM_ALL ON

        URL ${fftw3_SOURCE_LOCATION}/${fftw3_SOURCE_NAME}
        URL_MD5 ${fftw3_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ${fftw3_SOURCE_DIR}/configure
            CC=${CMAKE_C_COMPILER}
            CXX=${CMAKE_CXX_COMPILER}
            CXXFLAGS=${CMAKE_CXX_FLAGS}
            --prefix=${fftw3_INSTALL_DIR}
            --libdir=${fftw3_INSTALL_DIR}/lib
        BUILD_COMMAND make -j4
        TEST_COMMAND ""
        INSTALL_COMMAND make install

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )

    set(FFTW3_INCLUDE_PATH ${fftw3_INSTALL_DIR}/include)
    set(FFTW3_LIBRARY_DIR ${fftw3_INSTALL_DIR}/lib)
    add_library(fftw3_lib STATIC IMPORTED GLOBAL)
    set_property(TARGET fftw3_lib
        PROPERTY IMPORTED_LOCATION ${FFTW3_LIBRARY_DIR}/libfftw3.a
    )
    add_dependencies(fftw3_lib fftw3)
    set(FFTW3_LIBRARIES fftw3_lib)
    message(STATUS "  include path: ${FFTW3_INCLUDE_PATH}")
    message(STATUS "  library path: ${FFTW3_LIBRARY_DIR}")
endif()
