set(Boost_COMPONENTS "program_options")

if(pfasst_SYSTEM_BOOST)
    message(STATUS "trying to use system's Boost")
    set(Boost_USE_MULTITHREADED ON)

    find_package(Boost 1.53.0 QUIET COMPONENTS ${Boost_COMPONENTS} REQUIRED)

    message(STATUS "  found version: ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
    message(STATUS "  include path: ${Boost_INCLUDE_DIRS}")
    message(STATUS "  library path: ${Boost_LIBRARY_DIRS}")
else()
    set(Boost_SOURCE_NAME "boost_1_60_0.tar.bz2")
    set(Boost_SOURCE_MD5 "65a840e1a0b13a558ff19eeb2c4f0cbe")
    message(STATUS "  using version: 1.60.0  (file: ${Boost_SOURCE_NAME})")

    if(NOT Boost_SOURCES)
        set(Boost_SOURCE_LOCATION "https://downloads.sourceforge.net/project/boost/boost/1.60.0")
        message(STATUS "  going to download from: ${Boost_SOURCE_LOCATION}")
        message(STATUS "    HINT: specify -DBoost_SOURCES=/local/directory/with/${Boost_SOURCE_NAME} to avoid redownloading")
    else()
        set(Boost_SOURCE_LOCATION ${Boost_SOURCES})
        message(STATUS "  using pre-downloaded sources at: ${Boost_SOURCE_LOCATION}")
    endif()

    set(Boost_PREFIX "${3rdparty_PREFIX}/Boost")
    set(Boost_TMP_PREFIX "${3rdparty_TMP_PREFIX}/Boost")

    set(Boost_TMP_DIR "${Boost_TMP_PREFIX}/tmp")
    set(Boost_STAMP_DIR "${Boost_TMP_PREFIX}/stamp")
    set(Boost_DOWNLOAD_DIR "${Boost_PREFIX}/download")
    set(Boost_SOURCE_DIR "${Boost_PREFIX}/src")
    set(Boost_BINARY_DIR "${Boost_TMP_PREFIX}/build")
    set(Boost_INSTALL_DIR "${Boost_PREFIX}/install")

    set(Boost_build_CMD "${Boost_SOURCE_DIR}/bootstrap.sh")
    list(APPEND Boost_build_CMD
            --prefix=${Boost_INSTALL_DIR}
            --with-libraries=${Boost_COMPONENTS}
            )
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
        list(APPEND Boost_build_CMD toolset=clang)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND Boost_build_CMD toolset=gcc)
    endif()

    ExternalProject_Add(boost_build
        LIST_SEPARATOR " "
        TMP_DIR "${Boost_TMP_DIR}"
        STAMP_DIR "${Boost_STAMP_DIR}"
        DOWNLOAD_DIR "${Boost_DOWNLOAD_DIR}"
        SOURCE_DIR "${Boost_SOURCE_DIR}"
        # BINARY_DIR "${Boost_BINARY_DIR}"  # not used as in-source build is specified
        INSTALL_DIR "${Boost_INSTALL_DIR}"

        EXCLUDE_FROM_ALL ON

        URL ${Boost_SOURCE_LOCATION}/${Boost_SOURCE_NAME}
        URL_MD5 ${Boost_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${Boost_build_CMD}
        TEST_COMMAND ""
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
    )

    set(Boost_LIBRARIES)
    set(Boost_INCLUDE_DIRS ${Boost_SOURCE_DIR})
    set(Boost_LIBRARY_DIR ${Boost_SOURCE_DIR}/stage/lib)

    message(STATUS "  using components:")
    foreach(component ${Boost_COMPONENTS})
        message(STATUS "    ${component}")

        set(Boost_${component}_BUILD_CMD ${Boost_SOURCE_DIR}/b2)
        list(APPEND Boost_${component}_BUILD_CMD
            --build-dir=${Boost_BINARY_DIR}
            --with-${component}
            address-model=64
            cxxflags=-std=c++11
            cxxflags=-fPIC
            link=static
            threading=multi
            runtime-link=shared
            variant=release
        )
        if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
            list(APPEND Boost_${component}_BUILD_CMD toolset=clang)
        elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            list(APPEND Boost_${component}_BUILD_CMD toolset=gcc)
        endif()

        ExternalProject_Add(boost_${component}
            LIST_SEPARATOR " "
            TMP_DIR "${Boost_TMP_DIR}"
            STAMP_DIR "${Boost_STAMP_DIR}"
            DOWNLOAD_DIR "${Boost_DOWNLOAD_DIR}"
            SOURCE_DIR "${Boost_SOURCE_DIR}"
            # BINARY_DIR "${Boost_BINARY_DIR}"  # not used as in-source build is specified
            INSTALL_DIR "${Boost_INSTALL_DIR}"

            EXCLUDE_FROM_ALL ON
            DEPENDS boost_build

            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PATCH_COMMAND ""
            BUILD_IN_SOURCE ON
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ${Boost_${component}_BUILD_CMD} -j4 stage
            TEST_COMMAND ""
            INSTALL_COMMAND ""

            LOG_DOWNLOAD ON
            LOG_CONFIGURE ON
            LOG_BUILD ON
            LOG_INSTALL ON
        )

        list(APPEND pfasst_DEPENDEND_TARGETS boost_${component})
        list(APPEND Boost_LIBRARIES ${Boost_LIBRARY_DIR}/libboost_${component}.a)
    endforeach()

    message(STATUS "  include path: ${Boost_INCLUDE_DIRS}")
    message(STATUS "  library path: ${Boost_LIBRARY_DIR}")
endif()
