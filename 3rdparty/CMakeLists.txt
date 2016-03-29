# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS})
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS})

if(3rdparty_SOURCES)
    message(STATUS "trying to pick up pre-downloaded sources from: ${3rdparty_SOURCES}")
    if(NOT Boost_SOURCES)
        if(EXISTS ${3rdparty_SOURCES}/boost)
            set(Boost_SOURCES ${3rdparty_SOURCES}/boost)
        else()
            message(STATUS "!! for Boost: <3rdparty_SOURCES>/boost does not exists")
            message(STATUS "   use -DBoost_SOURCES=/local/path to overwrite <3rdparty_SOURCES>")
        endif()
    endif()
    if(NOT Eigen3_SOURCES)
        if(EXISTS ${3rdparty_SOURCES}/eigen3)
            set(Eigen3_SOURCES ${3rdparty_SOURCES}/eigen3)
        else()
            message(STATUS "!! for Eigen3: <3rdparty_SOURCES>/eigen3 does not exists")
            message(STATUS "   use -DEigen3_SOURCES=/local/path to overwrite <3rdparty_SOURCES>")
        endif()
    endif()
    if(NOT fftw3_SOURCES)
        if(EXISTS ${3rdparty_SOURCES}/fftw3)
            set(fftw3_SOURCES ${3rdparty_SOURCES}/fftw3)
        else()
            message(STATUS "!! for FFTW3: <3rdparty_SOURCES>/fftw3 does not exists")
            message(STATUS "   use -Dfftw3_SOURCES=/local/path to overwrite <3rdparty_SOURCES>")
        endif()
    endif()
    if(NOT gtest_SOURCES)
        if(EXISTS ${3rdparty_SOURCES}/gtest)
            set(gtest_SOURCES ${3rdparty_SOURCES}/gtest)
        else()
            message(STATUS "!! for GTest: <3rdparty_SOURCES>/gtest does not exists")
            message(STATUS "   use -Dgtest_SOURCES=/local/path to overwrite <3rdparty_SOURCES>")
        endif()
    endif()
    if(NOT gmock_SOURCES)
        if(EXISTS ${3rdparty_SOURCES}/gmock)
            set(gmock_SOURCES ${3rdparty_SOURCES}/gmock)
        else()
            message(STATUS "!! for GMock: <3rdparty_SOURCES>/gmock does not exists")
            message(STATUS "   use -Dgmock_SOURCES=/local/path to overwrite <3rdparty_SOURCES>")
        endif()
    endif()
endif()

set(3rdparty_DOWNLOAD_TIMEOUT "60")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Boost")

include(Boost.cmake)

list(APPEND 3rdparty_INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND 3rdparty_DEPENDEND_LIBS ${Boost_LIBRARIES})


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Eigen3")

include(Eigen3.cmake)

list(APPEND 3rdparty_INCLUDES "${Eigen3_INCLUDE_PATH}")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "FFTW3")

include(FFTW3.cmake)

set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} PARENT_SCOPE)
set(FFTW3_INCLUDE_PATH ${FFTW3_INCLUDE_PATH} PARENT_SCOPE)


if(pfasst_BUILD_TESTS)
    set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES})
    set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS})

    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "Google Testing Framework")

    include(googletest.cmake)

    list(APPEND TESTS_3rdparty_INCLUDES ${GOOGLETEST_INCLUDE_DIRS})
    list(APPEND TESTS_3rdparty_DEPENDEND_LIBS gtest_lib gmock_lib ${CMAKE_THREAD_LIBS_INIT})
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS} PARENT_SCOPE)
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS} PARENT_SCOPE)
set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES} PARENT_SCOPE)
