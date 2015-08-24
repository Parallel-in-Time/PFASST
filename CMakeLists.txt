cmake_minimum_required(VERSION 2.8.6)
project(pfasst)

list(APPEND CMAKE_MODULE_PATH ${pfasst_SOURCE_DIR}/cmake)
include(cmake/utility_functions.cmake)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(CMakeDependentOption)
include(FeatureSummary)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rdparty)

list(APPEND CMAKE_FIND_ROOT_PATH ${pfasst_SOURCE_DIR}/stack)

option(pfasst_DISABLE_LIBCXX    "Disable use of LLVM's libc++ when compiling with Clang." ON )
option(pfasst_BUILD_SHARED_LIBS "Build shared libraries."                                 ON )
option(pfasst_BUILD_EXAMPLES    "Build example programs."                                 ON )
cmake_dependent_option(
       pfasst_INSTALL_EXAMPLES  "Install example programs."                               ON
           "pfasst_BUILD_EXAMPLES" OFF)
option(pfasst_BUILD_TESTS       "Build test suite for PFASST."                            ON )
option(pfasst_WITH_MPI          "Build with MPI enabled."                                 ON )
cmake_dependent_option(
       pfasst_WITH_MPIP         "enable to link against MPIP"                             OFF
           "pfasst_WITH_MPI" ON)
option(pfasst_WITH_GCC_PROF     "Enable excessive debugging & profiling output with GCC." OFF)
cmake_dependent_option(
       enable_LTO               "enable LinkTimeOptimization"                             OFF
           "CMAKE_BUILD_TYPE" Release)
cmake_dependent_option(
       pfasst_WITH_SCOREP       "enable scalable instrumentation with Score-P"            OFF
           "pfasst_WITH_MPI" ON)

set(3rdparty_INCLUDES)
set(3rdparty_DEPENDEND_LIBS)
set(pfasst_INCLUDES)
set(pfasst_DEPENDEND_LIBS)
set(pfasst_DEPENDEND_TARGETS)
set(pfasst_TESTS_DEPENDEND_TARGETS)
option(pfasst_WITH_EXTRA_WRAPPER "Build script to allow to set an additional wrapper via $PREP" OFF)

option(pfasst_WITH_LOGGING      "Enable logging output on console"                        ON)

message(STATUS "********************************************************************************")
message(STATUS "Further tests on the compiler")
if(${pfasst_BUILD_SHARED_LIBS})
    message(STATUS "build dynamic libraries and shared linked executables")
else()
    message(STATUS "build static libraries and staticaly linked executables")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()
message(STATUS "")

if(${pfasst_WITH_MPI})
    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "Detecting MPI")
    find_package(MPI REQUIRED)
    if(NOT "${CMAKE_CXX_COMPILER}" STREQUAL "${MPI_CXX_COMPILER}")
        message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
        message(STATUS "MPI C++ Compiler Wrapper: ${MPI_CXX_COMPILER}")
        message(WARNING "Please make sure to set CXX to the MPI compiler wrappers!")
        set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        message(STATUS "Set default compilers to MPI compilers.")
    endif()
    message(STATUS "Using MPI C++ Compiler Wrapper: ${MPI_CXX_COMPILER}")
    add_definitions(-DWITH_MPI)
endif()

if(pfasst_WITH_EXTRA_WRAPPER)
    configure_file(
        "${pfasst_SOURCE_DIR}/cmake/cxx_wrapper.sh.in"
        "${CMAKE_BINARY_DIR}/cxx_wrapper.sh"
    )
    configure_file(
        "${pfasst_SOURCE_DIR}/cmake/cc_wrapper.sh.in"
        "${CMAKE_BINARY_DIR}/cc_wrapper.sh"
    )
    set(CMAKE_CXX_COMPILER "${CMAKE_BINARY_DIR}/cxx_wrapper.sh")
    set(CMAKE_C_COMPILER "${CMAKE_BINARY_DIR}/cc_wrapper.sh")
endif()

# Check for C++11 support
message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Checking C++11 standard support")
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
    check_cxx_compiler_flag(-std=c++11 HAVE_STD11)
    if(HAVE_STD11)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-std=c++11")
    else()
        message(FATAL_ERROR "No advanced standard C++ support of your GCC (-std=c++11 not defined).")
    endif()
    if(${CMAKE_BUILD_TYPE} MATCHES Debug)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Og")
    endif()
    if(pfasst_WITH_GCC_PROF)
        message(STATUS "enabling profiling support")
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-ggdb3 -pg")
    endif(pfasst_WITH_GCC_PROF)
    if(${CMAKE_BUILD_TYPE} MATCHES "Release" AND ${enable_LTO})
        message(STATUS "enabling Link Time Optimization (LTO)")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto -fwhole-program")
        set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -flto")
        set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} -flto")
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    check_cxx_compiler_flag(-std=c++11 HAVE_STD11)
    if(HAVE_STD11)
        if(pfasst_DISABLE_LIBCXX AND NOT APPLE)
            add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-std=c++11")
            add_to_string_list("${CMAKE_CXX_LINK_FLAGS}" CMAKE_CXX_LINK_FLAGS "-std=c++11")
        else()
            if(APPLE)
                message(WARNING "You are on an Apple system. libc++ is forcly enabled.")
            endif(APPLE)
            add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-std=c++11" "-stdlib=libc++")
            add_to_string_list("${CMAKE_CXX_LINK_FLAGS}" CMAKE_CXX_LINK_FLAGS "-std=c++11" "-stdlib=libc++")
        endif()
    else()
        message(FATAL_ERROR "No C++11 support for Clang version. Please upgrade Clang to a version supporting C++11.")
    endif()
    add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-fdiagnostics-show-category=name -fdiagnostics-show-template-tree")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES XL)
    # NOTE: This branch is not tested yet (in theory it should work)
    message(WARNING "IMB XL C/C++ support is experimental and not yet tested.")
    check_cxx_compiler_flag(-qlanglvl=extended0x HAVE_STD11)
    if(HAVE_STD11)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-qlanglvl=extended0x -qwarn0x")
    else()
        message(FATAL_ERROR "No advanced standard C++ support of your IBM XL C/C++ compiler (-qlanglvl=extended0x not defined).")
    endif()
else()
    message(FATAL_ERROR "Don't know how to check C++11 compatibility with compiler '${CMAKE_CXX_COMPILER_ID}'")
endif()
message(STATUS "Your compiler has C++11 support. Hurray!")
message(STATUS "")

if(CMAKE_CXX_COMPILER_VERSION)
    string(REGEX MATCH "^[0-9]+" CMAKE_CXX_COMPILER_VERSION_MAJOR ${CMAKE_CXX_COMPILER_VERSION})
    string(REGEX MATCH "\\.[0-9]+\\." cmake_cxx_compiler_version_minor ${CMAKE_CXX_COMPILER_VERSION})
    string(REGEX MATCH "[0-9]+" CMAKE_CXX_COMPILER_VERSION_MINOR ${cmake_cxx_compiler_version_minor})
    string(REGEX MATCH "[0-9]+$" CMAKE_CXX_COMPILER_VERSION_PATCH ${CMAKE_CXX_COMPILER_VERSION})
    set(compiler_version_available TRUE)
else()
    message(WARNING "Compiler Version not available. Configuration might fail.")
    set(compiler_version_available FALSE)
endif()

# Enable all compiler warnings
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU OR ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
endif()


# Set features from options
add_feature_info(Examples
    pfasst_BUILD_EXAMPLES
    "build example programs"
)
add_feature_info(Tests
    pfasst_BUILD_TESTS
    "build test suite"
)
add_feature_info(MPI
    pfasst_WITH_MPI
    "build with MPI"
)
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
    add_feature_info(Profiling
        pfasst_WITH_GCC_PROF
        "build-in profiling support"
    )
endif()

if(pfasst_BUILD_TESTS)
    enable_testing()
    set(TESTS_3rdparty_INCLUDES)
    set(TESTS_3rdparty_DEPENDEND_LIBS)
endif(pfasst_BUILD_TESTS)

# Add / include 3rd-party libraries
message(STATUS "********************************************************************************")
message(STATUS "Configuring 3rd party libraries")
# makes available:
#  - Boost headers in 3rdparty_INCLUDES
#  - Boost libraries in Boost_LIBRARY_DIRS
#  - Eigen headers in 3rdparty_INCLUDES
#  - Google test and mock headers in 3rdparty_INCLUDES (if pfasst_BUILD_TESTS)
#  - FFTW_INCLUDE_PATH (if pfasst_BUILD_EXAMPLES)
#  - FFTW_LIBRARIES (if pfasst_BUILD_EXAMPLES)
add_subdirectory(3rdparty)
message(STATUS "")

message(STATUS "********************************************************************************")
message(STATUS "Configuring sources")

if(NOT PYTHON_EXECUTABLE)
    find_package(PythonInterp 2.7)
    if(NOT PYTHONINTERP_FOUND)
        find_package(PythonInterp 3 REQUIRED)
    endif()
else()
    message(STATUS "Using specified Python executable: ${PYTHON_EXECUTABLE}")
endif()
add_custom_target(set_version ALL
    COMMAND ${PYTHON_EXECUTABLE} ${pfasst_SOURCE_DIR}/tools/get_pfasst_version.py
    WORKING_DIRECTORY ${pfasst_SOURCE_DIR}
    COMMENT "Updating PFASST++ version number"
    USES_TERMINAL
)
set(pfasst_DEPENDEND_TARGETS set_version ${pfasst_DEPENDEND_TARGETS})

list(LENGTH pfasst_DEPENDEND_LIBS pfasst_NUM_DEPENDEND_LIBS)
list(LENGTH pfasst_DEPENDEND_TARGETS pfasst_NUM_DEPENDEND_TARGETS)
list(LENGTH pfasst_TESTS_DEPENDEND_TARGETS pfasst_TESTS_NUM_DEPENDEND_TARGETS)

add_subdirectory(include)
add_subdirectory(src)
message(STATUS "")

if(pfasst_BUILD_EXAMPLES)
    message(STATUS "********************************************************************************")
    message(STATUS "Configuring examples")
    set(examples_to_install)
    add_subdirectory(examples)
    message(STATUS "")
endif()

if(pfasst_BUILD_TESTS)
    message(STATUS "********************************************************************************")
    message(STATUS "Configuring tests")
    add_subdirectory(tests)
    message(STATUS "")
endif()
message(STATUS "********************************************************************************")

message(STATUS "Creating 'install' target")
message(STATUS "  prefix:              ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  header files:        <PREFIX>/include/pfasst")
install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)
install(DIRECTORY src/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

if(pfasst_BUILD_EXAMPLES AND pfasst_INSTALL_EXAMPLES)
    message(STATUS "  example programs:    <PREFIX>/bin")
    message(STATUS "  Example programms have been added to 'install' target")
    message(STATUS "    example programs:")
    foreach(example_program ${examples_to_install})
        message(STATUS "      - ${example_program}")
    endforeach()
endif()
message(STATUS "")

message(STATUS "********************************************************************************")
if(${CMAKE_VERBOSE_MAKEFILE})
  message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
  message(STATUS "C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
  message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS}")
  message(STATUS "C++ link flags: ${CMAKE_CXX_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
  message(STATUS "Dependend Libs: ${pfasst_DEPENDEND_LIBS}")
  message(STATUS "********************************************************************************")
endif()

feature_summary(WHAT ALL)
