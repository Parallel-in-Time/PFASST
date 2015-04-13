cmake_minimum_required(VERSION 2.8.6)
project(pfasst)

list(APPEND CMAKE_MODULE_PATH ${pfasst_SOURCE_DIR}/cmake)
include(cmake/utility_functions.cmake)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(CMakeDependentOption)
include(FeatureSummary)

list(APPEND CMAKE_FIND_ROOT_PATH ${pfasst_SOURCE_DIR}/stack)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rdparty)

option(pfasst_DISABLE_LIBCXX    "Disable use of LLVM's libc++ when compiling with Clang." ON )
option(pfasst_BUILD_SHARED_LIBS "Build shared libraries."                                 ON )
option(pfasst_BUILD_EXAMPLES    "Build example programs."                                 ON )
CMAKE_DEPENDENT_OPTION(pfasst_INSTALL_EXAMPLES "Install example programs." ON
    "pfasst_BUILD_EXAMPLES" OFF)
option(pfasst_BUILD_TESTS       "Build test suite for PFASST."                            ON )
option(pfasst_WITH_MPI          "Build with MPI enabled."                                 OFF)
option(pfasst_WITH_GCC_PROF     "Enable excessive debugging & profiling output with GCC." OFF)
option(pfasst_DEFAULT_RAND_SEED "Using a hardcoded random seed"                           ON )

if(${pfasst_WITH_MPI})
    find_package(MPI REQUIRED)
    if(NOT "${CMAKE_CXX_COMPILER}" STREQUAL "${MPI_CXX_COMPILER}" OR 
       NOT "${CMAKE_C_COMPILER}" STREQUAL "${MPI_C_COMPILER}")
        message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
        message(STATUS "MPI C++ Compiler Wrapper: ${MPI_CXX_COMPILER}")
        message(WARNING "Please make sure to set CXX and CC to the MPI compiler wrappers!")
        set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        set(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
        message(STATUS "Set default compilers to MPI compilers.")
    endif()
    message(STATUS "Using MPI C++ Compiler Wrapper: ${MPI_CXX_COMPILER}")
    message(STATUS "USing MPI C Compiler Wrapper: ${MPI_C_COMPILER}")
    add_definitions(-DWITH_MPI)
endif()

if(pfasst_DEFAULT_RAND_SEED)
    set(pfasst_RANDOM_SEED "42")
endif()

# Check for C++11 support
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
    check_cxx_compiler_flag(-std=c++11 HAVE_STD11)
    if(HAVE_STD11)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-std=c++11")
    else()
        message(FATAL_ERROR "No advanced standard C++ support of your GCC (-std=c++11 not defined).")
    endif()
    if(pfasst_WITH_GCC_PROF)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-ggdb3 -pg")
    endif(pfasst_WITH_GCC_PROF)
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
add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")


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
if(pfasst_RANDOM_SEED)
    add_feature_info(DefaultRandomSeed
        pfasst_RANDOM_SEED
        "using specific default random seed ${pfasst_RANDOM_SEED}"
    )
endif()

set(3rdparty_INCLUDES)
set(3rdparty_DEPENDEND_LIBS)
set(pfasst_INCLUDES)
set(pfasst_DEPENDEND_LIBS)
set(pfasst_DEPENDEND_TARGETS)
set(pfasst_TESTS_DEPENDEND_TARGETS)

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

list(LENGTH pfasst_DEPENDEND_LIBS pfasst_NUM_DEPENDEND_LIBS)
list(LENGTH pfasst_DEPENDEND_TARGETS pfasst_NUM_DEPENDEND_TARGETS)
list(LENGTH pfasst_TESTS_DEPENDEND_TARGETS pfasst_TESTS_NUM_DEPENDEND_TARGETS)

message(STATUS "********************************************************************************")
message(STATUS "Configuring sources")
set(WARNING_COMMENT "/*\n * DO NOT ALTER THIS FILE\n *\n * It will get rewritten on CMake's next run\n *\n */")
if(pfasst_RANDOM_SEED)
    add_definitions(-DPFASST_DEFAULT_RANDOM_SEED)
endif()
configure_file(
    "${pfasst_SOURCE_DIR}/cmake/site_config.hpp.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/pfasst/site_config.hpp"
)
add_subdirectory(include)
add_subdirectory(src)

if(pfasst_BUILD_EXAMPLES)
    message(STATUS "********************************************************************************")
    message(STATUS "Configuring examples")
    set(examples_to_install)
    add_subdirectory(examples)
endif()

if(pfasst_BUILD_TESTS)
    message(STATUS "********************************************************************************")
    message(STATUS "Configuring tests")
    add_subdirectory(tests)
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

message(STATUS "********************************************************************************")
if(${CMAKE_VERBOSE_MAKEFILE})
  message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
  message(STATUS "C++ Compiler Names: ${CMAKE_CXX_COMPILER_NAMES}")
  message(STATUS "C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
  message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS}")
  message(STATUS "C++ link flags: ${CMAKE_CXX_LINK_FLAGS}")
  message(STATUS "********************************************************************************")
endif()

feature_summary(WHAT ALL)
