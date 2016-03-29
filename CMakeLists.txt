cmake_minimum_required(VERSION 3.0)
project(pfasst)

list(APPEND CMAKE_MODULE_PATH ${pfasst_SOURCE_DIR}/cmake)
include(cmake/utility_functions.cmake)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(CMakeDependentOption)
include(FeatureSummary)

list(APPEND CMAKE_FIND_ROOT_PATH ${pfasst_SOURCE_DIR}/stack)

set(3rdparty_PREFIX "${CMAKE_SOURCE_DIR}/3rdparty")
set(3rdparty_TMP_PREFIX "${CMAKE_BINARY_DIR}/3rdparty")
# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX "${3rdparty_PREFIX}")


option(pfasst_DISABLE_LIBCXX    "Disable use of LLVM's libc++ when compiling with Clang." ON )
option(BUILD_SHARED_LIBS        "Build shared libraries."                                 OFF)

option(pfasst_BUILD_EXAMPLES    "Build example programs."                                 ON )
cmake_dependent_option(
       pfasst_INSTALL_EXAMPLES  "Install example programs."                               ON
            "pfasst_BUILD_EXAMPLES" OFF)

option(pfasst_BUILD_TESTS       "Build test suite for PFASST."                            ON )

option(pfasst_WITH_MPI          "Build with MPI enabled."                                 ON )
option(pfasst_SYSTEM_BOOST      "Use system's Boost library; fail if not found"           OFF)
option(pfasst_SYSTEM_FFTW3      "Use system's FFTW library; fail if not found"            OFF)
option(pfasst_SYSTEM_EIGEN3     "Use system's Eigen3 library; fail if not found"          OFF)

option(pfasst_WITH_GCC_PROF     "Enable excessive debugging & profiling output with GCC." OFF)
cmake_dependent_option(
       enable_LTO               "enable LinkTimeOptimization"                             OFF
            "CMAKE_BUILD_TYPE" Release)

option(pfasst_WITH_EXTRA_WRAPPER "Build script to allow to set an additional wrapper via $PREP" OFF)

option(pfasst_WITH_LOGGING      "Enable logging output on console"                        ON)

if(NOT ${pfasst_WITH_LOGGING})
    add_definitions(-DPFASST_NO_LOGGING)
endif()

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
    add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Wdeprecated -fdiagnostics-show-category=name -fdiagnostics-show-template-tree")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
    check_cxx_compiler_flag(-std=c++11 HAVE_STD11)
    if(HAVE_STD11)
        add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-std=c++11")
    else()
        message(FATAL_ERROR "No advanced standard C++ support for your Intel compiler (-std=c++11 not defined).")
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES XL)
    message(FATAL_ERROR "IMB XL C/C++ does not support C++11 features PFASST++ is utilizing.")
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
if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
    add_to_string_list("${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS "-Wall -Wextra")
else()
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
#  - required headers in 3rdparty_INCLUDES
#  - required libraries in 3rdparty_DEPENDEND_LIBS and pfasst_DEPENDEND_LIBS
#  - Google test and mock headers in TESTS_3rdparty_INCLUDES (if pfasst_BUILD_TESTS)
#  - Google test and mock libraries in TESTS_3rdparty_DEPENDEND_LIBS (if pfasst_BUILD_TESTS)
#  - FFTW3_INCLUDE_PATH
#  - FFTW3_LIBRARIES
add_subdirectory(3rdparty)

message(STATUS "********************************************************************************")
message(STATUS "Configuring sources")
update_site_config()

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
