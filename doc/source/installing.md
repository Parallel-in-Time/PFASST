# Building and Installing                                                {#page_building_installing}

## Prerequesites

* A recent C++ compiler that supports (most of) the C++11 standard.
  This either implies at least GCC 4.7 (we recommend 4.8 or later) or Clang 3.2 (we recommend 3.3 or later).

* The [Boost](https://boost.org/) libraries at least of version 1.53.0.
  Especially `boost::program_options` is required.

* The [CMake] build tool (at least version 2.8).

* The [Eigen3] library is required for everything using the \em PFASST++ library.
  It will be automatically downloaded if it is not found on the system.

* The [FFTW3] library is required for (some) of the examples.
  It will be automatically downloaded and built as well if it is not found on the system.

* The [Google Testing] and [Mocking Framework] are required for the unit tests, and are automatically downloaded and 
  built.


## Obtaining the Sources

Either use _Git_ to clone the repository (https://github.com/Parallel-in-Time/PFASST.git) or download the latest 
release from [GitHub][github_releases].


## Building

1. Create a directory for out-of-source builds (usually, people are creative and name it `build`) in the source tree 
   (or elsewhere).

2. Decide on what and how you want to build:

   * __General__

     * By default everything gets built, i.e. tests and examples.

     * Use `-DCMAKE_BUILD_TYPE=<VALUE>` to specify general compiler flags.

           | `CMAKE_BUILD_TYPE` | implied compiler flags |
           |--------------------|------------------------|
           | `Debug`            | `-g`                   |
           | `Release`          | `-O3 -DNDEBUG`         |
           | `RelWithDebInfo`   | `-O2 -g -DNDEBUG`      |
           | `MinSizeRel`       | `-Os -DNDEBUG`         |

     * To enable profiling support, you need to specify `-Dpfasst_WITH_GCC_PROF=ON` and have the GNU compiler selected.
       When compiling with _Clang_ this option is obsolete as profiling with _Clang_ is not supported.

     * Users on Linux systems with the _Clang_ compiler and a working installation of LLVM's \em libc++ library may
       want to activate the usage of that by specifying `-Dpfasst_DISABLE_LIBCXX=OFF`.

   * __MPI__

     * To enable MPI, please specify `-Dpfasst_WITH_MPI=ON`.

   * __Test Suite__

     * Deactivate building of the test suite by passing `-Dpfasst_BUILD_TESTS=OFF` to the _CMake_ command line (not 
       recommended).

   * __Examples__

     * Deactivate building of the example programms by passing `-Dpfasst_BUILD_EXAMPLES=OFF` to the _CMake_ command
       line.

     * Add compiled example programs to the `install` target by passing `-Dpfasst_INSTALL_EXAMPLES=ON` to the _CMake_
       command line (by default, they will not get installed).

   * __Installing (optional)__

     * Use `-DCMAKE_INSTALL_PREFIX=<PREFIX>` to specify the prefix for installing the headers and optional compiled
       examples.

       By default, `<PREFIX>` is `/usr/local` on Linux systems.

       The headers will go into `<PREFIX>/include` and compiled binaries in `<PREFIX>/bin`.

     * In case you also want to install the compiled example programs on your system, specify
       `-Dpfasst_INSTALL_EXAMPLES=ON`.

   * For example, to build a release version with _Clang_ on a Linux system without the examples, the call to 
     _CMake_ looks like:

         cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -Dpfasst_BUILD_EXAMPLES=OFF ..

   * A full debug build on a Linux system with GCC and the desire to install everything to `/home/<USER>`, call _CMake_
     with:

         cmake -DCMAKE_BUILD_TYPE=Debug -Dpfasst_WITH_GCC_PROF=ON -Dpfasst_INSTALL_EXAMPLES=ON -DCMAKE_INSTALL_PREFIX=$HOME ..

3. Compile:

       make

4. Run the unit tests (if unit tests have been built):

       make test

   In case any of the tests do not pass, please [open an issue on GitHub][github_new_issue] and provide the full log of 
   your _CMake_ and _make_ invocations.

5. (optional) Install \em PFASST++ headers (and compiled examples if specified so):

       make install


[CMake]: http://cmake.org/
[Eigen3]: http://eigen.tuxfamily.org/
[Google Testing]: https://code.google.com/p/googletest/
[Mocking Framework]: https://code.google.com/p/googlemock/
[FFTW3]: http://fftw.org/
[github_releases]: https://github.com/Parallel-in-Time/PFASST/releases
[github_new_issue]: https://github.com/Parallel-in-Time/PFASST/issues/new
