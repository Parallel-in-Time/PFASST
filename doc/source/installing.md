# Building and Installing                                                {#page_building_installing}

## Prerequesites

If you use [CMake] to compile, most of the prerequesites are downloaded automatically if they aren't
already found on your system.

### Required

* A recent C++ compiler that supports (most of) the C++11 standard.  This implies either: GCC 4.7
  (we recommend 4.8 or later) or Clang 3.2 (we recommend 3.3 or later).

* The [Boost] `boost::program_options` library (at least version 1.53.0).  The [Boost]
  multiprecision component may also be of interest to \em PFASST++ users.

* The [Eigen3] library.

### Optional

* The [CMake] build tool (at least version 2.8).  This is especially useful for developers.

* The [FFTW3] library is required for some of the examples.

* The [Google Testing] and [Google Mocking] are required for the unit tests.

* An MPI implementation.

### Using [HashDist] to obtain prerequesites

The [HashDist] environment management system provides an easy, robust, and version controlled way
for our developers (and users) to download and install all prerequesites *locally* (ie, you don't
need special privileges on your system to use [HashDist], and it won't modify your system).

To use [HashDist]

1. Install [HashDist]

        git clone git@github.com:hashdist/hashdist.git
        export PATH=$PWD/hashdist/bin:$PATH

2. Install prerequesites

        cd PFASST
        hit build -v tools/hashdist/pfasst-stack.debian.yaml
        mv pfasst-stack.debian stack

This may take a while (ie, do this before going for lunch) the first time through.  However, once
complete, subsequent builds (using [CMake] or otherwise) will bypass prerequesite installation, even
if you change [CMake] build options or start a fresh out-of-source build directory.


## Obtaining the Sources

Use _Git_ to clone the repository (https://github.com/Parallel-in-Time/PFASST.git) or download the
latest release from [GitHub][github_releases].


## Building with CMake

1. Create a directory for out-of-source builds (usually, people are creative and name it `build`) in
   the source tree (or elsewhere).

2. Decide on what and how you want to build:

   * __General__

     * By default everything gets built, i.e. tests, examples and all dependencies which are not
       found on the system.

     * Use `-DCMAKE_BUILD_TYPE=<VALUE>` to specify general compiler flags.

           | `CMAKE_BUILD_TYPE` | implied compiler flags |
           |--------------------|------------------------|
           | `Debug`            | `-g`                   |
           | `Release`          | `-O3 -DNDEBUG`         |
           | `RelWithDebInfo`   | `-O2 -g -DNDEBUG`      |
           | `MinSizeRel`       | `-Os -DNDEBUG`         |

     * To enable profiling support, you need to specify `-Dpfasst_WITH_GCC_PROF=ON` and have the GNU
       compiler selected.
       When compiling with _Clang_ this option is obsolete as profiling with _Clang_ is not supported.

     * Users on Linux systems with the _Clang_ compiler and a working installation of LLVM's
       \em libc++ library may want to activate the usage of that by specifying 
       `-Dpfasst_DISABLE_LIBCXX=OFF`.
       As libc++ is highly experimental on non-Darwin systems, this is a very exotic option.

   * __Dependencies__

     * __Boost__
       To specify the root path to the Boost installation to be used, add `-DBOOST_ROOT=<PATH>`.
       The root path should contain the directory `lib` or `lib64` with the compiled Boost libraries
       as well as `include/boost` with the Boost header files.

     * __GMock__
       In case there is the Google Testing and Mocking framework installed on the system in a
       non-standard path, add `-DGMOCK_ROOT=<PATH>` to specify the root path to the framework.

   * __MPI__

     * To enable MPI, please specify `-Dpfasst_WITH_MPI=ON`.
       To avoid a warning and potential undefined behaviour, also set `-DCMAKE_C_COMPILER` and
       `-DCMAKE_CXX_COMPILER` to the MPI compiler wrappers.

   * __Test Suite__

     * Deactivate building of the test suite by passing `-Dpfasst_BUILD_TESTS=OFF` to the _CMake_
       command line (not recommended).

   * __Examples__

     * Deactivate building of the example programms by passing `-Dpfasst_BUILD_EXAMPLES=OFF` to the
       _CMake_ command line.

     * Add compiled example programs to the `install` target by passing `-Dpfasst_INSTALL_EXAMPLES=ON`
       to the _CMake_ command line (by default, they will not get installed).

   * __Installing (optional)__

     * Use `-DCMAKE_INSTALL_PREFIX=<PREFIX>` to specify the prefix for installing the headers and
       optional compiled examples.

       By default, `<PREFIX>` is `/usr/local` on Linux systems.

       The headers will go into `<PREFIX>/include` and compiled binaries in `<PREFIX>/bin`.

     * In case you also want to install the compiled example programs on your system, specify
       `-Dpfasst_INSTALL_EXAMPLES=ON`.

   * For example, to build a release version with _Clang_ on a Linux system without the examples
     (i.e. only the core test suite and dependencies therefor), the call to _CMake_ looks like:

         cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang \
           -DCMAKE_CXX_COMPILER=clang++ -Dpfasst_BUILD_EXAMPLES=OFF ..

   * For a full debug build with enabled profiling on a Linux system with GCC and the desire to
     install everything to `/home/<USER>`, call _CMake_ with:

         cmake -DCMAKE_BUILD_TYPE=Debug -Dpfasst_WITH_GCC_PROF=ON \
           -Dpfasst_INSTALL_EXAMPLES=ON -DCMAKE_INSTALL_PREFIX=$HOME ..

3. Compile:

       make

4. Run the unit tests (if unit tests have been built):

       make test

   In case any of the tests do not pass, please [open an issue on GitHub][github_new_issue] and
   provide the full log of your _CMake_ and _make_ invocations.

5. (optional) Install \em PFASST++ headers (and compiled examples if specified so):

       make install

### Supercomputers

There are a few tweaks and special care require to get up and running on supercomputers.
We have compiled a few working walk-throughs at \subpage page_supercomputers .

### Score-P and other compiler wrappers

It is possible to use extra compiler wrappers like Score-P when compiling. See \subpage page_scorep.

## Building with vanilla make

A sample `Makefile` is included in the `advection_diffusion` example.  This may be of particular
interested to advanced users wishing to incorporate \em PFASST++ into their existing code bases.


[Boost]: https://boost.org/
[CMake]: http://cmake.org/
[Eigen3]: http://eigen.tuxfamily.org/
[Google Testing]: https://code.google.com/p/googletest/
[Google Mocking]: https://code.google.com/p/googlemock/
[FFTW3]: http://fftw.org/
[github_releases]: https://github.com/Parallel-in-Time/PFASST/releases
[github_new_issue]: https://github.com/Parallel-in-Time/PFASST/issues/new
[HashDist]: http://hashdist.readthedocs.org/en/latest/
