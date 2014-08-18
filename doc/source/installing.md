# Building and Installing                                                {#page_building_installing}

## Prerequesites

* A recent C++ compiler that supports (most of) the C++11 standard.
* The [CMake](http://cmake.org/) build tool.


## Obtaining the Sources

Just use _Git_ to clone the repository.


## Building

1. Create a directory for out-of-source builds (usually, people are
   creative and name it `build`) in the source tree (or elsewhere).

2. Decide on what you want to build:

   * By default everything gets built, i.e. tests and examples.

   * To deactivate the tests, pass `-Dpfasst_BUILD_TESTS=OFF` to the _CMake_ command line.

   * To deactivate the examples, pass `-Dpfasst_BUILD_EXAMPLES=OFF` to the _CMake_ command line.

   * To enable profiling support, you need to specify `-Dpfasst_WITH_GCC_PROF` and have the GNU 
     compiler selected.

   * If you are not on a _Darvin_ system (i.e. MacOS, OSX, ...) and you are using _Clang_ as the 
     compiler, you want to deactivate the use of LLVM's libc++ by passing 
     `-Dpfasst_DISABLE_LIBCXX=ON` to the _CMake_ command line. (this is the default)

   * For example, to build with _Clang_ on a Linux system without the examples, the call to 
     _CMake_ looks like:

       cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ [-Dpfasst_DISABLE_LIBCXX=ON] -Dpfasst_BUILD_EXAMPLES=OFF ..

   The [Google Testing] and [Mocking Framework] are required for the unit tests, and are 
   automatically downloaded and built.

   The [FFTW3] library is required for (some) of the examples. It will be automatically downloaded 
   and built as well if it is not found on the system.

3. Run the unit tests (if unit tests have been built):

       make test

4. Run the examples (if they have been built):

       make run_example_all

   Or go to `./dist/examples` and look for the exetuable you want to run.

[Google Testing]: https://code.google.com/p/googletest/
[Mocking Framework]: https://code.google.com/p/googlemock/
[FFTW3]: http://fftw.org/
