# Building and Installing

## Prerequesites

* a recent C++ compiler, which supports most of the C++11 features
* the [CMake](http://cmake.org/) build tool


## Obtaining the Sources

Just use _Git_ to clone the repository.


## Building

1. create a directory for out-of-source builds (usually, people are creative and name it `build`) in the source tree
2. decide on what you want to build
   
   By default everythin gets build, i.e. tests and examples.
   
   To deactive the tests pass `-Dpfasst_BUILD_TESTS=OFF` to the _CMake_ command line.
   
   To deactive the examples, pass `-Dpfasst_BUILD_EXAMPLES=OFF` to the _CMake_ command line.
   
   If you are not on a _Darvin_ system (i.e. MacOS, OSX, ...) and you are using _Clang_ as the
   compiler, you want to deactive the use of LLVM's libc++ by passing `-Dpfasst_DISABLE_LIBCXX=ON`
   to the _CMake_ command line.
   
   In the end, for a build with _Clang_ on a Linux system without the examples, the call to _CMake_
   looks like:
   
       cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -Dpfasst_DISABLE_LIBCXX=ON -Dpfasst_BUILD_EXAMPLES=OFF ..

   As the [Google Testing] and [Mocking Framework] are required for the unit tests, those are automatically
   downloaded and build.
   
   In case, you enable the examples, [FFTW3] gets downloaded and build as well.
3. run the unit tests (if unit tests have been built)

       make test
4. run the examples (if those have been built). e.g.:

       ./examples/advection

[Google Testing]: https://code.google.com/p/googletest/
[Mocking Framework]: https://code.google.com/p/googlemock/
[FFTW3]: http://fftw.org/
