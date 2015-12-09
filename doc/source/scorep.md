# Instrumentation with Score-P and usage of extra compiler wrappers           {#page_scorep}

The _CMake_ build system of _PFASST++_ offers the option `pfasst_WITH_EXTRA_WRAPPER` that if enabled will create to shell scripts called `cc_wrapper.sh` and `cxx_wrapper.sh` in the build directory.These can be used to specify a compiler wrapper. Therefore the environment variable `PREP` is used. All compiler calls of the build system will result in calls to `$PREP <original-compiler>`. This will be demonstrated with the example of Score-P a scalable performance and communication analysis tool.

1. Call cmake as usual but pass the argument `-Dpfasst_WITH_EXTRA_WRAPPER=ON` option. There should now be two shell scripts in the build directory.

2. When compiling use `PREP=scorep make` instead of `make` otherwise you will compile as usual. The binaries are now instrumented using Score-P.

Note on linking: On some architectures you only want to add Score-P during the link step. Therefore first compile without setting `PREP` which will result in an unistrumented object and binary. Than delete the binary and compile with `PREP=scorep`. Because the object files have not been deleted only the binaries will be linked (now with Score-P support).
