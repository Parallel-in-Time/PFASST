#!/bin/sh
basepath="`pwd`"

function print_help {
  echo "#######################################################################"
  echo "###               Generation of Test Coverage Report                ###"
  echo "#                                                                     #"
  echo "# This only works for builds made with GCC and the following CMake    #"
  echo "# variables:                                                          #"
  echo "#   -Dpfasst_WITH_GCC_PROF=ON -Dpfasst_BUILD_TESTS=ON                 #"
  echo "#                                                                     #"
  echo "# First (and only) parameter must be the name of the build directory  #"
  echo "#                                                                     #"
  echo "# Example:                                                            #"
  echo "#   ./generate_coverage.sh build_gcc                                  #"
  echo "#                                                                     #"
  echo "#######################################################################"
  return 0
}


if [[ $# -ne 1 ]]
then
  print_help
  echo "ERROR: Please name the build directory as the first parameter."
  exit -1
fi

builddir="${basepath}/${1}"
if [[ ! -d "${builddir}" ]]
then
  echo "ERROR: Given build directory '${bilddir}' does not exist."
  exit -1
fi
cd "${builddir}"

coveragedir="${basepath}/coverage"
rm -rf "${coveragedir}"
mkdir -p "${coveragedir}"

for testdir in `find ${builddir} -type d | grep -o 'tests/.*/.*\dir'`
do
  testname=`expr "${testdir}" : '.*\(test_[a-zA-Z\-_]*\)\.dir'`
  echo "# ${testname}"
  cd "${testdir}"
  echo "    Deleting old tracing data"
  lcov --zerocounters  --directory . \
       > "${coveragedir}/log_${testname}_zero.log" 2&>1

  cd "${builddir}"
  echo "    Running test and generating tracing data"
  ctest -R ${testname} \
        > "${coveragedir}/log_${testname}_test.log"
  if [[ $? -ne 0 ]]
  then
    echo "ERROR: test '${testname}' failed."
    echo "       See '${coveragedir}/log_${testname}_test.log' for details."
    cd ${basepath}
    exit -1
  fi

  cd "${testdir}"
  echo "    Capturing all tracing data"
  lcov --directory . --capture --output-file "${testname}.info.tmp" \
       > "${coveragedir}/log_${testname}_capture.log"

  echo "    Extracting interesting tracing data"
  lcov --extract "${testname}.info.tmp" "*${basepath}/include/**/*" \
       --output-file "${testname}.info" \
       > "${coveragedir}/log_${testname}_extract.log"
  rm "${testname}.info.tmp"

  cd "${builddir}"
  echo "    Aggegrating coverage data"
  if [[ -e "${coveragedir}/all_tests.info" ]]
  then
    lcov --add-tracefile "${coveragedir}/all_tests.info" \
         --add-tracefile "${testdir}/${testname}.info" \
         --output-file "${coveragedir}/all_tests.info" \
         > "${coveragedir}/log_${testname}_aggregate.log"
  else
    lcov --add-tracefile "${testdir}/${testname}.info" \
         --output-file "${coveragedir}/all_tests.info" \
         > "${coveragedir}/log_${testname}_aggregate.log"
  fi
done

cd ${basepath}
echo "# Generating HTML report"
genhtml --output-directory "${coveragedir}" \
        --demangle-cpp --num-spaces 2 --sort \
        --title "PFASST++ Test Coverage" --prefix "${basepath}/include" \
        --function-coverage --branch-coverage --legend \
        "${coveragedir}/all_tests.info" \
        > "${coveragedir}/log_generate.log"

echo "# Coverage report can be found in ${coveragedir}/index.html"
echo "# Everything done."
