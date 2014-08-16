#!/bin/sh
basepath=`pwd`

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

builddir=${1}
cd ${builddir}

rm -rf ${basepath}/coverage
mkdir -p ${basepath}/coverage

for testdir in `find ${basepath}/${builddir} -type d | grep -o 'tests/.*/.*\dir'`
do
  testname=`expr "$testdir" : '.*\(test_[a-zA-Z\-_]*\)\.dir'`
  echo "Gathering Coverage for ${testname}"
  cd $testdir
  lcov --zerocounters  --directory .
  cd ${basepath}/${builddir}
  ctest -R $testname
  cd $testdir
  lcov --directory . --capture --output-file ${testname}.info.tmp
  lcov --extract ${testname}.info.tmp "*${basepath}/include/**/*" --output-file ${testname}.info
  rm ${testname}.info.tmp
  cd ${basepath}/${builddir}
  if [[ -e all_tests.info ]]
  then
    lcov --add-tracefile all_tests.info --add-tracefile ${testdir}/${testname}.info --output-file all_tests.info
  else
    lcov --add-tracefile ${testdir}/${testname}.info --output-file all_tests.info
  fi
done

cd ${basepath}
genhtml --output-directory ./coverage \
  --demangle-cpp --num-spaces 2 --sort \
  --title "PFASST++ Test Coverage" --prefix ${basepath}/include \
  --function-coverage --branch-coverage --legend ${basepath}/${builddir}/all_tests.info
