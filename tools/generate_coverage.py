#!/usr/bin/env python
# coding=utf-8
"""
This script will aid in generating a test coverage report for PFASST++
including its examples.

A standard CPython 3.3 compatible Python interpreter with standard library
support is required. No additional modules.

Run it with argument `-h` for usage instructions.

.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""
from sys import version_info
# require at least Python 3.3
#  (because subprocess.DEVNULL)
assert(version_info[0] >= 3 and version_info[1] >= 3)


import argparse
import os
import os.path
from os.path import join
import shutil
import subprocess as sp
import re
import logging
from logging.config import dictConfig


dictConfig(
    {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'default': {
                'style': '{',
                'format': '[{levelname!s:<8s}] {message!s}'
            }
        },
        'handlers': {
            'console': {
                'class': 'logging.StreamHandler',
                'formatter': 'default'
            }
        },
        'root': {
            'handlers': ['console'],
            'level': 'DEBUG'
        }
    }
)


class Options(object):
    coverage_dir = ""
    build_dir = ""
    base_dir = ""
    with_examples = True
    tests = []
    example_tests = []
    tracefiles = []
    final_tracefile = ""


options = Options()
options.base_dir = ""


def is_lcov_available():
    """
    Check whether lcov is available.

    Returns:
        bool: True if lcov found, False if not
    """
    try:
        sp.check_call('lcov --version', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    except sp.CalledProcessError:
        logging.critical("lcov command not available. It is required.")
        return False

    try:
        sp.check_call('genhtml --version', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    except sp.CalledProcessError:
        logging.critical("genhtml command not available. It is required.")
        return False

    return True


def get_project_root():
    """
    Try to get the project root path and sets working directory accordingly
    """
    logging.info("Determine project root directory")
    curr_dir = os.path.abspath(os.path.curdir)
    logging.debug("Trying current path: %s" % curr_dir)

    include_access = os.access(join(curr_dir, "include"), os.R_OK)
    examples_access = os.access(join(curr_dir, "examples"), os.R_OK)
    if include_access and examples_access:
        logging.debug("Project root is: %s" % curr_dir)
        options.base_dir = curr_dir
    else:
        logging.warning("Probably called from within the tools dir. "
                        "This should work but is not recommended. "
                        "Trying parent directory as project root.")
        os.chdir("..")
        get_project_root()


def setup_and_init_options():
    help_string = "Note:\n" \
                  "This only works for builds made with GCC and the following CMake variables:\n" \
                  "    -Dpfasst_WITH_GCC_PROF=ON -Dpfasst_BUILD_TESTS=ON"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=help_string)
    parser.add_argument('-d', '--build-dir', required=True,
                        help="name of build directory containing a debug build with GCC and enabled profiling")
    parser.add_argument('-o', '--output', default='coverage',
                        help="output directory for generated coverage report")
    parser.add_argument('--no-examples', default=False, action='store_true',
                        help="whether to not run and include tests from the examples")
    parser.add_argument('--debug', default=False, action='store_true',
                        help="enables more verbose debugging output")

    _args = parser.parse_args()

    if _args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Debug mode enabled.")

    get_project_root()

    if not os.access(_args.build_dir, os.W_OK):
        error = "Given build path could not be found: %s" % _args.build_dir
        logging.critical(error)
        raise ValueError(error)
    options.build_dir = os.path.abspath(_args.build_dir)

    cache_path = join(options.build_dir, "CMakeCache.txt")
    if os.access(cache_path, os.W_OK):
        with_gcc_prof = False
        with_mpi = False
        with open(cache_path, 'r') as cache:
            for line in cache:
                if "pfasst_WITH_GCC_PROF:BOOL=ON" in line:
                    with_gcc_prof = True
                if "pfasst_WITH_MPI:BOOL=ON" in line:
                    with_mpi = True
    if not with_gcc_prof:
        error = "PFASST++ must be built with 'pfasst_WITH_GCC_PROF=ON'"
        logging.critical(error)
        raise RuntimeError(error)
    if with_mpi:
        logging.warning("Coverage analysis only functional for non-MPI builds")
        exit(0)

    if not os.access(_args.output, os.W_OK):
        logging.info("Output directory not found. Creating: %s" % _args.output)
        os.mkdir(_args.output)
    else:
        logging.warning("Clearing out output directory: %s" % _args.output)
        shutil.rmtree(_args.output)
        os.mkdir(_args.output)
    options.coverage_dir = os.path.abspath(_args.output)

    options.with_examples = not _args.no_examples
    if not options.with_examples:
        logging.debug("Not running and tracing tests from examples.")


def get_test_directories():
    """
    Find tests in the build_dir
    """
    logging.info("Looking for tests ...")

    test_dir = join(options.build_dir, 'tests')
    regex_name = re.compile('^.*\/(?P<test_name>test_[a-zA-Z\d\-_]+)\.dir$')
    regex_example = re.compile('^.*\/tests\/examples\/.*$')

    for root, dirs, files in os.walk(test_dir):
        match_name = regex_name.search(root)
        match_is_example = regex_example.search(root)
        is_example = match_is_example is not None
        if match_name is not None:
            testname = match_name.groupdict()['test_name']
            if  is_example:
                test_case = {
                    'path': root,
                    'name': testname,
                    'is_example': is_example
                }
                options.example_tests.append(test_case)
            else:
                test_case = {
                    'path': root,
                    'name': testname,
                    'is_example': is_example
                }
                options.tests.append(test_case)

    num_tests = len(options.tests)
    num_example_tests = len(options.example_tests)
    num_total = num_tests + num_example_tests
    logging.info("%d tests found" % num_total)
    logging.info("  %d general tests" % num_tests)
    if options.with_examples:
        logging.info("  %d tests for examples" % num_example_tests)


def run_test(path, name, is_example):
    logging.info("- %s" % name)
    logging.debug("Found in %s" % path)
    output_file = open(join(options.coverage_dir, name) + ".log", mode='a')
    logging.debug("Output log: %s" % output_file.name)

    def _print(s):
        print(s, file=output_file, flush=True)

    def _call(cmd):
        return sp.check_call(cmd, shell=True,
                             stdout=output_file, stderr=output_file)

    os.chdir(os.path.abspath(path))
    logging.debug("Deleting old tracing data ...")
    _print('### deleting old tracing data ...')
    _call('lcov --zerocounters --directory .')
    _print('### done.')

    os.chdir(options.build_dir)
    logging.debug("Running test ...")
    _print('### running test ...')
    _call('ctest -R %s' % name)
    _print('### done.')

    os.chdir(os.path.abspath(path))
    logging.debug("Capturing all tracing data ...")
    _print('### capturing all tracing data ...')
    _call('lcov --capture --directory . --output-file "%s.info.complete"' % name)
    _print('### done.')

    logging.debug("Removing unnecessary data ...")
    _print('### removing unnecessary data ...')
    try:
        _call('lcov --remove "%s.info.complete" "%s/include/pfasst/easylogging++.h" --output-file %s.info.prelim' % (name, options.base_dir, name))
    except sp.CalledProcessError as e:
        logging.warning(e)
    _print('### done.')

    logging.debug("Extracting interesting tracing data ...")
    _print('### extracting interesting tracing data ...')
    try:
        glob = join(options.base_dir, "src", "pfasst", "**", "*")
        _call('lcov --extract "%s.info.prelim" "%s" --output-file %s.info' % (name, glob, name))
        options.tracefiles.append("%s/%s.info" % (os.path.abspath(path), name))
    except sp.CalledProcessError as e:
        logging.warning(e)

    if is_example:
        logging.debug("This test belongs to an example, thus also covering examples code")
        try:
            glob = join(options.base_dir, "examples", "**", "*")
            _call('lcov --extract "%s.info.prelim" "%s" --output-file %s.info.example' % (name, glob, name))
            options.tracefiles.append("%s/%s.info.example" % (os.path.abspath(path), name))
        except sp.CalledProcessError as e:
            logging.warning(e)

    os.chdir(options.base_dir)
    output_file.close()


def run_tests():
    """
    Executes all tests with run_test
    """
    logging.info("Running general tests ...")
    for test in options.tests:
        run_test(**test)
    if options.with_examples:
        logging.info("Running tests for examples ...")
        for example in options.example_tests:
            run_test(**example)


def aggregate_tracefiles():
    logging.info("Aggregating %d tracefiles ..." % len(options.tracefiles))
    output_file = open('%s/aggegrating.log' % (options.coverage_dir,), mode='a')
    logging.debug("Output log: %s" % output_file.name)
    options.final_tracefile = "%s/all_tests.info" % options.coverage_dir
    for tracefile in options.tracefiles:
        logging.debug("- %s" % (tracefile))
        # skip empty tracefiles
        if not os.path.getsize(tracefile):
            continue
        print("### adding tracefile: %s" % (tracefile,), file=output_file, flush=True)
        if os.access(options.final_tracefile, os.W_OK):
            sp.check_call('lcov --add-tracefile "%s" --add-tracefile "%s" --output-file "%s"'
                          % (options.final_tracefile, tracefile, options.final_tracefile),
                          shell=True, stdout=output_file, stderr=output_file)
        else:
            sp.check_call('lcov --add-tracefile "%s" --output-file "%s"'
                          % (tracefile, options.final_tracefile),
                          shell=True, stdout=output_file, stderr=output_file)
        print("### done.", file=output_file, flush=True)
    output_file.close()


def generate_html():
    logging.info("Generating HTML report ...")
    output_path = '%s/generate_html.log' % (options.coverage_dir)
    output_file = open(output_path, mode='a')
    cmd = (
        'genhtml --output-directory {} --demangle-cpp --num-spaces 2 --sort '
        '--title "PFASST++ Test Coverage" --prefix "{}" --function-coverage '
        '--legend {}'
    ).format(options.coverage_dir, options.base_dir, options.final_tracefile)
    sp.check_call(cmd, shell=True, stdout=output_file, stderr=output_file)
    output_file.close()
    logging.info("Coverage report can be found in: file://%s/index.html" % options.coverage_dir)


if __name__ == "__main__":
    if not is_lcov_available():
        raise RuntimeError("Required commands could not be found.")
    setup_and_init_options()
    get_test_directories()
    run_tests()
    aggregate_tracefiles()
    generate_html()
