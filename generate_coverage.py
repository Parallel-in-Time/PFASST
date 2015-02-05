# coding=utf-8
"""
.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""

from sys import version_info
# require at least Python 3.3
assert(version_info[0] >= 3 and version_info[1] >= 3)


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
            'level': 'INFO'
        }
    }
)


import subprocess as sp
try:
    sp.check_call('lcov --version', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
except sp.CalledProcessError as err:
    logging.critical("lcov command not found. It is required.")
    raise err

try:
    sp.check_call('genhtml --version', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
except sp.CalledProcessError as err:
    logging.critical("genhtml command not found. It is required.")
    raise err


import argparse
import os
import os.path
import shutil
import re


class Options(object):
    coverage_dir = ""
    build_dir = ""
    base_dir = ""
    tests = []
    tracefiles = []
    final_tracefile = ""


options = Options()
options.base_dir = os.path.abspath(os.path.curdir)


def setup_and_init_options():
    help_string = "Note:\n" \
                  "This only works for builds made with GCC and the following CMake variables:\n" \
                  "    -Dpfasst_WITH_GCC_PROF=ON -Dpfasst_BUILD_TESTS=ON"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=help_string)
    parser.add_argument('-d', '--build_dir', required=True,
                        help="name of build directory containing a debug build with GCC and enabled profiling")
    parser.add_argument('-o', '--output', default='coverage',
                        help="output directory for generated coverage report")

    _args = parser.parse_args()
    if not os.access(_args.build_dir, os.W_OK):
        logging.critical("Given build path could not be found: %s" % _args.build_dir)
        raise ValueError("Given build path could not be found: %s" % _args.build_dir)
    options.build_dir = os.path.abspath(_args.build_dir)

    if not os.access(_args.output, os.W_OK):
        logging.info("output directory not found. creating: %s" % _args.output)
        os.mkdir(_args.output)
    else:
        logging.warning("clearing out output directory: %s" % _args.output)
        shutil.rmtree(_args.output)
        os.mkdir(_args.output)
    options.coverage_dir = os.path.abspath(_args.output)


def get_test_directories():
    logging.info("looking for tests ...")
    for root, dirs, files in os.walk(options.build_dir + '/tests'):
        match_name = re.search('^.*/(?P<test_name>test_[a-zA-Z\-_]+)\.dir$', root)
        match_is_example = re.search('^.*/tests/examples/.*$', root)
        is_example = match_is_example is not None
        if match_name is not None:
            testname = match_name.groupdict()['test_name']
            options.tests.append({'path': root, 'name': testname, 'is_example': is_example})
    logging.info("%d tests found" % len(options.tests))


def run_test(path, name, is_example):
    logging.info("running %s" % name)
    logging.debug("in %s" % path)
    output_file = open('%s/%s.log' % (options.coverage_dir, name), mode='a')

    os.chdir(os.path.abspath(path))
    logging.debug("deleting old tracing data ...")
    print('### deleting old tracing data ...', file=output_file, flush=True)
    sp.check_call('lcov --zerocounters --directory .', shell=True, stdout=output_file, stderr=output_file)
    print('### done.', file=output_file, flush=True)

    os.chdir(options.build_dir)
    logging.debug("running test ...")
    print('### running test ...', file=output_file, flush=True)
    sp.check_call('ctest -R %s' % name, shell=True, stdout=output_file, stderr=output_file)
    print('### done.', file=output_file, flush=True)

    os.chdir(os.path.abspath(path))
    logging.debug("capturing all tracing data ...")
    print('### capturing all tracing data ...', file=output_file, flush=True)
    sp.check_call('lcov --capture --directory . --output-file "%s.info.complete"' % name,
                  shell=True, stdout=output_file, stderr=output_file)
    print('### done.', file=output_file, flush=True)

    logging.debug("removing unnecessary data ...")
    print('### removing unnecessary data ...', file=output_file, flush=True)
    try:
        sp.check_call('lcov --remove "%s.info.complete" "%s/include/pfasst/easylogging++.h" --output-file %s.info.prelim'
                      % (name, options.base_dir, name),
                      shell=True, stdout=output_file, stderr=output_file)
    except sp.CalledProcessError as e:
        logging.warning(e)
    print('### done.', file=output_file, flush=True)

    logging.debug("extracting interesting tracing data ...")
    print('### extracting interesting tracing data ...', file=output_file, flush=True)
    try:
        sp.check_call('lcov --extract "%s.info.prelim" "*%s/include/**/*" --output-file %s.info'
                      % (name, options.base_dir, name),
                      shell=True, stdout=output_file, stderr=output_file)
        options.tracefiles.append("%s/%s.info" % (os.path.abspath(path), name))
    except sp.CalledProcessError as e:
        logging.warning(e)
    if is_example:
        logging.debug("this test belongs to an example, thus also covering examples code")
        try:
            sp.check_call('lcov --extract "%s.info.prelim" "*%s/examples/**/*" --output-file %s.info.example'
                          % (name, options.base_dir, name),
                          shell=True, stdout=output_file, stderr=output_file)
            options.tracefiles.append("%s/%s.info.example" % (os.path.abspath(path), name))
        except sp.CalledProcessError as e:
            logging.warning(e)
    print('### done.', file=output_file, flush=True)

    os.chdir(options.base_dir)
    output_file.close()


def aggregate_tracefiles():
    logging.info("aggregating %d tracefiles ..." % len(options.tracefiles))
    output_file = open('%s/aggegrating.log' % (options.coverage_dir,), mode='a')
    options.final_tracefile = "%s/all_tests.info" % options.coverage_dir
    for tracefile in options.tracefiles:
        logging.debug("adding tracefile: %s" % (tracefile))
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
    logging.info("generating HTML report ...")
    output_file = open('%s/generate_html.log' % (options.coverage_dir,), mode='a')
    sp.check_call('genhtml --output-directory %s --demangle-cpp --num-spaces 2 --sort '
                  '--title "PFASST++ Test Coverage" --prefix "%s" --function-coverage --legend "%s"'
                  % (options.coverage_dir, options.base_dir, options.final_tracefile),
                  shell=True, stdout=output_file, stderr=output_file)
    output_file.close()
    logging.info("coverage report can be found in: %s" % options.coverage_dir)


if __name__ == "__main__":
    setup_and_init_options()
    get_test_directories()
    for test in options.tests:
        run_test(**test)
    aggregate_tracefiles()
    generate_html()
