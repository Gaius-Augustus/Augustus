#!/usr/bin/env python3

import argparse
import sys
import os
import importlib


# This script executes AUGUSTUS test cases based on the examples
# folder and compares the current results with reference results
# if the option --compare is set. It is expected that both results
# are identical for a successful test.
# This script must be called from "tests/short"!
# Python version 3.6 or higher is required for execution.

# Dict of available test modules.
# [0]: Argument to select the test case to be executed.
# [1]: Modulename of the selected test case.
testcases = {
    'examples' : 'examples.test_examples',
    'bam2hints' : 'auxprogs.bam2hints.test_bam2hints',
    'bam2wig' : 'auxprogs.bam2wig.test_bam2wig',
    'homGeneMapping' : 'auxprogs.homgenemapping.test_homgenemapping'
}

parser = argparse.ArgumentParser(description='Execute AUGUSTUS test cases.')
parser.add_argument('testcase',
                    action='store',
                    choices=list(testcases.keys()),
                    help='Testcase to execute.')
parser.add_argument('--mysql',
                    action='store_true',
                    help='cgp test cases are also executed with a MySQL database.')
parser.add_argument('--compare',
                    action='store_true',
                    help='Compare generated results with reference results.')
parser.add_argument('--html',
                    action='store_true',
                    help='Save diff results in html file.')
parser.add_argument('--clean',
                    action='store_true',
                    help='Remove all files created during the tests. If this option is set, no tests are executed.')
args = parser.parse_args()


def check_working_dir(clean):
    bindir = '../../bin/'
    wd = os.getcwd()
    if not (wd.endswith('tests/short')):
        errstr = 'Wrong working directory!' + '\n'
        errstr += 'This script must be called from "tests/short"!'
        sys.exit(errstr)
    if not clean and not (os.path.exists(f'{bindir}augustus')):
        errstr = 'Missing augustus binaries!' + '\n'
        errstr += f'The augustus binaries must be accessible in this path: "{bindir}"!'
        sys.exit(errstr)


if __name__ == '__main__':
    check_working_dir(args.clean)

    test_module = importlib.import_module(testcases[args.testcase])
    test_was_successful = True

    if args.clean:
        test_module.clean()
        sys.exit()
    else:
        test_was_successful = test_module.execute(args.compare, args.html, args.mysql)

    if test_was_successful:
        sys.exit(0)
    else:
        sys.exit(1)
