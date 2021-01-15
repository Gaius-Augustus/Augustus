#!/usr/bin/env python3

import argparse
import sys
import os

from examples import *

# This script executes AUGUSTUS test cases based on the examples
# folder and compares the current results with reference results
# if the option --compare is set. It is expected that both results
# are identical for a successful test.
# This script must be called from "tests/short"!
# Python version 3.6 or higher is required for execution.

testcases = [ 'examples', 'bam2wig', 'homeGeneMapping' ]

parser = argparse.ArgumentParser(description='Execute AUGUSTUS test cases.')
parser.add_argument('testcase',
                    action='store',
                    choices=testcases,
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

    # Remove only generated test files and do not execute test
    # cases if option --clean is set.
    if args.clean:
        clean_examples()
        sys.exit()

    test_was_successful = True
    if args.testcase == 'examples':
        test_was_successful = execute_examples(args.compare, args.html, args.mysql)
    
    if test_was_successful:
        sys.exit(0)
    else:
        sys.exit(1)
