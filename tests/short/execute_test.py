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


if __name__ == '__main__':
    #TODO: clean

    if args.testcase == 'examples':
        execute_examples(args.compare, args.html, args.mysql)
