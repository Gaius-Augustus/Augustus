#!/usr/bin/env python3

import subprocess
import argparse
import datetime
import sys

# import util script from parent directory
sys.path.append('..')
import lr_util as util


parser = argparse.ArgumentParser(
    description='Python wrapper to execute single genome test case.')
parser.add_argument('-g', '--pathToGitRepo',
                    help='path to the Augustus Git repository.')
parser.add_argument('-e', '--evalDir', help='path to Eval script.')
args = parser.parse_args()


def execute_test():
    print('Executing test... ')

    start = datetime.datetime.now()
    subprocess.call(['bash', 'test_single.sh', '-e' + args.evalDir])
    end = datetime.datetime.now()

    return (end - start).total_seconds() / 60.0


def analyze_commit():
    info = util.commit_info(args.pathToGitRepo)
    exec_minutes = execute_test()
    util.store_additional_data(
        info[1], info[0], exec_minutes, 'output/additional_information.json')


if __name__ == '__main__':
    if args.pathToGitRepo is None:
        print('The path to the Augustus Git repository is required, please make use of --pathToGitRepo to pass the path...')
        sys.exit()

    if args.evalDir is None:
        print('The path eval script collection, please make use of --evalDir to pass the path...')
        sys.exit()

    analyze_commit()
