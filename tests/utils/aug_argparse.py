#!/usr/bin/env python3

import argparse

def getDefaultArgParser():
    parser = argparse.ArgumentParser(description='Execute Augustus test cases.')
    parser.add_argument('--compare',
                        action='store_true',
                        help='Compare generated results with reference results.')
    parser.add_argument('--html',
                        action='store_true',
                        help='Save diff results in html file.')
    parser.add_argument('--clean',
                        action='store_true',
                        help='Remove all files created during the tests.')
    return parser
