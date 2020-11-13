#!/usr/bin/env python3

import os
import sys


def pred(sourcepath, targetpath):
    startwith = '# ----- prediction'
    ignores = []
    ignoredlines = search_for_lines_to_ignore(sourcepath, startwith, ignores)
    filter(sourcepath, targetpath, ignoredlines)


def cgp(sourcepath, targetpath):
    startwith = '#----- prediction on'  # this case could be removed by a small change in augustus
    ignores = []
    ignoredlines = search_for_lines_to_ignore(sourcepath, startwith, ignores)
    filter(sourcepath, targetpath, ignoredlines)


def eval(sourcepath, targetpath):
    startwith = '# ----- sequence number'
    ignores = ['# total time:']
    ignoredlines = search_for_lines_to_ignore(sourcepath, startwith, ignores)
    filter(sourcepath, targetpath, ignoredlines)


def cgp_out(sourcepath, targetpath):
    startwith = ''
    ignores = ['# total time:']
    ignoredlines = search_for_lines_to_ignore(sourcepath, startwith, ignores, startblock=False)
    filter(sourcepath, targetpath, ignoredlines)


def search_for_lines_to_ignore(sourcepath, startwith, ignores, startblock=True):
    ignoredlines = []

    # startblock
    if startblock:
        with open(sourcepath) as fp:
            line = fp.readline()
            lno = 1
            while line:
                if startwith in line.strip():
                    break
                line = fp.readline()
                ignoredlines.append(lno)
                lno += 1

    # other lines to ignore
    if len(ignores) > 0:
        with open(sourcepath) as fp:
            line = fp.readline()
            lno = 1
            while line:
                if any(s in line.strip() for s in ignores):
                    ignoredlines.append(lno)

                line = fp.readline()
                lno += 1

    return ignoredlines


def filter(sourcepath, targetpath, ignoredlines):
    if os.path.exists(targetpath):
        os.remove(targetpath)

    with open(targetpath, "w") as file_out:
        with open(sourcepath) as fp:
            line = fp.readline()
            lno = 1
            while line:
                if lno not in ignoredlines:
                    file_out.write(line.strip() + '\n')
                line = fp.readline()
                lno += 1


if __name__ == '__main__':
    source = sys.argv[1]
    target = sys.argv[2]
    pred(source, target)

#TODO: add option to remove source after target is created??
