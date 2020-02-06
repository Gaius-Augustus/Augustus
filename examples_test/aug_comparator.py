#!/usr/bin/env python3

import sys
import difflib


def compare_files(reffile, currentfile, html=False):
    with open(reffile) as ff:
        reflines = ff.readlines()
    with open(currentfile) as tf:
        currentlines = tf.readlines()

    diff = ''.join(
        difflib.context_diff(reflines, currentlines, reffile, currentfile,
                             n=0))

    if html and diff:
        generate_html(reflines, currentlines, reffile, currentfile)

    return diff


def generate_html(reflines, currentlines, reffile, currentfile):
    diff = difflib.HtmlDiff().make_file(reflines, currentlines, reffile,
                                        currentfile)

    with open('diff.html', 'w') as html:
        html.write(diff)


if __name__ == '__main__':
    reference = sys.argv[1]
    curent = sys.argv[2]
    res = compare_files(reference, curent)
    print(res)
