#!/usr/bin/env python3

import sys
import os
import difflib


def compare_folder(reffolder, currentfolder, html=False):
    if not os.path.isdir(reffolder):
        return 'Reference folder not found: ' + reffolder
    if not os.path.isdir(currentfolder):
        return 'Current folder not found: ' + currentfolder

    res = ''

    for subdir, dirs, files in os.walk(reffolder):
        for file in files:
            if not subdir.endswith('/'):
                subdir += '/'
            currentdir = subdir.replace(reffolder, currentfolder)

            diff = compare_files(subdir + file, currentdir + file, html=html)
            res += diff

    return res


def compare_files(reffile, currentfile, html=False):
    if not os.path.isfile(reffile):
        return 'Reference file not found: ' + reffile
    if not os.path.isfile(currentfile):
        return 'Current file not found: ' + currentfile

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

    filename = create_html_filename(reffile)
    with open(filename, 'w') as html:
        html.write(diff)


def create_html_filename(name):
    html_name = ''
    elements = []
    parts = name.split('/')

    for p in reversed(parts):
        elements.append(p)
        if 'test_' in p:
            break

    for e in reversed(elements):
        html_name += e + '_'

    return html_name + 'diff.html'


if __name__ == '__main__':
    reference = sys.argv[1]
    curent = sys.argv[2]
    #res = compare_files(reference, curent)
    res = compare_folder(reference, curent, html=True)
    print(res)
