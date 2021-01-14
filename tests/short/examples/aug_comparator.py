#!/usr/bin/env python3

import sys
import os
import difflib


def compare_folder(reffolder, currentfolder, html=False, outputfolder='output_html/'):
    if not os.path.isdir(reffolder):
        return 'Reference folder not found: ' + reffolder
    if not os.path.isdir(currentfolder):
        return 'Current folder not found: ' + currentfolder

    res = ''

    for subdir, dirs, files in os.walk(reffolder):
        for file in files:
            currentdir = subdir.replace(reffolder, currentfolder)

            diff = compare_files(os.path.join(subdir, file),
                                 os.path.join(currentdir, file),
                                 html=html, outputfolder=outputfolder)
            res += diff

    return res


def compare_files(reffile, currentfile, html=False, outputfolder='output_html/'):
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
        generate_html(reflines, currentlines, reffile,
                      currentfile, outputfolder)

    return diff


def generate_html(reflines, currentlines, reffile, currentfile, outputfolder):
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    diff = difflib.HtmlDiff().make_file(reflines, currentlines, reffile,
                                        currentfile)

    filename = create_html_filename(reffile)
    with open(os.path.join(outputfolder, filename), 'w') as html:
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
    current = sys.argv[2]
    print(compare_files(reference, current))
