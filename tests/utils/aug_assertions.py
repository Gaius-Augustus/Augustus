#!/usr/bin/env python3

from tests.examples import aug_comparator as comp

def assertEqualFiles(testcase, reffile, resfile, html=False, htmloutputfolder=None):
    """
    Assert that the specified files are equal
    :param testcase: the calling unittest.TestCase
    :param reffile: the reference file
    :param resfile: the file created by the test
    :param html: if true and the files differ a report is written in html format
    :param htmloutputfolder: path for the html diff file, if None a default
           path is used
    :return:
    """
    diff = comp.compare_files(reffile,
                              resfile,
                              html=html or htmloutputfolder is not None,
                              outputfolder=htmloutputfolder)
    testcase.assertEqual(diff, '')
