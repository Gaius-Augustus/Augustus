#!/usr/bin/env python3

import unittest
import os
import shutil

from utils import aug_assertions
from utils import aug_process
from utils import aug_path

pathname = 'auxprogs/bam2hints/'
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
htmldir = os.path.join('output_html/')
tmpdir = os.path.join(pathname, 'tmp')
exampledir = os.path.join(pathname, 'test_files')
testfilename = 'test.s.sorted.bam'
reffilename = 'test.intron.gff'

bindir = '../../bin/'
bam2hintsbin = f'{bindir}bam2hints'
default_wd = os.getcwd()


def clean(force_tmp_dir=True, force_html_dir=True, force_result_dir=True):
    """Remove empty directories or if forced"""
    aug_path.rmtree_if_exists(resultdir, force_result_dir)
    aug_path.rmtree_if_exists(tmpdir, force_tmp_dir)
    aug_path.rmtree_if_exists(htmldir, force_html_dir)


class TestBam2Hints(unittest.TestCase):

    opt_compare = False
    opt_html = False

    @classmethod
    def setUpClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=True, force_result_dir=True)
        aug_path.mkdir_if_not_exists(resultdir)
        aug_path.mkdir_if_not_exists(tmpdir)

        inputfile = os.path.join(exampledir, testfilename)
        testfile = os.path.join(tmpdir, testfilename)
        shutil.copyfile(inputfile, testfile)
        

    @classmethod
    def tearDownClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=False, force_result_dir=False)

    def test_bam2hints(self):
        os.chdir(default_wd)
        testfile = os.path.join(tmpdir, testfilename)
        reffile = os.path.join(referencedir, reffilename)
        resfile = os.path.join(resultdir, reffilename)

        aug_process.execute(
            self, f'{bam2hintsbin} --in={testfile} --out={resfile}')

        # compare results
        if TestBam2Hints.opt_compare:
            aug_assertions.assertEqualFiles(
                self, reffile, resfile, TestBam2Hints.opt_html, htmldir)
            os.remove(resfile)

    def test_bam2hints_pipes(self):
        os.chdir(default_wd)
        testfile = os.path.join(tmpdir, testfilename)
        reffile = os.path.join(referencedir, reffilename)
        resfile = os.path.join(resultdir, reffilename)

        aug_process.execute(
            self, f'cat {testfile} | {bam2hintsbin}', resfile)
        
        # compare results
        if TestBam2Hints.opt_compare:
            aug_assertions.assertEqualFiles(
                self, reffile, resfile, TestBam2Hints.opt_html, htmldir)
            os.remove(resfile)


def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestBam2Hints('test_bam2hints'))
    suite.addTest(TestBam2Hints('test_bam2hints_pipes'))

    return suite


def execute(compare, html, mysql):
    TestBam2Hints.opt_compare = compare
    TestBam2Hints.opt_html = html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    return result.wasSuccessful()
