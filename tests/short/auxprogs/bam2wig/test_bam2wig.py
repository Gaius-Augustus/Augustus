#!/usr/bin/env python3

import unittest
import os
import shutil

from utils import aug_assertions
from utils import aug_process
from utils import aug_path

pathname = 'auxprogs/bam2wig/'
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
htmldir = os.path.join('output_html/')
tmpdir = os.path.join(pathname, 'tmp')
exampledir = os.path.join(pathname, 'test_files')
exampletestfile = 'test.s.bam'
bindir = '../../bin/'
bam2wigbin = bindir + 'bam2wig'
default_wd = os.getcwd()


def clean(force_tmp_dir=True, force_html_dir=True, force_result_dir=True):
    """Remove empty directories or if forced"""
    aug_path.rmtree_if_exists(resultdir, force_result_dir)
    aug_path.rmtree_if_exists(tmpdir, force_tmp_dir)
    aug_path.rmtree_if_exists(htmldir, force_html_dir)


class TestBam2Wig(unittest.TestCase):

    opt_compare = False
    opt_html = False

    @classmethod
    def setUpClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=True, force_result_dir=True)
        aug_path.mkdir_if_not_exists(resultdir)
        aug_path.mkdir_if_not_exists(tmpdir)

        inputfile = os.path.join(exampledir, exampletestfile)
        testfile = os.path.join(tmpdir, exampletestfile)
        shutil.copyfile(inputfile, testfile)

        os.chdir(tmpdir)
        aug_process.execute(None, 'samtools index ' + exampletestfile)

        if os.path.isfile(os.path.join(tmpdir, exampletestfile+'.bai')):
            assert False, 'Test file not indexed.'

    @classmethod
    def tearDownClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=False, force_result_dir=False)

    def test_bam2wig(self):
        os.chdir(default_wd)
        testfile = os.path.join(tmpdir, exampletestfile)
        reffile = os.path.join(referencedir, 'test.s.wig')
        resfile = os.path.join(resultdir, 'test.s.wig')

        aug_process.execute(
            self, bam2wigbin + " " + testfile + ' | grep -v track',  resfile)

        # compare results
        if TestBam2Wig.opt_compare:
            aug_assertions.assertEqualFiles(
                self, reffile, resfile, TestBam2Wig.opt_html, htmldir)
            os.remove(resfile)

    def test_bam2wig_region(self):
        os.chdir(default_wd)
        testfile = os.path.join(tmpdir, exampletestfile)
        reffile = os.path.join(referencedir, 'test.s.chr3L.wig')
        resfile = os.path.join(resultdir, 'test.s.chr3L.wig')

        aug_process.execute(self,
                            bam2wigbin + ' -t "my_specified_track" -r chr3L ' + testfile + ' | grep -v track',
                            resfile)

        # compare results
        if TestBam2Wig.opt_compare:
            aug_assertions.assertEqualFiles(
                self, reffile, resfile, TestBam2Wig.opt_html, htmldir)
            os.remove(resfile)


def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestBam2Wig('test_bam2wig'))
    suite.addTest(TestBam2Wig('test_bam2wig_region'))

    return suite


def execute(compare, html, mysql):
    TestBam2Wig.opt_compare = compare
    TestBam2Wig.opt_html = html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    return result.wasSuccessful()
