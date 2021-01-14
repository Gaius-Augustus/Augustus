#!/usr/bin/env python3

import unittest
import os
import sys
import shutil

if not (os.path.exists('tests')):
    sys.exit('Wrong working directory!\nTest files must be accessible in ./tests folder.')
sys.path.insert(0, '.') # make package tests.utils accessible

from tests.utils import aug_assertions
from tests.utils import aug_process
from tests.utils import aug_argparse
from tests.utils import aug_path


args = aug_argparse.getDefaultArgParser().parse_args()

pathname = os.path.dirname(sys.argv[0])
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
htmldir = os.path.join(pathname, 'output_html')
tmpdir = os.path.join(pathname, 'tmp')
exampledir = os.path.join(pathname, 'test_files')
exampletestfile = 'test.s.bam'
bindir = 'bin/'
bam2wigbin = f'{bindir}bam2wig'
default_wd = os.getcwd()


def clean(force_tmp_dir=True, force_html_dir=True, force_result_dir=True):
    """Remove empty directories or if forced"""
    aug_path.rmtree_if_exists(resultdir, force_result_dir)
    aug_path.rmtree_if_exists(tmpdir, force_tmp_dir)
    aug_path.rmtree_if_exists(htmldir, force_html_dir)

def check_working_dir(clean):
    errstr = ''
    if not (os.path.exists(exampledir)):
        errstr += 'Wrong example directory!' + '\n'
        errstr += f'The example files must be accessible in this path: "{exampledir}"!'
    if not clean and not (os.path.exists(bam2wigbin)):
        errstr += 'bam2wig binary not yet created!' + '\n'
        errstr += f'The augustus binaries must be accessible in this path: "{bindir}"!' + '\n'
        errstr += 'Please run make.'
    if (errstr):
        sys.exit(errstr)


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

        if os.path.isfile(os.path.join(tmpdir, f'{exampletestfile}.bai')):
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

        aug_process.execute(self, f'{bam2wigbin} {testfile} | grep -v track',  resfile)

        # compare results
        if TestBam2Wig.opt_compare:
            aug_assertions.assertEqualFiles(self, reffile, resfile, TestBam2Wig.opt_html, htmldir)
            os.remove(resfile)

    def test_bam2wig_region(self):
        os.chdir(default_wd)
        testfile = os.path.join(tmpdir, exampletestfile)
        reffile = os.path.join(referencedir, 'test.s.chr3L.wig')
        resfile = os.path.join(resultdir, 'test.s.chr3L.wig')

        aug_process.execute(self,
                      f'{bam2wigbin} -t "my_specified_track" -r chr3L {testfile} | grep -v track',
                      resfile)

        # compare results
        if TestBam2Wig.opt_compare:
            aug_assertions.assertEqualFiles(self, reffile, resfile, TestBam2Wig.opt_html, htmldir)
            os.remove(resfile)

def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestBam2Wig('test_bam2wig'))
    suite.addTest(TestBam2Wig('test_bam2wig_region'))

    return suite


if __name__ == '__main__':
    check_working_dir(args.clean)

    if args.clean:
        clean()
        sys.exit()

    TestBam2Wig.opt_compare = args.compare
    TestBam2Wig.opt_html = args.html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    if result.wasSuccessful():
        sys.exit()
    else:
        sys.exit(1)
