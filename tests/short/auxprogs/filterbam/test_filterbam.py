#!/usr/bin/env python3

import unittest
import os
import shutil

from utils import aug_assertions
from utils import aug_process
from utils import aug_path

pathname = 'auxprogs/filterbam/'
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
htmldir = os.path.join('output_html/')
tmpdir = os.path.join(pathname, 'tmp')
exampledir = os.path.join(pathname, 'test_files')
bindir = '../../bin/'
binfilterbam = f'{bindir}filterBam'
default_wd = os.getcwd()


def clean(force_tmp_dir=True, force_html_dir=True, force_result_dir=True):
    """Remove empty directories or if forced"""
    aug_path.rmtree_if_exists(resultdir, force_result_dir)
    aug_path.rmtree_if_exists(tmpdir, force_tmp_dir)
    aug_path.rmtree_if_exists(htmldir, force_html_dir)


class TestFilterBam(unittest.TestCase):

    opt_compare = False
    opt_html = False

    @classmethod
    def setUpClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=True, force_result_dir=True)
        aug_path.mkdir_if_not_exists(resultdir)
        aug_path.mkdir_if_not_exists(tmpdir)


    @classmethod
    def tearDownClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=False, force_result_dir=False)

    def test_compare(self, refsamfile, resbamfile, ressamfile):
        os.chdir(default_wd)
        
        # compare results
        if TestFilterBam.opt_compare:
            aug_assertions.assertEqualFiles(self, refsamfile, ressamfile, TestFilterBam.opt_html, htmldir)
            os.remove(resbamfile)
            os.remove(ressamfile)

    def test_case(self, testfilename, refsamfilename, filterbamoptions):
        os.chdir(default_wd)
        testfile = os.path.join(exampledir, testfilename)
        refsamfile = os.path.join(referencedir, refsamfilename)
        resbamfile = os.path.join(resultdir, f'{refsamfilename}.bam')
        ressamfile = os.path.join(resultdir, refsamfilename)
        
        aug_process.execute(self, f'{binfilterbam} --in {testfile} --out {resbamfile} {filterbamoptions} ')
        
        if "--no-PG" in aug_process.execute(None, 'samtools view -?'):
            aug_process.execute(None, f'samtools view -h --no-PG {resbamfile} -o {ressamfile} ')
        else: 
            aug_process.execute(None, f'samtools view -h {resbamfile} -o {ressamfile} ')
        
        self.test_compare(refsamfile, resbamfile, ressamfile)


    def test_no_options(self):
        # unmapped fragments without coordinates are removed, default: minCover = 80, minId = 92
        self.test_case('example_single.bam', 'example_single.no_options.sam', '')

    def test_no_minId_no_minCover(self):
        # unmapped fragments without coordinates are removed
        self.test_case('example_single.bam', 'example_single.no_minId_no_minCover.sam', '--minId=0 --minCover=0')

    def test_no_minId(self):
        # unmapped fragments without coordinates are removed
        self.test_case('example_single.bam', 'example_single.no_minId.sam', '--minId=0')

    def test_no_minCover(self):
        # unmapped fragments without coordinates are removed
        self.test_case('example_single.bam', 'example_single.no_minCover.sam', '--minCover=0')

    def test_minCover99(self):
        # unmapped fragments without coordinates are removed
        self.test_case('example_single.bam', 'example_single.minCover99.sam', '--minCover=99')

    def test_uniq(self):
        self.test_case('example_single.bam', 'example_single.uniq.sam', '--uniq ')

    def test_uniqThresh(self):
        self.test_case('example_single.bam', 'example_single.uniqThresh1_1.sam', '--uniq --uniqThresh=1.1')

    def test_best(self):
        self.test_case('example_single.bam', 'example_single.best.sam', '--best')

    def test_noIntrons(self):
        self.test_case('example_single.bam', 'example_single.noIntrons.sam', '--noIntrons')

    def test_noIntrons_insertLimit(self):
        self.test_case('example_single.bam', 'example_single.noIntrons_insertLimit15.sam', '--noIntrons --insertLimit=15')

    def test_not_paired(self):
        self.test_case('example_single.bam', 'example_single.not_paired.sam', '--paired')

    def test_paired(self):
        self.test_case('example_paired.bam', 'example_paired.paired.sam', '--paired --minId=0 --minCover=0')

    def test_paired_pairwiseAlignments(self):
        self.test_case('example_paired.bam', 'example_paired.paired_pairwiseAlignments.sam', '--paired --pairwiseAlignments --minId=0 --minCover=0')

    def test_pairwiseAlignments_2(self):
        self.test_case('example_pairwise.bam', 'example_pairwise.paired_pairwiseAlignments.sam', '--paired --pairwiseAlignments --minId=97 --minCover=98')

    def test_paired_maxIntronLen(self):
        self.test_case('example_paired.bam', 'example_paired.paired_maxIntronLen100.sam', '--paired --maxIntronLen=100 --minId=0 --minCover=0')

    def test_paired_pairBedFile(self):
        refbedfile = os.path.join(referencedir, 'pairBedFile.bed')
        resbedfile = os.path.join(resultdir, 'pairBedFile.bed')
        self.test_case('example_paired.bam', 'example_paired.paired.sam', f'--paired --minId=0 --minCover=0 --pairBedFile={resbedfile}')
        if TestFilterBam.opt_compare:
            aug_assertions.assertEqualFiles(self, refbedfile, resbedfile, TestFilterBam.opt_html, htmldir)
            os.remove(resbedfile)

    def test_commonGeneFile(self):
        reffile = os.path.join(referencedir, 'commonGeneFile.txt')
        resfile = os.path.join(resultdir, 'commonGeneFile.txt')
        self.test_case('example_single.bam', 'example_single.best.sam', f'--best --commonGeneFile={resfile}')
        if TestFilterBam.opt_compare:
            aug_assertions.assertEqualFiles(
                self, reffile, resfile, TestFilterBam.opt_html, htmldir)
            os.remove(resfile)

    def test_calmd_seq_with_equalsigns_no_options(self):
        # unmapped fragments without coordinates are removed, default: minCover = 80, minId = 92
        self.test_case('example_pairwise_calmd_seq_with_equalsigns.bam', 'example_pairwise_calmd_seq_with_equalsigns.no_options.sam', '')

    def test_calmd_no_options(self):
        # unmapped fragments without coordinates are removed, default: minCover = 80, minId = 92
        self.test_case('example_pairwise_calmd.bam', 'example_pairwise_calmd.no_options.sam', '')
        
    def test_calmd_seq_with_equalsigns_uniq(self):
        self.test_case('example_pairwise_calmd_seq_with_equalsigns.bam', 'example_pairwise_calmd_seq_with_equalsigns.uniq.sam', '--uniq --')

    def test_calmd_uniq(self):
        self.test_case('example_pairwise_calmd.bam', 'example_pairwise_calmd.uniq.sam', '--uniq ')

    def test_calmd_seq_with_equalsigns_best(self):
        self.test_case('example_pairwise_calmd_seq_with_equalsigns.bam', 'example_pairwise_calmd_seq_with_equalsigns.best.sam', '--best')

    def test_calmd_best(self):
        self.test_case('example_pairwise_calmd.bam', 'example_pairwise_calmd.best.sam', '--best')

    def test_calmd_seq_with_equalsigns_paired_pairwiseAlignments(self):
        self.test_case('example_pairwise_calmd_seq_with_equalsigns.bam', 'example_pairwise_calmd_seq_with_equalsigns.paired_pairwiseAlignments.sam', '--paired --pairwiseAlignments --minId=0 --minCover=0')

    def test_calmd_paired_pairwiseAlignments(self):
        self.test_case('example_pairwise_calmd.bam', 'example_pairwise_calmd.paired_pairwiseAlignments.sam', '--paired --pairwiseAlignments --minId=0 --minCover=0')

def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestFilterBam('test_no_options'))
    suite.addTest(TestFilterBam('test_no_minId_no_minCover'))
    suite.addTest(TestFilterBam('test_no_minId'))
    suite.addTest(TestFilterBam('test_no_minCover'))
    suite.addTest(TestFilterBam('test_minCover99'))
    suite.addTest(TestFilterBam('test_uniq'))
    suite.addTest(TestFilterBam('test_uniqThresh'))
    suite.addTest(TestFilterBam('test_best'))
    suite.addTest(TestFilterBam('test_noIntrons'))
    suite.addTest(TestFilterBam('test_noIntrons_insertLimit'))
    suite.addTest(TestFilterBam('test_not_paired'))
    suite.addTest(TestFilterBam('test_paired'))
    suite.addTest(TestFilterBam('test_paired_pairwiseAlignments'))
    suite.addTest(TestFilterBam('test_pairwiseAlignments_2'))
    suite.addTest(TestFilterBam('test_paired_maxIntronLen'))
    suite.addTest(TestFilterBam('test_paired_pairBedFile'))
    suite.addTest(TestFilterBam('test_commonGeneFile'))
    suite.addTest(TestFilterBam('test_calmd_seq_with_equalsigns_no_options'))
    suite.addTest(TestFilterBam('test_calmd_no_options'))
    suite.addTest(TestFilterBam('test_calmd_seq_with_equalsigns_uniq'))
    suite.addTest(TestFilterBam('test_calmd_uniq'))
    suite.addTest(TestFilterBam('test_calmd_seq_with_equalsigns_best'))
    suite.addTest(TestFilterBam('test_calmd_best'))
    suite.addTest(TestFilterBam('test_calmd_seq_with_equalsigns_paired_pairwiseAlignments'))
    suite.addTest(TestFilterBam('test_calmd_paired_pairwiseAlignments'))

    return suite


def execute(compare, html, mysql):
    TestFilterBam.opt_compare = compare
    TestFilterBam.opt_html = html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    return result.wasSuccessful()
