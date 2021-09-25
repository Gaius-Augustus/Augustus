#!/usr/bin/env python3

import unittest
import os

from utils import aug_assertions
from utils import aug_process
from utils import aug_path


default_wd = os.getcwd()
pathname = os.path.join(default_wd, 'auxprogs/homgenemapping')
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
examplesdir = os.path.join('../../', 'examples')
htmldir = os.path.join(default_wd, 'output_html')
tmpdir = os.path.join(pathname, 'tmp')
testdir = os.path.join(pathname, 'test_files')
bindir = os.path.join('../../../../../', 'bin/')
hgmbin = bindir+ 'homGeneMapping'
load2sqlbin = '../../bin/load2sqlitedb'

gtffilenames_with_hints = os.path.join(
    os.path.join(default_wd, testdir), 'gtffilenames.tbl')
gtffilenames_without_hints = os.path.join(
    tmpdir, 'gtffilenames_without_hints.tbl')
sqlitedb = os.path.join(tmpdir, 'homGeneMapping_hints.db')


def clean(force_tmp_dir=True, force_html_dir=True, force_result_dir=True):
    """Remove empty directories or if forced"""
    aug_path.rmtree_if_exists(resultdir, force_result_dir)
    aug_path.rmtree_if_exists(tmpdir, force_tmp_dir)
    aug_path.rmtree_if_exists(htmldir, force_html_dir)


class TestHomGeneMapping(unittest.TestCase):

    opt_compare = False
    opt_html = False

    @classmethod
    def setUpClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=True, force_result_dir=True)
        aug_path.mkdir_if_not_exists(resultdir)
        aug_path.mkdir_if_not_exists(tmpdir)

        # create file with gtf filenames without hint file references
        aug_process.execute(None,
                            'cat ' + gtffilenames_with_hints + ' | sed "s/gtf\\t.*/gtf/g"',
                            gtffilenames_without_hints)

        # create sqllite database
        out_load2sqlitedb = os.path.join(tmpdir, 'out_load2sqlitedb.txt')
        if os.path.exists(sqlitedb):
            os.remove(sqlitedb)

        aug_process.execute(None,
                            load2sqlbin +' --species=hg19 --dbaccess=' + sqlitedb + ' ' + examplesdir + '/cgp/human.fa',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            load2sqlbin + ' --species=mm9 --dbaccess=' +sqlitedb +' '+examplesdir+'/cgp/mouse.fa',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            load2sqlbin +' --noIdx --species=hg19 --dbaccess=' + sqlitedb + ' ' + testdir + '/gtfs/human.hints.gff',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            load2sqlbin + ' --noIdx --species=mm9 --dbaccess=' + sqlitedb + ' ' + testdir + '/gtfs/mouse.hints.gff',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            load2sqlbin + ' --makeIdx --dbaccess=' + sqlitedb,
                            out_load2sqlitedb)

    @classmethod
    def tearDownClass(cls):
        os.chdir(default_wd)
        clean(force_tmp_dir=False, force_html_dir=False, force_result_dir=False)

    def test_homGeneMapping_without_hints(self):
        '''
        test without hints
        gtffilenames_without_hints.tbl - is gtffilenames.tbl with hints removed
        .../bin/homGeneMapping --noDupes --gtfs=.../gtffilenames_without_hints.tbl
            --halfile=aln.hal --tmpdir=.../tmp --outdir=...//outdir
            --printHomologs=.../homologs.txt > .../homGeneMapping.out)
        '''
        os.chdir(os.path.join(default_wd, testdir))
        refdir = os.path.join(referencedir, 'without_hints')
        outputdir = os.path.join(resultdir, 'without_hints')
        aug_path.mkdir_if_not_exists(refdir)
        aug_path.mkdir_if_not_exists(outputdir)
        resfile = os.path.join(outputdir, 'homGeneMapping.out')

        args = '--noDupes --gtfs=' + gtffilenames_without_hints + ' --halfile=aln.hal \
            --tmpdir=' + outputdir +'/tmp --outdir=' + outputdir + '/outdir \
            --printHomologs=' + outputdir + '/homologs.txt'

        aug_process.execute(self, hgmbin + ' ' +args, resfile)

        # compare results
        if TestHomGeneMapping.opt_compare:
            aug_assertions.assertEqualDirectory(self, refdir, outputdir,
                                                TestHomGeneMapping.opt_html,
                                                htmldir)

    def test_homGeneMapping_with_file_hints(self):
        '''
        test with hints provided by file
        .../bin/homGeneMapping --noDupes --gtfs=gtffilenames.tbl --halfile=aln.hal
            --tmpdir=.../tmp --outdir=.../outdir
            --printHomologs=.../homologs.txt > .../homGeneMapping.out
        '''

        os.chdir(os.path.join(default_wd, testdir))
        refdir = os.path.join(referencedir, 'with_hints')
        outputdir = os.path.join(resultdir, 'with_hints')
        aug_path.mkdir_if_not_exists(refdir)
        aug_path.mkdir_if_not_exists(outputdir)
        resfile = os.path.join(outputdir, 'homGeneMapping.out')
        args = '--noDupes --gtfs=' + gtffilenames_with_hints + ' --halfile=aln.hal \
            --tmpdir=' + outputdir + '/tmp --outdir=' + outputdir + '/outdir \
            --printHomologs=' + outputdir + '/homologs.txt'

        aug_process.execute(self, hgmbin + ' ' + args,  resfile)

        # compare results
        if TestHomGeneMapping.opt_compare:
            aug_assertions.assertEqualDirectory(
                self, refdir, outputdir, TestHomGeneMapping.opt_html, htmldir)

    def test_homGeneMapping_with_sql_hints(self):
        '''
        test with hints provided by SQLite database
        .../bin/homGeneMapping --noDupes --gtfs=.../gtffilenames_without_hints.tbl
            --dbaccess=.../homGeneMapping_hints.db --halfile=aln.hal
            --tmpdir=.../tmp --outdir=.../outdir
            --printHomologs=.../homologs.txt > .../homGeneMapping.out
        '''
        os.chdir(os.path.join(default_wd, testdir))
        refdir = os.path.join(referencedir, 'with_hints')
        outputdir = os.path.join(resultdir, 'with_sql_hints')
        aug_path.mkdir_if_not_exists(refdir)
        aug_path.mkdir_if_not_exists(outputdir)
        resfile = os.path.join(outputdir, 'homGeneMapping.out')

        args = '--noDupes --gtfs=' + gtffilenames_without_hints + \
               ' --dbaccess=' + sqlitedb + ' --halfile=aln.hal \
            --tmpdir=' + outputdir + '/tmp --outdir=' + outputdir + '/outdir \
            --printHomologs=' + outputdir + '/homologs.txt'

        aug_process.execute(self, hgmbin + ' ' + args,  resfile)

        # compare results
        if TestHomGeneMapping.opt_compare:
            aug_assertions.assertEqualDirectory(
                self, refdir, outputdir, TestHomGeneMapping.opt_html, htmldir)


def test_suite():
    print('Start Test...')
    suite = unittest.TestSuite()
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_without_hints'))
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_with_file_hints'))
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_with_sql_hints'))

    return suite


def execute(compare, html, mysql):
    TestHomGeneMapping.opt_compare = compare
    TestHomGeneMapping.opt_html = html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    return result.wasSuccessful()
