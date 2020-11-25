#!/usr/bin/env python3

import unittest
from pathlib import Path
import os
import sys
import shutil

if not (os.path.exists('tests')):
    sys.exit('Wrong working directory!\nTest files must be accessible in ./tests folder.')
sys.path.insert(0, '.') # make modul tests.utils accessible

from tests.utils import aug_assertions
from tests.utils import aug_process
from tests.utils import aug_argparse
from tests.utils import aug_path


args = aug_argparse.getDefaultArgParser().parse_args()

default_wd = os.getcwd()
pathname = os.path.join(default_wd, os.path.dirname(sys.argv[0]))
referencedir = os.path.join(pathname, 'expected_results')
resultdir = os.path.join(pathname, 'result_files')
examplesdir = os.path.join(default_wd, 'examples')
htmldir = os.path.join(pathname, 'output_html')
tmpdir = os.path.join(pathname, 'tmp')
testdir = os.path.join(pathname, 'test_files')
bindir = os.path.join(default_wd, 'bin/')
hgmbin = f'{bindir}homGeneMapping'
load2sqlbin = f'{bindir}load2sqlitedb'

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

def check_working_dir(clean):
    errstr = ''
    if not (os.path.exists(testdir)):
        errstr += 'Wrong example directory!' + '\n'
        errstr += f'The example files must be accessible in this path: "{testdir}"!'
    if not clean and not (os.path.exists(hgmbin)):
        errstr += 'Wrong working directory!' + '\n'
        errstr += f'The augustus binaries must be accessible in this path: "{bindir}"!'
    if (errstr):
        sys.exit(errstr)


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
                            f'cat {gtffilenames_with_hints} | sed "s/gtf\\t.*/gtf/g"',
                            gtffilenames_without_hints)

        # create sqllite database
        out_load2sqlitedb = os.path.join(tmpdir, 'out_load2sqlitedb.txt')
        if os.path.exists(sqlitedb):
            os.remove(sqlitedb)

        aug_process.execute(None,
                            f'{load2sqlbin} --species=hg19 --dbaccess={sqlitedb} {examplesdir}/cgp/human.fa',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            f'{load2sqlbin} --species=mm9 --dbaccess={sqlitedb} {examplesdir}/cgp/mouse.fa',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            f'{load2sqlbin} --noIdx --species=hg19 --dbaccess={sqlitedb} {testdir}/gtfs/human.hints.gff',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            f'{load2sqlbin} --noIdx --species=mm9 --dbaccess={sqlitedb} {testdir}/gtfs/mouse.hints.gff',
                            out_load2sqlitedb)
        aug_process.execute(None,
                            f'{load2sqlbin} --makeIdx --dbaccess={sqlitedb}',
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

        args = f'--noDupes --gtfs={gtffilenames_without_hints} --halfile=aln.hal \
            --tmpdir={outputdir}/tmp --outdir={outputdir}/outdir \
            --printHomologs={outputdir}/homologs.txt'

        aug_process.execute(self, f'{hgmbin} {args}', resfile)

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
        args = f'--noDupes --gtfs={gtffilenames_with_hints} --halfile=aln.hal \
            --tmpdir={outputdir}/tmp --outdir={outputdir}/outdir \
            --printHomologs={outputdir}/homologs.txt'

        aug_process.execute(self, f'{hgmbin} {args}',  resfile)

        # compare results
        if TestHomGeneMapping.opt_compare:
            aug_assertions.assertEqualDirectory(self, refdir, outputdir, TestHomGeneMapping.opt_html, htmldir)

    def test_homGeneMapping_with_sql_hints(self):
        '''
        test with hints provided by SQLite databaseX
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

        args = f'--noDupes --gtfs={gtffilenames_without_hints} \
            --dbaccess={sqlitedb} --halfile=aln.hal \
            --tmpdir={outputdir}/tmp --outdir={outputdir}/outdir \
            --printHomologs={outputdir}/homologs.txt'

        aug_process.execute(self, f'{hgmbin} {args}',  resfile)

        # compare results
        if TestHomGeneMapping.opt_compare:
            aug_assertions.assertEqualDirectory(self, refdir, outputdir, TestHomGeneMapping.opt_html, htmldir)


def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_without_hints'))
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_with_file_hints'))
    suite.addTest(TestHomGeneMapping('test_homGeneMapping_with_sql_hints'))

    return suite


if __name__ == '__main__':
    check_working_dir(args.clean)

    if args.clean:
        clean()
        sys.exit()

    TestHomGeneMapping.opt_compare = args.compare
    TestHomGeneMapping.opt_html = args.html

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite())

    if result.wasSuccessful():
        sys.exit()
    else:
        sys.exit(1)
