#!/usr/bin/env python3

import argparse
import unittest
import itertools
import json
import subprocess
import os
import sys
import shutil
import gzip
import aug_out_filter as afilter
import aug_comparator as comp

# This script executes AUGUSTUS test cases based on the examples
# folder and compares the current results with reference results
# if the option --compare is set. It is expected that both results
# are identical for a successful test.
# This script must be called from "tests/examples_test"!
# Python version 3.6 or higher is required for execution.

parser = argparse.ArgumentParser(description='Execute Augustus test cases.')
parser.add_argument('--mysql',
                    action='store_true',
                    help='cgp test cases are also executed with a MySQL database.')
parser.add_argument('--compare',
                    action='store_true',
                    help='Compare generated results with reference results.')
parser.add_argument('--html',
                    action='store_true',
                    help='Save diff results in html file.')
parser.add_argument('--clean',
                    action='store_true',
                    help='Remove all files created during the tests. If this option is set, no tests are executed.')
args = parser.parse_args()

# only import mysql connector if testcases using mysql should be executed
# MySQL Connector must be installed in this case
if args.mysql:
    import mysql.connector

resultdir = 'results/'
refdir = 'expected_results/'
htmldir = 'output_html/'
tmpdir = 'data/tmp/'
exampledir = '../../examples/'
bindir = '../../bin/'
augustusbin = f'{bindir}augustus'
datadir =  exampledir + 'chr2L/'
default_wd = os.getcwd()


def create_initial_resultdir():
    clean(False)
    os.mkdir(resultdir)


def clean(withtmpdir=True):
    print('Removing generated test files...')
    if os.path.exists(htmldir):
        shutil.rmtree(htmldir)

    if os.path.exists(resultdir):
        shutil.rmtree(resultdir)

    if withtmpdir and os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


def check_working_dir(clean):
    wd = os.getcwd()
    if not (wd.endswith('tests/examples')):
        errstr = 'Wrong working directory!' + '\n'
        errstr += 'This script must be called from "tests/examples"!'
        sys.exit(errstr)
    if not clean and not (os.path.exists(augustusbin)):
        errstr = 'Missing augustus binaries!' + '\n'
        errstr += f'The augustus binaries must be accessible in this path: "{bindir}"!'
        sys.exit(errstr)


class TestAugustus(unittest.TestCase):
    dbname = None
    dbhost = None
    dbuser = None
    dbpasswd = None
    cpuno = 2

    opt_compare = False
    opt_html = False
    opt_mysql = False

    @classmethod
    def read_config(cls):
        with open('testconfig.json', 'r') as file:
            config = json.load(file)

        cls.dbname = config['dbname']
        cls.dbhost = config['dbhost']
        cls.dbuser = config['dbuser']
        cls.dbpasswd = config['dbpasswd']
        cls.cpuno = int(config['cpuno'])

    @classmethod
    def init_test_data(cls):
        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)

        inputfile = os.path.join(tmpdir, 'chr2L.sm.fa.gz')
        testfile = os.path.join(tmpdir, 'chr2L.sm.fa')
        shutil.copyfile(os.path.join(datadir, 'chr2L.sm.fa.gz'), inputfile)
#            '../../docs/tutorial2015/data/chr2L.sm.fa.gz', inputfile)

        with gzip.open(inputfile, 'rb') as f_in:
            with open(testfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(inputfile)

    @classmethod
    def init_sqlite_db(cls):
        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)

        cmd_list = [[
            f'{bindir}load2sqlitedb', '--species=hg19',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/human.fa'
        ],
            [
            f'{bindir}load2sqlitedb', '--species=mm9',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/mouse.fa'
        ],
            [
            f'{bindir}load2sqlitedb', '--species=bosTau4',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/cow.fa'
        ],
            [
            f'{bindir}load2sqlitedb', '--species=galGal3',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/chicken.fa'
        ],
            [
            f'{bindir}load2sqlitedb', '--noIdx', '--species=hg19',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/human.hints.gff'
        ],
            [
            f'{bindir}load2sqlitedb', '--noIdx', '--species=mm9',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean',
            f'{exampledir}cgp/mouse.hints.gff'
        ],
            [
            f'{bindir}load2sqlitedb', '--makeIdx',
            f'--dbaccess={tmpdir}vertebrates.db', '--clean'
        ]]

        print('Creating SQLite database for cgp test cases...')

        cls.init_db(cmd_list)

    @classmethod
    def init_mysql_db(cls):
        cmd_list = [[
            f'{bindir}load2db', '--species=hg19', '--dbaccess=' + cls.dbname +
            ',' + cls.dbhost + ',' + cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/human.fa'
        ],
            [
            f'{bindir}load2db', '--species=mm9',
            '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
            cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/mouse.fa'
        ],
            [
            f'{bindir}load2db', '--species=bosTau4',
            '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
            cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/cow.fa'
        ],
            [
            f'{bindir}load2db', '--species=galGal3',
            '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
            cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/chicken.fa'
        ],
            [
            f'{bindir}load2db', '--species=hg19',
            '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
            cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/human.hints.gff'
        ],
            [
            f'{bindir}load2db', '--species=mm9',
            '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
            cls.dbuser + ',' + cls.dbpasswd,
            f'{exampledir}cgp/mouse.hints.gff'
        ]]

        print('  -' +
              'Inserting data into MySQL database for testing purposes...')

        cls.init_db(cmd_list)

    @classmethod
    def init_db(cls, cmd_list):
        for cmd in cmd_list:
            output = TestAugustus().process(cmd)
            # print(output)

    @classmethod
    def cleanup(cls):
        os.chdir(default_wd)
        # remove generated SQLite database
        if os.path.isfile(f'{tmpdir}vertebrates.db'):
            os.remove(f'{tmpdir}vertebrates.db')

        # remove copied/unzipped files
        if os.path.isfile(f'{tmpdir}chr2L.sm.fa'):
            os.remove(f'{tmpdir}chr2L.sm.fa')

    @classmethod
    def cleanup_mysqldb(cls):
        mysqldb = mysql.connector.connect(host=cls.dbhost,
                                          user=cls.dbuser,
                                          passwd=cls.dbpasswd,
                                          database=cls.dbname)

        print('\n' + '  -' + 'Clean up MySQL database...')
        augcursor = mysqldb.cursor()
        augcursor.execute('DROP TABLE IF EXISTS genomes;')
        augcursor.execute('DROP TABLE IF EXISTS speciesnames;')
        augcursor.execute('DROP TABLE IF EXISTS seqnames;')
        augcursor.execute('DROP TABLE IF EXISTS hints;')
        augcursor.execute('DROP TABLE IF EXISTS featuretypes;')

    @classmethod
    def setUpClass(cls):
        cls.read_config()

        # check config
        missing_arguments = False
        if (cls.opt_mysql):
            if TestAugustus.dbname is None:
                print('The database name is missing!')
                missing_arguments = True
            if TestAugustus.dbhost is None:
                print('The host name is missing!')
                missing_arguments = True
            if TestAugustus.dbuser is None:
                print('The db user name is missing!')
                missing_arguments = True
            if TestAugustus.dbpasswd is None:
                print('The db user passwd is missing!')
                missing_arguments = True
        if missing_arguments:
            assert False, 'Test case using MySQL are not executed.'

        cls.init_test_data()
        cls.init_sqlite_db()
        if (cls.opt_mysql):
            cls.cleanup_mysqldb()
            cls.init_mysql_db()

    @classmethod
    def tearDownClass(cls):
        cls.cleanup()
        if (cls.opt_mysql):
            cls.cleanup_mysqldb()

    def assertEqualFolders(self, reffolder, resfolder, html=None, outputfolder=None):
        if TestAugustus.opt_compare:
            if html is None:
                html = self.opt_html
            if outputfolder is None:
                diff = comp.compare_folder(reffolder,
                                           resfolder,
                                           html=html)
            else:
                diff = comp.compare_folder(reffolder,
                                           resfolder,
                                           html=html,
                                           outputfolder=outputfolder)
            self.assertEqual(diff, '', diff)

    def get_ref_folder(self, folder_name=None, path_to_wd=None):
        if folder_name is None:
            folder_name = self._testMethodName
        if path_to_wd is None:
            return os.path.join(refdir, folder_name)
        else:
            return os.path.join(path_to_wd, refdir, folder_name)

    def get_res_folder(self, folder_name=None, path_to_wd=None):
        if folder_name is None:
            folder_name = self._testMethodName
        if path_to_wd is None:
            return os.path.join(resultdir, folder_name)
        else:
            return os.path.join(path_to_wd, resultdir, folder_name)

    def process(self, cmd_list, out=subprocess.PIPE):
        isFile = isinstance(out, str)

        output = out
        if isFile:
            output = open(out, 'w')

        p = subprocess.Popen(cmd_list,
                             stdout=output,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        rc = p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)
        self.assertEqual(rc, 0, f'Returncode not 0! Error: {error}')

        if isFile:
            self.assertTrue(os.path.isfile(out),
                            'Output file was not created as expected!')
        else:
            stdout = p.stdout.read()
            p.stdout.close()
            return stdout

        return ''

    def test_utr_on(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        testtmpfile = os.path.join(resfolder, 'aug_utr_on_tmp.gff')
        testfile = os.path.join(resfolder, 'aug_utr_on.gff')
        os.mkdir(resfolder)

        self.process([
            augustusbin, '--species=human', '--UTR=on', '--softmasking=0',
            f'{exampledir}example.fa'
        ], testtmpfile)

        # filter output file
        afilter.pred(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_iterative_prediction(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        os.mkdir(resfolder)

        species_list = ['nasonia', 'zebrafish', 'tomato']

        # run augustus several times with different parameter sets
        for species in species_list:
            testtmpfile = os.path.join(
                resfolder, 'aug.' + species + '.1-1M_tmp.gff')
            self.process([
                augustusbin, '--species=' + species,
                f'{tmpdir}chr2L.sm.fa', '--softmasking=on',
                '--predictionEnd=1000000'
            ], testtmpfile)

            # filter output
            testfile = os.path.join(resfolder, 'aug.' + species + '.1-1M.gff')
            afilter.pred(testtmpfile, testfile)
            os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_iterative_prediction_with_hints(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        os.mkdir(resfolder)

        if not os.path.isfile('data/tmp/chr2L.sm.fa'):
            TestAugustus.init_test_data()

        for i in range(0, 3):
            testtmpfile = os.path.join(
                resfolder, f'aug.nasonia.hints.{str(i)}_tmp.gff')
            self.process([
                augustusbin, '--species=nasonia',
                f'{tmpdir}chr2L.sm.fa', '--softmasking=on',
                '--predictionStart=' + str(i * 2000000),
                '--predictionEnd=' + str((i + 1) * 2000000 + 50000),
                f'--hintsfile={datadir}/hints.gff',
                '--extrinsicCfgFile=extrinsic.M.RM.E.W.cfg'
            ], testtmpfile)

            # filter output
            testfile = os.path.join(
                resfolder, f'aug.nasonia.hints.{str(i)}.gff')
            afilter.pred(testtmpfile, testfile)
            os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_training_new_species(self):
        self.training_new_species(False)

    def test_training_new_species_crf(self):
        self.training_new_species(True)

    def training_new_species(self, crf):
        os.chdir(default_wd)
        speciesname = 'test_aug_dev_species'

        # Remove test species folder.
        # Just in case the deletion fails for whatever reason.
        if os.path.exists('../../config/species/' + speciesname):
            shutil.rmtree('../../config/species/' + speciesname)

        resfolder = self.get_res_folder()
        reffolder = self.get_ref_folder()
        testtmpfile = os.path.join(resfolder, 'test_tmp.out')
        testfile = os.path.join(resfolder, 'test.out')
        os.mkdir(resfolder)

        # call script to initialize new species
        self.process([
            'perl', '../../scripts/new_species.pl', '--species=' + speciesname,
            '--AUGUSTUS_CONFIG_PATH=../../config'
        ])

        # training
        self.process([
            f'{bindir}etraining', os.path.join(datadir, 'genes.gb.train'),
            '--species=' + speciesname
        ])

        # test
        cmd = [
            augustusbin, os.path.join(datadir, 'genes.gb.test'),
            '--species=' + speciesname, '--softmasking=0',
            '--AUGUSTUS_CONFIG_PATH=../../config'
        ]
        if (crf):
            cmd.append('--CRF=on')
            cmd.append('--CRF_N=2')
            cmd.append('--UTR=off')

        self.process(cmd, testtmpfile)

        # filter output file
        afilter.eval(testtmpfile, testfile)
        os.remove(testtmpfile)

        # move new species to result folder
        shutil.move('../../config/species/' + speciesname, resfolder)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_ab_initio_prediction(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        testtmpfile = os.path.join(resfolder, 'augustus_tmp.gff')
        testfile = os.path.join(resfolder, 'augustus.gff')
        os.mkdir(resfolder)

        self.process([
            augustusbin, f'{exampledir}autoAug/genome.fa', '--softmasking=1',
            '--species=caenorhabditis'
        ], testtmpfile)

        # filter output file
        afilter.pred(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_format_and_error_out(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        testtmpfile = os.path.join(resfolder, 'augustus_tmp.gff3')
        testfile = os.path.join(resfolder, 'augustus.gff3')
        os.mkdir(resfolder)

        cmd = [
            augustusbin, f'{exampledir}autoAug/genome.fa',
            '--species=caenorhabditis', '--gff3=on', '--softmasking=1',
            '--outfile=' + testtmpfile,
            '--errfile=' + resfolder + '/augustus.err'
        ]
        self.process(cmd)

        # filter output file
        self.assertTrue(os.path.isfile(testtmpfile),
                        'Output file was not created as expected!')
        afilter.pred(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_alternatives_from_sampling(self):
        os.chdir(default_wd)
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        testtmpfile = os.path.join(resfolder, 'augustus_tmp.gff')
        testfile = os.path.join(resfolder, 'augustus.gff')
        os.mkdir(resfolder)

        cmd = [
            augustusbin, f'{exampledir}autoAug/genome.fa',
            '--species=caenorhabditis', '--alternatives-from-sampling=on',
            '--minexonintronprob=0.08', '--minmeanexonintronprob=0.4',
            '--maxtracks=3'
        ]
        self.process(cmd, testtmpfile)

        # filter output file
        afilter.pred(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)

    def test_cgp(self):
        reffolder = self.get_ref_folder(path_to_wd='../../tests/examples')
        resfolder = self.get_res_folder(path_to_wd='../../tests/examples')
        testtmpfile = os.path.join(resfolder, 'output_tmp.txt')
        testfile = os.path.join(resfolder, 'output.txt')

        os.chdir(os.path.join(default_wd, f'{exampledir}cgp'))
        os.mkdir(resfolder)

        cmd = [
            augustusbin,
            '--species=human',
            '--speciesfilenames=genomes.tbl',
            '--treefile=tree.nwk',
            '--alnfile=aln.maf',
            '--softmasking=0',
            '--alternatives-from-evidence=0',  # removes warning
            '--/CompPred/outdir=' + resfolder + '/'
        ]
        self.process(cmd, testtmpfile)

        # filter output files
        for file in os.listdir(resfolder):
            filename = os.fsdecode(file)
            if filename.endswith('.gff'):
                afilter.cgp(os.path.join(resfolder, filename),
                            os.path.join(resfolder, filename.replace('.gff', '.filtered.gff')))
                os.remove(os.path.join(resfolder, filename))
        afilter.cgp_out(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder,
                                outputfolder=default_wd + '/output_html/')

    def test_cgp_sqlite(self):
        self.cgp_with_db_preparation(False, False)

    def test_cgp_sqlite_hints(self):
        self.cgp_with_db_preparation(True, False)

    def test_cgp_mysql(self):
        self.cgp_with_db_preparation(False, True)

    def test_cgp_mysql_hints(self):
        self.cgp_with_db_preparation(True, True)

    def cgp_with_db_execution(self, resfolder, reffolder, *args):
        os.mkdir(resfolder)
        proc_list = []

        # create groups according to the configured number of cpus
        grouped_args = [iter(args)] * TestAugustus.cpuno

        # parallel execution of the commands of each group
        for arg_list in itertools.zip_longest(*grouped_args):
            proc_list = []
            for cmd, filename in filter(None, arg_list):
                with open(filename, 'w') as file:
                    proc_list.append(
                        subprocess.Popen(cmd,
                                         stdout=file,
                                         stderr=subprocess.PIPE,
                                         universal_newlines=True))
            for p in proc_list:
                p.wait()
            for p in proc_list:
                error = p.stderr.read()
                p.stderr.close()
                self.assertEqual(error, '', error)

        # filter output prediction files
        for subdir, dirs, files in os.walk(resfolder):
            for file in files:
                filename = os.fsdecode(file)
                if filename.endswith('.gff'):
                    afilter.cgp(
                        subdir + '/' + filename, subdir + '/' +
                        filename.replace('.gff', '.filtered.gff'))
                    os.remove(subdir + '/' + filename)
                elif filename.endswith('.out'):
                    afilter.cgp_out(
                        subdir + '/' + filename,
                        subdir + '/' + filename.replace('_tmp', ''))
                    os.remove(subdir + '/' + filename)

        # compare results
        self.assertEqualFolders(reffolder, resfolder,
                                outputfolder=default_wd + '/output_html/')

    def cgp_with_db_preparation(self, hints, mysql):
        os.chdir(os.path.join(default_wd, f'{exampledir}cgp'))

        testname = 'test_cgp_with_db'
        if mysql:
            testname += '_mysql'
        if hints:
            testname += '_hints'
        resfolder = self.get_res_folder(testname, '../../tests/examples')
        reffolder = self.get_ref_folder(testname, '../../tests/examples')

        cmd = [
            augustusbin,
            '--species=human',
            '--speciesfilenames=genomes.tbl',
            '--treefile=tree.nwk',
            '--alnfile=aln.maf',
            '--softmasking=0',
            '--alternatives-from-evidence=0',  # removes warning
            '--/CompPred/outdir=' + resfolder + '/pred'
        ]

        if mysql:
            cmd.append('--dbaccess=' + TestAugustus.dbname + ',' +
                       TestAugustus.dbhost + ',' + TestAugustus.dbuser + ',' +
                       TestAugustus.dbpasswd)
        else:
            cmd.append(
                '--dbaccess=../../tests/examples/data/tmp/vertebrates.db')

        if hints:
            cmd.append('--dbhints=true')
            cmd.append('--extrinsicCfgFile=cgp.extrinsic.cfg')

        args = [[cmd, resfolder + '/aug_tmp.out']]

        self.cgp_with_db_execution(resfolder, reffolder, *args)

    def test_cgp_denovo_tutorial(self):
        os.chdir(default_wd)
        os.chdir('../../docs/tutorial-cgp/results/mafs')
        resfolder = self.get_res_folder('test_cgp_with_db')
        reffolder = self.get_ref_folder('test_cgp_with_db')
        args = []

        # create command list for all alignment files
        for idx, alin in enumerate(os.listdir(os.curdir), 1):
            args.append([
                [
                    '../../../' + augustusbin,
                    '--species=human',
                    '--softmasking=1',
                    '--speciesfilenames=../../../../examples_test/data/cgp_genomes.tbl',
                    '--treefile=../../data/tree.nwk',
                    '--alnfile=' + alin.__str__(),
                    '--alternatives-from-evidence=0',  # removes warning
                    '--dbaccess=../vertebrates.db',
                    '--/CompPred/outdir=' + resfolder + '/pred' + str(idx)
                ],
                resfolder + '/aug-' + str(idx) + '_tmp.out'
            ])

        self.cgp_with_db_execution(resfolder, reffolder, *args)

    def test_cgp_rna_hint_tutorial(self):
        os.chdir(default_wd)
        os.chdir('../../docs/tutorial-cgp/results/mafs')
        reffolder = self.get_ref_folder(path_to_wd='../../../../tests/examples')
        resfolder = self.get_res_folder(path_to_wd='../../../../tests/examples')
        args = []

        # create command list for all alignment files
        for idx, alin in enumerate(os.listdir(os.curdir), 1):
            args.append([
                [
                    '../../../' + augustusbin,
                    '--species=human',
                    '--softmasking=1',
                    '--speciesfilenames=../../../../tests/examples_test/data/cgp_genomes.tbl',
                    '--treefile=../../data/tree.nwk',
                    '--alnfile=' + alin.__str__(),
                    '--alternatives-from-evidence=0',  # removes warning
                    '--dbaccess=../vertebrates.db',
                    '--dbhints=1',
                    '--UTR=1',
                    '--allow_hinted_splicesites=atac',
                    '--extrinsicCfgFile=../extrinsic-rnaseq.cfg',
                    '--/CompPred/outdir=' + resfolder + '/pred' + str(idx)
                ],
                resfolder + '/aug-' + str(idx) + '_tmp.out'
            ])

        self.cgp_with_db_execution(resfolder, reffolder, *args)

    def test_hints_MPE(self):
        reffolder = self.get_ref_folder()
        resfolder = self.get_res_folder()
        testtmpfile = os.path.join(resfolder, 'aug_hints_MPE_tmp.gff')
        testfile = os.path.join(resfolder, 'aug_hints_MPE.gff')

        os.chdir(default_wd)
        os.mkdir(resfolder)
        self.process([
            augustusbin, '--species=human', f'--hintsfile={exampledir}hints.gff',
            '--extrinsicCfgFile=../../config/extrinsic/extrinsic.MPE.cfg',
            f'{exampledir}example.fa'
        ], testtmpfile)

        # filter output file
        afilter.pred(testtmpfile, testfile)
        os.remove(testtmpfile)

        # compare results
        self.assertEqualFolders(reffolder, resfolder)


def default_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_utr_on'))
    suite.addTest(TestAugustus('test_hints_MPE'))
    suite.addTest(TestAugustus('test_iterative_prediction'))
    suite.addTest(TestAugustus('test_iterative_prediction_with_hints'))
    suite.addTest(TestAugustus('test_training_new_species'))
    suite.addTest(TestAugustus('test_training_new_species_crf'))
    suite.addTest(TestAugustus('test_ab_initio_prediction'))
    suite.addTest(TestAugustus('test_format_and_error_out'))
    suite.addTest(TestAugustus('test_alternatives_from_sampling'))
    suite.addTest(TestAugustus('test_cgp'))
    suite.addTest(TestAugustus('test_cgp_sqlite'))
    suite.addTest(TestAugustus('test_cgp_sqlite_hints'))
    return suite


def small_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_utr_on'))
    suite.addTest(TestAugustus('test_hints_MPE'))
    suite.addTest(TestAugustus('test_training_new_species'))
    suite.addTest(TestAugustus('test_ab_initio_prediction'))
    suite.addTest(TestAugustus('test_format_and_error_out'))
    # suite.addTest(TestAugustus('test_alternatives_from_sampling'))
    suite.addTest(TestAugustus('test_cgp'))
    suite.addTest(TestAugustus('test_cgp_sqlite'))
    suite.addTest(TestAugustus('test_cgp_sqlite_hints'))
    return suite


def mysql_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_cgp_mysql'))
    suite.addTest(TestAugustus('test_cgp_mysql_hints'))
    return suite


def print_tc_header(tc_name):
    print(
        '----------------------------------------------------------------------'
    )
    print('Executing ' + tc_name)
    print(
        '----------------------------------------------------------------------'
    )


if __name__ == '__main__':
    check_working_dir(args.clean)
    default_wd = os.getcwd()

    # Remove only generated test files and do not execute test
    # cases if option --clean is set.
    if args.clean:
        clean()
        sys.exit()

    create_initial_resultdir()
    TestAugustus.opt_compare = args.compare
    TestAugustus.opt_html = args.html
    TestAugustus.opt_mysql = args.mysql
    runner = unittest.TextTestRunner(verbosity=2)
    #print_tc_header('default test suite')
    #result = runner.run(default_test_suite())
    print_tc_header('small test suite')
    result = runner.run(small_test_suite())

    mysql_was_successful = True
    if args.mysql:
        os.chdir(default_wd)
        print_tc_header('MySQL test suite')
        result_mysql = runner.run(mysql_test_suite())
        mysql_was_successful = result_mysql.wasSuccessful()

    if result.wasSuccessful() and mysql_was_successful:
        sys.exit()
    else:
        sys.exit(1)
