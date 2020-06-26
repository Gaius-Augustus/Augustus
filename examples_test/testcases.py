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

parser = argparse.ArgumentParser(description='Execute Augustus test cases.')
parser.add_argument('--mysql',
                    action='store_true',
                    help='Execute cgp test cases using a MySQL database.')
parser.add_argument('--compare',
                    action='store_true',
                    help='Compare generated results with reference results.')
parser.add_argument('--html',
                    action='store_true',
                    help='Save diff results in html file.')
parser.add_argument(
    '--set_default_wd',
    action='store_true',
    help=
    'Set the working directory to examples_test (if called from AUGUSTUS root).'
)
args = parser.parse_args()

# only import mysql connector if testcases using mysql should be executed
# MySQL Connector must be installed in this case
if args.mysql:
    import mysql.connector

resultdir = '../examples_test_results/'
refdir = '../examples_results/'
htmldir = 'output_html/'
augustusbin = '../bin/augustus'
default_wd = os.getcwd()

def create_initial_resultdir():
    if os.path.exists(htmldir):
        shutil.rmtree(htmldir)

    if os.path.exists(resultdir):
        shutil.rmtree(resultdir)
    os.mkdir(resultdir)


def check_working_dir():
    wd = os.getcwd()
    if not (wd.endswith('Augustus/examples_test')):
        errstr = 'Wrong working directory!' + '\n'
        errstr += 'This script must be called from "Augustus/examples_test"!'
        sys.exit(errstr)


class TestAugustus(unittest.TestCase):
    dbname = ''
    dbhost = ''
    dbuser = ''
    dbpasswd = ''
    cpuno = 2

    opt_compare = False
    opt_html = False

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
        if not os.path.exists('data/tmp'):
            os.mkdir('data/tmp')

        inputfile = 'data/tmp/chr2L.sm.fa.gz'
        shutil.copyfile('../docs/tutorial2015/data/chr2L.sm.fa.gz', inputfile)

        with gzip.open(inputfile, 'rb') as f_in:
            with open('data/tmp/chr2L.sm.fa', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(inputfile)

    @classmethod
    def init_sqlite_db(cls):
        if not os.path.exists('data/tmp'):
            os.mkdir('data/tmp')

        cmd_list = [[
            '../bin/load2sqlitedb', '--species=hg19',
            '--dbaccess=data/tmp/vertebrates.db', '--clean',
            '../examples/cgp/human.fa'
        ],
                    [
                        '../bin/load2sqlitedb', '--species=mm9',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean',
                        '../examples/cgp/mouse.fa'
                    ],
                    [
                        '../bin/load2sqlitedb', '--species=bosTau4',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean',
                        '../examples/cgp/cow.fa'
                    ],
                    [
                        '../bin/load2sqlitedb', '--species=galGal3',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean',
                        '../examples/cgp/chicken.fa'
                    ],
                    [
                        '../bin/load2sqlitedb', '--noIdx', '--species=hg19',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean',
                        '../examples/cgp/human.hints.gff'
                    ],
                    [
                        '../bin/load2sqlitedb', '--noIdx', '--species=mm9',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean',
                        '../examples/cgp/mouse.hints.gff'
                    ],
                    [
                        '../bin/load2sqlitedb', '--makeIdx',
                        '--dbaccess=data/tmp/vertebrates.db', '--clean'
                    ]]

        print('Creating SQLite database for cgp test cases...')

        for cmd in cmd_list:
            p = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
            p.wait()
            error = p.stderr.read()
            output = p.stdout.read()
            p.stdout.close()
            p.stderr.close()
            if error:
                print(error)
            #print(output)

    @classmethod
    def init_mysql_db(cls):
        cmd_list = [[
            '../bin/load2db', '--species=hg19', '--dbaccess=' + cls.dbname +
            ',' + cls.dbhost + ',' + cls.dbuser + ',' + cls.dbpasswd,
            '../examples/cgp/human.fa'
        ],
                    [
                        '../bin/load2db', '--species=mm9',
                        '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
                        cls.dbuser + ',' + cls.dbpasswd,
                        '../examples/cgp/mouse.fa'
                    ],
                    [
                        '../bin/load2db', '--species=bosTau4',
                        '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
                        cls.dbuser + ',' + cls.dbpasswd,
                        '../examples/cgp/cow.fa'
                    ],
                    [
                        '../bin/load2db', '--species=galGal3',
                        '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
                        cls.dbuser + ',' + cls.dbpasswd,
                        '../examples/cgp/chicken.fa'
                    ],
                    [
                        '../bin/load2db', '--species=hg19',
                        '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
                        cls.dbuser + ',' + cls.dbpasswd,
                        '../examples/cgp/human.hints.gff'
                    ],
                    [
                        '../bin/load2db', '--species=mm9',
                        '--dbaccess=' + cls.dbname + ',' + cls.dbhost + ',' +
                        cls.dbuser + ',' + cls.dbpasswd,
                        '../examples/cgp/mouse.hints.gff'
                    ]]

        print('  -' +
              'Inserting data into MySQL database for testing purposes...')

        for cmd in cmd_list:
            p = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
            p.wait()
            error = p.stderr.read()
            output = p.stdout.read()
            p.stdout.close()
            p.stderr.close()
            if error:
                print(error)
            #print(output)

    @classmethod
    def cleanup(cls):
        # remove generated SQLite database
        if os.path.isfile('data/tmp/vertebrates.db'):
            os.remove('data/tmp/vertebrates.db')

        # remove copied/unzipped files
        if os.path.isfile('data/tmp/chr2L.sm.fa'):
            os.remove('data/tmp/chr2L.sm.fa')

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
        cls.init_test_data()
        cls.init_sqlite_db()

    @classmethod
    def tearDownClass(cls):
        cls.cleanup()

    def test_utr_on(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_utr_on.__name__ + '/'

        os.mkdir(resfolder)
        with open(resfolder + 'aug_utr_on_tmp.gff', 'w') as file:
            p = subprocess.Popen([
                augustusbin, '--species=human', '--UTR=on', '--softmasking=0',
                '../examples/example.fa'
            ],
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        rc = p.wait()
        error = p.stderr.read()
        p.stderr.close()

        self.assertEqual(error, '', error)
        self.assertEqual(rc, 0, 'Returncode not 0!')
        self.assertTrue(os.path.isfile((resfolder + 'aug_utr_on_tmp.gff')),
                        'Output file was not created as expected!')

        # filter output file
        afilter.pred(resfolder + 'aug_utr_on_tmp.gff',
                     resfolder + 'aug_utr_on.gff')
        os.remove(resfolder + 'aug_utr_on_tmp.gff')

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(refdir + self.test_utr_on.__name__ +
                                       '/',
                                       resfolder,
                                       html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_iterative_prediction(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_iterative_prediction.__name__ + '/'
        species_list = ['nasonia', 'zebrafish', 'tomato']
        proc_list = []

        # run augustus several times with different parameter sets
        os.mkdir(resfolder)
        for species in species_list:
            with open(resfolder + 'aug.' + species + '.1-1M_tmp.gff',
                      'w') as file:
                proc_list.append(
                    subprocess.Popen([
                        augustusbin, '--species=' + species,
                        'data/tmp/chr2L.sm.fa', '--softmasking=on',
                        '--predictionEnd=1000000'
                    ],
                                     stdout=file,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True))
        for p in proc_list:
            p.wait()

        for p in proc_list:
            error = p.stderr.read()
            p.stderr.close()
            self.assertEqual(error, '', error)

        # filter output
        for species in species_list:
            source = resfolder + 'aug.' + species + '.1-1M_tmp.gff'
            self.assertTrue(os.path.isfile(source),
                            'Output file was not created as expected!')
            target = resfolder + 'aug.' + species + '.1-1M.gff'
            afilter.pred(source, target)
            os.remove(source)

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(
                refdir + self.test_iterative_prediction.__name__ + '/',
                resfolder,
                html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_iterative_prediction_with_hints(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_iterative_prediction_with_hints.__name__ + '/'
        proc_list = []
        if not os.path.isfile('data/tmp/chr2L.sm.fa'):
            init_test_data()

        os.mkdir(resfolder)
        for i in range(0, 3):
            with open(resfolder + 'aug.nasonia.hints.' + str(i) + '_tmp.gff',
                      'w') as file:
                proc_list.append(
                    subprocess.Popen([
                        augustusbin, '--species=nasonia',
                        'data/tmp/chr2L.sm.fa', '--softmasking=on',
                        '--predictionStart=' + str(i * 2000000),
                        '--predictionEnd=' + str((i + 1) * 2000000 + 50000),
                        '--hintsfile=../docs/tutorial2015/results/hints.gff',
                        '--extrinsicCfgFile=extrinsic.M.RM.E.W.cfg'
                    ],
                                     stdout=file,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True))

        for p in proc_list:
            p.wait()

        for p in proc_list:
            error = p.stderr.read()
            p.stderr.close()
            self.assertEqual(error, '', error)

        # filter output
        for i in range(0, 3):
            source = resfolder + 'aug.nasonia.hints.' + str(i) + '_tmp.gff'
            self.assertTrue(os.path.isfile(source),
                            'Output file was not created as expected!')
            target = resfolder + 'aug.nasonia.hints.' + str(i) + '.gff'
            afilter.pred(source, target)
            os.remove(source)

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(
                refdir + self.test_iterative_prediction_with_hints.__name__ +
                '/',
                resfolder,
                html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_training_new_species(self):
        os.chdir(default_wd)
        self.training_new_species(False)

    def test_training_new_species_crf(self):
        os.chdir(default_wd)
        self.training_new_species(True)

    def training_new_species(self, crf):
        speciesname = 'test_aug_dev_species'
        resfolder = (
            resultdir +
            self.test_training_new_species_crf.__name__) if crf else (
                resultdir + self.test_training_new_species.__name__) + '/'
        reffolder = (
            refdir + self.test_training_new_species_crf.__name__) if crf else (
                refdir + self.test_training_new_species.__name__) + '/'
        os.mkdir(resfolder)
        
        # call script to initialzie new species
        p = subprocess.Popen([
            'perl', '../scripts/new_species.pl', '--species=' + speciesname,
            '--AUGUSTUS_CONFIG_PATH=../config'
        ],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        stdout = p.stdout.read()
        p.stdout.close()
        p.stderr.close()
        #print(stdout)
        if error:
            print(error)

        # training
        p = subprocess.Popen([
            '../bin/etraining', '../docs/tutorial2015/results/genes.gb.train',
            '--species=' + speciesname
        ],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        stdout = p.stdout.read()
        p.stdout.close()
        p.stderr.close()
        #print(stdout)
        self.assertEqual(error, '', error)

        # test
        with open(resfolder + 'test_tmp.out', 'w') as file:
            cmd = [
                augustusbin, '../docs/tutorial2015/results/genes.gb.test',
                '--species=' + speciesname, '--softmasking=0', '--AUGUSTUS_CONFIG_PATH=../config'
            ]
            if (crf):
                cmd.append('--CRF=on')
                cmd.append('--CRF_N=2')
                cmd.append('--UTR=off')
            p = subprocess.Popen(cmd,
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + 'test_tmp.out'),
                        'Output file was not created as expected!')
        afilter.eval(resfolder + 'test_tmp.out', resfolder + 'test.out')
        os.remove(resfolder + 'test_tmp.out')

        # move new species to result folder
        shutil.move('../config/species/' + speciesname, resfolder)

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(reffolder,
                                       resfolder,
                                       html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_ab_initio_prediction(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_ab_initio_prediction.__name__ + '/'

        os.mkdir(resfolder)
        with open(resfolder + '/augustus_tmp.gff', 'w') as file:
            cmd = [
                augustusbin, '../examples/autoAug/genome.fa', '--softmasking=1',
                '--species=caenorhabditis'
            ]
            p = subprocess.Popen(cmd,
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + 'augustus_tmp.gff'),
                        'Output file was not created as expected!')
        afilter.pred(resfolder + 'augustus_tmp.gff',
                     resfolder + '/augustus.gff')
        os.remove(resfolder + 'augustus_tmp.gff')

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(
                refdir + self.test_ab_initio_prediction.__name__ + '/',
                resfolder,
                html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_format_and_error_out(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_format_and_error_out.__name__ + '/'

        os.mkdir(resfolder)
        cmd = [
            augustusbin, '../examples/autoAug/genome.fa',
            '--species=caenorhabditis', '--gff3=on', '--softmasking=1',
            '--outfile=' + resfolder + 'augustus_tmp.gff3',
            '--errfile=' + resfolder + 'augustus.err'
        ]
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stdout.close()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + 'augustus_tmp.gff3'),
                        'Output file was not created as expected!')
        afilter.pred(resfolder + 'augustus_tmp.gff3',
                     resfolder + 'augustus.gff3')
        os.remove(resfolder + 'augustus_tmp.gff3')

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(
                refdir + self.test_format_and_error_out.__name__ + '/',
                resfolder,
                html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_alternatives_from_sampling(self):
        os.chdir(default_wd)
        resfolder = resultdir + self.test_alternatives_from_sampling.__name__ + '/'

        os.mkdir(resfolder)
        with open(resfolder + 'augustus_tmp.gff', 'w') as file:
            cmd = [
                augustusbin, '../examples/autoAug/genome.fa',
                '--species=caenorhabditis', '--alternatives-from-sampling=on',
                '--minexonintronprob=0.08', '--minmeanexonintronprob=0.4',
                '--maxtracks=3'
            ]
            p = subprocess.Popen(cmd,
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + 'augustus_tmp.gff'),
                        'Output file was not created as expected!')
        afilter.pred(resfolder + 'augustus_tmp.gff',
                     resfolder + 'augustus.gff')
        os.remove(resfolder + 'augustus_tmp.gff')

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(
                refdir + self.test_alternatives_from_sampling.__name__ + '/',
                resfolder,
                html=TestAugustus.opt_html)
            self.assertEqual(diff, '', diff)

    def test_cgp(self):
        os.chdir(default_wd)
        os.chdir('../examples/cgp')
        resfolder = '../' + resultdir + self.test_cgp.__name__ + '/'
        reffolder = '../' + refdir + self.test_cgp.__name__ + '/'
        os.mkdir(resfolder)

        with open(resfolder + 'output_tmp.txt', 'w') as file:
            cmd = [
                '../' + augustusbin,
                '--species=human',
                '--speciesfilenames=genomes.tbl',
                '--treefile=tree.nwk',
                '--alnfile=aln.maf',
                '--softmasking=0',
                '--alternatives-from-evidence=0',  # removes warning
                '--/CompPred/outdir=' + resfolder
            ]
            p = subprocess.Popen(cmd,
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output files
        for file in os.listdir(resfolder):
            filename = os.fsdecode(file)
            if filename.endswith('.gff'):
                afilter.cgp(
                    resfolder + filename, resfolder + '/' +
                    filename.replace('.gff', '.filtered.gff'))
                os.remove(resfolder + filename)
        afilter.cgp_out(resfolder + 'output_tmp.txt', resfolder + 'output.txt')
        os.remove(resfolder + 'output_tmp.txt')

        # compare results
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(reffolder,
                                       resfolder,
                                       html=TestAugustus.opt_html,
                                       outputfolder=default_wd + '/output_html/')
            self.assertEqual(diff, '', diff)


    def test_cgp_sqlite(self):
        os.chdir(default_wd)
        self.cgp_with_db_preparation(False, False)

    def test_cgp_sqlite_hints(self):
        os.chdir(default_wd)
        self.cgp_with_db_preparation(True, False)

    def test_cgp_mysql(self):
        os.chdir(default_wd)
        self.cgp_with_db_preparation(False, True)

    def test_cgp_mysql_hints(self):
        os.chdir(default_wd)
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
        if TestAugustus.opt_compare:
            diff = comp.compare_folder(reffolder,
                                       resfolder,
                                       html=TestAugustus.opt_html,
                                       outputfolder=default_wd + '/output_html/')
            self.assertEqual(diff, '', diff)

    def cgp_with_db_preparation(self, hints, mysql):
        if mysql:
            missing_arguments = False
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
                print('Test case test_cgp_with_db was not executed.')
                return 1
            else:
                TestAugustus.cleanup_mysqldb()
                TestAugustus.init_mysql_db()

        os.chdir('../examples/cgp')

        resfolder = '../' + resultdir + 'test_cgp_with_db'
        reffolder = '../' + refdir + 'test_cgp_with_db'
        if mysql:
            resfolder += '_mysql'
            reffolder += '_mysql'
        if hints:
            resfolder += '_hints'
            reffolder += '_hints'
        resfolder += '/'
        reffolder += '/'

        cmd = [
            '../' + augustusbin,
            '--species=human',
            '--speciesfilenames=genomes.tbl',
            '--treefile=tree.nwk',
            '--alnfile=aln.maf',
            '--softmasking=0',
            '--alternatives-from-evidence=0',  # removes warning
            '--/CompPred/outdir=' + resfolder + 'pred'
        ]

        if mysql:
            cmd.append('--dbaccess=' + TestAugustus.dbname + ',' +
                       TestAugustus.dbhost + ',' + TestAugustus.dbuser + ',' +
                       TestAugustus.dbpasswd)
        else:
            cmd.append(
                '--dbaccess=../../examples_test/data/tmp/vertebrates.db')

        if hints:
            cmd.append('--dbhints=true')
            cmd.append('--extrinsicCfgFile=cgp.extrinsic.cfg')

        args = [[cmd, resfolder + 'aug_tmp.out']]

        self.cgp_with_db_execution(resfolder, reffolder, *args)


    def test_cgp_denovo_tutorial(self):
        os.chdir(default_wd)
        os.chdir('../docs/tutorial-cgp/results/mafs')
        resfolder = '../../../' + resultdir + self.test_cgp_denovo_tutorial.__name__ + '/'
        reffolder = '../../../' + refdir + self.test_cgp_denovo_tutorial.__name__ + '/'
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
                    '--/CompPred/outdir=' + resfolder + 'pred' + str(idx)
                ],
                resfolder + 'aug-' + str(idx) + '_tmp.out'
            ])

        self.cgp_with_db_execution(resfolder, reffolder, *args)


    def test_cgp_rna_hint_tutorial(self):
        os.chdir(default_wd)
        os.chdir('../docs/tutorial-cgp/results/mafs')
        resfolder = '../../../' + resultdir + self.test_cgp_rna_hint_tutorial.__name__ + '/'
        reffolder = '../../../' + refdir + self.test_cgp_rna_hint_tutorial.__name__ + '/'
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
                    '--dbhints=1',
                    '--UTR=1',
                    '--allow_hinted_splicesites=atac',
                    '--extrinsicCfgFile=../extrinsic-rnaseq.cfg',
                    '--/CompPred/outdir=' + resfolder + 'pred' + str(idx)
                ],
                resfolder + 'aug-' + str(idx) + '_tmp.out'
            ])

        self.cgp_with_db_execution(resfolder, reffolder, *args)


def default_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_utr_on'))
    suite.addTest(TestAugustus('test_iterative_prediction'))
    suite.addTest(TestAugustus('test_iterative_prediction_with_hints'))
    suite.addTest(TestAugustus('test_training_new_species'))
    suite.addTest(TestAugustus('test_training_new_species_crf'))
    suite.addTest(TestAugustus('test_ab_initio_prediction'))
    suite.addTest(TestAugustus('test_format_and_error_out'))
    suite.addTest(TestAugustus('test_alternatives_from_sampling'))
    suite.addTest(TestAugustus('test_cgp'))
    os.chdir(default_wd)
    suite.addTest(TestAugustus('test_cgp_sqlite'))
    os.chdir(default_wd)
    suite.addTest(TestAugustus('test_cgp_sqlite_hints'))
    os.chdir(default_wd)
    return suite


def small_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_utr_on'))
    suite.addTest(TestAugustus('test_training_new_species'))
    suite.addTest(TestAugustus('test_ab_initio_prediction'))
    suite.addTest(TestAugustus('test_format_and_error_out'))
    #suite.addTest(TestAugustus('test_alternatives_from_sampling'))
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
    if args.set_default_wd:
        os.chdir('examples_test/')

    check_working_dir()
    default_wd = os.getcwd()

    create_initial_resultdir()
    TestAugustus.opt_compare = args.compare
    TestAugustus.opt_html = args.html
    runner = unittest.TextTestRunner(verbosity=2)
    #print_tc_header('default test suite')
    #runner.run(default_test_suite())
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
