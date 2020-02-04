#!/usr/bin/env python3

import unittest
import itertools
import json
import subprocess
import os
import shutil
import gzip
import mysql.connector
import aug_out_filter as afilter


#TODO: generate output information while tests are running
#TODO: generate optional script arguments for mysql case and compare mode (argparse)
#TODO: use conditional mysql import
#TODO: use test suite for different configs (mysql)

testdir = '../examples_test_output/'
augustusbin = '../bin/augustus'

class TestAugustus(unittest.TestCase):
    dbname = ''
    dbhost = ''
    dbuser = ''
    dbpasswd = ''
    cpuno = 2


    @classmethod
    def read_config(cls):
        with open('testconfig.json', 'r') as file:
            config = json.load(file)

        cls.dbname = config['dbname']
        cls.dbhost = config['dbhost']
        cls.dbuser = config['dbuser']
        cls.dbpasswd = config['dbpasswd']
        cls.cpuno = config['cpuno']

    
    @classmethod
    def create_initial_testdir(cls):
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        os.mkdir(testdir)


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

        cmd_list = [
            ['../bin/load2sqlitedb', '--species=hg19', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/human.fa'],
            ['../bin/load2sqlitedb', '--species=mm9', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/mouse.fa'],
            ['../bin/load2sqlitedb', '--species=bosTau4', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/cow.fa'],
            ['../bin/load2sqlitedb', '--species=galGal3', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/chicken.fa'],
            ['../bin/load2sqlitedb', '--noIdx', '--species=hg19', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/human.hints.gff'],
            ['../bin/load2sqlitedb', '--noIdx', '--species=mm9', '--dbaccess=data/tmp/vertebrates.db', '--clean', '../examples/cgp/mouse.hints.gff'],
            ['../bin/load2sqlitedb', '--makeIdx', '--dbaccess=data/tmp/vertebrates.db', '--clean']]

        print('Creating SQLite database for testing purposes:' + '\n')

        for cmd in cmd_list:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            p.wait()   
            error = p.stderr.read()
            output = p.stdout.read()
            p.stdout.close()
            p.stderr.close()
            print(error)
            print(output)


    def init_mysql_db(dbname, host, user, passwd):
        cmd_list = [
            ['../bin/load2db', '--species=hg19', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/human.fa'],  
            ['../bin/load2db', '--species=mm9', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/mouse.fa'],  
            ['../bin/load2db', '--species=bosTau4', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/cow.fa'],  
            ['../bin/load2db', '--species=galGal3', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/chicken.fa'],  
            ['../bin/load2db', '--species=hg19', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/human.hints.gff'],  
            ['../bin/load2db', '--species=mm9', '--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd, '../examples/cgp/mouse.hints.gff']] 

        print('Inserting data into MySQL database for testing purposes:' + '\n')

        for cmd in cmd_list:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            p.wait()   
            error = p.stderr.read()
            output = p.stdout.read()
            p.stdout.close()
            p.stderr.close()            
            print(error)
            print(output)


    @classmethod
    def cleanup(cls):
        # remove generated SQLite database
        if os.path.isfile('data/tmp/vertebrates.db'):
            os.remove('data/tmp/vertebrates.db')

        # remove copied/unzipped files
        if os.path.isfile('data/tmp/chr2L.sm.fa'):
            os.remove('data/tmp/chr2L.sm.fa')


    def cleanup_mysqldb(dbname, host, user, passwd):
        mysqldb = mysql.connector.connect(
            host=host,
            user=user,
            passwd=passwd,
            database=dbname
        )

        print('Clean up MySQL database.' + '\n')
        augcursor = mysqldb.cursor()
        augcursor.execute('DROP TABLE IF EXISTS genomes;')    
        augcursor.execute('DROP TABLE IF EXISTS speciesnames;')    
        augcursor.execute('DROP TABLE IF EXISTS seqnames;')    
        augcursor.execute('DROP TABLE IF EXISTS hints;')    
        augcursor.execute('DROP TABLE IF EXISTS featuretypes;')    


    @classmethod
    def setUpClass(cls):
        cls.create_initial_testdir()
        cls.init_test_data()
        cls.init_sqlite_db()


    @classmethod   
    def tearDownClass(cls):
         cls.cleanup();   


    def test_utr_on(self):
        resfolder = testdir + self.test_utr_on.__name__

        os.mkdir(resfolder)
        with open(resfolder + '/aug_utr_on_tmp.gff', 'w') as file:
            p = subprocess.Popen([
                augustusbin, '--species=human', '--UTR=on',
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
        self.assertTrue(os.path.isfile((resfolder + '/aug_utr_on_tmp.gff')), 'Output file was not created as expected!')

        # filter output file
        afilter.pred(resfolder + '/aug_utr_on_tmp.gff',
                    resfolder + '/aug_utr_on.gff')
        os.remove(resfolder + '/aug_utr_on_tmp.gff')


    def test_iterative_prediction(self):
        resfolder = testdir + self.test_iterative_prediction.__name__
        species_list = ['nasonia', 'zebrafish', 'tomato']
        proc_list = []

        # run augustus several times with different parameter sets
        os.mkdir(resfolder)
        for species in species_list:
            with open(resfolder + '/aug.' + species + '.1-1M_tmp.gff',
                    'w') as file:
                proc_list.append(
                    subprocess.Popen([
                        augustusbin, '--species=' + species, 'data/tmp/chr2L.sm.fa',
                        '--softmasking=on', '--predictionEnd=1000000'
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
            source = resfolder + '/aug.' + species + '.1-1M_tmp.gff'
            self.assertTrue(os.path.isfile(source), 'Output file was not created as expected!')
            target = resfolder + '/aug.' + species + '.1-1M.gff'
            afilter.pred(source, target)
            os.remove(source)


    def test_iterative_prediction_with_hints(self):
        resfolder = testdir + self.test_iterative_prediction_with_hints.__name__
        proc_list = []
        if not os.path.isfile('data/tmp/chr2L.sm.fa'):
            init_test_data()

        os.mkdir(resfolder)
        for i in range(0, 3):
            with open(resfolder + '/aug.nasonia.hints.' + str(i) + '_tmp.gff',
                    'w') as file:
                proc_list.append(
                    subprocess.Popen([
                        augustusbin, '--species=nasonia', 'data/tmp/chr2L.sm.fa',
                        '--softmasking=on',
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
            source = resfolder + '/aug.nasonia.hints.' + str(i) + '_tmp.gff'
            self.assertTrue(os.path.isfile(source), 'Output file was not created as expected!')
            target = resfolder + '/aug.nasonia.hints.' + str(i) + '.gff'
            afilter.pred(source, target)
            os.remove(source)


    def test_training_new_species(self):
        self.training_new_species(False)


    def test_training_new_species_crf(self):
        self.training_new_species(True)


    def training_new_species(self, crf):
        speciesname = 'test_aug_dev_species'
        resfolder = (testdir + self.test_training_new_species_crf.__name__) if crf else (testdir +
                                        self.test_training_new_species.__name__)
        os.mkdir(resfolder)

        # call script to initialzie new species
        p = subprocess.Popen(
            ['perl', '../scripts/new_species.pl', '--species=' + speciesname, '--AUGUSTUS_CONFIG_PATH=../config'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        stdout = p.stdout.read()
        p.stdout.close()
        p.stderr.close()
        #print(stdout)
        print(error)

        # training
        p = subprocess.Popen([
            '../bin/etraining', '../docs/tutorial2015/results/genes.gb.train', '--species=' + speciesname
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
        with open(resfolder + '/test_tmp.out', 'w') as file:
            args = [augustusbin,
                    '../docs/tutorial2015/results/genes.gb.test',
                    '--species=' + speciesname,
                    '--AUGUSTUS_CONFIG_PATH=../config']
            if (crf):
                args.append('--CRF=on')
                args.append('--CRF_N=2')
                args.append('--UTR=off')
            p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + '/test_tmp.out'), 'Output file was not created as expected!')
        afilter.eval(resfolder + '/test_tmp.out', resfolder + '/test.out')
        os.remove(resfolder + '/test_tmp.out')

        # move new species to result folder
        shutil.move('../config/species/' + speciesname, resfolder)


    def test_ab_initio_prediction(self):
        resfolder = testdir + self.test_ab_initio_prediction.__name__

        os.mkdir(resfolder)
        with open(resfolder + '/augustus_tmp.gff', 'w') as file:
            args = [
                augustusbin, '../examples/autoAug/genome.fa',
                '--species=caenorhabditis'
            ]
            p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + '/augustus_tmp.gff'), 'Output file was not created as expected!')
        afilter.pred(resfolder + '/augustus_tmp.gff', resfolder + '/augustus.gff')
        os.remove(resfolder + '/augustus_tmp.gff')


    def test_format_and_error_out(self):
        resfolder = testdir + self.test_format_and_error_out.__name__

        os.mkdir(resfolder)
        args = [
            augustusbin, '../examples/autoAug/genome.fa',
            '--species=caenorhabditis', '--gff3=on',
            '--outfile=' + resfolder + '/augustus_tmp.gff3',
            '--errfile=' + resfolder + '/augustus.err'
        ]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stdout.close()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + '/augustus_tmp.gff3'), 'Output file was not created as expected!')
        afilter.pred(resfolder + '/augustus_tmp.gff3',
                    resfolder + '/augustus.gff3')
        os.remove(resfolder + '/augustus_tmp.gff3')


    def test_alternatives_from_sampling(self):
        resfolder = testdir + self.test_alternatives_from_sampling.__name__

        os.mkdir(resfolder)
        with open(resfolder + '/augustus_tmp.gff', 'w') as file:
            args = [
                augustusbin, '../examples/autoAug/genome.fa',
                '--species=caenorhabditis', '--alternatives-from-sampling=on',
                '--minexonintronprob=0.08', '--minmeanexonintronprob=0.4',
                '--maxtracks=3'
            ]
            p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output file
        self.assertTrue(os.path.isfile(resfolder + '/augustus_tmp.gff'), 'Output file was not created as expected!')
        afilter.pred(resfolder + '/augustus_tmp.gff', resfolder + '/augustus.gff')
        os.remove(resfolder + '/augustus_tmp.gff')


    def test_cgp(self):
        os.chdir('../examples/cgp')
        resfolder = '../' + testdir + self.test_cgp.__name__
        os.mkdir(resfolder)

        with open(resfolder + '/output.txt', 'w') as file:
            args = [
                '../' + augustusbin,
                '--species=human',
                '--speciesfilenames=genomes.tbl',
                '--treefile=tree.nwk',
                '--alnfile=aln.maf',
                '--alternatives-from-evidence=0',  # removes warning
                '--/CompPred/outdir=' + resfolder
            ]
            p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        error = p.stderr.read()
        p.stderr.close()
        self.assertEqual(error, '', error)

        # filter output files
        for file in os.listdir(resfolder):
            filename = os.fsdecode(file)
            if filename.endswith('.gff'):
                afilter.cgp(
                    resfolder + '/' + filename,
                    resfolder + '/' + filename.replace('.gff', '.filtered.gff'))
                os.remove(resfolder + '/' + filename)

        # set working directory back to base test directory
        os.chdir('../../examples_test')


    def test_cgp_sqlite(self):
        self.cgp_with_db_preparation(False, False)


    def test_cgp_sqlite_hints(self):
        self.cgp_with_db_preparation(True, False)


    def cgp_with_db_execution(self, max_sub_p, resfolder, *args):
        os.mkdir(resfolder)
        proc_list = []

        # create groups according to the number of maxumim subprocesses
        grouped_args = [iter(args)] * max_sub_p

        # parallel execution of the commands of each group
        for arg_list in itertools.zip_longest(*grouped_args):
            proc_list = []
            for cmd, filename in filter(None, arg_list):
                with open(filename, 'w') as file:
                    proc_list.append(
                        subprocess.Popen(cmd, stdout=file, stderr=subprocess.PIPE, universal_newlines=True))
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
                        subdir + '/' + filename,
                        subdir + '/' + filename.replace('.gff', '.filtered.gff'))
                    os.remove(subdir + '/' + filename)


    def cgp_with_db_preparation(self, hints, mysql):
        #TODO: use config file for mysql information 
        #dbname = kwargs.get('dbname', None)
        #host = kwargs.get('host', None)
        #user = kwargs.get('user', None)
        #passwd = kwargs.get('passwd', None)

        if mysql:
            missing_arguments = False
            if dbname is None:
                print('The database name is missing!')
                missing_arguments = True
            if host is None:
                print('The host name is missing!')
                missing_arguments = True
            if user is None:
                print('The db user name is missing!')
                missing_arguments = True
            if passwd is None:
                print('The db user passwd is missing!')
                missing_arguments = True
            
            if missing_arguments:
                print('Test case test_cgp_with_db was not executed.')
                return
            else:
                cleanup_mysqldb('aug_vertebrates', 'localhost', 'augustus', 'aug_passwd')
                init_mysql_db(dbname, host, user, passwd)

        os.chdir('../examples/cgp')

        resfolder = '../' + testdir + 'test_cgp_with_db'
        if mysql:
            resfolder += '_mysql'
        if hints:
            resfolder += '_hints'

        cmd = [
            '../' + augustusbin,
            '--species=human',
            '--speciesfilenames=genomes.tbl',
            '--treefile=tree.nwk',
            '--alnfile=aln.maf',
            '--alternatives-from-evidence=0',  # removes warning
            '--/CompPred/outdir=' + resfolder + '/pred' 
        ]

        if mysql:
            cmd.append('--dbaccess=' + dbname + ',' + host + ',' + user + ',' + passwd)
        else:
            cmd.append('--dbaccess=../../examples_test/data/tmp/vertebrates.db')

        if hints:
            cmd.append('--dbhints=true')
            cmd.append('--extrinsicCfgFile=cgp.extrinsic.cfg')

        args = [[cmd,  resfolder + '/aug.out']]

        self.cgp_with_db_execution(1, resfolder, *args)

        # set working directory back to base test directory
        os.chdir('../../examples_test')


    def test_cgp_denovo_tutorial(self):
        os.chdir('../docs/tutorial-cgp/results/mafs')
        resfolder = '../../../' + testdir + self.test_cgp_denovo_tutorial.__name__
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
                resfolder + '/aug-' + str(idx) + '.out'
            ])

        #TODO: use config for number of maximal sub processes
        self.cgp_with_db_execution(4, resfolder, *args)

        # set working directory back to base test directory
        os.chdir('../../../../examples_test/')


    def test_cgp_rna_hint_tutorial(self):
        os.chdir('../docs/tutorial-cgp/results/mafs')
        resfolder = '../../../' + testdir + self.test_cgp_rna_hint_tutorial.__name__
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
                    '--/CompPred/outdir=' + resfolder + '/pred' + str(idx)
                ],
                resfolder + '/aug-' + str(idx) + '.out'
            ])

        #TODO: use config for number of maximal sub processes
        self.cgp_with_db_execution(4, resfolder, *args)

        # set working directory back to base test directory
        os.chdir('../../../../examples_test/')


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
    suite.addTest(TestAugustus('test_cgp_sqlite'))
    suite.addTest(TestAugustus('test_cgp_sqlite_hints'))
    return suite


def small_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestAugustus('test_utr_on'))
    suite.addTest(TestAugustus('test_training_new_species'))
    suite.addTest(TestAugustus('test_ab_initio_prediction'))
    suite.addTest(TestAugustus('test_format_and_error_out'))
    suite.addTest(TestAugustus('test_alternatives_from_sampling'))
    suite.addTest(TestAugustus('test_cgp'))
    suite.addTest(TestAugustus('test_cgp_sqlite'))
    suite.addTest(TestAugustus('test_cgp_sqlite_hints'))
    return suite


if __name__ == '__main__':
    TestAugustus.read_config()
    runner = unittest.TextTestRunner(verbosity=2)
    #runner.run(default_test_suite())
    runner.run(small_test_suite())
