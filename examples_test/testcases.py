#!/usr/bin/env python3

import itertools
import subprocess
import os
import shutil
import aug_out_filter as afilter

#TODO: tmp copy and unzip required data for test cases:
#      - test_iterative_prediction()
#      - test_iterative_prediction_with_hints()
#      - test_training_new_species(True) # with crf
#      - test_training_new_species(False)

#TODO: generate output information while tests are running

testdir = '../examples_test_output/'
augustusbin = '../bin/augustus'


def create_initial_testdir():
    if os.path.exists(testdir):
        shutil.rmtree(testdir)
    os.mkdir(testdir)


def test_utr_on():
    resfolder = testdir + test_utr_on.__name__

    os.mkdir(resfolder)
    with open(resfolder + "/aug_utr_on_tmp.gff", "w") as file:
        subprocess.call([
            augustusbin, '--species=human', '--UTR=on',
            '../examples/example.fa'
        ],
                        stdout=file,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

    # filter output file
    afilter.pred(resfolder + "/aug_utr_on_tmp.gff",
                 resfolder + "/aug_utr_on.gff")
    os.remove(resfolder + "/aug_utr_on_tmp.gff")


def test_iterative_prediction():
    resfolder = testdir + test_iterative_prediction.__name__
    species_list = ['nasonia', 'zebrafish', 'tomato']
    proc_list = []

    # run augustus several times with different parameter sets
    os.mkdir(resfolder)
    for species in species_list:
        with open(resfolder + '/aug.' + species + '.1-1M_tmp.gff',
                  "w") as file:
            proc_list.append(
                subprocess.Popen([
                    augustusbin, '--species=' + species, 'data/chr2L.sm.fa',
                    '--softmasking=on', '--predictionEnd=1000000'
                ],
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True))

    for p in proc_list:
        p.wait()

    for p in proc_list:
        error = p.stderr.read()
        if (error):
            print(error)

    # filter output
    for species in species_list:
        source = resfolder + '/aug.' + species + '.1-1M_tmp.gff'
        target = resfolder + '/aug.' + species + '.1-1M.gff'
        afilter.pred(source, target)
        os.remove(source)


def test_iterative_prediction_with_hints():
    resfolder = testdir + test_iterative_prediction_with_hints.__name__
    proc_list = []

    os.mkdir(resfolder)
    for i in range(0, 3):
        with open(resfolder + '/aug.nasonia.hints.' + str(i) + '_tmp.gff',
                  "w") as file:
            proc_list.append(
                subprocess.Popen([
                    augustusbin, '--species=nasonia', 'data/chr2L.sm.fa',
                    '--softmasking=on',
                    '--predictionStart=' + str(i * 2000000),
                    '--predictionEnd=' + str((i + 1) * 2000000 + 50000),
                    '--hintsfile=data/hints.gff',
                    '--extrinsicCfgFile=extrinsic.M.RM.E.W.cfg'
                ],
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True))

    for p in proc_list:
        p.wait()

    for p in proc_list:
        error = p.stderr.read()
        if (error):
            print(error)

    # filter output
    for i in range(0, 3):
        source = resfolder + '/aug.nasonia.hints.' + str(i) + '_tmp.gff'
        target = resfolder + '/aug.nasonia.hints.' + str(i) + '.gff'
        afilter.pred(source, target)
        os.remove(source)


def test_training_new_species(crf):
    speciesname = 'test_aug_dev_species'
    resfolder = (testdir + test_training_new_species.__name__ +
                 '_crf') if crf else (testdir +
                                      test_training_new_species.__name__)
    os.mkdir(resfolder)

    # call script to initialzie new species
    p = subprocess.Popen(
        ['perl', '../scripts/new_species.pl', '--species=' + speciesname],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    error = p.stderr.read()
    if (error):
        print(error)

    # training
    p = subprocess.Popen([
        '../bin/etraining', '../docs/tutorial2015/results/genes.gb.train', '--species=' + speciesname
    ],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # test
    with open(resfolder + "/test_tmp.out", "w") as file:
        args = [augustusbin, '../docs/tutorial2015/results/genes.gb.test', '--species=' + speciesname]
        if (crf):
            args.append('--CRF=on')
            args.append('--CRF_N=2')
            args.append('--UTR=off')
        p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # filter output file
    afilter.eval(resfolder + "/test_tmp.out", resfolder + "/test.out")
    os.remove(resfolder + "/test_tmp.out")

    # move new species to result folder
    shutil.move('../config/species/' + speciesname, resfolder)


def test_ab_initio_prediction():
    resfolder = testdir + test_ab_initio_prediction.__name__

    os.mkdir(resfolder)
    with open(resfolder + "/augustus_tmp.gff", "w") as file:
        args = [
            augustusbin, '../examples/autoAug/genome.fa',
            '--species=caenorhabditis'
        ]
        p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # filter output file
    afilter.pred(resfolder + "/augustus_tmp.gff", resfolder + "/augustus.gff")
    os.remove(resfolder + "/augustus_tmp.gff")


def test_format_and_error_out():
    resfolder = testdir + test_format_and_error_out.__name__

    os.mkdir(resfolder)
    args = [
        augustusbin, '../examples/autoAug/genome.fa',
        '--species=caenorhabditis', '--gff3=on',
        '--outfile=' + resfolder + '/augustus_tmp.gff3',
        '--errfile=' + resfolder + '/augustus.err'
    ]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # filter output file
    afilter.pred(resfolder + "/augustus_tmp.gff3",
                 resfolder + "/augustus.gff3")
    os.remove(resfolder + "/augustus_tmp.gff3")


def test_alternatives_from_sampling():
    resfolder = testdir + test_alternatives_from_sampling.__name__

    os.mkdir(resfolder)
    with open(resfolder + "/augustus_tmp.gff", "w") as file:
        args = [
            augustusbin, '../examples/autoAug/genome.fa',
            '--species=caenorhabditis', '--alternatives-from-sampling=on',
            '--minexonintronprob=0.08', '--minmeanexonintronprob=0.4',
            '--maxtracks=3'
        ]
        p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # filter output file
    afilter.pred(resfolder + "/augustus_tmp.gff", resfolder + "/augustus.gff")
    os.remove(resfolder + "/augustus_tmp.gff")


def test_cgp():
    os.chdir('../examples/cgp')
    resfolder = '../' + testdir + test_cgp.__name__
    os.mkdir(resfolder)

    with open(resfolder + "/output.txt", "w") as file:
        args = [
            '../' + augustusbin,
            '--species=human',
            '--speciesfilenames=genomes.tbl',
            '--treefile=tree.nwk',
            '--alnfile=aln.maf',
            '--alternatives-from-evidence=0',  # removes warning
            '--/CompPred/outdir=' + resfolder
        ]
        p = subprocess.Popen(args, stdout=file, stderr=subprocess.PIPE)
    p.wait()
    error = p.stderr.read()
    if (error):
        print(error)

    # filter output files
    for file in os.listdir(resfolder):
        filename = os.fsdecode(file)
        if filename.endswith(".gff"):
            afilter.cgp(
                resfolder + "/" + filename,
                resfolder + "/" + filename.replace(".gff", ".filtered.gff"))
            os.remove(resfolder + "/" + filename)

    # set working directory back to base test directory
    os.chdir('../../examples_test')


def test_cgp_sqlite(max_sub_p, resfolder, *args):
    os.mkdir(resfolder)
    proc_list = []

    # create groups according to the number of maxumim subprocesses
    grouped_args = [iter(args)] * max_sub_p

    # parallel execution of the commands of each group
    for arg_list in itertools.zip_longest(*grouped_args):
        proc_list = []
        for cmd, filename in filter(None, arg_list):
            with open(filename, "w") as file:
                proc_list.append(
                    subprocess.Popen(cmd, stdout=file, stderr=subprocess.PIPE))
        for p in proc_list:
            p.wait()
        for p in proc_list:
            error = p.stderr.read()
            if (error):
                print(error)

    # filter output prediction files
    for subdir, dirs, files in os.walk(resfolder):
        for file in files:
            filename = os.fsdecode(file)
            if filename.endswith(".gff"):
                afilter.cgp(
                    subdir + "/" + filename,
                    subdir + "/" + filename.replace(".gff", ".filtered.gff"))
                os.remove(subdir + "/" + filename)

    # set working directory back to base test directory
    os.chdir('../../../../examples_test/')


def test_cgp_denovo(max_sub_p):
    os.chdir('../docs/tutorial-cgp/results/mafs')
    resfolder = '../../../' + testdir + test_cgp_denovo.__name__
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
            resfolder + "/aug-" + str(idx) + ".out"
        ])

    test_cgp_sqlite(max_sub_p, resfolder, *args)


def test_cgp_rna_hint(max_sub_p):
    os.chdir('../docs/tutorial-cgp/results/mafs')
    resfolder = '../../../' + testdir + test_cgp_rna_hint.__name__
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
            resfolder + "/aug-" + str(idx) + ".out"
        ])

    test_cgp_sqlite(max_sub_p, resfolder, *args)


if __name__ == '__main__':
    create_initial_testdir()
    test_utr_on()
    #    test_iterative_prediction()
    #    test_iterative_prediction_with_hints()
    test_training_new_species(True) # with crf
    test_training_new_species(False)
    test_ab_initio_prediction()
    test_format_and_error_out()
    test_alternatives_from_sampling()
    test_cgp()
    test_cgp_denovo(4)
    test_cgp_rna_hint(4)
