#!/usr/bin/env python3
import subprocess
import re
import json
import os
import shutil
import argparse
import fileinput
import filecmp
import sys
import random
from datetime import datetime
import numpy as np

# author Giovanna Migliorelli
# version beta 25.05.2020
# the code comes as a modification of original Augustus/longrunning_examples/execute_test.py script by Daniel Honsel

# todo 0 extend basic test about code correctness
# done 1 randomly pick chunks to build the test set avoiding any bias
# done 2 merge different chunks (join genes)
# done 3 parallelize
# todo 4 replace call to shell script with pybedtools

parser = argparse.ArgumentParser(description='Execute Augustus long running test cases.')
parser.add_argument('-p', '--predict', action='store_true',
                    help='to run original prediction.')
parser.add_argument('-e', '--prepare', action='store_true',
                    help='to build a new test set from scratch.')
parser.add_argument('-o', '--port', action='store_true',
                    help='to build a new test set from scratch.')
parser.add_argument('-r', '--run', action='store_true',
                    help='to run prediction using minimal data set.')
parser.add_argument('-v', '--eval', action='store_true',
                    help='to evaluate accuracy (respect to the last prediction obtained by launching the script with --run option).')
parser.add_argument('-c', '--chunks',
                    #required=True,
                    nargs='+',
                    help='a list of one or more positive integers indicating the chunk/s to be processed (refer to documentation for a list of chunks over hg38.chr1).') 
parser.add_argument('-t', '--test', action='store_true',
                    help='to run a basic test to assess the correctness in the creation of the minimal data set.')                   
parser.add_argument('-a', '--rand',
                    help='to pick a random subset of non overlapping chunks containing at least 300 genes.')  
parser.add_argument('-g', '--augustusDir',
                    help='path to comparative Augustus executable.')  
parser.add_argument('-l', '--evalDir',
                    help='path to Eval script.')  
parser.add_argument('-w', '--workingDir',
                    help='path to data set used in testing (link).')  
args = parser.parse_args()

# if not already existing, create dir to collect results for the current chunk
def make_dirs(paths_shared, paths, chunks):
    for chunk in chunks:
        if os.path.exists(paths[chunk]['result_dir']) == False:
            os.makedirs(paths[chunk]['result_dir'])
    if os.path.exists(paths_shared['accuracy']) == False:
            os.makedirs(paths_shared['accuracy'])
    if os.path.exists(paths_shared['joingenes_out_dir']) == False:
            os.makedirs(paths_shared['joingenes_out_dir'])

# extract minimal FASTAs on the base of gene ranges (BEDs) and prepare sqlitedb
def make_sqlitedb(paths_shared, paths, chunk):
    cleanup_db(paths, chunk)

    print('Creating minimal FASTAs and SQLite database for current chunk...', paths_shared['bedtool_bin'])

    # extract minimal FASTAs
    # todo replace the following with pybedtools
    subprocess.call([paths_shared['working_dir'] + 'extractFASTA.sh', paths_shared['bedtool_bin'], paths_shared['fasta_dir'], paths_shared['working_dir'], str(chunk), paths[chunk]['sqlitedb_dir']])
        
    for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
        # replace characters ':' and '-' with '_'
        for line in fileinput.input(paths[chunk]['sqlitedb_dir'] + species + '.MINIMAL.fasta', inplace=True):
            print(line.replace(':', "_").replace('-', "_"))

        # add sepcies to SQLite db
        cmd = [paths_shared['augustus_dir'] + 'bin/load2sqlitedb', '--dbaccess=' + paths[chunk]['sqlitedb_test_file'], '--species=' + species, paths[chunk]['sqlitedb_dir'] + species + '.MINIMAL.fasta']
        execute(cmd, paths[chunk]['result_dir'] + 'out.createMinimalFASTA', mode ='a+')

def port_sqlitedb(paths_shared, paths, chunk):
    cleanup_db(paths, chunk, False)
    
    for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
        # add sepcies to SQLite db
        cmd = [paths_shared['augustus_dir'] + 'bin/load2sqlitedb', '--dbaccess=' + paths[chunk]['sqlitedb_test_file'], '--species=' + species, paths[chunk]['sqlitedb_dir'] + species + '.MINIMAL.fasta']
        execute(cmd, paths[chunk]['result_dir'] + 'out.createMinimalFASTA', mode ='a+')

# create new genomes.tbl for this chunk
def make_genometbl_chunk(paths, chunk):
    cleanup_tbl_chunk(paths, chunk)
    print('Creating new genome tbl for current chunk...', paths[chunk]['tbl_test_file'])

    with open(paths[chunk]['tbl_test_file'], 'w') as f:
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            print(species, paths[chunk]['sqlitedb_dir'] + species + '.MINIMAL.fasta', sep='\t', file=f)

# create genomes.tbl for full size genomes
def make_genometbl(paths_shared):
    cleanup_tbl(paths_shared)
    print('Creating new genome tbl for full size genomes')

    with open(paths_shared['fasta_dir'], 'w') as f:
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            print(species, paths_shared['fasta_dir'] + species + '.fasta', sep='\t', file=f)
    
def cleanup_db(paths, chunk, removeFASTA = True):
    print('Cleaning up', paths[chunk]['sqlitedb_dir'])
    
    if removeFASTA:
        if os.path.exists(paths[chunk]['sqlitedb_dir']):
            shutil.rmtree(paths[chunk]['sqlitedb_dir'])
        os.makedirs(paths[chunk]['sqlitedb_dir'])

    print('Cleaning up', paths[chunk]['sqlitedb_test_file'])
    if os.path.exists(paths[chunk]['sqlitedb_test_file']):
        os.remove(paths[chunk]['sqlitedb_test_file'])

def cleanup_tbl_chunk(paths, chunk):
    print('Cleaning up', paths[chunk]['tbl_test_file'])
    if os.path.exists(paths[chunk]['tbl_test_file']):
        os.remove(paths[chunk]['tbl_test_file'])

def cleanup_tbl(paths_shared):
    print('Cleaning up', paths_shared['tbl_file'])
    if os.path.exists(paths_shared['tbl_file']):
        os.remove(paths_shared['tbl_file'])

# runs original prediction over full length genomes (used only for a basic test about correctness of the code)
def run_prediction(paths, chunk):
    print('Runnning prediction on chunk', chunk, '...')

    cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], '--alnfile=' + paths[chunk]['maf_file'],
    '--speciesfilenames=' + paths_shared['tbl_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + paths_shared['sqlitedb_file'],
    '--optCfgFile=../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
    '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction']
        
    execute(cmd, paths[chunk]['result_dir'] + 'out.runPrediction')

def run_prediction_parallel(paths_shared, paths, chunks):

    proc_list = []

    for chunk in chunks:
        print('Runnning prediction on chunk', chunk, '...')

        cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], 
        '--alnfile=' + paths[chunk]['maf_file'],
        '--speciesfilenames=' + paths_shared['tbl_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + paths_shared['sqlitedb_file'],
        '--optCfgFile=../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
        '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction']

        filename = paths[chunk]['result_dir'] + 'out.runPrediction'

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

def prepare_test(paths_shared, paths, chunks):
    for chunk in chunks:
        print('Preparing test for chunk', chunk, '...')
        cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], '--alnfile=' + paths[chunk]['maf_file'], '--dbaccess=' + paths_shared['sqlitedb_file'], 
        '--speciesfilenames=' + paths_shared['tbl_file'], '--softmasking=1', '--alternatives-from-evidence=0',
        '--optCfgFile=../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
        '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'prepare', '--/Testing/testMode=prepare']

        execute(cmd, paths[chunk]['result_dir'] + 'out.prepareTest')

        make_sqlitedb(paths_shared, paths, chunk)
        make_genometbl_chunk(paths, chunk)

def port_test(paths_shared, paths, chunks):
    for chunk in chunks:
        print('Porting test for chunk', chunk, '...')
        # port_sqlitedb(paths_shared, paths, chunk)
        make_genometbl_chunk(paths, chunk)
    #make_genometbl(paths_shared)

def run_test(paths, chunk):
    print('Runnning prediction on chunk', chunk, 'using the minimal data set...')
    cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], 
    '--alnfile=' + paths[chunk]['maf_file'],
    '--speciesfilenames=' + paths[chunk]['tbl_test_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + paths[chunk]['sqlitedb_test_file'],
    '--optCfgFile=../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
    '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'run', '--/Testing/testMode=run',
    '--/Testing/workingDir=' + paths_shared['working_dir'] + 'names']
    
    
    execute(cmd, paths[chunk]['result_dir'] + 'out.runTest')

# parallel execution : acknowldgement Daniel Honsel (revisited code from test_case.py)
def run_test_parallel(paths_shared, paths, chunks):

    proc_list = []

    for chunk in chunks:
        print('Runnning prediction on chunk', chunk, 'using the minimal data set...')

        cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], 
        '--alnfile=' + paths[chunk]['maf_file'],
        '--speciesfilenames=' + paths[chunk]['tbl_test_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + paths[chunk]['sqlitedb_test_file'],
        '--optCfgFile=../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
        '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'run', '--/Testing/testMode=run',
        '--/Testing/workingDir=' + paths_shared['working_dir']]

        output = paths[chunk]['result_dir'] + 'out.runTest'

        with open(output, 'w') as file:
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

def execute(cmd, output, mode='w'):
    with open(output, mode) as file:
        p = subprocess.Popen(cmd,
                             stdout=file,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)

    p.wait()
    error = p.stderr.read()
    p.stderr.close()
    if error:
        print(error)

# minimal test : prediction obtained working with minimal FASTAs is compared against original prediction, they should be identical for the test to succeed
def test_test(paths_shared, paths, chunks):
    for chunk in chunks:
        print('Runnning tests on chunk', chunk, '...')
        if os.path.exists(paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction'):
            for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
                pathToGFF = paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction/' + species + '.cgp.gff'
                if os.path.exists(pathToGFF) == False:
                    goahead = False
                    print('Cannot find', pathToGFF, 'no test will be run over the code because original prediction is missing.') 
                    return
        else:
            print('Cannot find', paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction', 'no test will be run over the code because original prediction is missing.') 
            return

        if os.path.exists(paths_shared['working_dir'] + 'out' + str(chunk) + 'run'):
            for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
                pathToGFF = paths_shared['working_dir'] + 'out' + str(chunk) + 'run/' + species + '.cgp.gff'
                if os.path.exists(pathToGFF) == False:
                    goahead = False
                    print('Cannot find', pathToGFF, 'no test will be run over the code because any new prediction over minimal data set is missing.') 
                    return
        else:
            print('Cannot find', paths_shared['working_dir'] + 'out' + str(chunk) + 'run', 'no test will be run over the code because any new prediction over minimal data set is missing.') 
            return
        
        print('Assessing correctness of code...')
        correct = True
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            pathToGFF_1 = paths_shared['working_dir'] + 'out' + str(chunk) + 'prediction/' + species + '.cgp.gff'
            pathToGFF_2 = paths_shared['working_dir'] + 'out' + str(chunk) + 'run/' + species + '.cgp.gff'
            if filecmp.cmp(pathToGFF_1, pathToGFF_2) == False:
                print('\tsomething went wrong: ', pathToGFF_1, 'differs from', pathToGFF_2)
                correct = False
        if correct:
            print('\ttest on code succeeded')

# currently not used (returns accuracy for single chunks) - not parallelized
def run_evaluate(paths, chunk):
    print('Runnning evaluation on chunk', chunk, '...')
    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
    execute(cmd, paths[chunk]['result_dir'] + 'out.eval')
    find_values(paths[chunk]['result_dir'] + 'out.eval')

# currently not used (returns accuracy for single chunks) - parallelized
def run_evaluate_parallel(paths, chunks):
    proc_list = []

    for chunk in chunks:
        print('Runnning evaluation on chunk', chunk, '...')
        cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
        execute(cmd, paths[chunk]['result_dir'] + 'out.eval')
        find_values(paths[chunk]['result_dir'] + 'out.eval')

        filename = paths[chunk]['result_dir'] + 'out.eval'

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


# currently in use (returns accuracy after merging the contributes from all chunks) - parallelized
def run_evaluate_global(paths_shared, paths, chunks):
    print('Runnning evaluation on chunks', chunks, '...')
    
    with open(paths_shared['joingenes_out_dir'] + 'jgGFFs', 'w') as f:
        for chunk in chunks:
            print(paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff', '1', sep='\t', file=f)  
        
    cmd = [paths_shared['joingenes_bin'], '-f' + paths_shared['joingenes_out_dir'] + 'jgGFFs', '-o' + paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['joingenes_out_dir'] + 'out.joingenes')

    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['accuracy'] + 'out.eval')
    find_values(paths_shared['accuracy'] + 'out.eval')
     

def find_values(file):
    d = init_dict()

    # Summary Stats
    with open(file, 'r') as input:
        input_str = input.read()

    d['gene_sensitivity'] = value_from_summary_stats(input_str,
                                                     'Gene Sensitivity')
    d['gene_specificity'] = value_from_summary_stats(input_str,
                                                     'Gene Specificity')
    d['transcript_sensitivity'] = value_from_summary_stats(
        input_str, 'Transcript Sensitivity')
    d['transcript_specificity'] = value_from_summary_stats(
        input_str, 'Transcript Specificity')
    d['exon_sensitivity'] = value_from_summary_stats(input_str,
                                                     'Exon Sensitivity')
    d['exon_specificity'] = value_from_summary_stats(input_str,
                                                     'Exon Specificity')
    d['nucleotide_sensitivity'] = value_from_summary_stats(
        input_str, 'Nucleotide Sensitivity')
    d['nucleotide_specificity'] = value_from_summary_stats(
        input_str, 'Nucleotide Specificity')

    # general stats
    with open(file, 'r') as input:
        general = re.search('\*\*General Stats\*\*.*\*\*Detailed Stats\*\*',
                            input.read(), re.DOTALL).group()

    b_gc = find_block(general, 'Gene', 'Transcript')
    d['gene_count'] = find_value(b_gc, 'Count')
    d['tx_count'] = find_value(b_gc, 'Total Transcripts')
    d['txs_per_gene'] = find_value(b_gc, 'Transcripts Per')

    b_tx = find_block(general, 'Transcript', 'Exon')
    b_tx_all = find_block(b_tx, 'All', 'Complete')
    d['avg_tx_length'] = find_value(b_tx_all, 'Average Length')
    d['median_tx_length'] = find_value(b_tx_all, 'Median Length')
    d['avg_coding_length'] = find_value(b_tx_all, 'Average Coding Length')
    d['median_coding_length'] = find_value(b_tx_all, 'Median Coding Length')
    d['avg_exons_per_tx'] = find_value(b_tx_all, 'Ave Exons Per')

    b_single_ex = find_block(b_tx, 'Single Exon', 'Exon')
    d['single_exon_count'] = find_value(b_single_ex, 'Count')

    b_ex = find_block(general, 'Exon', 'Nuc')
    b_ex_in = find_block(b_ex, 'Intron', 'InframeOptional')
    d['avg_intron_length'] = find_value(b_single_ex, 'Average Length')
    d['median_intron_length'] = find_value(b_single_ex, 'Median Length')

    # compute f1 scores
    d['gene_fscore'] = hmean(d['gene_sensitivity'], d['gene_specificity'])
    d['transcript_fscore'] = hmean(d['transcript_sensitivity'],
                                   d['transcript_specificity'])
    d['exon_fscore'] = hmean(d['exon_sensitivity'], d['exon_specificity'])
    d['nucleotide_fscore'] = hmean(d['nucleotide_sensitivity'],
                                   d['nucleotide_specificity'])

    # save results as JSON file
    with open(paths_shared['accuracy'] + 'eval.json', 'w') as file:
        json.dump(d, file, indent=4)

def hmean(v1, v2):
    if v1 > 0 and v2 > 0:
        return 2 * (v1 * v2 / (v1 + v2))
    else:
        return 0.0

def find_block(input, start, stop):
    match = re.search(start + '[ ]*\n\t.*' + stop, input, re.DOTALL)
    return match.group(0)

def find_value(input, name):
    match = re.search(
        '\t\t' + name + '[ ]*\t[0-9]+\.[0-9]+[ ]*\t[0-9]+\.[0-9]+', input)
    return float(match.group(0).rsplit('\t', 1)[1])

def value_from_summary_stats(input, name):
    match = re.search(name + '[ ]*\t[0-9]+\.[0-9]+', input).group(0)
    return float(match.split('\t')[1])

def init_dict():
    dictionary = {
        'gene_sensitivity': 0.0,
        'gene_specificity': 0.0,
        'gene_fscore': 0.0,
        'transcript_sensitivity': 0.0,
        'transcript_specificity': 0.0,
        'transcript_fscore': 0.0,
        'exon_sensitivity': 0.0,
        'exon_specificity': 0.0,
        'exon_fscore': 0.0,
        'nucleotide_sensitivity': 0.0,
        'nucleotide_specificity': 0.0,
        'nucleotide_fscore': 0.0,
        'gene_count': 0.0,
        'tx_count': 0.0,
        'txs_per_gene': 0.0,
        'avg_tx_length': 0.0,
        'median_tx_length': 0.0,
        'avg_coding_length': 0.0,
        'median_coding_length': 0.0,
        'avg_exons_per_tx': 0.0,
        'single_exon_count': 0.0,
        'avg_intron_length': 0.0,
        'median_intron_length': 0.0
    }
    return dictionary

def execute(cmd, output, mode='w'):
    with open(output, mode) as file:
        p = subprocess.Popen(cmd,
                             stdout=file,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)

    p.wait()
    error = p.stderr.read()
    p.stderr.close()
    if error:
        print(error)

def init_paths_shared(augustusDir, workingDir, evalDir):
    paths_shared = {
    'eval_dir' : evalDir,           #'/home/giovanna/Desktop/Alignment/eval-2.2.8/',    # path to eval
    'augustus_dir' : augustusDir,   # '../',                                    # path to augustus binaries
    'working_dir' : workingDir,     # '../examples/cgp12way/',                        # path to working directory (it contains tree, genome tbl, SQLite db and there results will be written)
    'anno_file' : workingDir + 'ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.FILTERED.gtf',  # path to annotation for hg38.chr1

    # the following three directories are required only if a new test set is to be built
    'maf_dir' : '../examples/cgp12way/MAF/',                                    # path to MAFs 
    'fasta_dir' : '/home/giovanna/Desktop/Alignment/DATA_UCSCSOFT/',            # path to original FASTA files for genomes of interest
    'bedtools_dir' : '/home/giovanna/Desktop/Alignment/bedtoolsBinaries/'       # path to bedtools (need to extract minimal FASTA from original FASTA, on the base of BED format)
    }

    paths_shared.update({'log' : paths_shared['working_dir'] + 'LOG/'})                             # todo : move here all outXresult
    paths_shared.update({'accuracy' : paths_shared['working_dir'] + 'ACCURACY/'})                   # results returned by eval
    # wrong paths_shared.update({'joingenes_bin' : paths_shared['augustus_dir'] + '/auxprogs/joingenes'})   # path to joingenes exec dir
    paths_shared.update({'joingenes_out_dir' : paths_shared['working_dir'] + 'JOINGENES/'})         # output from joingenes
    paths_shared.update({'sqlitedb_file' : paths_shared['working_dir'] + 'SQLITE/12way.db'})        # path to SQlite for full length genomes (only if a new data set is to be built)
    paths_shared.update({'tbl_file' : paths_shared['working_dir'] + 'GENOMETBL/genomes.tbl'})       # path to genomes.tbl for length genomes (only if a new data set is to be built)
    paths_shared.update({'tree_file' : paths_shared['working_dir'] + 'TREE/ucsc12way.nwk'})         # path to tree
    paths_shared.update({'augustus_bin' : paths_shared['augustus_dir'] + 'bin/augustus'})           # path to augustus exec
    paths_shared.update({'joingenes_bin' : paths_shared['augustus_dir'] + 'auxprogs/joingenes/joingenes'})  # path to joingenes exec
    paths_shared.update({'bedtool_bin' : paths_shared['bedtools_dir'] + 'bedtools'})                # path to bedtools (only if a new data set is to be built)
    paths_shared.update({'eval_bin' : paths_shared['eval_dir'] + 'evaluate_gtf.pl'})                # path to eval

    return paths_shared  

def init_paths(chunks):
    paths = {}
    for chunk in chunks:
        dictionary = {}

        # the following paths automatically reflect changes in paths_shared
        dictionary.update({'sqlitedb_dir' : paths_shared['working_dir'] + 'minimalFasta' + str(chunk) + '/'})               # will contain minimal fasta after their extraction
        dictionary.update({'sqlitedb_test_file' : paths_shared['working_dir'] + 'SQLITE/12wayTEST_' + str(chunk) + '.db'})
        dictionary.update({'tbl_test_file' : paths_shared['working_dir'] + 'GENOMETBL/genomesTEST_' + str(chunk) + '.tbl'})
        dictionary.update({'maf_file' : paths_shared['maf_dir'] + "chr1_chunk_" + str(chunk) + '.maf'})
        dictionary.update({'result_dir' : paths_shared['working_dir'] + 'out' + str(chunk) + 'result/'})

        paths[chunk] = dictionary

    return paths

# given a list of chunks (no header admitted, tabseparated) and, for each of them, the number of genes it contains, a random subset is picked according to the following:
# no chunks in the data set overlap (rule out contiguous ones)
# the sum of genes in greater than 299
def randomize_dataset(filename):
    random.seed(datetime.now())

    data = np.genfromtxt(filename)

    index = np.random.choice(range(0,len(data)), len(data), replace=False)

    dataset = []
    numgenes = 0
    for i in index:
        if data[i][1] > 0:
            res = [x for x in dataset if abs(data[i][0] - x)<2]
            if len(res) == 0:
                dataset.append(int(data[i][0]))
                numgenes = numgenes + data[i][1]
                if numgenes>=300:
                    break

    print('Sampled dataset contains', int(numgenes), 'genes from chunks:', [x for x in dataset])

def expandDir(path):
    tmp = path
    if len(path)>0 and path[len(path)-1] != '/':
        tmp += '/'
    return tmp

if __name__ == '__main__':
    if args.rand:
        randomize_dataset(args.rand)
        sys.exit()

    if args.chunks is None:
        print('No chunks specified, please make use of --chunks to pass a non empty list of positive integers...')
        sys.exit()

    chunks = [int(x) for x in list(dict.fromkeys(args.chunks))]
    chunks = [x for x in chunks if x>0 and x<126]          # range valid for chr1 chunk size 2.5 Mb, chunk overlap 0.5 Mb

    if len(chunks) == 0:
        print('No valid chunks specified...')
        sys.exit()

    if args.augustusDir is None:
        print('Path to comparative augustus executable required, please make use of --augustusDir to pass the path...')
        sys.exit()
    augustusDir = str(expandDir(args.augustusDir))  

    if args.workingDir is None:
        print('Path to data set used in testing required, please make use of --workingDir to pass the path...')
        sys.exit()
    workingDir = str(expandDir(args.workingDir))

    evalDir = ''
    if args.eval:
        if args.eval is None:
            print('Path to Eval script required, please make use of --evalDir to pass the path...')
            sys.exit()
        evalDir = str(expandDir(args.evalDir))

    paths_shared = init_paths_shared(augustusDir, workingDir, evalDir)
    paths = init_paths(chunks)

    make_dirs(paths_shared, paths, chunks)

    if args.predict:
        run_prediction_parallel(paths_shared, paths, chunks)
    if args.prepare:
        prepare_test(paths_shared, paths, chunks)
    if args.run:
        port_test(paths_shared, paths, chunks)
        run_test_parallel(paths_shared, paths, chunks)
    if args.eval:
        run_evaluate_global(paths_shared, paths, chunks)
    if args.test:
        test_test(paths_shared, paths, chunks)
    if args.port:
        port_test(paths_shared, paths, chunks)
