#!/usr/bin/env python3
import subprocess
import re
import json
import os
import shutil
import argparse
import fileinput
import filecmp

# author Giovanna Migliorelli
# version beta 25.05.2020
# the code comes as a modification of original Augustus/longrunning_examples/execute_test.py script by Daniel Honsel

# todo 0 extend basic test about code correctness
# todo 1 randomly pick chunks to build the test set avoiding any bias
# todo 2 merge different chunks (join genes)
# todo 3 parallelize
# todo 4 replace call to shell script with pybedtools

parser = argparse.ArgumentParser(description='Execute Augustus long running test cases.')
parser.add_argument('-p', '--predict', action='store_true',
                    help='to run original prediction.')
parser.add_argument('-e', '--prepare', action='store_true',
                    help='to build a new test set from scratch.')
parser.add_argument('-r', '--run', action='store_true',
                    help='to run prediction using minimal data set.')
parser.add_argument('-v', '--eval', action='store_true',
                    help='to evaluate accuracy (respect to the last prediction obtained by launching the script with --run option).')
parser.add_argument('-c-', '--chunk',
                    required=True,
                    help='a positive integer indicating the chunk to be processed (refer to documentation for a list of chunks over hg38.chr1).') 
parser.add_argument('-t', '--test', action='store_true',
                    help='to run a basic test to assess the correctness in the creation of the minimal data set.')                   
args = parser.parse_args()

# when moving on a different machine or changing paths to data, adjust the following 5+2 variables
augustus_dir = '../'
#'/home/giovanna/Desktop/githubAUGUSTUS/Augustus-TESTSET/'                        # path to augustus binaries
working_dir = '../examples/cgp12way/'
#'/home/giovanna/Desktop/githubAUGUSTUS/Augustus-TESTSET/examples/cgp12way/'       # path to working directory (it contains tree, genome tbl, SQLite db and there results will be written)
maf_dir = '../examples/cgp12way/MAF/'    # path to MAFs for chunks of interest
eval_dir = '/home/giovanna/Desktop/Alignment/eval-2.2.8/'                                       # path to eval
anno_file = '../examples/cgp12way/ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.gtf'  # path to annotation for hg38.chr1

# the following three directories are required only if a new test set is to be built
fasta_dir = '/home/giovanna/Desktop/Alignment/DATA_UCSCSOFT/'           # path to original FASTA files for genomes of interest
bedtools_dir = '/home/giovanna/Desktop/Alignment/bedtoolsBinaries/'     # path to bedtools (need to extract minimal FASTA from original FASTA, on the base of BED format)

# the following paths automatically reflect changes in the 6+3 paths above
sqlitedb_dir = working_dir + 'minimalFasta' + str(args.chunk) + '/'    # will contain minimal fasta after their extraction
sqlitedb_file = working_dir + 'SQLITE/12way.db'
sqlitedb_test_file = working_dir + 'SQLITE/12wayTEST_' + str(args.chunk) + '.db'
tbl_file = working_dir + 'GENOMETBL/genomes.tbl'
tbl_test_file = working_dir + 'GENOMETBL/genomesTEST_' + str(args.chunk) + '.tbl'
tree_file = working_dir + 'TREE/ucsc12way.nwk'
maf_file = maf_dir + "chr1_chunk_" + str(args.chunk) + '.maf' 
augustus_bin = augustus_dir + 'bin/augustus'
bedtool_bin = bedtools_dir + 'bedtools' 
eval_bin = eval_dir + 'evaluate_gtf.pl'
result_dir = '../examples/cgp12way/out' + str(args.chunk) + 'result/'


# if not already existing, create dir to collect results for the current chunk
def init():
    if os.path.exists(result_dir) == False:
        os.makedirs(result_dir)

# extract minimal FASTAs on the base of gene ranges (BEDs) and prepare sqlitedb
def make_sqlitedb():
    cleanup_db()

    print('Creating minimal FASTAs and SQLite database for current chunk...', bedtool_bin)

    # extract minimal FASTAs
    # todo replace the following with pybedtools
    subprocess.call([working_dir + 'extractFASTA.sh', bedtool_bin, fasta_dir, working_dir, args.chunk, sqlitedb_dir])
        
    for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
        # replace characters ':' and '-' with '_'
        for line in fileinput.input(sqlitedb_dir + species + '.MINIMAL.fasta', inplace=True):
            print(line.replace(':', "_").replace('-', "_"))

        # add sepcies to SQLite db
        cmd = [augustus_dir + 'bin/load2sqlitedb', '--dbaccess=' + sqlitedb_test_file, '--species=' + species, sqlitedb_dir + species + '.MINIMAL.fasta']
        execute(cmd, result_dir + 'out.createMinimalFASTA', mode ='a+')

# creating new genomes.tbl for this chunk
def make_genometbl():
    cleanup_tbl()
    print('Creating new genome tbl for current chunk...', tbl_test_file)

    with open(tbl_test_file, 'w') as f:
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            print(species, sqlitedb_dir + species + '.MINIMAL.fasta', sep='\t', file=f)
    
def cleanup_db():
    print('Cleaning up', sqlitedb_dir)
    if os.path.exists(sqlitedb_dir):
        shutil.rmtree(sqlitedb_dir)
    os.makedirs(sqlitedb_dir)

    print('Cleaning up', sqlitedb_test_file)
    if os.path.exists(sqlitedb_test_file):
        os.remove(sqlitedb_test_file)

def cleanup_tbl():
    print('Cleaning up', tbl_test_file)
    if os.path.exists(tbl_test_file):
        os.remove(tbl_test_file)

def run_prediction():
    print('Runnning prediction...')
    cmd = [augustus_bin, '--species=human', '--treefile=' + tree_file, '--alnfile=' + maf_file,
    '--speciesfilenames=' + tbl_file, '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + sqlitedb_file,
    '--stopCodonExcludedFromCDS=true', '--/CompPred/outdir=' + working_dir + 'out' + str(args.chunk) + 'prediction']

    execute(cmd, result_dir + 'out.runPrediction')

def prepare_test():
    print('Preparing test...')
    cmd = [augustus_bin, '--species=human', '--treefile=' + tree_file, '--alnfile=' + maf_file, '--dbaccess=' + sqlitedb_file, 
    '--speciesfilenames=' + tbl_file, '--softmasking=1', '--alternatives-from-evidence=0',
    '--stopCodonExcludedFromCDS=true', '--/CompPred/outdir=' + working_dir + 'out' + str(args.chunk) + 'prepare', '--/Testing/testMode=prepare']

    execute(cmd, result_dir + 'out.prepareTest')

    make_sqlitedb()
    make_genometbl()

def run_test():
    print('Runnning prediction using the minimal data set...')
    cmd = [augustus_bin, '--species=human', '--treefile=' + tree_file, '--alnfile=' + maf_file,
    '--speciesfilenames=' + tbl_test_file, '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + sqlitedb_test_file,
    '--stopCodonExcludedFromCDS=true', '--/CompPred/outdir=' + working_dir + 'out' + str(args.chunk) + 'run', '--/Testing/testMode=run']
    execute(cmd, result_dir + 'out.runTest')

def test_test():
    # test : prediction obtained working with minimal FASTAs is compared against original prediction, they should be identical for the test to succeed
    print('Runnning tests...')
    if os.path.exists(working_dir + 'out' + str(args.chunk) + 'prediction'):
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            pathToGFF = working_dir + 'out' + str(args.chunk) + 'prediction/' + species + '.cgp.gff'
            if os.path.exists(pathToGFF) == False:
                goahead = False
                print('Cannot find', pathToGFF, 'no test will be run over the code because original prediction is missing.') 
                return
    else:
        print('Cannot find', working_dir + 'out' + str(args.chunk) + 'prediction', 'no test will be run over the code because original prediction is missing.') 
        return

    if os.path.exists(working_dir + 'out' + str(args.chunk) + 'run'):
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            pathToGFF = working_dir + 'out' + str(args.chunk) + 'run/' + species + '.cgp.gff'
            if os.path.exists(pathToGFF) == False:
                goahead = False
                print('Cannot find', pathToGFF, 'no test will be run over the code because any new prediction over minimal data set is missing.') 
                return
    else:
        print('Cannot find', working_dir + 'out' + str(args.chunk) + 'run', 'no test will be run over the code because any new prediction over minimal data set is missing.') 
        return
    
    print('Assessing correctness of code...')
    correct = True
    for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
        pathToGFF_1 = working_dir + 'out' + str(args.chunk) + 'prediction/' + species + '.cgp.gff'
        pathToGFF_2 = working_dir + 'out' + str(args.chunk) + 'run/' + species + '.cgp.gff'
        if filecmp.cmp(pathToGFF_1, pathToGFF_2) == False:
            print('\tsomething went wrong: ', pathToGFF_1, 'differs from', pathToGFF_2)
            correct = False
    if correct:
        print('\ttest on code succeeded')

def run_evaluate():
    print('Runnning evaluation...')
    cmd = [eval_bin, anno_file, working_dir + 'out' + str(args.chunk) + 'run/hg38.cgp.gff']
    execute(cmd, result_dir + 'out.eval')
    find_values(result_dir + 'out.eval')

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
    with open(result_dir + 'eval.json', 'w') as file:
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

if __name__ == '__main__':
    init()
    if args.predict:
        run_prediction()
    if args.prepare:
        prepare_test()
    if args.run:
        run_test()
    if args.eval:
        run_evaluate()
    if args.test:
        test_test()