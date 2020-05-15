#!/usr/bin/env python3

import subprocess
import re
import json
import os
import shutil
import argparse


#TODO: script arguments for eval path and wd
#TODO: define input data to use for lr test


parser = argparse.ArgumentParser(
    description='Execute Augustus long running test cases.')
parser.add_argument('--data',
                    required=True,
                    help='The directory where the data is stored.')
args = parser.parse_args()

augustusbin = '../bin/augustus'
tmp_data_folder = 'tmp/'
result_folder = 'output/'


def init():
    if os.path.exists(tmp_data_folder):
        shutil.rmtree(tmp_data_folder)

    if os.path.exists(result_folder):
        shutil.rmtree(result_folder)

    os.makedirs(tmp_data_folder)
    os.makedirs(result_folder)


def create_test_files():
    cmd = [augustusbin, '--species=human', args.data + '/example.fa']
    execute(cmd, tmp_data_folder + 'file1.gff')

    cmd = [augustusbin, '--species=chicken', args.data + '/example.fa']
    execute(cmd, tmp_data_folder + 'file2.gff')

    create_seqlist(tmp_data_folder + 'file2.gff')


def create_seqlist(file):
    with open(file, 'r') as input:
        match_list = re.findall('name = [0-9A-Za-z_]*', input.read())
        with open(tmp_data_folder + 'seqlist', "w") as file_out:
            for m in match_list:
                file_out.write(str(m).replace('name =', '').strip() + '\n')


def evaluate():
    create_test_files()
    cmd = [
        '../scripts/eval_multi_gtf.pl', tmp_data_folder + 'seqlist',
        tmp_data_folder + 'file1.gff', tmp_data_folder + 'file2.gff',
        '--EVAL_PATH=../../../eval-2.2.8'
    ]
    execute(cmd, tmp_data_folder + 'eval.out')


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
    with open(result_folder + 'eval.json', 'w') as file:
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


def execute(cmd, output):
    with open(output, 'w') as file:
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
    evaluate()
    find_values(tmp_data_folder + 'eval.out')
