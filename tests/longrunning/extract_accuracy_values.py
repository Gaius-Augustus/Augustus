#!/usr/bin/env python3

import re
import json
import os
import sys
import argparse


parser = argparse.ArgumentParser(
    description='Extract desired accuracy values from eval file(s).')
parser.add_argument('--path',
                    required=True,
                    help='The path to the folder or the file to consider.')         
args = parser.parse_args()



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
    resultfile = file.replace('.eval', '.json') 
    with open(resultfile, 'w') as file:
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


def analyze_folder(folder):
    if not os.path.isdir(args.path):
        print('Folder not found: ' + folder)
        sys.exit()

    for subdir, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.eval'):
                find_values(subdir + '/' + file)


if __name__ == '__main__':
    if args.path is None:
        print('Path to data folder containing eval files required, please make use of --path to pass the path...')
        sys.exit()

    if os.path.isdir(args.path):
        analyze_folder(args.path)
    else:
        find_values(args.path)
