#!/usr/bin/env python3

# This script can be used to identify predicted proteins that might contain
# repetitive sequence elements. Do not use it for filtering out genes
# without individual inspection! Some genes (e.g. spidroins) are repetitive
# and they are real genes!

import argparse
import itertools
import collections
from Bio import SeqIO
import pandas as pd

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2022. All rights reserved."
__license__ = "Artistic License"
__credits__ = "wwii from stackoverflow"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

parser = argparse.ArgumentParser(
    description='Screen a protein file in fasta format for amino acid ' + \
                'sequences that contain repetitive k-mers.')
parser.add_argument('-p', '--protein_file', required=True, type=str,
                    help='File with protein sequences in fasta format')
parser.add_argument('-o', '--out_file', required=True, type=str,
                    help='Output file: list of transcripts with high ' + \
                    'abundance of repetitive k-mers (tab separated).')
parser.add_argument('-k', '--max_k', required=False, type=int, default = 15,
                    help="Script looks for k-mers in protein sequences that " + \
                    "occur >=t times, starts with the largest possible k, " + \
                    "max_k, default value 15.")
parser.add_argument('-t', '--times', required=False, type=int, default = 4,
                    help="Threshold of replicates of a k-mer that must be " + \
                    "reached for reporting the k-mer, default value 10.")
args = parser.parse_args()

''' ******************* BEGIN FUNCTIONS *************************************'''

def k_freq(text, k, t):
    '''outputs frequency of overlapping kmers that occur more than t times 
    in text'''
    if len(text) >= k:
        k_mers = itertools.tee(text,k)
        # advance iterators
        n = 0
        for i in range(0,len(k_mers)):
            for j in range(1, n+1):
                next(k_mers[j])
            n = n+1

        aa_strs = enumerate(''.join(y) for y in zip(*k_mers))
        d = collections.defaultdict(list)
        for i,aa_str in aa_strs:
            d[aa_str].append(i)
        all_kmers = [(k,v) for k,v in d.items()]
        # count how many kmers occur in total in the sequence
        n_all = 0
        for kmer in all_kmers:
            n_all = n_all + len(kmer[1])
        # count how many kmers are repeats
        repeats = [(k,v) for k,v in d.items() if len(v)>t]
        n_rep = 0
        for kmer in repeats:
            n_rep = n_rep + len(kmer[1])
        return round(100*n_rep/n_all, 2)
    else:
        return 0

''' ******************* END FUNCTIONS ***************************************'''

## Create 3 lists for dataframe
txids = []
freqs = []
ks = []


## Read protein sequences one by one

with open(args.protein_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        #print("Processing record " + record.id)
        # screen for repetitive k-mers, break when significant repeat fround
        # start with largest k, go to k = 2
        # decided to apply diffrent thresholds without much testing
        for k in range(args.max_k, 1, -1):
            prop_reps = k_freq(record.seq, k, args.times)
            if prop_reps >= 50 and k > 5:
                txids.append(record.id)
                freqs.append(prop_reps)
                ks.append(k)
                #print(record.id + "\t" + str(prop_reps) + "\t" + str(k))
                break
            elif prop_reps >= 75 and k > 2 and k <= 5:
                txids.append(record.id)
                freqs.append(prop_reps)
                ks.append(k)
                #print(record.id + "\t" + str(prop_reps) + "\t" + str(k))
                break
            elif prop_reps >= 98 and k == 2:
                txids.append(record.id)
                freqs.append(prop_reps)
                ks.append(k)
                #print(record.id + "\t" + str(prop_reps) + "\t" + str(k))
                break

df = pd.DataFrame({'txid': txids, 'repeatcontent' : freqs, 'k': ks})

df = df.sort_values(['repeatcontent', 'k'], ascending = [False, False])

df.to_csv(args.out_file, sep = "\t", index = False)