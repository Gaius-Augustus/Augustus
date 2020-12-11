#!/usr/bin/env python3
import argparse
import fileinput
import sys
import os
import numpy as np

# author Giovanna Migliorelli
# version beta 02.06.2020

parser = argparse.ArgumentParser(description='Extract transcripts from some Ensembl annoation on the base of given chunks.')
parser.add_argument('-c', '--chunks',
                    #required=True,
                    nargs='+',
                    help='a list of one or more positive integers indicating the chunk/s for which we want to build a minimal annotation from given Ensembl.') 
args = parser.parse_args()


# helper function to extract values for some key (e.g. gene_id, transcript_id) from ENSEMBL like 8th field in some gff line
def extractTransId(str, key):
    value = ''
    for x in str.split(';'):
        pos = x.find(key)
        if pos>-1:
            value = x.split("\"")[1]
    return value

# reads in some annotation file and filters in those transcripts which fall in the interval [left,right]
def extractAnno(infile, outfile, mode, left, right, overlap, subset):
    transcripts = {}
    with open(infile, 'r') as input:
        while True:
            line = input.readline()
            if not line:
                break
            if line.find('#')==-1 and line!='\n':
                split = [x for x in line.split("\t") if x != ""]

                if split[0] == '1':
                    trans_id = extractTransId(split[8], 'transcript_id')
                    if trans_id != '':
                        if trans_id not in subset:
                            # just in case gff contains some unordered pair
                            start = min(int(split[3]), int(split[4]))
                            end = max(int(split[3]), int(split[4]))

                            if trans_id in transcripts:
                                transcripts[trans_id][0] = min(transcripts[trans_id][0], start)
                                transcripts[trans_id][1] = max(transcripts[trans_id][1], end)
                            else:
                                transcripts[trans_id] = [start, end]
                                                
              
    # read once again and filter (not smart but helps saving mem)
    with open(infile, 'r') as input, open(outfile, mode) as output:
        while True:
            line = input.readline()
            if not line:
                break
            if line.find('#')==-1 and line!='\n':                       # here comments and blank lines are purged, we can change it though
                split = [x for x in line.split("\t") if x != ""]
                trans_id = extractTransId(split[8], 'transcript_id')
                if trans_id != '':
                    if trans_id in transcripts:
                        if overlap == 'full':
                            if transcripts[trans_id][0]>=left and transcripts[trans_id][0]<=right:
                                output.write(line)
                                subset[trans_id] = True
                        elif overlap == 'partial':
                            if not(transcripts[trans_id][0]>=right or transcripts[trans_id][0]<=left):
                                output.write(line)
                                subset[trans_id] = True
                            
if __name__ == '__main__':
    if args.chunks is None:
        print('No chunks specified, please make use of --chunks to pass a non empty list of positive integers...')
        sys.exit()

    chunks = [int(x) for x in list(dict.fromkeys(args.chunks))]
    chunks = [x for x in chunks if x>0 and x<126]          # range valid for chr1 chunk size 2.5 Mb, chunk overlap 0.5 Mb

    if len(chunks) == 0:
        sys.exit()

    # should be synced with chunk_size, chunk_overlap, dir names used in the project see executeTestCGP.py for details
    chunk_size = 2500000
    chunk_overlap = 500000
    infile = '../examples/cgp12way/ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.gtf'
    outfile = '../examples/cgp12way/ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.FILTERED.gtf'
    
    subset = {} # used to keep track of those transcripts which have been already included and avoid duplicates
    
    print('Cleaning up old', outfile)
    if os.path.exists(outfile):
        os.remove(outfile)

    for chunk in chunks:
        left = (chunk_size-chunk_overlap)*(chunk-1)
        right = left + chunk_size - 1

        print('Including transcripts entirely contained within chunk', chunk, '[', left, ',', right, ']')
        extractAnno(infile, outfile, 'a+', left, right, 'full', subset)
