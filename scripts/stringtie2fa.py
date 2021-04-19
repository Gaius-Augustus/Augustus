#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on April 19th 2021
#
# This Python script extracts exon features from a GTF file, excises
# corresponding sequence windows from a genome FASTA file, stitches the
# codingseq parts together, makes reverse complement
# if required
# Output file is:
#    * file with mRNA squences in FASTA format
# Beware: the script assumes that the gtf input file is sorted by coordinates!

try:
    import argparse
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install argparse\"')

try:
    import re
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install re\"')

try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError(
        'Failed to import biophython modules. Try installing with \"pip3 install biopython\"')


parser = argparse.ArgumentParser(
    description='Generate *.codingseq and *.aa FASTA-format files from genes \
                 in a GTF-file produced by AUGUSTUS auxprogs tool joingenes \
                 and a corresponding genomic FASTA-file.')
parser.add_argument('-g', '--genome', required=True,
                    type=str, help='genome sequence file (FASTA-format)')
parser.add_argument('-o', '--out', required=True, type=str,
                    help="name stem pf output file with coding sequences and \
                    protein sequences (FASTA-format); will be extended by \
                    .codingseq/.aa")
parser.add_argument('-p', '--print_format_examples', required=False, action='store_true',
                    help="Print gtf input format examples, do not perform analysis")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', '--gtf',
                    type=str, help='file with CDS coordinates (GTF-format)')
args = parser.parse_args()


if args.print_format_examples:
    print('This script requires an annotation of a transcript with exon ' +
        'features in GTF format. We here provide a compatible ' +
        'example. ' +
        'This script will only process the exon lines. The input data may contain other feature lines that ' +
        'are ignored.')
    print('\nGTF format example:\n')
    print('1\tStringTie\ttranscript\t1555206\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; cov "5.737374"; FPKM "5.261884"; TPM "18.775906";\n' +
         '1\tStringTie\texon\t1555206\t1555441\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "1"; cov "6.080509";\n' +
         '1\tStringTie\texon\t1565008\t1565346\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "2"; cov "5.917404";\n' +
         '1\tStringTie\texon\t1571901\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "3"; cov "4.533898";\n')
    print('\nThis script has successfully been tested with GTF format produced by Stringtie2.')
    exit(0)

# output file names:
mrnaFile = args.out + ".mrna"

# Read GTF file exon entries for transcripts
tx2seq = {}
tx2str = {}
mrna = {}

if args.gtf:
    try:
        with open(args.gtf, "r") as gtf_handle:
            for line in gtf_handle:
                if re.match(
                        r"\S+\t[^\t]+\texon\t\d+\t\d+\t\S+\t\S+\t\d\t.*transcript_id \"(\S+)\"", line):
                    seq_id, st, en, stx, tx_id = re.match(
                        r"(\S+)\t[^\t]+\texon\t(\d+)\t(\d+)\t\S+\t(\S+)\t\d\t.*transcript_id \"(\S+)\"", line).groups()
                    if seq_id not in mrna:
                        mrna[seq_id] = {}
                    if tx_id not in mrna[seq_id]:
                        mrna[seq_id][tx_id] = []
                    mrna[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx})
                    if not tx_id in tx2seq:
                        tx2seq[tx_id] = seq_id
                        tx2str[tx_id] = stx
    except IOError:
        print("Error: Failed to open file " + args.gtf + "!")
        exit(1)
else:
    print("Error: No annotation file in GTF format was provided!")
    exit(1)

# Read genome file (single FASTA entries are held in memory, only), extract
# exon sequence windows
seq_len = {}
mrnaseq = {}
try:
    with open(args.genome, "r") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            seq_len[record.id] = len(record.seq)
            if record.id in mrna:
                for tx in mrna[record.id]:
                    if tx not in mrnaseq:
                        if mrna[record.id][tx][0]['strand'] == '.':
                            descr = tx + ' strand_unknown'
                        else:
                            descr = tx
                        mrnaseq[tx] = SeqRecord(
                            Seq(""), id=tx, description=descr)
                    nExons = len(mrna[record.id][tx])
                    for i in range(0, nExons):
                        mrna_line = mrna[record.id][tx][i]
                        mrnaseq[tx].seq += record.seq[mrna_line['start'] -
                                                        1:mrna_line['end']]
                        if i == (nExons - 1) and mrna_line['strand'] == '-':
                            mrnaseq[tx].seq = mrnaseq[
                                tx].seq.reverse_complement()
except IOError:
    print("Error: Failed to open file " + args.genome + "!")
    exit(1)

# Print mRNA sequences to file
try:
    with open(mrnaFile, "w") as mrna_handle:
        for tx_id, seq_rec in mrnaseq.items():
            SeqIO.write(seq_rec, mrna_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + mrnaFile + "!")
    exit(1)
