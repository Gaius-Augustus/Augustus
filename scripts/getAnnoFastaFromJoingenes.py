#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on August 26th 2019
#
# This Python script extracts CDS features from a GTF file, excises
# corresponding sequence windows from a genome FASTA file, stitches the
# codingseq parts together, adds letters N at the ends if bases are
# annotated as missing by frame in the GTF file, makes reverse complement
# if required, and translates to protein sequence.
# Output files are:
#    * file with protein sequences in FASTA format,
#    * file with coding squences in FASTA format
# The script automatically checks for in-frame stop codons and prints a
# warning to STDOUT if such genes are in the GTF-file. The IDs of bad genes
# are printed to a file bad_genes.lst. Option -s allows to exclude bad genes
# from the FASTA output file, automatically.
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
parser.add_argument('-t', '--table', dest='translation_table', default=1,
                    type=int, help='Translational code table number (INT)')
parser.add_argument('-s', '--filter_out_invalid_stops', dest='filter',
                    type=bool, help='Suppress output of protein sequences \
                    that contain internal stop codons.', default=False)
parser.add_argument('-p', '--print_format_examples', required=False, action='store_true',
                    help="Print gtf/gff3 input format examples, do not perform analysis")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', '--gtf',
                    type=str, help='file with CDS coordinates (GTF-format)')
group.add_argument('-3', '--gff3',
                    type=str, help='file with CDS coordinates (GFF3 format)')
args = parser.parse_args()


if args.print_format_examples:
    print('This script requires an annotation of protein coding genes with CDS ' +
        'features either in GTF or GFF3 format. Since both formats are sometimes ' +
        'inhomogeneous accross different tools, we here provide compatible ' +
        'examples. ' +
        'This script will only process the CDS lines. The input data may contain other feature lines that ' +
        'are ignored.')
    print('\nGTF format example:\n')
    print('NW_018027262_1\tAUGUSTUS\tCDS\t347\t525\t0.5\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";\n' +
          'NW_018027262_1\tAUGUSTUS\tCDS\t1791\t2242\t0.59\t-\t2\ttranscript_id "g1.t1"; gene_id "g1";\n' +
          'NW_018027262_1\tAUGUSTUS\tCDS\t2451\t2682\t0.24\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";\n')
    print('GFF3 format example:\n')
    print('NW_018027262_1\tAUGUSTUS\tCDS\t347\t525\t0.5\t-\t0\tID=g1.t1.CDS1;Parent=g1.t1;\n' +
          'NW_018027262_1\tAUGUSTUS\tCDS\t1791\t2242\t0.59\t-\t2\tID=g1.t1.CDS2;Parent=g1.t1;\n' +
          'NW_018027262_1\tAUGUSTUS\tCDS\t2451\t2682\t0.24\t-\t0\tID=g1.t1.CDS3;Parent=g1.t1;\n')
    print('\nThis script has successfully been tested with GTF format produced by BRAKER, and with' +
        'the GFF3 format produced by gtf2gff3.pl (Augustus/scripts).')
    exit(0)

# output file names:
codingseqFile = args.out + ".codingseq"
proteinFile = args.out + ".aa"

# Read GTF file CDS entries for transcripts
tx2seq = {}
tx2str = {}
cds = {}

if args.gtf:
    try:
        with open(args.gtf, "r") as gtf_handle:
            for line in gtf_handle:
                if re.match(
                        r"\S+\t[^\t]+\tCDS\t\d+\t\d+\t\S+\t\S+\t\d\t.*transcript_id (\S+)", line):
                    seq_id, st, en, stx, fr, tx_id = re.match(
                        r"(\S+)\t[^\t]+\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\d)\t.*transcript_id (\S+)", line).groups()
                    tx_id = re.sub(r'\"(\S+)\"', r'\1', tx_id)
                    tx_id = re.sub(r';', r'', tx_id)
                    if seq_id not in cds:
                        cds[seq_id] = {}
                    if tx_id not in cds[seq_id]:
                        cds[seq_id][tx_id] = []
                    cds[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx,
                         'frame': int(fr)})
                    if not tx_id in tx2seq:
                        tx2seq[tx_id] = seq_id
                        tx2str[tx_id] = stx
    except IOError:
        print("Error: Failed to open file " + args.gtf + "!")
        exit(1)
elif args.gff3:
    try:
        with open(args.gff3, "r") as gff3_handle:
            for line in gff3_handle:
                if re.match(r"\S+\t\S+\tCDS\t\d+\t\d+\t\S+\t\S+\t\d\tID=[^;]+;Parent=[^;]+;", line):
                    seq_id, st, en, stx, fr, tx_id = re.match(
                        r"(\S+)\t\S+\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\d)\tID=[^;]+;Parent=([^;]+);", line).groups()
                    if seq_id not in cds:
                        cds[seq_id] = {}
                    if tx_id not in cds[seq_id]:
                        cds[seq_id][tx_id] = []
                    cds[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx,
                         'frame': int(fr)})
                    if not tx_id in tx2seq:
                        tx2seq[tx_id] = seq_id
                        tx2str[tx_id] = stx
    except IOError:
        print("Error: Failed to open file " + args.gtf + "!")
        exit(1)
else:
    print("Error: Neither annotation file in GTF, nor in GFF3 was provided!")
    exit(1)

# Read genome file (single FASTA entries are held in memory, only), extract
# CDS sequence windows, add N when frame information states missing nucleotides
seq_len = {}
codingseq = {}
try:
    with open(args.genome, "r") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            seq_len[record.id] = len(record.seq)
            if record.id in cds:
                for tx in cds[record.id]:
                    if tx not in codingseq:
                        codingseq[tx] = SeqRecord(
                            Seq(""), id=tx, description=tx)
                    nCDS = len(cds[record.id][tx])
                    for i in range(0, nCDS):
                        cds_line = cds[record.id][tx][i]
                        if i == 0 and cds_line['strand'] == '+' and cds_line['frame'] != 0:
                            codingseq[tx].seq += Seq((3 - cds_line['frame'])
                                                     * 'N')
                        codingseq[tx].seq += record.seq[cds_line['start'] -
                                                        1:cds_line['end']]
                        if i == (nCDS - 1):
                            if cds_line['strand'] == '+':
                                if (len(codingseq[tx].seq) % 3) != 0:
                                    codingseq[tx].seq += Seq(
                                        (3 - (len(codingseq[tx]) % 3)) * 'N')
                        if i == (nCDS - 1) and cds_line['strand'] == '-':
                            if cds_line['frame'] != 0:
                                codingseq[tx].seq += Seq((3 - cds_line['frame'])
                                                         * 'N')
                            if(len(codingseq[tx].seq) % 3) != 0:
                                codingseq[tx].seq = Seq(
                                    (3 - (len(codingseq[tx].seq) % 3)) * 'N') + codingseq[tx].seq
                            codingseq[tx].seq = codingseq[
                                tx].seq.reverse_complement()
except IOError:
    print("Error: Failed to open file " + args.genome + "!")
    exit(1)

# Print coding sequences to file
try:
    with open(codingseqFile, "w") as codingseq_handle:
        for tx_id, seq_rec in codingseq.items():
            SeqIO.write(seq_rec, codingseq_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + codingseqFile + "!")
    exit(1)


# Translate coding sequences, identify sequences with in-frame stop codons,
# print proteins in FASTA format
bad_tx = {}
try:
    with open(proteinFile, "w") as protein_handle:
        for tx_id, seq_rec in codingseq.items():
            is_good = True
            seq_rec.seq = seq_rec.seq.translate(table=args.translation_table)
            match = re.search(r"(\*)(\w|\*)", str(seq_rec.seq))
            if match:
                # 1-based amino acid index of stop codon
                bad_tx[tx_id] = match.start(0) + 1
            if ((args.filter == True) and (tx_id not in bad_tx)) or args.filter == False:
                SeqIO.write(seq_rec, protein_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + proteinFile + "!")
    exit(1)

# Print IDs of genes with in-frame stop codons if any were found
if len(bad_tx) > 0:
    try:
        with open("bad_genes.lst", "w") as bad_handle:
            bad_handle.write(
                "# tx_id\tstop_in_amino_acid\tstop_in_mRNA_start\tstop_in_mRNA_end\tgenomic_transcript_start\tgenomic_transcript_end\tstrand\n")
            for tx in bad_tx:
                bad_handle.write(tx + "\t" + str(bad_tx[tx]) + "\t" + str(
                    bad_tx[tx] * 3 - 2) + "\t" + str(bad_tx[tx] * 3) + "\t")
                left_most = -1
                right_most = -1
                for i in range(0, len(cds[tx2seq[tx]][tx])):
                    if (left_most == -1) or (cds[tx2seq[tx]][tx][i]['start'] <= left_most):
                        left_most = cds[tx2seq[tx]][tx][i]['start']
                    if (right_most == -1) or (cds[tx2seq[tx]][tx][i]['end'] >= right_most):
                        right_most = cds[tx2seq[tx]][tx][i]['end']
                bad_handle.write(str(left_most) + '\t' + str(right_most) + '\t' + tx2str[tx] +
                                 '\n')
    except IOError:
        print("Error: Failed to open file bad_genes.lst!")
        exit(1)
    print("WARNING: The GTF file contained " + str(len(bad_tx)) + " gene(s) " +
          "with internal Stop codons (see file bad_genes.lst).")
    if args.filter:
        print("Note: Translations of the bad genes have already been excluded " +
              "from the output file " + args.gtf + ".")
    else:
        print("Note: Translations of the bad genes have been included in the " +
              "output file. If you wish to exclude them, run " +
              "getAnnoFastaFromJoingenes.py with option --filter!")
