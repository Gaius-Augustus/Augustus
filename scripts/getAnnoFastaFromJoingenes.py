#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on September 10th 2018
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

try:
    import argparse
except ImportError:
    raise ImportError('Failed to import argparse. Try installing with \"pip3 install argparse\"')

try:
    import re
except ImportError:
    raise ImportError('Failed to import argparse. Try installing with \"pip3 install re\"')

try:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna, generic_protein
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError('Failed to import biophython modules. Try installing with \"pip3 install biopython\"')


parser = argparse.ArgumentParser(
    description='Generate *.codingseq and *.aa FASTA-format files from genes \
                 in a GTF-file produced by AUGUSTUS auxprogs tool joingenes \
                 and a corresponding genomic FASTA-file.')
parser.add_argument('-g', '--genome', required=True,
                    type=str, help='genome sequence file (FASTA-format)')
parser.add_argument('-f', '--gtf', required=True,
                    type=str, help='file with CDS coordinates (GTF-format)')
parser.add_argument('-o', '--out', required=True, type=str,
                    help="name stem pf output file with coding sequences and \
                    protein sequences (FASTA-format); will be extended by \
                    .codingseq/.aa")
parser.add_argument('-t', '--table', dest='translation_table', default=1,
                    type=int, help='Translational code table number (INT)')
parser.add_argument('-s', '--filter_out_invalid_stops', dest='filter',
                    type=bool, help='Suppress output of protein sequences \
                    that contain internal stop codons.', default=False)
args = parser.parse_args()

# output file names:
codingseqFile = args.out + ".codingseq"
proteinFile = args.out + ".aa"

# Read GTF file CDS entries for transcripts
cds = {}
try:
    with open(args.gtf, "r") as gtf_handle:
        for line in gtf_handle:
            if re.match(
                    r"\S+\t\S+\tCDS\t\d+\t\d+\t\S+\t\S+\t\d\t.*transcript_id \"(\S+)\";", line):
                seq_id, st, en, stx, fr, tx_id = re.match(
                    r"(\S+)\t\S+\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\d)\t.*transcript_id \"(\S+)\";", line).groups()
                if seq_id not in cds:
                    cds[seq_id] = {}
                if tx_id not in cds[seq_id]:
                    cds[seq_id][tx_id] = []
                cds[seq_id][tx_id].append(
                    {'start': int(st), 'end': int(en), 'strand': stx,
                     'frame': int(fr)})
except IOError:
    print("Error: Failed to open file " + args.gtf + "!")

# Read genome file (single FASTA entries are held in memory, only), extract
# CDS sequence windows, add N when frame information states missing nucleotides
codingseq = {}
try:
    with open(args.genome, "rU") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            if record.id in cds:
                for tx in cds[record.id]:
                    if tx not in codingseq:
                        codingseq[tx] = SeqRecord(
                            Seq("", generic_dna), id=tx, description=tx)
                    nCDS = len(cds[record.id][tx])
                    for i in range(0, nCDS):
                        cds_line = cds[record.id][tx][i]
                        if i == 0 and cds_line['strand'] == '+' and cds_line['frame'] != 0:
                            codingseq[tx].seq += Seq((3 - cds_line['frame'])
                                                     * 'N', generic_dna)
                        codingseq[tx].seq += record.seq[cds_line['start'] -
                                                        1:cds_line['end']]
                        if i == (nCDS-1):
                            if cds_line['strand'] == '+':
                                if (len(codingseq[tx].seq) % 3) != 0:
                                    codingseq[tx].seq += Seq(
                                        (3 - (len(codingseq[tx]) % 3)) * 'N',
                                        generic_dna)
                        if i == (nCDS-1) and cds_line['strand'] == '-':
                            if cds_line['frame'] != 0:
                                codingseq[tx].seq += Seq((3 - cds_line['frame'])
                                                         * 'N', generic_dna)
                            if(len(codingseq[tx].seq) % 3) != 0:
                                codingseq[tx].seq = Seq(
                                    (3 - (len(codingseq[tx].seq) % 3)) * 'N',
                                    generic_dna) + codingseq[tx].seq
                            codingseq[tx].seq = codingseq[tx].seq.reverse_complement()
except IOError:
    print("Error: Failed to open file " + args.genome + "!")

# Print coding sequences to file
try:
    with open(codingseqFile, "w") as codingseq_handle:
        for tx_id, seq_rec in codingseq.items():
            SeqIO.write(seq_rec, codingseq_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + codingseqFile + "!")


# Translate coding sequences, identify sequences with in-frame stop codons,
# print proteins in FASTA format
bad_tx = {}
try:
    with open(proteinFile, "w") as protein_handle:
        for tx_id, seq_rec in codingseq.items():
            is_good = True
            seq_rec.seq = seq_rec.seq.translate(table=args.translation_table)
            if re.search(r"\*\w", str(seq_rec.seq)):
                bad_tx[tx_id] = 1
            if ((args.filter == True) and (tx_id not in bad_tx)) or args.filter == False:
                SeqIO.write(seq_rec, protein_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + proteinFile + "!")

# Print IDs of genes with in-frame stop codons if any were found
if len(bad_tx) > 0:
    try:
        with open("bad_genes.lst", "w") as bad_handle:
            for tx in bad_tx:
                bad_handle.write(tx + "\n")
    except IOError:
        print("Error: Failed to open file bad_genes.lst!")
    print("WARNING: The GTF file contained " + str(len(bad_tx)) + " gene(s) " +
          "with internal Stop codons (see file bad_genes.lst).")
    if args.filter:
        print("Note: Translations of the bad genes have already been exclude " +
              "from the output file " + args.gtf + ".")
    else:
        print("Note: Translations of the bad genes have been included in the " +
              "output file. If you wish to exclude them, run " +
              "getAnnoFastaFromJoingenes.py with option --filter!")
