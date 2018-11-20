#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on November 15th 2018
#
# WARNINGS:
#
# 1) This script has been adapted to BRAKER output files.
#    It might fail e.g. at producing hints tracks from hints files that have
#    not been produced within BRAKER!
# 2) This script retrieves the binaries for linux 64 bit systems. It will
#    not work on other architectures and systems unless the required
#    UCSC tool binaries are already present.


import os
import errno
import os.path
import argparse
import re
import shutil
import platform
import urllib.request
import subprocess
from difflib import SequenceMatcher
try:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna, generic_protein
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError(
        'Failed to import biophython modules. Try installing with \"pip3 install biopython\"')

ucsc_tools = {'bedToBigBed': '', 'genePredCheck': '', 'faToTwoBit': '',
              'gtfToGenePred': '', 'hgGcPercent': '', 'ixIxx': '',
              'twoBitInfo': '', 'wigToBigWig': '', 'genePredToBed': ''}

parser = argparse.ArgumentParser(
    description='Generate UCSC assembly hub from BRAKER output.')
parser.add_argument('-e', '--email', required=True, type=str,
                    help='Contact e-mail adress for assembly hub')
parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('--no_repeats', required=False, type=bool, default=False,
                    help="Bool flag to disable repeat track generation from softmasked genome sequence")
parser.add_argument('--long_label', required=True, type=str,
                    help='Long label for hub, e.g. species name in english ' +
                    'and latin, pass in single quotation marks, e.g. ' +
                    '--long_label \'Dorosphila melanogster (fruit fly)\'')
parser.add_argument('--short_label', required=True, type=str,
                    help='Short label for hub, will also be used as directory ' +
                    'name for hub, should not contain spaces or special ' +
                    'characters, e.g. --short_label fly')
parser.add_argument('--annot', required=False, type=str,
                    help='GTF file with reference annotation')  # todo
parser.add_argument('--hints', required=False, type=str,
                    help='GFF file with AUGUSTUS hints')  # todo
parser.add_argument('--bam', required=False, type=str, nargs="+",
                    help='BAM file with RNA-Seq information')
parser.add_argument('--genemark', required=False, type=str,
                    help='GTF file with GeneMark predictions')  # todo
parser.add_argument('--aug_ab_initio', required=False, type=str,
                    help='GTF file with ab initio AUGUSTUS predictions')  # todo
parser.add_argument('--aug_hints', required=False, type=str,
                    help='GTF file with AUGUSTUS predictions with hints')  # todo
parser.add_argument('--aug_ab_initio_utr', required=False, type=str,
                    help='GTF file with ab initio AUGUSTUS predictions with UTRs')  # todo
parser.add_argument('--aug_hints_utr', required=False, type=str,
                    help='GTF file with AUGUSTUS predictions with hints with UTRs')  # todo
parser.add_argument('--traingenes', required=False, type=str,
                    help='GTF file with training genes')  # todo
parser.add_argument('--SAMTOOLS_PATH', required=False, type=str,
                    help="Path to samtools executable")
parser.add_argument('--outdir', required=False, type=str, default='.',
                    help="output directory to write hub to")
args = parser.parse_args()

''' Find samtools (if bam file provided) '''
samtools = ""
if args.bam:
    print("Searching for samtools:")
if args.bam and args.SAMTOOLS_PATH:
    samtools = args.SAMTOOLS_PATH + "/samtools"
    if not(os.access(samtools, os.X_OK)):
        print("Error: " + samtools + " is not executable!")
        exit(1)
    else:
        print("Will use " + samtools)
elif args.bam:
    if shutil.which("samtools") is not None:
        samtools = shutil.which("samtools")
        print("Will use " + samtools)
    else:
        print("Error: unable to locate samtools binary!")
        exit(1)


''' Find or obitain UCSC tools '''
# the URLs of UCSC tool download are hardcoded for linux.x84_64
arch = platform.machine()

print("Searching for required UCSC tools:")
for key, val in ucsc_tools.items():
    if shutil.which(key) is not None:
        ucsc_tools[key] = shutil.which(key)
    elif os.path.isfile(os.getcwd() + "/" + key):
        ucsc_tools[key] = os.getcwd() + "/" + key
        if not(os.access(ucsc_tools[key], os.X_OK)):
            os.chmod(ucsc_tools[key], 0o777)
    else:
        if not(arch == 'x86_64'):
            print(arch)
            print("Error: This script was implemented for linux.x84_64 " +
                  " architecture and will not be able to locate the exectuables " +
                  "for your system. Please download " + key + ", manually. UCSC" +
                  "tools are generally available at " +
                  "http://hgdownload.soe.ucsc.edu/admin/exe")
            exit(1)
        else:
            tool_url = "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/" + key
            print("Was unable to locate " + key +
                  ", will try to download it from " + tool_url)
            print("This may take a while...")
            with urllib.request.urlopen(tool_url) as response, open(key, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            ucsc_tools[key] = os.getcwd() + "/" + key
            os.chmod(ucsc_tools[key], 0o777)
    print("Will use " + ucsc_tools[key])


''' track color defintion '''

col_idx = 0
rgb_cols = ['0,0,0', '255,0,0', '0,255,0', '0,0,255', '176,196,222', '0,255,255', '255,0,255',
            '192,192,192', '128,128,128', '128,0,0', '128,128,0', '0,128,0', '128,0,128',
            '0,128,128', '0,0,128', '139,0,0', '220,20,60', '233,150,122', '255,140,0',
            '218,165,32', '154,205,50', '34,139,34', '152,251,152', '72,209,203', '176,224,230',
            '138,43,226']  # 0 - 25

''' Create hub directory structure '''

try:
    os.makedirs(args.outdir + "/" + args.short_label + "/" + args.short_label)
    os.makedirs(args.outdir + "/" + args.short_label + "/tmp")
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

''' Function that runs a subprocess with arguments '''


def run_simple_process(args_lst):
    try:
        print("Trying to execute the following command:")
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            print("Error: return code of subprocess was " +
                  str(result.returncode))
    except subprocess.CalledProcessError as grepexc:
        print("Failed executing: ", " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


''' Function that runs a subprocess with input from STDIN '''


def run_process_stdinput(args_lst, byte_obj):
    try:
        print("Trying to execute the following command with input from STDIN:")
        print(" ".join(args_lst))
        result = subprocess.run(args_lst, input=byte_obj,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            print("Error: return code of subprocess was "
                  + str(result.returncode))
    except subprocess.CalledProcessError as grepexc:
        print("Failed executing: ", " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


''' Function that sanity checks and fixes a gtf file, e.g. for strand problems '''


def make_gtf_sane(annot_file, ucsc_file):
    try:
        with open(annot_file, "r") as annot_handle:
            txs = {}
            for line in annot_handle:
                seq, first_part, strand, second_part, txid = re.search(r'(\S+)(\t\S+\t\S+\t\d+\t\d+\t\S+\t)(\S+)(\t\S+\t).*transcript_id\s\"(\S+)\"', line).groups()
                if txid not in txs:
                    txs[txid] = {'strand' : strand, 'lines' : [], 'sane' : True}
                else:
                    if not(strand == txs[txid]['strand']):
                        txs[txid]['sane'] = False
                txs[txid]['lines'].append(line)
    except IOError:
        print("Error: failed to open file " + annot_file + " for reading!")
    try:
        with open(ucsc_file, "w") as ucsc_handle:
            for key, value in txs.items():
                if value['sane'] == True:
                    for line in value['lines']:
                        ucsc_handle.write(line)
                else:
                    for line in value['lines']:
                        seq, first_part, strand, second_part, txid = re.search(r'^(\S+)(\t\S+\t\S+\t\d+\t\d+\t\S+\t)(\S+)(\t\S+\t).*transcript_id\s\"(\S+)\"', line).groups()
                        new_gid = txid + "_" + seq + "_" + strand
                        new_txid = new_gid + ".t1"
                        ucsc_handle.write(seq + first_part + strand + second_part + "gene_id \"" + new_gid + "\"; transcript_id \"" + new_txid + "\";\n")
    except IOError:
        print("Error: failed to open file " + ucsc_file + " for writing!")



''' Function that reformats AUGUSTUS gtf format to UCSC gtf format '''


def aug2ucsc_gtf(augustus_file, ucsc_file):
    try:
        with open(augustus_file, "r") as aug_handle:
            try:
                with open(ucsc_file, "w") as ucsc_handle:
                    for line in aug_handle:
                        if re.search(r'\tAUGUSTUS\tCDS\t', line) or re.search(r'\tAUGUSTUS\texon\t', line) or re.search(r'\tAUGUSTUS\tstart_codon\t', line) or re.search(r'\tAUGUSTUS\t\d+\'-UTR\t', line):
                            if re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\ttranscript_id\s\"\S+\";\sgene_id\s\"\S+\";', line):
                                first_part, feature, second_part, gid, txid = re.search(r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(transcript_id\s\"\S+\";)\s(gene_id\s\"\S+\";)', line).groups()
                            elif re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\tgene_id\s\"\S+\";\stranscript_id\s\"\S+\";', line):
                                first_part, feature, second_part, txid, gid = re.search(r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(gene_id\s\"\S+\";)\s(transcript_id\s\"\S+\";)', line).groups()
                            if re.search(r'5\'-UTR', feature):
                                feature = "5UTR"
                            elif re.search(r'3\'-UTR', feature):
                                feature = "3UTR"
                            ucsc_handle.write(first_part + feature + second_part + gid + " " + txid + "\n")
            except IOError:
                print("Error: failed to open file " +
                      ucsc_file + " for writing!")
    except IOError:
        print("Error: failed to open file " + augustus_file + " for reading!")


''' Function that writes subprocess byte object to flat file '''


def write_byteobj(byte_obj, outfile):
    try:
        with open(outfile, 'w') as byteobj_handle:
            byteobj_handle.write(byte_obj.decode('utf-8'))
    except IOError:
        print("Error: failed to open file " + outfile + " for writing!")


''' Function that determines which regions in genome are softmasked '''


def find_masked_intervals(genome_file, bed3_file):
    try:
        masked_intervals = []
        with open(genome_file, "rU") as genome_handle:
            for record in SeqIO.parse(genome_handle, "fasta"):
                masked_seq = str(record.seq)
                inMasked = False
                mmEnd = 0
                for x in range(0, len(masked_seq)):
                    if masked_seq[x].islower() and inMasked == False:
                        mmStart = x + 1
                        inMasked = True
                    elif not(masked_seq[x].islower()) and inMasked == True:
                        mmEnd = x
                        inMasked = False
                        if x > 0 and mmEnd > 0:
                            masked_intervals.append(
                                {'id': record.id, 'start': mmStart, 'end': mmEnd})
                if inMasked == True:
                    masked_intervals.append(
                        {'id': record.id, 'start': mmStart, 'end': len(masked_seq)})
    except IOError:
        print("Error: failed to open file " + genome_file + " for reading!")

    try:
        with open(bed3_file, "w+") as bed2_handle:
            for x in masked_intervals:
                bed2_handle.write(
                    x['id'] + '\t' + str(x['start']) + '\t' + str(x['end']) +
                    '\n')
    except IOError:
        print("Error: failed to open file " + bed3_file + " for writing!")


''' Function that sorts a bed3 file with LC_COLLATE=C (using a bash script) '''


def sort_bed3(bed3_file, bed3_sorted_file):
    script_file = args.outdir + '/sort_bed3.sh'
    try:
        with open(script_file, 'w') as script_handle:
            script_handle.write(
                '#!/bin/bash\nLC_COLLATE=C\nsort -k1,1 -k2,2n ' + bed3_file +
                ' > ' + bed3_sorted_file)
    except IOError:
        print("Error: failed to open file " + script_file + " for writing!")
    os.chmod(script_file, 0o777)
    subprcs_args = [script_file]
    run_simple_process(subprcs_args)


''' Function that converts bed to bigBed '''

def bed2bigBed(btype, bed_file, chrom_sizes_file, bigBedFile):
    print('Generating bigBed file for ' + bed_file + '...')
    subprcs_args = [ucsc_tools['bedToBigBed'], '-type=bed' +
                    str(btype), bed_file, chrom_sizes_file, bigBedFile]
    run_simple_process(subprcs_args)
    print('Done.')


''' Function that converts gtf to bb format '''

def gtf2bb(gtf_file, gp_file, bed_file, bb_file, info_out_file, chrom_size_file):
    subprcs_args = [ucsc_tools['gtfToGenePred'], '-infoOut=' +
                    info_out_file, '-genePredExt', gtf_file, gp_file]
    run_simple_process(subprcs_args)
    subprcs_args = [ucsc_tools['genePredCheck'], gp_file]
    result = run_simple_process(subprcs_args)
    # parse result for failed annotations
    annotation_validation_result = result.stderr.decode('utf-8')
    print(annotation_validation_result)
    regex_result = re.match(
        r'checked: \d+ failed: (\d+)', annotation_validation_result)
    if int(regex_result.group(1)) > 0:
        print("Warning: " + regex_result.group(1) +
              " annotations did not pass validation!")
    subprcs_args = [ucsc_tools['genePredToBed'], gp_file, 'stdout']
    result = run_simple_process(subprcs_args)
    sort_tool = shutil.which('sort')
    if sort_tool is None:
        print("Error: unable to locate bash tool 'sort'")
    subprcs_args = [sort_tool, '-k1,1', '-k2,2n']
    result = run_process_stdinput(subprcs_args, result.stdout)
    write_byteobj(result.stdout, bed_file)
    bed2bigBed(12, bed_file, chrom_size_file, bb_file)


''' Function that writes info about gene pred or hints to trackDb file '''


def info_to_trackDB(trackDb_file, short_label, long_label, rgb_color, group, bed_no):
    try:
        with open(trackDb_file, "a") as trackDb_handle:
            trackDb_handle.write("track " + short_label + "\n" +
                                 "longLabel " + long_label + "\n" +
                                 "shortLabel " + short_label + "\n" +
                                 "group " + group + "\ntype bigBed " + str(bed_no) + " .\n" +
                                 "bigDataUrl " + short_label + ".bb\n" +
                                 "color " + rgb_color + "\n" +
                                 "visibility 4\n\n" +
                                 "group " + group + "\ntype bigBed " + str(bed_no) + " .\n" +
                                 "bigDataUrl " + short_label + ".bb\n" +
                                 "color " + rgb_color + "\n\n")
    except IOError:
        print("Error: failed to open file " +
              trackDb_file + " for writing!")


''' Function that converts gtf to gene prediction track '''


def make_gtf_track(trackDb_file, gtf_file, chrom_size_file, short_label, long_label, rgb_color):
    gp_file = args.outdir + "/" + args.short_label + "/tmp/" + short_label + ".gp"
    info_out_file = args.outdir + "/" + args.short_label + \
        "/tmp/" + short_label + ".infoOut.txt"
    bed_file = args.outdir + "/" + args.short_label + "/tmp/" + short_label + ".bed"
    bb_file = args.outdir + "/" + args.short_label + \
        "/" + args.short_label + "/" + short_label + ".bb"
    gtf2bb(gtf_file, gp_file, bed_file, bb_file,
           info_out_file, chrom_size_file)
    info_to_trackDB(trackDb_file, short_label, long_label,
                    rgb_color, "genePreds", 12)
    # parse info_out_file to produce txt file for creating nameIndex files
    name_index_txt_file = args.outdir + "/" + args.short_label + \
        "/tmp/" + short_label + ".nameIndex.txt"
    try:
        with open(name_index_txt_file, "w") as name_index_handle:
            try:
                with open(info_out_file, "r") as info_out_handle:
                    for line in info_out_handle:
                        if not(re.match(r'^\#', line)) and (re.match(r'(\S+)\s+(\S+)\s+(\S+)\s+', line)):
                            tx_id, g_id, src = re.match(
                                r'(\S+)\s+(\S+)\s+(\S+)\s+', line).groups()
                            name_index_handle.write(
                                tx_id + "\t" + g_id + "," + src + ",,,\n")
            except IOError:
                print("Error: failed to open file " +
                      info_out_file + " for reading!")
    except IOError:
        print("Error: failed to open file " +
              name_index_txt_file + " for writing!")
    ix_file = args.outdir + "/" + args.short_label + \
        "/" + args.short_label + "/" + short_label + ".nameIndex.ix"
    ixx_file = args.outdir + "/" + args.short_label + \
        "/" + args.short_label + "/" + short_label + ".nameIndex.ixx"
    subprcs_args = [ucsc_tools['ixIxx'],
                    name_index_txt_file, ix_file, ixx_file]
    run_simple_process(subprcs_args)


''' Generate essential files for genome display '''


TwoBit_file = args.outdir + "/" + args.short_label + "/" + \
    args.short_label + "/" + args.short_label + ".2bit"
print('Generating genome 2bit file ' + TwoBit_file + '...')
subprcs_args = [ucsc_tools['faToTwoBit'], args.genome, TwoBit_file]
run_simple_process(subprcs_args)
print('Done.')

ChromSizes_file = args.outdir + "/" + args.short_label + \
    "/tmp/" + args.short_label + ".chrom.sizes"
print('Generating chromsome size info file ' + ChromSizes_file + '...')
subprcs_args = [ucsc_tools['twoBitInfo'], TwoBit_file, 'stdout']
result = run_simple_process(subprcs_args)
sort_tool = shutil.which('sort')
if sort_tool is None:
    print("Error: unable to locate bash tool 'sort'")
subprcs_args = [sort_tool, '-k2rn']
result = run_process_stdinput(subprcs_args, result.stdout)
write_byteobj(result.stdout, ChromSizes_file)
print('Done.')

WigVarStep_file = args.outdir + "/" + args.short_label + \
    "/tmp/" + args.short_label + ".gc5Base.wigVarStep"
WigVarStep_file_compr = WigVarStep_file + ".gz"
print('Generating variable step wiggle file for GC content ' +
      WigVarStep_file_compr + '...')
subprcs_args = [ucsc_tools['hgGcPercent'], '-wigOut', '-doGaps',
                '-file=stdout', '-win=5', '-verbose=0', args.short_label,
                TwoBit_file]
result = run_simple_process(subprcs_args)
write_byteobj(result.stdout, WigVarStep_file)
gzip_tool = shutil.which('gzip')
if gzip_tool is None:
    print("Error: unable to locate bash tool 'gzip'")
if os.path.isfile(WigVarStep_file_compr):
    os.unlink(WigVarStep_file_compr)
subprcs_args = [gzip_tool, WigVarStep_file]
run_simple_process(subprcs_args)
print('Done.')

BigWigGC_file = args.outdir + "/" + args.short_label + "/" + \
    args.short_label + "/" + args.short_label + ".gc5Base.bw"
print('Generating bigWig file for GC content ' + BigWigGC_file + '...')
subprcs_args = [ucsc_tools['wigToBigWig'],
                WigVarStep_file_compr, ChromSizes_file, BigWigGC_file]
run_simple_process(subprcs_args)
print('Done.')

hub_txt_file = args.outdir + "/" + args.short_label + "/hub.txt"
try:
    with open(hub_txt_file, "w+") as hub_txt_handle:
        hub_txt_handle.write("hub " + args.short_label + "\n")
        hub_txt_handle.write("shortLabel " + args.short_label + "\n")
        hub_txt_handle.write("longLabel " + args.long_label + "\n")
        hub_txt_handle.write("genomesFile genomes.txt\n")
        hub_txt_handle.write("email " + args.email + "\n")
        hub_txt_handle.write("descriptionUrl aboutHub.html\n")
except IOError:
    print("Error: failed to open file " + hub_txt_file + " for writing!")

default_seq_id = ""
default_seq_end = 0
try:
    with open(args.genome, "rU") as genome_handle:
        nSeq = 0
        for record in SeqIO.parse(genome_handle, "fasta"):
            default_seq_id = record.id
            if len(record.seq) > 15000:
                default_seq_end = 15000
            else:
                default_seq_end = len(record.seq)
            nSeq = nSeq + 1
            if nSeq > 0:
                continue
except IOError:
    print("Error: failed to open file " + args.genome + " for reading!")

genomes_txt_file = args.outdir + "/" + args.short_label + "/genomes.txt"
try:
    with open(genomes_txt_file, "w+") as genomes_txt_handle:
        genomes_txt_handle.write("genome " + args.short_label + "\n")
        genomes_txt_handle.write(
            "trackDb " + args.short_label + "/trackDb.txt\n")
        genomes_txt_handle.write(
            "groups " + args.short_label + "/groups.txt\n")
        genomes_txt_handle.write(
            "description Automatically generated Hub\n")
        genomes_txt_handle.write(
            "twoBitPath " + args.short_label + "/" + args.short_label +
            ".2bit\n")
        genomes_txt_handle.write("organism " + args.long_label + "\n")
        genomes_txt_handle.write(
            "defaultPos " + default_seq_id + ":1-" + str(default_seq_end) +
            "\n")
        genomes_txt_handle.write("orderKey 4800\n")
        genomes_txt_handle.write("scientificName " + args.long_label + "\n")
        genomes_txt_handle.write("htmlPath hmi/description.html\n")
except IOError:
    print("Error: failed to open file " + genomes_txt_file + " for writing!")

trackDb_file = args.outdir + "/" + args.short_label + \
    "/" + args.short_label + "/trackDb.txt"
try:
    with open(trackDb_file, "w+") as trackDb_handle:
        trackDb_handle.write("track gcPercent\nlongLabel GC Percent in 5-base " +
                             "Window\nshortLabel GC Percent\n" +
                             "type bigWig 0 100\ngroup map\nvisibility dense" +
                             "\nwindowingFunction Mean\nbigDataUrl " +
                             args.short_label + ".gc5Base.bw\npriority 2\nautoScale Off\n" +
                             "maxHeightPixels 128:36:16\ngraphTypeDefault Bar\ngridDefault OFF\n" +
                             "ncolor 0,0,0\naltColor 128,128,128\nviewLimits 30:70\nhtml ../documentation/gcPercent\n\n")
except IOError:
    print("Error: failed to open file " + trackDb_file + " for writing!")

groups_txt_file = args.outdir + "/" + args.short_label + \
    "/" + args.short_label + "/groups.txt"
try:
    with open(groups_txt_file, "w+") as groups_handle:
        groups_handle.write(
            "name genePreds\nlabel Gene Predictions\npriority 2\ndefaultIsClosed 0\n\n")
        groups_handle.write(
            "name reps\nlabel Repeats\npriority 2\ndefaultIsClosed 0\n\n")
        groups_handle.write(
            "name hints\nlabel Hints\npriority 2\ndefaultIsClosed 0\n\n")
except IOError:
    print("Error: failed to open file " + groups_txt_file + " for writing!")


''' Generate repeat masking track '''

if not args.no_repeats:
    softmaskedBed3_file = args.outdir + "/" + args.short_label + \
        "/tmp/" + args.short_label + ".RMsoft.bed3"
    print('Generating softmasking information bed3 file ' +
          softmaskedBed3_file + " from genome data (this may take a while)...")
    find_masked_intervals(args.genome, softmaskedBed3_file)
    print('Done.')
    softmaskedBed3_sorted_file = args.outdir + "/" + \
        args.short_label + "/tmp/" + args.short_label + ".RMsoft.s.bed3"
    print('Sorting file ' + softmaskedBed3_file + '...')
    sort_bed3(softmaskedBed3_file, softmaskedBed3_sorted_file)
    print('Done.')
    ChromSizes_file = args.outdir + "/" + args.short_label + \
        "/tmp/" + args.short_label + ".chrom.sizes"
    softmaskedBigBed_file = args.outdir + "/" + args.short_label + "/" + \
        args.short_label + "/" + args.short_label + "_softmasking.bb"
    bed2bigBed(3, softmaskedBed3_sorted_file,
               ChromSizes_file, softmaskedBigBed_file)
    try:
        with open(trackDb_file, "a") as trackDb_handle:
            trackDb_handle.write("track RMsoft\nlongLabel Softmaked Repeats\n" +
                                 "shortLabel Repeats\ngroup reps\ntype bigBed 3 .\n" +
                                 "bigDataUrl " + args.short_label + "_softmasking.bb\n" +
                                 "color " + rgb_cols[col_idx] + "\n\ngroup reps\ntype bigBed 3 .\n" +
                                 "bigDataUrl " + args.short_label + "_softmasking.bb\n" +
                                 "color " + rgb_cols[col_idx] + "\n\n")
            col_idx = col_idx + 1
            if col_idx > 25:
                col_idx = 0
    except IOError:
        print("Error: failed to open file " + trackDb_file + " for writing!")


''' Creating RNA-Seq bam track(s) '''


if args.bam:
    print('Generating BAM track(s)...')
    bam_index = 1
    for bam_file in args.bam:
        bam_sorted_file = args.outdir + "/" + args.short_label + "/" + \
            args.short_label + "/" + "rnaseq_" + str(bam_index) + \
            ".s.bam"
        subprcs_args = [samtools, "sort", bam_file, "-o", bam_sorted_file]
        run_simple_process(subprcs_args)
        bam_index_file = args.outdir + "/" + args.short_label + \
            "/" + args.short_label + "/" + \
            "rnaseq_" + str(bam_index) + ".s.bam.bai"
        subprcs_args = [samtools, "index",
                        bam_sorted_file, bam_index_file]
        run_simple_process(subprcs_args)
        try:
            with open(trackDb_file, "a") as trackDb_handle:
                trackDb_handle.write("track RNASeq_" + str(bam_index) + "\n" +
                                     "bigDataUrl rnaseq_" + str(bam_index) + ".s.bam\n" +
                                     "shortLabel RNASeq_" + str(bam_index) + "\n" +
                                     "longLabel RNASeq bam file " + str(bam_index) + " from file " +
                                     bam_file + "\ngroup hints\nvisibility 4\ntype bam\n\n" +
                                     "group hints\nbigDataUrl rnaseq_" + str(bam_index) + ".s.bam\n")
        except IOError:
            print("Error: failed to open file " +
                  trackDb_file + " for writing!")
        bam_index = bam_index + 1
    print('Done.')


''' Creating reference annotation gene track '''

if args.annot:
    print('Generating reference annotation gene track...')
    ucsc_file = args.outdir + "/" + args.short_label + "/tmp/annot_ucsc.gtf"
    make_gtf_sane(args.annot, ucsc_file)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "annot",
                   "Reference Annotation", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating genemark prediction track '''

if args.genemark:
    print('Generating GeneMark prediction track...')
    make_gtf_track(trackDb_file, args.genemark, ChromSizes_file, "genemark",
                   "GeneMark predictions", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating AUGUSTUS ab initio track '''

if args.aug_ab_initio:
    print('Generating AUGUSTUS ab initio prediction (no UTRs) track...')
    ucsc_file = args.outdir + "/" + args.short_label + "/tmp/aug_ab_initio_ucsc.gtf"
    aug2ucsc_gtf(args.aug_ab_initio, ucsc_file)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "aug_ab_initio_no_utr",
                   "AUGUSTUS ab initio predictions without UTRs", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating AUGUSTUS with hints track '''

if args.aug_hints:
    print('Generating AUGUSTUS prediction with hints (no UTRs) track...')
    ucsc_file = args.outdir + "/" + args.short_label + "/tmp/aug_hints_ucsc.gtf"
    aug2ucsc_gtf(args.aug_hints, ucsc_file)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "aug_hints_no_utr",
                   "AUGUSTUS predictions with hints without UTRs", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating AUGUSTUS ab initio with UTRs track '''

if args.aug_ab_initio_utr:
    print('Generating AUGUSTUS ab initio prediction with UTRs track...')
    ucsc_file = args.outdir + "/" + args.short_label + "/tmp/aug_ab_initio_utr_ucsc.gtf"
    aug2ucsc_gtf(args.aug_ab_initio_utr, ucsc_file)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "aug_ab_initio_utr",
                   "AUGUSTUS ab initio predictions with UTRs", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating AUGUSTUS with hints and UTRs track '''

if args.aug_hints_utr:
    print('Generating AUGUSTUS prediction with hints and with UTRs track...')
    ucsc_file = args.outdir + "/" + args.short_label + "/tmp/aug_hints_utr_ucsc.gtf"
    aug2ucsc_gtf(args.aug_hints_utr, ucsc_file)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "aug_hints_utr",
                   "AUGUSTUS predictions with hints and UTRs", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating traingenes track '''

if args.traingenes:
    print('Generating training gene track...')
    make_gtf_track(trackDb_file, args.traingenes, ChromSizes_file, "traingenes",
                   "Training genes", rgb_cols[col_idx])
    col_idx = col_idx + 1
    if col_idx > 25:
        col_idx = 0
    print('Done.')


''' Creating hints tracks '''

if args.hints:
    print('Separating hints into separate tracks...')
    hint_categ = {}
    try:
        with open(args.hints, "r") as hints_handle:
            for line in hints_handle:
                h_src, h_type = re.match(
                    r'^\S+\t(\S+)\t(\S+)\t', line).groups()
                if h_src not in hint_categ:
                    hint_categ[h_src] = {}
                    hint_categ[h_src][h_type] = []
                else:
                    if h_type not in hint_categ[h_src]:
                        hint_categ[h_src][h_type] = []
                hint_categ[h_src][h_type].append(line)
    except IOError:
        print("Error: failed to open file " + args.hints + " for reading!")
    for h_src in hint_categ:
        for h_type in hint_categ[h_src]:
            this_hints_file = args.outdir + "/" + args.short_label + "/tmp/" + \
                h_src + "_" + h_type + ".hints"
            try:
                with open(this_hints_file, "w") as hints_handle:
                    hint_no = 0
                    for hint in hint_categ[h_src][h_type]:
                        hint_no = hint_no + 1
                        if (h_type == 'CDS') or (h_type == 'CDSpart') or (h_type == 'exon'):
                            # GTF like hints (all treated as exons)
                            if re.match(r'(\S+\t\S+\t)\S+(\t\S+\t\S+\t\S+\t\S+\t\S+\t)', hint):
                                first_part, second_part = re.match(
                                    r'(\S+\t\S+\t)\S+(\t\S+\t\S+\t\S+\t\S+\t\S+\t)', hint).groups()
                                if re.search(r'grp=', hint):
                                    seq_id, strand, grp = re.search(
                                        r'(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t\S+\t\S*grp=([^;]+)', hint).groups()
                                    tx_id = seq_id + "_" + grp + "_" + strand
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id" +
                                                       "\" " + tx_id + "\"; transcript_id \"" + tx_id + ".t1\"\n")
                                elif re.search(r'group=', hint):
                                    seq_id, strand, grp = re.search(
                                        r'(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t\S+\t\S*group=([^;]+)', hint).groups()
                                    tx_id = seq_id + "_" + grp + "_" + strand
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id" +
                                                       "\" " + grp + "\"; transcript_id \"" + grp + ".t1\"\n")
                                elif re.search(r'mult=', hint):
                                    mult = re.search(
                                        r'mult=(\d+)', hint).group(1)
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id \"hint_" + str(hint_no) +
                                                       "_mult_" + mult + "\"; transcript_id \"hint_" + str(hint_no) + "_mult_" + mult + ".t1\"\n")
                                else:
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id \"hint_" + str(hint_no) +
                                                       "_mult_1\"; transcript_id \"hint_" + str(hint_no) + "_mult_1.t1\"\n")

                        else:
                            # segmental hints, i.e. introns
                            if re.match(
                                    r'(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t', hint):
                                seq, start, end, strand = re.match(
                                    r'^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t', hint).groups()
                                if re.search(r'mult=', hint):
                                    mult = re.search(
                                        r'mult=(\d+)', hint).group(1)
                                else:
                                    mult = 1
                                hints_handle.write(seq + "\t" + start + "\t" + end + "\t" + h_type + "_" + str(hint_no) + "_mult_" + str(mult) + "_strand_" +
                                                   strand + "\t1\t" + strand + "\t" + start + "\t" + end + "\t0\t1\t" + str(int(end) - int(start)) + "\t0\n")

            except IOError:
                print("Error: failed top open file " +
                      this_hints_file + " for writing!")
            bb_file = args.outdir + "/" + args.short_label + \
                "/" + args.short_label + "/"  + h_src + "_" + h_type + ".hints.bb"
            if (h_type == 'CDS') or (h_type == 'CDSpart') or (h_type == 'exon'):
                gp_file = this_hints_file + ".gp"
                bed_file = this_hints_file + ".bed"
                info_out_file = this_hints_file + ".infoOut.txt"
                gtf2bb(this_hints_file, gp_file, bed_file,
                       bb_file, info_out_file, ChromSizes_file)
            else:
                this_sorted_hints_file = this_hints_file + ".sorted"
                sort_bed3(this_hints_file, this_sorted_hints_file)
                bed2bigBed(12, this_sorted_hints_file,
                           ChromSizes_file, bb_file)
            info_to_trackDB(trackDb_file, h_src + "_" + h_type + ".hints",
                            "Hints of type " + h_type + " from source " + h_src, rgb_cols[col_idx], "hints", 12)
            col_idx = col_idx + 1
            if col_idx > 25:
                col_idx = 0

    print('Done')

''' Delete temporary directory '''

print("Deleting temporary files...")
#shutil.rmtree(args.outdir + "/" + args.short_label + "/tmp/")
print("Done.")

print('\nHub is ready, please copy to a web server, e.g.')
print('\"scp -r ' + args.outdir + "/" + args.short_label + " user@server:/target/location\"")
print('Feed the link to hub.txt to UCSC genome browser:')
print('\tMy Data -> Track Hubs -> MyHubs')
