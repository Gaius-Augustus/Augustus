#!/usr/bin/env python3

# Author: Anica Hoppe
# Last modified: March 25th 2020

# Given a list of transcript IDs, this python script identifies the
# boundaries of corresponding genes in a GTF or GFF3 file, makes
# new gene predictions in these regions with AUGUSTUS using mea and
# replaces genes with in frame stop codons in the GTF or GFF3 file
# by newly predicted genes.

import os
import os.path
import shutil
import subprocess
from inspect import currentframe, getframeinfo
import logging
import random
import string
import multiprocessing
from functools import partial

try:
    import argparse
except ImportError:
    raise ImportError(
        'Failed to import argparse. ' +
        'Try installing with \"pip3 install argparse\"')

try:
    import re
except ImportError:
    raise ImportError(
        'Failed to import re. Try installing with \"pip3 install re\"')

try:
    from Bio.Seq import Seq
    from Bio import SeqIO
except ImportError:
    frameinfo = getframeinfo(currentframe())
    raise ImportError('In file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      'Failed to import biophython modules. ' +
                      'Try installing with \"pip3 install biopython\"')

parser = argparse.ArgumentParser(description='Replaces genes with in-frame \
        stop codons (IFS) in a GTF or GFF3 file by genes without IFS \
        that are newly predicted with AUGUSTUS using mea.')
parser.add_argument('-g', '--genome', required=True, type=str,
                    help='genome sequence file (FASTA format)')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-t', '--gtf', type=str, help='GTF input file')
group.add_argument('-3', '--gff3', type=str, help='GFF3 input file')
parser.add_argument('-b', '--badGenes', required=True, type=str,
                    help='File with list of transcript IDs of genes with \
                    in-frame stop codons')
parser.add_argument('-H', '--hintsfile', required=False, type=str,
                    help='File with hints in gff format used for \
                    (re-)predicting genes with AUGUSTUS')
parser.add_argument('-o', '--out', required=True, type=str,
                    help="Name stem of the output file; will be extended with \
                    .gtf or .gff3 depending on the input format.")
parser.add_argument('-s', '--species', required=True, type=str,
                    help='Set the species to be used for running AUGUSTUS')
parser.add_argument('-e', '--extrinsicCfgFile', required=False, type=str,
                    help='Set extrinsic config file for AUGUSTUS')
parser.add_argument('-m', '--softmasking', required=False,
                    choices=['on', 'off'], default='off', type=str,
                    help='Choose \'on\' if the genome file is softmasked')
parser.add_argument('-u', '--UTR', required=False, choices=['on', 'off'],
                    default='off', type=str,
                    help='Predict the untranslated regions in addition \
                    to the coding sequence. If UTR=on was used in the original \
                    AUGUSTUS run, use \'--UTR on\' here, otherwise not')
parser.add_argument('-U', '--print_utr', required=False, choices=['on', 'off'],
                    default='off', type=str,
                    help='Choose \'on\' if --print-utr=on was used in the \
                    original AUGUSTUS run')
parser.add_argument('--additional_aug_args', required=False, type=str, 
                    help='One or several command line arguments to be passed \
                    to AUGUSTUS (which can not be given with another specific \
                    command line argument here). The list of arguments has to \
                    be given in quotes. If several arguments are given, \
                    they have to be separated by whitespace, i.e. \
                    \"--first_arg=sth --second_arg=sth\". If only one argument is \
                    given, the argument still has to contain a whitespace, i.e. \
                    \"--first_arg=sth \". Beware: Do not choose \
                    --alternatives-from-evidence=true because mea can not use \
                    this. Also do not set --exonnames (this parameter will be \
                    set automatically depending on the input gtf/gff3 file).')
parser.add_argument('-a', '--augustus_config_path', required=False, type=str,
                    help='Set path to the config directory of AUGUSTUS. \
                    If not given, will try to set augustus_config_path to \
                    environment variable AUGUSTUS_CONFIG_PATH. If this does not \
                    work, will try to set augustus_config_path to \
                    augustus_scripts_path/../config/. \
                    The commandline argument --AUGUSTUS_CONFIG_PATH has higher \
                    priority than the environment variable with the same name.')
parser.add_argument('-A', '--augustus_bin_path', required=False, type=str,
                    help='Set path to the AUGUSTUS directory that contains \
                    augustus binary. If not given, will try to locate the path \
                    with which(augustus)')
parser.add_argument('-S', '--augustus_scripts_path', required=False, type=str,
                    help='Set path to the AUGUSTUS scripts directory. If not \
                    given, will try to locate the path with which(gtf2gff.pl). \
                    If this does not work, will try to set the path relative \
                    to the augustus_bin_path (augustus_bin_path/../scripts/).')
parser.add_argument('-n', '--noCleanUp', required=False, action='store_true',
                    help='Unless chosen, temporary files created while running \
                    this script will be deleted at the end')
parser.add_argument('-p', '--print_format_examples', required=False,
                    action='store_true', help="Print gtf/gff3 input format \
                    examples, do not perform analysis")
parser.add_argument('-C', '--cdbtools_path', required=False, type=str, 
                    help = "Set path to cdbfasta/cdbyank. If not given, \
                    will try to locate the path with which(cdbfasta).")
parser.add_argument('-c', '--cores', required=False, default=1, type=int, 
                    help = "Set the number of cores used. Default will be 1.")
args = parser.parse_args()

### As args.hintsfile and args.extrinsicCfgFile have to be given together: ###
### Check if only one of args.hintsfile and args.extrinsicCfgFile is given ###
if (args.hintsfile is not None) ^ (args.extrinsicCfgFile is not None):
    parser.error("--hintsfile and --extrinsicCfgFile must be given together")

if args.print_format_examples:
    print('This script requires an annotation of protein coding genes ' +
          ' either in GTF or GFF3 format, containing genes with in-frame ' +
          'stop codon (IFS). The file must contain a gene line and a ' +
          'transcript/mRNA line. \n')
    print('GTF format example:\n')
    print('ctg1\tAUGUSTUS\tgene\t7816\t12187\t0.25\t+\t.\tg1\n' +
          'ctg1\tAUGUSTUS\ttranscript\t7816\t12187\t0.25\t+\t.\tg1.t1\n' +
          'ctg1\tAUGUSTUS\ttss\t7816\t7816\t.\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\t5\'-UTR\t7816\t8109\t0.63\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\texon\t7816\t8143\t.\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tstart_codon\t8110\t8112\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tCDS\t8110\t8143\t0.71\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tintron\t8144\t11630\t0.74\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tCDS\t11631\t11896\t0.89\t+\t2\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\texon\t11631\t11896\t.\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tintron\t11897\t12070\t0.8\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tCDS\t12071\t12139\t0.8\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\texon\t12071\t12187\t.\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\tstop_codon\t12137\t12139\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\t3\'-UTR\t12140\t12187\t0.47\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n' +
          'ctg1\tAUGUSTUS\ttts\t12187\t12187\t.\t+\t.\ttranscript_id \"g1.t1\"; gene_id \"g1\";)\n')
    print('GFF3 format example:\n')
    print('ctg1\tAUGUSTUS\tgene\t7816\t12187\t0.25\t+\t.\tID=g1\n' +
          'ctg1\tAUGUSTUS\ttranscript\t7816\t12187\t0.25\t+\t.\tID=g1.t1;Parent=g1\n' +
          'ctg1\tAUGUSTUS\ttranscription_start_site\t7816\t7816\t.\t+\t.\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tfive_prime_utr\t7816\t8109\t0.63\t+\t.\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tstart_codon\t8110\t8112\t.\t+\t0\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tintron\t8144\t11630\t0.74\t+\t.\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tintron\t11897\t12070\t0.8\t+\t.\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tCDS\t8110\t8143\t0.71\t+\t0\tID=g1.t1.cds;Parent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tCDS\t11631\t11896\t0.89\t+\t2\tID=g1.t1.cds;Parent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tCDS\t12071\t12139\t0.8\t+\t0\tID=g1.t1.cds;Parent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tstop_codon\t12137\t12139\t.\t+\t0\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\tthree_prime_utr\t12140\t12187\t0.47\t+\t.\tParent=g1.t1\n' +
          'ctg1\tAUGUSTUS\t\ttranscription_end_site\t12187\t12187\t.\t+\t.\tParent=g1.t1\n')
    print('Further more, this script requires a list of IDs of IFS affected transcripts. This list must be ' +
          'a file with transcript ids in the first column, that may also contain more ' +
          'columns and may contain lines beginning with \'#\' that will not be considered.\n')
    print('Examples:')
    print('1.\n')
    print('# tx_id\tstop_in_amino_acid\tstop_in_mRNA_start\tstop_in_mRNA_end\tgenomic_transcript_start\tgenomic_transcript_end\tstrand\n' +
          'g11.t1\t10\t28\t30\t376803\t386044\t+\n' +
          'g772.t1\t410\t1228\t1230\t731145\t752055\t+\n')
    print('or 2.\n')
    print('g11.t1\ng722.t1')
    exit(0)


''' ******************* BEGIN FUNCTIONS *************************************'''


def create_random_string():
    """ Funtion that creates a random string added to the logfile name 
        and tmp dir name """
    letters = string.ascii_lowercase
    randomString = ''.join(random.choice(letters) for i in range(8))
    tmp = "tmp_" + randomString + "/"
    log = "fix_IFS_log_" + randomString
    # if directory  or log_file exists, create a new random string
    while(os.path.exists(tmp) or os.path.exists(log)):
        randomString = ''.join(random.choice(letters) for i in range(8))
        tmp = "tmp" + randomString + "/"
        log = "fix_IFS_log_" + randomString
    return(randomString)


def create_log_file_name(randomString):
    """ Function that creates a log file with a random name """
    log = "fix_IFS_log_" + randomString
    return(log)


def create_tmp_dir(randomString):
    """ Function that creates a directory for temporary files with a random name """
    tmp = "tmp_" + randomString + "/"
    os.mkdir(tmp)
    logger.info("Creating directory " + tmp + ".")
    return(tmp)


def find_tool(toolname):
    """ Function that tries to locate a tool with which """
    if shutil.which(toolname) is not None:
        toolbinary = shutil.which(toolname)
    else:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': '
                + "Unable to locate binary of " + toolname + "!")
        exit(1)
    return(toolbinary)


def check_tool_in_given_path(given_path, toolname):
    """ If TOOL_PATH is provided as command line option, check whether the
    required binary is executable in that directory; return tool binary 
    with full path. """
    if not re.search(r"/$", given_path):
        given_path = given_path + "/"
    toolbinary = given_path + toolname
    toolbinary = os.path.abspath(toolbinary)
    if not(os.access(toolbinary, os.X_OK)):
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + toolbinary + 
                    " is not executable!")
        exit(1)
    return toolbinary


def run_process(args_lst, prc_out, prc_err):
	''' Function that runs a subprocess with arguments and specified STDOUT and STDERR '''
	try:
		logger.info("Trying to execute the following command:")
		logger.info(" ".join(args_lst))
		result = subprocess.run(args_lst, stdout=prc_out, stderr=prc_err)
		logger.info("Suceeded in executing command.")
		if(result.returncode == 0):
			return(result)
		else:
			frameinfo = getframeinfo(currentframe())
			logger.info('Error in file ' + frameinfo.filename + ' at line ' +
				  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " + 
		        str(result.returncode) + str(result.args))
			quit(1)
	except subprocess.CalledProcessError as grepexc:
		frameinfo = getframeinfo(currentframe())
		print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
		print("Error code: ", grepexc.returncode, grepexc.output)
		quit(1)


def run_grep_process(args_lst, prc_out, prc_err):
	''' Function that runs a subprocess with arguments and specified STDOUT and STDERR '''
	try:
		logger.info("Trying to execute the following command:")
		logger.info(" ".join(args_lst))
		result = subprocess.run(args_lst, stdout=prc_out, stderr=prc_err)
		logger.info("Suceeded in executing command.")
		if(result.returncode == 0 or result.returncode == 1):
			return(result)
		else:
			frameinfo = getframeinfo(currentframe())
			logger.info('Error in file ' + frameinfo.filename + ' at line ' +
				  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " + 
		        str(result.returncode) + str(result.args))
			quit(1)
	except subprocess.CalledProcessError as grepexc:
		frameinfo = getframeinfo(currentframe())
		print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
		print("Error code: ", grepexc.returncode, grepexc.output)
		quit(1)


def run_process_stdinput(args_lst, prc_in):
	''' Function that runs a subprocess with arguments and input from STDIN '''
	try:
		logger.info("Trying to execute the following command with input from STDIN:")
		logger.info(" ".join(args_lst))
		result = subprocess.run(args_lst, stdin=prc_in, stderr=subprocess.PIPE)
		logger.info("Suceeded in executing command.")
		if(result.returncode == 0):
			return(result)
		else:
			frameinfo = getframeinfo(currentframe())
			logger.info('Error in file ' + frameinfo.filename + ' at line ' +
				  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " + 
		        str(result.returncode) + str(result.args))
			quit(1)
	except subprocess.CalledProcessError as grepexc:
		frameinfo = getframeinfo(currentframe())
		print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
		print("Error code: ", grepexc.returncode, grepexc.output)
		quit(1)

def run_augustus_process(lst_entry):
    """ Function that runs an AUGUSTUS subprocess """
    subprcs_args = lst_entry[0]
    augustus_out = lst_entry[1]
    augustus_err = lst_entry[2]
    try:
        with open(augustus_out, "w") as augustus_out_handle:
            try:
                with open(augustus_err, "w") as augustus_err_handle:
                    run_process(subprcs_args, augustus_out_handle, augustus_err_handle)
            except IOError:
                frameinfo = getframeinfo(currentframe())
                logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                            str(frameinfo.lineno) + ': ' + "Could not open file " +
                            augustus_err + " for writing!")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Could not open file " +
                    augustus_out + " for writing!")

def sortFirst(val):
    """  """
    return val[0]

def sortSecond(val):
    """  """
    return val[1]

''' ******************* END FUNCTIONS *************************************'''

### Check whether sufficient options have been provided ###

### Create log file and tmp directory for saving files that can be removed afterwards ###
rString = create_random_string()
log = create_log_file_name(rString)
logger = logging.getLogger("")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(log)
fh.setLevel(logging.INFO)
formatter = logging.Formatter("%(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)
tmp = create_tmp_dir(rString)

### Find required binaries on system ###
grep = find_tool("grep")
perl = find_tool("perl")
if args.cdbtools_path:
    cdbfasta = check_tool_in_given_path(args.cdbtools_path, "cdbfasta")
    cdbyank = check_tool_in_given_path(args.cdbtools_path, "cdbyank")
else:
    cdbfasta = find_tool("cdbfasta")
    cdbyank = find_tool("cdbyank")

### Find augustus_bin_path ###
if args.augustus_bin_path:
    augustus = check_tool_in_given_path(args.augustus_bin_path, "augustus")
else:
    augustus = find_tool("augustus")

### Find augustus_scripts_path ###
if args.augustus_scripts_path:
    gtf2gff = check_tool_in_given_path(args.augustus_scripts_path, "gtf2gff.pl")
elif shutil.which("gtf2gff.pl") is not None:
    gtf2gff = find_tool("gtf2gff.pl")
else:
    logger.info(
        "Trying to set augustus_scripts_path to augustus_bin_path/../scripts/")
    test_augustus_scripts_path = augustus_bin_path + "../scripts/"
    if os.path.existis(test_augustus_scripts_path):
        augustus_scripts_path = test_augustus_scripts_path
        gtf2gff = augustus_scripts_path + "gtf2gff.pl"
    else:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + 
                    "Unable to set augustus_scripts_path!")
        exit(1)

### Find augustus_config_path ###
augustus_config_path = ""
if args.augustus_config_path:
    augustus_config_path = args.augustus_config_path
else:
    logger.info("Trying to find environment variable " + 
                "AUGUSTUS_CONFIG_PATH")
    if os.environ.get('AUGUSTUS_CONFIG_PATH') is not None:
        test_augustus_config_path = os.environ.get('AUGUSTUS_CONFIG_PATH')
        if os.path.exists(test_augustus_config_path):
            augustus_config_path = test_augustus_config_path
            logger.info("Found environment variable AUGUSTUS_CONFIG_PATH " + 
                        "and set augustus_config_path to environment " +
                        "variable AUGUSTUS_CONFIG_PATH")
if augustus_config_path == "":
    logger.info(
        "Did not find environment variable AUGUSTUS_CONFIG_PATH " +
        "(either variable does not exist, or the path given in variable " +
        "does not exist). Will try to set this variable in a different " +
        "way. \n" +
        "Trying to set augustus_config_path to " +
        "augustus_scripts_path/../config/")
    augustus_scripts_path = gtf2gff
    augustus_scripts_path = re.sub(r'gtf2gff\.pl', '', augustus_scripts_path)
    test_augustus_config_path = augustus_scripts_path + "../config/"
    if os.path.exists(test_augustus_config_path):
        augustus_config_path = test_augustus_config_path
    else:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Unable to set " +
                    "augustus_config_path!")
        exit(1)

### Save gtf file and create output file name ###
if args.gtf is not None:
    gtf_file = args.gtf # setting gtf_file because identically named variable
                        # is also used if a gff3 file was originally
                        # provided
    out_file = args.out + ".gtf"
else:
    gtf_file = tmp + "converted_augustus.gtf"
    out_file = args.out + ".gff3"


### If input file is in gff3 format, convert input file to gtf format ###
if args.gff3 is not None:
    g_id = ""
    tx_id = ""
    try:
        with open(args.gff3, "r") as gff3_handle:
            try:
                with open(gtf_file, "w") as gtf_handle:
                    for line in gff3_handle:
                        if not re.match(r"#.+", line):
                            if re.match(r".+\tgene\t.+", line):
                                match = re.match(r"(.+)ID=(\S+)", line)
                                g_id = match.group(2)
                                gtf_handle.write(match.group(1) + g_id + "\n")
                                transcript_seen = False
                                mRNA_seen = False
                            elif re.match(r".+\ttranscript\t.+", line):
                                if not mRNA_seen:
                                    match = re.match(
                                        r"(.+)ID=(\S+);Parent=.+", line)
                                    tx_id = match.group(2)
                                    gtf_handle.write(
                                        match.group(1) + tx_id + "\n")
                                transcript_seen = True
                            elif re.match(r".+\tmRNA\t.+", line):
                                if not transcript_seen:
                                    match = re.match(
                                        r"(.+)mRNA(.+)ID=(\S+);Parent=.+", line)
                                    tx_id = match.group(3)
                                    gtf_handle.write(match.group(
                                        1) + "transcript" + match.group(2) + tx_id + "\n")
                                mRNA_seen = True
                            elif re.match(r".+\ttranscription_start_site\t.+", line):
                                match = re.match(
                                    r"(.+)transcription_start_site(\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                                gtf_handle.write(match.group(1) + "tss" + match.group(
                                    2) + "transcript_id \"" + tx_id + "\"; gene_id \"" + g_id + "\";\n")
                            elif re.match(r".+\ttranscription_end_site\t.+", line):
                                match = re.match(
                                    r"(.+)transcription_end_site(\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                                gtf_handle.write(match.group(1) + "tts" + match.group(
                                    2) + "transcript_id \"" + tx_id + "\"; gene_id \"" + g_id + "\";\n")
                            elif re.match(r".+\tfive_prime_utr\t.+", line):
                                match = re.match(
                                    r"(.+)five_prime_utr(\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                                gtf_handle.write(match.group(1) + "5\'-UTR" + match.group(
                                    2) + "transcript_id \"" + tx_id + "\"; gene_id \"" + g_id + "\";\n")
                            elif re.match(r".+\tthree_prime_utr\t.+", line):
                                match = re.match(
                                    r"(.+)three_prime_utr(\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                                gtf_handle.write(match.group(1) + "3\'-UTR" + match.group(
                                    2) + "transcript_id \"" + tx_id + "\"; gene_id \"" + g_id + "\";\n")
                            elif re.match(r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line):
                                match = re.match(
                                    r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                                gtf_handle.write(match.group(
                                    1) + "transcript_id \"" + tx_id + "\"; gene_id \"" + g_id + "\";\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                             str(frameinfo.lineno) + ': ' + "Could not open file " +
                             gtf_file + " for writing!")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Could not open file " +
                    args.gff3 + " for reading!")


### Read file with list of transcript IDs of genes with IFS and extract transcript IDs ###
bad_tx = []
try:
    with open(args.badGenes, "r") as bad_handle:
        for line in bad_handle:
            if not re.match(r"#.*", line):
                match = re.search(r"(\S+).*", line)
                tx_id = match.group(1)
                bad_tx.append(tx_id)
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + "Could not open file " +
                args.badGenes + " for reading!")


### Read gtf file and extract seq id, gene start and end ###
genes = {}
genesStart_ofSeqID = {}
genesEnd_ofSeqID = {}
g_id_of_tx = {}
g_id = ""
exonnames = "off"
try:
    with open(gtf_file, "r") as gtf_handle:
        for line in gtf_handle:
            if re.match(r".+\tgene\t.+", line):
                seq_id, st, en, stx, g_id = re.match(
                    r"(\S+)\t\S+\tgene\t(\d+)\t(\d+)\t\S+\t(\S)\t\S+\t+[ID=]*(\S+)", line).groups()
                genes[g_id] = {}
                genes[g_id] = {'seq_id': seq_id, 'start': int(
                    st), 'end': int(en), 'strand': stx}
                if not(seq_id in genesStart_ofSeqID):
                    genesStart_ofSeqID[seq_id] = []
                genesStart_ofSeqID[seq_id].append((g_id, int(st)))
                if not(seq_id in genesEnd_ofSeqID):
                    genesEnd_ofSeqID[seq_id] = []
                genesEnd_ofSeqID[seq_id].append((g_id, int(en)))
            elif re.match(r".+\ttranscript\t.+", line):
                match = re.match(
                    r".+transcript\t\d+\t\d+\t\S+\t\S+\t\S+\t[ID=]*(\S+)[;Parent=\S+]*", line)
                tx_id = match.group(1)
                g_id_of_tx[tx_id] = g_id
            elif re.match(r".+internal.+", line):
                exonnames = "on"
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + "Could not open file " +
                gtf_file + " for reading!")

### Create list of "bad genes" = genes with in-frame stop codon ###
bad_genes = []
for tx_id in bad_tx:
    bad_genes.append(g_id_of_tx[tx_id])


### Create list of scaffold lengths ###
seq_len = {}
try:
    with open(args.genome, "r") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            seq_len[record.id] = len(record.seq)
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + "Could not open file " +
                args.genome + " for reading!")


### Compute positions for prediction start and prediction end ###
regions = {}
genesStart = []
genesEnd = []
for g_id in bad_genes:
    start = genes[g_id]['start']
    end = genes[g_id]['end']
    seq_id = genes[g_id]['seq_id']
    genesStart = genesStart_ofSeqID[seq_id]
    genesStart.sort(key=sortSecond)
    geneIndexStartList = genesStart.index((g_id, start))
    genesEnd = genesEnd_ofSeqID[seq_id]
    genesEnd.sort(key=sortSecond)
    geneIndexEndList = genesEnd.index((g_id, end))
    # compute start position of the region between the neighbour genes
    pgIndexEndList = geneIndexEndList
    while genesEnd[pgIndexEndList][1] > start and pgIndexEndList != 0:
        pgIndexEndList = pgIndexEndList - 1
    if genesEnd[pgIndexEndList][1] < start:
        pg_id = genesEnd[pgIndexEndList][0]
        region_start = genes[pg_id]['end'] + 1
    else:
        region_start = 1
    if (start - region_start) > 50000:
        region_start = (start - 50000)
    # compute end position of the region between the neighbour genes
    ngIndexStartList = geneIndexStartList
    while genesStart[ngIndexStartList][1] < end and ngIndexStartList != (len(genesStart) - 1):
        ngIndexStartList = ngIndexStartList + 1
    if genesStart[ngIndexStartList][1] > end:
        ng_id = genesStart[ngIndexStartList][0]
        region_end = genes[ng_id]['start'] - 1
    else:
        region_end = seq_len[genes[g_id]['seq_id']]
    if (region_end - end) > 50000:
        region_end = (end + 50000)
    # add region to list regions
    if not(seq_id in regions):
        regions[seq_id] = []
    if not((region_start, region_end) in regions[seq_id]):
        regions[seq_id].append((region_start, region_end))

for seq_id in regions:
    regions[seq_id].sort(key=sortFirst)


### Generate FASTA files for each scaffold containing a gene with in-frame stop codon ###
# Run cdbfasta
genome_cidx = args.genome + ".cidx"
try:
    with open(genome_cidx, "w") as cidx_handle:
        subprcs_args = [cdbfasta, args.genome]
        run_process(subprcs_args, cidx_handle, subprocess.PIPE)
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                 str(frameinfo.lineno) + ': ' + "Could not open file " +
                 genome_cidx + " for writing!")

# Run cdbyank
for seq_id in regions:
    fasta_file = tmp + "genome." + seq_id + ".fa"
    try:
        with open(fasta_file, "w") as fasta_tmp_handle:
            subprcs_args = [cdbyank, "-a", seq_id, genome_cidx]
            run_process(subprcs_args, fasta_tmp_handle, subprocess.PIPE)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Could not open file " +
                    fasta_file + " for writing!")


### If a hintsfile is given, create hintsfiles for each scaffold ### 
### containing genes with in-frame stop codon ###
if args.hintsfile is not None:
    for seq_id in regions:
        new_hintsfile = tmp + "hints." + seq_id + ".gff"
        try:
            with open(new_hintsfile, "w") as new_hints_handle:
                subprcs_args = [grep, "-w", seq_id, args.hintsfile]
                run_grep_process(subprcs_args, new_hints_handle, subprocess.PIPE)
        except IOError:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                        str(frameinfo.lineno) + ': ' + "Could not open file " +
                        new_hintsfile + " for writing!")


### Run AUGUSTUS for each gene with in-frame stop codon ###
### (with predictionstart=region_start and predictionend=region_end) ###
# Directory for AUGUSTUS output:
out = tmp + "out/"
if not os.path.exists(out):
    os.mkdir(out)
else:
    logger.info("Directory out already exists.")

augustus_lst = []
for seq_id in regions:
    i = 1
    genome = tmp + "genome." + seq_id + ".fa"
    if args.hintsfile is not None:
        hintsfile = tmp + "hints." + seq_id + ".gff"
    for reg in regions[seq_id]:
        augustus_out = out + "augustus." + seq_id + "." + str(i) + ".out"
        augustus_err = out + "augustus." + seq_id + "." + str(i) + ".err"
        start = str(reg[0])
        end = str(reg[1])

        subprcs_args = [augustus, "--mea=1", "--species="+args.species,
                        "--softmasking="+args.softmasking, "--UTR="+args.UTR, 
                        "--print_utr="+args.print_utr, "--genemodel=complete",
                        "--alternatives-from-evidence=0", "--exonnames="+exonnames,
                        "--predictionStart="+start, "--predictionEnd="+end,
                        "--AUGUSTUS_CONFIG_PATH="+augustus_config_path, genome]
        if args.hintsfile is not None:
            subprcs_args.append("--hintsfile="+hintsfile)
            subprcs_args.append("--extrinsicCfgFile="+args.extrinsicCfgFile)
        if args.additional_aug_args is not None:
            additional_aug_args = args.additional_aug_args.split()
            subprcs_args.extend(additional_aug_args)

        augustus_lst_entry = (subprcs_args, augustus_out, augustus_err)
        augustus_lst.append(augustus_lst_entry)
        
        i += 1

if __name__ == "__main__":
    with multiprocessing.Pool(processes = args.cores) as pool:
        pool.map(run_augustus_process, augustus_lst)


### Replace genes with in-frame stop codon by new genes  ###
augustus_tmp_file = tmp + "augustus.tmp"
gene_id = ""
gene_num = 0
new_gene_id = ""
seq_id = ""
ifs_gene_num = 0
tx_num = 0
try:
    with open(gtf_file, "r") as gtf_handle:
        try:
            with open(augustus_tmp_file, "w") as tmp_handle:
                for line in gtf_handle:
                    if re.match(r".+\tgene\t.+", line):
                        match = re.match(
                            r"(\S+)(\t\S+\tgene\t)(\d+)(\t)(\d+)(\t\S+\t\S+\t\S+\t)(\S+)", line)
                        if match.group(1) != seq_id:
                            ifs_gene_num = 0
                        seq_id = match.group(1)
                        gene_id = match.group(7)
                        if not gene_id in bad_genes:
                            gene_num += 1
                            new_gene_id = "g" + str(gene_num)
                            tx_num = 0
                            tmp_handle.write(match.group(1) + match.group(2) + match.group(
                                3) + match.group(4) + match.group(5) + match.group(6) + new_gene_id + "\n")
                        else:
                            ifs_gene_num += 1
                            ifs_gene_file = out + "augustus." + \
                                seq_id + "." + str(ifs_gene_num) + ".out"
                            try:
                                with open(ifs_gene_file, "r") as ifs_handle:
                                    for ifs_line in ifs_handle:
                                        if re.match(r".+\tgene\t.+", ifs_line):
                                            m = re.match(
                                                r"(\S+\t\S+\tgene\t\d+\t\d+\t\S+\t\S+\t\S+\t)(\S+)", ifs_line)
                                            gene_num += 1
                                            new_gene_id = "g" + str(gene_num)
                                            tx_num = 0
                                            tmp_handle.write(
                                                m.group(1) + new_gene_id + "\n")
                                        elif re.match(r".+\ttranscript\t.+", ifs_line):
                                            tx_num += 1
                                            new_tx_id = new_gene_id + \
                                                ".t" + str(tx_num)
                                            m = re.match(
                                                r"(\S+\t\S+\ttranscript\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", ifs_line)
                                            tmp_handle.write(
                                                m.group(1) + new_tx_id + "\n")
                                        elif re.match(r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t).+", ifs_line):
                                            m = re.match(
                                                r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t).+", ifs_line)
                                            tmp_handle.write(m.group(
                                                1) + "transcript_id \"" + new_tx_id + "\"; gene_id \"" + 
                                                new_gene_id + "\";\n")
                            except IOError:
                                frameinfo = getframeinfo(currentframe())
                                logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                                            str(frameinfo.lineno) + ': ' + "Could not open file " +
                                            ifs_gene_file + " for reading!")
                    elif re.match(r".+\ttranscript\t.+", line):
                        tx_num += 1
                        new_tx_id = new_gene_id + ".t" + str(tx_num)
                        match = re.match(
                            r"(\S+\t\S+\ttranscript\t\d+\t\d+\t\S+\t\S+\t\S+\t)\S+", line)
                        if not gene_id in bad_genes:
                            tmp_handle.write(match.group(1) + new_tx_id + "\n")
                    elif re.match(r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t).+", line):
                        match = re.match(
                            r"(\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\t).+", line)
                        if not gene_id in bad_genes:
                            tmp_handle.write(match.group(
                                1) + "transcript_id \"" + new_tx_id + "\"; gene_id \"" + new_gene_id + "\";\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                        str(frameinfo.lineno) + ': ' + "Could not open file " +
                        augustus_tmp_file + " for writing!")

except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + "Could not open file " +
                gtf_file + " for reading!")

### Sort gtf file or convert gtf file to gff3 format (sorted) depending on the input format ###
try:
    with open(augustus_tmp_file, "r") as augustus_tmp_handle:
        subprcs_args = [perl, gtf2gff, "--out="+out_file, "--printExon"]
        if args.gff3 is not None:
            subprcs_args.append("--gff3")
        run_process_stdinput(subprcs_args, augustus_tmp_handle)
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + "Could not open file " +
                augustus_tmp_file + " for reading!")


### If not 'noCleanUp' is chosen, remove tmp directory ###
if not args.noCleanUp:
    shutil.rmtree(tmp)
