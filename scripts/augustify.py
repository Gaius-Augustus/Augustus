#!/usr/bin/env python3

import os
import errno
import os.path
import argparse
import re
import shutil
import subprocess
import logging
import random
import string
import json
import numpy as np
from difflib import SequenceMatcher
from inspect import currentframe, getframeinfo
from itertools import compress
import multiprocessing
from datetime import datetime
import uuid
import ast
import glob



__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2022. All rights reserved."
__credits__ = "Mario Stanke, Anica Hoppe, Marnix Medema"
__license__ = "Artistic License"
__version__ = "1.0.3"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

parser = argparse.ArgumentParser(
    description='augustify.py will either assign a majority vote Augustus ' +
                'parameter set to a given set of input sequences, or it ' +
                'will sort input sequences according to their respective ' +
                'best parameter sets and perform ab initio gene prediction ' +
                'with these parameter sets.')
parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-p', '--parameter_list', required=True, type=str,
                    help='File that lists parameter sets to be tried on ' +
                    'input sequence (one parameter set name per line)')
parser.add_argument('-m', '--metagenomic_classification_outfile', type=str,
                    help='Output a tabulator separated text file that ' +
                    'assigns sequences to parameter sets (last column ' +
                    'contains probability).')
parser.add_argument('-P', '--prediction_file', type=str,
                    help='GFF file with gene predictions (only compatible ' +
                    'with option --metagenomic_classification_outfile/-m).')
parser.add_argument('-s', '--species', action='store_true',
                    help='Output the dominant and most suitable ' +
                    'parameter set name across all input sequences.')
parser.add_argument('-a', '--augustus_config_path', required=False, type=str,
                    help='Set path to the config directory of AUGUSTUS. \
                    If not given, will try to set augustus_config_path to \
                    environment variable AUGUSTUS_CONFIG_PATH. If this does not \
                    work, will try to set augustus_config_path to \
                    augustus_bin_path/../config/. \
                    The command-line argument --augustus_config_path has higher \
                    priority than the environment variable with the same name.')
parser.add_argument('-A', '--augustus_bin_path', required=False, type=str,
                    help='Set path to the AUGUSTUS directory that contains \
                    augustus binary. If not given, will try to locate the path \
                    with which(augustus)')
parser.add_argument('-t', '--threads', required=False, type=int, default=1,
                    help='Number of threads for running augustus. The number ' +
                    'of threads should not be greater than the number of ' +
                    'species parameter sets.')
parser.add_argument('-c', '--use_coordinates', required=False, type=str, #default='selected_coordinates.txt',
                  help='This is the file with the selected coordinates for contigs. that Augustify will use for predictions.' +
                  'Augustify currently takes the n nucleotides of a contig RANDOMLY to pick a parameter set.' +
                  'Prediction is performed on the full set of sequence after picking.')
args = parser.parse_args()

if ( (args.metagenomic_classification_outfile) and (args.species) ):
    print("Incompatible options selected: you must either specify " +
          "--metagenomic/-m OR --species/-s. You cannot run " +
          "augustify.py with both options!")
    exit(1)
elif ( (args.metagenomic_classification_outfile) and (args.species) ):
    print("Missing opion: you must either specify " +
          "--metagenomic/-m OR --species/-s.")

if ( (args.prediction_file is True) and (args.metagenomic_classification_outfile is not True)):
    print("Incompatible options selected: " +
          "--prediction_file/-P requires argument " +
          "--metagenomic_classification_outfile/-m!")
    exit(1)

if args.metagenomic_classification_outfile:
    if os.path.isfile(args.metagenomic_classification_outfile):
        print('File ' + args.metagenomic_classification_outfile +
            ' already exists. Please specify different file or delete file.')
        exit(1)

if args.prediction_file:
    if os.path.isfile(args.prediction_file):
        print('File ' + args.prediction_file +
            ' already exists. Please specify different file or delete file.')
        exit(1)

if args.use_coordinates:
    if not os.path.isfile(args.use_coordinates):
        print('File ' + args.use_coordinates +
            ' does not exists. Please specify different file correctly.')
        exit(1)

''' ******************* BEGIN FUNCTIONS *************************************'''

def normalize_float(manti, expo):
    """ Function that normalizes a float consisting of mantisse and exponent
        to a mantisse that has one digit in front of comma. Assumption is 
        that input manisse is either 0.xxxx or 1.xxx, i.e. the number of 
        digits in front of comma is already one in all cases."""
    if manti < 1:
        manti = manti*10
        expo = expo + 1
    return manti,expo


def create_random_string():
    """ Function that creates a random string added to the logfile name
        and tmp dir name """
    letters = string.ascii_lowercase
    randomString = ''.join(random.choice(letters) for i in range(8))
    tmp = "tmp_" + randomString + "/"
    log = "augustify_log_" + randomString
    # if directory or log_file exists, create a new random string
    while(os.path.exists(tmp) or os.path.exists(log)):
        randomString = ''.join(random.choice(letters) for i in range(8))
        tmp = "tmp" + randomString + "/"
        log = "augustify_log_" + randomString
    return(randomString)


def create_log_file_name(randomString):
    """ Function that creates a log file with a random name """
    log = "augustify_log_" + randomString
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
    return toolbinary


def check_tool_in_given_path(given_path, toolname):
    """ If TOOL_PATH is provided as command-line option, check whether the
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
                    ' is not executable!')
        exit(1)
    return toolbinary


def run_process(args_lst, prc_out, prc_err):
    ''' Function that runs a subprocess with arguments and specified STDOUT and STDERR '''
    try:
        logger.info("Trying to execute the following command:")
        logger.info(" ".join(args_lst))
        result = subprocess.run(args_lst, stdout=prc_out, stderr=prc_err)
        logger.info("Succeeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + 'Return code of subprocess was ' +
                str(result.returncode) + str(result.args))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + 'Failed executing: ',
              " ".join(grepexc.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)

def run_process_stdinput(args_lst, prc_in):
    ''' Function that runs a subprocess with arguments and input from STDIN '''
    try:
        logger.info("Trying to execute the following command with input from STDIN:")
        logger.info(" ".join(args_lst))
        result = subprocess.run(args_lst, stdin=prc_in, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logger.info("Succeeded in executing command.")
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
              str(frameinfo.lineno) + ': ' + 'Failed executing: ',
              " ".join(grepexc.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


def work_augustus(cmd_ext_lst):
    ''' Function that keeps running AUGUSTUS to compute emiprobs until it
    returns a result or until the entire sequence has been tried; walk
    the input in steps 4000 nt windows, overlapping 200 nt, if possible;
    we do this because AUGUSTUS sometimes crashes when computing the desired
    probability; problem should ultimately be fixed in AUGUSTUS '''

    #global seqs, starts

    segmlen = 15000  # These parameters need to be evaluated, later!
    stepwidth = 7500  # These parameters need to be evaluated, later!
    cmd1 = cmd_ext_lst[0:7]
    cmd2 = cmd_ext_lst[7:9]
    sub = re.search(r"seq(.+?)\.", str(cmd2[0]))[0][:-1]
    seqlen = cmd_ext_lst[-1]

    if seqlen <= segmlen:
        cmd = cmd1 + cmd2
        print("short sequence")
        return subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=False), sub, np.NaN
    else:
        if args.use_coordinates:
            with open(args.use_coordinates) as f:
                data = f.read()
                dict = ast.literal_eval(data)
                curr_start = dict.get(sub)
        else:
            # if the selected_coordinates.txt wasn't provided
            curr_start = random.randint(0, seqlen - segmlen)
        start = curr_start

        print("Randomly selected start coordinate for ", sub, "is: ", curr_start)

        curr_end = segmlen + curr_start
        returncode = 1
        while (curr_start < curr_end) and (curr_end - curr_start + 1 >= segmlen) and not (returncode == 0):
            cmd = cmd1 + ["--predictionStart=" + str(curr_start), "--predictionEnd=" + str(curr_end)] + cmd2
            result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, shell=False)
            returncode = result.returncode
            curr_start = curr_start + stepwidth
            if (curr_end + stepwidth) <= seqlen:
                curr_end = curr_end + stepwidth
            else:
                curr_end = seqlen
        print("run work_augustus for long seq")
        return result, sub, start

def augustify_seq(hindex, header, seqs, tmp, params, id):
    ''' Function that runs a subprocess with arguments and specified STDOUT and STDERR '''
    logger.info("Processing sequence: " + header)
    # store length for loop (emiprobs bug)
    currlen = len(seqs)
    # sequence files for prediction
    try:
        with open(tmp + "seq" + str(hindex) + ".fa", "w") as seq_handle:
            seq_handle.write(header)
            seq_handle.write(seqs)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Could not open file " +
                    tmp + "seq" + str(hindex) + ".fa" + " for writing!")

    # construct augustus calls
    calls = []
    for species in params:
        curr_call =  [augustus, '--AUGUSTUS_CONFIG_PATH=' + augustus_config_path, 
                     '--species=' + species.rstrip(), '--genemodel=complete', 
                     '--emiprobs=on',
                     '--softmasking=0', tmp + "seq" + str(hindex) + ".fa",
                     '--outfile=' + tmp + "seq" + str(hindex) + "_" + species.rstrip() + ".gff",
                     '--errfile=' + tmp + "seq" + str(hindex) + "_" + species.rstrip() + ".err", currlen]
        calls.append(curr_call)

# execute processes in parallel

    if __name__ == '__main__':

        if args.use_coordinates:

            with multiprocessing.Pool(processes=args.threads) as pool:
                results = pool.map(work_augustus, calls)

        else:
            time = datetime.now().strftime("%H%M%S")
            filename = f'selected_coordinates_{id}_{time}.txt'

            with open(filename, 'w') as file:
                with multiprocessing.Pool(processes=args.threads) as pool:
                    print("multiprocessing")
                    res = pool.map_async(work_augustus, calls)
                    res.wait()
                    results_dict = {}
                    for result in res.get():
                        print(f'Got result: {result}', flush=True)
                        results_dict[result[1]] = result[2]
                    print(res)
                    print("Get all_results list")
                file.write(json.dumps(results_dict))
            logger.info("Saved randomly selected coordinates for contigs to file.")

    logger.info("Finished parallel execution!")

    # parse results and find max prob for this sequence
    results = {}

    for species in params:
        try:
            with open(tmp + "seq" + str(hindex) + "_" + species.rstrip() + ".gff", "r") as spec_result_handle:
                for line in spec_result_handle:
                    thismatch1 = re.search(r'\# joint probability of gene structure and sequence in \S+ model: (\d+\.\d+)$', line)
                    thismatch2 = re.search(r'\# joint probability of gene structure and sequence in \S+ model: (\d+\.\d+)e(-\d+)', line)
                    thismatch3 = re.search(r'\# joint probability of gene structure and sequence in \S+ model: 0$', line)
                    thismatch4 = re.search(r'\# joint probability of gene structure and sequence in \S+ model: (\d+)e(-\d+)', line)
                    
                    if thismatch1:
                        mantisse,exponent = normalize_float(float(thismatch1.group(1)), int(0))
                        results[species.rstrip()] = {'mantisse' : mantisse,
                                                     'exponent' : exponent,
                                                     'original' : thismatch1.group(1)}
                    elif thismatch2:
                        mantisse,exponent = normalize_float(float(thismatch2.group(1)), int(thismatch2.group(2)))
                        results[species.rstrip()] = {'mantisse' : mantisse,
                                                     'exponent' : exponent,
                                                     'original' : thismatch2.group(1) + 'e' + thismatch2.group(2)}
                    elif thismatch3:
                        results[species.rstrip()] = {'mantisse' : float(0),
                                                     'exponent' : int(0),
                                                     'original' : 0}
                    elif thismatch4:
                        mantisse,exponent = normalize_float(float(thismatch4.group(1)), int(thismatch4.group(2)))
                        results[species.rstrip()] = {'mantisse' : mantisse,
                                                     'exponent' : exponent,
                                                     'original' : thismatch4.group(1) + 'e' + thismatch4.group(2)}
            
                if species.rstrip() not in results:
                    # AUGUSTUS may not predict any gene at all, hence no probability available
                    results[species.rstrip()] = {'mantisse' : float(0),
                                                 'exponent' : int(0),
                                                 'original' : 0}

        except IOError:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + "Could not open file " +
                    tmp + "seq" + str(hindex) + "_" + species.rstrip() + ".gff" + " for reading!")

        try:
            os.remove(tmp + "seq" + str(hindex) + '_' + species.rstrip() + '.gff')
        except OSError:
            pass

        try:
            os.remove(tmp + 'seq' + str(hindex) + '_' + species.rstrip() + '.err')
        except OSError:
            pass

    # determine maximum exponent
    max_exp = float('-inf')
    max_mant = 0.0
    max_species = "undef"
    for species in results:
        if not(results[species]['mantisse'] == 0):
            if (results[species]['exponent'] == max_exp) and (results[species]['mantisse'] > max_mant):
                max_mant = results[species]['mantisse']
                max_species = species
            elif results[species]['exponent'] > max_exp:
                max_exp = results[species]['exponent']
                max_mant = results[species]['mantisse']
                max_species = species

    thismatch = re.search(r'>(\S+)', header.rstrip())
    # write result if appropriate
    if args.metagenomic_classification_outfile:
        try:
            with open(args.metagenomic_classification_outfile, "a+") as classify_handle:
                if not (max_mant == 0):
                    classify_handle.write(thismatch.group(1) + "\t" + max_species +
                                        "\t" + results[max_species]['original'] + "\n")
                else: # this case should not occur anymore unless you use only completely unfitting parameter sets
                    classify_handle.write(thismatch.group(1) + "\t" + "undef" +
                                        "\t" + str(0) + "\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                        str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                        args.metagenomic_classification_outfile + ' for writing!')
            
    # run augustus with all gene models for the selected species
    if args.prediction_file and not(max_mant == 0):
        curr_call = [augustus, '--AUGUSTUS_CONFIG_PATH=' + augustus_config_path,
                     '--species=' + max_species,
                     '--softmasking=1', tmp + 'seq' + str(hindex) + '.fa',
                     '--outfile=' + tmp + 'seq' + str(hindex) + '_' + max_species + '_max.gff',
                     '--errfile=' + tmp + 'seq' + str(hindex) + '_' + max_species + '_max.err']
        aug_res = subprocess.run(curr_call, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, shell=False)
        try:
            with open(tmp + "seq" + str(hindex) + "_" + max_species + "_max.gff", "r") as aug_handle:
                try:
                    with open(tmp + "all_preds.gff", "a+") as gff_handle:
                        for line in aug_handle:
                            gff_handle.write(line)
                except IOError:
                    frameinfo = getframeinfo(currentframe())
                    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                                str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                                tmp + 'all_preds.gff' + ' for writing!')
        except IOError:
            frameinfo = getframeinfo(currentframe())
            logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                        str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                        tmp + 'seq' + str(hindex) + '_' + max_species + 
                        '_max.gff' + ' for reading!')

        try:
            os.remove(tmp + 'seq' + str(hindex) + '_' + max_species + '_max.gff')
        except OSError:
            pass

        try:
            os.remove(tmp + 'seq' + str(hindex) + '_' + max_species + '_max.err')
        except OSError:
            pass

        try:
            os.remove(tmp + 'seq' + str(hindex) + '.fa')
        except OSError:
            pass

    if max_mant == 0:
        return "undef"
    else:
        return max_species

''' ******************* END FUNCTIONS *************************************'''

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

### Check whether provided threads are available
count = len(os.sched_getaffinity(0))
if args.threads > count:
    frameinfo = getframeinfo(currentframe())
    logger.info('Warning in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + 'You provided ' +
                str(args.threads) + ' threads but your system has only ' +
                str(count) + ' cores. Decreasing to that number.')
    args.threads = count

### Find augustus_bin_path ###
if args.augustus_bin_path:
    augustus = check_tool_in_given_path(args.augustus_bin_path, "augustus")
else:
    augustus = find_tool("augustus")

if args.prediction_file:
    ### Find join_aug_preds.pl - used only if option prediction_file is set
    perl = find_tool("perl")
    if args.augustus_bin_path:
        join_aug_pred = check_tool_in_given_path(args.augustus_bin_path + "../scripts/", "join_aug_pred.pl")
    else:
        join_aug_pred = find_tool("join_aug_pred.pl")


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
        "augustus_bin_path/../config/")
    augustus_bin_path = augustus
    augustus_bin_path = re.sub(r'/augustus', '/', augustus_bin_path)
    test_augustus_config_path = augustus_bin_path + "../config/"
    if os.path.exists(test_augustus_config_path):
        augustus_config_path = test_augustus_config_path
    else:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + 'Unable to set ' +
                    'augustus_config_path!')
        exit(1)

### Read species parameter names
params = [];
try:
    with open(args.parameter_list, "r") as params_handle:
        for line in params_handle:
            params.append(line)
except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                args.parameter_list + ' for reading!')


### Check whether given parameter sets are valid
must_delete = [];
for species_set in params:
    species_set = species_set.rstrip()
    curr_spec_set_file = augustus_config_path +  '/species/' + species_set + '/' + species_set + '_exon_probs.pbl'
    # print("Checking for " + curr_spec_set_file)
    if os.path.isfile(curr_spec_set_file):
        must_delete.append(True)
    else:
        must_delete.append(False)
        frameinfo = getframeinfo(currentframe())
        logger.info('Warning in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + 'Species parameter set ' + species_set +
                    ' does not exist in ' + augustus_config_path +
                    '. Will ignore this parameter set!')
params = list(compress(params, must_delete))

### Check whether more than 1 set remains
if len(params) < 2:
    frameinfo = getframeinfo(currentframe())
    if len(params) == 0:
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' +
                'Remaining number of parameter sets is zero. ' +
                'Please run augustus outside of augustify.py with ' +
                'a parameter set!')
    else:
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                   str(frameinfo.lineno) + ': ' +
                    'Remaining number of parameter sets is one. ' +
                    'Please run augustus outside of augustify.py with ' +
                    'only one parameter set (here: ' + params[0] +')!')
    exit(1)

### Check whether number of threads is appropriate
if args.threads > len(params):
    frameinfo = getframeinfo(currentframe())
    logger.info('Warning in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' +
                'Number of threads is with ' + str(args.threads) +
                ' greater than number of species parameter sets (here: ' +
                str(len(params)) + '). ' +
                'Decreasing number of threads to ' + str(len(params)) +
                '.')
    args.threads = len(params)


### Open genome file, loop over single sequences
seq_to_spec = {}
try:
    with open(args.genome, "r") as genome_handle:
        id = str(uuid.uuid4().hex)[:10]

        curr_seq = ""
        curr_header = ""
        hindex = 0;
        for line in genome_handle:
            if re.match(r'^>', line):
                hindex = hindex + 1
                if len(curr_seq) > 0:
                    fitting_species = augustify_seq(hindex, curr_header, curr_seq, tmp, params, id)
                    seq_to_spec[curr_header] = fitting_species
                curr_header = line;
                curr_seq = "";
            else:
                curr_seq += line
        # process the last sequence
        if hindex > 0 and len(curr_seq) > 0:
            fitting_species = augustify_seq(hindex, curr_header, curr_seq, tmp, params, id)
            seq_to_spec[curr_header] = fitting_species
        else:
            logger.info('Wrong formatted genome fasta file ' + args.genome + ' !')

        filename = f'selected_coordinates_{id}_*.txt'
        list_of_files = glob.glob(filename)

        with open("selected_coordinates.txt", "w") as fout:
            results_dict = {}
            for fileName in list_of_files:
                with open(fileName, 'r') as finp:
                    data = finp.read()
                    dict_i = ast.literal_eval(data)
                    results_dict = results_dict | dict_i
            fout.write(json.dumps(results_dict))

        with multiprocessing.Pool(processes=args.threads) as pool:
            list(pool.map(os.remove, glob.glob(filename)))

except IOError:
    frameinfo = getframeinfo(currentframe())
    logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                args.genome + ' for reading!')

### Compute majority vote: most often selected parameter set wins
if(args.species):
    count_specs = {}
    for seq in seq_to_spec:
        if seq_to_spec[seq] not in count_specs:
            count_specs[seq_to_spec[seq]] = 1
        else:
            count_specs[seq_to_spec[seq]] = count_specs[seq_to_spec[seq]] + 1

    max_setcount = 0
    max_name = ""
    for spec in count_specs:
        if count_specs[spec] > max_setcount:
            max_setcount = count_specs[spec]
            max_name = spec

    print("Best species: " + max_name)

### Merge augustus predictions
if args.prediction_file:
    try:
        with open(tmp + "all_preds.gff", "r") as aug_single_handle:
            subprcs_args = [perl, join_aug_pred]
            result = run_process_stdinput(subprcs_args, aug_single_handle)
            try:
                with open(args.prediction_file, "w") as aug_merged_handle:
                    aug_merged_handle.write(result.stdout.decode('utf-8'))
            except IOError:
                frameinfo = getframeinfo(currentframe())
                logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                            str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                            args.prediction_file + ' for writing!')
    except IOError:
        frameinfo = getframeinfo(currentframe())
        logger.info('Error in file ' + frameinfo.filename + ' at line ' +
                    str(frameinfo.lineno) + ': ' + 'Could not open file ' +
                    tmp + 'all_preds.gff' + ' for reading!')

### Cleanup
shutil.rmtree(tmp)
