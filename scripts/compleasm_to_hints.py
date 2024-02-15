#!/usr/bin/env python3

"""
Description: This script runs compleasm on a genome with a given BUSCO partition
             and parses the output table for complete BUSCOs to generate a hints
             file for AUGUSTUS/BRAKER/GALBA.
Author: Katharina J. Hoff
Email: katharina.hoff@uni-greifswald.de
Date: November 27th, 2023

Copyright (C) 2023, Katharina J. Hoff, University of Greifswald

This program is free software; you can redistribute it and/or modify
it under the terms of the Artistic License.
"""

# install instructions for sepp that will be needed for containers:
# git clone https://github.com/smirarab/sepp.git
# sudo apt-get install default-jre
# python3 setup.py config -c
# sudo python3 setup.py install
# other than that, the script depends on compleasm.py and its dependencies


import argparse
import re
import shutil
import os
import csv
import subprocess
from inspect import currentframe, getframeinfo

''' Function that runs a subprocess with arguments '''


def run_simple_process(args_lst):
    try:
        # bed files need sorting with LC_COLLATE=C
        myenv = os.environ.copy()
        myenv['LC_COLLATE'] = 'C'
        print("Trying to execute the following command:")
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=myenv)
        print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
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


def extract_tx_ids_from_tsv(tsv_file):
    busco_ids = {}
    with open(tsv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if (row['Status'] == 'Single' or row['Status'] == 'Duplicated') and (row['Frameshift events'] == '0'):
                busco_ids[row['Best gene']] = True
    return busco_ids

def read_and_filter_gff(gff_file, id_dict):
    gff_list = []
    try:
        with open(gff_file, 'r') as gff_handle:
            for line in gff_handle:
                if not line.startswith("#"):
                    if re.search(r'Target=', line):
                        bid = re.search(r'Target=(\S+)', line).group(1)
                        if bid in id_dict.keys():
                            gff_list.append(line)
    except IOError:
        print("Could not open file ", gff_file)
    return gff_list

def miniprot_to_hints(gff_lines):
    hints_lines = []
    for line in gff_lines:
        line = line.rstrip()
        line_fields = line.split("\t")
        if re.search(r'\tCDS\t.*Parent=([^;]+)', line):
            if (int(line_fields[4])-3) - (int(line_fields[3])+3) + 1 > 0:
                regex_result = re.search(r'\tCDS\t.*Parent=([^;]+)', line)
                hint = line_fields[0] + "\t" + "c2h" + "\tCDSpart\t" + str(int(line_fields[3])+3) + "\t" + str(int(line_fields[4])-3) + "\t" + "1" + "\t" + line_fields[6] + "\t" + line_fields[7] + "\t" + "src=M;grp=" + regex_result.group(1) + ";pri=4\n"
                hints_lines.append(hint)
    return hints_lines



def main():
    parser = argparse.ArgumentParser(description="Run compleasm and generate a gtf file with complete BUSCO genes.")
    
    # Mandatory input arguments
    parser.add_argument("-g", "--genome", required=True, help="Genome file in fasta format")
    parser.add_argument("-d", "--database", required=True, help="BUSCO database to use")
    parser.add_argument("-p", "--compleasm", required=False, help="Location of compleasm binary incl. binary name")
    parser.add_argument("-t", "--threads", required=False, help="Number of threads to use, default is 1")
    parser.add_argument("-o", "--output", required=True, help="Output file name, file is in GTF format.")
    parser.add_argument("-s", "--scratch_dir", required=False, help="Temporary directory for compleasm output, default is current directory, must be writable.")
    args = parser.parse_args()

    if args.compleasm is None and shutil.which('compleasm.py'):
        # Try to find 'compleasm' in the system PATH
        args.compleasm = shutil.which('compleasm')
    elif args.compleasm is None and os.environ.get('COMPLEASM_PATH'):
        # If 'compleasm' is not in PATH, check for an environment variable
        args.compleasm = str(os.environ.get('COMPLEASM_PATH')) + '/compleasm.py'
    elif args.compleasm is None:
        # If still not found, raise an error
        raise FileNotFoundError("compleasm is not in PATH and COMPLEASM_PATH is not set")
    if args.compleasm is not None:
        # check whether provided compleasm is executable
        if not os.access(args.compleasm, os.X_OK):
            raise FileNotFoundError("compleasm is not executable")
        
    # check whether the database has the ending odb_10, if not, add it
    if not args.database.endswith("_odb10"):
        args.database = args.database + "_odb10"

    # apply compleasm to genome file with run_subprocess and the database
    if args.scratch_dir is None:
        args.scratch_dir = "compleasm_genome_out"
    run_simple_process([args.compleasm, 'run', '-l', args.database, '-a', args.genome, '-t', args.threads, '-o', args.scratch_dir])

    # parse compleasm output table for complete BUSCOs without frameshifts
    busco_ids = extract_tx_ids_from_tsv(args.scratch_dir + '/' + args.database + '/full_table.tsv')
    
    # filter the miniprot alignments for those that have no frame shifts
    gff_lines = read_and_filter_gff(args.scratch_dir + '/' + args.database + '/miniprot_output.gff', busco_ids)

    # convert the miniprot lines to CDSpart hints
    hints_lines = miniprot_to_hints(gff_lines)

    # gff_lines is not compatible with getAnnoFastaFromJoingenes, we need to fix this!
    try:
        with open(args.output, "w") as out_handle:
            for line in hints_lines:
                out_handle.write(line)
    except IOError:
        print("Failed to open file", args.out)

    # print the BUSCO scores for genome level statistics to STDOUT
    print("The following BUSCOs were found in the genome:")
    try:
        with open(args.scratch_dir + '/' + 'summary.txt', 'r') as busco_summary:
            for line in busco_summary:
                line = line.strip()
                print(line)
    except IOError:
        print("Failed to open file", args.scratch_dir + '/' + 'summary.txt')

    # delete the temporary directory
    shutil.rmtree(args.scratch_dir)

if __name__ == "__main__":
    main()