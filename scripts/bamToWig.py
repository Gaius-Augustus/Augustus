#!/usr/bin/env python3

# Convert bam files to wiggle files for usage with AUGUSTUS (exonpart hints)
# Alternative to Augustus/auxprogs/bam2wig/bam2wig, which often causes
# compilation problems. However, this script also has a number of dependencies:
# UCSC tools twoBitInfo, fatToTwoBit,  Unix sort, samtools.
# Unix sort and samtools are usually installed on AUGUSTUS user PCs,
# UCSC tools are automatically obtained if not present on machine.
# In comparison to bam2wig, this script requires the genome file
# as additional input argument. In return, it ensures that aligments
# do not exceed sequence boundaries (which sometimes happens and is not caught
# by bam2wig).
# This script is a re-assembly of functions from MakeHub
# at https://github.com/Gaius-Augustus/MakeHub

import string
import random
import os
import os.path
import argparse
import re
import urllib.request
import subprocess
import platform
from inspect import currentframe, getframeinfo
import shutil


''' ******************* BEGIN HEAD ******************************************'''
''' Contains arg parser & software requirement checks            '''


__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2019. All rights reserved."
__credits__ = "MakeHub project"
__license__ = "Artistic Licsense"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"


ucsc_tools = {'twoBitInfo': '', 'faToTwoBit': ''}

parser = argparse.ArgumentParser(
    description='Convert bam file to wiggle format for usage with AUGUSTUS ' +
    'as exonpart hints.')
parser.add_argument('-b', '--bamFile', required=True, type=str,
                    help="Input file in Bam format.")
parser.add_argument('-g', '--genomeFile', required=True, type=str,
                    help="Input genome file in FASTA format.")
parser.add_argument('-o', '--outFile', required=True, type=str,
                    help="Output file in wiggle format.")
parser.add_argument('-s', '--SAMTOOLS_PATH', required=False, type=str,
                    help="Path to samtools executable, e.g. \'/usr/bin\'.")
args = parser.parse_args()

''' This script contains bash calls for samtools. It has not been tested
    on OS X or Windows. Therefore check platform prior further actions. '''

plat_sys = platform.system()
if plat_sys != "Linux":
    frameinfo = getframeinfo(currentframe())
    print('Warning in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': '
          + "This script has been developed and tested on Linux, " +
          "only. Your operating system is " + plat_sys +
          ". The script " +
          "might not work properly on your system!")

''' This script depends on samtools pileup. Check whether samtools is
    available. '''

samtools = ""
if args.SAMTOOLS_PATH:
    samtools = args.SAMTOOLS_PATH + "/samtools"
    if not(os.access(samtools, os.X_OK)):
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + samtools + " is not executable!")
        exit(1)
else:
    if shutil.which("samtools") is not None:
        samtools = shutil.which("samtools")
    else:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': '
              + "Unable to locate samtools binary!")
        print("samtools is available as package in many Linux distributions.")
        print("For example, on Ubuntu, try installing with:")
        print("\"sudo apt install samtools\"")
        print("If samtools is unavailable as a package, you can obtain it " +
              "from github at:")
        print("https://github.com/samtools/samtools")
        exit(1)

''' Find bash sort '''

sort_tool = shutil.which('sort')
if sort_tool is None:
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': ' + "Unable to locate bash tool 'sort'")
    print('sort is part of most Linux distributions. On Ubuntu, it is ' +
          'part of the package coreutils. Try re-installing your bash if sort' +
          ' is missing on your system.')
    quit(1)


''' Find or obtain UCSC tools  (here, tool for chromsome sizes file).
    script has never been tested on Darwin, tool download is implemented
    for Darwin though. '''

arch = platform.machine()
for key, val in ucsc_tools.items():
    if shutil.which(key) is not None:
        ucsc_tools[key] = shutil.which(key)
    elif os.path.isfile(os.getcwd() + "/" + key):
        ucsc_tools[key] = os.getcwd() + "/" + key
        if not(os.access(ucsc_tools[key], os.X_OK)):
            os.chmod(ucsc_tools[key], 0o777)
    else:
        if not(arch == 'x86_64'):
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "This script depends on binaries that are available for " +
                  "x86_64 architecture for linux and MacOX. " +
                  "We have determined that your system architecture is " +
                  arch + "." +
                  " Please try downloading " + key +
                  " for your architecture from: " +
                  "http://hgdownload.soe.ucsc.edu/admin/exe")
            exit(1)
        elif not(plat_sys == 'Linux') and not(plat_sys == 'Darwin'):
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "This script depends on binaries that are available for " +
                  "x86_64 architecture for linux and MacOX. " +
                  "We have determined that your system is " +
                  plat_sys + "." +
                  " Please try downloading " + key +
                  " for your operating system from: " +
                  "http://hgdownload.soe.ucsc.edu/admin/exe")
            exit(1)
        else:
            if plat_sys == 'Linux':
                tool_url = "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/" + key
            elif plat_sys == 'Darwin':
                tool_url = "http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/" + key
            print("Was unable to locate " + key +
                  " on your system, will try to download it from " + tool_url + "...")
            with urllib.request.urlopen(tool_url) as response, open(key, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            ucsc_tools[key] = os.getcwd() + "/" + key
            os.chmod(ucsc_tools[key], 0o777)

''' ******************* END HEAD ********************************************'''

''' ******************* BEGIN FUNCTIONS *************************************'''

''' Function that runs a subprocess with arguments '''


def run_simple_process(args_lst):
    try:
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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


''' Function that runs a subprocess with input from STDIN '''


def run_process_stdinput(args_lst, byte_obj):
    try:
        result = subprocess.run(args_lst, input=byte_obj,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "run_process_stdinput: return code of subprocess was "
                  + str(result.returncode))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Failed executing: ", " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


''' Function that writes subprocess byte object to flat file '''


def write_byteobj(byte_obj, outfile):
    try:
        with open(outfile, 'w') as byteobj_handle:
            byteobj_handle.write(byte_obj.decode('utf-8'))
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that converts bam file to wig file with RSeQC '''


def bamToWig(bam_file, wig_file, size_file, mpileup_file):
    subprcs_args = [samtools, "mpileup", "-o", mpileup_file, bam_file]
    run_simple_process(subprcs_args)
    # observed that in rare cases a wig file might contain coverage for one
    # more base than present in the sequence; seems to be an alignment
    # error, not a bam2wig error, because the same problem arises if I
    # convert bam to wig in python in a different way
    # therefore check wiggle file for sanity and modify if required
    # (i.e. cleave coverage at sequence end)
    chrom_sizes = {}
    try:
        with open(size_file, "r") as size_handle:
            for line in size_handle:
                split_list = re.split(r'\t', line.rstrip('\n'))
                chrom_sizes[split_list[0]] = int(split_list[1])
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              size_file + " for reading!")
    try:
        with open(mpileup_file, "r") as pileup_handle:
            try:
                with open(wig_file, "w") as wig_handle:
                    wig_handle.write(
                        "track name=" + bam_file + " type=wiggle_0\n")
                    lastSeq = ""
                    lastStart = 0
                    for line in pileup_handle:
                        seq, start, t1, depth, t2, t3 = line.split()
                        if (seq != lastSeq) and (start != lastStart):
                            wig_handle.write(
                                "variableStep chrom=" + seq + "\n")
                        if(int(start) <= chrom_sizes[seq]):
                            wig_handle.write(start + " " + depth + "\n")
                        lastSeq = seq
                        lastStart = start
            except IOError:
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      "Could not open file " + wig_file + " for writing!")
    except IOError:
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Could not open file " + pilup_file + " for reading!")


''' Function that generates a random string of fixed length, default length is
    10.'''


def randomString(stringLength=10):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


''' ******************* END FUNCTIONS ***************************************'''


''' Generate random file prefix for temporary files '''

prefix = randomString()
twoBitFile = "twoBitFile_" + prefix
ChromSizes_file = "chromSizes_" + prefix
pileup_file = "pileup_" + prefix

''' Generate chrosome sizes file '''

subprcs_args = [ucsc_tools['faToTwoBit'], args.genomeFile, twoBitFile]
run_simple_process(subprcs_args)
subprcs_args = [ucsc_tools['twoBitInfo'], twoBitFile, 'stdout']
result = run_simple_process(subprcs_args)
subprcs_args = [sort_tool, '-k2rn']
result = run_process_stdinput(subprcs_args, result.stdout)
write_byteobj(result.stdout, ChromSizes_file)


''' Generate Wiggle file '''

bamToWig(args.bamFile, args.outFile, ChromSizes_file, pileup_file)

''' Delete temporary files '''

os.remove(twoBitFile)
os.remove(ChromSizes_file)
os.remove(pileup_file)
