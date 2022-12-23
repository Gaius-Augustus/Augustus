#!/usr/bin/env python3

__author__ = "Katharina J. Hoff with a little help of ChatGPT"
__copyright__ = "Copyright 2022. All rights reserved."
__credits__ = "Katharina Hoff"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

# Import the necessary libraries
import csv
import re
import argparse

parser = argparse.ArgumentParser(description='Identify single exon transcripts in a gtf file')
parser.add_argument('-g', '--gtf', required=True, type=str, help="Name of (sorted) gtf file")
parser.add_argument('-o', '--output_file', required=True, type=str, help="Output plain text file")
args = parser.parse_args()

# First, we will read the GTF file into a list of lines
with open(args.gtf, 'r') as gtf_file:
    gtf_lines = gtf_file.readlines()

# Next, we will create a dictionary to store the transcript IDs of single exon genes
single_exon_transcripts = {}

# Then, we will iterate through each line in the GTF file
for line in gtf_lines:
    # Split the line by tabs to get the different fields
    fields = line.split('\t')
    if fields[2] == "CDS":
      # Get the transcript_id field
      transcript_id = fields[8]
      # Check if the transcript_id is already in the dictionary
      if transcript_id not in single_exon_transcripts:
        # If it is not, add it to the dictionary and set the value to 1
        single_exon_transcripts[transcript_id] = 1
      else:
        # If it is, increment the value by 1
        single_exon_transcripts[transcript_id] += 1

with open(args.output_file, "w") as out:
  # Finally, we will iterate through the dictionary and print the transcript IDs of single exon genes
  for transcript_id, exon_count in single_exon_transcripts.items():
    if exon_count == 1:
        match = re.search(r'([^0-9;"]+\d+\.t\d+)', transcript_id)
        if(match):
            out.write(match.group(1) + "\n")

      
