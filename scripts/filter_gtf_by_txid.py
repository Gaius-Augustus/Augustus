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

parser = argparse.ArgumentParser(description='Remove transcripts form a gtf file')
parser.add_argument('-g', '--gtf', required=True, type=str, help="Name of (sorted) gtf file")
parser.add_argument('-l', '--exclude_lst', required=True, type=str, help="Text file with transcript IDs")
parser.add_argument('-o', '--output_file', required=True, type=str, help="Output file in gtf format")
args = parser.parse_args()

# Initialize a set to store the IDs of the transcripts to remove
transcripts_to_remove = set()

# Read the transcript IDs from the file and add them to the set
with open(args.exclude_lst, "r") as transcript_ids_file:
  for line in transcript_ids_file:
    transcripts_to_remove.add(line.strip())

# Initialize the dictionary to store the gene and transcript data
gene_data = {}

# Open the GTF file and use the csv reader to read through the rows
with open(args.gtf, "r") as gtf_file:
  reader = csv.reader(gtf_file, delimiter="\t")
  for row in reader:
    if row[2] == "CDS" or row[2]=="exon" or row[2]=="intron" or row[2]=="start_codon" or row[2]=="stop_codon":
      # Extract the transcript ID from the row  
      transcript_id = row[8].split("transcript_id ")[1].split(";")[0].strip('"')
      gene_id = row[8].split("gene_id ")[1].split(";")[0].strip('"')
    elif row[2]=="transcript":
      transcript_id = row[8].strip()
    elif row[2]=="gene":
      gene_id = row[8].strip()
    if row[2]=="gene":
      gene_data[gene_id] = {"gene_line": row}
    elif row[2]=="transcript":
      gene_data[gene_id][transcript_id]=[row]
    else:
      gene_data[gene_id][transcript_id].append(row)

for tx in transcripts_to_remove:
  match = re.match(r'([^g]*g\d+)\.t\d+', tx)
  gid = match.group(1)
  del gene_data[gid][tx]
  if len(gene_data[gid]) == 1:
    del gene_data[gid]
  
# Open the GTF file and use the csv reader to read through the rows
with open(args.output_file, "w") as safe:
  for key, value in gene_data.items():
    safe.write("\t".join(value['gene_line']) + "\n")
    del value['gene_line']
    for innerkey in value.keys():
      for ele in value[innerkey]:
        safe.write("\t".join(ele) + "\n")


      
