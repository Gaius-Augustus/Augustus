#!/usr/bin/perl

# Jonas Behr's tool Transcript Skimmer will report two ORFs in assembled transcripts if strand unspecific RNA-Seq data was used for assembly. This script selects and reports the longer of those two transcripts.

# Katharina Hoff, 11.10.2012

my $usage = "selectLongestORFs.pl input.gff > output.gff";

if (@ARGV != 1) {
    print $usage;
    exit;
}

my 

open(INPUT, "<", $ARGV[0]) or die("Could not open input file $ARGV[0]!\n");
while(<INPUT>){

}
close(INPUT) or die("Could not close input file $ARGV[0]!\n");
