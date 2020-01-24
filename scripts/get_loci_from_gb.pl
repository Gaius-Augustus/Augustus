#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# get_loci_from_gb.pl                                                                              #
#                                                                                                  #
# Author: Katharina Hoff                                                                           # 
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
####################################################################################################

use strict;
use warnings;
use Getopt::Long;

my $usage = <<'ENDUSAGE';

get_loci_from_gb.pl    script for extracting a loci/transcript list from genbank format file

SYNOPSIS

get_loci_from_gb.pl --gb=train.gb --out=loci.lst

INPUT FILE OPTIONS

--gb=train.gb                       file with gene structures in genbank format
--out=loci.lst                      tabulator separated output file
--help                              print this help message

ENDUSAGE

my $gb;
my $outfile;
my $help;

GetOptions(
    'gb=s' => \$gb,
    'out=s'=> \$outfile,
    'help!'=> \$help);

if($help){
	print $usage;
	exit(1);
}

if(not(defined($gb)) or not (defined($outfile))){
	print("Error: input argument missing!");
	print $usage;
	exit(1);
}

my $txLocus;
my %txInGb3;

open(GB, "<", $gb) or die("Failed to open file $gb for reading!\n");
while(<GB>){
	if ($_ =~ m/LOCUS\s+(\S+)\s/) {
		$txLocus = $1;
	} elsif ($_ =~ m/\/gene="(\S+)"/) {
		$txInGb3{$1} = $txLocus
	}
}
close(GB) or die ("Failed to close file $gb !\n");

open(OUT, ">", $outfile) or die ("Failed to open file $outfile for writing!\n");
foreach (keys %txInGb3) {
	print "$_\t$txInGb3{$_}\n";
}
close(OUT) or die ("Failed to close file $outfile !\n");