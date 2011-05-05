#!/usr/bin/perl

# Katharina J. Hoff, May 5th 2011
#
# PASA sometimes has problems with weird fasta headers.
# This script converts the headers of a fasta file to "simple headers" without any special characters or spaces

my $USAGE="simpleFastaHeaders.pl in.fa prefix out.fa mapping.txt\n\n in.fa - the file to be reformatted\n prefix - the prefix of every new header, e.g. contig\n out.fa - the reformatted fasta file\n mapping.txt - a tab-separated mapping table (newName \t oldName)\n\n";

if(@ARGV!=4){print STDERR $USAGE; exit -1}

my $inFile = $ARGV[0];
my $prefix = $ARGV[1];
my $outFa = $ARGV[2];
my $mapping = $ARGV[3];
my $preC = 1; # counter for prefix

open(INFILE, "<", $inFile) or die("Could not open file $inFile!\n");
open(OUTFA, ">", $outFa) or die ("Could not open file $outFa!\n");
open(MAP, ">", $mapping) or die ("Could not open file $mapping!\n");

while(<INFILE>){
    if($_=~m/^>/){
	chomp;
	print MAP ">".$prefix.$preC."\t$_\n";
	print OUTFA ">".$prefix.$preC."\n";
	$preC = $preC+1;
    }else{
	print OUTFA $_;
    }
}



close(INFILE) or die("Could not close file $inFile!\n");
close(OUTFA) or die("Could not close file $outFa!\n");
close(MAP) or die("Could not close file $mapping!\n");
