#!/usr/bin/env perl

# Author: Katharina J. Hoff
# Date: Oct. 27 2017

# This script is part of WebAUGUSTUS
# It check whether scaffold names in a given gff file are present in a related genome fasta file as headers

use strict;
use Getopt::Long;

my $usage = "$0 -- check whether scaffold names in a given gff file are present in a related genome fasta file as headers\n\n";
$usage .= "Usage: $0 --gff=gffFile --genome=genomeFile --out=outFile\n";
$usage .= "An integer 1 will be casted to STDOUT if the gff file contains scaffold names that are NOT present in the genome file, otherwise 0 will be casted.\nAnother integer (0/1) will be cast to outFile concerning the number of columns; if 9, then 0; if different from 9, then 1.";

my $gffFile;
my $genomeFile;
my $outFile;

GetOptions('gff=s'=>\$gffFile,
           'genome=s'=>\$genomeFile,
           'out=s'=>\$outFile);

my %headers;
open(GENOME, "<", $genomeFile) or die("Could not open file $genomeFile!\n");
while(<GENOME>){
    if(m/^>/){
	chomp;
	$_ =~ s/>//;
	$headers{$_} = 1;
    }
}
close(GENOME) or die("Could not close file $genomeFile!\n");

my @line;
my $error1 = 0;
my $error2 = 0;
my $nCols;
open(GFF, "<", $gffFile) or die("Could not open file $gffFile!\n");
while(<GFF>){
    @line = split(/\t/);
    if(not(defined($headers{$line[0]}))){
	$error1 = 1;
    }
    $nCols = @line;
    if(not($nCols==9) && not($nCols==0)){
	$error2 = 1;
    }

    
    
}
close(GFF) or die("Could not close file $gffFile!\n");

print STDOUT "$error1";

open(OUT, ">", $outFile) or die ("Could not open file $outFile!\n");
print OUT $error2;
close(OUT) or die ("Could not close file $outFile!\n");
