#!/usr/bin/env perl

# Katharina J. Hoff
# November 20th 2018

use strict;
use warnings;
use Getopt::Long;

my $usage = << 'ENDUSAGE';
merge_masking.pl	take two softmasked (same) genome assemblies and mask all bases that are
                    masked in one of the two files (softmasking).

SYNOPSIS

compare_masking file1.fa file2.fa

	file1.fa  softmasked fasta file
	file2.fa  softmasked second fasta file
	file3.fa  output softmasked third fasta file

OPTIONS

    --help    output this help message

WARNING: this script keeps two assemblies in memory, i.e. it is not suitable for large genomes!

ENDUSAGE

my ($help);

GetOptions('help' => \$help);

if($help){
	print $usage;
	exit(1);
}


my %masking1;
my $key;
open(FILE1, "<", $ARGV[0]) or die ("Could not open file $ARGV[0]!\n");
	while(<FILE1>){
		chomp;
		if(m/^>/){
			$masking1{$_} = "";
			$key = $_;
		}else{
			$masking1{$key} .= $_;
		}
	}
close(FILE1) or die ("Could not close file $ARGV[0]!\n");

my %masking2;
open(FILE2, "<", $ARGV[1]) or die ("Could not open file $ARGV[1]!\n");
	while(<FILE2>){
		chomp;
		if(m/^>/){
			$masking2{$_} = "";
			$key = $_;
		}else{
			$masking2{$key} .= $_;
		}
	}
close(FILE2) or die ("Could not close file $ARGV[1]!\n");


open(FILE3, ">", $ARGV[2]) or die ("Could not open file $ARGV[2]!\n");
while( my($key, $value) = each(%masking1)){
	my @arr1 = split(//, $value);
	my @arr2 = split(//, $masking2{$key});
	my $counter = 0;
	print(FILE3 $key."\n");
	foreach(@arr1){
		if( ($_ =~ m/\p{Lowercase}/) or ($arr2[$counter] =~ m/\p{Lowercase}/)){
			print(FILE3 lc($_));
		}else{
			print(FILE3 $_);
		}
		$counter++;
	}
	print(FILE3 "\n");
}
close(FILE3) or die ("Could not close file $ARGV[2]!\n");
