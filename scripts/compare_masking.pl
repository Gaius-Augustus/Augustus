#!/usr/bin/env perl

# Katharina J. Hoff
# November 20th 2018


use strict;
use warnings;
use Getopt::Long;

my $usage = << 'ENDUSAGE';
compare_masking	compare the repeat masking content of differently masked (same) assemblies

SYNOPSIS

compare_masking file1.fa file2.fa

	file1.fa  softmasked fasta file
	file2.fa  softmasked second fasta file

OPTIONS

    --help    output this help message

WARNING: This script keeps two assemblies in memory, i.e. it is not suitable for large genomes!

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


my $only1 = 0;
my $only2 = 0;
my $both = 0;
while( my($key, $value) = each(%masking1)){
	my @arr1 = split(//, $value);
	my @arr2 = split(//, $masking2{$key});
	my $counter = 0;
	foreach(@arr1){
		if( ($_ =~ m/\p{Lowercase}/) and ($arr2[$counter] =~ m/\p{Lowercase}/)){
			$both++;
		}elsif($_ =~ m/\p{Lowercase}/){
			$only1++;
		}elsif($arr2[$counter] =~ m/\p{Lowercase}/){
			$only2++;
		}
		$counter++;
	}
}

print("Masked in both files: $both\n");
print("Masked in File1, only: $only1\n");
print("Masked in File2, only: $only2\n");
