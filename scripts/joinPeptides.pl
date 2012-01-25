#!/usr/bin/perl

# Creates a fasta file with unique peptides from two peptide fasta input files
# Deletes redundant entries. Multiplicity information is taken from the first
# of two fasta input files.

my $usage = "joinPeptides.pl in1.fa in2.fa > out.fa\n\nRead script head comments for futher documentation!\n\n";

if(@ARGV!=2){print STDERR $usage; exit -1;}

my $in1 = $ARGV[0];
my $in2 = $ARGV[1];

my %peptides = ();

my $header;
my $maxLen = 100;

open(IN, "<", $in1) or die("Could not open file $in1!\n");

while ( <IN> ) { 
	chomp;
	if($_=~m/^>/){
		$header = $_;
	}else{
		if(length($_) > $maxLen){
			print STDERR "The peptide $_ is longer $maxLen. Aborting.\n";
			exit -1;
		}
		while($_=~m/-/){ # delete weird dashes in sequences
			$_ = s/-//;
		}
		if(not(exists($locHash{$_}))){
			$peptides{$_} = "$header";			
		}
		
		
	}
}
close IN or die("Could not close file $in1!\n");

open(IN2, "<", $in2) or die("Could not open file $in2!\n");
while( <IN2>){
		chomp;
		if($_=~m/^>/){
		$header = $_;
	}else{
		if(length($_) > $maxLen){
			print STDERR "The peptide $_ is longer $maxLen. Aborting.\n";
			exit -1;
		}
		while($_=~m/-/){ # delete weird dashes in sequences
			$_ = s/-//;
		}
		if(not(exists($locHash{$_}))){
			$peptides{$_} = "$header";			
		}
		
		
	}
}
close IN2 or die("Could not close file $in2!\n");
# print unique fasta file
while ( ($k,$v) = each %peptides ) {
	print "$v\n";
	print "$k\n";
}

