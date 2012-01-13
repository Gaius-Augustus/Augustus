#!/usr/bin/perl

# Creates a fasta file with unique peptides from a fasta input file
# New fasta headers are random
# Creates two output files: unique.fa header.mapping
# Do not use for long protein sequences, only for short peptides!!!

my $usage = "uniquePeptides.pl in.fa unique.fa header.mapping\n";

if(@ARGV!=3){print STDERR $usage; exit -1;}

my $in = $ARGV[0];
my $uniquefa = $ARGV[1];
my $mapping = $ARGV[2];

my %mapHash = ();
my %seqHash = ();
my $maxLen = 100;

my $header;
my $c = 1;

open(IN, "<", $in) or die("Could not open file $in!\n");

while ( <IN> ) { 
	chomp;
	if($_=~m/^>/){
		$_=~s/>//;
		$header = $_;
	}else{
		if(length($_) > $maxLen){
			print STDERR "The peptide $_ is longer $maxLen. Aborting.\n";
			exit -1;
		}
		if(not(exists($seqHash{$_}))){
			$seqHash{$_} = ">pep$c";
			$mapKey = ">pep$c";
			$c++;
			
		}
		push (@{$mapHash{">pep$c"}}, $header);
	}
}

close IN or die("Could not close file $in!\n");

# print unique fasta file
open(OUT1, ">", $uniquefa) or die("Could not open file $uniquefa!\n");
while ( ($k,$v) = each %seqHash ) {
    print OUT1 "$v\n$k\n";
}
close(OUT1) or die("Could not close file $uniquefa!\n");

# print mapping file
open(OUT2, ">", $mapping) or die("Could not close file $mapping!\n");
for $k ( keys %mapHash ) {
	print OUT2 "$k:\n";
	foreach(@{$mapHash{$k}}){
		print OUT2 "$_\n";
	}
}
close(OUT2) or die("Could not close file $mapping!\n");

