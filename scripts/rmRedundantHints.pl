#!/usr/bin/perl

# Katharina J. Hoff, Jan 30th 2012
#
# Scenario:
# E.g. when peptide hints are generated from alternative transcript AUGUSTUS predictions and a six-frame translation, a certain amount of redundant hints is expected. This script discards redundant hint copies (multiplicity is already accounted for). The last column is ignored when checking redundancy.
# Input file format: AUGUSTUS hints format
# Output file format: AUGUSTUS hints format

my $usage = "rmRedundantHints.pl list1.hints list2.hints > nonredundant.hints\n";

if(@ARGV!=2){print STDERR $usage; exit -1;}

my $set1 = $ARGV[0];
my $set2 = $ARGV[1];
my %hints = ();
my @t1;

open(SET1, "<", $set1) or die "Could not open hints file 1 $set1!\n";
while(<SET1>){
	@t1 = split(/\t/);
	if(not(exists($hints{"$t1[0]\t$t1[1]\t$t1[2]\t$t1[3]\t$t1[4]\t$t1[5]\t$t1[6]\t$t1[7]"}))){
		$hints{"$t1[0]\t$t1[1]\t$t1[2]\t$t1[3]\t$t1[4]\t$t1[5]\t$t1[6]\t$t1[7]"} = "$t1[8]";
	}
}
close(SET1) or die "Could not close hints file 1 $set1!\n";

open(SET2, "<", $set2) or die "Could not open hints file 2 $set2!\n";
while(<SET2>){
	@t1 = split(/\t/);
	if(not(exists($hints{"$t1[0]\t$t1[1]\t$t1[2]\t$t1[3]\t$t1[4]\t$t1[5]\t$t1[6]\t$t1[7]"}))){
		$hints{"$t1[0]\t$t1[1]\t$t1[2]\t$t1[3]\t$t1[4]\t$t1[5]\t$t1[6]\t$t1[7]"} = "$t1[8]";
	}
}
close(SET2) or die "Could not close hints file 2 $set2!\n";

while(($k, $v) = each %hints){
	print "$k\t$v";
}
