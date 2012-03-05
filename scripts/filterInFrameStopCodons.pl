#!/usr/bin/perl

# Katharina J. Hoff
# March 5th 2012
#
# Find all predicted genes that do not contain in-frame stop codons in an amino acid file and return the gene identifier (fasta header).
#
# Input format: multiple fasta with fasta headers that contain only the gene identifier.
# Output format: List of gene identifiers
#
# Note: consider running AUGUSTUS with the option --noInFrameStop=true if you want to avoid stop codons in the first place.

my $usage = "filterInFrameStopCodons.pl protein.fa > no-stop.lst\nConsider running AUGUSTUS with the option --noInFrameStop=true if you want to avoid stop codons in the first place.\n";

if (@ARGV != 1) {
    print $usage;
    exit;
}

my $protein = $ARGV[0];
my %hasStop;
my $currentID;

open(PROT, "<", $protein) or die "Could not open protein file $protein!\n";
while(<PROT>){
    if($_=~m/^>/){
	$_=~s/>//;
	chomp;
	$currentID = $_;
	$hasStop{$currentID} = 0;
    }elsif($_=~m/X/){
	$hasStop{$currentID} = 1;
    }
}
close(PROT) or die "Could not close protein file $protein!\n";

while ( ($id,$stop) = each %hasStop ) {
    if($stop == 0){
	print "$id\n";
    }
}
