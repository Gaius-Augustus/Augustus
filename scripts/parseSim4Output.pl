#!/usr/bin/perl
#
# parse the output of sim4, fiter the alignments and
# output them to stdout
#
# Mario Stanke, August 17th, 2005

use SplicedAlignment;
use strict;

my $usage = "$0 -- filter the output of sim4\n";
$usage .= "Usage: $0 sim4.out seqfile.fa\n";
$usage .= "output to stdout. first messages. then complete genes. then genes complete only at the 5' end\n";
  
if ($#ARGV != 1) {
    die "Unknown option\n\n$usage";
}
my $simfilename = $ARGV[0];
my $seqfilename = $ARGV[1];
my @alignmentList=();

###########################################################################################
#
# read in the sim4 alignments
#
###########################################################################################


open(SIM, "<$simfilename") or die ("Could not open $simfilename");

my ($sa, $estName, $estLength, $genomicName, $genomicLength, $estExonBegin, $estExonEnd, $exonBegin, $exonEnd, $intronFlag, $quality, $complement);
my @genomicMatches;
$/="--------\n";
while (<SIM>) {
    next unless /seq1 =/;
    />(.*)\n/;
    $estName = $1;
    /seq1 = .*, (\d+) bp/;
    $estLength = $1;
    @genomicMatches = split(/seq2 = /, $_);
    shift @genomicMatches;
    foreach my $match (@genomicMatches) {
	$match =~ /\((.*)\), (\d+) bp/;
	$genomicName = $1;
	$genomicLength = $2;
	$complement = ($match =~ /complement/);
	$sa = new SplicedAlignment($estName, $estLength, $genomicName, $genomicLength, $complement);
	foreach (split (/\n/, $match)){
	    if (/(\d+)-(\d+) +\((\d+)-(\d+)\) +(\d+)%(.*)/) {
		$estExonBegin = $1;
		$estExonEnd = $2;
		$exonBegin = $3;
		$exonEnd = $4;
		$quality = $5 / 100;
		$intronFlag = $6;
		$intronFlag =~ s/ //;
		$sa->addExon($estExonBegin-1, $estExonEnd-1, $exonBegin-1, $exonEnd-1 , $quality, $intronFlag);
	    }
	}
	push @alignmentList, $sa;
    }
}

my @geneList=();

###########################################################################################
#
# filter and enrich the data
#
###########################################################################################

foreach $sa (@alignmentList){
    $sa->checkAlignment();
    print STDERR $sa->get_status(), "\n";
    if ($sa->get_status() eq "OK"){
	my $seq = getSequence($sa->get_contigname);
	$sa->makeGene($seq);
	$sa->findSplicedUTR();
	print STDERR "\t\t", $sa->get_status(), "\n"; 
	if ($sa->get_status() eq "OK"){
	    push @geneList, $sa;
	}
    }
}

###########################################################################################
#
# output the data
#
###########################################################################################

print "### complete genes\n";
foreach $sa (@geneList){
    if ($sa->get_complete5prime() && $sa->get_complete3prime()){
	print $sa->output;
    }
}
print "### genes incomplete at the 3' end and complete at the 5' end\n";
foreach $sa (@geneList){
    if ($sa->get_complete5prime() && !$sa->get_complete3prime()){
	print $sa->output;
    }
}



###########################################################################################
#
# subroutines
#
###########################################################################################

sub getSequence {
    my $seqname = shift;
    my $seq;
    # let cdbyank use the index to search the sequence
    system("echo '$seqname' | cdbyank ${seqfilename}.cidx -d $seqfilename > genomic.fa");
    system ("perl -i.orig -p -e 's/^Incorrectly.*fixed.\n//' genomic.fa");
    open (SEQ, "<genomic.fa") or die ("Could not open genomic.fa");
    $/="\n";
    my @lines=<SEQ>;
    shift @lines;
    $seq = join ("", @lines);
    $seq =~ s/\n//g;
    return $seq;
}
