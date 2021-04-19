#!/usr/bin/env perl
#
# parse psl format, fiter the alignments and
# output them to stdout
#
# Katharina J. Hoff, April 14th 2021, katharina.hoff@uni-greifswald.de

use File::Basename;
use lib dirname (__FILE__);
use SplicedAlignment;
use strict;
use Getopt::Long;
use File::Which qw(which);

my $usage = "$0 -- parse psl format to orf\n";
$usage .= "Usage: $0 in.psl seqfile.fa\n";
$usage .= "output to stdout. first messages. then complete genes. then genes complete only at the 5' end\n";
  
if ($#ARGV != 1) {
    die "Unknown option\n\n$usage";
}
my $pslFile = $ARGV[0];
my $seqFile = $ARGV[1];
my @alignmentList=();

###########################################################################################
#
# read genome into memory
#
###########################################################################################

my %genome;
my $header;
open(GENOME, "<$seqFile") or die("Could not open $seqFile!\n");
while(<GENOME>){
    chomp;
    if(m/^>(\S+)/){
        $header = $1;
    }else{
        $genome{$header} .= $_;
    }
}
close(GENOME) or die("Could not close $seqFile!\n");

###########################################################################################
#
# read in the psl alignments
#
###########################################################################################

open(PSL, "<$pslFile") or die ("Could not open $pslFile\n");

my ($estName, $estLength, $genomicName, $genomicLength, $estExonBegin, $estExonEnd, $exonBegin, $exonEnd, $intronFlag, $quality, $complement);
my @genomicMatches;

while (<PSL>) {
    my @pslLine = split(/\t/, $_);
    $estName = $pslLine[9];
    $estLength = $pslLine[12];
    $genomicName = $pslLine[13];
    $genomicLength = $pslLine[14]; # it is unclear whether the genome sequence size is meant or the range of the alignment on the genome sequence 
    $complement = ($pslLine[8] =~ m/-/);
    my $sa = new SplicedAlignment($estName, $estLength, $genomicName, $genomicLength, $complement);
    my $nBlocks = $pslLine[17];
    my @blockSizes = split(/,/, $pslLine[18]);
    my @qStarts = split(/,/, $pslLine[19]); # zero based, we want one based in the end, I think
    my @tStarts = split(/,/, $pslLine[20]); # zero based
    for(my $i = 0; $i < $nBlocks; $i++){
        $estExonBegin = $qStarts[$i] + 1;
        $estExonEnd = $estExonBegin + $blockSizes[$i] - 1;
        $exonBegin = $tStarts[$i] + 1;
        $exonEnd = $exonBegin + $blockSizes[$i] - 1;;
        $quality = 1; # don't know how to compute quality for single exons from psl
        if($i < ($nBlocks-1)){
            $intronFlag = 1;
        }else{
            $intronFlag = 0;
        }
        $sa->addExon($estExonBegin-1, $estExonEnd-1, $exonBegin-1, $exonEnd-1 , $quality, $intronFlag);
    }
	push @alignmentList, $sa;
}

close(PSL) or die ("Could not close $pslFile!\n");

my @geneList=();

###########################################################################################
#
# filter and enrich the data
#
###########################################################################################

foreach my $sa (@alignmentList){
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
foreach my $sa (@geneList){
    if ($sa->get_complete5prime() && $sa->get_complete3prime()){
	   print "This is a complete gene: ".$sa->output;
    }
}
print "### genes incomplete at the 3' end and complete at the 5' end\n";
foreach my $sa (@geneList){
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
    return $genome{$seqname};
}
