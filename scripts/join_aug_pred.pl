#!/usr/bin/perl
#
# join_aug_pred.pl
# Join AUGUSTUS predictions from separate runs of possibly overlapping sequence segments.
# Mario Stanke, 25.10.2006
#
use strict;
use Getopt::Long;

my $usage = "$0 -- join AUGUSTUS predictions from separate runs of possibly overlapping sequence segments.\n\n";
$usage .= "Usage: $0 < augustus.concat > augustus.joined\n\n";
$usage .= "Reads from standard input a file with AUGUSTUS outputs. This file will typically be the concatenation\n";
$usage .= "of the outputs of several simulaneous runs of AUGUSTUS of different sequence segments on a computer cluster.\n";
$usage .= "The sequence segments can be overlapping, in which case this script intelligently omits one version of\n";
$usage .= "overlapping genes from different runs. The outputs from the runs must have been concatenated IN ORDER.\n";
$usage .= "The lines '# This output was generated with AUGUSTUS' must still be present in the input file so the\n";
$usage .= "boundaries between runs can be determined.\n";
$usage .= "The output file contains a nonredundant set of renumbered genes.\n";
$usage .= "Options:\n";
$usage .= "   droplist=s   file with list of gene identifiers that should not be included in output.\n";

my $geneNumber = 0;
my $droplist;
my %drop = ();

GetOptions('droplist:s'=>\$droplist);

if ($#ARGV != -1) {
    die "Unknown option\n\n$usage";
}

if (defined($droplist)){
    open (DROP, "<$droplist") or die ("Could not open file $droplist.");
    while (<DROP>){
	chomp;
	s/\.t\d+//;
	$drop{$_}=1;
    }
    print STDERR "Excluding " . scalar(keys %drop) . " gene ids from joining.\n";
}


# Run: the results for one segment
# I use a hash with two elements:
# header: string with the head part preceeding the first gene
# genes: reference to an array of Genes
# seqname: sequence name

# Gene:
# Hash with the following elements
# begin: begin position of the transcribed part
# end: end position of the transcribed part
# gff: string with the complete AUGUSTUS output about this gene

my $printAllHeaders = 0; # if 0, only the first header is printed
my $lastRun=0;
my $currentRun=0;
my $line;
my $runSeparator = "\# This output was generated with AUGUSTUS";
my $geneid = 1;
my $headerPrinted = 0;
my $gff3flag = 0;

# start main loop that gets run by run
do {
    if (($currentRun = getNextRun())) {
	if ($lastRun) {
	    cleanRedundant($lastRun, $currentRun);
	    if ($gff3flag && $geneid == 1) {
		print "##gff-version 3\n";
	    }
	    printRun($lastRun);
	}
	$lastRun = $currentRun;
    }
} until ($currentRun == 0);

if ($lastRun) {
    printRun($lastRun);
}

sub getNextRun {
    my %run;
    my $numSeps = 0;
    my $geneNr = 0;
    my $inGene = 0;
    my $endHeader = 0;
    do {
	if (defined $line) {
	    if ($line =~ /\#\#gff-version 3/){
		$gff3flag = 1;
	    }
	    if ($line =~ /^$runSeparator/) {
		$numSeps++;
	    }
	    if ($numSeps == 1) {
		if ($line =~ /^\#\#\# gene g/ || $line =~ /^\# start gene g/ ) { # the first option for back compatibility with a previous AUGUSTUS version
		    if ($line =~ /gene (g\d+)/ && $drop{$1}){
			print STDERR "dropping $1\n";
		    }  else {
			$geneNr++;
			$inGene=1;
		    }
		}
		if ($geneNr == 0) {
		    if ( $line=~ /prediction on sequence number/ || 
			 $line =~ /Looks like .* is in .* format/) { # this excludes the hints and hint messages for the first run from the header
			$endHeader = 1;
		    }
		    if (!$endHeader){
			$run{"header"}.= $line;
		    }
		}
		if ($geneNr > 0 && $inGene) {
		    if ($line =~ /^(\S+)\t.*\tgene\t(\d+)\t(\d+)\t/) {
			my $seqname = $1;
			my $begin = $2;
			my $end = $3;
			if (!defined $run{"seqname"}){
			    $run{"seqname"} = $seqname;
			}
			$run{"genes"}->[$geneNr-1]->{"begin"} = $begin;
			$run{"genes"}->[$geneNr-1]->{"end"} = $end;
		    }
		    $run{"genes"}->[$geneNr-1]->{"gff"} .= $line;
		}
		if ($line =~ /^\#\#\# end gene g/ || $line =~ /^\# end gene g/) { # the first option for back compatibility with a previous AUGUSTUS version
		    $inGene = 0;
		}
	    }
	}
	if ($numSeps<2) {
	    $line = <>;
	}
    } until (!(defined $line) || $numSeps >1);
    if (defined $line || $geneNr>0){
	return \%run;
    } else {
	return 0;
    }
};


sub printRun {
    my $run = shift @_;
    if (!$headerPrinted || $printAllHeaders) {
#	print "================ header ==================\n";
	print $run->{"header"};
	$headerPrinted = 1;
    }
#    print "================ genes ==================\n";
    if (defined $run->{"genes"}) {
	foreach my $gene (@{$run->{"genes"}}) {
#	    print "============= gene ======= begin = " .$gene->{"begin"} . " end = " . $gene->{"end"} . " =======\n";
	    foreach my $gffline (split /\n/, $gene->{"gff"}){
		$gffline =~ s/\bg(\d+)\b/g$geneid/g;
		print "$gffline\n";
	    }
	    $geneid++;
	}
    }
}

#
# Clean redundant genes of two possibly overlapping neighboring runs.
# Run1 starts left of run2.
#
sub cleanRedundant {
    my $run1 = shift @_;
    my $run2 = shift @_;
    if (!(defined $run1->{"seqname"}) || !(defined $run2->{"seqname"}) || $run1->{"seqname"} ne $run2->{"seqname"}){
	return;
    }
    if (!(defined $run1->{"genes"}) || !(defined $run2->{"genes"}) || @{$run1->{"genes"}} == 0 || @{$run2->{"genes"}} == 0) {
	if (!(defined $run1->{"genes"})){
	    return;
	}
	if (!(defined $run2->{"genes"})){
	    return;
	}
	return;
    }
 #   print "Anzahl Gene in Run1: " . (scalar @{$run1->{"genes"}}) . "\n";
 #   print "Anzahl Gene in Run2: " . (scalar @{$run2->{"genes"}}) . "\n";
    my $firstgene2 = $run2->{"genes"}->[0];
    my $firstbegin2 = $firstgene2->{"begin"};
    my $lastend1 = -1;
    my $lastgene1;
    foreach my $gene1 (@{$run1->{"genes"}}) {
	if ($gene1->{"end"} > $lastend1) {
	    $lastend1 = $gene1->{"end"};
	    $lastgene1 = $gene1;
	}
    } 
    my $lastend2 = -1;
    foreach my $gene2 (@{$run2->{"genes"}}) {
	if ($lastend2 == -1 || $gene2->{"end"} > $lastend2) {
	    $lastend2 = $gene2->{"end"};
	}
    }
    if ($lastend2 < $run1->{"genes"}->[0]->{"begin"}){
	print STDERR "Prediction runs are not in the right order in sequence "  . $run1->{"seqname"} . ". Please sort along the chromosome first.\n";
    }

 #   print "rightmost end  of gene in run1 : " . $lastend1 . "\n";
 #   print "leftmost begin of gene in run2 : " . $firstbegin2 . "\n";
    if ($lastend1 < $firstbegin2) {
	return; # no overlap between any genes
    }
    # the two runs really have overlapping genes
    # find a good breakposition and then remove all genes 
    # in run1 ending to the right of the breakposition and all genes
    # in run2 beginning at or to the left of breakposition.
    #                                                     lastgene1
    # run1  ------              -------       ------      ---------
    #          -------                                      ---
    # run2                               -----------      ------------------------        -----------         -----
    #                                    firstgene2                                                        ---------------
    #                                    |-d2-|                   |-     d1     -|
    #                                                    |breakpoint
    my $breakpoint;
    my $d1=-1; 
    my $d2=-1;
    
    foreach my $gene2 (@{$run2->{"genes"}}) {
	if ($gene2->{"begin"} <= $lastend1 && $gene2->{"end"} - $lastend1 > $d1) {
	    $d1 = $gene2->{"end"} - $lastend1;
	}
    }
#    print "d1=$d1\n";
    foreach my $gene1 (@{$run1->{"genes"}}) {
	if ($gene1->{"end"} >= $firstbegin2 && $firstbegin2 - $gene1->{"begin"} > $d2) {
	    $d2 = $firstbegin2 - $gene1->{"begin"};
	}
    }
#    print "d2=$d2\n";
    if ($d1 >= $d2) {
	# set the breakpoint directly left to the leftmost begin of any gene in run2 that
	# ends at or after lastend1.
	$breakpoint = $lastend1;
	foreach my $gene2 (@{$run2->{"genes"}}) {
	    if ($gene2->{"end"} >= $lastend1 && $gene2->{"begin"}-1 < $breakpoint) {
		$breakpoint = $gene2->{"begin"}-1;
	    }
	}
    } else {
	# set the breakpoint to the rightmost end of any gene in run1 that
	# begins at or before firstbegin2.
	$breakpoint = $firstbegin2;
	foreach my $gene1 (@{$run1->{"genes"}}) {
	    if ($gene1->{"begin"} <= $firstbegin2 && $gene1->{"end"} > $breakpoint) {
		$breakpoint = $gene1->{"end"};
	    }
	}
    }
#    print "breakpoint = $breakpoint\n";
    #now delete the redundant genes from both lists
    my @newgenes1 = ();
    my @newgenes2 = ();
    foreach my $gene1 (@{$run1->{"genes"}}) {
	if ($gene1->{"end"} <= $breakpoint) {
	    push @newgenes1, $gene1;
	}
    }
    $run1->{"genes"} = \@newgenes1;
    foreach my $gene2 (@{$run2->{"genes"}}) {
	if ($gene2->{"begin"} > $breakpoint) {
	    push @newgenes2, $gene2;
	}
    }
    $run2->{"genes"} = \@newgenes2;
}
