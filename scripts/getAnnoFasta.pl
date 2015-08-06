#!/usr/bin/perl
# getAnnoFast.pl
# Creates fasta sequence files from the AUGUSTUS output.
# 
# Mario Stanke, 10.05.2007
#

use strict;
use Getopt::Long;

my $usage = "getAnnoFasta.pl augustus.gff\n";
$usage .= "   Makes a fasta file with protein sequences (augustus.aa)\n";
$usage .= "   and one with coding sequences (augustus.codingseq)\n";
$usage .= "   from the sequences provided in the comments of the AUGUSTUS output.\n";
$usage .= "   These sequence comments are turned on with --protein=on and --codingseq=on, respectively\n";
$usage .= "Options:\n";
$usage .= "   --seqfile=s  Input a fasta file with the genomic sequences that AUGUSTUS was run on.\n";
$usage .= "                When this option is given, an additional file with the individual\n";
$usage .= "                coding exon sequences (augustus.cdsexons) is output.\n";
$usage .= "                and a file with the complete mRNA including UTRs (augustus.mrna) is output.\n";
$usage .= "   --chop_cds   for incomplete transcripts: cut off bases before first codon.\n";

my ($seqname, $trid, $status, $haveCod, $haveAA, $haveCDS, $haveRNA, $seq, $seqfile, $chop_cds);

GetOptions('seqfile=s'=>\$seqfile, 'chop_cds!'=>\$chop_cds);

if ($#ARGV != 0) {
    print $usage;
    exit;
} 

my $separator = ";";
my $augustusfilename = $ARGV[0];
open(AUG, "<$augustusfilename") || die "Couldn't open $augustusfilename\n";
my $stemfilename = $augustusfilename;
$stemfilename =~ s/(\.gff|\.gtf|.gff3|\.txt)//;
my %sequence = (); # Hash with the DNA sequence. The sequence names are the keys.

# Read in the sequence file in one chunk.
# And sort it in the sequence hash.
# Yes, this requires a lot of memory for large genomes.
if ($seqfile){
    open (SEQ, "<$seqfile") or die ("Could not open sequence file $seqfile\n");
    $/=">";
    while (<SEQ>){
	s/>$//;
	next unless /\S+/;
	/(.*)\n/;
	
	$seqname = $1;
	my $sequencepart = $'; #'
	$seqname =~ s/\s.*//; # seqname only up to first white space
	$sequencepart =~ s/\s//g;
	$sequence{$seqname} = $sequencepart;
    }
    print "Read in " . (scalar keys %sequence) . " sequence(s) from $seqfile.\n";
}
$/="\n";

#
# Go through the augustus output transcript by transcript.
#
$haveCod = $haveAA = $haveCDS = $haveRNA = 0;
$status = 0;
my $exonUTRFormat = 0; # UTR implicitly given by exon features
my $UTRFormat = 0; # UTR explicitly given by *UTR features
my $cdsSeq = "";
my $aaSeq = "";
my %cdsnr = ();
my %cdsTx = (); #keys transcript id, values concatenated coding exon sequences
my %mrnaTx = (); #keys transcript id, values concatenated exon sequences (including UTR)
my %strandTx = (); #keys transcript id, values strands
my %frameTx = (); #keys transcript id, values frames

while(<AUG>) {
    if ($seqfile && (/^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\ttranscript_id "([^"]*)"; gene_id "([^"]*)";$/
                 || /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*Parent=([^;]+)/)){
	my $feat = $2;
	$seqname = $1;
	my $start = $3;
	my $end = $4;
	my $strand = $5;
	my $frame = $6;
	$trid=$7;
	$trid =~ s/\s$//;
	next unless ($feat eq "CDS" || $feat =~ /UTR/ || $feat eq "exon");
	# decide whether to use exon or UTR format for mRNA by whether we see UTR or exon first
	$UTRFormat = 1 if (!$exonUTRFormat && $feat =~ /UTR/);
	$exonUTRFormat = 1 if (!$UTRFormat && $feat eq "exon");
	$cdsnr{$trid}++ if ($feat eq "CDS");
	my $seqpart = lc(substr($sequence{$seqname}, $start-1, $end - $start + 1));
	#print "$seqname $trid CDS $cdsnr{$trid} $start -> $end $seqpart\n";
	if ($seqpart ne "") {
	    # add mRNA if applicable
	    if (($exonUTRFormat && $feat eq "exon") ||
		($UTRFormat && ($feat eq "CDS" || $feat =~ /UTR/))) {
		$mrnaTx{$trid} = "" if (!defined($mrnaTx{$trid}));
		$mrnaTx{$trid} .= $seqpart;
	    }
	    if ($feat eq "CDS"){
		if (!$haveCDS) {
		    open (CDSEXON, ">$stemfilename.cdsexons");
		    $haveCDS++;
		}
		$cdsTx{$trid} = "" if (!defined($cdsTx{$trid}));
		$cdsTx{$trid} .= $seqpart;
		push(@{$frameTx{$trid}}, $frame);
		if ($strand eq '-') {
		    $seqpart = rc($seqpart);
		    $strandTx{$trid} = $strand;
		}
		print CDSEXON ">$trid.cds" . $cdsnr{$trid} . "\n$seqpart\n";
	    }
	}
    }
    if (/^(\S+)\t.*\ttranscript_id "([^"]*)"; gene_id "([^"]*)";$/ ||
	/^(\S+)\t.*Parent=([^;]+)/){
	$seqname=$1;
	$trid=$2;
	$trid =~ s/\s$//;
	$status=1;
    } elsif (/coding sequence = \[(.*)/ && $status == 1){
	if ($haveCod == 0) {
	    open (COD, ">$stemfilename.codingseq");
	}
	$haveCod++;
	$seq = $1;
	$seq =~ s/\]$//;
	print COD ">$seqname.$trid\n$seq\n";
	$status=2;
    } elsif ($status == 2 && /^\# ([\w\]]*)$/){
	$seq = $1;
	$seq =~ s/\]$//;
	print COD "$seq\n";
	$status=2;
    } elsif (/protein sequence = \[(.*)/ && $status >= 1){
	if ($haveAA == 0) {
	    open (AA, ">$stemfilename.aa");
	}
	$haveAA++;
	$seq = $1;
	$seq =~ s/\]$//;
#	print AA ">$seqname$separator$trid\n$seq\n";
#	print AA ">$trid\n";
	$aaSeq .= $seq;
	if (!/\]/){
	    $status=3;
	} else {
	    if ($aaSeq ne ""){
 		print AA ">$trid\n";
		print AA getFa($aaSeq, 100);
	    }
	    $aaSeq = "";
	    $status=1;
	}
    } elsif ($status == 3 && /^\# (.*)/){
	$seq = $1;
	$seq =~ s/\]$//;
#	print AA "$seq\n";
	$aaSeq .= $seq;
	if (!/\]/){
	    $status=3;
	} else {
	    if ($aaSeq ne ""){
		print AA ">$trid\n";
		print AA getFa($aaSeq, 100);
	    }
	    $aaSeq = "";
	    $status=1;
	}
    }
}

#
# print coding sequences, if not already done (because included in output)
#
if (!$haveCod && scalar(keys %cdsTx)>0){
    open (COD, ">$stemfilename.codingseq") or die ("Could not open $stemfilename.codingseq for writing.");
    foreach my $trid (sort by_id keys %cdsTx){
	print COD ">$trid\n";
	my $codingseq = $cdsTx{$trid};
	$codingseq = rc($codingseq)  if ($strandTx{$trid} eq "-");
	if($chop_cds){
	    if ($strandTx{$trid} eq "-"){
		if($frameTx{$trid}[-1] != 0){
		    $codingseq = substr($codingseq, $frameTx{$trid}[-1]);
		}
	    }else{
		if($frameTx{$trid}[0] != 0){
		    $codingseq = substr($codingseq, $frameTx{$trid}[0]);
		}
	    } 
	}
	print COD getFa($codingseq);
    }
}

#
# print mRNA sequences
#
if (scalar(keys %mrnaTx)>0){
    open (MRNA, ">$stemfilename.mrna") or die ("Could not open $stemfilename.mrna for writing.");
    foreach my $trid (sort by_id keys %mrnaTx){
	print MRNA ">$trid\n";
	my $mrnaseq = $mrnaTx{$trid};
	$mrnaseq = rc($mrnaseq) if ($strandTx{$trid} eq "-");
	print MRNA getFa($mrnaseq);
    }
}


#
# sort by increasing transcript id
#
sub by_id{
    $a =~ /g(\d+)\.t(\d+)/;
    my ($ag,$at)=($1,$2);
    $b =~ /g(\d+)\.t(\d+)/;
    my ($bg,$bt)=($1,$2);
    if ($ag>$bg){
	return 1;
    } elsif ($bg>$ag){
	return -1;
    } else {
	return $at <=> $bt;
    }
}


# reverse complement
sub rc{
    my $s = shift;
    $s = reverse $s;
    $s =~ s/a/T/g;
    $s =~ s/c/G/g;
    $s =~ s/g/C/g;
    $s =~ s/t/A/g;
    $s = lc $s;
    return $s;
}


sub getFa{
    my $seq = shift;
    my $cols = 100;
    $cols = shift if (@_);
    
    my $start = 0;
    my $ret = "";
    while (length($seq)-$start >= $cols) {
	my $shortline = substr($seq, $start, $cols);
	$ret .= "$shortline\n";
	$start += $cols;
    }
    $ret .= substr($seq, $start, $cols) . "\n" if ($start<length($seq));
    return $ret;
}
