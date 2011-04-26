#!/usr/bin/perl
#
# Convert a polyA polyT position table to tts hints
# This is for the polyA/T tables from Chun Liang
#
# Mario Stanke, 11.7.2009

use strict;
use Getopt::Long;

my $usage = "$0 -- polyA/polyT table to tts hints for AUGUSTUS\n";
$usage .= "\n";
$usage .= "Usage: $0 --in=polyA.txt --out=hintsfile\n";
$usage .= "  options:\n";
$usage .= "  hintradius=n       a 3' transcript end at position p gives rise to a tts hint from p-n to p+n (default 10)\n";
$usage .= "  format=new|newer   new: input file in newer format with 12 columns\n";
$usage .= "                     newer: input file in format with these 13 columns\n";
$usage .= "                     SeqName EstDir  Chromsome  ChrStrand ChrStart ChrEnd  EstStrand EstStart EstEnd  EstMapLen ChrSite Type PolyLen PolyStr\n";
$usage .= "  swapstrand         swap strand, currently use this for polyT files\n";


my $infilename;
my $hintsfilename;
my $hintradius = 10;
my $source="E";
my $format = "old";
my $swapstrand = 0;
my $priority = 4;
my $prgsrc = "polyA";
my ($qname, $strand, $targetname, $ttspos, $start, $end);

if ($#ARGV < 1 ) {
    print "$usage";
    exit;
}

GetOptions(
    'in=s'=>\$infilename,
    'out=s'=>\$hintsfilename,
    'format=s'=>\$format,
    'swapstrand!'=>\$swapstrand,
    'hintradius:i'=>\$hintradius);

open(POLYA, "<$infilename") || die "Couldn't open $infilename\n";
open(HINTS, ">$hintsfilename") || die "Could not open $hintsfilename";

while (<POLYA>) {
    next if (/Chromosome/ && /SeqName/); # skip header line
    my @f = split /\t/, $_;
    if (($format eq "old" && @f < 26) || ($format eq "new" && @f<12) || ($format eq "newer" && @f<13)) { warn "Not poly(A) format"; next } 
    
    if ($format eq "new"){
	$targetname  = $f[1];
	$ttspos      = $f[9];
	$strand      = $f[2];
    } elsif ($format eq "newer"){
        $targetname  = $f[2];
        $ttspos      = $f[10];
	$strand      = $f[3];
    } else {
	$targetname  = $f[3];
	$ttspos      = $f[18];
	$strand      = $f[2];
    }
    $qname       = $f[0];

    $start       = $ttspos - $hintradius;
    $end         = $ttspos + $hintradius;
    if ($swapstrand) {
	if ($strand eq "+") {
	    $strand = "-";
	} else {
	    $strand = "+";
	}
    }
    print HINTS "$targetname\t$prgsrc\ttts\t$start\t$end\t.\t$strand\t.\tpri=$priority;src=$source;est=$qname\n";
}
