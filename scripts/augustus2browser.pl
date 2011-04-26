#!/usr/bin/perl
#
# This takes the AUGUSTUS output in the standard input
# and outputs to standard output a file with UCSC browser gtf format
# The display includes tracks for those types of hints that are present in the file.
# Mario 01.08.2006
#
use strict;

my @augustus = <STDIN>;

my $firstseqname="";
my $firstseqbegin=-1;
my $lastseqend=-1;
my ($seqname, $begin, $end, $type, $prgname);
#
# determine the browser position entry
# display the first sequence in the full range of annotation
#
my %hints = ('ass' => 0, 'dss' => 0, 'tss' => 0, 'start' => 0, 'stop' => 0, 'tts' => 0, 'intronpart' => 0, 'intron' => 0,
	     'UTRpart' => 0, 'CDSpart' => 0, 'nonexonpart' => 0, 'irpart' => 0, 'CDS' => 0, 'UTR' => 0, 'exon' => 0, 'exonpart' => 0);

foreach my $line (@augustus) {
    if ($line !~/^#/) {
	# expand shortcuts for exonpart, intronpart, CDSpart, ..
	$line =~ s/\tep\t/\texonpart\t/;
	$line =~ s/\tip\t/\tintronpart\t/;
	$line =~ s/\tcp\t/\tCDSpart\t/;
	$line =~ s/\tup\t/\tUTRpart\t/;
	$line =~ s/\tnep\t/\tnonexonpart\t/;
	
	if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/) {
	    $seqname = $1;
	    $prgname = $2;
	    $type = $3;
	    $begin = $4;
	    $end = $5;
	    if (exists $hints{$type} && !($prgname =~ /^AUGUSTUS/)) {
		$hints{$type} += 1;
	    }
	    if ($firstseqname eq "") {
		$firstseqname = $seqname;
	    }
	    if ($seqname eq $firstseqname) {
		if ($firstseqbegin < 0 || $lastseqend < 0) {
		    $firstseqbegin = $begin;
		    $lastseqend = $end;
		}
		if ($begin < $firstseqbegin) {
		    $firstseqbegin = $begin;
		}
		if ($end > $lastseqend) {
		    $lastseqend = $end;
		}
	    }
	}
    }
}

my $browserbegin = $firstseqbegin - 1000;
my $browserend =   $lastseqend + 1000;
if ($browserbegin < 1) {
    $browserbegin = 1;
}

#print "browser position $firstseqname:$browserbegin-$browserend\n";
my $grpNameCount = 1;
my $grpName;
foreach my $feature (keys %hints) {
    next unless $hints{$feature}>0;
    print "track name=\"$feature\" description=\"$feature hints\" visibility=3\n";
    foreach my $line (@augustus) {
	if ($line !~/^#/) {
	    if ($line =~ /\t$feature\t/ && !($line =~/\tAUGUSTUS/)) {
		$line =~ s/;?source=[^;]*;?//;
		$line =~ s/;?src=[^;]*;?//;
		$line =~ s/grp=//;
		$line =~ s/pri=\d+;?//;
		$line =~ s/priority=\d+//;
		$line =~ s/\t([^\t]*)$//;
		$grpName = $1;
		$grpName =~ s/\s*\n//;
		if ($grpName eq "" || ($grpName =~ /mult=\d+;/)){ # no name, use counter as name so the browser does not complain about same group different strand
		    $grpName = "" . ($grpNameCount++) . "$grpName";
		}
		print $line, "\t" , $grpName, "\n";
	    }
	}
    }	
}

my $haveAUG=0;
foreach my $line (@augustus) {
#    $line =~ /\tAUGUSTUS(\d+)\t/;
#    if ($1 ne $id) {
#	$id = $1;
#	print "track name=AUGUSTUS$id description=\"AUGUSTUS $id\" visibility=3\n";
#    }
    if ($line =~ /\tAUGUSTUS\S*\texon\t/ || $line =~ /\tAUGUSTUS\S*\tCDS\t/) {
	print "track name=AUGUSTUS description=\"AUGUSTUS\" visibility=3\n" if (!$haveAUG);
	print $line;
	$haveAUG = 1;
    }
}

