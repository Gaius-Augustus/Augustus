#!/usr/bin/perl
#
# generate hints from an Exonerate gff dump
# Mario Stanke June, 2007

use strict;
use Getopt::Long;

my $usage = "$0 -- convert exonerate output file to hints file for AUGUSTUS\n";
$usage .= "\n";
$usage .= "Usage: $0 --in=exfile --out=hintsfile\n";
$usage .= "  Exonerate run like this: exonerate --model protein2genome --showtargetgff T ... > exfile\n";
$usage .= "  options:\n";
$usage .= "  --priority=n       priority of hint group (default 4)\n";
$usage .= "  --minintronlen=n   alignments with gaps shorter than this and longer than maxgaplen are discarded (default 41)\n";
$usage .= "  --maxintronlen=n   alignments with longer gaps are discarded (default 350000\n";
$usage .= "  --CDSpart_cutoff=n this many bp are cut off of each CDSpart hint w.r.t. the exonerate cds (default 9)\n";
$usage .= "  --source=s         source identifier (default 'P')\n";
#$usage .= "  --ssOn             include splice site (dss, ass) hints in output (default false)\n";

my $exfile;
my $hintsfilename;
my $minintronlen = 41;
my $maxintronlen = 350000;
my $CDSpart_cutoff = 15;
my $source="XNT";
my $priority = 4;
my $prgsrc = "xnt2h";
my $line=0;
my $coloffset=0;
my $ssOn=0;
my $CDSpartid = "CDSpart"; # abbreviate to decrease file size
my $prot="";

if ($#ARGV < 1 ) {
    print "$usage";
    exit;
}

GetOptions(
	   'in=s'=>\$exfile,
	   'out=s'=>\$hintsfilename,
	   'minintronlen:i'=>\$minintronlen,
	   'maxintronlen:i'=>\$maxintronlen,
	   'CDSpart_cutoff:i'=>\$CDSpart_cutoff,
	   'source:s'=>\$source,
	   'priority:i'=>\$priority,
	   'ssOn!'=>\$ssOn);

open(XNT, "<$exfile") || die "Couldn't open $exfile\n";
open(HINTS, ">$hintsfilename") || die "Could not open $hintsfilename";


while (<XNT>) {
    s/#.*//;
    next unless /\S/;
    next unless /\texonerate:protein2genome:local\t/; # modified by Katharina on September 8th 2010 because exonerate output format apparently changed regarding the string "local".q
   # print "I am in the file!\n";
    my @f = split /\t/, $_, 9;
    if (@f < 8) { warn "Not gff format"; next }
    my $seqname = $f[0];
    my $type = $f[2];
    my $start = $f[3];
    my $end = $f[4];
    my $score = $f[5];
    my $strand = $f[6];
#    if ($strand eq '-') {
#	$start++; # exonerate has a weird coordinate system
#	$end++;
#    } # apparently this has changed with exonerate versions and it is now standard
    if ($end < $start) {
	my $tmp = $start;
	$start = $end;
	$end = $tmp;
    }
    if ($type eq "intron") {
	# length in range?
	if ($end - $start + 1 >= $minintronlen && $end - $start + 1 <= $maxintronlen){
	    print HINTS "$seqname\t$prgsrc\tintron\t$start\t$end\t$score\t$strand\t.\tsrc=$source;grp=$prot;pri=$priority\n";
	}
    } elsif ($type eq "cds") {
	$start += $CDSpart_cutoff;
	$end -= $CDSpart_cutoff;
	if ($start > $end) {
	    $start = $end = int(($start+$end)/2);
	}	
	print HINTS "$seqname\t$prgsrc\t$CDSpartid\t$start\t$end\t$score\t$strand\t.\tsrc=$source;grp=$prot;pri=$priority\n";
    } elsif ($type eq "gene") {
	/sequence (\S+) ; /;
	$prot = $1;
    }
}
