#!/usr/bin/perl
#
# convert a blat output file (tab format) for transcript alignments
# to a GBrowse GFF file and the multiple fasta file with the matching ests

use strict;
use Getopt::Long;

my $usage = "$0 -- convert blat file to gbrowse file\n";
$usage .= "\n";
$usage .= "Usage: $0 blat.psl gbrowse.gff\n";
$usage .= "Options:\n";
$usage .= "    --oldformat        output format for old GBrowse (before 2.0)\n";
$usage .= "    --estnames=file    output file with the names of the ESTs\n";
$usage .= "    --source=name      identifyier in the source column\n";
$usage .= "\n";

my $coloffset=0;
my $oldformat=0;
my $source;
my $estfilename;

if ($#ARGV <1) {
    die "Unknown option \n\n$usage";
}

GetOptions('estnames:s'=>\$estfilename,
	   'oldformat!'=>\$oldformat,
	   'source:s'=>\$source);

$source = "bl2gbr" unless (defined($source));

my $blatfilename = $ARGV[0];
my $hintsfilename = $ARGV[1];


open(BLAT, "<$blatfilename") || die "Couldn't open $blatfilename\n";
open(GFF, ">$hintsfilename") || die "Could not open $hintsfilename";
if (defined $estfilename){
    open(EST, ">$estfilename") || die "Could not open $estfilename";
}
my ($i, $j, $mstart, $mend);
my (@dsshints, @asshints, @exonhints, @exonparthints, @intronhints);
my ($match,$TgapCount,$strand,$qname,$qsize,$blockSizes,$tStarts, $qStarts);
my (@f,@b,@t,@q);
my (@blockbegins, @blockends);
my $numBlocks;
my $id=0;

# hint lists are sorted by by increasing begin position
my @hint; # (begin, end, strand, tname, qname)
my $hintref;
my $targetname;
my $skiplines=0;
my %estnames;

while (<BLAT>) {
    if (/psLayout/){
	$skiplines=5;
    }
    if ($skiplines>0) {
	$skiplines--;
	next;
    }
    s/#.*//;
    next unless /\S/;

    @f = split /\t/, $_, $coloffset+21;
    if (@f < $coloffset+20) { warn "Not BLAT format"; next } # blat format from the GenomeBrowser has an additional first column
    
    $match = $f[$coloffset+0];
    $TgapCount = $f[$coloffset+6];
    $strand = $f[$coloffset+8];
    $qname = $f[$coloffset+9];
    $qsize = $f[$coloffset+10];
    $targetname = $f[$coloffset+13];
    $blockSizes = $f[$coloffset+18];
    $qStarts = $f[$coloffset+19];
    $tStarts = $f[$coloffset+20];

    
    #print "match=", $match;
    #print " TgapCount=", $TgapCount;
    #print " strand=", $strand;
    #print " qname=", $qname;
    #print " blockSizes=", $blockSizes;
    #print " tStarts=", $tStarts, "\n";

    $blockSizes =~ s/[, ]$//;
    $tStarts =~ s/[, ]$//;
    @b = split /,/, $blockSizes;
    @t = split /,/, $tStarts;
    @q = split /,/, $qStarts;
    
    #print "blocksizes ", (join ", ", @b), " blockbegins ", (join ", ", @t) , "\n"; 
    
    $numBlocks = scalar @t;
    # Go throught the line
    #
    @blockbegins=();
    @blockends=();
    for ($i=0; $i<$numBlocks; $i++) {
	$mstart = $t[$i]+1; # blat is 0-based
	$mend = $mstart + $b[$i] - 1;

	push @blockbegins, $mstart;
	push @blockends, $mend;
    }
    $numBlocks = scalar @blockbegins;
    
    $id++;
    my $a, $b;
    if ($oldformat){
	if ($strand ne "-") {
	    print GFF "$targetname\t$source\tmatch\t",($t[0]+1),"\t", ($t[$numBlocks-1]+$b[$numBlocks-1]), "\t0\t$strand\t.\tTarget $source:$qname ", ($q[0]+1)," ",($q[$numBlocks-1]+$b[$numBlocks-1]),"\n";
	} else {
	    print GFF "$targetname\t$source\tmatch\t",($t[0]+1),"\t", ($t[$numBlocks-1]+$b[$numBlocks-1]), "\t0\t$strand\t.\tTarget $source:$qname ", ($qsize+1-($q[$numBlocks-1]+$b[$numBlocks-1]))," ",($qsize-$q[0]),"\n";
	}
    }
	
    for ($i=0; $i<$numBlocks; $i++) {
	if ($strand ne "-") {
	    if ($oldformat){
		print GFF "$targetname\t$source\tHSP\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t0\t$strand\t.\tTarget $source:$qname ", ($q[$i]+1), " ", ($q[$i]+$b[$i]),"\n";
	    } else {
		print GFF "$targetname\t$source\tHSP\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t0\t$strand\t.\tID=m$id;Name=$qname;Target=$qname ", ($q[$i]+1), " ", ($q[$i]+$b[$i]),"\n";	
	    }
	} else {
	    if ($oldformat){
		print GFF "$targetname\t$source\tHSP\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t0\t$strand\t.\tTarget $source:$qname ", ($qsize+1-($q[$i]+$b[$i])), " ",($qsize-$q[$i]) ,"\n";
	    } else {
		print GFF "$targetname\t$source\tHSP\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t0\t$strand\t.\tID=m$id;Name=$qname;Target=$qname ", ($qsize+1-($q[$i]+$b[$i])), " ",($qsize-$q[$i]) ,"\n";
	    }
	}
    }
    $estnames{$qname}++;
}

if (defined $estfilename) {
    foreach (keys %estnames) {
	print EST;
	print EST "\n";
    }
}
