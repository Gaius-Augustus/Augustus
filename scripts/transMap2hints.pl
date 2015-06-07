#!/usr/bin/perl
#
# Mario Stanke, 26.01.2007
# convert a transmap output file (tab format) to a hints file
# include only transmap hints from the list of ids in idfile
# alignment file format:
# name   chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds 
# this script is quicker when you sort the input file by chromosomes and by gene position
# cat transmap.gp | sort -n -k 4,4 | sort -s -k 2,2 > transmap.sorted.gp

use strict;
use Getopt::Long;

my $usage = "$0 -- convert transmap alignments to hints file\n";
$usage .= "  include only transmap hints from the list of ids in idfile\n";
$usage .= "  alignment file format:\n";
$usage .= "  name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\n";
$usage .= "  this script is quicker when you sort the input file by chromosomes and by gene position, e.g. with\n";
$usage .= "  cat transmap.gp | sort -n -k 4,4 | sort -s -k 2,2 > transmap.sorted.gp\n";
$usage .= "\n";
$usage .= "Usage: $0 --in=transmapfile.psl --out=hints.gff\n";
$usage .= "  Options:\n";
$usage .= "  --keepids=idfile       If idfile is given, then only those genes are used.\n";
$usage .= "  --priority=n           larger number = higher priority (default 4)\n";
$usage .= "  --ep_cutoff=n          (default 1)\n";
$usage .= "  --ep_margin=n          (default 18)\n";
$usage .= "                         each exon suggested by the input yields 3 exonpart hints:\n";
$usage .= "                          |             exon as suggested by transmap          |\n";
$usage .= "                                     |  ep   |    ep       |   ep  |\n";
$usage .= "                          |ep_cutoff |       |   exon      |       | ep_cutoff |\n";
$usage .= "                          |    ep_margin     |   exon      |    ep_margin      |\n";
$usage .= "  --ip_cutoff=n          the end of the intronpart interval are shorter by n bp than the intron (default 0)\n";
$usage .= "  --utrend_cutoff=n      UTRpart hint are cut off at the end by n bp (default 15)\n";
$usage .= "  --min_intron_len=n     minimal intron len (default 50)\n";
$usage .= "  --min_intron_len_utr=n a gap in the utr must be at least this long to be counted as an intron (default 80)\n";
$usage .= "  --start_stop_radius=n  (start and stop codon hint interval size)/2 (default 15)\n";
$usage .= "  --tss_tts_radius=n     (tss and tts interval size)/2 (default 100)\n";
$usage .= "  --source=s             source identifier in output (default 'T')\n";


my $ep_cutoff = 1;        # | ep_cutoff |       |   exon      |       | ep_cutoff |
my $ep_margin = 18;#was 12# |     ep_margin     |   exon      |    ep_margin      |
my $utrend_cutoff = 15;   # in addition to $ep_cutoff
my $ip_cutoff = 0;        # this is cut off at both ends of a suggested intron
my $min_intron_len = 50;
my $max_intronpart_len = 200000; # some transmap introns are insanely long and mess everything up
my $min_intron_len_utr = 80;
my $min_exon_len = 3;
my $priority = 4;
my $start_stop_radius = 15;
my $tss_tts_radius = 100;
my $source = 'T';
my $line=0;
my $have_idfile=0;
my ($idfilename, $hintsfilename, $transmapfilename);
my $prgname='t2h';

if ($#ARGV < 1 ) {
    print "$usage";
    exit;
}

GetOptions('in=s'=>\$transmapfilename,
	   'out=s'=>\$hintsfilename,
	   'keepids:s'=>\$idfilename,
	   'priority:i'=>\$priority,
	   'ep_cutoff:i'=>\$ep_cutoff,
	   'ep_margin:i'=>\$ep_margin,
	   'ip_cutoff:i'=>\$ip_cutoff,
	   'utrend_cutoff:i'=>\$utrend_cutoff,
	   'min_intron_len:i'=>\$min_intron_len,
	   'min_intron_len_utr:i'=>\$min_intron_len_utr,
	   'start_stop_radius:i'=>\$start_stop_radius,
	   'tss_tts_radius:i'=>\$tss_tts_radius,
	   'source:s'=>\$source);

if (defined $idfilename) {
    $have_idfile = 1;
}

open(TRANSMAP, "<$transmapfilename") || die "Couldn't open $transmapfilename\n";
if ($have_idfile) {
    open(IDS, "<$idfilename") || die "Couldn't open $idfilename\n";
}
open(HINTS, ">$hintsfilename") || die "Could not open $hintsfilename";

my (@starthints, @stophints, @tsshints, @ttshints, @dsshints, @asshints, @exonhints, @exonparthints, @intronhints, @intronparthints, @CDSparthints, @UTRparthints);
my ($strand,$qname,$txStart, $txEnd, $cdsStart, $cdsEnd, $exonstarts, $exonends, $unknown, $leftcmpl, $rightcmpl, $frames, $hasIntrons);
my (@f,@s,@e,@in);
#my (@blockbegins, @blockends);
my $numBlocks;
my ($from, $to, $ifrom, $ito, $i);
# hint lists are sorted by increasing begin position
my @hint; # (begin, end, strand, tname, qname)
my $hintref;
my ($targetname, $oldtargetname);
$oldtargetname = "no name yet";
my %goodids=();
my $id;
my $ephintbegin;
my $ephintend;

if ($have_idfile) {
    while(<IDS>){
	$id = $_;
	chop $id;
	$goodids{$id}=1;
    }
}

while (<TRANSMAP>) {
    $line++;
    s/#.*//;
    next unless /\S/;
    if ($line % 10000 == 0){
	print "line $line\n";
    }

    @f = split /\t/, $_, 16;
    if (@f < 15) { warn "Incorrect format"; next } # alignment of mRNA

    $qname = $f[0];      # mRNA       
    $targetname = $f[1]; # chromosome
    $strand = $f[2];
    $txStart = $f[3];
    $txEnd = $f[4];
    $cdsStart = $f[5];
    $cdsEnd = $f[6];
    # $f[7] not used (number of exons?)
    $exonstarts = $f[8];
    $exonends = $f[9];
    $unknown = $f[10];
    $leftcmpl = $f[12];
    $rightcmpl = $f[13];
    $frames = $f[14];
    $hasIntrons = $f[15];

    # filter alignments
    my $qnameuniq = $qname;
    $qnameuniq =~ s/-\d+$//; # transMap adds -1 -2 -3 etc to the diffent alignments of the same query
    if ($have_idfile && !exists $goodids{$qnameuniq}) { 
	next;
    }
    if ($targetname ne $oldtargetname) {
	printHints();
    }

    # now add the hints
    $exonstarts =~ s/[, ]$//;
    $exonends =~ s/[, ]$//;
    $hasIntrons =~ s/[, ]$//;
    @s = split /,/, $exonstarts;
    @e = split /,/, $exonends;
    if ($hasIntrons ne ""){
	@in = split /,/, $hasIntrons;
    } else {
	# not all transMap output has the intron information about gaps
	# assume that all of the gaps within size are introns
	@in = (1) x (@s-1);
    }

    # print "exonstarts ", (join ", ", @s), " exonends ", (join ", ", @e) , " hasintrons ", (join ", ", @in) , "\n";

    $numBlocks = scalar @s;
    # start and stop hint
    if ($strand eq '+') {
	if ($txStart != $cdsStart && $cdsStart > 0) {
	    @hint = ($cdsStart+1-$start_stop_radius, $cdsStart+3+$start_stop_radius, '+', $qname);
	    addSignalHint(\@starthints, [@hint]);
	}
	if ($txEnd != $cdsEnd && $cdsEnd > 0) {
	    @hint = ($cdsEnd-2-$start_stop_radius, $cdsEnd+$start_stop_radius, '+', $qname);
	    addSignalHint(\@stophints, [@hint]);
	}
    } else{
	if ($txStart != $cdsStart && $cdsStart > 0) {
	    @hint = ($cdsStart+1-$start_stop_radius, $cdsStart+3+$start_stop_radius, '-', $qname);
	    addSignalHint(\@stophints, [@hint]);
	}
	if ($txEnd != $cdsEnd && $cdsEnd > 0) {
	    @hint = ($cdsEnd-2-$start_stop_radius, $cdsEnd+$start_stop_radius, '-', $qname);
	    addSignalHint(\@starthints, [@hint]);
	}
    }
    # transcription start and transcription stop hint (tss and tts)
    if ($strand eq '+') {
	if ($txStart != $cdsStart && $txStart > 0) {
	    @hint = ($txStart+1-$tss_tts_radius, $txStart+1+$tss_tts_radius, '+', $qname);
	    addSignalHint(\@tsshints, [@hint]);
	}
	if ($txEnd != $cdsEnd && $txEnd > 0) {
	    @hint = ($txEnd-$tss_tts_radius, $txEnd+$tss_tts_radius, '+', $qname);
	    addSignalHint(\@ttshints, [@hint]);
	}
    } else{
	if ($txStart != $cdsStart && $txStart > 0) {
	    @hint = ($txStart+1-$tss_tts_radius, $txStart+1+$tss_tts_radius, '-', $qname);
	    addSignalHint(\@ttshints, [@hint]);
	}
	if ($txEnd != $cdsEnd && $txEnd > 0) {
	    @hint = ($txEnd-$tss_tts_radius, $txEnd+$tss_tts_radius, '-', $qname);
	    addSignalHint(\@tsshints, [@hint]);
	}
    }

    $ephintbegin=-1;
    $ephintend=-1;
    for ($i=0; $i<$numBlocks; $i++) {
	$from = $s[$i]+1;
	$to   = $e[$i];
	if ($i==0) {
	    if ($from + $utrend_cutoff <= $to) { 
		$from += $utrend_cutoff;
	    } else {
		$from = $to;
	    } 
	    if ($from > $cdsStart && $cdsStart >= $txStart ) {
                $from = $cdsStart;
            }
	}
	if ($i == $numBlocks-1) {
	    if ($to - $utrend_cutoff >= $from) {
		$to -= $utrend_cutoff;
	    } else {
		$to = $from;
	    } 
	    if ($to < $cdsEnd && $cdsEnd <= $txEnd) {
                $to = $cdsEnd;
            }
	}
	if ($ephintbegin<0 || $ephintend <0) {
	    $ephintbegin = $from;
	    $ephintend = $to;
	} elsif ((($ephintend < $cdsStart || $ephintbegin>$cdsEnd) && ($ephintend + $min_intron_len_utr + 1 >= $from))||
		 ($ephintend + $min_intron_len + 1 >= $from)){
	    $ephintend = $to;
	} else { # large gap 
	    $ifrom = $ephintend+1;
	    $ito = $from-1;
	    if ($ito-$ifrom+1 >= $min_intron_len && $in[$i-1]) {
		@hint = ($ifrom, $ito, $strand, $qname);
		addIntervalHint(\@intronhints, [@hint]);
		# also add dss and ass hints in case it is an utr intron. those intron hints don't help
		if ($ifrom < $cdsStart || $ifrom  > $cdsEnd ) {
		    @hint = ($ifrom, $ifrom, $strand, $qname);
		    if ($strand eq '+') {
			addSignalHint(\@dsshints, [@hint]);
		    } else {
			addSignalHint(\@asshints, [@hint]);
		    }
		}
		if ($ito < $cdsStart || $ito  > $cdsEnd ) {
		    @hint = ($ito, $ito, $strand, $qname);
		    if ($strand eq '+') {
			addSignalHint(\@asshints, [@hint]);
		    } else {
			addSignalHint(\@dsshints, [@hint]);
		    }
		}
		$ifrom += $ip_cutoff;
		$ito -= $ip_cutoff;
		if ($ifrom < $ito && $ifrom > $cdsStart && $ito < $cdsEnd && $ito-$ifrom+1 <= $max_intronpart_len) {
		    @hint = ($ifrom, $ito, $strand, $qname);
		    addIntervalHint(\@intronparthints, [@hint]);
		}
	    }
	    addExonpartFuzzyHint($ephintbegin, $ephintend, $strand, $qname);
	    $ephintbegin = $from;
	    $ephintend = $to;
	}
    }
    addExonpartFuzzyHint($ephintbegin, $ephintend, $strand, $qname);
    $ephintbegin = $from;
    $ephintend = $to;

    $oldtargetname = $targetname;
}

printHints();


###########################################################################################
#
# subroutines
#
###########################################################################################

#
# printHints
# print and delete the hints
#
sub printHints {
    # finished computing the hints for the old sequence output them
    # delete all exonpart hints that are contained in an exon hint
    my $startidx=0;
    my $curidx;
    foreach my $exon (@exonhints){
	my $start = $exon->[0];
	my $end = $exon->[1];
	my $strand = $exon->[2];
	# increase $startidx until the exonpart does not start to the left of start
	while ($startidx <= $#exonparthints && @exonparthints[$startidx]->[0] < $start){
	    $startidx++;
	}
	$curidx = $startidx;
	while ($curidx <= $#exonparthints && @exonparthints[$curidx]->[0] <= $end){
	    if (@exonparthints[$curidx]->[0] >= $start && @exonparthints[$curidx]->[1] <= $end && @exonparthints[$curidx]->[2] eq $strand ) {
		#redundant, delete it
		#print "deleting " , (join " ", @{$exonparthints[$curidx]}), " as it is contained in " , (join " ", @{$exon}), "\n";
		splice @exonparthints, $curidx, 1;
	    } else {
		$curidx++;
	    }
	}
    }	
    foreach $hintref (@tsshints) {
	print HINTS "$oldtargetname\t$prgname\ttss\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t0\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
     foreach $hintref (@starthints) {
	print HINTS "$oldtargetname\t$prgname\tstart\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t0\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
    foreach $hintref (@stophints) {
	print HINTS "$oldtargetname\t$prgname\tstop\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t0\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
    foreach $hintref (@ttshints) {
	print HINTS "$oldtargetname\t$prgname\ttts\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t0\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    }
    foreach $hintref (@asshints) {
	print HINTS "$oldtargetname\t$prgname\tass\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
    foreach $hintref (@dsshints) {
	print HINTS "$oldtargetname\t$prgname\tdss\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
    foreach $hintref (@exonparthints) {
	print HINTS "$oldtargetname\t$prgname\texonpart\t$hintref->[0]\t$hintref->[1]\t$hintref->[4]\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    }
    foreach $hintref (@exonhints) {
	print HINTS "$oldtargetname\t$prgname\texon\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    }
    foreach $hintref (@intronhints) {
	print HINTS "$oldtargetname\t$prgname\tintron\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    } 
#    foreach $hintref (@intronparthints) {
#	print HINTS "$oldtargetname\t$prgname\tintronpart\t$hintref->[0]\t$hintref->[1]\t0\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
#   }
    foreach $hintref (@CDSparthints) {
	print HINTS "$oldtargetname\t$prgname\tCDSpart\t$hintref->[0]\t$hintref->[1]\t$hintref->[4]\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    }
    foreach $hintref (@UTRparthints) {
	print HINTS "$oldtargetname\t$prgname\tUTRpart\t$hintref->[0]\t$hintref->[1]\t$hintref->[4]\t$hintref->[2]\t.\tgrp=$hintref->[3];src=$source;pri=$priority\n";
    }
    
    # delete all hints as the new sequence starts
    @starthints = ();
    @stophints = ();
    @tsshints = ();
    @ttshints = ();
    @dsshints = ();
    @asshints = ();
    @exonhints = ();
    @exonparthints = ();
    @intronhints = ();
    @intronparthints = ();
    @CDSparthints = ();
    @UTRparthints = ();
}

# addExonpartFuzzyHint(begin, end, strand, qname)
# split the interval in case it covers the CDS start or stop
# and then add the more specific CDS/UTR hints
#
sub addExonpartFuzzyHint {
    my $begin = shift;
    my $end = shift;
    my $strand = shift;
    my $name = shift;
   
    if ($end-$begin+1 < $min_exon_len) {
	return;
    }
    if ($cdsStart > $begin && $cdsStart <= $end){
	addFuzzyHint($begin, $cdsStart-1, $strand, $name);
	$begin = $cdsStart;
    }
    if ($cdsEnd >= $begin && $cdsEnd < $end){
	addFuzzyHint($begin, $cdsEnd, $strand, $name);
	$begin = $cdsEnd+1;
    }
    addFuzzyHint($begin, $end, $strand, $name);
}

# addFuzzyHint(begin, end, strand, qname)
#
#      | ep_cutoff |       |       exon      |       | ep_cutoff |
#      |     ep_margin     |       exon      |    ep_margin      |
#               fuzzybegin                         fuzzyend
#                       corebegin        coreend
#
sub addFuzzyHint {
    my $begin = shift;
    my $end = shift;
    my $strand = shift;
    my $name = shift;
    my ($corebegin, $coreend, $fuzzybegin, $fuzzyend);

    if ($begin>$end) {return;}
    $fuzzybegin = $begin + $ep_cutoff;
    $fuzzyend   = $end - $ep_cutoff;
    $corebegin  = $begin + $ep_margin;
    $coreend    = $end - $ep_margin;
    if ($corebegin > $coreend) {
	$corebegin = $coreend = int (($corebegin + $coreend)/2);
    }
    @hint = ($corebegin, $coreend, $strand, $qname, 2); # main part in the middle, score 2
    addExonpartHint([@hint]);
    if ($fuzzybegin < $corebegin) {
	@hint = ($fuzzybegin, $corebegin-1, $strand, $qname, 1); # left end piece, score 1
	addExonpartHint([@hint]);
    } 
    if ($fuzzyend > $coreend) {
	@hint = ($coreend+1, $fuzzyend, $strand, $qname, 1); # right end piece, score 1
	addExonpartHint([@hint]);
    }
}

#
# addExonpartHint(hintref)
# search in the list of existing exonpart hints for an including or included one
# if no such hint exists, sort the parameter hint into this list
# if there is a stronger one, do nothing, if there are weaker ones, replace them with this one
#

sub addExonpartHint {
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    my $k;
    my $typeref;
    
    if ($end < $cdsStart || $begin > $cdsEnd) {
	$typeref = \@UTRparthints;
    } elsif ($begin >= $cdsStart && $end <= $cdsEnd) {
	$typeref = \@CDSparthints;
    } else {
	$typeref = \@exonparthints;
    }

    # insert hint at the right position	
    $k = scalar(@{$typeref})-1;
    while ($k>=0 && $typeref->[$k]->[0] > $begin) {
	$k--;
    }
    my @temparray = ($href);
    if ($k == scalar(@{$typeref})-1) {
	push @{$typeref}, @temparray;
    } else {
	splice (@{$typeref}, $k+1, 0, @temparray);
    }
}

#
# addSignalHint(hintref)
# add the hint if it is not already there

sub addSignalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $strand = $href->[2];
    #print (join " ", @{$href});

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1) {
	push @{$hintlistref}, $href;
    } else {
	my $k = @{$hintlistref}-1;
	while ($k>=0 && $hintlistref->[$k]->[0]>= $begin) {
	    $k--;
	}
	my @temparray = ($href);
	if (!(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0]== $begin && $hintlistref->[$k+1]->[2] eq $strand) ||
	    ($k+2 <= @{$hintlistref}-1 && $hintlistref->[$k+2]->[0]== $begin && $hintlistref->[$k+2]->[2] eq $strand))) {
	    # hint does not previously exist, insert it
	    if ($k== @{$hintlistref}-1) {
		push @{$hintlistref}, @temparray;
	    } else {
		splice (@{$hintlistref}, $k+1, 0, @temparray);
	    }
	}
    }
}

#
# addIntervalHint(hintref)
# for exon and intron hints (not exonpart)
# add the hint if it is not already there

sub addIntervalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    #print (join " ", @{$href});
    #print "\n";

    push @{$hintlistref}, $href;
    return;

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1) {
	push @{$hintlistref}, $href;
    } else {
	my $k = @{$hintlistref}-1;
	while ($k>=0 && ($hintlistref->[$k]->[0]> $begin || ($hintlistref->[$k]->[0]== $begin && $hintlistref->[$k]->[1]>= $end))) {
	    $k--;
	}
	my @temparray = ($href);
	if ( 1
	    #!(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0]== $begin && $hintlistref->[$k+1]->[1]== $end && $hintlistref->[$k+1]->[2] eq $strand) ||
	    #($k+2 <= @{$hintlistref}-1 && $hintlistref->[$k+2]->[0]== $begin && $hintlistref->[$k+2]->[1]== $end && $hintlistref->[$k+2]->[2] eq $strand))
	    ) {
	    # hint does not previously exist, insert it
	    if ($k== @{$hintlistref}-1) {
		push @{$hintlistref}, @temparray;
	    } else {
		splice (@{$hintlistref}, $k+1, 0, @temparray);
	    }
	}
    }
}

