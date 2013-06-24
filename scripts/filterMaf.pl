#!/usr/bin/perl
# filter blocks from a MAF alignment that intersect with a given genomic interval
# or filter sequences for a given subset of species
# Stefanie Koenig 23.06.2013

use strict;
use warnings;

use Getopt::Long; # for parameter specification on the command line

my $usage = <<'ENDUSAGE';
filterMaf.pl	filter blocks from a MAF alignment that intersect with a given genomic interval
		or filter sequences for a given subset of species
SYNOPSIS

filterMaf.pl < alignment.maf > filter.maf

OPTIONS

    --help             			output this help message
    --species=species1,species2,...	a comma separated list of the species which shall be kept in the filtered alignment (at least 2 species,
                                        by default all species are kept)		
    --min-seq N				only blocks that contain at list N sequences are taken (N=2 by default)
    --interval=seqID:start-end	        a genomic interval (seqID is equal to the first field of an 's' line in the alignment)

DESCRIPTION
      
  Example:

    filterMaf.pl --species=hg19,mm9,rheMac2,bosTau4 < alignment.maf > filter.maf
    filterMaf.pl --interval=hg19.chr21:43490000-43560000 < alignment.maf > filter.maf
    filterMaf.pl --min-seq 5 < alignment.maf > filter.maf
    filterMaf.pl --species=hg19,mm9,rheMac2,bosTau4 --min-seq 4 --interval=mm9.chr17:43490000-43560000 < alignment.maf > filter.maf

ENDUSAGE

my ($specieslist, $interval, $help); # options
my $minseq=2;

GetOptions('interval=s'=>\$interval,
	   'species=s'=>\$specieslist,
	   'min-seq:f' =>\$minseq,
           'help!'=>\$help);

if ($help){
    	print $usage;
    	exit(1);
}

my @slist=();
if (defined($specieslist)){
    @slist=split(/,/,$specieslist);
    if(@slist < 2){
	print "$specieslist has to be a comma separated list of at least two species.\n$usage";
    }
}


my ($s,$b,$e);
if (defined($interval)){
    if($interval=~m/(.+):(\d+)\-(\d+)/){
	$s=$1;
	$b=$2;
	$e=$3;
    }
    else{
	print "$interval has the wrong syntax.\n$usage";
    }
}

my @block=();
my $overlap = 0;

while(<>){
    if(/##maf/){
	print $_;
    }
    elsif(/^a\s/){
	if ((@block >= $minseq+1) && $overlap){
	    for (my $j=0; $j<@block; $j++){
		print $block[$j];
	    }
	    print "\n";
    	}
	@block=();
	push @block, $_;
	$overlap=0;
    }
    elsif (/^s\s/){
	my @f = split(/\s+/,$_);
	my ($seqID, $start, $alen, $strand, $slen) = ($f[1], $f[2], $f[3], $f[4], $f[5]);
	if($strand eq '-'){
	    $start= $slen - $start - $alen;
	}
	my $end=$start + $alen;
	$start++;
	if (!defined($specieslist)){
	    push @block, $_;
	}else{
	    foreach my $species (@slist){
		if($seqID=~m/$species/){
		    push @block, $_;
		    last;
		}
	    }
	}
	if (!defined($interval)){
	    $overlap=1;
	}
	else{
	    if(!($end<$b || $start>$e) && ($s eq $seqID)){
		$overlap=1;
	    }
	}
    }
}

if ((@block >= $minseq+1) && $overlap){
    for (my $j=0; $j<@block; $j++){
      print $block[$j];
  }
}
