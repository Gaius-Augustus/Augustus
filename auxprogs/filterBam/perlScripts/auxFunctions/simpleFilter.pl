#!/usr/bin/perl
#
# Auxiliary variable estimates the "basic" filters that follow the following criteria:
# a) Coverage
# b) Percentage identity
# c) Intron size
#
# Created:  8-November-2011
# Last modified: 9-November-2011
#

# Package definition
use strict;

sub simpleFilter
{
	my $fRef = shift;
	my $paramsRef = shift;
	# De-reference the array list
	my @f = @$fRef;
	my @params = @$paramsRef;
	my $minintronlen = $params[0];
	my $maxintronlen = $params[1];
	my $minId        = $params[2];
	my $outMinId     = $params[3];
	my $minCover     = $params[4];
	my $outMinCover  = $params[5];
	my $nointrons    = $params[6];
	my $outIntrons   = $params[7];
    my $match       = $f[0];  # No. of bases that match that aren't repeats
    my $qBaseInsert = $f[5];  # No. of bases inserted in query
    my $TgapCount   = $f[6];  # No. of inserts in target
    my $tBaseInsert = $f[7];  # No. of bases inserted in target
    my $strand      = $f[8];  # +, - for query strand. 2nd +,- for genomic strand
    my $qname       = $f[9];  # Query name, e.g. SRR027108.21/1
    my $qsize       = $f[10]; # Query sequence size
    my $qstart      = $f[11]; # Alignment start position in query
    my $qend        = $f[12]; # Alignment end position in query
    my $targetname  = $f[13]; # Target sequence name
    my $tstart      = $f[15]; # Alignment start position in target
    my $tend        = $f[16]; # Alignment end position in target
    my $blockSizes  = $f[18]; # List of sizes of each block
    my $qStarts     = $f[19]; # List of starting positions of each block in query
    my $tStarts     = $f[20]; # List of starting positions of each block in target
    $blockSizes =~ s/[, ]$//;    # take away commas at the end of $BLOCKSIZES
    $tStarts =~ s/[, ]$//;       # takes away commas at the end of $TSTARTS
    my @b = split /,/, $blockSizes; # split $BLOCKSIZES, if several of them
    my @t = split /,/, $tStarts;    # split $TSTARTS, if several of them
    my @q = split /,/, $qStarts;    # split $QSTARTS, if several of them

	# filter for minimum percentage of identity
    my $tgap = 0; # inserted bases in target sequence, excluding introns
    my $qgap = 0; # inserted bases in query sequence
    my $gaps = 0; # total amount of gaps
    
    # @b contains $BLOCKSIZES
	# @t contains $TSTARTS
	# @q contains $QSTARTS

    for (my $i=0; $i<@b-1; $i++)
	  {
		$tgap = $t[$i+1]-$t[$i]-$b[$i];
		$qgap = $q[$i+1]-$q[$i]-$b[$i];
		$tgap = 0 if ($qgap ==0 && $tgap >= $minintronlen && $tgap <= $maxintronlen); #target gap is intron
		$gaps += ($tgap>$qgap)?  $tgap : $qgap; # count the larger gap if both seqs happen to have a gap
      }
	
	 # filter for percent identity
     my $percid = sprintf("%.1f", 100*$match/((($qend - $qstart + $gaps)))); 
     # filter for minimum coverage
     my $coverage =  sprintf("%.1f", 100*($qend - $qstart)/$qsize);

	return ($percid, $coverage);
	
} # end of function

return 1;
