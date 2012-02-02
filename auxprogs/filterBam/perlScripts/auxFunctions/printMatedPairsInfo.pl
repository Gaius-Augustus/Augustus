#!/usr/bin/perl
# 
# Prints mated pairs
#
# Created: 14-November-2011
# Last modified: 19-January-2012
#

# Package definition
use strict;
use 5.010; # so we can use say()

# Prints contents of MatePairs
sub printMatedPairsInfo
{
	#printMatedPairsInfo(\@matepairs, \@qali);
	my ($refMatePairs, $refQali) = @_; 
	my $counter = 0;
	my @matepairs = @$refMatePairs;
	my @qali = @$refQali;
    my @mate1 = (); 
	my @mate2 = (); 

	print STDERR "Printing Mate-pairs info:\n";
	for my $mp (@matepairs) 
	 {
		@mate1 = split /\t/, $qali[@$mp[0]]->[0], 21;
		@mate2 = split /\t/, $qali[@$mp[1]]->[0], 21;
	 	print STDERR "(@$mp[0],@$mp[1])=($mate1[9],$mate2[9]),scoreMate=@$mp[2]\n";
	 	$counter++
	 }

}

return 1;
