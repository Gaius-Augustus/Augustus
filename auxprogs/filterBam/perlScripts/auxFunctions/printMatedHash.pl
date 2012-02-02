#!/usr/bin/perl
# 
# Prints hash with mates found
#
# Created: 14-November-2011
# Last modified: 19-January-2012
#

# Package definition
use strict;
use 5.010; # so we can use say()

# Prints contents of Hash with mates found
sub printMatedHash()
{
	my $refmated = shift();
	print STDERR "Printing mated summary:\n";
	my %mated = %{$refmated};
	 foreach my $key (sort { $a <=> $b } keys %mated) 
	 {
	 	print STDERR "mate:$key, mated: $mated{$key} times\n";
	 }

}

return 1;
