#!/usr/bin/perl
# 
# Prints contents of hash array
#
# Created: 22-January-2012
# Last modified: 22-January-2012
#

# Package definition
use strict;
use 5.010; # so we can use say()
require '/home/tonatiuh/web-home/scripts/auxFunctions/printAoA.pl';

# Prints contents of Hash with mates found
sub printPairCovSteps()
{
	my $ref = shift();
	print STDERR "Printing pairCovSteps contents: key=>[value(s)]\n";
	my %hash = %{$ref};
	 foreach my $key (sort { $a <=> $b } keys %hash) 
	 {
	 	print STDERR "$key=>"; 
		printAoA(@{$hash{$key}});
	 }

}

return 1;
