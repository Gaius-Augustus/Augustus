#!/usr/bin/perl
# 
# Prints contents of the hash RSCORES, which has as key value a QALI array
#
# Created: 6-January-2012
# Last modified: 9-January-2012
#

# Package definition
use strict;
use 5.010; # so we can use say()

# Prints contents of Hash with mates found
sub printRscores()
{
	my $ref = shift();
	print "key=>value\n";
	my %hash = %{$ref};
	 foreach my $key (sort { $a <=> $b } keys %hash) 
	 {
	 	print "$key=>$hash{$key}\n";
	 }

}

return 1;


