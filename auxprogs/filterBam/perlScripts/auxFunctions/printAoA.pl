#!/usr/bin/perl
# 
# Auxiliary that prints an Array of Arrays
#
# Created: 14-October-2011
# Last modified: 16-January-2012
#

# Package definition
use strict;
use 5.010; # so we can use say()

# Prints contents of ArrayOfArrays
sub printAoA
{
	for my $aref (@_ ) 
	{
		print STDERR "\t [ @$aref ],\n";
		# say "\t [ @$aref ],";
	}
}

return 1;
