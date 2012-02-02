#!/usr/bin/perl
# 
# Prints contents of hash array
#
# Created: 23-December-2011
# Last modified: 23-December-2011
#

# Package definition
use strict;
use 5.010; # so we can use say()

# Prints contents of Hash with mates found
sub printHash()
{
	my $ref = shift();
	print "key=>value\n";
	my %hash = %{$ref};
	 foreach my $key (sort { $a <=> $b } keys %hash) 
	 {
	 	print STDERR "$key=>$hash{$key}\n";
	 }

}

return 1;
