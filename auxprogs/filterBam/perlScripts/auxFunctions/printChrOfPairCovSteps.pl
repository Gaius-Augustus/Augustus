# Auxiliary function that prints the pairedness coverage, 
# in terms of the coordinates covered by the reads, wrt 
# the targets
#
# Created: 28-January-2012
# Last modified: 28-January-2012

sub printChrOfPairCovSteps()
{
	my $ref = shift();
	my $refChr = shift();
	my $chr = $$refChr;
	print STDERR "Printing coverage of chr=$chr\n";
	my %hash = %{$ref};
	 foreach my $key (@{$hash{$chr}}) 
	 {
	 	print STDERR "@$key\n";
	 }
}

return 1;
