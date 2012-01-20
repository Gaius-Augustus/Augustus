#!/usr/bin/perl

# Creates a fasta file with unique peptides from a fasta input file
# where identical peptide match multiple positions. The first found 
# peptide occurence is used for counting multiplicity.
# Input format fasta header examples:
# >scaffold4083_450 [15866 - 17011]  (65.81)
# >au14.g511.t1 (46.17)
# ... there is a score in the end of the header in round brackets!
# Output format fasta header:
# >peptideX mult=x
# Do not use for long protein sequences, only for short peptides!!!
# prefix

my $usage = "uniquePeptides.pl in.fa prefix > out.fa\n\nRead script head comments for futher documentation!\n\n";

if(@ARGV!=2){print STDERR $usage; exit -1;}

my $in = $ARGV[0];
my $prefix = $ARGV[1];

# two hashes for the same peptide (hash of hashes would be more elegant, but I am too lazy for that), key is always the peptide sequence, trimmed of "-" characters that sometimes appear to occur in Harald's peptide files.
my %locHash = ();
my %multHash = ();
my $maxLen = 100;

my $header;

open(IN, "<", $in) or die("Could not open file $in!\n");

while ( <IN> ) { 
	chomp;
	if($_=~m/^>/){
		$_=~s/>//;
		$_=~s/\(\d+\.\d+\)//; # delete score
		$header = $_;
	}else{
		if(length($_) > $maxLen){
			print STDERR "The peptide $_ is longer $maxLen. Aborting.\n";
			exit -1;
		}
		while($_=~m/-/){ # delete weird dashes in sequences
			$_ = s/-//;
		}
		if(not(exists($locHash{$_}))){
			$locHash{$_} = "$header";
			$multHash{$_} = 1;			
		}elsif($locHash{$_} eq $header){
			$multHash{$_} = $multHash{$_} + 1;
		}
		
		
	}
}

close IN or die("Could not close file $in!\n");

# print unique fasta file
my $c = 1;
while ( ($k,$v) = each %multHash ) {
	print ">$prefix$c mult=$v\n$k\n";
	$c++;
}


# print mapping file
#open(OUT2, ">", $mapping) or die("Could not close file $mapping!\n");
#for $k ( keys %mapHash ) {
#	print OUT2 "$k:\n";
#	foreach(@{$mapHash{$k}}){
#		print OUT2 "$_\n";
#	}
#}
#close(OUT2) or die("Could not close file $mapping!\n");

