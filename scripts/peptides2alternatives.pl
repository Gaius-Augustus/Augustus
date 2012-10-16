#!/usr/bin/perl

# Katharina J. Hoff, Jan 30th 2012
#
# Scenario:
# A list of peptides was identified by using a library of
# alternatively spliced gene predictions from AUGUSTUS.
# You want to map those peptides back to the genome via
# mapping to the gene amino acid sequences and placing
# their positions in the genome. Since alternative 
# transcripts are redundant, you cannot simply map against
# the complete library of predicted genes, but you have
# to map iteratively, i.e. remove peptides that have already 
# a hit in one of the alternative transcript isoforms from
# the list before proceeding to the next alternative
# transcript isoform. 
#
# Input files must be a in multiple fasta format.
# Output is a psl file.

my $usage = "peptides2alternatives.pl peptides.aa augustus-genes.aa > iteratively-mapped.psl\n";

if(@ARGV!=2){print STDERR $usage; exit -1;}

my $peptides = $ARGV[0];
my $allGenes = $ARGV[1];
my $toBeMapped = "tobemapped.fa"; # file with peptides that are still to be mapped, must be writable!
my $target = "target.fa"; # file with proteins to be mapped against, must be writable!
my $tmpPsl = "tmp.psl"; # temporary psl file, must be writable!
my $blat = "/usr/local/bin/blat"; # BLAT version on your system
my @t1; # temporary variable 1

print STDERR "Find number of transcripts\n";
# identify max number of transcripts
my $maxT = 0;
open(GENES, "<", $allGenes) or die "Could not open file $allGenes!\n";
while(<GENES>){
	if($_ =~ m/^>/){
		chomp;
		$_ =~ m/^>g\d+\.t(\d+)/;
		if($1 > $maxT){$maxT = $1;}
	}
}
close(GENES) or die "Could not close file $allGenes!\n";

print STEDERR "Read peptides for mapping\n";
# initialize toBeMapped Hash
my %mapHash = ();
my $curPep;
open(PEPTIDES, "<", $peptides) or die "Could not open peptide file $peptides!\n";
while(<PEPTIDES>){
	chomp;
	if($_=~m/^>/){
		$_=~s/^>//;
		@t1 = split(/ /, $_);
		$curPep = $t1[0];
	}else{
		$mapHash{$curPep} = $_;
	}
}
close(PEPTIDES) or die "Could not close peptide file $peptides!\n";

# for each transcriptNo Do
my $printFlag = 0;
my $count;
print STDERR "Looping from 1 to $maxT with mapping\n";
for ($count=1; $count<=$maxT; $count++){
	# get transcripts that are going to be used as mapping target
	open(GENES, "<", $allGenes) or die "Could not open file $allGenes!\n";
	open(TARGET, ">", $target) or die "Could not open file $target!\n";
	while(<GENES>){
		if($_=~m/^>/){
			if($_=~m/^>g\d+\.t$count\n/){$printFlag = 1;}else{$printFlag = 0}
		}
		if($printFlag == 1){
			print TARGET $_;
		}
	}
	close(TARGET) or die "Could not close file $target!\n";
	close(GENES) or die "Could not close file $allGenes!\n";
	# write toBeMapped peptides
	open(TOBEMAPPED, ">", $toBeMapped) or die "Could not open to be mapped file $toBeMapped!\n";
	while (($k, $v) = each %mapHash){
		print TOBEMAPPED ">$k\n$v\n";
	}
	close(TOBEMAPPED) or die "Could not close to be mapped file $toBeMapped!\n";
	print STDERR "blat process no $count... ";
	## map pulled transcripts
	`$blat -t=prot -q=prot -noHead $target $toBeMapped $tmpPsl`;
	## read in psl file, delete already mapped peptides, print psl
	open(PSL, "<", $tmpPsl) or die "Could not open psl file $tmpPsl!\n";
	while(<PSL>){
		@t1 = split(/\t/);
		delete $mapHash{$t1[9]};
		print $_;
	}
	close(PSL) or die "Could not close psl file $tmpPsl!\n";
	`rm $toBeMapped`;
}
print STDERR "Done!\n";
`rm $target $tmpPsl`;
