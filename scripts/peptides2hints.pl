#!/usr/bin/perl

# Evaluate a psl-file that contains peptide 2 protein mappings 
# in context with a gff-file of the protein coding genes 
# or a six frame translation (gff!) that
# the peptides were mapped against.
#
# Assumptions:
# The gtf file was sorted by startCoord column, e.g.:
# cat in.gtf | sort -n -k4 > out.gtf
# The gtf file contains in its last column an AUGUSTUS transcript
# descriptor with the following format: 
# transcript_id "au(\d*).g(\d*).t(\d*)
#
# In a gene gtf file, only the features CDS and intron are used.
# In a six-frame translation gtf file, the feature frame is used.
#
# Example call:
# perl peptides2hints.pl tiny.psl tiny.gtf E
#
# A six-frame fasta file produced with getorf can be converted to gtf format
# with the following command:
# cat six-frame.fa | perl -ne 'if(m/^>/){$_=~s/>//; @t1 = split(/_/); print $t1[0]."\tGETORF\tframe\t"; $t1[1]=~m/\d+ \[(\d+) - (\d+)\]/; $staC = $1; $stoC = $2; if($_ =~ m/REVERSE/){print $stoC."\t".$staC."\t.\t-\t";}else{print $staC."\t".$stoC."\t.\t+\t";} @t2 = split(/ /, $_); print "0\t$t2[0]\n";}' > six-frame.gtf

my $usage = "peptides2hints.pl psl-file gff-file src > hint-file\n";

if(@ARGV!=3){print STDERR $usage; exit -1;}

my $pslFile = $ARGV[0];
my $gffFile = $ARGV[1];
my $src = $ARGV[2];

my $gffL;
my %gffHash = ();
my $geneName;
my @gffLine = ();

# read gff file into hash of arrays
open(GFF, "<", $gffFile) or die ("Could not open gff-file $gffFile!\n");
while(<GFF>){
	$gffL = $_;
	if(($_=~m/CDS/) or($_=~m/intron/)){
		$gffL =~ m/transcript_id "au(\d*).g(\d*).t(\d*)"/;
		$geneName = "au$1.g$2.t$3";
		push (@{$gffHash{$geneName}}, $gffL);
	}elsif($_=~m/frame/){
		# muss noch genename für six-Frame einlesen implementieren!
		@gffLine = split(/\t/, $gffL);
		$geneName = $gffLine[8];
		push (@{$gffHash{$geneName}}, $gffL);
	}
}
close(GFF) or die ("Could not close gff-file $gffFile!\n");

# compute gene lengths and delete genes whose length is not a multiple of 3
my $thisGeneLen = 0;
my @geneGff = ();
for $geneName ( keys %gffHash ) {
	$thisGeneLen = 0;
	@geneGff = @{$gffHash{$geneName}};
	foreach(@{geneGff}){
		@gffLine = split(/\t/);
		if($gffLine[2]=~m/CDS/){
			$thisGeneLen = $thisGeneLen + ($gffLine[4] - $gffLine[3] + 1);
		}
	}
	if(not($thisGeneLen%3 == 0)){
		delete $gffHash{$geneName};
	}
}

my @pslLine = ();
my $strand;
my $globalPslStart;
my $globalPslStop;
my @tmpPslLine;
open(PSL, "<", $pslFile) or die("Could not open psl-file $pslFile\n");

# each line of PSL is processed
while(<PSL>){
	chomp;
	@pslLine = split(/\t/);
	# continue working with perfect matches, only
	if($pslLine[0] == $pslLine[10]){
		$geneName = $pslLine[13];
		# get gff-entry for gene
		if(exists $gffHash{$geneName}){
			@geneGff = @{$gffHash{$geneName}};
			if(@geneGff[0] =~ m/\t\+\t/){
				$strand = "+";
			}else{
				$strand="-";
				# convert psl file to the opposite strand
				@tmpPslLine = @pslLine;
				$pslLine[15] = $pslLine[14] - $tmpPslLine[16] + 1;
				$pslLine[16] = $pslLine[14] - $tmpPslLine[15] + 1;
				# pslLine[14] = protein length/target length; pslLine[15] = psl start; pslLine[16] = psl stop
			}
			$exonC = 1;
			$intronLen = 0;
			#print "# $geneName $strand\n";
			foreach(@geneGff){
				chomp;
				@gffLine = split(/\t/);
				if($gffLine[2]=~m/CDS/){
					if($exonC==1){
						$firstExonStart = $gffLine[3]; # $gffLine[3] = gffStart; $gffStop = $gffLine[4];
					}
					$globalPslStart = $pslLine[15]*3-2 + $firstExonStart - 1 + $intronLen;
					$globalPslStop = $pslLine[16]*3 + $firstExonStart -1 + $intronLen;
					if($globalPslStart>=$gffLine[3] and $globalPslStart<=$gffLine[4] and $globalPslStop>=$gffLine[3] and $globalPslStop<=$gffLine[4]){	
						print $gffLine[0]."\tpep2hints\tCDSpart\t".$globalPslStart."\t".$globalPslStop."\t.\t$strand\t0\tsrc=$src\n";	
						$printIntron = 0;
					}elsif($globalPslStart>=$gffLine[3] and $globalPslStart<=$gffLine[4] and $globalPslStop>=$gffLine[4]){
						print $gffLine[0]."\tpep2hints\tCDSpart\t".$globalPslStart."\t".$gffLine[4]."\t.\t$strand\t";
						if($strand eq "+"){
							print "0";
						}else{
							print "$gffLine[7]";
						}
						print "\tsrc=$src\n";
						$printIntron = 1;
					}elsif($globaPslStart<=$gffLine[3] and $globalPslStop>=$pslLine[15] and $globalPslStop<=$gffLine[4] and $printIntron==1)
						print $gffLine[0]."\tpep2hints\tCDSpart\t".$gffLine[3]."\t".$globalPslStop."\t.\t$strand\t";
						if($strand eq "+"){
							print "$gffLine[7]";
						}else{
							print "0";
						}
						print "\tsrc=$src\n";
						$printIntron = 0;
					}
					$exonC++;
				}elsif($gffLine[2]=~m/intron/){
					$intronLen = $intronLen + ($gffLine[4] -$gffLine[3] + 1);
					if($printIntron==1){
						print $gffLine[0]."\tpep2hints\tintron\t$$gffLine[3]\t$gffLine[4]\t.\t$strand\t.\tsrc=$src\n";
					}
				}elsif($gffLine[2]=~m/frame/){
					$globalPslStart = $pslLine[15]*3-2 + $gffLine[3] - 1;
					$globalPslStop = $pslLine[16]*3 + $gffLine[3] -1;
					if($globalPslStart>=$gffLine[3] and $globalPslStart<=$gffLine[4] and $globalPslStop>=$gffLine[3] and $globalPslStop<=$gffLine[4]){	
						print $gffLine[0]."\tpep2hints\tCDSpart\t".$globalPslStart."\t".$globalPslStop."\t.\t$strand\t0\tsrc=$src\n";
					}
				}
			}
		}
	}
}
close(PSL) or die("Could not close psl-file $pslFile\n");
