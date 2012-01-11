#!/usr/bin/perl

# Evaluate a psl-file that contains peptide 2 protein mappings 
# in context with a gff-file of the protein coding genes 
# the peptides were mapped against.
#
# Assumptions:
# The gtf file was sorted by startCoord column, e.g.:
# cat in.gtf | sort -n -k4 > out.gtf
# The gtf file contains in its last column an AUGUSTUS transcript
# descriptor with the following format: 
# transcript_id "au(\d*).g(\d*).t(\d*)
#
# Example call:
# perl peptides2hints.pl tiny.psl tiny.gtf E

my $usage = "peptides2hints.pl psl-file gff-file src > hint-file\n";

if(@ARGV!=3){print STDERR $usage; exit -1;}

my $pslFile = $ARGV[0];
my $gffFile = $ARGV[1];
my $src = $ARGV[2];

my $gffL;
my %gffHash = ();
my $geneName;

# read gff file into hash of arrays
open(GFF, "<", $gffFile) or die ("Could not open gff-file $gffFile!\n");
while(<GFF>){
	if(($_=~m/CDS/) or( $_=~m/intron/)){
		$gffL = $_;
		$gffL =~ m/transcript_id "au(\d*).g(\d*).t(\d*)"/;
		$geneName = "au$1.g$2.t$3";
		push (@{$gffHash{$geneName}}, $gffL);

	}

}
close(GFF) or die ("Could not close gff-file $gffFile!\n");

# compute gene lengths and delete genes whose length is not a multiple of 3
my %geneLengths = ();
my $thisGeneLen = 0;
my @gffLine = ();
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
	if($thisGeneLen%3 == 0){
		$geneLengths{$geneName} = $thisGeneLen;
	}else{
		delete $gffHash{$geneName};
	}
}

my @pslLine = ();
my $pslStart;
my $pslStop;
my $gffStart;
my $gffStop;
my $strand;
my $globalPslStart;
my $globalPslStop;
my $protLen;
my @tmpPslLine;
open(PSL, "<", $pslFile) or die("Could not open psl-file $pslFile\n");

# each line of PSL is processed
while(<PSL>){
	chomp;
	@pslLine = split(/\t/);
	# continue working with perfect matches, only
	if($pslLine[0] == $pslLine[10]){
		$protLen = $pslLine[14];
		#print "Local $pslStart $pslStop\n";
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
				$pslLine[15] = $protLen - $tmpPslLine[16] + 1;
				$pslLine[16] = $protLen - $tmpPslLine[15] + 1;
			}
			$pslStart = $pslLine[15];
			$pslStop = $pslLine[16];


			$exonC = 1;
			$intronLen = 0;
			print "# $geneName $strand\n";
			foreach(@geneGff){
				chomp;
				@gffLine = split(/\t/);
				$gffStart = $gffLine[3];
				$gffStop = $gffLine[4];
				if($strand eq "+"){
					if($gffLine[2]=~m/CDS/){
						if($exonC==1){
							$firstExonStart = $gffStart;
						}
						#print "-----\n".$_."\n";
						$globalPslStart = $pslStart*3-2 + $firstExonStart - 1 + $intronLen;
						$globalPslStop = $pslStop*3 + $firstExonStart -1 + $intronLen;
						#print $globalPslStart." ".$globalPslStop."\n";
						if($globalPslStart>=$gffStart and $globalPslStart<=$gffStop and $globalPslStop>=$gffStart and $globalPslStop<=$gffStop){	
							print $gffLine[0]."\tpep2hints\tCDSpart\t".$globalPslStart."\t".$globalPslStop."\t.\t$strand\t0\tsrc=$src\n";	
							#print "Exon no. $exonC: GlobalPslStart: $globalPslStart GlobalPslStop: $globalPslStop vs. gffStart $gffStart and gffStop $gffStop\n";
						}elsif($globalPslStart>=$gffStart and $globalPslStart<=$gffStop and $globalPslStop>=$gffStop){
							#print "Partial hint found!\n";
							print $gffLine[0]."\tpep2hints\tCDSpart\t".$globalPslStart."\t".$gffStop."\t.\t$strand\t0\tsrc=$src\n";
							#print "Exon no. $exonC: LocalPslStart: $pslStart GlobalPslStart: $globalPslStart vs. gffStart $gffStart and gffStop $gffStop\n";
							$printIntron = 1;
						}elsif($globaPslStart<=$gffStart and $globalPslStop>=$pslStart and $globalPslStop<=$gffStop and $printIntron==1){
							#print "Partial hint end found!\n";
							print $gffLine[0]."\tpep2hints\tCDSpart\t".$gffStart."\t".$globalPslStop."\t.\t$strand\t".$gffLine[7]."\tsrc=$src end of intronspanning crap\n";
							$printIntron = 0;
						}
						$exonC++;
	
					}elsif($gffLine[2]=~m/intron/){
						$intronLen = $intronLen + ($gffStop -$gffStart + 1);
						if($printIntron==1){
							print $gffLine[0]."\tpep2hints\tintron\t$gffStart\t$gffStop\t.\t$strand\t.\tsrc=$src\n";
						}
						
					}
	
				}else{
	
				}
			}
		}
	}
}
close(PSL) or die("Could not close psl-file $pslFile\n");
