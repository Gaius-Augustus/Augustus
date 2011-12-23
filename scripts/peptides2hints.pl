#!/usr/bin/perl

# Evaluate a psl-file that contains peptide 2 protein mappings 
# in context with a gff-file of the protein coding genes 
# the peptides were mapped against.
# it is assumed that the gtf file was sorted by startCoord column, e.g.:
# cat in.gtf | sort -n -k4 > out.gtf

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


my @pslLine = ();
my $pslStart;
my $pslStop;
my @geneGff = ();
my $gffStart;
my $gffStop;
my $strand;
my $globalPslStart;
my $globalPslStop;
my @gffLine = ();
open(PSL, "<", $pslFile) or die("Could not open psl-file $pslFile\n");
while(<PSL>){
	chomp;
	@pslLine = split(/\t/);
	# continue working with perfect matches, only
	if($pslLine[0] == $pslLine[10]){
		$pslStart = $pslLine[15];
		$pslStop = $pslLine[16];
		#print "Local $pslStart $pslStop\n";
		$geneName = $pslLine[13];
		# get gff-entry for gene
		@geneGff = @{$gffHash{$geneName}};
		if(@geneGff[0] =~ m/\t\+\t/){
			$strand = "+";
		}else{$strand="-";}
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
close(PSL) or die("Could not close psl-file $pslFile\n");
