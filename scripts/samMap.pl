#!/usr/bin/perl
#
# Map a sam file that was created by aligning against an exon-exon junction file back to genome level.
# Filters out non-intron-containing alignments.
#
# Katharina J. Hoff, May 3rd 2012

#use strict;
#use Getopt::Long;



my $usage = "\n$0 -- Map a sam file that was created by aligning against an exon-exon junction file back to genome level\n\n";
$usage .= "Usage: $0 exex.sam map.psl > out.sam\n\n";
$usage .= "exex.sam was created e.g. with GSNAP in the following way:\n";
$usage .= "gmap_build -d exex_db exex.fa\n";
$usage .= "gsnap --format=sam --nofails -d exex_db rnaseq.fastq > exex.sam\n\n";
$usage .= "map.psl was created in the following way:\n";
$usage .= "intron2exex.pl --introns=introns.lst --seq=genome.masked.fa --exex=exex.fa --map=map.psl\n\n";

if(@ARGV != 2){
    print $usage;
    exit;
}

my $exexFile = $ARGV[0];
my $mapFile = $ARGV[1];

my %mapHash = ();
my @t;
# read map file into hash, lookup by exexID
open(MAP, "<", $mapFile) or die("Could not open mapFile $mapFile!\n");
while(<MAP>){
	@t = split(/\t/);
	if(not(exists($mapHash{$t[9]}))){
		$mapHash{$t[9]} = $_;
	}else{
		print STDERR "Hash element $t[9] in mapHash already exists!\n";
	}
}
close(MAP) or die("Could not close mapFile $mapFile!\n");

my @tSam;
my $readLen;
my $exexName;
my $mapLine;
my @tMap;
my $junctionStart;
my $beforeJunction;
my $globalStart;
my $cigarLen;
my $cigPos;
my $safeForAfterJunctionNum;
my $safeForAfterJunctionDesc;
my $seenJunction;
my $intronLen;
my $samLen;
my $letter;
my $number;
my $newCigar;
my $startCoord;
# process each line of sam format
open(SAM, "<", $exexFile) or die("Could not open exexFile $exexFile!\n");
while(<SAM>){
    if(not(m/^@/)){       
	@tSam = split(/\t/);
	$readLen = length($tSam[9]);
	$exexName = $tSam[2];
	$mapLine = $mapHash{$exexName};
	@tMap = split(/\t/, $mapLine);
	$junctionStart = $tMap[15]+$tMap[0]/2+1; # first base after junction
	$globalStart = $tSam[3]+$tMap[15];
	$intronLen = $tMap[7];
	# process the cigar string
	my $softclipped = 0;	
	if($tSam[5]=~m/^(\d+)S/){
		$softclipped = $1;		
	}
	$beforeJunction = $junctionStart - $globalStart + $softclipped;
	if($beforeJunction > 0){
		$startCoord = $tSam[3]+$tMap[15]  ;

	}else{
		$startCoord = ($tSam[3]+$tMap[15]+$tMap[7]);
	}
	$cigarLen = length($tSam[5]);
	$safeForAfterJunction = 0;
	$cigPos = 1;
	$seenJunction = 0;
	$newCigar = "";
	while($cigarLen > 0){
		$tSam[5] =~ m/^(\d+)(\w)/;
		$letter = $2; 
		$number = $1;
		if($tSam[3]>($tMap[0]/2)){
		    $seenJunction = 1;
		}
		if($number < $beforeJunction and $seenJunction == 0){
			$newCigar = $newCigar.$number.$letter;		
			$beforeJunction = $beforeJunction - $number;
		}elsif($seenJunction == 0 and $beforeJunction => 1 ){
			$newCigar = $newCigar.($number - ($number - $beforeJunction));	
			$newCigar = $newCigar.$letter;
			$newCigar = $newCigar.$intronLen."N";
			$newCigar = $newCigar.($number-$beforeJunction);
			$newCigar = $newCigar.$letter;
			$seenJunction = 1;
		}else{
			$newCigar = $newCigar.$number.$letter;
		}
	
		$cigPos = $cigPos + $number;
		$tSam[5]=~ s/$number$letter//;
		$cigarLen = length($tSam[5]);
	}
	if($newCigar =~ m/N/){
		print $tSam[0]."\t".$tSam[1]."\t".$tMap[13]."\t".$startCoord."\t".$tSam[4]."\t$newCigar";
		$samLen = @tSam;
		for ($count = 6; $count < $samLen; $count ++){
			print "\t$tSam[$count]";
		}
	}
    }
}
close(SAM) or die("Could not close exexFile $exexFile!\n");
