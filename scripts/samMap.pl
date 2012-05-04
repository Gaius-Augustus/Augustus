#!/usr/bin/perl
#
# Map a sam file that was created by aligning against an exon-exon junction file back to genome level
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
	@tSam = split(/\t/);
	$readLen = length($tSam[9]);
	$exexName = $tSam[2];
	# print $exexName."\t";
	$mapLine = $mapHash{$exexName};
	@tMap = split(/\t/, $mapLine);
	$intronLen = $tMap[7];
	$junctionStart = $tMap[15]+$tMap[0]/2; # last base before junction
	
	$globalStart = ($tSam[3]+$tMap[15]-1);
	#
	# process the cigar string

	$beforeJunction = $junctionStart - $globalStart +1;
	if($beforeJunction > 0){
		$startCoord = ($tSam[3]+$tMap[15]-1+2);
	}else{
		$startCoord = ($tSam[3]+$tMap[15]-1+$tMap[7]);
	}
	$cigarLen = length($tSam[5]);
	$safeForAfterJunction = 0;
	$cigPos = 1;
	$seenJunction = 0;
	$newCigar = "";
	# print "\n\noldCigar: $tSam[5]\nnewCigar computation:\nThe junction is at $beforeJunction in the read\n";
	while($cigarLen > 0){
		#print "Cigar Len at the beginning is $cigarLen\n";
		$tSam[5] =~ m/^(\d+)(\w)/;
		$letter = $2; 
		$number = $1;
		# print "We are now dealing with $number $letter\n";
		if($number < $beforeJunction and $seenJunction == 0){
			#print "This claims that $number + $cigPos -1 is smaller or equal to $beforeJunction.\n";
			# print "Truely before Junction at $beforeJunction: ";
			# print "$number$letter\n";
			$newCigar = $newCigar.$number.$letter;		
			# print "Remaining positions before the junction: ".($beforeJunction - $number)." Reset beforeJunction to this value.\n";
			$beforeJunction = $beforeJunction - $number;
		}elsif($seenJunction == 0 and $beforeJunction > 1 ){
			# print "Partially before Junction: ";
			# print ($number - ($number - $beforeJunction)-1);
			$newCigar = $newCigar.($number - ($number - $beforeJunction)-1);	
			# print "$letter\n";
			$newCigar = $newCigar.$letter;
			# print " at junction: ";
			# print "$intronLen"."N\n";
			$newCigar = $newCigar.$intronLen."N";
			# print " partially after Junction:";
			# print ($number-$beforeJunction+1);
			$newCigar = $newCigar.($number-$beforeJunction+1);
			# print "$letter\n";
			$newCigar = $newCigar.$letter;
			$seenJunction = 1;
		}else{
			# print " truely after junction: ";
			# print "$number$letter\n";
			$newCigar = $newCigar.$number.$letter;
		}
	
		$cigPos = $cigPos + $number;
		#print "\ncigPos reset to $cigPos\n";
		$tSam[5]=~ s/$number$letter//;
		$cigarLen = length($tSam[5]);
		#print "Cigar Len at the end is $cigarLen\n";


	}
# 	print "New Cigar: $newCigar\n";
	if($newCigar =~ m/N/){
		print $tSam[0]."\t".$tSam[1]."\t".$tMap[13]."\t".$startCoord."\t".$tSam[4]."\t$newCigar";
		$samLen = @tSam;
		for ($count = 6; $count < $samLen; $count ++){
			print "\t$tSam[$count]";
		}
	}
}
close(SAM) or die("Could not close exexFile $exexFile!\n");
