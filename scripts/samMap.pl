#!/usr/bin/perl
#
# Map a sam file that was created by aligning against an exon-exon junction file back to genome level.
# Filters out non-intron-containing alignments.
#
# Katharina J. Hoff, May 3rd 2012

#use strict;
#use Getopt::Long;



my $usage = "\n$0 -- Map a sam file that was created by aligning against an exon-exon junction file back to genome level\n\n";
$usage .= "Usage: $0 exex.sam map.psl [flank] > out.sam\n\n";
$usage .= "exex.sam was created e.g. with GSNAP in the following way:\n";
$usage .= "gmap_build -d exex_db exex.fa\n";
$usage .= "gsnap --format=sam --nofails -d exex_db rnaseq.fastq > exex.sam\n\n";
$usage .= "map.psl was created in the following way:\n";
$usage .= "intron2exex.pl --introns=introns.lst --seq=genome.masked.fa --exex=exex.fa --map=map.psl\n\n";
$usage .= "Please note that mapSam.pl assumes that the flanking region of introns is 75+75 bp! Specify other value if needed.\n";
if(! defined $ARGV[1]){
    print $usage;
    exit;
}

my $exexFile = $ARGV[0];
my $mapFile = $ARGV[1];
my $flank = 75;
if(defined $ARGV[2]){
    $flank = $ARGV[2];
}

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
my $cigarNumber;
my $badCigar = 0;
# process each line of sam format
open(SAM, "<", $exexFile) or die("Could not open exexFile $exexFile!\n");
while(<SAM>){
    if(not(m/^@/)){       
	@tSam = split(/\t/);
	$readLen = length($tSam[9]);
	$exexName = $tSam[2];
	$mapLine = $mapHash{$exexName};
	@tMap = split(/\t/, $mapLine);
	$junctionStart = $tMap[15]+$tMap[0]/2+1; # first base after
	if(($tMap[0]/2)==$flank){ # results from fragment border introns!
	    $globalStart = $tSam[3]+$tMap[15];
	    $intronLen = $tMap[7];
	    # process the cigar string
	    my $softclipped = 0;	
	    if($tSam[5]=~m/^(\d+)S/){
		$softclipped = $1;		
	    }
	    $beforeJunction = $junctionStart - $globalStart + $softclipped;
	    if($beforeJunction > 0){
		$startCoord = $tSam[3]+$tMap[15];
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
		if($number <= $beforeJunction and $seenJunction == 0){
		    $cigarNumber = $number;
		    if(($cigarNumber==0 and $letter eq 'M') or $cigarNumber=~m/\D/){
				$badCigar = 1;
		    }
			$newCigar = $newCigar.$cigarNumber.$letter;		
			$beforeJunction = $beforeJunction - $number;
		}elsif($seenJunction == 0 and $beforeJunction >= 1 ){
		    $cigarNumber = ($number - ($number - $beforeJunction));
		    if(($cigarNumber==0 and $letter eq 'M') or $cigarNumber=~m/\D/){
				$badCigar = 1;
		    }
			$newCigar = $newCigar.$cigarNumber;	
			$newCigar = $newCigar.$letter;
			$newCigar = $newCigar.$intronLen."N";
		    $cigarNumber = ($number-$beforeJunction);
                    if(($cigarNumber==0 and $letter eq 'M') or $cigarNumber=~m/\D/){
                        $badCigar = 1;
                    }
			$newCigar = $newCigar.$cigarNumber;
			$newCigar = $newCigar.$letter;
			$seenJunction = 1;
		}else{
		    $cigarNumber = $number;
            if(($cigarNumber==0 and $letter eq 'M') or $cigarNumber=~m/\D/){
				$badCigar = 1;
            }
		    $newCigar = $newCigar.$cigarNumber.$letter;
		}
	
		$cigPos = $cigPos + $number;
		$tSam[5]=~ s/$number$letter//;
		$cigarLen = length($tSam[5]);
	}
	if($newCigar =~ m/N/ and $badCigar==0 and not($newCigar=~m/\d+H\d+N\d+H/)){
		print $tSam[0]."\t".$tSam[1]."\t".$tMap[13]."\t".$startCoord."\t".$tSam[4]."\t$newCigar";
		$samLen = @tSam;
		for ($count = 6; $count < $samLen; $count ++){
			print "\t$tSam[$count]";
		}
	}elsif($badCigar==1 or $newCigar=~m/\d+H\d+N\d+H/){
	    print STDERR "Warning - Bad CIGAR (this alignment is not included in the output at stdout):\n".$tSam[0]."\t".$tSam[1]."\t".$tMap[13]."\t".$startCoord."\t".$tSam[4]."\t$newCigar\n";
                $samLen = @tSam;
                for ($count = 6; $count < $samLen; $count ++){
                        print STDERR "\t$tSam[$count]";
                }
	    print STDERR "Sam:\n";
	    print STDERR $_;
	    
	    print STDERR "Map:\n";
	    foreach(@tMap){
	    print STDERR "\t".$_;
	    }
	}
	$badCigar = 0;
    }
}
}
close(SAM) or die("Could not close exexFile $exexFile!\n");
