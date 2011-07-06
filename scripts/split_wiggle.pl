#!/usr/bin/perl

# split a wiggle file in single contig files

my $usage = "split_wiggle.pl in.wig outputPath prefix\n";

if (@ARGV != 2) {
    print $usage;
    exit;
}

$inFile = $ARGV[0];
$outPath = $ARGV[1];
if(defined($ARGV[2])){
  $prefix = $ARGV[2];
}else{$prefix="split_";}
$outFile = $outPath."/".$prefix;

open(IN, "<", $inFile) or die "Could not open file $inFile!\n";
$counter = 0;
while(<IN>){
  if($_=~m/track name/){$headerLine = $_;
  }elsif($_=~m/variableStep/){
     if($counter>=1){
	close(OUT) or die "Could not close File $thisFile!\n";
     }
	$counter = $counter + 1;
	$thisFile = $outFile.$counter;
	open(OUT, ">", $thisFile) or die "Could not open File $thisFile!\n";
	print OUT $headerLine;
	print OUT $_;
  }else{print OUT $_;}
}
close(OUT) or die "Could not close file $outFile\n";
close(IN) or die "Could not close file $inFile\n";
