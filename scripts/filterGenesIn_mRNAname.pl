#!/usr/bin/perl

#############################################################
# filterGenes
# filter genes from a genbank flat file database
# usage: fileterGenes namefile dbfile
#
#
# Mario Stanke, 13.08.2002
#############################################################

if ($#ARGV != 1) {
    print "usage:filterGenes namefile dbfile\n";
    print "names of the loci to be kept come from\n";
    print "the first parameter. Only the the first of identical loci is kept\n";
    exit;
} 
$origfilename = $ARGV[1];
$goodfilename = $ARGV[0];

open(origfile, "$origfilename") || die "Couldn't open dbfile\n";
@data = <origfile>;
close(origfile);


$/="\n//\n";

my $head;
my $mRNAflag = 0;
my $cdsFlag = 0;
my $genename;
my $printFlag = 0;
my $firstPrintFlag = 0;


foreach(@data) {
    if($_=~m/^LOCUS/){
	$head = "";
	$printFlag = 0;
	$genename = "";
	$head = $head.$_;
    }
    if($_=~m/FEATURES/){
	$head = $head.$_;
    }
    if($_=~m/source/){
	$head = $head.$_;
    }
    if($mRNAflag==1 and not($_=~m/CDS/)){
	$head = $head.$_;
    }
    if($_=~m/mRNA/){
	$mRNAflag = 1;
	$head = $head.$_;
    }
    if($cdsFlag==1){
        if($_=~m/gene="/){
		@tmp = split(/\"/);
		$genename = $tmp[1];
		$cdsFlag = 0;
		$firstPrintFlag = 1;
        }else{
		$head = $head.$_;
        }
    }
    if($_=~m/CDS/){
	$mRNAflag = 0;
        $head = $head.$_;
	$cdsFlag = 1;
    }

    if($firstPrintFlag==1 and length($head)>=2){
	# correction by Simone Lange on November 6th 2014
	@line = `grep -w $genename $goodfilename`;
	if(length($line[0])>=2){
		#print "Got into head printing $genename\n";
		print $head;
		$head = "";
		$printFlag = 1;
        }
	$firstPrintFlag = 0;
	@line = ();
    }
    if($printFlag==1){
	#print "Got into body printing $genename\n";
	print $_;
    }
}




	        
		

