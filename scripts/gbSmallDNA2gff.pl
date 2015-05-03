#!/usr/bin/perl
#
# create a genome-scale gff file from a gb file that created for AUGUSTUS training (has a flanking region, coordinates of the origin of sequence are given in the LOCUS tag)
#
# Katharina J. Hoff, Oct 15th 2012

#use strict;
#use Getopt::Long;



my $usage = "\nUsage: $0 file.gb > file.gff\n";

if(@ARGV != 1){
    print $usage;
    exit;
}

my $chr;
my $begChr;
my $endChr;
my $begSrc;
my $endSrc;
my $mRNAstr;
my $inmRNA = 0;
my $geneName;
my $mRNAbeg = 0;
my $mRNAend;
my $CDSstr;
my $inCDS = 0;
my $CDSbeg = 0;
my $CDSend;
my $haveSeenGene = 0;

open(GB, "<", $ARGV[0]) or die("Could not open genbank input file $ARGV[0]!\n");
while(<GB>){
    if(m/LOCUS/){
	m/\s(\S+)_(\d+)-(\d+)\s/;
	$chr = $1;
	$begChr = $2;
	$endChr = $3;
	$haveSeenGene = 0;
	$stopReading = 0;
    }elsif(m/source/){
	m/\s(\d+)\.\.(\d+)\s/;
	$begSrc = $1;
	$endSrc = $2;
	
    }elsif(m/mRNA/ && $haveSeenGene == 0){
	$haveSeenGene = 1;
	chomp;
	s/\s+mRNA\s+//;
	if(m/complement/){
	    $strand = "-";
	    s/\s+complement\(//;
	}else{
	    $strand = "+";
	}
	s/.*join\(//;
	$mRNAstr = $_;
	$inmRNA = 1;
    }elsif(m/mRNA/ && $haveSeenGene == 1){
	$stopReading = 1;
    }elsif(m/CDS/){
	chomp;
	s/\s+CDS\s+//;
        if(m/complement/){
            $strand = "-";
            s/\s+complement\(//;
        }else{
            $strand = "+";
        }
        s/.*join\(//;
        $CDSstr = $_;
        $inCDS = 1;
    }elsif($inmRNA == 1 && not(m/gene/) && $stopReading==0){
	chomp;
	s/\s+//;
	s/\)//g;
	$mRNAstr = $mRNAstr.$_;
    }elsif($inCDS == 1 && not(m/gene/)){
        chomp;
        s/\s+//;
        s/\)//g;
        $CDSstr = $CDSstr.$_;
    }elsif(m/gene/ && $inmRNA==1 && $stopReading == 0){
	$inmRNA = 0;
	m/gene=\"(.+)\"/;
	$geneName = $1;
	@t = split(/,/, $mRNAstr);
	foreach(@t){
	    m/(\d+)\.\.(\d+)/;
	    if($mRNAbeg == 0){
		$mRNAbeg = $begChr+$1-1;
	    }
	    $mRNAend = ($begChr+$2-1);
	    print "$chr\tsomeSource\texon\t".($begChr+$1-1)."\t".($begChr+$2-1)."\t.\t$strand\t.\tgene_id \"$geneName\"; transcript_id \"$geneName\";\n";
	}
	if($begChr < ($mRNAbeg-1)){
	    print "$chr\tsomeSource\tflanking_region\t".$begChr."\t".($mRNAbeg-1)."\t.\t$strand\t.\tgene_id \"$geneName\"; transcript_id \"$geneName\";\n";
	}
	print "$chr\tsomeSource\tmRNA\t$mRNAbeg\t$mRNAend\t.\t$strand\t.\tgene_id \"$geneName\"; transcript_id \"$geneName\";\n";
        print "$chr\tsomeSource\tgene\t$mRNAbeg\t$mRNAend\t.\t$strand\t.\tgene_id \"$geneName\";\n";
	if(($mRNAend+1) < $endChr){
	    print "$chr\tsomeSource\tflanking_region\t".($mRNAend+1)."\t".($endChr)."\t.\t$strand\t.\tgene_id \"$geneName\"; transcript_id \"$geneName\";\n";
	}
	$mRNAbeg = 0;     
    }elsif(m/gene/ && $inCDS==1){
        $inCDS = 0;
        m/gene=\"(.+)\"/;
        $geneName = $1;
        @t = split(/,/, $CDSstr);
        foreach(@t){
            m/(\d+)\.\.(\d+)/;
            $CDSbeg = $begChr+$1-1;
            $CDSend = ($begChr+$2-1);
	    print "$chr\tsomeSource\tCDS\t$CDSbeg\t$CDSend\t.\t$strand\t.\tgene_id \"$geneName\"; transcript_id \"$geneName\";\n";
        }
        $CDSbeg = 0;
    }
}

close(GB)
