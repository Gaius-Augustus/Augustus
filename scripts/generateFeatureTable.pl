#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $exonCandFile;
my $OEfile;
my $sampledECfile;
my $outfile;
my $reffile;
my $printGFF = "";
my $help = 0;
my $introns = 0;

GetOptions(
    'out=s'=>\$outfile,
    'ec=s'=>\$exonCandFile,
    'oe=s'=>\$OEfile,
    'samp_ec=s'=>\$sampledECfile,
    'ref=s'=>\$reffile,
    'printGFF=s'=>\$printGFF,
    'classifyIntrons!'=>\$introns,
    'help!'=>\$help);

exec("perldoc $0") if ($help || !defined($reffile));

sub classify;
sub srtcmp;
sub cannot_overlap;
    


print STDERR "reading Ortho Exon list\n";

my %allExons = ();
my $numRedundant = 0; 
if(defined $OEfile){
    open(OE,"$OEfile") or die $!;
    
# OEs with same coordinates exist since we have overlapping gene ranges. Only keep the one with highest number of species involved
    while(<OE>){
	next if(m/#/);
	chomp;
	my $line = $_;
	my @oe = split(/\t/);
	my $key = "$oe[0]\t$oe[3]\t$oe[4]\t$oe[6]\t$oe[7]";
	my $numSpecies = 0;
	if ($oe[8] =~ m/;(n|numSpecies)=(\d+);/){
	    $numSpecies = $2;
	} else {
	    print STDERR "ERROR: could not detect number of species in ortho exon!\n";
	}
	if ($allExons{$key}){
	    $numRedundant++; # all but one exons with the same key are redundant
	}
	if (!$allExons{$key} || $allExons{$key}{"num"} < $numSpecies){
	    $allExons{$key} = { "attr" => $oe[8], "num" => $numSpecies, "score" => $oe[5] };
	}
    }
    close(OE);
}    

foreach my $key(keys %allExons){
    my $attr = $allExons{$key}{"attr"};
    my ($MLomega, $omega, $varOmega, $cons, $div, $contain, $length, $numSpecies, $type, $score, $oescore, $postProb, $mbp, $id, $gapProp, $baseFrameProp, $ec_pp, $fs1, $fs2, $fs3, $fs4, $fs5, $leftBfEx, $leftBfIn, $rightBfEx, $rightBfIn, $leftCons, $rightCons, $numSubs);
    $cons = $div = $contain = $length = $numSpecies = $type = $score = $oescore = $postProb = $mbp = $id = $gapProp = $baseFrameProp = $ec_pp = $fs1 = $fs2 = $fs3 = $fs4 = $fs5 = $leftCons = $rightCons = 0;

    $MLomega = $omega = $varOmega = $leftBfEx = $leftBfIn = $rightBfEx = $rightBfIn = $numSubs  = -1;
   
    if($attr =~ m/MLomega=([^;]+)/){
        $MLomega = $1;
    }
    if($attr =~ m/(E|PM)omega=([^;]+)/){
	$omega = $2;
    }
    if($attr =~ m/VarOmega=([^;]+)/i){
	$varOmega = $1;
    }
    if($attr =~ m/cons=([^;]+)/){
	$cons = $1;
    }
    if($attr =~ m/div=([^;]+)/){
	$div = $1;
    }
    if($attr =~ m/contain(ment)?=([^;]+)/){
	$contain = $2;
    }
    my @k = split(/\t/,$key);
    $length = $k[2] - $k[1] + 1;
    
    if($attr =~ m/(Note|type)=([^;]+)/){
	$type = $2;
    }
    if($attr =~ m/oescore=([^;]+)/){
	$oescore = $1;
    }
    if($attr =~ m/postProb=([^;]+)/){
        $postProb = $1;
    }
    if($attr =~ m/mbp=([^;]+)/){
        $mbp = $1;
    }
    if($attr =~ m/ID=([^;]+)/){
        $id = $1;
    }
    if($attr =~ m/gapProp=([^;]+)/){
        $gapProp = $1;
    }
    if($attr =~ m/baseFrameProp=([^;]+)/){
        $baseFrameProp = $1;
	chomp($baseFrameProp);
    }

    if($attr =~ m/ec_postProbs=([^;]+)/){
	$ec_pp = $1;
	chomp($ec_pp);
    }
    if($attr =~ m/numGaps=([^;]+)/){
        $fs1 = $1;
        chomp($fs1);
    }
    if($attr =~ m/numFrameShift=([^;]+)/){
        $fs2 = $1;
        chomp($fs2);
    }
    if($attr =~ m/numFrameShiftGap=([^;]+)/){
        $fs3 = $1;
        chomp($fs3);
    }
    if($attr =~ m/numWrongAli=([^;]+)/){
        $fs4 = $1;
        chomp($fs4);
    }
    if($attr =~ m/numWrongAliGap=([^;]+)/){
        $fs5 = $1;
        chomp($fs5);
    }
    if($attr =~ m/leftBoundaryExtOmega=([^;]+)/){
	$leftBfEx = $1;
	chomp($leftBfEx);
    }
    if($attr =~ m/rightBoundaryExtOmega=([^;]+)/){
        $rightBfEx = $1;
	chomp($rightBfEx);
    }
    if($attr =~ m/leftBoundaryIntOmega=([^;]+)/){
        $leftBfIn = $1;
        chomp($leftBfIn);
    }
    if($attr =~ m/rightBoundaryIntOmega=([^;]+)/){
        $rightBfIn = $1;
        chomp($rightBfIn);
    }

    if($attr =~ m/subst=([^;]+)/){
        $numSubs = $1;
	chomp($numSubs);
    }
    if($attr =~ m/LeftCons=([^;]+)/){
        $leftCons = $1;
        chomp($leftCons);
    }
    if($attr =~ m/rightCons=([^;]+)/){
        $rightCons = $1;
        chomp($rightCons);
    }
    

    $numSpecies = $allExons{$key}{"num"};
    $score = $allExons{$key}{"score"};
    $allExons{$key} = { "MLomega"=>$MLomega, "omega"=>$omega, "varOmega"=>$varOmega, "cons"=>$cons, "div"=>$div, "contain"=>$contain, "len"=>$length, "isOE"=>1, "type"=>$type, "num"=>$numSpecies, "score"=>$score, "oescore"=>$oescore, "postProb"=>$postProb, "mbProb"=>$mbp, "id"=>$id , "upSig"=>0, "downSig"=>0, "gapProp"=>$gapProp, "baseFrameProp"=>$baseFrameProp, "ec_pp"=>$ec_pp, "fs1"=>$fs1, "fs2"=>$fs2, "fs3"=>$fs3, "fs4"=>$fs4, "fs5"=>$fs5, "leftBfEx"=>$leftBfEx, "rightBfEx"=>$rightBfEx, "leftBfIn"=>$leftBfIn, "rightBfIn"=>$rightBfIn, "leftCons"=>$leftCons, "rightCons"=>$rightCons, "numSubs"=>$numSubs};
}

if(defined $sampledECfile){
    print STDERR "reading list of sampled exon candidates\n";
    open(SEC,"$sampledECfile") or die("Could not find file with sampled exon candidates: $sampledECfile");
    while(<SEC>){
	next if(m/#/ or (!$introns and !m/CDS/) or ($introns and !m/intron/));
	chomp;
	my $line = $_;
	my @sec = split(/\t/);
	$sec[7] = 0 if($introns);
	my $key = "$sec[0]\t$sec[3]\t$sec[4]\t$sec[6]\t$sec[7]";
	if($allExons{$key}){
	    if(!defined $allExons{$key}{"postProb"} or $allExons{$key}{"postProb"} == 0){
		if($sec[8] =~ m/postProb=([^;]+)/){
		    $allExons{$key}{"postProb"} = $1;
		}
		if($sec[8] =~ m/(avgBaseProb|mbp)=([^;]+)/){
		    $allExons{$key}{"mbProb"} = $2;
		}
		if($sec[8] =~ m/(Name|type)=([^;]+)/){
		    $allExons{$key}{"type"} = uc($2);
		}
		$allExons{$key}{"isOE"} = 1;
	    }
	}else{
	    $allExons{$key} = {};
	    if($sec[8] =~ m/postProb=([^;]+)/){
		$allExons{$key}{"postProb"} = $1;
	    }
	    if($sec[8] =~ m/(avgBaseProb|mbp)=([^;]+)/){
		$allExons{$key}{"mbProb"} = $2;
	    }
	    if($sec[8] =~ m/(Name|type)=([^;]+)/){
		$allExons{$key}{"type"} = uc($2);
	    }
	    $allExons{$key}{"len"} = $sec[4] - $sec[3] + 1;
	    $allExons{$key}{"isOE"} = 0;
	    $allExons{$key}{"score"} = $sec[5];
	}
    }
    close(SEC);
}


if(defined $exonCandFile){
    print STDERR "read list of exon candicates\n";
    open(EC,"$exonCandFile") or die $!;
     while(<EC>){
	next if(m/#/);
	chomp;
	my $line = $_;
	my @ec = split(/\t/);
	my $key = "$ec[0]\t$ec[3]\t$ec[4]\t$ec[6]\t$ec[7]";
	if(! $allExons{$key}){
	    $allExons{$key} = {};
	    $allExons{$key}{"len"} = $ec[4] - $ec[3] + 1;
	    $allExons{$key}{"isOE"} = 0;
	    $allExons{$key}{"score"} = $ec[5];
	    if($ec[8] =~ m/(Name|type)=([^;]+)/){
		$allExons{$key}{"type"} = $2;
	    }
	}
	$allExons{$key}{"upSig"} = $1 if ($ec[8] =~ m/upSig=([^;]+)/); 
	$allExons{$key}{"downSig"} = $1	if($ec[8] =~ m/downSig=([^;]+)/);
    }
    close(EC);
}


#print Dumper(\%allExons);

my @ex = ();
foreach my $key(keys %allExons){
    
    my @l = split(/\t/,$key);
    push(@ex, { "seq" => $l[0], "start" => $l[1], "end" => $l[2], "strand" => $l[3], "frame" => $l[4], "attr" => $allExons{$key} });
    delete $allExons{$key};
}


###### classify all exons #######

print STDERR "classify exon list\n";

classify(\@ex);

###### generate feature table ######

print STDERR "print feature table\n";

my $header = "class\tomega\tvarOmega\tMLomega\tleftBoundaryExtOmega\trightBoundaryExtOmega\tleftBoundaryIntOmega\trightBoundaryIntOmega\tconservation\tleftBoundaryCons\trightBoundaryCons\tdiversity\tcontainment\tlength\tnumSpecies\tisOrthoExon\tapostProb\tmeanBaseProb\ttype\tid\tupSig\tdownSig\tgapProp\tbaseFrameProp\tecPostProb\tnumGaps\tnumFrameShift\tnumFrameShiftGap\tnumWrongAli\tnumWrongAliGap\tnumSubst\tleftBoundFailCode\trightBoundFailCode\n";
if (defined $outfile){
    open(OUT,">$outfile") or die $!;
    print OUT $header;
} else {
    print $header;
}

if($printGFF ne ""){
    open(GFF,">$printGFF") or die $!;
}

foreach my $id(@ex){
    my $attr = $id->{"attr"};
#    print "$id->{seq}:$id->{start}..$id->{end}($id->{strand}):\t";

    my ($class, $omega, $varOmega, $MLomega, $cons, $div, $contain, $length, $numSpecies, $isOE, $postProb, $mbp, $type, $exon_id, $upSig, $downSig, $gapProp, $baseFrameProp, $ecPostProb, $fs1, $fs2, $fs3, $fs4, $fs5, $leftBfEx, $rightBfEx, $leftBfIn, $rightBfIn, $leftCons, $rightCons, $numSubs, $leftBoundFailCode, $rightBoundFailCode);
    $cons = $div = $contain = $length = $numSpecies = $isOE = $postProb = $mbp = $type = $exon_id = 
	$upSig = $downSig # upstream/downstream signal score (e.g. splice site scores
	= $gapProp = $baseFrameProp
	= $ecPostProb
	= $fs1 = $fs2 = $fs3 = $fs4 = $fs5 
	= $leftCons = $rightCons
	= 0;
    $omega = $varOmega = $MLomega = $leftBfEx = $rightBfEx = $leftBfIn = $rightBfIn = $numSubs = $leftBoundFailCode = $rightBoundFailCode = -1;
    $class=-5;

    $length = $attr->{len} if($attr->{len});
    $isOE = $attr->{isOE} if($attr->{isOE});
    $upSig = $attr->{upSig} if($attr->{upSig}); 
    $downSig = $attr->{downSig} if($attr->{downSig}); 
    if($isOE){
	$omega = $attr->{omega} if(defined $attr->{omega});
	$varOmega = $attr->{varOmega} if(defined $attr->{varOmega});
	$MLomega = $attr->{MLomega} if(defined $attr->{MLomega});
	$cons = $attr->{cons} if($attr->{cons});
	$div = $attr->{div} if($attr->{div});
	$contain = $attr->{contain} if($attr->{contain});
	$numSpecies = $attr->{num} if($attr->{num});
	$exon_id = $attr->{id} if($attr->{id});
	$gapProp = $attr->{gapProp} if($attr->{gapProp});
	$baseFrameProp = $attr->{baseFrameProp} if($attr->{baseFrameProp});
	$ecPostProb = $attr->{ec_pp} if($attr->{ec_pp});
	$fs1 = $attr->{fs1} if($attr->{fs1});
	$fs2 = $attr->{fs2} if($attr->{fs2});
	$fs3 = $attr->{fs3} if($attr->{fs3});
	$fs4 = $attr->{fs4} if($attr->{fs4});
	$fs5 = $attr->{fs5} if($attr->{fs5});
	$leftBfEx = $attr->{leftBfEx} if(defined $attr->{leftBfEx});
	$rightBfEx = $attr->{rightBfEx} if(defined $attr->{rightBfEx});
	$leftBfIn = $attr->{leftBfIn} if(defined $attr->{leftBfIn});
        $rightBfIn = $attr->{rightBfIn} if(defined $attr->{rightBfIn}); 
	$leftCons = $attr->{leftCons} if($attr->{leftCons});
        $rightCons = $attr->{rightCons} if($attr->{rightCons});
	$numSubs = $attr->{numSubs} if($attr->{numSubs} >= 0);
    }
    
    $postProb = $attr->{postProb} if($attr->{postProb});
    $mbp = $attr->{mbProb} if($attr->{mbProb});
    $type = $attr->{type} if($attr->{type});
    $class = $attr->{class} if(defined $attr->{class});
    $leftBoundFailCode = $attr->{leftBoundFailCode} if(defined $attr->{leftBoundFailCode});
    $rightBoundFailCode = $attr->{rightBoundFailCode} if(defined $attr->{rightBoundFailCode});

    my $line = "$class\t$omega\t$varOmega\t$MLomega\t$leftBfEx\t$rightBfEx\t$leftBfIn\t$rightBfIn\t$cons\t$leftCons\t$rightCons\t$div\t$contain\t$length\t$numSpecies\t$isOE\t$postProb\t$mbp\t$type\t$exon_id\t$upSig\t$downSig\t$gapProp\t$baseFrameProp\t$ecPostProb\t$fs1\t$fs2\t$fs3\t$fs4\t$fs5\t$numSubs\t$leftBoundFailCode\t$rightBoundFailCode\n";

    if (defined $outfile){
        print OUT $line;
    } else {
        print $line;
    }

    # print gff line
    if($printGFF ne ""){
        print GFF "$id->{seq}\tGFT\texon\t$id->{start}\t$id->{end}\t$attr->{score}\t$id->{strand}\t$id->{frame}\t";
	if($isOE){
	    print GFF "ID=$exon_id;PMomega=$omega;varOmega=$varOmega;MLomega=$MLomega;leftBoundaryExtOmega=$leftBfEx;rightBoundaryExtOmega=$rightBfEx;leftBoundaryIntOmega=$leftBfIn;rightBoundaryIntOmega=$rightBfIn;cons=$cons;div=$div;contain=$contain;numSpecies=$numSpecies;oescore=$attr->{oescore};gapProp=$gapProp;baseFrameProp=$baseFrameProp;ecPostProb=$ecPostProb;numGaps=$fs1;numFrameShift=$fs2;numFrameShiftGap=$fs3;numWrongAli=$fs4;numWrongAliGap=$fs5;subst=$numSubs;leftCons=$leftCons;rightCons=$rightCons;";
	}
	if($postProb>0){
	    print GFF "postProb=$postProb;mbp=$mbp;";
	}
	print GFF "type=$type;class=$class;leftBoundFailCode=$leftBoundFailCode;rightBoundFailCode=$rightBoundFailCode;\n";
    }

}


sub classify {

    my $protr = 0; # maximum allowed protrusion
    my %supportedRef = ();
    my %halfrightOE = ();
    my $numIdentical = 0;
    my @r = (); # reference coding exons
    my $oe_ref = shift; # orthologous exons

   
    open(REF,"$reffile") or die("Could not open ref file $reffile.");
    
    while(<REF>){
	if ((!$introns and m/\tCDS\t/) or ($introns and m/\tintron\t/)){
	    my @l = split /\t/;
	    $l[7] = 0 if($introns);
	    push(@r, {"seq" => $l[0], "start" => $l[3], "end" => $l[4], "strand" => $l[6], "frame" => $l[7], "line" => $_ });
	}
    }
    close REF;
#print Dumper(\@r);
    @r = sort {srtcmp($a,$b)} @r;
    my $numRef = @r;
    
    @{$oe_ref} = sort {srtcmp($a,$b)} @{$oe_ref};
    
    foreach my $oexon (@{$oe_ref}){
	my $real = 0;
	my $halfright = 0;
	my $leftBound = 0;
	my $rightBound = 0;
	my $false = 1;
	while (@r && cannot_overlap($oexon, $r[0])){
	    shift @r;
	}
	foreach my $ref (@r){
	    if($ref->{seq} eq $oexon->{seq} && !( $ref->{end} < $oexon->{start} || $oexon->{end} < $ref->{start}) ){
		$false = 0;
		if($oexon->{start} < $ref->{start}){
		    $leftBound = 2;
		}elsif($oexon->{start} == $ref->{start}){
		    $leftBound = 1;
		}else{
		    $leftBound = 3;
		}
		if($oexon->{end} > $ref->{end}){
		    $rightBound = 2;
		}elsif($oexon->{end} == $ref->{end}){
		    $rightBound = 1;
		}else{
		    $rightBound = 3;
		}
	    }
	    last if ($ref->{seq} ne $oexon->{seq} || $ref->{start} > $oexon->{end} + $protr);
	    # check for compatibility
	    if ($ref->{strand} eq $oexon->{strand} # same strand
		&& !($ref->{end} < $oexon->{start} || $oexon->{end} < $ref->{start}) # overlap of interval
		&& (($ref->{strand} eq "+"
		     && ($oexon->{frame} - $ref->{frame} + $oexon->{start} - $ref->{start})%3 == 0) # +, frame compatible
		    || ($ref->{strand} eq "-"
			&& ($ref->{end} - $oexon->{end} + $oexon->{frame} - $ref->{frame})%3 == 0)) # -, frame compatible
		){
		#print "included in exon: " . $ref->{line};
		$supportedRef{$ref->{line}}++;
		if ($oexon->{start} == $ref->{start} && $oexon->{end} == $ref->{end} && $oexon->{frame} == $ref->{frame}) {
		    $numIdentical++;
		    $real = 1;
		}else{
		    $halfrightOE{"$oexon->{seq}\t$oexon->{start}\t$oexon->{end}\t$oexon->{strand}\t$oexon->{frame}"} = 1;
		    $halfright = 1;
		}
	    }
	}
	if($real){
	    $oexon->{"attr"}{"class"} = 1;
	    $oexon->{"attr"}{"leftBoundFailCode"} = 1;
	    $oexon->{"attr"}{"rightBoundFailCode"} = 1;
	}elsif($halfright){
	    $oexon->{"attr"}{"class"} = -1;
	    if($leftBound != 0 and $rightBound != 0){
	    $oexon->{"attr"}{"leftBoundFailCode"} = $leftBound;
	    $oexon->{"attr"}{"rightBoundFailCode"} = $rightBound;
	    }else{
		print STDERR "no boundFailCode was assigned!\n";
	    }
	}elsif($false){
	    $oexon->{"attr"}{"class"} = 0;
	}else{
	    $oexon->{"attr"}{"class"} = -2; 
	    $oexon->{"attr"}{"leftBoundFailCode"} = $leftBound;
            $oexon->{"attr"}{"rightBoundFailCode"} = $rightBound;
	}
	if($introns and $oexon->{"attr"}{"class"} != 1){
	    $oexon->{"attr"}{"class"} = 0;
	}
    }
    my $sensitivity = $numIdentical/$numRef;
    my $specificity = $numIdentical/ @{$oe_ref};
    print STDERR "$numRedundant exons were ignored as an identical other exons exists with an alignment at least as deep (numSpecies)\n";
    print STDERR scalar(keys %supportedRef) . " of $numRef reference CDS are supported by exon candidates (frame-compatibly included).\n";
    print STDERR "$numIdentical of " . @{$oe_ref} . " exons are identical to a reference CDS (sensitivity of " . sprintf("%.2f", $sensitivity)*100 . "%)\n";
    print STDERR scalar(keys %halfrightOE) . " exons are frame-compatibly contained (up to $protr bp) in a reference CDS but not identical to a reference CDS.\n";
    print STDERR "specificity (with respect to identical exons): " . sprintf("%.4f", $specificity)*100 . "%\n";
}

# sort increasingly by sequence name, start, end
sub srtcmp {
    my $a = shift;
    my $b = shift;
    if ($a->{seq} lt $b->{seq}){
        return -1;
    } elsif ($a->{seq} gt $b->{seq}){
        return 1;
    }
    if ($a->{start} < $b->{start}){
        return -1;
    } elsif ($a->{start} > $b->{start}){
        return 1;
    }
    if ($a->{end} < $b->{end}){
        return -1;
    } elsif ($a->{end} > $b->{end}){
        return 1;
    }
    return 0;
}

sub cannot_overlap {
    my $a = shift;
    my $b = shift;

    # assume that no exon is longer than 10kb 
    if ($b->{seq} lt $a->{seq} or ( $a->{seq} eq $b->{seq} and $b->{start} < $a->{start} - 10000 ) ){
        return 1;
    }
    return 0;
}






__END__

=pod

=head1 NAME

generateFeatureTable.pl      convert exon gff files into feature table file (readable in R)

=head1 SYNOPSIS

generateFeatureTable.pl --oe=orthoExons.gff --ref=reffile.gtf
    
=head1 OPTIONS

    --ec=exonCandidates.gff   provide file with exon candidates (not neccessary ortho exons)                    
    --samp_ec=sampledEC.gff   provide file with sampled exon candidates
    --out=outfile.tbl         print table to outputfile.tbl
    --classifyIntrons         classify introncandidates instead of exoncandidates (DEFAULT false)
    --printGFF=oe.gff         prints gff lines with all feature information in the attribute column to file oe.gff                
    --help                    print this help message

=head1 DESCRIPTION

    example input:    
   
    input orthoExon file has to be in the following gff3 format:

CP000253        OE1     exon    555273  556745  0       +       0       ID=74;Name=74;Note=SINGLE;n=11;omega=0.177;Eomega=0.177;VarOmega=0.000788;cons=0.948;div=0.266;containment=0;labelpattern=111111-11111:111111-11111
CP000253        OE1     exon    555278  555367  0       -       0       ID=75;Name=75;Note=RSINGLE;n=11;omega=0.776;Eomega=0.776;VarOmega=0.224;cons=0.972;div=0.266;containment=0;labelpattern=000000-00000:000000-00000
CP000253        OE1     exon    555385  555546  0       +       0       ID=76;Name=76;Note=SINGLE;n=11;omega=1.42;Eomega=1.42;VarOmega=0.692;cons=0.965;div=0.266;containment=0;labelpattern=000000-00000:000000-00000
CP000253        OE1     exon    555415  555546  0       +       0       ID=77;Name=77;Note=SINGLE;n=11;omega=1.43;Eomega=1.43;VarOmega=0.701;cons=0.964;div=0.266;containment=30;labelpattern=000000-00000:000000-00000
CP000253        OE1     exon    555439  555546  0       +       0       ID=78;Name=78;Note=SINGLE;n=11;omega=1.35;Eomega=1.35;VarOmega=0.651;cons=0.961;div=0.266;containment=54;labelpattern=000000-00000:000000-00000

    input exon candidate file has to be in gff3 format

    CP000253        EC      exon    549445  549519  0       -       0       ID=1;Name=RSINGLE
    CP000253        EC      exon    549445  549534  0       -       0       ID=2;Name=RSINGLE
    CP000253        EC      exon    549445  549537  0       -       0       ID=3;Name=RSINGLE
    CP000253        EC      exon    549536  550738  0       +       0       ID=4;Name=SINGLE
    CP000253        EC      exon    549588  549686  0       +       0       ID=5;Name=SINGLE

    input sampled exon candidates in gff3 format

    CP000253        SAMPLED_ECs     exon    561684  562316  10.5    .       .       Name=single;score=10.5;postProb=1;avgBprob=0.9
    CP000253        SAMPLED_ECs     exon    562318  562866  4.27    .       .       Name=single;score=4.27;postProb=0.69;avgBprob=0.4
    CP000253        SAMPLED_ECs     exon    562328  562408  -10.3   .       .       Name=rsingle;score=-10.3;postProb=0.01;avgBprob=0.12
    CP000253        SAMPLED_ECs     exon    562338  562412  -10.1   .       .       Name=rsingle;score=-10.1;postProb=0.02;avgBprob=0.06

    input reference file in gtf format

    CP000253        protein_coding  start_codon     517     519     .       +       0       transcript_id "g1.t1"; gene_id "g1";
    CP000253        protein_coding  CDS     517     1878    .       +       0       transcript_id "g1.t1"; gene_id "g1";
    CP000253        protein_coding  stop_codon      1876    1878    .       +       0       transcript_id "g1.t1"; gene_id "g1";
    CP000253        protein_coding  start_codon     2156    2158    .       +       0       transcript_id "g2.t1"; gene_id "g2";

    example output: (features that do not occur in an exon are set to zero for that exon)

    class    omega   varOmega        conservation    diversity       containment     length  numSpecies      isOrthoExon     apostProb
    meanBaseProb    type
    0       1.51    0.74    0.406   0.131   0       81      5       1       0       0       SINGLE
    -1      0.898   0.264   0.416   0.116   0       549     4       1       0.01    0.978   RSINGLE
    0       0       0       0       0       0       66      0       0       0.01    0.01    rsingle
    -1      0.167   1.01e-06        0.629   0.297   948     543     7       1       0       0       RSINGLE
    1       0.177   0.000888        0.495   0.169   0       1359    6       1       1       1       RSINGLE
    -2      0       0       0       0       0       198     0       0       0       0       SINGLE
    -1      0.199   0.00392 0.8     0.24    84      135     10      1       0       0       RSINGLE


    values of "class":
     1  idential to a reference exon 
     0  not overlapping any reference exon
    -1  frame and strand compatible and overlapping but not identical
    -2  overlapping but not frame or strand compatible
    -5  if non of these classes are assigned (should not happen)

    values of "leftBoundFailCode" and "rightBoundFailCode":
     1     boundary correct
     2     boundary too long
     3     boundary too short

=cut
