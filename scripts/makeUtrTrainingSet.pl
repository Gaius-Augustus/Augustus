#!/usr/bin/perl
#
# Construct likely UTR regions from EST alignments and codon region boundaries
# e.g. to make a training set to train the AUGUSTUS UTR models.
#
# Mario Stanke, 1.4.2008

use strict;
use Getopt::Long;

my $usage = "$0 -- make a genbank or gff training file with 3' UTR and/or 5' UTR regions.\n\n";
$usage .= "EST alignments are used to guess the UTR and its end point.\n";
$usage .= "Usage: $0 codons.gff seq.fa ests.psl trainfile\n\n";
$usage .= "codons.gff is a file with (predicted) stop and/or start codons (containing stop_codon/start_codon lines\n";
$usage .= "seq.fa is the (multiple) fasta file with the assembly\n";
$usage .= "ests.psl is the (filtered) blat output of the ests against seq.fa\n";
$usage .= "train is the prefix-name of output files with the annotation for augustus training.\n";
$usage .= "two files will be created: train.gb (with genbank format) and train.gff (with gff format).\n";
$usage .= "options:\n";
$usage .= "--onlybest output for each stop/start codon only the 3'/5'-UTR from the most frequent splice variant\n";
$usage .= "--dist=n   consider all alignments that start at most n bp downstream of stop codon/upstream of the start codon (default 0)\n";

my $dist = 0;
my $minintronlen = 40;
my $estContaminMargin = 30; # This much of the EST is allowed to mismatch at the 3' end (3'UTR) or 5'end (5'UTR)
                            # Apart from this allowed margin, the EST must align up to its end.
my $radius = 20;            # Radius for smoothing/consensus building at the transcription end
my $thresh = 2*$radius+3;   # More than two exactly matching EST ends needed
my $gff = 0;
my $onlybest = 0;
my @f;
my @psl;
my %seqs =();
my %ali; # keys: genomic seqnames, value: (values: arrays: index is block, values: arrays of psl lines)
my $blocksize = 500000; # for indexing which bucket/block an alignment belongs to
my @hits;
my $kind; # stop or start meaning 3'UTR or 5'UTR
my ($codon, $strand, $name, $seq, $gname, $trname, $median);

my %c; # hash: Keys: genomic sequence names, values: rerefence to (sorted) array of start and stop positions


GetOptions('gff!'=>\$gff,
	   'onlybest!'=>\$onlybest,
           'dist=i'=>\$dist);

if ($#ARGV != 3) {
    die "Unknown option\n\n$usage";
}

my $codonfname = $ARGV[0];
my $seqfname = $ARGV[1];
my $pslfname = $ARGV[2];
my $trainfname = $ARGV[3];



open(CODONS, "<$codonfname") or die ("Could not open $codonfname");
open(PSL, "<$pslfname") or die ("Could not open $pslfname");
open(SEQ, "<$seqfname") or die ("Could not open $seqfname");
@psl = <PSL>;
close PSL;

# store the EST alignments in a hash indexed by the sequence name
foreach my $line (@psl) {
    my @f = split /\t/, $line;
    $ali{$f[13]} = [] unless exists $ali{$f[13]};
    my $blockstart = int($f[15]/$blocksize);
    my $blockend = int($f[16]/$blocksize);
    for (my $block = $blockstart; $block <= $blockend; $block++){ # usually just a single block
	$ali{$f[13]}->[$block] = [] if (!defined($ali{$f[13]}->[$block]));
	push @{$ali{$f[13]}->[$block]}, $line;
    }
}

# read in all sequences and store them in hash
$/="\n>";
while(<SEQ>){
    />?(.*)\n/;
    $name = $1;
    $seq = $'; #'
    $name =~ s/\s+.*//;
    $seq =~ s/>//;
    $seq =~ s/\n//g;
    if (exists $seqs{$name}) {
	warn "duplicate sequence name: $name\n";
    }
    $seqs{$name} = $seq;
}
close SEQ;
$/="\n";

open(GB, ">$trainfname.gb") or die ("Could not write $trainfname.gb.");
open(GFF, ">$trainfname.gff") or die ("Could not write $trainfname.gff."); 

PrepareCodonIdx();

my $nOfoverlap = 0;            # number of filtered hints by overlapsGene
my $nOfhints = 0;              # number of compatible hints if one wouldn't use the overlapsGene

while(<CODONS>) {

    my %splicevar;
    @f = split /\t/, $_, 9;
    if (@f < 8) { next }
    
    if ($f[2] eq 'stop_codon') {
	$kind = "stop";
    } elsif ($f[2] eq 'start_codon') {
	$kind = "start";
    } else {
	next;
    }
    $name = $f[0];
    $codon = $f[4];
    $f[8] =~ /gene_id."(.*)"/;
    $gname = $1;
    $f[8] =~ /transcript_id."([^"]*)"/; # "
    $trname = $1;
    $strand = $f[6];
    if (($strand eq '-' && $kind eq 'stop') || ($strand eq '+' && $kind eq 'start')) {
	$codon -=2;
    }
    #print "$kind codon $codon $strand\n";
    # search the est alignments for ones that constitute UTRs like this:
    # 
    # stop codon             >
    #        
    #                      -----------    ----       
    #                      -----------    -----
    #                      -----------    ---
    #test
    #print "codon: $codon\n" 
    #print "strand: $strand\n";
    #test
    if (defined $ali{$name}) { # sequence has EST alignments at all
	undef $seq;
	my $block = int($codon/$blocksize);
	if (defined($ali{$name}->[$block])){
	    @hits = @{$ali{$name}->[$block]};
	} else {
	    @hits = ();
	}
			
	foreach my $line (@hits) {
	    @f = split /\t/, $line;
	    my $aliStrand = $f[8];
	    my $queryLength = $f[10];
	    my $queryBegin = $f[11];
	    my $queryEnd = $f[12];
	    my $targetBegin = $f[15];
	    my $targetEnd = $f[16];
	    my $isEstEnd = 0;

	    if (((($kind eq 'stop' && $strand eq $aliStrand) || ($kind eq 'start' && $strand ne $aliStrand)) && $queryEnd >= $queryLength - $estContaminMargin) ||
		((($kind eq 'stop' && $strand ne $aliStrand) || ($kind eq 'start' && $strand eq $aliStrand)) && $queryBegin <= $estContaminMargin)){
		$isEstEnd = 1;
	    }
	    if ($isEstEnd && # the EST/cDNA aligns almost up to its end
		((($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-')) && $targetBegin <= $codon + $dist && $targetEnd >= $codon) ||
		((($kind eq 'stop' && $strand eq '-') || ($kind eq 'start' && $strand eq '+')) && $targetBegin <= $codon && $targetEnd >= $codon - $dist)) { # alignment overlaps start/stop codon or is not too far away from it
		if (!defined $seq) {
		    die ("Sequence $name not in sequence file.") unless (defined $seqs{$name});
		    $seq = $seqs{$name};
		}
		# close small gaps and create list of segments @segs
		my @segs;
		my @blocklengths = split /,/, $f[18];
		chomp $f[20];
		my @targetbegins = split /,/, $f[20];
		

		# fill the gap to the stop/start codon in case the alignment is (at most $dist) downstream of stop codon
		if ((($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-'))  && $targetbegins[0] > $codon) {
		    @blocklengths = ($targetbegins[0]-$codon+1, @blocklengths);
		    @targetbegins = ($codon - 1, @targetbegins);
		} elsif ((($kind eq 'stop' && $strand eq '-') || ($kind eq 'start' && $strand eq '+')) && $targetbegins[-1] < $codon) {
		    push @blocklengths , ($codon - $targetbegins[-1]+1);
		    push @targetbegins , $targetbegins[-1]+1;
		}


		push @segs, $targetbegins[0]+1;
		for (my $i=0; $i < @blocklengths-1; $i++) {
		    if ($targetbegins[$i] + $blocklengths[$i] + $minintronlen < $targetbegins[$i+1]){
			push @segs, $targetbegins[$i] + $blocklengths[$i];
			push @segs, $targetbegins[$i+1]+1;
		    }
		}
		push @segs, $targetbegins[@targetbegins-1] + $blocklengths[@blocklengths-1];
		#print "segs= ", (join " ", @segs);
		#print "\n";
		# check whether the codon position is contained in a block at all
		my $ok = 0;
		for (my $i=0; $i < @segs && !$ok; $i += 2) {
		    if ($segs[$i] <= $codon && $segs[$i+1] >= $codon) {
			$ok = 1;
		    }
		}
		#print "alignment covers $kind codon: $ok\n";
		# check whether each remaining gap could be an intron
		if (@segs >1 && $ok) {
		    for (my $i=1; $i < @segs-1 && $ok; $i += 2) {
			my $intronbegin = $segs[$i]+1;
			my $intronend = $segs[$i+1]-1;
			#print "ibegin=$intronbegin, $splicevariantiend=$intronend\t";
			#print substr($seq, $intronbegin-1, 2) . "  " . substr($seq, $intronend-2, 2) . "\n"; 
			if (($strand eq '+' && (uc((substr($seq, $intronbegin-1, 2)) ne 'GT' && uc(substr($seq, $intronbegin-1, 2)) ne 'GC') ||
					       (uc(substr($seq, $intronend-2, 2)) ne 'AG'))) ||
			    ($strand eq '-' && (uc((substr($seq, $intronend-2, 2)) ne 'AC' && uc(substr($seq, $intronend-2, 2)) ne 'GC') ||
						(uc(substr($seq, $intronbegin-1, 2)) ne 'CT')))) {
			    $ok = 0;
			    #print "intron not ok\n";
			}
		    }
		}
		
		my ($a,$b); # fill in a and b, Z.B. a=alignment start, b= Start codon-3
		if (($strand eq "+" && $kind eq "start") || ($strand eq "-" && $kind eq "stop")) {
			$a = $targetBegin; #anfang des alignments
			$b = $codon - 3;
		} else {
			$a = $codon +3;
			$b = $targetEnd; # ende des alignments
		}
		
		
		my $geneOverlap=overlapsGene($name, $a, $b);
				
		if ($ok && !$geneOverlap) {	#falls alle Gaps Intron sind...
		    
		    $nOfhints++;		   
		    
		    if (($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-')) {
			while (@segs >= 4 && $segs[1] < $codon) {
			    shift @segs;	#entfernt das erste Element von @segs
			    shift @segs;
			}
		    } else {
			while (@segs >= 4 && $segs[-2] > $codon) {
			    pop @segs;	#entfernt das letzte Element von @segs
			    pop @segs;
			}
		    }
		    # build equivalence classes of alignments, the very start and very end do not matter
		    # make a hash with the string consisting of all coordinates but the first and the last as key
		    my @copy = @segs;
		    shift @copy;
		    pop @copy;
		    my $key = join " ", @copy; #join the array elements by " "
		    #print "key=$key\n";
		    $splicevar{$key} = [] unless exists $splicevar{$key};	
		    push @{$splicevar{$key}}, \@segs;
		} 
		else {
		    
			if ($ok && $geneOverlap) {
				$nOfoverlap++;
				$nOfhints++;
			    }
		    }
	    }
	}
			
	# now go through the equivalence classes of alignments and determine for each one the tts/tss as a cluster point
	my @ends;
	my @segs;
	my $numpeak;
	my $num;
	my $numalis;
	my $maxalis=0;
	foreach my $splicevariant (keys %splicevar) {
	    if (scalar(@{$splicevar{$splicevariant}}) > $maxalis) {
		$maxalis =  scalar(@{$splicevar{$splicevariant}});
	    }
	}
	my $haveOne=0;
	foreach my $splicevariant (keys %splicevar) {
	    $numalis = scalar(@{$splicevar{$splicevariant}});
	    next if ($onlybest && ($numalis < $maxalis || $haveOne));
	    #print "splice variant (" . $numalis . " alignments): $splicevariant\n";
	    @ends = ();
	    foreach my $segref (@{$splicevar{$splicevariant}}){ 
		@segs = @$segref;
		if (($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-')) {
		    push @ends, $segs[@segs-1];
		} else {
		    push @ends, $segs[0];
		}
	    }
	    sort {$a <=> $b } @ends;

	    # compute the likely transcription termination site (tts)  or transcripton start site (tss)
	    # take that end point that has the most other end points within a range of [-$radius,$radius] (weighting: 1 2 3 4 5 ... $radius+1 ... 5 4 3 2 1)
	    $numpeak = 0;
	    $median = -1;
	    foreach my $p (@ends){
		$num = 0;
		foreach my $q (@ends){
		    if(abs($q-$p) <= $radius) {
			$num += $radius+1-abs($q-$p);
		    }
		}
		if ($num > $numpeak) {
		    $numpeak = $num;
		    $median = $p;
		}
	    }
	    #print "ends: " . join(" ", @ends) . " median=$median\n";
	    if ($numpeak >= $thresh){	
		printGFF($codon, $splicevariant, $median, $strand, $numpeak);
		printGB($codon, $splicevariant, $median, $strand);
	    }
	    $haveOne=1;
	}
	
    }
}

print "$nOfoverlap hints were filtered because of gene overlap.\n";

print "$nOfhints hints would be compatible if the hints with gene-overlap wouldn't be filtered.\n"; 

close GB;

close GFF;

sub printGFF {
    my $codon = shift;
    my $splicevariant = shift;
    my $tend = shift;
    my $strand = shift;
    my $score = shift;
    my @exons;
    if (($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-')) {
	@exons = ($codon+1, split (/ /, $splicevariant), $tend);
    } else {
	@exons = ($tend, split (/ /, $splicevariant), $codon-1);
    }
    for (my $i=0; $i<@exons; $i+=2) {
	print GFF "$name\tmakeUtr\t" . ($kind eq 'stop'? "3'-UTR" : "5'-UTR") . "\t$exons[$i]\t$exons[$i+1]\t$score\t$strand\t.\ttranscript_id \"$trname\"; gene_id \"$gname\";\n";
    }
}

sub printGB {
    my $codon = shift;
    my $splicevariant = shift;
    my $tend = shift;
    my $strand = shift;
    my $length;
    my $locusname;
    my $beginpos;
    my $sequence;
    my ($a,$c,$g,$t,$n,$rest,$pos,$i,$j,$ten,$zahlzeile);
    my @join = split / /, $splicevariant;
    if (($kind eq 'stop' && $strand eq '+') || ($kind eq 'start' && $strand eq '-')) {
	die ("negative length") unless $tend >= $codon;
	return unless ($codon >= 100 && $tend + 400 < length $seq);
	$beginpos = $codon-99;
	$sequence = substr($seq, $beginpos, $tend - $codon + 100 + 400);
	@join =  ($codon, @join, $tend);
    } else {
	die ("negative length") unless $codon >= $tend;
	return unless ($tend >= 401 && $codon + 100 < length $seq);
	$beginpos = $tend-1-400;
	$sequence = substr($seq, $beginpos, $codon - $tend + 100 + 400);
	@join = ($tend, @join, $codon);
    }
    foreach my $p (@join){
	$p -= $beginpos;
    }
    $length = length $sequence;
    $locusname = $name . "_" . $gname;
    if ($splicevariant ne "") {
	my $variant = $splicevariant;
	$variant =~ s/ /_/g;
	$locusname .= ":" . $variant;
    }
    
    
    # print "codon=$codon\n";
    # print "splicevariant=$splicevariant\n";
    # print "tend=$tend\n";
    # print "strand=$strand\n";
   

    print GB  "LOCUS       $locusname   $length bp  DNA\n";
    print GB  "FEATURES             Location/Qualifiers\n";
    print GB  "     source          1..$length\n";    
    print GB  "     mRNA            ";
    my $joinstr;
    # create mRNA join string
    if ($strand eq '-') {
	$joinstr .= "complement(";
    }
    if (@join > 2) {
	$joinstr .= "join(";
    }
    for ($i=0; $i<@join-1; $i+=2) {
	if ($strand eq '+' && $i==0) {
	    $joinstr .= "" . $join[$i] . ".." . $join[$i+1];
	} elsif ($strand eq '-' && $i == @join-2){
	    $joinstr .= $join[$i] . ".." . $join[$i+1];
	} else {
	    $joinstr .= $join[$i] . ".." . $join[$i+1];
	}
	if ($i<@join-2) {
	    $joinstr .= ",";
	}
    }
    if (@join > 2) {
	$joinstr .= ")";
    } 
    if ($strand eq '-') {
	$joinstr .= ")";
    }
    
    while (length $joinstr > 0) {
	print GB  substr $joinstr,0,59;
	$joinstr = substr $joinstr, 59;
	if (length $joinstr > 0) {
	    print GB  "\n                     ";
	}
    }
    print GB  "\n";
    # create CDS join string
    if ($strand eq '+') {
	$joinstr = "" . ($codon-$beginpos) . ".." . ($codon-$beginpos);
    } else {
	$joinstr = "complement(" . ($codon-$beginpos) . ".." . ($codon-$beginpos) . ")";
    }
    print GB "     CDS             " . $joinstr . "\n";
    # print sequence
    $sequence =~ s/A/a/g;
    $sequence =~ s/C/c/g;
    $sequence =~ s/G/g/g;
    $sequence =~ s/T/t/g;
    $sequence =~ s/N/n/g;

    $a = $sequence =~ s/a/a/g;
    $c = $sequence =~ s/c/c/g;
    $g = $sequence =~ s/g/g/g;
    $t = $sequence =~ s/t/t/g;
    $n = $sequence =~ s/n/n/g;
    $rest = $length - $a - $c - $g - $t -$n;
    print GB  "BASE COUNT     $a a   $c c  $g g   $t t";
    if ($n>0) {
        print GB  "   $n n";
    }
    if ($rest > 0) {
          print GB  "   $rest ?";
    }
    print GB  "\nORIGIN\n";
    $i = 1;
    $pos = 0;
    while ($pos <= length $sequence) {
        $zahlzeile = "";
        for ($j=0; $j < 9-length "$i"; $j=$j+1) {
            print GB  " ";
        }
        print GB  "$i";
        for ($j=0; $j < 6; $j=$j+1) {
             $ten = substr $sequence, $pos, 10;
             if (length $ten > 0) {
                 print GB  " $ten";
             }
             $pos = $pos + 10;
        }
        print GB  "\n";
        $i += 60;
    }
    print GB  "//\n";
}


# prepare an index of where codons lie in the genome for fast access later

sub PrepareCodonIdx{
	open(CODONS_PRE, "<$codonfname") or die ("Could not open $codonfname"); 
	while(<CODONS_PRE>) {
   		@f = split /\t/, $_, 9;	
    		if (@f < 8) { next }
    		if ($f[2] eq 'stop_codon') {$kind = "stop";} elsif ($f[2] eq 'start_codon') {$kind = "start";} else {next;}
    		$name = $f[0];
    		$codon = $f[4];
    		$strand = $f[6];
    		if (($strand eq '-' && $kind eq 'stop') || ($strand eq '+' && $kind eq 'start')) {
			$codon -=2;	#nv
    		}
		# add $codon to array 
		$c{$name} = [] unless (exists($c{$name}));
		push @{$c{$name}}, $codon;
	};
	# now sort all lists of codons 
	foreach my $seq  (keys %c){
		@{$c{$seq}}=sort{$a <=> $b}@{$c{$seq}}	#sort numerically
	}	
	close CODONS_PRE;
}


sub overlapsGene{
	# binary sort in $c{seq}
	my $seq = shift;
	my $a = shift;
	my $b = shift;
	my ($l, $u) = (0, $#{$c{$seq}} );	## lower, upper end of search interval
	my $i;	
	while ($l<=$u){
		$i = int(($l + $u)/2);
		if ($c{$seq}->[$i]<$a){$l=$i+1;}
		elsif ($c{$seq}->[$i]>$b){$u=$i-1;}
		else {return $c{$seq}->[$i]} #the position number for overlap
	}
	return 0;	#not found
}
