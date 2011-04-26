# Mario Stanke, August 17th, 2005

package SplicedAlignment;
use strict;
use locale;

#
# constants/thresholds
#

my $min_alilen      = 300;   # minimum total length of the aligned regions
my $min_cdslen      = 300;   # minimum total length of the inferred coding region
my $min_exonquality = 0.80;  # minimum percentage identity for each aligned region
my $min_avrgquality = 0.90;  # minimum percentage identity of complete aligned regions
my $min_intron_len  = 30;    # minimum length of all introns

#
# class SplicedAlignment
#

sub new {
    my $class = shift;
    my $self = {};
    bless ($self, $class);
    $self->{numexons}=0;
    $self->{estExonBegins}=();
    $self->{estExonEnds}=();
    $self->{exonBegins}=();
    $self->{exonEnds}=();
    $self->{qualities}=();
    $self->{intronFlag}=();
    $self->{estName} = shift;
    $self->{estLength} = shift;
    $self->{strand} = ".";
    $self->{genomicName} = shift;
    $self->{genomicLength} = shift;
    $self->{complement} = shift;
    $self->{status} = "OK";
    $self->{hasCDS} = 0;
    $self->{c5p} = 0;
    $self->{c3p} = 0;
    $self->{cdsBegins} = ();
    $self->{cdsEnds} = ();
    $self->{aliLen} = 0;
    $self->{codinglen} = 0;
    $self->{spliced5pUTR} = 0;
    $self->{spliced3pUTR} = 0;
    return $self;
}


#
# member functions
#


sub get_status {
    my $self = shift;
    return $self->{status};
}

sub get_contigname {
    my $self = shift;
    return $self->{genomicName};
}

sub get_complete5prime {
    my $self = shift;
    return $self->{c5p};
}

sub get_complete3prime {
    my $self = shift;
    return $self->{c3p};
}

sub addExon {
    my $self = shift;
    $self->{numexons}++;
    push @{$self->{estExonBegins}}, shift;
    push @{$self->{estExonEnds}}, shift;
    push @{$self->{exonBegins}}, shift;
    push @{$self->{exonEnds}}, shift;
    push @{$self->{qualities}}, shift; 
    push @{$self->{intronFlag}}, shift; 
}

sub checkAlignment {
    my $self = shift;
    
    # check whether there is an alignment at all
    if ($self->{numexons} == 0){
	$self->{status} = "no alignment found";
    }
    
    # check if there are no gaps in the aligned EST sequence
    # set strand if indicated by introns 
    foreach my $intronstatus (@{$self->{intronFlag}}){
	if ($intronstatus eq "=="){
	    $self->{status} = "not spliced alignment";
	    return;
	}
	if ($intronstatus eq "->"){
	    if ($self->{strand} eq '.') {
		$self->{strand} = '+';
	    } elsif ($self->{strand} eq '-') {
		$self->{status} = "inconsistent splice sites";
		return;
	    }
	}
	if ($intronstatus eq "<-"){
	    if ($self->{strand} eq '.') {
		$self->{strand} = '-';
	    } elsif ($self->{strand} eq '+') {
		$self->{status} = "inconsistent splice sites";
		return;
	    }
	}
    }
    # check intron lengths
    for (my $i=0; $i<$self->{numexons}-1; $i++){
	if($self->{exonBegins}[$i+1] - $self->{exonEnds}[$i] - 1 < $min_intron_len){
	   $self->{status} = "intron too short";
	   return;
	} 
    }

    # determine alignment length
    $self->{aliLen} = 0;
    for (my $i=0; $i<$self->{numexons}; $i++){
	$self->{aliLen} += $self->{exonEnds}[$i] - $self->{exonBegins}[$i] + 1; 
    }
    
    print STDERR "$self->{genomicName}:$self->{estName}\talignment length = $self->{aliLen}\n";
    if ($self->{aliLen} < $min_alilen) {
	$self->{status} = "alignment too short";
	return;
    }
    
    # determine average alignment quality
    # dicard spliced alignment if any exon quality is below 80% 
    # or the average quality is below 90%
    my $avqual=0;
    for (my $i=0; $i<$self->{numexons}; $i++){
	$avqual += ($self->{estExonEnds}[$i] - $self->{estExonBegins}[$i] + 1) * $self->{qualities}[$i];
	if ($self->{qualities}[$i] < $min_exonquality){
	    $self->{status} = "bad alignment quality";
	}
    }
    $avqual /= $self->{aliLen};
    print STDERR "average quality = $avqual\n";
    if ($avqual < $min_avrgquality) {
	$self->{status} = "bad alignment quality";
    }
    if ($self->{status} eq "bad alignment quality"){
	return;
    }
    
    $self->{status} = "OK";
}

sub makeGene {
    my $self = shift;
    my $seq = shift;
    if (length $seq != $self->{genomicLength}) {
	alert("Error in sequence $self->{genomicName}:$self->{estName}. Sequence length " . (length $seq) . " != $self->{genomicLength}");
    }
    #print STDERR "seq =$seq\n";    
    #construct the presumable mrna sequence
    my $mrna="";
    for (my $i=0; $i<$self->{numexons}; $i++) {
	$mrna .= substr $seq, $self->{exonBegins}[$i], $self->{exonEnds}[$i] - $self->{exonBegins}[$i]+1;
    } 
    my ($ORFbegin, $ORFend, $complete5prime, $start, $complete3prime) = findLongestORF($mrna, $self->{strand});
    print STDERR "findLongestORF: " , $ORFbegin, ", ", $ORFend, ", ", $complete5prime, ", ", $start, ", ", $complete3prime, "\n"; 

    $self->{strand} = ($ORFbegin < $ORFend) ? '+' : '-';
    $self->{c5p}=$complete5prime;
    $self->{c3p}=$complete3prime;
    print STDERR "ORF: $ORFbegin-$ORFend\n";
    # construct the CDS exons by mapping the mRNA back to the genomic region
    my ($mrnaseen, $cdsseen, $offset, $i);
    if ($self->{strand} eq '+') { # forward strand gene
	$mrnaseen = $self->{exonEnds}[0] - $self->{exonBegins}[0]+1;
	$cdsseen = 0;
	$offset = $complete5prime? $start : $ORFbegin;
	$self->{codinglen} = 1 + $ORFend - ($complete5prime? $start : $ORFbegin);
	$i = 0;
	while ($mrnaseen < $offset) {
	    $i++;
	    $mrnaseen += $self->{exonEnds}[$i] - $self->{exonBegins}[$i] + 1;
	}
	push @{$self->{cdsBegins}},  $self->{exonEnds}[$i] - ($mrnaseen - $offset) + 1;
	$cdsseen = $mrnaseen - $offset;
	while ($cdsseen < $self->{codinglen}){
	    push @{$self->{cdsEnds}},  $self->{exonEnds}[$i];
	    $i++;
	    push @{$self->{cdsBegins}}, $self->{exonBegins}[$i];
	    $cdsseen += $self->{exonEnds}[$i]-$self->{exonBegins}[$i]+1;
	}
	push @{$self->{cdsEnds}},  $self->{exonEnds}[$i]-($cdsseen-$self->{codinglen});
    } else { # reverse strand gene
	$mrnaseen = $self->{exonEnds}[$self->{numexons}-1] - $self->{exonBegins}[$self->{numexons}-1] + 1;
	$cdsseen = 0;
	$offset = (length $mrna) - ($complete5prime? $start : $ORFbegin) - 1;
	$self->{codinglen} = 1 + ($complete5prime? $start : $ORFbegin) - $ORFend;
	$i = $self->{numexons}-1;
	while ($mrnaseen < $offset) {
	    $i--;
	    $mrnaseen += $self->{exonEnds}[$i] - $self->{exonBegins}[$i] + 1;
	}
	unshift @{$self->{cdsEnds}}, $self->{exonBegins}[$i] + ($mrnaseen - $offset) - 1;
	$cdsseen = $mrnaseen - $offset;
	while ($cdsseen < $self->{codinglen}){
	    unshift @{$self->{cdsBegins}},  $self->{exonBegins}[$i];
	    $i--;
	    unshift @{$self->{cdsEnds}}, $self->{exonEnds}[$i];
	    $cdsseen += $self->{exonEnds}[$i]-$self->{exonBegins}[$i]+1;
	}
	unshift @{$self->{cdsBegins}},  $self->{exonBegins}[$i] + ($cdsseen-$self->{codinglen});
    }
    $self->{hasCDS} = 1;
    if ($self->{codinglen} < $min_cdslen){
	$self->{status} = "short CDS";
    }
}

sub findSplicedUTR {
    my $self = shift;
    if($self->{hasCDS}){
	# find spliced 5' UTR
	if ($self->{c5p}){
	    if ((($self->{strand} eq "+") && ($self->{exonEnds}[0] < $self->{cdsBegins}[0])) ||
		(($self->{strand} eq "-") && ($self->{exonBegins}[$#{$self->{exonBegins}}] > $self->{cdsEnds}[$#{$self->{cdsEnds}}]))){
		$self->{spliced5pUTR}=1;
	    }
	}
	# find spliced 3' UTR
	if ($self->{c3p}){
	    if ((($self->{strand} eq "+") && $self->{exonBegins}[$#{$self->{exonBegins}}] > $self->{cdsEnds}[$#{$self->{cdsEnds}}])||
		(($self->{strand} eq "-") && ($self->{exonEnds}[0] < $self->{cdsBegins}[0]))){
		$self->{spliced3pUTR}=1;
	    }
	}
    }
}

sub output {
    my $self = shift;
    my $outstring="";
    $outstring .= "# sequence $self->{genomicName}, len=$self->{genomicLength}\n";
    $outstring .= "# gene constructed from EST: $self->{estName}, len=$self->{estLength}, ";
    if ($self->{complement}){
	$outstring .= "(other strand)";
    } else {
	$outstring .= "(same strand)";
    }
    if ($self->{hasCDS}){
	if ($self->{c5p} && $self->{c3p}) {
	    $outstring .= ", CDS complete";
	} elsif ($self->{c5p} && !$self->{c3p}) {
	    $outstring .= ", CDS incomplete at 3'";
	} elsif (!$self->{c5p} && $self->{c3p}) {
	    $outstring .= ", CDS incomplete at 5'";
	} else {
	    $outstring .= ", CDS incomplete at both ends";
	}
    }
    if ($self->{spliced5pUTR}) {
	$outstring .= ", spliced 5'UTR";
    }
    if ($self->{spliced3pUTR}) {
	$outstring .= ", spliced 3'UTR";
    }

    $outstring .= "\n# alignment length = " . $self->{aliLen}. ", coding length = " . $self->{codinglen} ."\n";
    for (my $i=0; $i<$self->{numexons}; $i++) {
	$outstring .= "$self->{genomicName}\tmario\tmrna\t";
	$outstring .= (1+$self->{exonBegins}[$i]) . "\t";
	$outstring .= (1+$self->{exonEnds}[$i]) . "\t";
	$outstring .= "$self->{qualities}[$i]\t";
	$outstring .= "$self->{strand}\t";
	$outstring .= ".\t";
	$outstring .= "gene_id=\"$self->{genomicName}:$self->{estName}\";";
	$outstring .= "EST=" . ($self->{estExonBegins}[$i]+1). "-" . ($self->{estExonEnds}[$i]+1) . ";";
	if ($self->{intronFlag}[$i] ne ""){
	    $outstring .= "intronflag=$self->{intronFlag}[$i];";
	}
	$outstring .= "\n";
    }
    if ($self->{hasCDS}){
	for (my $i=0; $i<=$#{$self->{cdsBegins}}; $i++) {
	    $outstring .= "$self->{genomicName}\tmario\tCDS\t";
	    $outstring .= (1+$self->{cdsBegins}[$i]) . "\t";
	    $outstring .= (1+$self->{cdsEnds}[$i]) . "\t";
	    $outstring .= ".\t";
	    $outstring .= "$self->{strand}\t";
	    $outstring .= ".\t";
	    $outstring .= "gene_id=\"$self->{genomicName}:$self->{estName}\";";
	    $outstring .= "\n";
	}
    }
    return $outstring;
}

#
# non-member functions
#

#
# find longes open reading frame in sequence
#
# if begin < end it is on the forward strand
# if begin > end it is on the reverse strand

sub findLongestORF {
    my $mrna = shift;
    my $strand = shift;
    my $n = length $mrna;
    my ($begin, $end, $complete5prime, $start, $complete3prime) = findLongestORFOnPlusstrand($mrna);
    # print "plus  longest orf $begin, $end\n";
    my $rcmrna = reverseComplement($mrna);
    my ($b, $e, $c5p, $s, $c3p) = findLongestORFOnPlusstrand($rcmrna);
    # print "minus longest orf $b, $e\n";
    if (($strand eq '-' || ($strand eq '.' && $e-$b > $end-$begin))) {
	$begin = $n-1 - $b;
	$end = $n-1 - $e;
	$complete5prime=$c5p;
	$start=$n-1 - $s;
	$complete3prime=$c3p;
    }
    return ($begin, $end, $complete5prime, $start, $complete3prime);
}

sub findLongestORFOnPlusstrand {
    my $mrna = shift;
    my $begin=0;
    my $end=0;
    my $complete5prime=0; # boolean
    my $complete3prime=0; # boolean
    my $start =-1; # position of the start codon
    my ($curbegin, $curend);
    my $n = length $mrna;

    my @stppos = ((),(),());
    while($mrna =~ /taa|tga|tag/icg){
	my $m = $-[0];
	push @{$stppos[$m % 3]}, $m;
    }
    # print STDERR "mrna=$mrna\n";
    
    for (my $rf=0; $rf<3; $rf++){
	my $last = $n - ($n % 3) + $rf;
	if ($last>$n) {
	    $last -= 3;
	}
	push @{$stppos[$rf]}, $last;
	#print STDERR "rf=$rf, ", (join " ", @{$stppos[$rf]}), "\n";

	$curbegin = $rf;
	$curend = $curbegin;
	while($curend = shift @{$stppos[$rf]}){
	    next unless ($curend % 3 == $rf);
	    if ($curend-$curbegin > $end-$begin) {
		$end = $curend;
		$begin = $curbegin;
	    }
	    $curbegin=$curend;
	}
    }
    # check whether the ORF is complete at the 3 prime end
    if ($end < $n-2) {
	$complete3prime=1;
	$end += 3;
    }
    # check whether the ORF is likely to be complete at the 5 prime end
    # by checking whether there is another in-frame stop codon upstream of the start codon.
    if ($begin>2){
	$complete5prime=1;
	$begin += 3; # skip the stop codon
    }
    $start = $begin;
    while ($start<$end && (uc(substr $mrna, $start, 3) ne "ATG")){
	$start += 3;
    }
    if ($start >= $end){
	$complete5prime=0;
	$start = $begin;
	print STDERR "no ATG in the largest ORF\n";
    }

    return ($begin, $end-1, $complete5prime, $start, $complete3prime);
}

sub reverseComplement {
    my $seq = shift;
    my $rc = "";
    my %rcmap = ('A' => 'T', 'a' => 'T',
		 'C' => 'G', 'c' => 'G',
		 'G' => 'C', 'g' => 'C',
		 'T' => 'A', 't' => 'T');
    foreach my $char (split "", reverse $seq) {
	if ($rcmap{$char}) {
	    $rc .= $rcmap{$char};
	} else {
	    $rc .= 'n';
	}
    }
    return $rc;
}

1;

