#!/usr/bin/perl -w
################################################################################
#
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2010-07-07
#   Version: accompanying Scipio 1.4
#
#
#   yaml2gff.pl
#
#
#   This script transforms the yaml output of scipio into easily readable GFF 
#   
#   The output produces is conforming to the GFF3 standard; however the term "Target" used
#   in GFF3 is replaced by "Query".
#   In the context of Scipio, the "target" sequence is the referenced genomic sequence while the
#   "query" is the aligned protein sequence. In other contexts (as in the GFF standard too), the 
#   "target" sequence is what Scipio is referring to as the "query". To avoid confusion, in the
#   GFF produces by yaml2gff the term "target" is omitted completely. 
#
#   the protein sequence shown contains only matched amino acids, together with the following 
#   extra symbols:
#   "." for an unmatched amino acid in the query (gaps, case 5)
#   "x" for a mismatch / undetermined amino acid (cases 1-3, 10-13, 15, 16)
#   "-" for an extra codon (or pair or single nucleotide) in the target (cases 4, 6-9, 17)
#   "*" for a matched stop codon (case 14)
#
#   Counting of bps/aas: is transforming coordinates from YAML file (see usage.html)
#                        into low-high coordinates according to GFF standard
#



use strict;
use List::Util('sum', 'max');
use YAML;
use Getopt::Long;


my $PROT_WIDTH=98;
my $DNA_ORDER=0;

my $QPFX="prot";
my $TPFX="dna";
my $TRPFX="trans";
my $QNPFX="nucl";

my $filterStatus; # filter this status out
my $help=0;       # whether to print usage information
my %PARAMETER = ("filterstatus:s" => \$filterStatus,
		 "help!" => \$help);

sub usage
{   
    print STDERR "                                                                                                                                                                                                                           
    usage:                                                                                                                                                                                                                                   
    yaml2gff.pl [<options>] < scipio.yaml > scipio.gff                                                                                                                                                                                       
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             
    Options:    --help                     print this help message                                                                                                                                                                           
                --filterstatus=<value>     filter out alignments with given status, e.g. 'incomplete'                                                                                                                                        
";  
    exit 1;
}

sub nu_to_aa {
    my $arg = shift;
    my $residue = ($arg+1) % 3 -1;   # -1, 0, 1
    return ($arg-$residue)/3;
}

sub escape {
    foreach(@_) {
	s/([^a-zA-Z0-9.:^*$@!+_?-|])/sprintf("%%%02X",ord($1))/ge;
    }
}

sub with_number {
    my ($n, $s, $alt) = @_;

    $alt = "${s}s" unless (defined $alt);
    return ($n==1) ? "$n$s" : "$n$alt";
}

sub get_gapstr {
    my ($frameshifts, $qfrom, $qto) = @_;
    return "" unless(defined $frameshifts && @$frameshifts);

    my $result="";
    my @matches = map { &nu_to_aa($_->{"${QNPFX}_start"}) } @$frameshifts;
    $matches[0] -= $qfrom;
    push @matches, $qto;
    for my $i (1..@$frameshifts) {
	$matches[$i] -= &nu_to_aa($frameshifts->[$i-1]{"${QNPFX}_end"});
    }
    foreach (@$frameshifts) {
	my $base_count = ($_->{"${TPFX}_end"}) - ($_->{"${TPFX}_start"});
	my $aa_count = (($_->{"${QNPFX}_end"}) - ($_->{"${QNPFX}_start"})) / 3;
	$result.="M".(shift @matches)." ";
	if ($aa_count == 0) {
	    my $fs = $base_count % 3;
	    $base_count -= $fs;
	    my $deletes = $base_count / 3;
	    $result.="D$deletes " if ($deletes);
	    $result.="F$fs " if ($fs);
	} else {
	    $result.="I$aa_count ";
	    $result.="F$base_count " if ($base_count);
	}
    }	
    return ";Gap=${result}M$matches[0]";
}

sub output {
    my ($queryname, $hitlist) = @_;
    next if ($filterStatus && $hitlist->[0]{"status"} eq $filterStatus);
    escape($queryname);
    print "##query\t$queryname 1 ".$hitlist->[0]{"${QPFX}_len"};
    my ($pstart, $pend, $plen) = ($hitlist->[0]{"${QPFX}_start"}, @{$hitlist->[-1]}{"${QPFX}_end", "${QPFX}_len"});
    my $prev_end = 0;
    my $matched_prot = "." x $plen;
    my @minuses = ();
    foreach my $hitref (@$hitlist) {
	my ($hit_id, $matchings, $strand, $targetname, $score, 
	    $gapcount,  $stopcodon, $status, $protseq, 
	    $querystart, $queryend, $total_length) = 
	    @$hitref{"ID", "matchings", "strand", "target", "score", 
		     "gaps", "stopcodon", "status", "${QPFX}_seq",
		     "${QPFX}_start", "${QPFX}_end", "${QPFX}_len"};
	$targetname =~ s/\s.*//;
	escape($targetname);
	$gapcount = 0 unless (defined $gapcount);
	$stopcodon = "" unless (defined $stopcodon);
	
	if ($strand eq "-" || $strand =~ /minus/ || $strand =~ /back/) {
	    $strand = "-";
	} elsif ($strand eq "+" || $strand =~ /plus/ || $strand =~ /forward/) {
	    $strand = "+";
	} else {
	    my $dnapos = $hitref->{"${TPFX}_start"};
	    $dnapos = $hitref->{matchings}[0]{"${TPFX}_start"} unless defined $dnapos;
	    if (defined $dnapos) {
		$strand = $dnapos < 0 ? "-" : "+";
	    } else {
		$strand = ".";
	    }
	}
	my $complement = $strand eq "-";

	my @exons = grep { $_->{type} eq "exon" } @$matchings;
	my @mismatchcounts = map  { scalar @{$_->{mismatchlist}} } @exons;
	my @undeterminedcounts = map  { scalar @{$_->{undeterminedlist}} } @exons;
	my $total_mismatches = sum(@mismatchcounts,@undeterminedcounts,0);
	
	$querystart=0 unless (defined $querystart);
	$queryend=$total_length unless (defined $queryend);
	if (defined $protseq) {
	    substr($matched_prot, $querystart, $queryend-$querystart) = $protseq;
	}
	if ($prev_end < $querystart) {
	    print "\n# gap ".($prev_end+1)." $querystart";
	}

	my $matchsize = $queryend-$querystart;
	my $matchcount = $matchsize - $gapcount - $total_mismatches;

	my $has_bad_introns = 
	    grep ({ /intron\?/ } map { $_->{type} } @$matchings);

	my @gff_lines = ();
	foreach (@$matchings) {
	    my ($type, $qnfrom,$qnto) = @$_{"type", "${QNPFX}_start","${QNPFX}_end"};
	    if ($type eq "exon") {
		my ($qfrom, $qto) = map { &nu_to_aa($_); } ($qnfrom, $qnto);
		my ($frameshifts,$mismatches,$undetermined,$addprotseq, $trans,
		    $tfrom, $tto) 
		    = @$_{"seqshifts","mismatchlist","undeterminedlist","${QPFX}_seq","translation",
			  "${TPFX}_start", "${TPFX}_end"};
		if (defined $addprotseq) {
		    substr($matched_prot, $qfrom, $qto-$qfrom)=$addprotseq;
		}
		my ($first, $last) = $complement? (-$tto+1, -$tfrom+0) : ($tfrom+1, $tto+0);
		my @xpos = sort { $a <=> $b } (@$mismatches,@$undetermined);
		my $misstr = (@xpos)? ";Mismatches=".join(" ", @xpos) : "";
		map { substr($matched_prot,$_-1,1)="x"; } @xpos;
		my $gapstr = &get_gapstr($frameshifts, $qfrom, $qto);
		
		foreach (@$frameshifts) {
		    my ($fsstart, $fsend) = map { &nu_to_aa($_) } @{$_}{"${QNPFX}_start", "${QNPFX}_end"};
		    my $excount = int(($_->{"${TPFX}_end"} - $_->{"${TPFX}_start"} + 2) / 3);
		    if ($fsend > $fsstart) {
			substr($matched_prot, $fsstart, $fsend-$fsstart) = ("x" x $excount).("." x ($fsend-$fsstart-$excount));
		    } else {
			push @minuses, ($fsstart) x $excount;
		    }
		}

		push (@gff_lines, 
		      "$targetname\tScipio\tprotein_match\t$first\t$last\t$score\t$strand\t".
		      (-$qnfrom % 3)."\tID=$hit_id;Query=$queryname ".($qfrom+1).
		      " $qto$misstr$gapstr");
	    } elsif ($type eq "gap") {
		my ($qfrom, $qto) = map { &nu_to_aa($_) } ($qnfrom, $qnto);
		substr($matched_prot, $qfrom, $qto-$qfrom) = "." x ($qto-$qfrom);
		push @gff_lines, "# gap ".($qfrom+1)." $qto";
	    }
	}
	@gff_lines = reverse @gff_lines if ($complement && $DNA_ORDER);
	print "\n".join("\n", @gff_lines);
	$prev_end = $queryend;
    }   
    if (defined $pend && $pend < $plen) {
	print "\n# gap ".($pend+1)." $plen";
    }
    foreach (reverse @minuses) {
	substr($matched_prot, $_, 0) = "-";
    }
    $matched_prot = "protein sequence = [$matched_prot";
    do {
	print "\n# ".substr($matched_prot,0,$PROT_WIDTH,"");
    } while ($matched_prot);
    print "]\n#\n";
} 

### Command Line                                                                                                                                                                                                                             
&GetOptions (%PARAMETER) or &usage;
&usage if ($help);

my @input  = split("\n---",join("",<>));
shift @input if (@input && $input[0]!~/^---/);

print "##gff-version\t3\n";
print STDERR "Found ".&with_number(scalar @input, " YAML entry", " YAML entries").".\n";
foreach (@input) {
    s/^(---)?/---/;
    s/$/\n/;
    my $hits = YAML::Load($_);
    print STDERR "[";
    foreach (sort { $a cmp $b } keys %$hits) {
	&output($_, $hits->{$_});
	print STDERR " $_ ";
    }
    print STDERR "]\n";
}



