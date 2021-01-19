#!/usr/bin/env perl
#
# This script is used by the braker.pl pipeline.
# Please be extremely careful when changing this script because the braker.pl
# pipeline may fail upon custom modification of this script.
# In case of doubt, contact katharina.hoff@uni-greifswald.de
#
# convert a gff/gtf file with CDS and mrna/exon annotations together with a 
# fasta file containing the sequences to a short genbank
# file that can be read by augustus and etraining containing
# the DNA only flanking a gene and not the whole sequence
#
# Mario Stanke, 18.09.2006
# Last modified by Katharina J. Hoff on Feb 21st 2018

use strict;
use warnings;
use Getopt::Long;

my $usage .= "$0 -- convert GFF file and sequence fasta file to minimal genbank format\n";
$usage .= "\n";
$usage .= "Usage: $0 gff-file seq-file max-size-of-gene-flanking-DNA output-file [options]\n";
$usage .= "options:\n";
$usage .= "--bad=badfile    Specify a file with gene names. All except these are included in the output.\n";
$usage .= "--good=goodfile  Specify a file with gene names. Only these genes are considered\n";
$usage .= "                 from the input, also for overlap detection.\n";
$usage .= "--overlap        Overlap filtering turned off.\n";
$usage .= "--connected      Do not cut a sequence into gene-pieces anymore.\n";
$usage .= "--softmasked     Keep softmasking information from input sequence\n";
$usage .= "--hardmask       Convert softmasking information from input sequence to hardmasking in output file\n";
$usage .= "\n";

if ($#ARGV < 3) {
    die "Unknown option $ARGV\n\n$usage";
}

my $badfilename = "";
my $exceptionfilename = "";
my $exceptiontype = "";
my ($overlap, $good, $bad, $connected, $softmasked, $hardmask);
my %exceptionlist=();
my $num_ambig_utr_mgs = 0;

GetOptions( 'good=s' => \$good, 'bad=s' => \$bad, 'overlap!' => \$overlap, 'connected!' => \$connected, 'softmasked!' => \$softmasked, 'hardmask' => \$hardmask);

my $gfffilename = $ARGV[0];
my $seqfilename = $ARGV[1];
my $flank = $ARGV[2];
my $outputfilename = $ARGV[3];
my $seqnameerrcount = 0;

$overlap = 0 if (!defined($overlap));
$connected = 0 if (!defined($connected));

if (defined($good) && defined($bad)){
  die "Good and bad cannot both be specified.\n";
}
if (defined($bad)){
    $exceptionfilename = $bad;
    $exceptiontype = "bad";
    print "Not using genes in $exceptionfilename.\n";
}
if (defined($good)){
    $exceptionfilename = $good;
    $exceptiontype = "good";
    print "Using only genes in $exceptionfilename.\n";
}

my $minContigLen = 3*$flank;

open(GFFFILE, "<$gfffilename") || die "Couldn't open $gfffilename.\n";

#
# If exception file exists, read in the names of the genes that should be skipped.
#
if ($exceptionfilename ne "") {
    open(EXCEPTIONFILE, "<$exceptionfilename") || die "Couldn't open $exceptionfilename.\n";
    foreach my $name (<EXCEPTIONFILE>){
	chomp $name;
	$exceptionlist{$name}=1;
    }
}


#
# data structure:
# annos: hash of annotations
# annotation: hash of gbfkeys, keys: CDS:genename or mRNA:genename
# gbfkey: list:
#         {name, strand, exon1start, exon1end, ..., exonXstart, exonXend}
#         name is "mRNA" or "CDS"
#
my %annos; # keys: sequences, elements: annotations
my $annotation; # reference to a hash
my ($fkey, $seqname, $type, $begin, $end, $grp, $strand, $genename);

my $cdsNumber = 0;
my $mrnaNumber = 0;
my $utrNumber = 0;

#
# Read in and store all annotations
#

while (<GFFFILE>) {
    s/#.*//;
    next unless /\S/;
    
    my @f = split /\t/, $_, 9;
    if (@f < 8) { warn @f,"Not GFF format"; next }
    $seqname = $f[0];
    $type = $f[2];
    $begin = $f[3];
    $strand = $f[6];
    $end = $f[4];
    $grp = $f[8];
    $grp =~ s/\n//;
    if ($grp =~ m/(transcript_id|Transcript)."([^"]+)"/){ #"
        $genename = $2;
    } elsif ($grp =~ m/gene_id."([^"]+)"/){ #"
        $genename = $1;
    } elsif ($grp =~ m/Parent=Transcript:([^;]+)/){
        $genename = $1;
    } elsif ($grp =~ m/Parent=([^;]+)/){
        $genename = $1;
    } else {
        $genename = $grp;
    }
    if ($type =~ /mrna/i || $type =~ /^exon$/i){
        $fkey = "mRNA";
        $mrnaNumber++;
    } elsif ($type =~ /cds/i || $type =~ /coding_exon/i){
	$fkey = "CDS";
        $cdsNumber++;
    } elsif ($type =~ /5'-UTR/i || $type =~ /3'-UTR/i || $type =~ /UTR/i || $type =~ /.*_UTR/i){
        $fkey = "UTR";
        $utrNumber++;
    } else {
        $fkey = "other";
    }
    foreach my $genename (split(",", $genename)){ # if an exon belongs to several transcripts ...
	if (($fkey eq "CDS" || $fkey eq "mRNA" || $fkey eq "UTR") && 
	    ($exceptiontype eq "" ||
	     ($exceptiontype eq "good" && (exists $exceptionlist{$genename})) ||
	     ($exceptiontype eq "bad" && !exists $exceptionlist{$genename}))) {
	    if (!exists($annos{$seqname})){
		my %newanno=();
		$annos{$seqname}=\%newanno;
	    }
	    $annotation = $annos{$seqname};
	    if (exists($annotation->{$fkey.$genename})){
		# just add the exon
		insertExon($annotation->{$fkey.$genename}, $begin, $end);
	    } else {
		my @newgbfkey;
		push @newgbfkey, ($fkey, $genename, $strand, $begin, $end);
		$annotation->{$fkey.$genename}=\@newgbfkey;
	    }
	}
    }
}

#
# insert exon at right place into sorted array
#
sub insertExon {
    my ($annoref, $a, $b) = @_;
    if ($a > $annoref->[-1]) {# sorted, just add at end of list
	push @{$annoref}, ($a,$b);
    } else {# not sorted
	my $i = 3;
	while ($i<@{$annoref} && $a > $annoref->[$i]){
	    $i += 2;
	}
	splice @{$annoref}, $i, 0, ($a,$b);
    }
}

# escape characters according GFF3 Format restrictions for Column 1: "seqid"  
# see http://gmod.org/wiki/GFF3
sub escape_seqid {
    my ($seqID) = @_;
    $seqID =~ s/([^a-zA-Z0-9\.:\^\*\$@!\+_\?\-\|])/sprintf("%%%02X",ord($1))/ge;
    return $seqID;
}

#
# Now write the data
#
open (OUTPUT, ">$outputfilename") || die "Could not open output file $outputfilename\n";
open(FASTA, "<$seqfilename") || die "Couldn't open $seqfilename\n";
my ($seq, $length, $an, $cn, $gn, $tn, $nn, $rest, $i, $pos, $zahlzeile, $join, $j, $ten);

$/="\n>";
while(<FASTA>) {
    /[>]*(.*)\n/;
    $seqname = $1;
    $seq = $';#'
	$seq =~ s/>//; 
    $seq =~ s/\n//g;
    
    $length = length $seq;
    $annotation = $annos{$seqname};
    my $annotatedseqname = $seqname;
    
    if (!defined($annotation)) {
        my $shortseqname = $seqname;
        $shortseqname =~ s/\s.*//;
        $annotatedseqname = $shortseqname;
        $annotation = $annos{$annotatedseqname}; # try the short sequence name (split at first whitespace)
        my $usedshortseqname = defined($annotation);
        
        if (!defined($annotation)) { # check if seq id in gff-file contains escaped characters
            $annotatedseqname = escape_seqid($seqname);
            $annotation = $annos{$annotatedseqname}; # try the escaped sequence name (as demanded by GFF3 format)
            
            if (!defined($annotation)) {
                $annotatedseqname = escape_seqid($shortseqname);
                $annotation = $annos{$annotatedseqname}; # try the escaped short sequence name
                $usedshortseqname = defined($annotation);
            }
        }
        
        if ($usedshortseqname) {
            $seqnameerrcount++;
            if ($seqnameerrcount <= 10) {
                print STDERR "Sequence $seqname has no annotation but $shortseqname has. ";
                print STDERR "Assuming that space truncates name.\n";
                if ($seqnameerrcount == 10) {
                    print STDERR "Supressing this error message from now on.\n";
                }
            }
            $seqname = $shortseqname;
        }
    }
    # For each UTR annotation check whether there is a CDS annoation with the same name and matching boundaries.
    # If yes, make an mRNA from the UTR and discard the UTR annotation.
    # This is necessary for gff files which have a UTR annotated separately instead of 'exon' for all biological exons.
    foreach my $gbfkeyref (values %$annotation) {
	if ($gbfkeyref->[0] eq "UTR") {
	    my $fehler = 0;
            my $genename = $gbfkeyref->[1];
            if (exists($annotation->{"CDS$genename"})) {
                 my $cdsref = $annotation->{"CDS$genename"};
                 if (exists($annotation->{"mRNA$genename"})) {
		     $num_ambig_utr_mgs++;
		     if ($num_ambig_utr_mgs < 11){
			 print "Have mRNA and UTR for gene $genename in sequence $seqname. ";
			 print "Ignoring UTR annotation and using mRNA annotation only.\n";
		     }
		     if ($num_ambig_utr_mgs == 10) {
			 print "Suppressing this error message from now on.\n";
		     }
                 } else {
		     # print "UTR " , (join " ", @$gbfkeyref), "\n";
		     # print "CDS " , (join " ", @$cdsref), "\n";
                       my @newgbfkey = ();
	               push @newgbfkey, ("mRNA", $genename, $gbfkeyref->[2]);
                       my $utrindex = 4;		    
                       while ($utrindex < @$gbfkeyref && $fehler == 0) {
			   if ($gbfkeyref->[$utrindex] < ($cdsref->[3]-1)) {
			       # UTR interval completely before begin of CDS
			       push @newgbfkey, ($gbfkeyref->[$utrindex-1], $gbfkeyref->[$utrindex]);
			       $utrindex += 2;
			   } elsif ($gbfkeyref->[$utrindex] == ($cdsref->[3]-1)) {
			       # UTR interval ends exactly before begin of CDS
			       push @newgbfkey, $gbfkeyref->[$utrindex-1];
			       for (my $cdsindex = 4; $cdsindex < @$cdsref-1; $cdsindex++) { # multi CDS gene
				   push @newgbfkey, $cdsref->[$cdsindex];
			       }
			       if ($utrindex+2 < @$gbfkeyref && (($cdsref->[-1] + 1) == $gbfkeyref->[$utrindex+1])) {
				   push @newgbfkey, $gbfkeyref->[$utrindex+2];
                                   $utrindex += 4;
			       } else {
				  push @newgbfkey, $cdsref->[-1];
				  $utrindex += 2;
			       }
			   } elsif ($gbfkeyref->[$utrindex-1] > ($cdsref->[-1] + 1)) {
			       # UTR interval ends completely after end of CDS
			       push @newgbfkey, ($gbfkeyref->[$utrindex-1], $gbfkeyref->[$utrindex]);
			       $utrindex += 2;
			   } elsif (($gbfkeyref->[$utrindex-1] == ($cdsref->[-1] + 1))) { # removed $utrindex==4 condition; Katharina Hoff Dec 11th 2017
			       # first UTR interval ends begins directly after end of CDS
			       for (my $cdsindex = 3; $cdsindex < @$cdsref-1; $cdsindex++) {
				   push @newgbfkey, $cdsref->[$cdsindex];
			       }
			       push @newgbfkey, $gbfkeyref->[$utrindex];
			       $utrindex += 2;
			   } else {
			       print "Error: UTR interval in CDS range of gene $genename on sequence $seqname. Ignoring UTR in this case.\n";
			       $fehler = 1;
   			   }
                       }
                       if ($fehler == 0) {
			   #print "add mRNA " , (join " ", @newgbfkey), "\n";
			   $annotation->{"mRNA$genename"}=\@newgbfkey;
		       }
                 }
            }
            delete $annotation->{"UTR$genename"};
	}
    }

    my @nr_sort_fkeyarray = nrsort(values %{$annos{$annotatedseqname}});

    # extract only the cds annotations
    my @cdsfkeys = ();
    my @gene_boundaries = ();
    foreach my $gbfkeyref (@nr_sort_fkeyarray){
	if (@{$gbfkeyref}[0] eq "CDS"){
	    push @cdsfkeys, $gbfkeyref; 
	}
    }

    # check whether there is at least one coding gene
    # if not then output nothing and go to the next sequence
    if (@cdsfkeys <= 0){ next;}
    my $geneNumber = scalar(@cdsfkeys);

    my $currentSeqBegin;
    my $currentSeqEnd;
    my $lastSeqEnd = 0;
    my $genebegin=$length+1;
    my $geneend=0;
    my $count = 0;
    my $i=0;
    my $j;
    my $newseq;
    my $newlength;
    my $offset;
    my @shiftedfeatures;
    my $newname;
    my @theseMRNA = ();
    &printhead($seqname, $length) if ($connected);
    while ($i<@nr_sort_fkeyarray) {
        #print "$i:  " , join (", ", @{$nr_sort_fkeyarray[$i]}), "\n";
        
        # Read up to the next CDS (possibly several mRNA and 1 CDS.)
        # Then determine whether there are more genes and output the LOCUS.
        if ($nr_sort_fkeyarray[$i]->[-1] > $geneend) {
            $geneend = $nr_sort_fkeyarray[$i]->[-1];
        }
        if ($nr_sort_fkeyarray[$i]->[3] < $genebegin) {
	    $genebegin = $nr_sort_fkeyarray[$i]->[3];
	}

        if ($nr_sort_fkeyarray[$i]->[0] eq "CDS") {
	    $currentSeqEnd = $geneend + $flank;
	    if (!$overlap){
		if ($i+1<@nr_sort_fkeyarray && $nr_sort_fkeyarray[$i+1]->[3] < $geneend + 2*$flank && $nr_sort_fkeyarray[$i+1]->[3] > $geneend) {
		    # next gene close but non-overlapping
		    $currentSeqEnd = int(($geneend + $nr_sort_fkeyarray[$i+1]->[3])/2);
		} elsif ($i+1<@nr_sort_fkeyarray && $nr_sort_fkeyarray[$i+1]->[3] <= $geneend ) {
		    # next gene overlapping
		    $currentSeqEnd = $geneend;
		}
	    }
            if ($currentSeqEnd > $length) {
                $currentSeqEnd = $length;
            }
	    $currentSeqBegin = $genebegin - $flank;
	    if ($currentSeqBegin < 1) {
		$currentSeqBegin = 1;
	    }
	    if (!$overlap && $currentSeqBegin < $lastSeqEnd + 1) {
		$currentSeqBegin = $lastSeqEnd + 1;
	    }
	    if ($overlap || $currentSeqBegin < $genebegin) {# UTR not overlapping with previous gene
		$newname = $seqname . "_${currentSeqBegin}-$currentSeqEnd";
		&printhead($newname, $currentSeqEnd-$currentSeqBegin+1) if (!$connected);
		$offset = $currentSeqBegin-1;
		foreach my $theseMRNAref (@theseMRNA) {
		    @shiftedfeatures = @{$theseMRNAref};
		    for($j = 3; $j < @shiftedfeatures; $j++){
			$shiftedfeatures[$j] -= $offset if (!$connected);
		    }
		    &printdata(@shiftedfeatures);
		}
		
		@shiftedfeatures = @{$nr_sort_fkeyarray[$i]};
		for($j = 3; $j < @shiftedfeatures; $j++){
		    $shiftedfeatures[$j] -= $offset if (!$connected);
		}
		$newseq = substr($seq, $currentSeqBegin - 1, $currentSeqEnd - $currentSeqBegin + 1 );
		&printdata(@shiftedfeatures);
		&printseq($newseq, $newlength) if (!$connected);
	    }
	    #reset variables
	    $lastSeqEnd = $currentSeqEnd;
	    $currentSeqBegin = $length+1; #infinity
	    $currentSeqEnd = -1;
	    $count++;
	    $geneend=0;
	    $genebegin=$length+1;
	    @theseMRNA = ();
	} else { # mRNA
	    push @theseMRNA, $nr_sort_fkeyarray[$i];
	}
        $i++;
    }
    &printseq($seq, $length) if ($connected);
} # while <FASTA>

if ($num_ambig_utr_mgs > 0) {
    print "Warning: Had redundant UTR exon information for $num_ambig_utr_mgs genes.\n";
}
if ($seqnameerrcount > 0){
    print "Warning: I assumed $seqnameerrcount times that sequence names end at first space.\n";
}

sub nrsort {
    my @ret = ();    
    my @sorted = sort { 
        # sort by begin position. if two begin positions are equal then 
        # the longer one comes first
        # if two cds are identical then take the one where the corresponding
        # mRNA starts more upstream
	if (@{$a}[3] < @{$b}[3] || (@{$a}[3] == @{$b}[3] && @{$a}[-1]>@{$b}[-1])){
	    return -1;
	} elsif (@{$a}[3] == @{$b}[3] && @{$a}[-1] == @{$b}[-1] && @{$a}[0] eq "CDS" && @{$b}[0] eq "CDS" ) {
	    # CDS range identical, check the mRNA
	    my ($aname, $aTSS, $bname, $bTSS);
	    $aTSS=$bTSS=0;
	    $aname= @{$a}[1];
	    $bname= @{$b}[1];
	    for ($i=0; $i<=$#_; $i++){
		if (@{$_[$i]}[0] eq "mRNA" && @{$_[$i]}[1] eq $aname){
		    if (@{$_[$i]}[2] eq "+"){
			$aTSS=@{$_[$i]}[3];
		    } else {
			$aTSS=-@{$_[$i]}[-1];
		    }
		}
		if (@{$_[$i]}[0] eq "mRNA" && @{$_[$i]}[1] eq $bname){
		    if (@{$_[$i]}[2] eq "+"){
			$bTSS=@{$_[$i]}[3];
		    } else {
			$bTSS=-@{$_[$i]}[-1];
		    }
		}
	    }
	    if (($aTSS>0 && ($aTSS<$bTSS)) || ($aTSS<0 && $aTSS<$bTSS)){
		return -1;
	    }
	    return 0;
	} elsif (@{$a}[3] == @{$b}[3] && @{$a}[-1] == @{$b}[-1] && @{$a}[0] eq "CDS" && @{$b}[0] eq "mRNA" ) {
	    return 1;
	} elsif (@{$a}[3] == @{$b}[3] && @{$a}[-1] == @{$b}[-1] && @{$a}[0] eq "mRNA" && @{$b}[0] eq "CDS" ) {
	    return -1;
	} elsif (@{$a}[3] == @{$b}[3] && @{$a}[-1] == @{$b}[-1] && @{$a}[0] eq "mRNA" && @{$b}[0] eq "mRNA" ) {
	    return 0;
	}
	return 1;
    } @_;
    # now delete overlapping CDS
    my $lastend=0;
    my $end;
    my %geneNames = ();
    $i=0;
    while ($i <= $#sorted){
	if (@{$sorted[$i]}[0] eq "CDS"){
	    if ($overlap || @{$sorted[$i]}[3] > $lastend){
		push @ret, $sorted[$i];
		$lastend = @{$sorted[$i]}[-1];
		$geneNames{@{$sorted[$i]}[1]} = 1;
	    }
	} else {
	    push @ret, $sorted[$i];  
	}
	$i++;
    }
    # remove all mRNAs without CDS
    $i=0;
    while($i < @ret) {
	if (!exists($geneNames{$ret[$i][1]})){
	    splice @ret, $i, 1;
	} else {
	    $i++;
	}
    }

    return @ret;
}

sub printhead {
    my $seqname = shift;
    my $length = shift;
    print OUTPUT "LOCUS       $seqname   $length bp  DNA\n";
    print OUTPUT "FEATURES             Location/Qualifiers\n";
    print OUTPUT "     source          1..$length\n";
}

    #  write dataset to stdout
sub printdata {
    my @gbfkey = @_;

    $fkey = shift @gbfkey;
    $genename = shift @gbfkey;
    $strand = shift @gbfkey;
    if ($strand eq "-") {
	$join = "complement(";
    } else {
	$join = "";
    }
    if (@gbfkey>2) {
	$join .= "join(";
    }
    for ($i=0; $i<@gbfkey; $i+=2) {
	if ($i != 0) {
	    $join .= ",";
	}
	$join .= "$gbfkey[$i]..$gbfkey[$i+1]";
    }
    if (@gbfkey>2) {
	$join .= ")";
    }
    if ($strand eq "-") {
	$join .= ")";
    }
    if ($fkey eq "mRNA"){
	$genename = $genename;
	$genename =~ s/^mRNA//;
	print OUTPUT "     mRNA            ";
    } else {
	$genename = $genename;
	$genename =~ s/^CDS//;
	print OUTPUT "     CDS             ";
    }
    while (length $join > 0) {
	print OUTPUT substr $join,0,59;
	$join = substr $join, 59; 
	if (length $join > 0) {
	    print OUTPUT "\n                     ";
	}
    }
    if (length $genename > 0){
	print OUTPUT "\n                     /gene=\"$genename\"";
    }
    print OUTPUT "\n";
}

sub printseq {

    my $seq = shift;
    my $length = shift;
    if($hardmask){
	$seq =~ s/[actg]/N/g;
    }
    if(not($softmasked)){
	$seq =~ s/A/a/g;
	$seq =~ s/C/c/g;
	$seq =~ s/G/g/g;
	$seq =~ s/T/t/g;
	$seq =~ s/N/n/g;
	$an = $seq =~ s/a/a/g;
	$cn = $seq =~ s/c/c/g;
	$gn = $seq =~ s/g/g/g;
	$tn = $seq =~ s/t/t/g;
	$nn = $seq =~ s/n/n/g;
    }else{
	$an = $seq =~ s/a/a/g;
	$an += $seq =~ s/A/A/g;
	$cn = $seq =~ s/c/c/g;
	$cn += $seq =~ s/C/C/g;
	$gn = $seq =~ s/g/g/g;
	$gn += $seq =~ s/G/G/g;
	$tn = $seq =~ s/t/t/g;
	$tn += $seq =~ s/T/T/g;
	$nn = $seq =~ s/n/n/g;
	$nn += $seq =~ s/N/N/g;
    }
    
    $rest = $length - $an - $cn - $gn - $tn -$nn;
    print OUTPUT "BASE COUNT     $an a   $cn c  $gn g   $tn t";
    if ($nn>0) {
        print OUTPUT "   $nn n";
    }
    if ($rest > 0) {
          print OUTPUT "   $rest ?";
    } 
    print OUTPUT "\nORIGIN\n";
    $i = 1;
    $pos = 0;
    while ($pos <= length $seq) {
        $zahlzeile = "";
        for ($j=0; $j < 9-length "$i"; $j=$j+1) {
            print OUTPUT " ";
        }
        print OUTPUT "$i";
        for ($j=0; $j < 6; $j=$j+1) {
             $ten = substr $seq, $pos, 10;
             if (length $ten > 0) {
                 print OUTPUT " $ten";
             }
             $pos = $pos + 10; #$seq = substr $seq, 10;
        }
        print OUTPUT "\n";
        $i += 60;
    }
    print OUTPUT "//\n";
}
