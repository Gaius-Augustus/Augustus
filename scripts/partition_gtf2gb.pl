#!/usr/bin/perl

# Partition a gff file into gene clusters of a specified minimal number of genes,
# clusters shall not end in another gene.
# Convert the partitioned gtf and a fasta file into a genbank file.
# Katharina Hoff, 21.3.2013

use strict;
use warnings;
use Getopt::Long; # for parameter specification on the command line

my $usage = <<'ENDUSAGE';
parition_gtf2gb.pl      Partition a gtf file into gene clusters of a specified 
    maximal number of genes. Clusters shall not end in 
    another gene. Convert the partitioned gff and a fasta 
    file into a genbank file where one cluster corresponds 
    to one LOCUS. Sequences in the genbank file begin and
    end in the middle of intergenic regions between the 
    clusters. Gene features in the genbank file may 
    overlap. Only single exon genes without UTR are 
    currently supported.

    SYNOPSIS

    parition_gtf2gb.pl --genome=genome.fa --gtf=annotation.gtf --out.gb

    genome.fa is a genome file in (multiple) fasta format. The fasta headers 
    must correspond to the sequence names in the gff-file. Example:
    
    >entry1
    ATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATA
    ATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATA
    >entry2
    ATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATATATNATATATCGNATA
    ATCGNATATATNATATATCGNATATATNATAT
    ...

    annotation.gtf is a file that specifies the coordinates of single-exon 
    genes without UTR in tabulator separated gtf format. Example:
    
    entry1   NCBI   CDS 1   6   0 + . gene_id "g1";transcript_id "g1.t1"
    entry1   NCBI   CDS 19  30  0 - . gene_id "g2"; transcript_id "g2.t1"
    entry2   NCBI   CDS 7   28  0 + . gene_id "g3"; transcript_id "g3.t1"
    ...

    out.gb is the output file with clusters of CDS per LOCUS in genbank
    format. Example: 

    LOCUS       NC_010473.1_5128-6494   1367 bp  DNA
    FEATURES             Location/Qualifiers
             source          1..1367
             CDS             107..402
                             /gene="170079668"
	     CDS             complement(556..1332)
	                     /gene="170079669"
	BASE COUNT     346 a   349 c  336 g   336 t
	ORIGIN
	        1 aatgcggtaa cttagagatt aggattgcgg agaataacaa ccgccgttct catcgagtaa
	       61 tctccggata tcgacccata acgggcaatg ataaaaggag taacctgtga aaaagatgca
	     ...
	     1261 cagcgtatag cgcgtggtgg tcaacgggct ttggtaatca agcgttttcg caggtgaaat
	     1321 aagaatcagc atatccagtc cttgcaggaa atttatgccg
	//


OPTIONS

    --help                     output this help message
    --minClusterSize=n         minimal number of genes in a one cluster, 
                               default n=20
    --maxFlankingRegionSize=n  maximal size of flanking intergenic region to be 
                               excised around clusters of genes,
                               default n=2000
    --verbose		       print verbose warning messages

ENDUSAGE

my $help = 0;
my $verbose = 0;
my ($genomeFname, $gtfFname, $outName);
my $minClusterSize = 20;
my $maxFlankSize = 2000;
my $msg;

GetOptions('genome=s'=>\$genomeFname,
           'gtf=s'=>\$gtfFname,
	   'out=s'=>\$outName,
           'minClusterSize:n'=>\$minClusterSize,
	   'maxFlankingRegionSize:n'=>\$maxFlankSize,
	   'verbose!'=>\$verbose,
           'help!'=>\$help);

if($help){
	print $usage;
	exit(1);
}

if(!defined($genomeFname)){
	print "\nMissing input genome file\n\n$usage";
	exit(1);
}

if(!defined($gtfFname)){
	print "\nMissing input gtf file\n\n$usage";
	exit(1);
}

if(!defined($outName)){
	$msg = "\nWARNING: Missing output file name. Results will be printed to out.gb\n";
	print STDERR $msg;
	$outName = "out.gb";
}

# main body of code starts

# read sequence and annotation;
my %fasta = read_genome($genomeFname, $verbose);
my %gtf = read_gtf($gtfFname, \%fasta, $verbose); # note how the reference of the hash is passed to the function!
# partition genes into clusters
my %partitions = identify_cluster_boundaries(\%fasta, \%gtf, $minClusterSize, $maxFlankSize, $verbose);

# write genbank format file
open(my($OUT), ">", $outName) or die("Could not open output file $outName!\n");
write_gb_format(\%fasta, \%gtf, \%partitions, $OUT, $verbose);
close($OUT) or die("Could not close output file $outName!\n");

#----------------   subroutines     -------------------

# read genome file into a hash
# remove new lines and empty newlines
# die if a header is used twice
# issue a warning if the sequence does not consist of valid nucleotide letters
# returns a hash with sequence IDs as keys and another hash as values. The inner
# hash contains a field sequence and a field seqLen.
sub read_genome {
    my $fname = shift; # genome file name
    my $verbose = shift; # verbosity flag for printing warnings
    my $msg; # print message variable
    $msg = "Reading genome file...\n";
    print STDERR $msg if($verbose);
    my %seq = (); # hash that will be returned
    open(FASTA, "<", $fname) or die ("Could not open genome input file $fname.\n");
    my $seqName; # temporary variable to store the sequence name
    while(<FASTA>){
	if($_=~m/^>/){
	    chomp $_; # cleave off the new line and maybe whitespaces
	    $_=~s/^>//; # remove the leading header sign
	    $seqName=$_;
	    if(defined($seq{$seqName})){
		die("ERROR: the same sequence name occured twice in the genome file $fname!\n");
	    }
	}else{	
	    chomp $_;
	    if($_=~m/([^AaTtCcGgNn])/){ 
		# other valid nucleotides would have been XxUuYyRrWwSsKkMmBbDdHhVvPpZz, 
		#but augustus does not like those, anyway, and they are rare
		$msg = "WARNING: Genomic nucleotide sequence contains non nucleotide character $1!\n";
		$msg .= "         All non ATGC nucleotides will be replaced by n!\n";
		print STDERR $msg if($verbose);
		$_=~s/[XxUuYyRrWwSsKkMmBbDdHhVvPpZz]/n/g;
	    }
	    if(not(defined($seq{$seqName}))){
		$seq{$seqName}{sequence} = $_;
	    }else{
		$seq{$seqName}{sequence} = $seq{$seqName}{sequence}.$_;
	    }
	}
    }
    close(FASTA) or die ("Could not close genome input file $fname.\n");
    while ( $seqName = each(%seq) ) {
	$seq{$seqName}{seqLen} = length($seq{$seqName}{sequence});
    }
    $msg = "Done!\n";
    print STDERR $msg if($verbose);
    return %seq;
}

# read gtf file into a hash of arrays of hashes
# currently implementation for single-exon genes without UTR, only!
# one would use a different data structure for multi-exon genes.
# issue a warning if no sequence was provided for a gtf line
# issue a warning if the file does not have 8 columns
# issue a warning if any of the important columns are not formatted correctly
# die if input file was not sorted by start coordinate
# returns a hash where the sequence names serve as keys, and arrays as entries.
# the inner arrays have one entry per gene, and each of those entrys is
# another hash with the fields feature, start, stop, strand, and id.
sub read_gtf {
    my $fname = shift; # genome file name
    my $genomeRef = shift; # reference to genome fasta hash!
    my $verbose = shift; # verbosity flag for printing warnings
    my $msg; # print message variable
    $msg = "Reading gtf file...\n";
    print STDERR $msg if($verbose);
    my %gtf = (); # hash of arrays that will be returned
    open(GTF, "<", $fname) or die ("Could not open gtf input file $fname.\n");
    my @t; # temporary split array
    my $line; # temporary line variable
    my $w1 = 0; # counters to suppress printing of different warnings after 5 examples
    my $w2 = 0;
    my $w3 = 0;
    my $w4 = 0;
    my $w5 = 0;
    my $seqName;
    while(<GTF>){
	$line = $_;
	@t = split(/\t/, $_); # obtain features that will be required
	$seqName = $t[0];
	# if the gtf line is formatted correctly, included in the genomic sequence, and 
	# not a comment or an empty newline	
	if($genomeRef->{$seqName}{seqLen} > 1 && $#t == 8 && not($_=~m/\#/) && not($_=~m/^\n/) && ($_=~m/\tCDS\t\d+\t\d+\t[\d+\.]\t[+-]\t[\d+\.]\ttranscript_id "\w+\"/) && ($t[4]<= $genomeRef->{$seqName}{seqLen})){ 
            # note how to work with the hash reference
	    # store gtf line for further usage
	    push (@{$gtf{$t[0]}}, {feature => $t[2], start => $t[3], stop => $t[4], strand => $t[6], id => $t[8]});
	    # if the gtf line refers to a nonexisting sequence, issue a warning and ignore entry
	}elsif(not(defined($genomeRef->{$t[0]}))){
	    $msg = "WARNING: gtf file $fname contains entries that refer to sequences that were not included\n";
	    $msg .= "        in the provided fasta file! Entry\n        ".$line."        will be ignored!\n";
	    print STDERR $msg if($w1<6 && $verbose);
	    $w1++;
	    supress($w1) if($verbose);
	    # if the number of columns is not 9, something must be very wrong. Ignore line, issue a warning.
	}elsif($#t != 9){
	    $msg = "WARNING: entry in gtf file\n        $_        does not have 9 columns. ";
	    $msg .= "Entry will be ignored!\n";
	    print STDERR $msg if($w2<6 && $verbose);
	    $w2++;
	    supress($w2) if($verbose);
	    # if something else is weird in the content of the columns, issue warning, print weird line, ignore.
	}elsif(not($_=~m/\tCDS\t\d+\t\d+\t[\d+\.]\t[+-]\t[\d+\.]\ttranscript_id "\w+\"/)){
	    $msg = "WARNING: entry in gtf file $fname\n".$line."is not compliant with gtf format\n";
	    $msg .= "        as it is required for this script. Correctly formatted,\n";
	    $msg .= "        the file contains the features CDS in column 3, numeric values in columns 4\n";
	    $msg .= "        and 5, a + or - in column 7, and a transcript_id with\n";
	    $msg .= "        quotation marks in the last column! This entry will be ignored!\n";
	    print STDERR $msg if($w3 < 6 && $verbose);
	    $w3++;
	    supress($w3) if($verbose);
	    # if the gtf feature exceeds the range of sequence length, issue warning, ignore.
	}elsif($t[4] > $genomeRef->{$t[0]}{seqLen}){
	    $msg = "WARNING: gtf file $fname contains entries that do not correspond to the genome sequence\n";
	    $msg .= "        because the genomic sequences are shorter!\n        $line";
	    $msg .= "        This entry will be ignored!\n";
	    print STDERR $msg if($w5 < 6 && $verbose);
	    $w5++;
	    supress($w5) if($verbose);
	}else{
	    # if some other case might occur that I forgot to implement:
	    $msg = "WARNING: Entry\n        $line"."        was not imported due to an unforeseen format error!\n";
	    print STDERR $msg if($w4 < 6 && $verbose);
	    $w4++;
	    supress($w4) if($verbose);
	}
    }
    # check whether gtf file was sorted
    my $prev; # temporary variable for sortedness check
    for my $contig ( keys %gtf ) {
	$prev = 1;
	@t = @{ $gtf{$contig} };
	foreach(@{ $gtf{$contig} }){ # note how data access in hash of arrays of hashes!
	    $msg = "ERROR: gtf file $fname must be sorted by start coordinate!\n";
	    $msg .= "       Run for example \"cat unsorted.gtf | sort -n -k4,4 > sorted.gtf\"\n";
	    die($msg) if(not($$_{start} >= $prev));
	    $prev = $$_{start}
	}
    }
    $msg = "Done!\n";
    print STDERR $msg if($verbose);
    return %gtf;
}

# issue message about supress printing of warning message after 5 examples
sub supress {
    my $w = shift;
    if($w == 6){
	print STDERR "This warning will not be shown for further examples.\n";
    }
}

# identify cluster boundaries for each LOCUS to be printed later
# returns a hash of arrays
# positions start counting at 1 for the first position
sub identify_cluster_boundaries {
    my $genomeRef = shift;
    my $gtfRef = shift;
    my $clusterSize = shift;
    my $flankSize = shift;
    my $verbose = shift;
    my $msg;
    $msg = "Identifying gene cluster borders...\n";
    print $msg if($verbose);
    my %boundaries; # array with partition boundaries to be returned
    my ($lower, $upper);
    my $nEle;
    my $gC; # gtf cluster counter
    # loop over all fasta entries that have an annotation
    for my $contig ( keys %$gtfRef) {
	$nEle = @{$$gtfRef{$contig}};
	# set first lower boundary
	if(${$$gtfRef{$contig}}[0]{start} < $flankSize){
	    $lower=1;
	}else{
	    $lower = ${$$gtfRef{$contig}}[0]{start} - $flankSize;
	}
	# loop over all remaining gtf entries of one fasta entry except for the last to identify cluster borders
	$gC = 0;
	for(my $i = 1; $i <= $nEle-2; $i++){
	    $gC++;
	    if($gC >= $clusterSize){
		# if genes do not overlap, and or not adjacent without spacer and if the overlap is smaller than 
		# twice the max flanking region, cassette genes shall not be placed in different clusters
		if((${$$gtfRef{$contig}}[$i+1]{start} - ${$$gtfRef{$contig}}[$i]{stop} > 1) && (${$$gtfRef{$contig}}[$i+1]{start} - ${$$gtfRef{$contig}}[$i]{stop} < 2*$flankSize) && (${$$gtfRef{$contig}}[$i]{stop} > ${$$gtfRef{$contig}}[$i-1]{stop})){
		    $upper = int( ${$$gtfRef{$contig}}[$i]{stop} + (${$$gtfRef{$contig}}[$i+1]{start} - ${$$gtfRef{$contig}}[$i]{stop})/2 + 0.5);
		    push(@{$boundaries{$contig}}, {lower => $lower, upper => $upper});
		    $lower = $upper + 1;
		    $gC = 0;
		    # if genes do not overlap but have a very large intergenic region, cassette genes shall not be 
		    # placed in different clusters
		}elsif(${$$gtfRef{$contig}}[$i+1]{start} - ${$$gtfRef{$contig}}[$i]{stop} > 1 && (${$$gtfRef{$contig}}[$i]{stop} > ${$$gtfRef{$contig}}[$i-1]{stop})){
		    $upper = ${$$gtfRef{$contig}}[$i]{stop} + $flankSize;

		    push(@{$boundaries{$contig}}, {lower => $lower, upper => $upper});
		    $lower = ${$$gtfRef{$contig}}[$i+1]{start} - $flankSize;
		    $gC = 0;
		}
	    }
	}
	# set last upper boundary
	if(${$$gtfRef{$contig}}[$nEle-1]{stop} < $genomeRef->{$contig}{seqLen} - $flankSize){
	    $upper = ${$$gtfRef{$contig}}[$nEle-1]{stop} + $flankSize;
	}else{
	    $upper = $genomeRef->{$contig}{seqLen};
	}
	push(@{$boundaries{$contig}}, {lower => $lower, upper => $upper});					
    }
    $msg = "Done!\n";
    print $msg if($verbose);
    return %boundaries;
}

# write genbank format
# use fasta hash reference to retrieve sequence
# use paritions reference to define borders of LOCI
# use gtf reference to identify suitable CDS entries for each LOCUS
# no return value
sub write_gb_format{
    my $genomeRef = shift;
    my $gtfRef = shift;
    my $partitionsRef = shift;
    my $OUT = shift;
    my $verbose = shift;
    my $msg;
    $msg = "Writing genbank format...\n";
    print $msg if($verbose);
    my $gC; # gene Counter for array
    my $gtfSize;
    for my $contig ( keys %$gtfRef) {
	$gC = 0;
	$gtfSize = @{$$gtfRef{$contig}};
	foreach(@{$$partitionsRef{$contig}}){
	    print $OUT "LOCUS       $contig"."_$$_{lower}-$$_{upper}   ".($$_{upper}-$$_{lower}+1)." bp  DNA\n";
	    print $OUT "FEATURES             Location/Qualifiers\n";
	    print $OUT "     source          1..".($$_{upper}-$$_{lower}+1)."\n";
	    # print single gene CDS features
	    for (my($i) = $gC; $i <= ($gtfSize-1); $i++) {
		print  $OUT "     CDS             ";
		if(${$$gtfRef{$contig}}[$i]{strand} eq "+"){
		    print $OUT "".(${$$gtfRef{$contig}}[$i]{start}-$$_{lower}+1)."..".(${$$gtfRef{$contig}}[$i]{stop}-$$_{lower}+1)."\n";
		}else{
		    print $OUT "complement(".(${$$gtfRef{$contig}}[$i]{start}-$$_{lower}+1)."..".(${$$gtfRef{$contig}}[$i]{stop}-$$_{lower}+1).")\n";
		}
		${$$gtfRef{$contig}}[$i]{id}=~m/transcript_id \"(\w+)\"/;
		print $OUT "                     /gene=\"$1"."\"\n";
		if(defined(${$$gtfRef{$contig}}[$i+1]{start})){ # don't continue after the last entry...
		    if(${$$gtfRef{$contig}}[$i+1]{start} >= $$_{upper}){$gC = $i+1; last;}
		}
	    }
	    printseq(substr(${$$genomeRef{$contig}}{sequence}, ($$_{lower}-1), ($$_{upper}-$$_{lower}+1)), $OUT);

	}
    }
    $msg = "Done!\n";
    print $msg if($verbose);
}

# print DNA sequence
sub printseq {
    my $seq = shift;
    my $seqLen = length($seq);
    my $OUT = shift;
    # make all DNA letters lower case
    $seq =~ s/A/a/g;
    $seq =~ s/C/c/g;
    $seq =~ s/G/g/g;
    $seq =~ s/T/t/g;
    $seq =~ s/N/n/g;
    # count nucleotide quantities
    my $an = $seq =~ s/a/a/g;
    my $cn = $seq =~ s/c/c/g;
    my $gn = $seq =~ s/g/g/g;
    my $tn = $seq =~ s/t/t/g;
    my $nn = $seq =~ s/n/n/g;
    my $rest = $seqLen - $an - $cn - $gn - $tn -$nn;
    # print nucleotide counts
    print $OUT "BASE COUNT     $an a   $cn c  $gn g   $tn t";
    if ($nn>0) {
        print $OUT "   $nn n";
    }
    # remaining nucleotides should actually not occur at all
    if ($rest > 0) {
	# this is an error that we always want to see!
	print STDERR "ERROR: Number of persisting unknown nucleotides is $rest\n";
    } 
    print $OUT "\nORIGIN\n";
    my $i = 1;
    my $pos = 0;
    my ($j, $ten, $zahlzeile);
    while ($pos <= $seqLen && $i <= $seqLen) {
        $zahlzeile = "";
        for ($j=0; $j < 9-length "$i"; $j=$j+1) {
            print $OUT " ";
        }
        print $OUT "$i";
        for ($j=0; $j < 6; $j=$j+1) {
	    if(($pos+9) <= $seqLen){
		$ten = substr $seq, $pos, 10;
	    }elsif(($seqLen-$pos+1) > 0){
		$ten = substr $seq, $pos, ($seqLen-$pos+1);
		last;
	    }
	    print $OUT " $ten";
	    $pos = $pos + 10; # $seq = substr $seq, 10;
        }
        print $OUT "\n";
        $i += 60;
    }
    print $OUT "//\n";
}

