#!/usr/bin/perl
# convert a utr.gff from makeUtrTrainingSet.pl to a gbrowse file
# Mario Stanke, 9.4.2008
#
# example input:
#
# chromosome1     makeUtr 5'-UTR  43947   44004   128     -       .       transcript_id "g6.t1"; gene_id "g6";
# chromosome1     makeUtr 5'-UTR  44317   44476   128     -       .       transcript_id "g6.t1"; gene_id "g6";
# chromosome1     makeUtr 3'-UTR  43022   43173   46      -       .       transcript_id "g7.t1"; gene_id "g7";
# chromosome1     makeUtr 3'-UTR  43436   43609   46      -       .       transcript_id "g7.t1"; gene_id "g7";
# chromosome1     makeUtr 3'-UTR  43837   43995   46      -       .       transcript_id "g7.t1"; gene_id "g7";
#
# example output:
#
# chromosome1     makeUtr 5UTREX  43947   44004   128     -       .       5UTR g6.t1
# chromosome1     makeUtr 5UTREX  44317   44476   128     -       .       5UTR g6.t1
# chromosome1     makeUtr 5UTR    43947   44476   128     -       .       5UTR g6.t1
# chromosome1     makeUtr 3UTREX  43022   43173   46      -       .       3UTR g7.t1
# chromosome1     makeUtr 3UTREX  43436   43609   46      -       .       3UTR g7.t1
# chromosome1     makeUtr 3UTREX  43837   43995   46      -       .       3UTR g7.t1
# chromosome1     makeUtr 3UTR    43022   43995   46      -       .       3UTR g7.t1

my ($gene, $oldgene);
my ($type, $newtype);
my $start = -1;
my $end = -1;
my $score = -1;
my $seqname;
my $source;
my $strand;

while(<>){
    @f = split /\t/;
    next unless @f>7;
    if ($f[8] =~ /transcript_id "(\S+)"/) {
	$newgene = $1;
	if ($f[2] =~ /^5/) {
	    $newtype = 5;
	} else {
	    $newtype = 3;
	}
     
	if ((defined($gene) && defined($type)) && ($gene ne $newgene  || $type ne $newtype)) {
	    if ($type == 5) {
		print "$seqname\t$source\t5UTR\t$start\t$end\t$score\t$strand\t.\t5UTR $gene\n";
	    } else {
		print "$seqname\t$source\t3UTR\t$start\t$end\t$score\t$strand\t.\t3UTR $gene\n";
	    }
	    $start = $end = -1;
	    
	}	
	if ($newtype == 5) {
	    print "$f[0]\t$f[1]\t5UTREX\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t5UTR $newgene\n";
	} else {
	    print "$f[0]\t$f[1]\t3UTREX\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t3UTR $newgene\n";
	}
	$seqname = $f[0];
	$source = $f[1];
	$start = $f[3] unless ($start != -1 && $start < $f[3]);
	$end   = $f[4] unless ($end != -1   && $end   > $f[4]);
	$score = $f[5];
	$strand = $f[6];
	$gene = $newgene;
	$type = $newtype;
    }
}
