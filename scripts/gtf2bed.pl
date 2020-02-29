#!/usr/bin/env perl
#
# convert gtf/gff/gff3 to BED format
# 23.9.2016 Stefanie Koenig

use strict;
use Getopt::Long;


my $help = 0;
my $verbose = 0;
my $itemRgb="0,0,255";
my $includeStopInCDS = 0;

GetOptions(
    'itemRgb=s'=>\$itemRgb,
    'help!'=>\$help,
    'verbose!'=>\$verbose,
    'includeStopInCDS!'=>\$includeStopInCDS);

exec("perldoc $0") if ($help);

my @txorder = (); # tx ids in the order of the input file
my %geneOf = ();  # keys tx ids, values: gene ids
my %geneLine =  (); # keys gene ids, values: array refs for the gene GTF line (if exists)
# hash of transcripts
#    keys: transcript ids
#    values: hash reference
#            keys: txstart
#                  txend
#                  codingstart
#                  codingend
#                  strand, chr, source
#                  txline  array of columns if a transcript/mRNA line exists
#                  CDS  array of arrays (lines and columns) for coding parts of exons
#                  UTR  array of arrays for UTR exons
#                  exon array of arrays for complete exons
#                  intron array of arrays for introns
#                  rest array of arrays all other features, like tts,tss start_codon, stop_codon

my %txs = ();

parseAndStoreGTF();
convert();
printBed();

sub parseAndStoreGTF{
    my %seen = ();
    my ($txid, $geneid, $chr, $start, $end, $feature, $strand, $source, $stop_codon);
    foreach my $line (<STDIN>){
	my @f = split /\t/, $line;
	next if (@f<8);
	($chr,$source,$feature,$start,$end,$strand) = ($f[0],$f[1],$f[2],$f[3],$f[4],$f[6]);
	# check whether it is a line with 'gene' feature
	if ($f[2] eq "gene" && ($f[8] =~ /ID=([^;]+)/ || $f[8] =~ /gene_id."?([^";]+)"?/ || $f[8] =~ /^(\S+)$/)){
	    $geneid = $1;
	    $geneLine{$geneid} = \@f;
	    next;
	} 
	# check whether it is a line with 'transcript' feature
	if ($f[2] =~ "(transcript|mRNA)" && ($f[8] =~ /ID=([^;]+)/ || $f[8] =~ /transcript_id."?([^";]+)"?/ || $f[8] =~ /^(\S+)$/)){
	    $txid = $1;
	    $txs{$txid} = {"strand"=>$strand, "chr"=>$chr, "source"=>$source, "CDS"=>[], "UTR"=>[], "exon"=>[], "intron"=>[], "rest"=>[]} if (!exists($txs{$txid}));
	    $txs{$txid}{"txline"} = \@f;
	    if($f[8] =~ /Parent=([^;]+)/ || $f[8] =~ /gene_id."?([^";]+)"?/){
		$geneOf{$txid} = $1;
	    }
	    next;
	}
	
	$txs{$txid}{"CDS"} = [] if (!defined($txs{$txid}{"CDS"}));
	$txs{$txid}{"UTR"} = [] if (!defined($txs{$txid}{"UTR"}));
	$txs{$txid}{"exon"} = [] if (!defined($txs{$txid}{"exon"}));
	$txs{$txid}{"rest"} = [] if (!defined($txs{$txid}{"rest"}));

	# all other lines must belong to a transcript and a gene
	if ($f[8] =~ /(transcript_id|Transcript)."?([^";]+)"?/ ){
	    $txid = $2;
	} else {
	    if($f[8] =~ /Parent=([^;]+)/){
		$txid = $1;
	    }else{
		die ("Neither GTF nor GFF format in the following line:\n$line\ntranscript_id not found.\n");
	    }
	}
	if ($f[8] =~ /gene_id."?([^";]+)"?/){
	    $geneid = $1;
	} else {
	    if($f[8] =~ /Parent=([^;]+)/){
		$geneid = $geneOf{$1};
	    }else{
		die ("Neither GTF nor GFF format in the following line:\n$line\ngene_id not found.\n");
	    }
	}
	$txs{$txid} = {"strand"=>$strand, "chr"=>$chr, "source"=>$source, "CDS"=>[], "UTR"=>[], "exon"=>[], "intron"=>[], "rest"=>[]} if (!exists($txs{$txid}));

	if (!$seen{$txid}){
	    push @txorder, $txid; # remember the input order for transcripts for the output
	    $seen{$txid} = 1;
	}
	# assign parental gene id to tx id
	die ("transcript $txid has conflicting gene parents: and $geneid. Remember: In GTF txids need to be overall unique.")
	    if (defined($geneOf{$txid}) && $geneOf{$txid} ne $geneid);
	
	if ($feature eq "CDS" || $feature eq "coding_exon" || $feature eq "exon" || $feature =~ /UTR/i){
	    $txs{$txid} = {"strand"=>$strand, "chr"=>$chr, "source"=>$source, "CDS"=>[], "UTR"=>[], "exon"=>[], "intron"=>[], "rest"=>[]} if (!exists($txs{$txid}));
	    $txs{$txid}{"txstart"} = $start if (!defined($txs{$txid}{"txstart"}) || $txs{$txid}{"txstart"} > $start);
	    $txs{$txid}{"txend"} = $end if (!defined($txs{$txid}{"txend"}) || $txs{$txid}{"txend"} < $end);
	}
	if ($feature eq "CDS" || $feature eq "coding_exon"){
	    $txs{$txid}{"codingstart"} = $start if (!defined($txs{$txid}{"codingstart"}) || $txs{$txid}{"codingstart"} > $start);
	    $txs{$txid}{"codingend"} = $end if (!defined($txs{$txid}{"codingend"}) || $txs{$txid}{"codingend"} < $end);
	    push @{$txs{$txid}{"CDS"}}, \@f;
	} elsif ($feature =~ /UTR/i){
	    push @{$txs{$txid}{"UTR"}}, \@f;
	} elsif ($feature eq "exon"){
	    push @{$txs{$txid}{"exon"}}, \@f;
	} elsif ($feature eq "intron") {
	    $txs{$txid}{"intron"} = [] if (!defined($txs{$txid}{"intron"}));
            push @{$txs{$txid}{"intron"}}, \@f;
	} else {
	    push @{$txs{$txid}{"rest"}}, \@f;
	}
	if ($feature eq "stop_codon"){
	    $txs{$txid}{"stop_codon"} = $start
	}
    }
}


sub convert{
    my @f;
    foreach my $txid (keys %txs){
	
        # optionally, include stop codon in CDS
	if($includeStopInCDS){
	    my @cdslines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} @{$txs{$txid}{"CDS"}};
	    if( $txs{$txid}{"strand"} eq '-' && $txs{$txid}{"stop_codon"}+3 == $txs{$txid}{"codingstart"}){
		$cdslines[0]->[3]-=3;
	    }
	    if( $txs{$txid}{"strand"} eq '+' && $txs{$txid}{"stop_codon"} == $txs{$txid}{"codingend"}+1){
		$cdslines[$#cdslines]->[4]+=3;
	    }
	    # TODO: if stop_codon is a separate CDS exon, then insert new CDS exon
	}
	
	# remember whether exons were not in the input file
	my $exonArrayWasEmpty;
	if(@{$txs{$txid}{"exon"}} == 0){
	    $exonArrayWasEmpty = 1;
	}else{
	    $exonArrayWasEmpty = 0;
	}
	# add exon lines if not already present and if desired in output
	if (@{$txs{$txid}{"exon"}} == 0){
	    print "Creating exon lines for $txid\n" if ($verbose);
	    # sort UTR and CDS lines by coordinates
	    my @exonpartlines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} (@{$txs{$txid}{"CDS"}}, @{$txs{$txid}{"UTR"}});
	    next if (@exonpartlines == 0);
	    @f = @{$exonpartlines[0]};
	    shift @exonpartlines;
	    ($f[2], $f[5], $f[7]) = ("exon", '.', '.'); # score and frame are not defined
	    foreach my $g (@exonpartlines){
		if ($f[4] >= $g->[3]){ # check for non-overlappingness
		    die ("In transcript $txid two UTR/CDS features are overlapping. Not allowed by definition.");
		} elsif ($f[4] + 1 == $g->[3]){ # exactly adjacent
		    # join two UTR/CDS features to one
		    $f[4] = $g->[4];
		} else {
		    # push exon
		    my @ff = @f; # deep copy array
		    push @{$txs{$txid}{"exon"}}, \@ff;
		    @f = @$g;
		    ($f[2], $f[5], $f[7]) = ("exon", '.', '.'); # score and frame are not defined
		}
	    }
	    # push remaining, last exon
	    my @ff = @f;
	    push @{$txs{$txid}{"exon"}}, \@ff;
	}

    }
}

sub printBed {

    foreach my $txid (@txorder){

	my @lines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} @{$txs{$txid}{"exon"}};
	my $blockCount = @lines;
	my @blockSizes = ();
	my @blockStarts = ();

	foreach my $line (@lines){
	    push @blockSizes, $line->[4] - $line->[3] + 1;
	    push @blockStarts,  $line->[3] - $txs{$txid}{"txstart"};
	}
	
	print $txs{$txid}{"chr"} . "\t";
	print $txs{$txid}{"txstart"}-1 . "\t";
	print $txs{$txid}{"txend"} . "\t";
	print $txid . "\t0\t";
	print $txs{$txid}{"strand"} . "\t";
	if(@{$txs{$txid}{"CDS"}} == 0){
	    print $txs{$txid}{"txstart"}-1 . "\t";
	    print $txs{$txid}{"txend"} . "\t";
	}
	else{
	    print $txs{$txid}{"codingstart"}-1 . "\t";
	    print $txs{$txid}{"codingend"} . "\t";
	}
	print "$itemRgb\t";
	print "$blockCount\t";
	print join(',',@blockSizes) ."\t";
	print join(',',@blockStarts) ."\n";
    }
}

__END__

=pod

=head1 NAME

gtf2bed.pl      convert gtf/gff/gff3 to BED format

=head1 SYNOPSIS

gtf2bed.pl <in.gtf >out.bed
    
=head1 OPTIONS

  --itemRgb=s              a string s encoding the RGB value of the form R,G,B (default 0,0,225).
  --includeStopInCDS       include stop codon into the coding sequence (default off)

=head1 DESCRIPTION
    
    example input:
    chr16   AUGUSTUS        transcript      100472  160062  .       -       .       jg7.t1
    chr16   AUGUSTUS        tts     100472  100472  .       -       .       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        UTR     100472  101158  0       -       .       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        stop_codon      101159  101161  .       -       0       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     101159  101285  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     102477  102644  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     104227  104334  0.9     -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     106525  106755  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     110093  110263  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     110759  111288  1       -       0       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     117341  117478  1       -       0       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     123010  123106  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     127580  127720  1       -       1       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        CDS     159185  160062  0       -       0       transcript_id "jg7.t1"; gene_id "jg7";
    chr16   AUGUSTUS        start_codon     160060  160062  .       -       0       transcript_id "jg7.t1"; gene_id "jg7";

    example output :

    chr16   100471  160062  jg7.t1  0       -       101158  160062  0,0,255 10      814,168,108,231,171,530,138,97,141,878  100471,102476,104226,106524,110092,110758,117340,123009,127579,159184

=cut
