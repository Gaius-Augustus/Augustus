#!/usr/bin/perl

###################################################################################################
# evalCGP
# evaluates a prediction in GTF format against an annotation
# using the external evaluation package Eval by Evan Keibler and Michael R. Brent¹
# and returns accuracy values (SN and SP on gene, exon and nucleotide level)
# evalCGP only compares gene features on given genomic intervals that
# are parsed from the prediction files.
#
# ¹(Eval: A software package for analysis of genome annotations. BMC Bioinformatics 4: 50 (2003))
#
# usage
#
# evalCGP.pl --anno=annotation.gtf --pred=prediction.gtf
#
# Stefanie Koenig, 11.08.2014
###################################################################################################


use strict;
use IO::File;

my %cmdpars = ( 'pred'              => '',
		'anno'              => '',
		'joingenes'         => '',
		'wholeGenome'       => '',
                'alternatives'      => '',
		'noselection'       => '',
		'nojoin'            => '',
		'jg_exec_dir'       => '',
		'eval_exec_dir'     => '');


my $usage = <<'ENDUSAGE';

evalCGP.pl      evaluates a prediction in GTF format against an annotation
                using the external evaluation package Eval by Evan Keibler and Michael R. Brent
                and returns accuracy values (SN and SP on gene, exon and nucleotide level)
                evalCGP only compares gene features on given genomic intervals that
                are parsed from the prediction files.

USAGE

evalCGP.pl --anno=annotation.gtf --pred=prediction.gtf

      annotation.gtf           Annotation file in GTF format.
      prediction.gtf           Prediction file in GTF format.

OPTIONS

    --eval_exec_dir=d          Directory that contains the executable evaluate_gtf.pl from the eval package.
                               If not specified it must be in \$PATH environment variable.
    --joingenes=1              Use this option to merge genes in the prediction set and filter out duplicates (default: 0)
    --wholeGenome=1            If this flag is set evaluation is on the whole genome. Per default, evaluation
                               is restricted to the gene ranges
    --alternatives=1           Parameter of joingenes. If this flag is set, joingenes keeps alternative splice forms of a gene, otherwise
                               it only keeps the best splicing form. Per definition, alternative splice forms are either transcripts
                               with the same gene ID or the same coding start AND end coordinates (default: 0).
    --noselection=1            Parameter of joingenes. If this flag is set, joingenes does NOT select a single best transcripts
                               among multiple conflicting transcripts. Two transcripts are confliciting if they overlap
                               each other and are no alternative splice forms.
                               considered as conflicting.
    --nojoin=1                 Parameter of joingenes. If this flag is set, joingenes does NOT create new
                               transcripts by merging input transcripts, f.i. it does NOT combine two
                               incomplete transcripts to a single complete transcript, where possible.
                               


ENDUSAGE

##############################################################
# Check the command line
##############################################################

if ($#ARGV<0) {
    print "$usage";
    exit;
}

foreach (@ARGV) {
    if (/--(\w+)=(.*)/){
	if (!exists($cmdpars{$1})){
	    print "unknown parameter: " . $1 . "\n$usage";
	    exit;
	}
	$cmdpars{$1}=$2;
    } 
}
if ($cmdpars{"pred"} eq ""){
    print "prediction file missing\n$usage";
    exit;
}
if ($cmdpars{"anno"} eq ""){
    print "annotation file missing\n$usage";
    exit;
}
if ($cmdpars{'eval_exec_dir'} =~ /.[^\/]$/) {
    $cmdpars{'eval_exec_dir'} .= '/';
}

if ($cmdpars{'jg_exec_dir'} =~ /.[^\/]$/) {
    $cmdpars{'jg_exec_dir'} .= '/';
}

my $joingenes=0;
if ($cmdpars{'joingenes'} eq '1'){
    $joingenes=1;
}

my $wholegenome=0;
if ($cmdpars{'wholeGenome'} eq '1'){
    $wholegenome=1;
}

my $jg_pars="";
if ($cmdpars{'alternatives'} eq '1'){
    $jg_pars.=" -a";
}
if ($cmdpars{'noselection'} eq '1'){
    $jg_pars.=" -l";
}
if ($cmdpars{'nojoin'} eq '1'){
    $jg_pars.=" -j";
}

# check whether joingenes is properly installed
if ($joingenes && qx(which "$cmdpars{'jg_exec_dir'}joingenes") !~ /joingenes$/){
    die ("joingenes is not executable. Please add the directory which contains the executable joingenes to the PATH environment variable or specify the path with --jg_exec_dir.");
} 

# check whether the eval package is properly installed
if (qx(which "$cmdpars{'eval_exec_dir'}evaluate_gtf.pl") !~ /evaluate_gtf.pl$/){
    die ("eval is not executable. Please add the directory which contains the executable evaluate_gtf.pl to the PATH environment variable or specify the path with --eval_exec_dir.");
}
if (qx("$cmdpars{'eval_exec_dir'}evaluate_gtf.pl" 2>&1) =~ /^Can\'t\slocate\s(\w+\.pm)/ ){
    die ("eval is not executable. The perl libary " . $1 . " cannot be located.\n" . 
	 "Please add the directory which contains " . $1 . " to the PERL5LIB environment variable, e.g. add the following line to your .bashrc file:\n\n" . 
	 "export PERL5LIB=\$PERL5LIB:/path/to/" . $1 . "\n\n");
}

my @gfflines = ();

# reading in annotation file and checking if it is in a valid GTF format
print STDERR "reading in annotation file $cmdpars{'anno'}...\n";
open (ANNO, <$cmdpars{"anno"}>) or die ("Could not open $cmdpars{'anno'} for reading: $!");
while(<ANNO>){
    if (/^\s*\#.*/ || /^\s*$/){ # skip comment lines and empty lines
	next;
    }
    my @line = split(/\t/,$_);
    if(@line < 9){
	die ("Not GTF format in the following line:\n$_\n");
    }
    if($line[2] eq "CDS" || $line[2] eq "stop_codon" || $line[2] eq "start_codon"){
	if ($line[8] !~ /transcript_id\s"?[^";]+"?;/){
	    die ("Not GTF format in the following line:\n$_\ntranscript_id not found.\n");
	}
	if ($line[8] !~ /gene_id\s"?[^";]+"?;/){
	    die ("Not GTF format in the following line:\n$_\ngene_id not found.\n");
	}
	push @gfflines, $_;
    }	
}
close(ANNO);

# temporary directory that contains the prediction and the annotation split by seq
my $gffDir = "tempGFF";
system ("rm -rf $gffDir; mkdir $gffDir");
    
my @intervals=(); # hash of genomic intervals
my %seqlist=(); # hash of sequences (only keys, no values)
    
# make gene IDs unique
#my $geneID = 0;
   
open (PRED, <$cmdpars{"pred"}>) or die ("Could not open $cmdpars{'pred'} for reading: $!");
open (JOINPRED, ">$gffDir/pred.gtf") or die("Could not open $gffDir/pred.gtf for writing: $!");

while(<PRED>){
    if(/prediction on sequence range (\w+):(\d+)\-(\d+)/){
	push @intervals,[$1, $2, $3]; # store genomic intervals on which gene prediction is executed
	print JOINPRED $_;
    }
 #   if(/\tgene\t/){
  #      $geneID++;
   # }
    #s/g\d+/g$geneID/g; # replace gene ID with new unique gene ID
    if(/\t(CDS|stop_codon|start_codon)\t/){
	print JOINPRED $_;
    }	
}
close(PRED);
close(JOINPRED);

# join overlapping genomic intervals
# sort intervals by 1. chromosome, 2. start
@intervals = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @intervals;

my @joined=(); # array of joined genomic intervals
my ($chr, $start, $end);
foreach (@intervals){
    if(defined($chr) && $_->[0] eq $chr && $_->[1] - $end <= 0 ) { # overlap between the last and the current interval
	if($end < $_->[2]){
	    $end = $_->[2];
	}
    } else {
	if (defined($chr)){
	    push @joined,[$chr, $start, $end];
	    if (!$seqlist{$chr}){ # add new sequences to seqlist
		$seqlist{$chr} = 1;
	    }
	}
	($chr, $start, $end) = ($_->[0], $_->[1], $_->[2]);
    }
}
push @joined,[$chr, $start, $end];
if (!$seqlist{$chr}){ # add new sequences to seqlist
    $seqlist{$chr} = 1;
}

# make a new annotation file that only contains features from the training set that are completely contained in one of the intervals
# (if necessary, this can be done faster with a single loop over the intervals, requires presorting of @gfflines)
open(ANNO, '>', "$gffDir/anno.gtf") or die ("Could not open $gffDir/anno.gtf for writing: $!");
if($wholegenome){
    foreach my $line (@gfflines){
	print ANNO $line;
    }
}
else{
    foreach my $line (@gfflines){
	my @gffline = split(/\t/,$line);
	my ($chr, $start, $end)=($gffline[0], $gffline[3], $gffline[4]);
	foreach my $i (@joined){
	    if($chr eq $i->[0] && !($start > $i->[2]) && !($end < $i->[1]) ){
		print ANNO $line;
		last;
	    }
	}
    }
}
close(ANNO);

# join genes
if($joingenes){
    system("mv $gffDir/pred.gtf $gffDir/pred.unfiltered.gtf");
    system("$cmdpars{'jg_exec_dir'}joingenes $jg_pars -g $gffDir/pred.unfiltered.gtf -o $gffDir/pred.gtf");
}

# split annotation and prediction file by seqs and prepare
# list files that contain the directories of the GTF files being compared (required by eval)
system ("rm -f $gffDir/annotation_list");
system ("rm -f $gffDir/prediction_list");
foreach my $seq (keys %seqlist) {
    system ("echo '$gffDir/$seq.anno.gtf' >> $gffDir/annotation_list");
    system ("echo '$gffDir/$seq.pred.gtf' >> $gffDir/prediction_list");
    system ("grep \"^$seq\\b\" $gffDir/pred.gtf > $gffDir/$seq.pred.gtf");
    system ("grep \"^$seq\\b\" $gffDir/anno.gtf > $gffDir/$seq.anno.gtf");
}

# call evaluate_gtf
system ("evaluate_gtf.pl $gffDir/annotation_list $gffDir/prediction_list");
