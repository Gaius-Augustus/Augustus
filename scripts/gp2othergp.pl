#!/usr/bin/perl

# convert genepred format from UCSC gtfToGenePred -genePredExt
# table genePredExt
#"A gene prediction with some additional info."
#    (
#     string name;        "Name of gene (usually transcript_id from GTF)"
#     string chrom;       "Chromosome name"
#     char[1] strand;     "+ or - for strand"
#     uint txStart;       "Transcription start position"
#     uint txEnd;         "Transcription end position"
#     uint cdsStart;      "Coding region start"
#     uint cdsEnd;        "Coding region end"
#     uint exonCount;     "Number of exons"
#     uint[exonCount] exonStarts; "Exon start positions"
#     uint[exonCount] exonEnds;   "Exon end positions"
#     uint id;            "Unique identifier"
#     string name2;       "Alternate name (e.g. gene_id from GTF)"
#     string cdsStartStat; "enum('none','unk','incmpl','cmpl')"
#     string cdsEndStat;   "enum('none','unk','incmpl','cmpl')"
#     lstring exonFrames; "Exon frame offsets {0,1,2}"
#    )
# 
# to the genepred format that can be converted to bigGenePred:
#
#   string chrom;       "Reference sequence chromosome or scaffold"
#   uint   chromStart;  "Start position in chromosome"
#   uint   chromEnd;    "End position in chromosome"
#   string name;        "Name or ID of item, ideally both human readable and unique"
#   uint score;         "Score (0-1000)"
#   char[1] strand;     "+ or - for strand"
#   uint thickStart;    "Start of where display should be thick (start codon)"
#   uint thickEnd;      "End of where display should be thick (stop codon)"
#   uint reserved;       "RGB value (use R,G,B string in input file)"
#   int blockCount;     "Number of blocks"
#   int[blockCount] blockSizes; "Comma separated list of block sizes"
#   int[blockCount] chromStarts; "Start positions relative to chromStart"
#   string name2;       "Alternative/human readable name"
#   string cdsStartStat; "enum('none','unk','incmpl','cmpl')"
#   string cdsEndStat;   "enum('none','unk','incmpl','cmpl')"
#   int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
#   string type;        "Transcript type"
#   string geneName;    "Primary identifier for gene"
#   string geneName2;   "Alternative/human readable gene name"
#   string geneType;    "Gene type"

while(<>){
    @f = split;
    $name = shift @f;
    $chrom = shift @f;
    $strand = shift @f;
    $txStart = shift @f;
    $txEnd = shift @f;
    $cdsStart = shift @f;
    $cdsEnd = shift @f;
    $exonCount = shift @f;
    $exonStarts = shift @f;
    $exonEnds = shift @f;
    $id = shift @f;
    $name2 = shift @f;
    $cdsStartStat = shift @f;
    $cdsEndStat = shift @f;
    $exonFrames = shift @f;

    @s = split /,/, $exonStarts;
    @e = split /,/, $exonEnds;
    @len = ();
    @relStarts = ();
    while (@s){
	my $st = shift @s;
	push @relStarts, $st - $txStart;
	push @len, (shift @e) - $st;
    }
    
    print "$chrom\t$txStart\t$txEnd\t$name\t0\t$strand\t$cdsStart\t$cdsEnd\t255,0,0\t$exonCount\t" . join(",", @len), "\t" 
	. join(",", @relStarts) . "\t$name\t$cdsStartStat\t$cdsEndStat\t$exonFrames\tcoding\t?\t?\t?\n";
}
