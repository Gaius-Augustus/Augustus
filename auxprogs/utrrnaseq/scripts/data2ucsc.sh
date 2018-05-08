#!/bin/bash
# parse GFF files such that they are accecepted by the UCSC Browser
# Ingo Bulla
# 16 May 2013

# replace 9th field of a GFF file by a unique ID
# parameters: the source and the target GFF file
make_id_for_gff() {
  line_num_file=`mktemp`
  num_lines=`wc -l $1 | cut -f1 -d' ' `
  for ((i=1; i<=num_lines; i++)); do
    echo $i >> $line_num_file
  done

  trunc_gff_file=`mktemp`
  cut -f1-8 $1 > $trunc_gff_file

  paste $trunc_gff_file $line_num_file > $2
}


cd ~/workspace/utrrnaseq


# generate GFF files containing UTRs, introns of UTRs, and start & stop codons
# "gene_id" not accepted in ninth column
echo 'track name=utr description="UTRs" visibility=2 color=255,0,0' > human19_utrs.gff
sed 's/gene_id "\([^"]*\)";.*/\1/' utrs.gff | grep -v intron >> human19_utrs.gff
echo 'track name=utr_introns description="UTR introns" visibility=2 color=127,0,0' > human19_utr_introns.gff
sed 's/gene_id "\([^"]*\)";.*/\1/' utrs.gff | grep intron >> human19_utr_introns.gff

echo 'track name=start_stop description="start and stop codons" visibility=2' > human19_start_stop.gff
sed 's/gene_id "\([^"]*\)";.*/\1/' input/human-chr19/start_stop.gff >> human19_start_stop.gff

# generate GFF files containing introns and repeats
make_id_for_gff input/human-chr19/lib1.introns.ff.gff human19_introns.gff
sed -i '1itrack name=intron description="introns" visibility=2 color=0,255,0' human19_introns.gff

make_id_for_gff input/human-chr19/repeats.gff human19_repeats.gff
sed -i '1itrack name=repeats description="repeats" visibility=1 color=0,0,255' human19_repeats.gff

# generate WIG file
tail -n +2 input/human-chr19/trimmed-wig/trimmed.sf.wig > human19_coverage.gff
sed -i '1itrack name=coverage description="RNA-seq coverage" visibility=2 type=wiggle_0' human19_coverage.gff



