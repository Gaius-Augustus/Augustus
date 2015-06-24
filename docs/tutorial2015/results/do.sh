# these commands will run the whole tutorial automatically
# input in this directory:
# chr2L.sm.fa  rnaseq1.fq  rnaseq2.fq

mkdir gindex
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir gindex --genomeFastaFiles chr2L.sm.fa # ~1m
STAR --runThreadN 4 --genomeDir gindex --readFilesIn rnaseq1.fq rnaseq2.fq --alignIntronMax 100000 \
  --outSAMtype BAM  SortedByCoordinate --outWigType wiggle --outWigStrand Unstranded # ~2m

for s in nasonia zebrafish tomato; do
  augustus --species=$s chr2L.sm.fa --softmasking=on --predictionEnd=1000000 > aug.$s.1-1M.gff &
done

sleep 2m

echo "chr2L   23513712" > chrom.sizes # size of chromosomes
mv Signal.Unique.str1.out.wig rnaseq.wig    # rename for convenience
mv Aligned.sortedByCoord.out.bam rnaseq.bam 
wigToBigWig rnaseq.wig chrom.sizes rnaseq.bw
bamtools index -in rnaseq.bam # index required by UCSC browser
scp rnaseq.bw rnaseq.bam rnaseq.bam.bai mario@hgwdev:~/public_html/tutorial2015/results/ # copy to web space


echo "browser hide all" > customtrack
echo "track name=\"STAR RNA-Seq alignments\" type=bam visibility=4 bigDataUrl=http://hgwdev.cse.ucsc.edu/~mario/tutorial2015/results/rnaseq.bam" >> customtrack
echo "track name=\"STAR coverage\" type=bigWig visibility=2 bigDataUrl=http://hgwdev.cse.ucsc.edu/~mario/tutorial2015/results/rnaseq.bw" >> customtrack

for s in nasonia zebrafish tomato; do
   echo "track name=\"AUGUSTUS $s abinitio\"  db=dm6 visibility=3" >> customtrack
   grep -P "AUGUSTUS\t(CDS|exon)\t" aug.$s.1-1M.gff >> customtrack # use only the relevant coordinate lines
done

# create hints from rnaseq
bam2hints --intronsonly --in=rnaseq.bam --out=hints.intron.gff
cat rnaseq.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep \
  --radius=4.5 --pri=4 --strand="." > hints.ep.gff
cat hints.intron.gff hints.ep.gff > hints.gff

# Predict genes using hints:
for i in {0..3}; do
    augustus --species=nasonia chr2L.sm.fa --softmasking=on --predictionStart=$((i*2000000)) --predictionEnd=$(((i+1)*2000000+50000)) \
    --hintsfile=hints.gff --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg  > aug.nasonia.hints.$i.gff &
done

sleep 7m

cat aug.nasonia.hints.{0..3}.gff | join_aug_pred.pl > aug.nasonia.hints.gff

#  Add the browser track of genes with evidence:

echo "track name=\"AUGUSTUS nasonia hints\"  db=dm6 visibility=3" >> customtrack
grep -P "AUGUSTUS\t(CDS|exon)\t" aug.nasonia.hints.gff  >> customtrack
gzip -c customtrack > customtrack.whints.gz
scp customtrack.whints.gz mario@hgwdev:~/public_html/tutorial2015/results/ # copy to web space

# subset of 100 supported tx
cat aug.nasonia.hints.gff | perl -ne 'if (/\ttranscript\t.*\t(\S+)/){$tx=$1;} if (/transcript supported.*100/) {print "$tx\n";}' | tee supported.lst | wc -l

gff2gbSmallDNA.pl --good=supported.lst aug.nasonia.hints.gff chr2L.sm.fa 5000 genes.gb

##################### training.html

randomSplit.pl genes.gb 100

grep -c LOCUS genes.gb*
new_species.pl --species=bug

etraining --species=bug genes.gb.train

ls -ort $AUGUSTUS_CONFIG_PATH/species/bug/

augustus --species=bug genes.gb.test | tee firsttest.out # takes ~1m

grep -A 22 Evaluation firsttest.out

