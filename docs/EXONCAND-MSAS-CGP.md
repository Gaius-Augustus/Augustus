# Comparatively Finding Exon Candidates (conserved ORFs) and their MSA

This funcionality of comparative AUGUSTUS (AUGUSTUS-CGP) can be used to find conserved exon candidates in aligned genomes.
It is a preliminary step for comparative gene prediction, but can be useful on its own right to find candidate conserved coding exons, e.g. to score them with an external tool like aladdin or PhyloCSF.

A candidate exon (**Exon Candidate**, EC) is a region in one sequence, together with frame and strand information that could be a coding exon (CDS). In particular,
 - it does not have a stop codon in the given reading frame and
 - the signals at the boundaries (e.g. splice sites) exceed a weak minimum threshold.

We call an "**Ortho Exon**" (OE) a tuple of candidate exons (ExonCandidate, EC) from multiple species 
 - whose left and right boundaries are aligned with eath other
 - which are of the same type (single, initial, internal, terminal) - they have the same signal pair at the boundaries
 - which share the same phase (reading frame) at both boundaries
 - which are on the same "strand the of the genome alignment". They can be on mixed strands in their respective genomes, though.
 
Example call
```
augustus --/CompPred/printExonCandsMSA=1 \
  --treefile=tree.nwk \
  --alnfile=test2.maf \
  --referenceFile=refanno.gtf --refSpecies=dmel \
  --species=fly \
  --dbaccess=flies.db \
  --speciesfilenames=genomes.tbl \
  --exoncands=1 \
  --/CompPred/outdir=out --outfile=out/train.out --errfile=out/train.err \
  --alternatives-from-evidence=0 --printOEs=1 --/CompPred/omega=off \
```

  -  `--/CompPred/printExonCandsMSA=1` is the specific mode of AUGUSTUS described here.
          With this flag, no gene predictions are made and no comparative parameter are trained.  
  - `--dbaccess` an *sqlite3* database than contains an index for fasta access to the genomes
  - `--referenceFile` is the path to a GFF file with a reference annotation of CDS lines. They are used to classify the OrthoExons as true (`y=1`, identical to a reference CDS) or false (`y=0`), which can be used for training. If this is not required the file can be empty.
  - `--printOEs=1` prints the "OrthoExons" in a GFF file
  - `--exoncands` Prints the exon candidates into a GFF file if set to true
  - `--speciesfilenames` a text file with the paths to the genomes
  - `-treefile=tree.nwk` contains the names of species, the structure is disregarded here
  - `--species=fly` contains a parameter set for predictions, not used here
  
  
  Forall other input formats and parameters see [RUNNING-AUGUSTUS-IN-CGP-MODE.md](RUNNING-AUGUSTUS-IN-CGP-MODE.md).


Example output:
```
-------Speciesnames:--------
species 0       dmoj
species 1       dvir
species 2       dgri
species 3       dmel
species 4       dsim
species 5       dsec
species 6       dere
species 7       dyak
species 8       dana
species 9       dpse
species 10      dper
species 11      dwil
...
y=1     OE2030: CDS     2L      6012179 6012236 -       0
0       tgctgcgtcta--tttgtgg--c---at--tccgcgagtcgagttgacgcgtttgcacattacccat
1       tgctgcgtcta--tttgtcg--c---at--tgcgcgactccagctgacgcgtttgcacattgcccat
2       tacttcgtcta--tttgtcg--c---at--tgcgagactcgaattgacgcgtttgcacattgcccat
3       tgctcctccga--t-----t--cggcacttggcggttctccagttgtcgggcttgggtgttgcccat
4       tgctcctccga--t-----t--cggcacttgggggttctccagcggtcgggcatgggtgttgcccat
5       tgctcctccga--t-----t--cggctcttggcggttcaccagcggtcgggcctgggtgttgcccat
6       tgctcctccga--t-----t--cggcacctggcggttctcaatcggacgggcttgggtgttgcccat
7       tgctcctccga--t-----t--cggcacattgcggttctccagctgccgggcttgggtgttgcccat
8       tacttcttcga--g-----aagcagaac--cacgattctccaatggacgcgtctgtgtatttcccat
9       tactgcgccggcttgcc--a--c---at--tccgagtctcaaatggacgtggctgcacattgcccat
10      tactgcgccggcttgcc--a--c---at--tccgagtctcaaatggacgtggctgcacattgcccat
11      tgctacgtcga--tttccgg--c---at--tgcgagcctccaattgacgcggttgagcattgcccat

y=0     OE2252: CDS     2L      6013809 6013837 +       0
3       atgtggcag-cggcggcggcggcgggcgag
4       atgtggcag-cggcggcggcggcgggcgag
5       atgtggcgg-cggcagcggcggcgggcgag
7       atgcggctg-cggct---gcggcgggcgag
9       atgtgccaccagcct---gcg-------at
10      atgtgccacccgcct---gcg-------at

```

The last two fields in the lines starting with "y=" are the strand in the reference genome and the frame of all exon candidates.

Example file content of `orthoExons.dmel.gff3`:

```
2L  OE1 exon  6000131 6000136 -6.88 + 0 ID=12;Name=12;Note=TERMINAL;n=3;cons=0.417;LeftCons=0.36;rightCons=0.383;....
2L  OE1 exon  6000131 6000137 -6.83 + 1 ID=13;Name=13;Note=INTERNAL0;n=3;cons=0.395;LeftCons=0.36;rightCons=0.381;...

```
