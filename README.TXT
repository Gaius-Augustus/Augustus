%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%           Gene Prediction with the Command-Line Version of AUGUSTUS                            %
%                                                                                                %
%           Authors: Lizzy Gerischer, Oliver Keller, Stefanie König, Lars Romoth, Mario Stanke   %
%             Mario Stanke                                                                       %
%             Universitaet Greifswald                                                            %
%             Institut fuer Mathematik und Informatik                                            %
%             17497 Greifswald                                                                   %
%             Fon: +49 3834 864642                                                               %
%             mario.stanke@uni-greifswald.de                                                     %
%                                                                                                %
%             Date: October 23rd, 2015                                                           %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  1.  INTRODUCTION
  2.  INSTALLATION
  3.  RUNNING AUGUSTUS
  4.  SAMPLING: ALTERNATIVE TRANSCRIPTS AND POSTERIOR PROBABILITIES
  5.  USING HINTS (AUGUSTUS+)
  6.  PREDICTIONS USING CDNA
  7.  AUGUSTUS-PPX: PREDICTIONS USING PROTEIN PROFILES
  8.  RETRAINING AUGUSTUS
  9.  WEB-SERVER
  10. MEA: USING THE MAXIMUM EXPECTED ACCURACY APPROACH
  11. CONTACT ME
  12. REFERENCES
  13. LICENCES

              1. INTRODUCTION
              ---------------
AUGUSTUS is a gene prediction program for eukaryotes written by Mario Stanke, Oliver Keller, Stefanie König
and Lizzy Gerischer. It can be used as an ab initio program, which means it bases its prediction purely on the
sequence. AUGUSTUS may also incorporate hints on the gene structure coming from extrinsic sources
such as EST, MS/MS, protein alignments and syntenic genomic alignments. Since version 3.0 AUGUSTUS
can also predict the genes simultaneously in several aligned genomes (see README-cgp.txt).
Currently, AUGUSTUS has been trained for predicting genes in: 

  - Homo sapiens (human), 
  - Drosophila melanogaster (fruit fly), 
  - Arabidopsis thaliana (plant),
  - Brugia malayi (nematode),
  - Aedes aegypti (mosquito),
  - Coprinus cinereus (fungus),
  - Tribolium castaneum (bug)
  - Schistosoma mansoni (worm)
  - Tetrahymena thermophila (ciliate)
  - Galdieria sulphuraria (red algae)
  - Zea mays (maize)
  - Toxoplasma gondii (parasitic protozoa)
  - Caenorhabditis elegans (worm)
  - Aspergillus fumigatus
  - Aspergillus nidulans
  - Aspergillus oryzae
  - Aspergillus terreus
  - Botrytis cinerea
  - Callorhinchus milii
  - Candida albicans
  - Candida guilliermondii
  - Candida tropicalis
  - Chaetomium globosum
  - Coccidioides immitis
  - Cryptococcus neoformans gattii
  - Cryptococcus neoformans neoformans
  - Danio rerio
  - Debaryomyces hansenii
  - Encephalitozoon cuniculi
  - Eremothecium gossypii
  - Fusarium graminearum
  - Gallus gallus
  - Histoplasma capsulatum
  - Kluyveromyces lactis
  - Laccaria bicolor
  - Lodderomyces elongisporus
  - Magnaporthe grisea
  - Neurospora crassa
  - Nicotiana attenuata (coyote tobacco)
  - Petromyzon marinus (sea lamprey)
  - Phanerochaete chrysosporium
  - Pichia stipitis
  - Rhizopus oryzae
  - Saccharomyces cerevisiae
  - Schizosaccharomyces pombe
  - Staphylococcus aureus
  - Ustilago maydis
  - Yarrowia lipolytica
  - Nasonia vitripennis (wasp)
  - Solanum lycopersicum (tomato)
  - Chlamydomonas reinhardtii (green algae)
  - Amphimedon queenslandica (sponge)
  - Acyrthosiphon pisum (pea aphid)
  - Leishmania tarentolae (protozoa, intronless)
  - Trichinella spiralis
  - Theobroma cacao (cacao)
  - Escherichia coli
  - Thermoanaerobacter tengcongensis (a bacterium)
  - Triticum aestivum (wheat)
  - Ancylostoma ceylanicum

The training annotation files of the following species are a courtesy of Jason Stajich (see also http://fungal.genome.duke.edu/): 
Aspergillus fumigatus, Aspergillus nidulans, Aspergillus oryzae, Aspergillus terreus, Botrytis cinerea, Candida albicans
, Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Cryptococcus neoformans gattii
, Cryptococcus neoformans neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum
, Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa
, Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Ustilago maydis
, Yarrowia lipolytica.

The training for lamprey (Petromyzon marinus) was performed by Falk Hildebrand and Shigehiro Kuraku, based on the
genome assembly (PMAR3.0) provided by the Genome Sequencing Center at Washington University School of Medicine (WUGSC)
in St. Louis.

The training for elephant shark (Callorhinchus milii) was performed by Tereza Manousaki and Shigehiro Kuraku, based on the genome assembly
(made up of 1.4x whole genome shotgun reads) available at http://esharkgenome.imcb.a-star.edu.sg/resources.html.

The training for Pneumocystis jirovecii was performed by Marco Pagni, Philippe Hauser et al as described in Hauser PM, Burdet FX, Cissé OH,
Keller L, Taffé P, Sanglard D, Pagni M., Comparative Genomics Suggests that the Fungal Pathogen Pneumocystis Is an Obligate Parasite
Scavenging Amino Acids from Its Host's Lungs. PLoS One. 2010, Dec 20;5(12):e15152. PubMed PMID: 21188143; PubMed Central PMCID: PMC3004796. 

The parameter training for cacao was done by Don Gilbert (gilbertd at indiana.edu).

Parameters for Ancylostoma ceylanicum were trained by Erich Schwarz (Cornell University).

Parameters for camponotus_floridanus (ant) were contributed by Shishir K Gupta.

Parameters for Apis dorsata were contributed by Francisco Camara Ferreira.

              2. INSTALLATION
              ---------------

1. unpack

>tar -xzf augustus-3.2.1.tar.gz

The tar-archive contains one directory 'augustus-3.1' with the following sub-directories:
bin
src
include
config
examples
scripts
auxprogs
docs

2. compile (if not already compiled)

Install dependencies:
	
   - To turn on the optional support of gzip compressed input files you
    1) Edit common.mk and uncomment the line ZIPINPUT = true
    2) Install these dependencies:
       - Boost C++ Libraries: libboost-iostreams-dev (on Ubuntu: sudo apt-get install libboost-iostreams-dev)
       - zlib library for compression methods. Download from http://www.zlib.net/ or install via package manager.
   
If comparative (multi-species, CGP) AUGUSTUS is enabled via uncommenting 
   COMPGENEPRED = true
   and 
   SQLITE = true or MYSQL = true
in common.mk then the following additional dependencies are required:

   - GNU scientific library (for eigen decompositions of matrices, sudo apt-get install libgsl0-dev)
   - libsqlite3-dev or mysql++ (sudo apt-get install libmysql++-dev)
   - Boost C++ Library: libboost-graph-dev (at least version 1.45 from Nov. 2010)

Compiling bam2hints and filterBam requires:

   - The packages bamtools and libbamtools-dev
     (Ubuntu: sudo apt-get install bamtools libbamtools-dev,
      Possibly, you will need to add the repository http://us.archive.ubuntu.com/ubuntu vivid main universe
      to your /etc/apt/sources.lst, first: 
      deb http://us.archive.ubuntu.com/ubuntu vivid main universe)

> make

After compilation has finished, the command bin/augustus should be executable and print a usage message.

3. As a normal user, add the directory of the executables to the PATH environment variable. E.g. issue

PATH=$PATH:~/augustus/bin:~/augustus/scripts

As an administrator, globally install augustus by typing

make install

or by executing similar commands to those in the "install" section of the top-level Makefile.

4. Optional: set environment variable AUGUSTUS_CONFIG_PATH

> export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/

If the environment variable AUGUSTUS_CONFIG_PATH is set, augustus and etraining will look there for the config directory that contains the
configuration and parameter files, e.g. '~/augustus/config'. You may want to add this line to a startup script (like ~/.bashrc).
If this environment variable is not set, then the programs will look in the path ../config relative to the directory in which the executable lies.
As a third alternative, you can specify this directory on the command line when you run augustus:
--AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/


              3. RUNNING AUGUSTUS
              -------------------

AUGUSTUS has 2 mandatory arguments. The query file and the species.
The query file contains the DNA input sequence and must be in uncompressed (multiple) fasta format,
e.g. the file may look like this

>name_of_sequence_1
agtgctgcatgctagctagct
>name_of_sequence_2
gtgctngcatgctagctagctggtgtnntgaaaaatt

Every letter other than a,c,g,t,A,C,G and T is interpreted as an unknown base. Digits and white spaces are ignored. The number of characters per
line is not restricted.

usage:
augustus [parameters] --species=SPECIES queryfilename

or if you want to output to be redirected to a file (see also outfile parameter):
augustus [parameters] --species=SPECIES queryfilename > output.gff

SPECIES is one of the following identifiers
The directory names under config/species constitutes the complete list.

identifier                               | species
-----------------------------------------|----------------------
human                                    | Homo sapiens
fly                                      | Drosophila melanogaster
arabidopsis                              | Arabidopsis thaliana
brugia                                   | Brugia malayi
aedes                                    | Aedes aegypti
tribolium                                | Tribolium castaneum
schistosoma                              | Schistosoma mansoni
tetrahymena                              | Tetrahymena thermophila
galdieria                                | Galdieria sulphuraria
maize                                    | Zea mays
toxoplasma                               | Toxoplasma gondii
caenorhabditis                           | Caenorhabditis elegans
(elegans)                                | Caenorhabditis elegans 
aspergillus_fumigatus                    | Aspergillus fumigatus
aspergillus_nidulans                     | Aspergillus nidulans
(anidulans)                              | Aspergillus nidulans
aspergillus_oryzae                       | Aspergillus oryzae
aspergillus_terreus                      | Aspergillus terreus
botrytis_cinerea                         | Botrytis cinerea
candida_albicans                         | Candida albicans
candida_guilliermondii                   | Candida guilliermondii
candida_tropicalis                       | Candida tropicalis
chaetomium_globosum                      | Chaetomium globosum
coccidioides_immitis                     | Coccidioides immitis
coprinus                                 | Coprinus cinereus
coprinus_cinereus                        | Coprinus cinereus
coyote_tobacco                           | Nicotiana attenuata
cryptococcus_neoformans_gattii           | Cryptococcus neoformans gattii
cryptococcus_neoformans_neoformans_B     | Cryptococcus neoformans neoformans
cryptococcus_neoformans_neoformans_JEC21 | Cryptococcus neoformans neoformans
(cryptococcus)                           | Cryptococcus neoformans
debaryomyces_hansenii                    | Debaryomyces hansenii
encephalitozoon_cuniculi_GB              | Encephalitozoon cuniculi
eremothecium_gossypii                    | Eremothecium gossypii
fusarium_graminearum                     | Fusarium graminearum
(fusarium)                               | Fusarium graminearium
histoplasma_capsulatum                   | Histoplasma capsulatum
(histoplasma)                            | Histoplasma capsulatum
kluyveromyces_lactis                     | Kluyveromyces lactis
laccaria_bicolor                         | Laccaria bicolor
lamprey                                  | Petromyzon marinus
leishmania_tarentolae                    | Leishmania tarentolae
lodderomyces_elongisporus                | Lodderomyces elongisporus
magnaporthe_grisea                       | Magnaporthe grisea
neurospora_crassa                        | Neurospora crassa
(neurospora)                             | Neurospora crassa
phanerochaete_chrysosporium              | Phanerochaete chrysosporium
(pchrysosporium)                         | Phanerochaete chrysosporium
pichia_stipitis                          | Pichia stipitis
rhizopus_oryzae                          | Rhizopus oryzae
saccharomyces_cerevisiae_S288C           | Saccharomyces cerevisiae
saccharomyces_cerevisiae_rm11-1a_1       | Saccharomyces cerevisiae
(saccharomyces)                          | Saccharomyces cerevisiae
schizosaccharomyces_pombe                | Schizosaccharomyces pombe
thermoanaerobacter_tengcongensis         | Thermoanaerobacter tengcongensis
trichinella                              | Trichinella spiralis
ustilago_maydis                          | Ustilago maydis
(ustilago)                               | Ustilago maydis
yarrowia_lipolytica                      | Yarrowia lipolytica
nasonia                                  | Nasonia vitripennis
tomato                                   | Solanum lycopersicum
chlamydomonas                            | Chlamydomonas reinhardtii
amphimedon                               | Amphimedon queenslandica
pneumocystis                             | Pneumocystis jirovecii
wheat                                    | Triticum aestivum
chicken                                  | Gallus gallus
zebrafish                                | Danio rerio
E_coli_K12                               | Escherichia coli
s_aureus                                 | Staphylococcus aureus


The identifiers in parentheses denote older versions for that species.

'queryfilename' is the filename (including relative path) to the file containing the query sequence(s) in fasta format.

important parameters:

--strand=both, --strand=forward or --strand=backward 
 report predicted genes on both strands, just the forward or just the backward strand. 
 default is 'both'

--genemodel=partial, --genemodel=intronless, --genemodel=complete, --genemodel=atleastone or --genemodel=exactlyone
  partial      : allow prediction of incomplete genes at the sequence boundaries (default)
  intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes
  complete     : only predict complete genes
  atleastone   : predict at least one complete gene
  exactlyone   : predict exactly one complete gene

--singlestrand=true
 predict genes independently on each strand, allow overlapping genes on opposite strands
 This option is turned off by default.

--hintsfile=hintsfilename
  When this option is used the prediction considering hints (extrinsic information) is turned on.
  hintsfilename contains the hints in gff format.

--extrinsicCfgFile=cfgfilename
  Optional. This file contains the list of used sources for the hints and their boni and mali.
  If not specified the file "extrinsic.cfg" in the config directory $AUGUSTUS_CONFIG_PATH is used.

--maxDNAPieceSize=n
  This value specifies the maximal length of the pieces that the sequence is cut into for the core algorithm (Viterbi) to be run. Default is --maxDNAPieceSize=200000.
  AUGUSTUS tries to place the boundaries of these pieces in the intergenic region, which 
  is inferred by a preliminary prediction. GC-content dependent parameters are chosen for each piece of DNA if /Constant/decomp_num_steps > 1 for that species.
  This is why this value should not be set very large, even if you have plenty of memory. 

--protein=on/off
--introns=on/off
--start=on/off
--stop=on/off
--cds=on/off
--codingseq=on/off
  Output options. Output predicted protein sequence, introns, start
  codons, stop codons. Or use 'cds' in addition to 'initial', 'internal',
  'terminal' and 'single' exon. The CDS excludes the stop codon (unless stopCodonExcludedFromCDS=false)
  whereas the terminal and single exon include the stop codon.
--AUGUSTUS_CONFIG_PATH=path
  path to config directory (if not specified as environment variable)
--alternatives-from-evidence=true/false
  report alternative transcripts when they are suggested by hints
--alternatives-from-sampling=true/false
  report alternative transcripts generated through probabilistic sampling
--sample=n
--minexonintronprob=p
--minmeanexonintronprob=p
--maxtracks=n
  For a description of these parameters see section 4 below.

--proteinprofile=filename
  Read a protein profile from file filename. See section 7 below.

--predictionStart=A, --predictionEnd=B
  A and B define the range of the sequence for which predictions should be found. 
  Quicker if you need predictions only for a small part.
--gff3=on/off
  output in gff3 format
--UTR=on/off
  predict the untranslated regions in addition to the coding sequence. This currently works only for
  human, galdieria, toxoplasma and caenorhabditis.

--outfile=filename
  print output to filename instead to standard output. 
  This is useful for computing environments, e.g. parasol jobs, which do not allow shell redirection.

--noInFrameStop=true/false
  Don't report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. Default: false

--noprediction=true/false
  If true and input is in genbank format, no prediction is made. Useful for getting the annotated protein sequences.

--contentmodels=on/off
  If 'off' the content models are disabled (all emissions uniformly 1/4). The content models are; coding region Markov chain
  (emiprobs), initial k-mers in coding region (Pls), intron and intergenic regin Markov chain. This option is intended for
  special applications that require judging gene structures from the signal models only, e.g. for predicting the effect of
  SNPs or mutations on splicing. For all typical gene predictions, this should be true. Default: on


For example you may type

>augustus --species=human --UTR=on ../examples/example.fa


The output format is gtf similar to General Feature Format (gff), see
http://www.sanger.ac.uk/Software/formats/GFF/. 
It contains one line per predicted exon. Example:

HS04636   AUGUSTUS   initial    966     1017    .       +       0       transcript_id "g1.1"; gene_id "g1";
HS04636   AUGUSTUS   internal   1818    1934    .       +       2       transcript_id "g1.1"; gene_id "g1";

The columns (fields) contain: 
seqname   source     feature    start   end   score   strand   frame    transcript and gene name


AUGUSTUS also accepts files in annotated GENBANK format as input. This is needed for training.
Also when predicting on a genbank file AUGUSTUS compares its prediction to the annotation and prints out a
statistic. Example genbank file format accepted by AUGUSTUS: 

LOCUS       HS04636   9453 bp  DNA
FEATURES             Location/Qualifiers
     source          1..9453
     CDS             join(966..1017,1818..1934,2055..2198,2852..2995,3426..3607,
                     4340..4423,4543..4789,5072..5358,5860..6007,6494..6903)
BASE COUNT     2937 a   1716 c  1710 g   3090 t
ORIGIN
        1 gagctcacat taactattta cagggtaact gcttaggacc agtattatga ggagaattta
       61 cctttcccgc ctctctttcc aagaaacaag gagggggtga aggtacggag aacagtattt
      121 cttctgttga aagcaactta gctacaaaga taaattacag ctatgtacac tgaaggtagc
      ...
     9421 aaaaaaaaaa aaaaatcgat gtcgactcga gtc	
//

Another example that is important for training the UTR models. The following genbank file will be interpreted 
as having three genes. One gene ('A') with both 5' and 3' UTR and two single UTRs without matching coding sequence.
Gene 'B' consists just of the 5'UTR, gene 'C' just of the 3' UTR.

LOCUS       example2   9453 bp  DNA
FEATURES             Location/Qualifiers
     source          1..9453
     mRNA            join(100..200,900..1017,1818..2000,2100..2200)
                     /gene="A"
     CDS             join(966..1017,1818..1934)
                     /gene="A"
     mRNA            join(3100..3200,3500..>3600)
                     /gene="B"
     mRNA            join(<4100..4200,4500..4600)
                     /gene="C"

BASE COUNT     2937 a   1716 c  1710 g   3090 t
ORIGIN
        1 gagctcacat taactattta cagggtaact gcttaggacc agtattatga ggagaattta
       61 cctttcccgc ctctctttcc aagaaacaag gagggggtga aggtacggag aacagtattt
      121 cttctgttga aagcaactta gctacaaaga taaattacag ctatgtacac tgaaggtagc
      ...
     9421 aaaaaaaaaa aaaaatcgat gtcgactcga gtc	
//

              4. SAMPLING: ALTERNATIVE TRANSCRIPTS AND POSTERIOR PROBABILITIES
              ----------------------------------------------------------------

Note that for the prediction of alternative splicing there is another method described in 5. below.

Alternative transcripts (from sampling)
---------------------------------------
When you say on the command line
--alternatives-from-sampling=true
or edit the appropriate line in the configuration file for your species to
alternatives true then AUGUSTUS may report multiple transcripts per gene. A gene is then defined as a set of transcripts, whose 
coding sequences (indirectly) overlap. The number of alternatives AUGUSTUS reports for a gene depends on which 
ones are likely alternatives. If just one transcript is likely in that region then also just one transcript 
is reported. The behavior of AUGUSTUS can be adjusted with the parameters
--minexonintronprob=p
--minmeanexonintronprob=p
--maxtracks=n (default -1, no limit)

The posterior probability of every exon and every intron in a transcript must be at least 'minexonintronprob', otherwise the transcript is not reported.
minexonintronprob=0.1 is a reasonable value. In addition the geometric mean of the probabilities of all exons and introns 
must be at least 'minmeanexonintronprob'. minmeanexonintronprob=0.4 is a reasonable value. 
The maximum number of tracks when displayed in a genome browser is 'maxtracks' (unless maxtracks=-1, then it is unbounded).
In cases where all transcripts of a gene overlap at some position this is also the maximal number of transcripts for that gene.
I recommend increasing the parameter 'maxtracks' for improving sensitivity and setting 'maxtracks' to 1 and increasing 
minmeanexonintronprob and/or minexonintronprob in order to improve specificity.

Posterior probabilities
-----------------------
AUGUSTUS reports the posterior probabilities of exons, introns, transcripts and genes. The posterior probability
of an exon is the conditional probability that the random gene structure has some exon with these coordinates on this strand 
given the input sequence. It not only depends on the sequence in the range of the exon itself like an exon score but 
is influenced for example by the possibilities of compatible neighboring exons. The intron score is similar.
The reported probability of a transcript is the probability that a splice variant is exactly like in the given transcript.
The reported probability of a gene is the probability that SOME coding sequence is in the reported range on the reported 
strand, regardless of the exact transcript.

The posterior probabilities are estimated using a sampling algorithm. The parameter --sample==n adjusts the number 
of sampling iterations. The higher 'n' is the more accurate is the estimation but it usually isn't important 
that the posterior probability is very accurate. Every 30 sample iterations take about the same time as one run without 
sampling, e.g. --sample=60 takes about 3 times as much time as --sample=0 (which was standard up to version 1.6).
The default is
--sample=100
If you do not need the posterior probabilities or alternative transcripts, say
--sample=0

There are 3 common scenarios for above parameters, depending on what you want:

Just output the most likely gene structure as in previous versions. No posterior probabilities, no alternatives:
--sample=0 --alternatives=false

Output the most likely gene structure and report posterior probabilities:
--sample=100 --alternatives=false

Output alternative transcripts and report posterior probabilities:
--sample=100 --alternatives=true

Be aware that sampling is pseudorandom and that the results may vary from machine to machine.

Heating
-------
The probabilistic model of AUGUSTUS can be seen as a rough approximation to reality. A consequence is that the posterior probabilities for the 
strong exons (e.g. the ones called by the Viterbi algorithm) tend to be larger than the actually measured precision (specificity) values. For
example, in human, only 94.5% of the exons with predicted posterior probability >= 98% (under the default --sample=100) are actually correct.
See docs/CDS.sp.{pdf,README} for more data and an explanation. If the aim of the sampling is to produce a diverse, sensitive (including) set of gene
structures, you can use this parameter

--temperature=t

where t is one of 0,1,2,3,4,5,6,7. All probabilities of the model are then taken to the power of (8-t)/8, i.e. t=0 (the default) does nothing. The larger t
the more alternatives are sampled. t=3 is a good compromise between getting a high sensitivity but not getting too many exons sampled in total.
For t=3, 96.1% of human exons with posterior probability >= 98% are correct.

              5. USING HINTS (AUGUSTUS+)
              --------------------------
AUGUSTUS can take hints on the gene structure. It currently accepts 16 types of hints:  
start, stop, tss, tts, ass, dss, exonpart, exon, intronpart, intron, CDSpart, CDS, UTRpart, UTR, irpart, nonexonpart.
The hints must be stored in a file in gff format containing one hint per line. 
Example of a hintsfile:

HS04636	mario	exonpart	500	506	.	-	.	source=M
HS04636	mario	exon	966	1017	.	+	0	source=P
HS04636	AGRIPPA	start	966	968	6.3e-239	+	0	group=gb|AAA35803.1;source=P
HS04636	AGRIPPA	dss	2199	2199	1.3e-216	+	.	group=gb|AAA35803.1;source=P
HS04636	mario	stop	7631	7633	.	+	0	source=M
HS08198	AGRIPPA	intron	2000	2000	0	+	.	group=ref|NP_000597.1;source=E
HS08198	AGRIPPA	ass	757	757	1.4e-52	+	.	group=ref|NP_000597.1;source=E

The fields must be separated by a tabulator. In the first column (field) the sequence name is given. 
In this case the hints are together about two sequences. The second field is the name of the program that produced the hint. 
It is ignored here. The third column specifies the type of the hint. The 4th and 5th column specify the begin and end position of the hint.
Positions start at 1. The 6th colum gives a score. The 7th the strand. The 8th the reading frame as defined in the GFF standard. 
The 9th column contains arbitrary extra information but it must contain a string 'source=X' where X
is the source identifier of the hint. Which values for X are possible is specified in the file augustus/config/extrinsic.cfg,
e.g. X=M, E, or P.

AUGUSTUS can follow a hint, i.e. predict a gene structure that is compatible with it, or 
AUGUSTUS can ignore a hint, i.e. predict a gene structure that is not compatible with it. The probability that AUGUSTUS ignores a 
hint is the smaller the more reliable the hints of this type are. 

Run AUGUSTUS using the --hintsfile option. Example:

augustus --species=human --hintsfile=../examples/hints.gff --extrinsicCfgFile=../config/extrinsic/extrinsic.MPE.cfg ../examples/example.fa

As an alternative to giving the option --extrinsicCfgFile you can replace augustus/config/extrinsic.cfg with the appropriate file, as this 
file is read by default when the option --extrinsicCfgFile is not given.

The preferable way to use repeat information is via softmasking in which the bases in repeat regions are lower case
(a,c,g,t instead of A,C,G,T) in the input. Running augustus like this

augustus --species=human --softmasking=1

will interpret masked regions as evidence against exons (nonexonparts hints with a default bonus of 1.15). This is slightly more accurate
than hard masking (with N's), which looses information. On human augustus is also more than twice as fast with softmasking=1 than on hard
masked sequence.

Explanation of the file format of the extrinsic.cfg file.

The gff/gtf file containint the hints must contain somewhere in the last
column an entry source=?, where ? is one of the source characters listed in
the line after [SOURCES] above. You can use different sources when you have
hints of different reliability of the same type, e.g. exon hints from ESTs
and exon hints from evolutionary conservation information.

In the [GENERAL] section the entries second column specify a bonus for obeying
a hint and the entry in the third column specify a malus (penalty) for
predicting a feature that is not supported by any hint. The bonus and the
malus is a factor that is multiplied to the posterior probability of gene
structures. 
Example: 
  CDS     1000  0.7  ....
means that, when AUGUSTUS is searching for the most likely gene structure,
every gene structure that has a CDS exactly as given in a hint gets
a bonus factor of 1000. Also, for every CDS that is not supported the
probability of the gene structure gets a malus of 0.7. Increase the bonus to
make AUGUSTUS obey more hints, decrease the malus to make AUGUSTUS predict few
features that are not supported by hints. The malus helps increasing
specificity, e.g. when the exons predicted by AUGUSTUS are suspicious because
there is no evidence from ESTs, mRNAs, protein databases, sequence
conservation, transMapped expressed sequences.
Setting the malus to 1.0 disables those penalties. Setting the bonus to 1.0
disables the boni. 

      start: translation start (start codon), specifies an interval that contains
             the start codon. The interval can be larger than 3bp, in which case
             every ATG in the interval gets a bonus. The highest bonus is given
             to ATGs in the middle of the interval, the bonus fades off towards the ends.
       stop: translation end  (stop codon), see 'start'
        tss: transcription start site, see 'start'
        tts: transcription termination site, see 'start'
        ass: acceptor (3') splice site, the last intron position, for only approximately known ass an interval can be specified
        dss: donor (5') splice site, the first intron position, for only approximately known dss an interval can be specified
   exonpart: part of an exon in the biological sense. The bonus applies only
             to exons that contain the interval from the hint. Just
             overlapping means no bonus at all. The malus applies to every
             base of an exon. Therefore the malus for an exon is exponential
             in the length of an exon: malus=exonpartmalus^length.
	     Therefore the malus should be close to 1, e.g. 0.99.
       exon: exon in the biological sense. Only exons that exactly match the
             hint get a bonus. Exception: The exons that contain the start
             codon and stop codon. This malus applies to a complete exon
             independent of its length.
 intronpart: introns both between coding and non-coding exons. The bonus
             applies to every intronic base in the interval of the hint.
     intron: An intron gets the bonus if and only if it is exactly as in the hint.
    CDSpart: part of the coding part of an exon. (CDS = coding sequence)
        CDS: coding part of an exon with exact boundaries. For internal exons
             of a multi exon gene this is identical to the biological
             boundaries of the exon. For the first and the last coding exon
             the boundaries are the boundaries of the coding sequence (start, stop).
        UTR: exact boundaries of a UTR exon or the untranslated part of a
             partially coding exon.
    UTRpart: The hint interval must be included in the UTR part of an exon.
     irpart: The bonus applies to every base of the intergenic region. If UTR
             prediction is turned on (--UTR=on) then UTR is considered
             genic. If you choose against the usual meaning the bonus of
             irparts to be much smaller than 1 in the configuration file you
             can force AUGUSTUS to not predict an intergenic region in the
             specified interval. This is useful if you want to tell AUGUSTUS
             that two distant exons belong to the same gene, when AUGUSTUS
             tends to split that gene into smaller genes.
nonexonpart: intergenic region or intron. The bonus applies to very non-exon
             base that overlaps with the interval from the hint. It is
             geometric in the length of that overlap, so choose it close to
             1.0. This is useful as a weak kind of masking, e.g. when it is
             unlikely that a retroposed gene contains a coding region but you
             do not want to completely forbid exons.
  genicpart: everything that is not intergenic region, i.e. intron or exon or UTR if
             applicable. The bonus applies to every genic base that overlaps with the
             interval from the hint. This can be used in particular to make Augustus
             predict one gene between positions a and b if a and b are experimentally
             confirmed to be part of the same gene, e.g. through ESTs from the same clone.
             alias: nonirpart

Any hints of types dss, intron, exon, CDS, UTR that (implicitly) suggest a donor splice
site allow AUGUSTUS to predict a donor splice site that has a GC instead of the much more common GT.
AUGUSTUS does not predict a GC donor splice site unless there is a hint for one.

Starting in column number 4 you can tell AUGUSTUS how to modify the bonus 
depending on the source of the hint and the score of the hint. 
The score of the hints is specified in the 6th column of the hint gff/gtf.
If the score is used at all, the score is not used directly through some
conversion formula but by distinguishing different classes of scores, e.g. low
score, medium score, high score. The format is the following:
First, you specify the source character, then the number of classes (say n), then you
specify the score boundaries that separate the classes (n-1 thresholds) and then you specify
for each score class the multiplicative modifier to the bonus (n factors). 

Examples:

M 1 1e+100
means for the manual hint there is only one score class, the bonus for this
type of hint is multiplied by 10^100. This practically forces AUGUSTUS to obey
all manual hints.

T    2       1.5 1 5e29
For the transMap hints distinguish 2 classes. Those with a score below 1.5 and
with a score above 1.5. The bonus if the lower score hints is unchanged and
the bonus of the higher score hints is multiplied by 5x10^29.

D    8     1.5  2.5  3.5  4.5  5.5  6.5  7.5  0.58  0.4  0.2  2.9  0.87  0.44 0.31  7.3
Use 8 score classes for the DIALIGN hints. DIALIGN hints give a score, a strand and
reading frame information for CDSpart hints. The strand and reading frame are often correct but not
often enough to rely on them. To account for that I generated hints for all
6 combinations of a strand and reading frame and then used 2x2x2=8 different
score classes:
{low score, high score} x {DIALIGN strand, opposite strand} x {DIALIGN reading frame, other reading frame}
This example shows that scores don't have to be monotonous. A higher score
does not have to mean a higher bonus. They are merely a way of classifying the
hints into categories as you wish. In particular, you could get the effect of
having different sources by having just hints of one source and then distinguishing
more scores classes.

Alternative Transcripts / Alternative Splicing (evidence based)
---------------------------------------------------------------

AUGUSTUS can predict alternative splicing or - more general - alternative transcripts that are suggested by evidence given in hints.
The method is very general. But to give an example: If two EST alignments to the same genomic area cannot be explained by a single
transcript then AUGUSTUS can predict a gene with two different splice forms, one splice form
compatible with each of the EST alignments.

Grouping hints:
Each hint can be given a group name, by specifying 'group=goupname;' or 'grp=goupname;' in the last column for the hint
in the gff file. This should be used to group all the hints coming from the alignment of the same
sequence to the genome. For example, if an EST with the name est_xyz aligns to the genome with one gap
suggesting an intron then the hints resulting from that alignment could look like this

HS04636	blat2hints	exonpart	500	599	.	+	.	group=est_xyz; source=E
HS04636	blat2hints	intron		600	700	.	+	.	group=est_xyz; source=E
HS04636	blat2hints	exonpart	701	900	.	+	.	group=est_xyz; source=E

Grouping tells AUGUSTUS that hints belong together. Ideally, all hints of a group are obeyed by a predicted
transcript or the whole group of hints is ignored when making the prediction.

Prioritizing groups:
Hints or hint groups can be given a priority by specifying 'priority=n;' or 'pri=n' in the last
column for the hint in the gff file. For example

HS04636	blat2hints	exonpart	500	599	.	+	.	priority=2; source=E
HS04636	blat2hints	intron		550	650	.	+	.	priority=5; source=mRNA

When two hints or hint groups contradict each other then the hints with the lower priority number
are ignored. This is especially usefull if for a genome several sources of hints are available,
where one source should be trusted when in doubt. For example, the rhesus macaque currently has few native ESTs
but human ESTs often also align to rhesus. Giving the hints from native ESTs a higher priority means
that AUGUSTUS uses only them for genes with support by native ESTs and uses the alien EST alignments
when native ESTs alignments are not available for a gene. When the priority is not specified, it is
internally set to -1.

When AUGUSTUS is run with --alternatives-from-evidence=false then all hints are given to AUGUSTUS at
the same time no whether they can be explained with a single transcript per gene. AUGUSTUS will then
choose the most likely transcript variant.

When AUGUSTUS is run with --alternatives-from-evidence=true then AUGUSTUS will predict alternative
transcripts based on the alternatives the hints suggest. This can be any form of alternative
splicing, including nested genes contained in introns of other genes, overlapping genes, alternative
translation starts and variation in UTR.


              6. PREDICTIONS USING CDNA
              -------------------------

Improving the predictions by integrating ESTs or mRNA data is fairly simple. Let cdna.fa be a fasta file with ESTs and/or mRNAs.
Here is the list of commands that will do the trick:

blat -minIdentity=92 genome.fa cdna.fa cdna.psl
blat2hints.pl --in=cdna.psl --out=hints.E.gff
augustus --species=human --hintsfile=hints.E.gff --extrinsicCfgFile=extrinsic.ME.cfg genome.fa


Explanation and possible improvements:
BLAT is a fast spliced alignment program from Jim Kent. blat2hints.pl is a script from the AUGUSTUS scripts directory. The file extrinsic.ME.cfg 
states the parameters for the inclusion of the hints. You can manually adjust the few parameters for your genome. I recommend adjusting the 
bonuses and maluses in extrinsic.ME.cfg after a visual inspection of the predictions. For example, if it looks as if AUGUSTUS is trying to 
fit too many spurious EST alignments then reduce the bonuses. 
From experience, some ESTs often align to very many places in the genome. Most of those matches do not correspond to real protein coding gene structures.
It is therefore better to add another step after the BLAT run. The command

pslCDnaFilter -maxAligns=1 cdna.psl cdna.f.psl

will filter the cDNA alignments and report only the highest-scoring spliced alignment(s) for each cDNA. Then use the filtered file cdna.f.psl to 
create hints. The program pslCDnaFilter is part of the Kent source tree (but not in the BLAT distribution).

For RNA-Seq integration please refer to the documentation in doc/readme.rnaseq.html.


              7. AUGUSTUS-PPX: PREDICTIONS USING PROTEIN PROFILES
              ---------------------------------------------------

AUGUSTUS can make its prediction based on a protein profile that can be generated from a 
Multiple Sequence Alignment. The protein profile is passed to AUGUSTUS by specifying the 
parameter --proteinprofile as in the following example:
> msa2prfl.pl fam.aln > fam.prfl 
> augustus --proteinprofile fam.prfl genome.fa

The profile consists of a set of position-specific frequency matrices that model conserved 
regions in a MSA, without deletions or insertions. When equipped with a profile, AUGUSTUS will 
make extra effort to predict genes that are similar to the profile, for example members of
a specific protein family of interest. Prediction accuracy for these genes is generally enhanced
by the extra information from the protein model, while other genes are predicted identically 
to the ab-initio version. 

Creating protein profiles from Multiple Sequence Alignments
-----------------------------------------------------------
The script msa2prfl.pl converts a Multiple Sequence Alignment given in FASTA or CLUSTAL format 
to a protein profile, by computing frequencies from all blocks of at least 6 gapless columns in
the alignment. Minimal block width can be changed with the parameter --width. The script
blocks2prfl.pl converts a flat file from the BLOCKS database to a protein profile

Preparing Core Alignments
-------------------------
Large alignments may not be represented by a block profile, if they do not have sufficiently
many gapless columns. It is then recommended to cluster the sequences according to similarity,
or to discard sequences from the alignment that do not cover most of the blocks. The program
prepareAlign can do that with an MSA in FASTA format. Usage:
prepareAlign < input.fa > output.fa

The environment variables PA_FULL_COL_WEIGHT, PA_SKIP_COL_WEIGHT, PA_MINSIZE, PA_MIN_COL_COUNT
control the behaviour of the program. See the source file for details.

Format of the protein profile inputfile
---------------------------------------
A section "[name]", followed by the name of the family.
Alternating sections "[dist]", and "[block]"; each "[dist]" section contains a line with minimum
and maximum distance between blocks. <max> can be specified as "*" to indicate unbounded distance.

[dist]
<min> <max>

Each "[block]" section contains a frequency matrix, one line in the section corresponding to a column 
in the alignment. Each line contains 21 tab-separated values, the first is the column index in the 
block (0,1,2,...), the other 20 values are the frequencies (adding up to 1), given in the order 
G,D,E,R,K,N,Q,S,T,A,V,L,I,F,Y,W,H,M,C,P 

Example protein profiles are in the directory examples/profile/

Running AUGUSTUS-PPX
--------------------
The running time of AUGUSTUS-PPX is proportional to the size of the profile; as a rule of thumb, the
factor compared to AUGUSTUS is approximately the number of blocks in the profile. For large profiles,
it is recommended to restrict the prediction with --predictionStart and --predictionEnd. On standard 
intel machines, running times of about one hour for a large profile on a region of 1 Mbps size, were 
observed. To determine regions where a profile is relevant, run a fastBlockSearch (see below). 
Important parameters for running AUGUSTUS-PPX are:

/ProteinModel/allow_truncated:                   Enables having profile hits in right-truncated genes 
                                                                                       (default: yes)

/ProteinModel/block_threshold_spec               These parameters control the specificity and 
/ProteinModel/block_threshold_sens               sensitivity when determining block hits. 
                                                                             (default: ..._spec=4.0,
                                                                                       ..._sens=0.4)

Increasing ..._sens and decreasing ..._spec will both result in more block hits found (and possibly 
more genes with profile hits), at the expense of more false positive hits. When the requirements cannot 
be both met, a block is discarded from the profile used for the prediction. Specificity and Sensitivity
are given in units of standard deviation from the expected block score (Percentages can be calculated
by applying the Gaussian distribution function; e.g., the default value of 2.5 corresponds to an estimated
specificity of 99.3%: 7 FP hits in 1000bps). Note that filtering block hits is mainly a performance issue, 
and it is very unlikely that a false positive block hit affects the prediction if its score is low.
To prevent blocks from being discarded from the profile, decrease either of the parameters. 

/ProteinModel/blockpart_threshold_spec          specificity for block prefixes or suffixes (4.5)
/ProteinModel/blockpart_threshold_sens          sensitivity for block prefixes or suffixes (2.0)

The same for the case that a block is disconnected by an intron.

/ProteinModel/weight                            influence of the protein model to the combined score,
                                                can be weighted (default: 1, equal contribution)
                                                A higher value will result in more gene structures closer
 						to the protein model if there are.

Output of AUGUSTUS-PPX
----------------------
If a gene is a profile hit, the following lines are added to the gff output:
- a protein_match feature for every block mapped to the DNA (or block part, if the block was disconnected
  by an intron). If --gff3=on is specified, then target block and protein location are given in the attributes
  column:
chr1    AUGUSTUS        protein_match   161494506       161494595       7.54    +       0       ID=pp.g2.t1.PF00012.13_A;Target=PF00012.13_A 1 30;Target_start=19;

- a interblock_region feature for every exon part in the gene that is not mapped to a block
chr1    AUGUSTUS        interblock_region       161494449       161494505       .       +       0       ID=pp.g2.t1.iBR0

- each translated protein sequence of a block match, as a comment (if --protein=on is specified)
# sequence of block PF00012.13_A   19 [VGVFQQGRVEILANDQGNRTTPSYVAFTDT] 49

Fast block search to determine regions for the gene prediction
--------------------------------------------------------------
If a protein profile and a genome are given, a preliminary search can be performed with the program
fastBlockSearch. It will output the locations of profile hits. An AUGUSTUS-PPX run can then be restricted
to regions containing these locations. It should be run with the same parameters as the AUGUSTUS-PPX run. 
In addition, a threshold can be specified with the parameter --cutoff that controls the number of 
profile hits shown. 
> fastBlockSearch --cutoff=0.5 genome.fa fam.prfl

The profile hits found by the fastBlockSearch may not contain always all blocks. In this
case, it may improve the prediction to modify the profile with the following command
> del_from_prfl.pl fam.prfl 2,3,5
where 2,3,5 is to be replaced with the list of blocks to be deleted from the profile.
              

              8. RETRAINING AUGUSTUS
              ----------------------

Please see the file README.autoAug for documentation for the automatic training script autoAug.pl.
See also the file retraining.html. Here is some background:

AUGUSTUS uses parameters which are species specific like the Markov chain transition probability of coding and non-coding regions.
These parameters can be trained on training sets of annotated genes in genbank format. They are stored in the config directory in 3 
files containing the parameters for the exon-related, intron-related and intergenic-region-related parameters, 
e.g. human_exon_probs.pbl, human_intron_probs.pbl, human_igenic_probs.pbl. For each species there are also parameters like the order 
of the markov chain or the size of the window used for the splice site models. Let's call these meta parameters. 
These meta parameters are stored in a separate file, e.g. human_parameters.cfg. Which meta parameters work best depends on 
the species and on the training set, in particular on the size of the training set. Using the meta parameters of another species 
or for another training set is likely to result in poor prediction performance. The meta parameters are not documented sufficiently.
However, when optimizing the meta parameters for a new species it helps to know their meaning. Please contact me in case you want me 
to do the training. 
The program 'etraining' reads the meta parameters from the .cfg file and a genbank file with annotated genes and writes the 
other species specific parameters into the 3 .pbl files.

usage:
etraining --species=SPECIES trainfilename

'trainfilename' is the filename (including relative path) to the file in genbank format containing the training sequences. 
These can be multi-gene sequences and genes on the reverse strand. However, the genes must not overlap.

              9. WEB-SERVER
              -------------

AUGUSTUS can also be run through a web-interface avaible on the AUGUSTUS home page: http://augustus.gobics.de.

	     10. MEA: USING THE MAXIMUM EXPECTED ACCURACY APPROACH
 	     -----------------------------------------------------

MEA is an alternative decoding approach and is described in the Poster found in docs/MEAposter.pdf.
For predictions using MEA the corresponding AUGUSTUS parameter has to be turned on: --mea=1

              11. CONTACT
              -----------

In case you find a bug, or miss a desireable feature or need executables for a different platform, need
detailed info about the model, please check this forum and post your question at:
http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Forum.Forum

In case you need to train AUGUSTUS on a different organism consider this training web service:
bioinf.uni-greifswald.de/trainaugustus


              12. REFERENCES
              --------------

Stefanie König, Lars Romoth, Lizzy Gerischer, and Mario Stanke (2015)
Simultaneous gene finding in multiple genomes. PeerJ PrePrints, e1594

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008)
"Using native and syntenically mapped cDNA alignments to improve de novo gene finding"
Bioinformatics, doi: 10.1093/bioinformatics/btn013

Mario Stanke, Ana Tzvetkova, Burkhard Morgenstern (2006)
"AUGUSTUS at EGASP: using EST, protein and genomic alignments for improved gene prediction in the human genome"
BMC Genome Biology, 7(Suppl 1):S11

Mario Stanke , Oliver Keller, Irfan Gunduz, Alec Hayes, Stephan Waack, Burkhard Morgenstern (2006)
"AUGUSTUS: ab initio prediction of alternative transcripts"
Nucleic Acids Research, 34: W435-W439. 

Mario Stanke, Oliver Schoeffmann, Burkhard Morgenstern and Stephan Waack 
"Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources",
BMC Bioinformatics, 7:62 (2006)

Mario Stanke and Burkhard Morgenstern (2005)
"AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints",
Nucleic Acids Research, 33, W465-W467

Mario Stanke, Rasmus Steinkamp, Stephan Waack and Burkhard Morgenstern, 
"AUGUSTUS: a web server for gene finding in eukaryotes" (2004), 
Nucleic Acids Research, Vol. 32, W309-W312

Mario Stanke (2003),
Gene Prediction with a Hidden-Markov Model. Ph.D. thesis, Universitaet Goettingen,
http://webdoc.sub.gwdg.de/diss/2004/stanke/

Mario Stanke and Stephan Waack (2003),
Gene Prediction with a Hidden-Markov Model and a new Intron Submodel. Bioinformatics, Vol. 19, Suppl. 2, pages ii215-ii225 


              13. LICENCES
              ------------

All source code, i.e.
  - the AUGUSTUS source code (src/*.cc, include/*.hh)
  - the scripts (scripts/*.pl)
  - the auxiliary programs (auxprogs/)
  - the tree-parser (src/scanner,src/parser)
are under the Artistic Licence (see src/LICENCE.TXT).
