## RUNNING AUGUSTUS

1. [BASIC USAGE](#basic-usage)
2. [OPTIONAL ARGUMENTS](#optional-arguments)
3. [DATABASE ACCESS](#database-access)
4. [USING HINTS](#using-hints)
5. [TRAINING OF CLADE-SPECIFIC PARAMETERS](#training-of-clade-specific-parameters)  
    (USUALLY NOT REQUIRED!!!)
6. [TRAINING CGP SCORE PARAMETERS](#training-cgp-score-parameters)

# BASIC USAGE
   
In order to call AUGUSTUS in the comparative gene prediction mode, 4 mandatory 
arguments need to be passed:

    --species=identifier
       a species for which model parameters are trained. 
       For a list of identifiers see README.md.
       Decide on one of the species in your set that you can find in the list of 
       identifiers or that comes closest to one of the identifiers.
       For instance, use the identifier 'human' for a mammal data set or the 
       identifier 'chicken' for a bird data set.

    --speciesfilenames=genomes.tbl
       a file containing for each species the path to its genome file.
       Each line in 'genomes.tbl' consists of two tab-separated fields.
       The first field is a genome or species identifier
       (does not correspond to the identifier in --species !!!).
       The second field is the directory and file name for the genome file, e.g.

         hg19      /dir/to/genome/human.fa
         mm9       /dir/to/genome/mouse.fa
         bosTau4   /dir/to/genome/cow.fa
         galGal3   /dir/to/genome/chicken.fa

       The genome files must be in fasta format and may contain the sequences of
       one or more chromosomes/scaffolds.
       The file 'mouse.fa' might look as follows

         >chr16
         AGCTCGCAGTGTTGATGCTTCAGTCTC
         >chr3
         ccagaggagacagttagtactaaatgcaccaa

       For running Augustus-cgp on a subset of genomes, simply delete all lines of
       non-target genomes in --speciesfilenames.
       The alignment and phylogenetic tree need no modification if only a subset 
       of genomes is used.

    --alnfile=aln.maf
       a file containing a multiple sequence alignment of the genomes in MAF format.
       The sequence names (first field in an 's' line) must be the species identifiers
       (as they appear in 'genomes.tbl')
       and the sequences identifiers (as they appear in  the genome files) delimited 
       by '.', e.g.  

         a score=628177.000000
         s	hg19.chr21    2032  36  +  32085  TAGG-----------TCTTGCTTCGCCGCAGGAGCGTGGCGGCGGGG
         s	mm9.chr16     1745  36  +  25968  TAGG-----------TCTTGCTGCGTCGGAGCAACGTGGCAGCAGAG
         s	bosTau4.chr1  1935  36  +  30875  TAGG-----------TCTTGCTCCGCCGGAGGAGCGTGGCGGCAGGA
         s	galGal3.chr1  1000  45  +  22283  CAGTAACTGAGCTATTGCTGCTCTGCTG--AGTGAGCCGGAGCAGGG

         a score=3843.000000
         s	hg19.chr21    2068   9  +  32085  CCATGGCCG
         s	bosTau4.chr1  1971   9  +  30875  CTGGGGCCG
         s	mm9.chr16     1781   9  +  25968  CTCTGGCCT
         
       Alignment rows of species that are not listed in --speciesfilenames are ignored.
         
    --treefile=tree.nwk
       a phylogenetic tree of the species in Newick format, e.g.   

         (((hg19:0.16268,mm9:0.352605):0.020666,bosTau4:0.219477):0.438357,galGal3:0.474279);

       All branch lengths are required and leaf nodes must be named after the 
       genome/species identifier (as in 'aln.maf' and 'genomes.tbl').
       Also a valid format (often output of phylogenetic tree reconstruction tools 
       such as MrBayes, PHYLIP, ...) is for instance

         begin trees;
                translate
                        1       hg19,
                        2       mm9,
                        3       bosTau4,
                        4       galGal3
                        ;
         tree con_50_majrule = [&U] (((1:0.16268,2:0.352605):0.020666,3:0.219477):0.438357,4:0.474279);
         end;

       In cases where the phylogeny is not known, a star-like tree with uniform branch lengths might be
       used instead, e.g.

         (hg19:0.01,mm9:0.01,bosTau4:0.01,galGal3:0.01);

       If --speciesfilenames only contains a subset of the species in --treefile, a subtree of is extracted.

   example usage:

    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --speciesfilenames=genomes.tbl

   a small data set for testing can be found in examples/cgp/


# OPTIONAL ARGUMENTS


a) General Options:
-------------------

    --/CompPred/onlySampledECs=on/off
       if on, only exons from the sampling of gene structures are taken as the set 
       of possible candidate exons.
       Otherwise additional candidate exons are determined by combining all possible 
       pairs of ASS/DSS, start/DSS, ASS/stop and start/stop that are within the 
       maximum length of exons (--max_exon_len, default: 12000).
       Turn this flag on, to reduce the overall runtime memory requirements at 
       the cost of a potential decrease in accuracy.
       (default: off)

    --/CompPred/liftover_all_ECs=on/off
       by default only likely exon candidates (the ones from sampling) are lifted 
       over to the other genomes. If this flag is turned on ALL exon candidates 
       are lifted over to the other genomes. This increases the runtime and memory 
       requirement, but is potentially more accurate 
       (default: off)

    --UTR=on/off
       predict the untranslated regions in addition to the coding sequence.
       Note that the 3'-UTR, 5'UTR or both can be absent in some predicted genes 
       if candidate UTRs perform poorly in the ab initio model and are not supported 
       by extrinsic evidence. Enforce the prediction of UTRs with 
       --/CompPred/genesWithoutUTRs=false
       This option requires that a UTR model was trained for the species specified with 
       --species=...

    --nc=on/off
       simultaneous prediction of coding genes and non-coding genes (mostly lincRNA) 
       (default: off)
       Non-coding genes are only predicted if they have RNA-Seq support.
       This option is experimental, as the scores of exons/introns of non-coding 
       genes in the gene structure graph still need to be defined.
       Usage only intended for developers!

    --genemodel=partial/complete
       partial : allow prediction of incomplete genes at the sequence boundaries (default)
       complete: only predict complete genes

    --/CompPred/genesWithoutUTRs=true/false
       if true, all predicted genes are flanked by a 5'- and 3'- untranslated region
       (with the exception of partial genes at the sequence boundaries).
       this option only makes sense together with --UTR=on.

    --noprediction=true/false
       If true, no prediction is made. Useful for getting the gene ranges and 
       homologous candidate exons.

    --/CompPred/outdir=path
       send all output files to this directory (default is the current directory)

    --printOEs=true/false
       print all homologous candidate exons to the file orthoExons.<species>.gff3 
       (default: off)

    --/CompPred/printOrthoExonAli=true/false
       prints codon alignments of all tuples of homologous candidate exons to 
       the file 'orthoexons_codonAlignment.maf'
       (default: false)

    --/CompPred/printConservationWig=true/false
       prints conservation tracks (in wiggle format) of all syntenic regions to 
       the file <species>.wig
       (default: off)

    --exoncands=true/false
       print all candidate exons to the file exonCands.<species>.gff3 (default: off)

    --softmasking=1
       adds regions with lowercase nucleotides as nonexonpart hints of source "RM".
       This is the preferrable way to deal with repeat (soft) masked genomes.
       If --extrinsicCfgFile is not given, it used the default extrinsic-cgp.cfg 
       with bonus 1.15, if another extrinsic config file is given, it must contain 
       the "RM" source.

    --temperature=t
       heat the posterior distribution for sampling, 0=cold, 7=hottest, take probs 
       to the power of (8-temperature)/8
       A higher temperature tends to include more suboptimal gene structures during sampling.
       (default: 3)

    --optCfgFile=cgp_parameters.cfg
       include parameter file from training. The training uses logistic regression
       to classify candidate exons based on cross-species features like selective 
       pressure (dN/dS), conservation, phylogenetic diversity, etc.
       (default: config/cgp/log_reg_parameters_default.cfg)

    --allow_hinted_splicesites=atac
       comma-separated list of non-canonical splice site pairs that enables the 
       prediction of rare introns with unusual splice sites in addition to the 
       GT-AG and GC-AG introns that are allowed by default.

b) Options to adjust the generation of geneRanges (syntenic regions) from the input alignment
---------------------------------------------------------------------------------------------

    --maxDNAPieceSize=n
       This value specifies the maximal length of a sequence chunk in a geneRange. 
       For longer sequence chunks, the geneRange is cut into several pieces that 
       are processed separately. 
       (default: 200000)

    --/CompPred/maxCov
       Decreasing maxCov punishes multiple overlapping geneRanges more.
       By default, maxCov = 3 penalizes the same region covered by more than 3 alignments

    --/CompPred/covPen
       Increasing the coverage penalty covPen punishes long overlaps between geneRanges more.
       By default, covPen = 0.2 punishes uncovered bases 5 times more than each 
       base covered too often

c) Options to adjust properties of splice sites, exons, introns and genes
-------------------------------------------------------------------------

    --max_exon_len=n
       maximum length of a candidate exon 
       (default: 12000)
       Typically, this needs not be changed.

    --min_intron_len=n
       minimum length of a candidate intron 
       (default: 39)

    --min_coding_len=n
       minimum length of a coding region 
       (default: 102)

    --/CompPred/mil_factor=f
       mean intron length factor (>=1), the higher the less are long introns penalized 
       (default: 1).
       A value of 100 roughly corresponds to not penalizing long introns.
       (does not concern explicit introns from sampling whose lengths are modeled 
       explicitly by a geometric-tail distribution. This only concerns implicit 
       introns that are constructed along auxiliary bars in the gene structure graph).

    --/CompPred/dssqthresh=q
       threshold for the inclusion of donor splice sites based on the pattern 
       probability (q in [0,1] )
       q=0.05 means that only dss are considered that have a pattern, such that 
       5% of true splice site patterns have lower probability.
       q=0 means that all splice site patterns are considered.

    --/CompPred/assqthresh=q --/CompPred/assmotifqthresh=q
       thresholds for the inclusion of acceptor splice sites
       (the inclusion of an acceptor splice site depends both on the ASS and the 
       ASS motif threshold)

d) Options to adjust the scoring function of candidate exons/introns:
---------------------------------------------------------------------

    --/CompPred/omega=on/off
       estimate selective pressure (non-synonymous to synonymous rate ration dN/dS) 
       for each codon alignment of homologous candidate exons
       (default: on)

    --/CompPred/conservation=on/off
       compute an average columnwise conservation score for each tuple of homologous 
       candidate exons  
       (default: on)

    --/CompPred/ec_thold=a
       parameter that is added to the scoring function of candidate exons 
       (default: 0)
       Enables the shifting of the scoring function such that sensitivity (SN) and 
       specificity (SP) are in balance.
       If t>0 is the threshold from logistic regression for which SN and SP are 
       in balance on the training set, then set
       a = log((1/t) - 1)

    --/CompPred/ic_thold=b
       parameter that is added to the scoring function of candidate introns
       (default: 0).
       Analogous to parameter --/CompPred/ec_thold above.

    --/CompPred/scale_codontree=f
       scaling factor to scale branch lengths in the codon tree to one codon 
       substitution per time unit.
       After applying this factor to each branch length for the input tree, the 
       tree should be scaled for the expected number of CODON substitutions.
       (default: 1)

e) Options to adjust the phylogenetic model:
--------------------------------------------

    --/CompPred/phylo_model=2,3 or 4
       number of states in the phylogenetic model (default: 2)
       model 2: state 1: EC present but not predicted
                state 2: EC present and predicted
                rate Matrix Q = [(-lambda, lambda), (mu, -mu)]
                depends on parameters --/CompPred/exon_gain (lambda) and --/CompPred/exon_loss (mu)
                (see next parameters)       
       model 3: model 2 + state 3: EC not present but alignment present
                additionally depends on --/CompPred/ali_error
                (see next parameters)
       model 4: model 3 + state 4: no alignment present
       models 3 and 4 are experimental as the rate matrices are not well defined.
       (usage only recommended for developers that want to play around with the rate matrices)

    --/CompPred/exon_loss=r
       rate r>0 of exon loss (parameter of the phylogenetic models, see above)
       (default: 0.0001)

    --/CompPred/exon_gain=r
       rate r>0 of exon gain (parameter of the phylogenetic models, see above)
       (default: 0.0001)

    --/CompPred/ali_error=r
       rate r of alignment errors (parameter of the phylogenetic model 3 and 4, see above)
       (default: 0.1)

    --/CompPred/phylo_factor=f
       specifies the influence of the phylogenetic model
       (default: 1).
       The higher f is chosen, the more weight is given to the phylogenetic model, 
       i.e. the more consistent the gene structures are across the species.

f) Options to adjust the DD algorithm:
--------------------------------------

    --/CompPred/dd_rounds=r
       the number of Dual Decomposition rounds (default: 5).

    --/CompPred/maxIterations=n
       the maximum number of Dual Decomposition iterations per round
       (default: 500)

    --/CompPred/dd_step_rule=harmonic/square_root/base_2/base_e/polyak/constant/mixed
       the step size function (default: mixed)
       - constant:       c
       - harmonic:       c / (v+1)
       - square_root:    c / sqrt(v+1)
       - base_2:         c / (2^v)
       - base_6:         c / (e^v)
       - polyak:         (d_t - p_best) / numInconsistencies
       - mixed:          1. round: polyak, all other rounds: square_root
       where c is the step size parameter (see next parameter) and 
       v is the number of iterations prior to the current iteration, in which the 
       value of the dual  problem increases. The polyak step size adjusts the step 
       size dynamically from quantities computed in previous iterations:
       d_t is the current dual value, 
       p_best is the best primal value, seen so far and
       numInconsistencies is the current number of inconsistencies between the 
       complicating variables.

    --/CompPred/dd_factor=a-b
       value range of the step size parameter c (default: 1-4). Only required for 
       step size rules other than "polyak".
       When only a single round of DD (--/CompPred/dd_rounds=1) is chosen, specify 
       a single value for the step size parameter,
       e.g. --/CompPred/dd_factor=a. For r>1 rounds of DD, the value range [a-b] 
       is split into equidistant values, e.g. for r=4, a=1 and b=4, the values 
       1,2,3 and 4 are used for the first, second, ... and fourth round of DD, 
       respectively.

# DATABASE ACCESS

   The flat-file option above reads in all genomes into RAM. This may require too 
   much memory, e.g. for a large number of vertebrate-sized genomes. Also, this is 
   inefficient when many parallel comparative AUGUSTUS runs are started on a compute 
   cluster. Therefore, another option allows to read only the required sequences 
   from a database.

## Option 1: SQLite

   Sequences and hints can be accessed using an SQLite database 
   (in our experience the SQLite access runs more stable than MySQL).  
   Other than the MySQL database that stores the full sequences, the SQLite database 
   only stores file offsets to achieve random access to the genome files.

a) create an SQLite database and populate it  
   Use the program 'load2sqlitedb' in the AUGUSTUS repository.
   Run load2sqlitedb with the parameter "--help" to view the usage instructions

    load2sqlitedb --help

   example code for loading genome files to the database vertebrates.db: 

    load2sqlitedb --species=hg19 --dbaccess=vertebrates.db human.fa
    load2sqlitedb --species=mm9 --dbaccess=vertebrates.db mouse.fa
    load2sqlitedb --species=bosTau4 --dbaccess=vertebrates.db cow.fa
    load2sqlitedb --species=galGal3 --dbaccess=vertebrates.db chicken.fa

b) running AUGUSTUS with SQLite db access:  
   call AUGUSTUS with parameters --dbaccess AND --speciesfilenames

    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --dbaccess=vertebrates.db --speciesfilenames=genomes.tbl
    
## Option 2: MySQL

This is an alternative to the SQLite database from above.

a) creating a MySQL database (example code) and a user:  

    apt install mysql-server
    /* The initial MySQL root account passwords are usally empty */
    mysql -u root -p
    > CREATE DATABASE vertebrates;
    > CREATE USER `cgp`@`%` IDENTIFIED BY 'db_passwd'; /* or any other password */
    > GRANT ALL PRIVILEGES ON vertebrates.* TO cgp@'%';
    > exit

b) loading sequences into the database:  

   Use the program 'load2db' in the AUGUSTUS repository.  
   Run load2db with the parameter "--help" to view the usage instructions

    load2db --help

   Call 'load2db' for each genome, double check that the correct species identifier is used, e.g.

    load2db --species=hg19 --dbaccess=vertebrates,localhost,cgp,db_passwd human.fa
    load2db --species=mm9 --dbaccess=vertebrates,localhost,cgp,db_passwd mouse.fa
    load2db --species=bosTau4 --dbaccess=vertebrates,localhost,cgp,db_passwd cow.fa
    load2db --species=galGal3 --dbaccess=vertebrates,localhost,cgp,db_passwd chicken.fa

c) running AUGUSTUS with database access:

    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=vertebrates,localhost,cgp,db_passwd

# USING HINTS

   Extrinsic evidence (or [hints](RUNNING-AUGUSTUS.md#using-hints) can be 
   integrated using a flat file or database access.  
   Note that you have to retrieve BOTH genomes and hints either from a flat file or
   from the database. Mixed combinations are not possible.  

   Let's assume we have extrinsic evidence for human and mouse and already prepared the 
   hints files for human and mouse in GFF format (just as you would do it in the single 
   species version of AUGUSTUS):  

   human.hints.gff contains hints from human RNA-seq and repeat masking

    chr21       b2h          intron       9908433    9909046         0  .  .  pri=4;src=E
    chr21       repmask      nonexonpart  10018268  10018612         0  .  .  src=RM
    chr21       w2h          ep           48084612  48084621    41.600  .  .  src=W;pri=4;mult=41;

   mouse.hints.gff contains hints from the mouse Refseq annotation

    chr10       mm9_refGene  CDS          50409921  50410055  0.000000  +  0  source=M
    chr10       mm9_refGene  intron       50410056  50419745  0.000000  +  .  source=M

## retrieving hints from a flat file

   First concatenate the hints files into a single file. Prepend the species identifier 
   to the sequence identifier (first column) in the hints files:

    cat human.hints.gff | perl -pe 's/(^chr\d+)/hg19\.$1/' >> hints.gff
    cat mouse.hints.gff | perl -pe 's/(^chr\d+)/mm9\.$1/'  >> hints.gff

   hints.gff now looks as follows

    hg19.chr21  b2h          intron        9908433   9909046         0  .  .  pri=4;src=E
    hg19.chr21  repmask      nonexonpart  10018268  10018612         0  .  .  src=RM
    hg19.chr21  w2h          ep           48084612  48084621    41.600  .  .  src=W;pri=4;mult=41;
    mm9.chr10   mm9_refGene  CDS          50409921  50410055  0.000000  +  0  source=M
    mm9.chr10   mm9_refGene  intron       50410056  50419745  0.000000  +  .  source=M

   prepare the extrinsic config file. Use config/extrinsic/extrinsic-cgp.cfg as template

   call AUGUSTUS (just as in the single species version) with the hints file and 
   the extrinsic config file

    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --speciesfilenames=genomes.tbl --hintsfile=hints.gff --extrinsicCfgFile=cgp.extrinsic.cfg

## retrieving hints from a SQLite database

   Loading hints into the SQLite database works exactly the same as loading genomes 
   into the database. After the genome files are loaded, call 'load2sqlitedb' to 
   load hints for a particular species. Use the same species identifier as for the genomes:

    load2sqlitedb --species=hg19 --dbaccess=vertebrates.db human.hints.gff
    load2sqlitedb --species=mm9 --dbaccess=vertebrates.db mouse.hints.gff

   prepare the extrinsic config file. Use config/extrinsic/extrinsic-cgp.cfg as template

   call AUGUSTUS with --dbhints enabled:
    
    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --dbaccess=vertebrates.db --speciesfilenames=genomes.tbl --dbhints=true \
      --extrinsicCfgFile=cgp.extrinsic.cfg

## retrieving hints from a MySQL database

   Loading hints into the MySQL database works exactly the same as loading genomes 
   into the database. After the genome files are loaded, call 'load2db' to load 
   hints for a particular species. Use the same species identifier as for the genomes:

    load2db --species=hg19 --dbaccess=vertebrates,localhost,cgp,db_passwd human.hints.gff
    load2db --species=mm9  --dbaccess=vertebrates,localhost,cgp,db_passwd mouse.hints.gff

   prepare the extrinsic config file. Use config/extrinsic/extrinsic-cgp.cfg as template

   call AUGUSTUS with --dbhints enabled:

    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --dbaccess=vertebrates,localhost,cgp,db_passwd --dbhints=true \
      --extrinsicCfgFile=cgp.extrinsic.cfg

# TRAINING OF CLADE-SPECIFIC PARAMETERS
   (USUALLY NOT REQUIRED!!!)

   Clade-specific parameters include the rates for exon gain and loss

    --/CompPred/exon_loss=r
    --/CompPred/exon_loss=r

   as well as the scaling factor

    --/CompPred/phylo_factor=f

   If necessary, these parameters can be optimized similar to the meta parameters 
   in single species gene prediction using the script 'optimize_augustus.pl'.  
   In short, a range of parameter values is specified for each parameter in a 
   config file with the extension _metapars.cgp.cfg (e.g. human_metapars.cgp.cfg).  
   Different values in these ranges are tried out in several rounds and values 
   giving highest accuracy are chosen.  
   In the evaluation step, the external program Eval¹ and a reference gene set are required.

a) Installation of Eval

   The software package eval by Keibler and Brent is required for retrieving accuracy 
   values of predictions.  
   It can be downloaded from

    wget http://mblab.wustl.edu/media/software/eval-2.2.8.tar.gz
    tar zxvf eval-2.2.8.tar.gz

   add following lines to your .bashrc file to include the perl executable evaluate_gtf.pl 
   to your $PATH environment variable (optional), and the perl modules EVAL.pm and GTF.pm 
   to your $PERL5LIB environment variable (mandatory)

    export PATH=$PATH:/path/to/eval-2.2.8
    export PERL5LIB=$PERL5LIB:/path/to/eval-2.2.8

   to check that the installation was successful, run following command

    evaluate_gtf.pl -v /path/to/eval-2.2.8/chr22.refseq.gtf /path/to/eval-2.2.8/chr22.twinscan.gtf \
      /path/to/eval-2.2.8/chr22.genscan.gtf

b) Running optimize_augustus.pl for cgp parameter training
   Run optimize_augustus.pl and read the instructions in USAGE 2 for further information

    optimize_augustus.pl

   example code:

    optimize_augustus.pl --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --dbaccess=db.vertebrates --speciesfilenames=genomes.tbl --eval_against=hg19 \
      --stopCodonExcludedFromCDS=1 eval.gtf

   The file eval.gtf contains a reference gene set for the human genome that is used 
   for evaluation.

   ¹Keibler, E. and M.R. Brent. 2003. "Eval: A software package for analysis of genome annotations."
   BMC Bioinformatics 4:50.

# TRAINING CGP SCORE PARAMETERS

  To train the parameters used to score exon and intron candidates you have two options:

1. if your training set (input alignment) is small, run AUGUSTUS as shown in the following example
   ```
    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --speciesfilenames=genomes.tbl --referenceFile=referenceFeature.gff \
      --refSpecies=hg19 --param_outfile=params.cfg
   ```
   This command will not predict genes. Specifying a reference file with --referenceFile 
   will make AUGUSTUS train feature parameters used to score exon and intron candidates 
   using reference coding exons (CDS) and introns provided by the reference file.  
   referenceFeature.gff needs to be in gtf, gff or gff3 format.  
   All lines with type "intron" or "CDS" are used for training.  
   Other lines will be ignored. Note, that the program is case sensitive.  
   Stop codons need to be included in terminal coding exons.  
   Specify the reference species with --refSpecies. The reference species must be 
   one of the clade species and is denoted by the identifier used in the tree or 
   alignment file. Intron and CDS features in referenceFeature.gff must be
   from the reference species. 

   The trained parameters are written to the file params.cfg.

   After training, run AUGUSTUS in cgp mode with params.cfg as optional config file
   ```
    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --speciesfilenames=genomes.tbl --optCfgFile=params.cfg 
   ```
   If --param_outfile is not specified parameters will be written to  
   $AUGUSTUS_CONFIG_PATH/cgp/log_reg_parameters_trained.cfg.  
   Of course the genomes can also be stored in MySQL or SQLite databases.  
   Adjust the commands accordingly.

2. if option 1. takes to long you can parallelize the training by splitting the 
   alignment file and run AUGUSTUS in parallel to collect all relevant attributes 
   of all exon and intron candidates in the training set as follows
   ```
    mkdir run1 run2 ...
    for dir in run*; do
      augustus --species=human --treefile=tree.nwk --alnfile=$dir/aln.maf \
        --speciesfilenames=genomes.tbl --exoncands=1 --/CompPred/outdir=$dir \
        --outfile=$dir/aug.out --errfile=$dir/aug.err &
    done
   ```
   concatenate the following files after all jobs are done
   ```
    cat run*/exonCands.refSpecies.gff3 run*/refSpecies.sampled_GFs.gff \
      run*/orthoExons.refSpecies.gff3 > all_exon_intron_candidates.gff
   ```
   run AUGUSTUS again, now in training mode using all_exon_intron_candidates.gff (this is fast)
   ```
    augustus --species=human --treefile=tree.nwk --referenceFile=referenceFeature.gff \
      --refSpecies=hg19 --trainFeatureFile=all_exon_intron_candidates.gff --param_outfile=params.cfg
   ```
   The trained parameters are written to the file params.cfg.

   After training, run AUGUSTUS in cgp mode with params.cfg as optional config file
   ```
    augustus --species=human --treefile=tree.nwk --alnfile=aln.maf \
      --speciesfilenames=genomes.tbl --optCfgFile=params.cfg
   ```
   If --param_outfile is not specified parameters will be written to  
   $AUGUSTUS_CONFIG_PATH/cgp/log_reg_parameters_trained.cfg.  
   Of course the genomes can also be stored in MySQL or SQLite databases.  
   Adjust the commands accordingly.