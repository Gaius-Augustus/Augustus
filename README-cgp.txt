# manual for AUGUSTUS in comparative gene prediction (cgp) mode
# genes are predicted simulteneously in several aligned genomes
# Stefanie Koenig, September 25th, 2015

 1. INTRODUCTION
 2. DEPENDENCIES
 3. INSTALLATION
 4. RUNNING AUGUSTUS IN CGP MODE
 5. OPTIONAL ARGUMENTS
 6. RETRIEVING GENOMES FROM A MYSQL DATABASE
 7. USING HINTS
 8. SQLITE ACCESS
 9. OPTIMIZING CGP PARAMETERS
10. BUILDING THE NEWICK PARSER FROM SCRATCH
    (not needed unless you run into compiler errors related to 'parse.cc' or 'lex.cc')

1. INTRODUCTION
----------------
The cgp mode is an extension to AUGUSTUS that takes an alignment of two or more genomes
and simultaneously predicts genes in all of them.
Beside the genomes and the alignment, a phylogenetic tree of the species is required input.
AUGUSTUS-cgp can either be used

- de novo, or
- with extrinsic evidence for any subset of species
  Such evidence includes for example already existing and trusted gene structures or hints from RNA-Seq alignments.


Both genomes and extrinsic evidence can either be read in from a flat file or 
alternatively retrieved from a MYSQL or SQLITE database.
All three approaches are described below in more detail.

This manual assumes that you are already familiar with AUGUSTUS
and that you know how to use AUGUSTUS for gene prediction in a single genome.

2. DEPENDENCIES
-----------------

The following programs need be installed in cgp mode:

- GSL (GNU Scientific Library)
- Boost C++ Libraries (>= V1_46_1)
- g++ compiler with C++0X support (>= V4.4)
- lpsolve (mixed linear integer programming)

3. INSTALLATION
----------------

a) install all dependencies

   GSL:      download source code from http://www.gnu.org/software/gsl/ and follow the installation instructions
   Boost:    install via package manager, on UBUNTU/Debian linux:
             > sudo apt-get install libboost-all-dev
   g++       install via package manager:
             > sudo apt-get install build-essential
   lpsolve   > sudo apt-get install liblpsolve55-dev	     

  optional (for gzipped input):
   zlib:     The compression library. Download from http://www.zlib.net/ or install via package manager.

b) recompile AUGUSTUS with cgp mode enabled

   open the file common.mk with a text editor and uncomment the following line to enable comparative gene prediction

#COMGENEPRED = true

   recompile AUGUSTUS

>  cd src
>  make clean all

  
4. RUNNING AUGUSTUS IN CGP MODE
--------------------------------  
   
In order to call AUGUSTUS in the comparative gene prediction mode, 4 mandatory arguments need to be passed:

--species=identifier
          a species for which model parameters are trained. For a list of identifiers see README.txt.
	  Decide on one of the species in your set that you can find in the list of identifiers or that
          comes closest to one of the identifiers. For instance, use the identifier 'human' for a mammal
          data set or the identifier 'chicken' for a bird data set
	  

--speciesfilenames=genomes.tbl
                   a file containing for each species the path to its genome file.
		   Each line in 'genomes.tbl' consists of two tab-separated fields.
		   The first field is a species identifier (does not correspond to the
		   identifier in --species !!!).
		   The second field is the directory and file name for the genome file, e.g.
                   
hg19	 /dir/to/genome/human.fa
rheMac2	 /dir/to/genome/rhesus.fa
mm9	 /dir/to/genome/mouse.fa
bosTau4	 /dir/to/genome/cow.fa
galGal3	 /dir/to/genome/chicken.fa

		   The genome files must be in fasta format and may contain the sequences of one or more chromosomes/scaffolds.
		   The file 'mouse.fa' might look as follows

>chr16
AGCTCGCAGTGTTGATGCTTCAGTCTC
>chr3
ccagaggagacagttagtactaaatgcaccaa


--alnfile=aln.maf
          a file containing a multiple sequence alignment of the genomes in MAF format.
          The sequence names (first field in an 's' line) must be the species identifiers (as they appear in 'genomes.tbl')
          and the sequences identifiers (as they appear in  the genome files) delimited by '.', e.g.  

a score=235085.000000
s hg19.chr21                        15725769 27 +  48129895 AGCTATTGCTGTTTATGTCTCAATTTC
s rheMac2.chr3                     163875558 27 - 196418989 AGCTCTTGCTGTTTACGTCTCGATTTC
s mm9.chr16                         75744509 27 +  98319150 AGCTCGCAGTGTTGATGCTTCAGTCTC
s bosTau4.chr1                     138520043 27 - 161106243 AGCTATTGATGTTTATGTCTTCATTTC
s galGal3.chr1                     101466793 21 + 200994015 AGCTCGAGAAG------AGCCATTATA

a score=128487.000000
s hg19.chr21                        15725796 32 +  48129895 CCAGAGGAGAGGGTTAGTACCAAATGCACCAA
s bosTau4.chr1                     138520070 30 - 161106243 CCAGAGGAGA--GTTCATATTGAGTGCACCAA
s mm9.chr16                         75744536 30 +  98319150 TCAGAGAAGA--ACTTGGACAAAGTGCACCCA
s rheMac2.chr3                     163875585 32 - 196418989 CCAGAGGAGACAGTTAGTACTAAATGCACCAA
  
         
--treefile=tree.nwk
           a phylogenetic tree of the species in Newick format, e.g.   

((((hg19:0.032973,rheMac2:0.036199):0.129706,mm9:0.352605):0.020666,bosTau4:0.219477):0.438357,galGal3:0.474279);

           All branch lengths are required and leaf nodes must be named after the species identifier (as
           in 'aln.maf' and 'genomes.tbl'). Also a valid format (often output of phylogenetic
           tree reconstruction tools such as MrBayes, PHYLIP, ...)  is f.i.

begin trees;
        translate
                1       hg19,
                2       rheMac2,
                3       mm9,
                4       bosTau4,
                5       galGal3
                ;
tree con_50_majrule = [&U] ((((1:0.032973,2:0.036199):0.129706,3:0.352605):0.020666,4:0.219477):0.438357,5:0.474279);
end;

           In cases where the phylogeny is not known, a star-like tree with uniform branch lengths might be used instead, e.g.

(hg19:1,rheMac2:1,mm9:1,bosTau4:1,galGal3:1);

example usage:

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --speciesfilenames=genomes.tbl

a small data set for testing can be found in examples/cgp/


5. OPTIONAL ARGUMENTS
------------------------

--/CompPred/dssqthresh=q
  threshold for the inclusion of donor splice sites based on the pattern probability (q in [0,1] )
  q=0.05 means that only dss are considered that have a pattern, such that 5% of true splice site patterns have lower probability.
  q=0 means that all splice site patterns are considered.

--/CompPred/assqthresh=q --/CompPred/assmotifqthresh=q
  thresholds for the inclusion of aceptor splice sites
  (the inclusion of an acceptor splice site depends both on the ASS and the ASS motif threshold)

--/CompPred/onlySampledECs=on/off
  if on, only exons from the sampling of gene structures are taken as the set of possible candidate exons.
  Otherwise additional candidate exons are determined by combining all possible pairs of ASS/DSS
  start/DSS, ASS/stop and start/stop that are within the maximum length of exons (--max_exon_len, default: 12000).
  Turn this flag on, to reduce the overall runtime memory requirements (default: off)

--/CompPred/exon_gain=q_gain --/CompPred/exon_loss=q_loss
  rate of exon gain and rate of exon loss (parameters of the phylogenetic model)
  q_gain and q_loss are positive real numbers

--/CompPred/maxIterations=n
  the maximum number of dual decomposition iterations (default 100).

--/CompPred/only_species=f
  only predict genes for a subset of the species in the phylogenetic tree.
  f is a file that contains the species identifiers, one per line, of the subset 

--UTR=on/off
  predict the untranslated regions in addition to the coding sequence.
  Note that the 3'-UTR, 5'UTR or both can be absent in some genes if candidate UTRs
  perform poorly in the ab initio model and are not supported by extrinsic evidence. Enforce the prediction
  of UTRs with --/CompPred/genesWithoutUTRs=false

--genemodel=partial/complete
   partial      : allow prediction of incomplete genes at the sequence boundaries (default)
   complete     : only predict complete genes

--/CompPred/genesWithoutUTRs=true/false
  if true, all predicted genes are flanked by a 5'- and 3'- untranslated region (with the exception of partial genes at the sequence boundaries).
  this option only makes sense together with --UTR=on.

--noprediction=true/false
  If true, no prediction is made. Useful for getting the gene ranges and homologous candidate exons.

--/CompPred/outdir=path
  send all output files to this directory (default is the current directory)

--softmasking=1
 adds regions with lowercase nucleotides as nonexonpart hints of source "RM"
 If --extrinsicCfgFile is not given, it used the default cgp.extrinsic.cfg with bonus 1.15, if 
 another extrinsic config file is given, it must contain the "RM" source.
  


6. RETRIEVING GENOMES FROM A MYSQL DATABASE
------------------------------------------------

The flat-file option above reads in all genomes into RAM. This may require too much memory, e.g. for a large number
of vertebrate-sized genomes. Also, this is inefficient when many parallel comparative AUGUSTUS runs are started on a
compute cluster. Therefore, another option allows to read only the required sequences from a MYSQL database:

a.) enabling mysql access:
    follow the instructions in docs/mysql.install.readme to install a mysql client and compile the mysql++ library
    
b.) creating a mysql database (example code) and a user:

> mysql -u root -p
> create database saeuger;
> select password('AVglssd8');
> create user `cgp`@`%` identified by password '*9D01B966C9648BD3B72A75CEB20A7BCCD41EDE5D'; /* or whatever the password code is*/
> grant all privileges on saeuger.* to cgp@'%';

c.) loading sequences into the database:
    Use the program 'load2db' in the AUGUSTUS repository.
    Run load2db with the parameter "--help" to view the usage instructions

> load2db --help
  
    Call 'load2db' for each genome, double check that the correct species identifier is used, e.g.

> load2db --species=hg19 --dbaccess=saeuger,localhost,cgp,AVglssd8 dir/to/genome/human.fa
> load2db --species=rheMac2 --dbaccess=saeuger,localhost,cgp,AVglssd8 dir/to/genome/rhesus.fa
> load2db --species=mm9 --dbaccess=saeuger,localhost,cgp,AVglssd8 dir/to/genome/mouse.fa
> load2db --species=bosTau4 --dbaccess=saeuger,localhost,cgp,AVglssd8 dir/to/genome/cow.fa
> load2db --species=galGal3 --dbaccess=saeuger,localhost,cgp,AVglssd8 dir/to/genome/chicken.fa    

d.) running AUGUSTUS with database access:

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=saeuger,localhost,cgp,AVglssd8


7. USING HINTS
---------------

Extrinsic evidence (or hints) can be integrated using a flat file or database access.
Note that you have to retrieve BOTH genomes and hints either from a flat file or
from the database. Mixed combinations are not possible.

Let's assume we have extrinsic evidence for human and mouse and already prepared the hints files for human and mouse in GFF format
(just as you would do it in the single species version of AUGUSTUS):

human.hints.gff contains hints from human RNA-seq and repeat masking

chr21   b2h     intron         	9908433	        9909046		0       .       .       pri=4;src=E
chr21   repmask nonexonpart     10018268        10018612        0       .       .       src=RM
chr21   w2h     ep      	48084612        48084621        41.600  .       .       src=W;pri=4;mult=41;


mouse.hints.gff contains hints from the mouse Refseq annotation

chr10   mm9_refGene     CDS     50409921        50410055        0.000000        +       0       source=M
chr10   mm9_refGene     intron  50410056        50419745        0.000000        +       .       source=M

a) retrieving hints from a flat file

   First concatenate the hints files into a single file. Prepend the species identifier to the sequence identifier (first column) in the hints files:

> cat human.hints.gff | perl -pe 's/(^chr\d+)/hg19\.$1/' >>hints.gff
> cat mouse.hints.gff | perl -pe 's/(^chr\d+)/mm9\.$1/' >>hints.gff

hints.gff now looks as follows

hg19.chr21   b2h     intron          9908433         9909046         0       .       .       pri=4;src=E
hg19.chr21   repmask nonexonpart     10018268        10018612        0       .       .       src=RM
hg19.chr21   w2h     ep              48084612        48084621        41.600  .       .       src=W;pri=4;mult=41;
mm9.chr10    mm9_refGene     CDS     50409921        50410055        0.000000        +       0       source=M
mm9.chr10    mm9_refGene     intron  50410056        50419745        0.000000        +       .       source=M

   prepare the extrinsic config file. Use config/extrinsic/cgp.extrinsic.cfg as template

   call AUGUSTUS (just as in the single species version) with the hints file and the extrinsic config file

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --speciesfilenames=genomes.tbl --hintsfile=hints.gff --extrinsicCfgFile=cgp.extrinsic.cfg

b) retrieving hints from database
   
   loading hints into the database works exactly the same as loading genomes into the database. Call 'load2db' to
   load hints for a particual species. Use the same species identifier as for the genomes:

> load2db --species=hg19 --dbaccess=saeuger,greifserv2,cgp,AVglssd8 human.hints.gff
> load2db --species=mm9 --dbaccess=saeuger,greifserv2,cgp,AVglssd8 mouse.hints.gff

   prepare the extrinsic config file. Use config/extrinsic/cgp.extrinsic.cfg as template

   call AUGUSTUS with --dbhints enabled:

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=saeuger,localhost,cgp,AVglssd8 --dbhints=true --extrinsicCfgFile=cgp.extrinsic.cfg


8. SQLITE ACCESS
------------------

Alternatively to Mysql, sequences and hints can also be accessed using an SQLite database
(in our experience the sqlite access runs more stabe than the mysql).
Other than the Mysql database that stores the full sequences, the SQLite database only stores
file offsets to achieve random access to the genome files.

a) Installation

   To enable access to an SQLITE database, install the package libsqlite3-dev with your package manager or 
   download the SQLite source code from http://www.sqlite.org/download.html/ 
   (tested with  SQLite 3.8.5 ) and install as follows:

   > tar zxvf sqlite-autoconf-3080500.tar.gz
   > cd sqlite-autoconf-3080500
   > ./configure
   > sudo make
   > sudo make install

   If you encounter an "SQLite header and source version mismatch" error, try

   > ./configure --disable-dynamic-extensions --enable-static --disable-shared

   Turn on the flag SQLITE in augustus/trunks/common.mk and recompile AUGUSTUS

b) create an SQLite database and populate it
   Use the program 'load2sqlitedb' in the AUGUSTUS repository.
   Run load2sqlitedb with the parameter "--help" to view the usage instructions

   > load2sqlitedb --help

   example code for loading a genome and a hints file to the database vertebrates.db
   (always load the genome file first, before loading hints):

   > load2sqlitedb --species=chicken --dbaccess=vertebrates.db genome.fa
   > load2sqlitedb --species=chicken --dbaccess=vertebrates.db hints.gff

c) running AUGUSUTS with SQLite db access:
   call AUGUSTUS with parameters --dbaccess AND --speciesfilenames

   > augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=vertebrates.db --speciesfilenames=genomes.tbl

   in order to retrieve hints from the database, enable --dbhints and pass an extrinsic config file

   > augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=vertebrates.db --speciesfilenames=genomes.tbl --dbhints=true --extrinsicCfgFile=cgp.extrinsic.cfg


9. OPTIMIZING CGP PARAMETERS
-------------------------------

The parameters specific to comparative gene prediction can be automatically optimized
similar to the meta parameters in single species gene prediction using the script 'optimize_augustus.pl'.
In short, a range of parameter values is specified for each parameter in a config file with the extension _metapars.cgp.cfg (e.g. human_metapars.cgp.cfg).
Different values in these ranges are tried out in several rounds and values giving highest accuracy are chosen.
In the evaluation step, the external program Eval¹ and a reference gene set are required.

a) Installation of Eval

   The software package eval by Keibler and Brent is required for retrieving accuracy values of predictions.
   It can be downloaded from

   > wget http://mblab.wustl.edu/media/software/eval-2.2.8.tar.gz
   > tar zxvf eval-2.2.8.tar.gz

   add following lines to your .bashrc file to include the perl executable evaluate_gtf.pl to your $PATH environment variable (optional),
   and the perl modules EVAL.pm and GTF.pm to your $PERL5LIB environment variable (mandatory)

   export PATH=$PATH:/path/to/eval-2.2.8
   export PERL5LIB=$PERL5LIB:/path/to/eval-2.2.8

   to check that the installation was successful, run following command

   >  evaluate_gtf.pl -v /path/to/eval-2.2.8/chr22.refseq.gtf /path/to/eval-2.2.8/chr22.twinscan.gtf /path/to/eval-2.2.8/chr22.genscan.gtf

b) Running optimize_augustus.pl for cgp parameter training
   Run optimize_augustus.pl and read the instructions in USAGE 2 for further information

   > optimize_augustus.pl

   exampe code: 

   > optimize_augustus.pl --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=vertebrates.db --speciesfilenames=genomes.tbl --eval_against=hg19 --stopCodonExcludedFromCDS=1 eval.gtf
   
   the file eval.gtf contains a reference gene set for the human genome that is used for evaluation

¹Keibler, E. and M.R. Brent. 2003. "Eval: A software package for analysis of genome annotations." BMC Bioinformatics 4:50.
   

10. BUILDING THE NEWICK PARSER FROM SCRATCH
    (not needed unless you run into compiler errors related to 'parse.cc' or 'lex.cc')
---------------------------------------------------------------------------------------

To parse phylogenetic trees in Newick format, AUGUSTUS-cgp uses a scanner and parser class, generated by Flexc++ and Bisonc++, respectively.
These classes are part of the AUGUSTUS package and can be used in most applications without the need for changing them.
However, if you have trouble compiling the Augustus source code and the compiler error is related to these classes (parse.cc' and 'lex.cc'), it is recommended to
rebuild the scanner and parser class from scratch.
The following will give you a step-by-step instruction, how to do this:

a) installation of Flexc++ and Bisonc++
b) creation of a scanner class with Flexc++
c) creation of a parser  class with Bisonc++
d) recompilation of AUGUSTUS-cgp

a) installation of Flexc++ and Bisonc++

- Flexc++ (lexical scanner generator, tested with V0.94.0) 

  download source code from http://flexcpp.sourceforge.net/ and follow the installation instructions. (Flexc++
  has several dependencies including the bobcat library. I recommend to use bobcat-3.10.00, if you do not want to install
  the latest g++ compiler)

-  Bisonc++ (parser generator, tested with V2.09.03)

   download source code from http://bisoncpp.sourceforge.net/ and follow the installation instructions or, alternatively,
   install bisonc++ via package manager:
   > sudo apt-get install bisonc++

b) creation of a scanner class with Flexc++:

   Switch to the directory src/scanner/ in your augustus folder and compile the file 'lexer' with Flexc++:

> flexc++ lexer
   
  The following files are generated: lex.cc scanner.h scanner.ih scannerbase.h
  Add the include directive '#include "../parser/parserbase.h"' to scanner.ih

> echo '#include "../parser/parserbase.h"' >>scanner.ih

c) creation of a parser class with Bisonc++

   Switch to the directory trunks/src/parser/ in your augustus folder and compile the file 'grammar' with Bisonc++

>  bisonc++ grammar

  The following files are generated: parse.cc parser.h  parser.ih  parserbase.h
  Edit the 'public' part of parser.h such that it looks as follows

class Parser: public ParserBase
{
    public:
        Parser(std::list<Treenode*> *tree, std::vector<std::string> *species, std::istream &in):d_scanner(in), ptree(tree), pspecies(species) {}
        Scanner d_scanner;
        std::list<Treenode*> *ptree;
        std::vector<std::string> *pspecies;

        int parse();

    private:
    ... // more code
};

d) recompilation of AUGUSTUS-cgp
   
   Switch to the directory 'trunks' in your augustus folder and type 

>  make clean all

