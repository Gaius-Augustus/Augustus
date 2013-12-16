# manual for AUGUSTUS in cgp (comparative gene prediction) mode
# Stefanie Koenig, 15.12.2013

1. INTRODUCTION
2. DEPENDENCIES
3. INSTALLATION
4. RUNNING AUGUSTUS IN CGP MODE
5. RETRIEVING GENOMES FROM A MYSQL DATABASE
6. USING HINTS


1. INTRODUCTION
----------------
The cgp mode is an extension to AUGUSTUS that takes a genome
alignment of two or more closely related species and simultaneously
predicts genes in all input genomes.
Beside the genomes and the alignment, a phylogenetic tree of the species is required input.
AUGUSTUS-cgp can either be used

- de novo, or
- with extrinsic evidence (such as rnaseq or EST data) for one or more species.
  This also includes the annotation of a well-studied genome to predict genes
  in newly sequenced, related genomes.

Both genomes and extrinsic evidence can either be read in from a flat file or 
alternatively retrieved from a mysql database.
Both approaches are described below in more detail. 

This manual assumes that you are already familiar with AUGUSTUS
and that you know how to use AUGUSTUS for gene prediction in a single genome.

2. DEPENDENCIES
-----------------

- GSL (GNU Scientific Library)
- Boost C++ Libraries (>= V1_46_1)
- Flexc++ (lexical scanner generator, tested with V0.94.0)
- Bisonc++ (parser generator, tested with V2.09.03)
- g++ compiler with C++0X support (>= V4.4)

3. INSTALLATION
----------------

a) install all dependencies

   GSL:      download source code from http://www.gnu.org/software/gsl/ and follow the installation instructions
   Boost:    install via package manager:
             > sudo apt-get install libboost-all-dev
   Flexc++:  download source code from http://flexcpp.sourceforge.net/ and follow the installation instructions. (Flexc++
             has several dependencies including the bobcat library. I recommend to use bobcat-3.10.00, if you do not want to install
             the latest g++ compiler)
   Bisonc++: download source code from http://bisoncpp.sourceforge.net/ and follow the installation instructions or, alternatively,
             install bisonc++ via package manager:
	     > sudo apt-get install bisonc++
   g++       install via package manager:
             > sudo apt-get install build-essential

b) generate a scanner and a parser class that parses files in Newick format

  generate a scanner class with Flexc++:
  Switch to the directory trunks/src/scanner/ in your augustus folder and compile the file 'lexer' with Flexc++:

>  flexc++ lexer
   
  The following files are generated: lex.cc scanner.h scanner.ih scannerbase.h
  Add the include directive '#include "../parser/parserbase.h"' to scanner.ih

> echo '#include "../parser/parserbase.h"' >>scanner.ih

  generate a parser class with Bisonc++:
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

c) recompile AUGUSTUS with cgp mode enabled

   open the file trunks/common.mk with a text editor and uncomment following line to enable comparative gene prediction

#COMGENEPRED = true

   recompile AUGUSTUS

>  cd src
>  make clean all

  
4. RUNNING AUGUSTUS IN CGP MODE
--------------------------------  
   
In order to call Augustus in the comparative gene prediction mode, 4 mandatory arguments need to be passed:

--species=identifier
          a species for which model parameters are trained. For a list of identifiers see README.txt.
	  Decide on one of the species in your set that you can find in the list of identifiers or that
          comes closest to one of the identifiers. For instance, use the identifier 'human' for a mammal
          data set or the identifier 'chicken' for a bird data set
	  

--speciesfilenames=genomeDirectories.tbl
                   a file containing for each species the directory to its genome file.
		   Each line in 'genomeDirectories.tbl' consists of two tab-separated fields.
		   The first field is a species identifier (has nothing to do with the
		   identifier in --species !!!).
		   The second field is the directory to the genome file, e.g.
                   
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
          The sequence names (first field in an 's' line) must be the species identifiers (as they appear in 'genomeDirectories.tbl')
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

           Note that all branch lengths are required and the leaf nodes must be named after the species identifier (as
           in 'aln.maf' and 'genomeDirectories.tbl'). Also a valid format (often output of phylogenetic
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

           In cases where the phylogeny is not known, a star like tree with uniform branch lengths might be used instead, e.g.

(hg19:1,rheMac2:1,mm9:1,bosTau4:1,galGal3:1);

example usage:

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --speciesfilenames=genomeDirectories.tbl

a small data set for testing can be found in examples/cgp/

5. RETRIEVING GENOMES FROM A MYSQL DATABASE
------------------------------------------------

Reading in the genome files from a flat file may require a lot of RAM. Instead,
sequences can be retrieved from a mysql database.

a.) enabling mysql access:
    follow the instructions in docs/mysql.install.readme to install a mysql client and compile the mysql++ library
    
b.) creating a mysql database (example code):

>mysql -u root -p
>create database saeuger;
>select password('AVglssd8');
>create user `cgp`@`%` identified by password '*9D01B966C9648BD3B72A75CEB20A7BCCD41EDE5D'; /* or whatever the password code is*/
>grant all privileges on saeuger.* to cgp@'%';

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


6. USING HINTS
---------------

Extrinsic evidence (or hints) can be integrated using a flat file or database access.
Note that you have to retrieve BOTH genomes and hints either from a flat file or
from the database. Mixed combinations are not possible.

Let's assume we have extrinsic evidence for human and mouse and already prepared the hints files for human and mouse in gff format
(just as you would do it in the single species version of AUGUSTUS):

human.hints.gff contains hints from human rnaseq data

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

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --speciesfilenames=genomeDirectories.tbl --hintsfile=hints.gff --extrinsicCfgFile=cgp.extrinsic.cfg

b) retrieving hints from database
   
   loading hints into the database works exactly the same as loading genomes into the database. Call 'load2db' to
   load hints for a particual species. Use the same species identifier as for the genomes:

> load2db --species=hg19 --dbaccess=saeuger,greifserv2,cgp,AVglssd8 human.hints.gff
> load2db --species=mm9 --dbaccess=saeuger,greifserv2,cgp,AVglssd8 mouse.hints.gff

   prepare the extrinsic config file. Use config/extrinsic/cgp.extrinsic.cfg as template

   call AUGUSTUS with --dbhints enabled:

> augustus --species=human --treefile=tree.nwk --alnfile=aln.maf --dbaccess=saeuger,localhost,cgp,AVglssd8 --dbhints=true --extrinsicCfgFile=cgp.extrinsic.cfg
