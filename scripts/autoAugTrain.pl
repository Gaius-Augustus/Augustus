#!/usr/bin/perl
#

###########################################################################################################################
#                                                                                                                         #
# autoAugTrain.pl                                                                                                         #
# train augustus automatically                                                                                            #
#                                                                                                                         #
# usage:                                                                                                                  #
# perl ./autoAugTrain.pl --trainingset=training.gb --species=sname [OPTIONS]                                              #
# perl ./autoAugTrain.pl --genome=genome.fa --trainingset=training.gff --species=sname [OPTIONS]                          #
# perl ./autoAugTrain.pl --species=sname --utr --est=cdna.f.psl --workingdir=dir --aug=augustus.gff [OPTIONS]             #
#                                                                                                                         #
###########################################################################################################################

use Getopt::Long;
use Cwd;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname);

BEGIN {
    my $path=rel2abs($0);
    $0=$path;
    our $directory = dirname($path);
}
use lib $directory;
use helpMod qw(find checkFile formatDetector relToAbs setParInConfig uptodate);
use Term::ANSIColor qw(:constants);
use strict;


my $genome = '';           # name of sequence file
my $trainingset='';        # name of training set file
my $species='';            # species name
my $utr='';                # default value: without "utr"
my $flanking_DNA=4000;     # length of flanking DNA, default value is 4000
my $cwd = cwd();           # current working directory
my $positionWD = $cwd;     # position where working directory placed, default: current working directory 
my $workDir;               # working directory = $positionWD/autoAugTrain
my $aug;                   # prediction file made by AUGUSTUS
my $estali;                # est file to train untranslated region
my $verbose=0;             # verbose level
my $CRFtrain=0;            # try CRF training in addition to normal HMM training
my $optrounds=1;           # optimization rounds
my $perlCmdString;         # to store perl commands
my $cmdString;             # to store shell commands
my $useexisting=0;         # start with and change existing config and parameter files

my $evalString;

# absolute path of this script
my $absPath=dirname(relToAbs($0));

# usage
my $usage ="\nName: $0\n\n\n"; 
$usage.="Function: train AUGUSTUS automatically\n\n\nUsage:\n\n"; 
$usage.="autoAugTrain.pl [OPTIONS] --trainingset=training.gb --species=sname\n";
$usage.="autoAugTrain.pl [OPTIONS] --genome=genome.fa --trainingset=training.gff --species=sname\n";
$usage.="autoAugTrain.pl [OPTIONS] --genome=genome.fa --trainingset=training.fa --species=sname\n";
$usage.="autoAugTrain.pl [OPTIONS] --genome=genome.fa --species=sname --utr --est=cdna.f.psl --aug=augustus.gff";
$usage.=" --workingdir=autoTrainNum\n\n";
$usage.="--trainingset=training.gb      training.gb is a file with training gene structures in Genbank format\n";
$usage.="--trainingset=training.gff     training.gff is a file with training genes in GFF format\n\n";
$usage.="--trainingset=training.fa      training.fa is a file with training protein sequences in FASTA format\n\n";
$usage.="--genome=genome.fa             FASTA file with DNA sequences for training\n";
$usage.="                               genome.fa should include the corresponding sequences in this case\n\n";
$usage.="--species=sname                species name as used by AUGUSTUS\n\n";
$usage.="--estali=cdna.f.psl            EST alignments are used to guess the UTR and its end point.\n\n";
$usage.="--utr                          Switch it on to train AUGUSTUS with untranslated regions. Off by default\n\n";
$usage.="options:\n\n";
$usage.="--flanking_DNA=length          flanking_DNA length, default:4000\n";
$usage.="--workingdir=/path/to/wd/      In the working directory results and temporary files are stored.\n";
$usage.="                               Default: current working directory\n";
$usage.="                               By case with \"utr\", the directory \"autoTrainRandomNumber\" which made by autoAugTrain.pl\n";
$usage.="                               without \"utr\" is expected.\n";
$usage.="--verbose                      increase the verbosity level of the program by 1 each, default: 1, max level: 3, e.g. say -v -v -v\n";
$usage.="--useexisting                  use and change the present config and parameter files if they exist for 'species'\n";
$usage.="--optrounds=n                  optimization rounds - each meta parameter is optimized this often (default 1)\n";
$usage.="--CRF                          try training as Conditional Random Field. Off by default\n";
$usage.="--aug=augustus.gff             Previous CDS predictions for constructing a training set of UTRs.\n";

# set options
GetOptions('genome=s' => \$genome,
	   'trainingset=s' => \$trainingset,
	   'species=s' => \$species,
	   'utr!' => \$utr,
	   'CRF!' => \$CRFtrain,
	   'flanking_DNA=i' => \$flanking_DNA,
           'workingdir=s' => \$positionWD,
	   'aug=s' => \$aug,
	   'estali=s' => \$estali,
	   'verbose+' => \$verbose,
	   'optrounds=i' => \$optrounds,
	   'useexisting!' => \$useexisting
	   );
if ($flanking_DNA > 10000){
    print "flanking_DNA larger than neccessary ($flanking_DNA). Resetting flanking_DNA to 10000.\n";
    $flanking_DNA = 10000;
}
# show error information and stop the program if $species not specified
die("\nError: Need to specify the species name!\n\n$usage\n") unless ($species);

# show error information and stop the program if $trainingset not specified
chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
$trainingset = checkFile($trainingset, "training set", $usage) if (!$utr);

# show error information and stop the program if the specified $positionWD cannot be found
$positionWD = relToAbs($positionWD);              # with absolute path
die("The working directory not found! Please specify a valid one! \n") unless (-d $positionWD);
chdir "$positionWD" or die ("Error: Could not change directory to $positionWD!\n");

# check the write permission of $positionWD before building of the work directory
if(!(-w $positionWD)){
    print "Don\'t have write permission at $positionWD!\n";
    die("Please use command \"chmod\" to permission or specify another working directory.\n");
}

# build working directory
$workDir = "$positionWD/autoAugTrain";
if (!$useexisting && -d $workDir){
    print STDERR "$workDir already exists. Please move or use --useexisting.\n";
    exit(1);
}
system("mkdir -p $workDir");
    
# check AUGUSTUS_CONFIG_PATH
my $AUGUSTUS_CONFIG_PATH = $ENV{'AUGUSTUS_CONFIG_PATH'};          # the environment varialbe AUGUSTUS_CONFIG_PATH
$AUGUSTUS_CONFIG_PATH=relToAbs($AUGUSTUS_CONFIG_PATH);
die("Error: The environment variable AUGUSTUS_CONFIG_PATH haven't be evaluated.\n") unless $AUGUSTUS_CONFIG_PATH;
die("Error: The environment variable AUGUSTUS_CONFIG_PATH coudn't be found.\n") unless (-d $AUGUSTUS_CONFIG_PATH);

my $configDir="$AUGUSTUS_CONFIG_PATH"."/species/$species";        # the config directory
my $matName = "$species".'_weightmatrix.txt';                     # species_weightmatrix.txt
my $paraName="$species".'_parameters.cfg';                        # species_parameters.cfg
my $metaName="$species".'_metapars.cfg';                          # species_metapars.cfg
my $metaUtrName="$species".'_metapars.utr.cfg';                   # species_metapars.cfg for UTR

print "\n";

train($genome, $trainingset, $species, $positionWD, $flanking_DNA) if (!$utr);
trainWithUTR($estali, $species, $aug) if ($utr);

print "\n\n";

##################################################################################################
#                                       Train without UTR                                        #
##################################################################################################

sub train{
    my $genome=shift;
    my $trainingset=shift;
    my $species=shift;
    my $positionWD=shift;
    my $flanking_DNA=shift;
    # other pars (???)
    my $string;          # temp string for perl-scripts, which will be called in this script
 
    # build all necessary directories
    chdir "$workDir" or die ("Could not change directory to $workDir!\n");
    for(("gbrowse", "training", "predictions", "predictions/abinitio.1",
	 "predictions/hints.E", "predictions/hints.UTR.E", "training/test", "training/utr")){
	mkdir "$_" if (! -d $_)
    }
    print "3 All necessary diretories have been created unter $workDir.\n" if ($verbose>=3);

    if (!uptodate(["$trainingset"], ["$workDir/training/training.gb"])){
	# detect format of the training set file and create training.gb
	my $format=formatDetector($trainingset);
	die("Error: training set file $trainingset has neither Genbank nor GFF nor FASTA format!\n") if(!$format);
	if($format eq "gb"){ # if it has Genbank format, simply make a link under directory $workDir/training/ 
	    print "2 The input training set file has Genbank format.\n" if ($verbose>=2);
	    print "3 Creating a link $workDir/training/training.gb to $trainingset.\n" if ($verbose>=3);
	    chdir "$workDir/training/" or die("Could not change directory to $workDir/training/.\n");
	    $string = "ln -s $trainingset training.gb";
	    print "3 $string\n" if ($verbose>=3);
	    system("ln -s $trainingset training.gb")==0 or die ("failed to execute: $string\n");
	} elsif ($format eq "gff"){ # if it has GFF format, convert it to Genbank format
	    print "2 The input training set file has GFF format.\n" if ($verbose>=2);
	    chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
	    $genome = checkFile($genome, "fasta", $usage);
	    check_fasta_headers($genome);
	    $string = find("gff2gbSmallDNA.pl");
	    chdir "$workDir/training/" or die("Error: could not change directory to $workDir/training/!\n");
	    # check whether gff file contains any entries
	    my $gffLines = 0;
	    open(GFFforCounting, "<", $trainingset) or die("Cannot open training gene structure gff file $trainingset!\n");
		$gffLines++ while (<GFFforCounting>);
	    close(GFFforCounting) or die("Cannot close training gene structure gff file $trainingset!\n");
	    if($gffLines==0){
		print STDERR "Number of lines in gff file was zero! This is likely to cause problems because no training gene genbank entries can be created! etraining will crash when the training gene genbank file is empty!\n";
	    }
	    # The GFF-file needs to be sorted such that for each gene or mRNA the exons are in increasing order 
	    $cmdString="cat $trainingset | perl -pe 's/\t\S*Parent=/\t/' | sort -n -k 4 | sort -s -k 9 | sort -s -k 1,1 > training.gff";
	    print "3 Running \"$cmdString\" ..." if ($verbose >2);
	    system("$cmdString")==0 or die("failed to execute $!\n");
	    print " Finished!\n" if ($verbose >2);
	    
	    # converting trainingset to Genbank format
	    print "3 converting $trainingset to Genbank format file ...\n" if ($verbose>=3);
	    $perlCmdString = "perl $string $trainingset $genome $flanking_DNA training.gb";
	    print "3 $perlCmdString\n" if ($verbose>2);
	    system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	    print "3 training.gb with genbank format has been created under $workDir/seq/training/.\n" if ($verbose>=3);
	} elsif ($format eq "fasta-prot"){
	    print "2 The input training set file has protein FASTA format. Converting to Genbank format...\n" if ($verbose>=2);
	    scipio_conversion($trainingset);
	} elsif ($format eq "fasta-dna"){
	    die ("Cannot train with a DNA FASTA file, only with protein FASTA.\n");
	} else {
	    die ("Unknown file format: $format\n");
	}
    } else {
	print "2 Reusing existing training set training.gb.\n" if ($verbose>=2);
    }

    
    # run new_species to create the parameter files necessary for training AUGUSTUS for a new species
    $string = find("new_species.pl");
    $perlCmdString = "perl $string --species=$species --silent ";
    $perlCmdString .= "--ignore " if ($useexisting);
    $perlCmdString .= "1>new_species.stdout";
    print "2 Running command  \"$perlCmdString\" to create AUGUSTUS config files for new species $species\n" if ($verbose>=2);
    system("$perlCmdString")==0 or die("Program aborted $!\n");
    open(NSPECIES, "new_species.stdout") or die ("Cannot open new_species.stdout");
    while(defined (my $i=<NSPECIES>)){
	my $string=$i;
	if ($string =~ /^creating/){
	    print "3 $i" if ($verbose>=3);
	}
	else{
	    print "2 $i" if ($verbose>=1 && !($string =~ /^\n/));
	}
    }
    close NSPECIES;
    system ("rm -f new_species.stdout new_species.stderr");

    chdir  "$workDir/training" or die ("Could not change directory to $workDir/training!\n");
    print "2 Now training AUGUSTUS parameters for $species.\n" if ($verbose>=1);

    # count how many squences and genes training.gb contains
    my $counter_seq=0;
    my $counter_gen=0;
    open(TS, "training.gb");
    while(<TS>){
	$counter_seq++ if(/^LOCUS/);
	$counter_gen++ if(/^     CDS             /);
    }
    close TS;
    if($counter_gen == 0){
	die("ERROR: training.gb is empty. Possible reasons:\n\ta) features in a provided training gene structure gff file were not compliant with the autoAug.pl pipeline (for instructions read at e.g. http://bioinf.uni-greifswald.de/augustus-training-0.1/help.gsp#structure\n\tb) Scipio failed to generate training gene structures\n\tThis will cause a crash of the autoAug.pl pipeline!\n");
    }
    my $ave=$counter_gen/$counter_seq;
    print "1 training.gb contains $counter_seq sequences and $counter_gen genes," if ($verbose>=1);
    print " each sequence contains $ave gene(s) on average.\n" if ($verbose>=1); 
    
    # stop Pipeline if the number of training genes is lower than 100
    if($counter_gen < 100){
	die("Number of training genes is with $counter_gen too low (at least 100 genes required)! Training aborted.\n");
    }

    # set $v to the smaller one of 200 and int(0.1*$counter_gen)
    my $v=200;     # the number, how many genes the file training.gb.test should contain 
    $v=int(0.1*$counter_gen) unless (200 <= int(0.1*$counter_gen));
    my $split_num_v=int($v/$ave);
    $split_num_v=1 if ($split_num_v==0);
    
    # make training.gb.test and training.gb.train, output information
    if (!uptodate(["training.gb"], ["training.gb.train", "training.gb.test"])){
	$string=find("randomSplit.pl"); # detect where the script randomSplit.pl is
	$perlCmdString="perl $string training.gb $split_num_v";
	print "2 Creating training.gb.train and training.gb.test with \"$perlCmdString\"...\n" if ($verbose>=2);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	print "2 Created training.gb.train and training.gb.test under $workDir.\n" if ($verbose>=2);
    }
    my $gen_num_v_test=`grep -c ^LOCUS training.gb.test`;
    $gen_num_v_test *= 1;          # in order to delete the "\n"
    my $gen_num_v_train = $counter_gen-$gen_num_v_test;
    print "1 test/evaluation set training.gb.test contains $split_num_v sequences and $gen_num_v_test genes.\n" if ($verbose>=1);
    my $t_b_t=$counter_seq-$split_num_v;  # sequence number in training.gb.train
    print "1 training set training.gb.train contains $t_b_t sequences and $gen_num_v_train genes.\n" if ($verbose>=1);

    # set $o to the smaller one of 300 and int(0.9*$counter)
    my $o = 300;     # the number of genes the file training.gb.train.test should contain
    $o = int(0.9*$counter_gen) unless (300 <= int(0.9 * $counter_gen));
    my $split_num_o = int($o/$ave);
    $split_num_o = 1 if ($split_num_o == 0);
    
    # create training.gb.train.test and training.gb.onlytrain, output information
    if (!uptodate(["training.gb.train"], ["training.gb.onlytrain", "training.gb.train.test"])){
	print "2 Creating training.gb.train.test and training.gb.onlytrain:\n" if ($verbose>=2);
	$perlCmdString="perl $string training.gb.train $split_num_o";
	print "1 randomly selecting $split_num_o genes from the training set training.gb.train..." if ($verbose>=1);
	print "3 $perlCmdString ...\n" if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	system("mv training.gb.train.train training.gb.onlytrain")==0 or die("failed to execute: $!\n");
	print "2 training.gb.train.test and training.gb.onlytrain have been made under $workDir.\n" if ($verbose>=2);
    }
    my $gen_num_o_test=`grep -c ^LOCUS training.gb.train.test`;
    $gen_num_o_test*=1;    # in order to delete the "\n"
    my $t_b_o=$t_b_t-$split_num_o;  # sequence number in training.gb.onlytrain
    my $gen_num_o_train=$t_b_t-$gen_num_o_test; # gene number in training.gb.onlytrain
    print "1 training.gb.train.test contains $split_num_o sequences and $gen_num_o_test genes.\n" if ($verbose>=1);
    print "1 training.gb.onlytrain contains $t_b_o sequences and $gen_num_o_train genes.\n" if ($verbose>=1);
      
    # update parameter file with respect to stop codons (included in CDS? frequencies of the three stop codons)
    if (!uptodate(["$workDir/training/training.gb.train"], ["$configDir/$paraName"]) ||
	!(-f "$configDir/${species}_exon_probs.pbl")){
	# set "stopCodonExcludedFromCDS" to true
	chdir "$configDir" or die ("Can not chdir to $configDir.\n");
	print "2 Seting value of \"stopCodonExcludedFromCDS\" in $paraName to \"true\"\n" if ($verbose>=2);
	setParInConfig($paraName, "stopCodonExcludedFromCDS", "true");
	
	# first try with etraining
	chdir "$workDir/training/" or die ("Can not change directory to $workDir/training.");
	print "2 First try with etraining: etraining --species=$species training.gb.train >train.out ..." if ($verbose>=2);
	system("etraining --species=$species training.gb.train 1>train.out 2>train.err")==0 or die("failed to execute: $!\n");
	print " Finished!\n" if ($verbose>=2);
	print "3 train.out and train.err have been made under $workDir/training.\n" if ($verbose>=3);
	
	# set "stopCodonExcludedFromCDS" to false and run etraining again if necessary
	my $err_stopCodonExcludedFromCDS=`grep -c "exon doesn't end in stop codon" train.err`;
	my $err_rate=$err_stopCodonExcludedFromCDS/$t_b_t;
	print "3 Error rate of missing stop codon is $err_rate\n" if ($verbose>=3);
	
	if($err_rate>=0.5){
	    print "3 The appropriate value for \"stopCodonExcludedFromCDS\" seems to be \"false\".\n" if ($verbose>=3);
	    chdir "$configDir" or die ("Can not chdir to $configDir.\n");
	    print "2 Seting value of \"stopCodonExcludedFromCDS\" in $paraName to \"false\"\n" if ($verbose>=2);
	    setParInConfig($paraName, "stopCodonExcludedFromCDS", "false");
	    print "3 Trying etraining again: etraining --species=$species training.gb.train >train.out ..." if ($verbose>=3);
	    chdir "$workDir/training/" or die ("Can not change directory to $workDir/training.");
	    system("etraining --species=$species training.gb.train 1>train.out 2>train.err")==0 or die("failed to execute: $!\n");
	    print " Finished!\n" if ($verbose>=3);
	    print "3 train.out and train.err have been made again under $workDir/training.\n" if ($verbose>=3);
	}
	
	# adjust the stop-codon frequency in species_parameters.cfg according to train.out
	my $freqOfTag;            
	my $freqOfTaa;
	my $freqOfTga;
	open(TRAIN, "train.out") or die ("Can not open file train.out\n");
	while(<TRAIN>){
	    if(/tag:\s*.*\((.*)\)/) {$freqOfTag=$1}
	    elsif(/taa:\s*.*\((.*)\)/){$freqOfTaa=$1}
	    elsif(/tga:\s*.*\((.*)\)/){$freqOfTga=$1}
	}
	close(TRAIN);
	chdir "$configDir" or die ("Can not chdir to $configDir.\n");
	print "2 Setting frequency of stop codons to tag=$freqOfTag, taa=$freqOfTaa, tga=$freqOfTga.\n" if ($verbose>=2);
	setParInConfig($paraName, "/Constant/amberprob", $freqOfTag);
	setParInConfig($paraName, "/Constant/ochreprob", $freqOfTaa);
	setParInConfig($paraName, "/Constant/opalprob", $freqOfTga);
    }
    if (!uptodate(["$workDir/training/training.gb.test"],
		  ["$workDir/training/test/augustus.1.out"])){
	# first test with augustus
	$cmdString = "augustus --species=$species $workDir/training/training.gb.test > $workDir/training/test/augustus.1.out";
	print "2 First evaluation of parameters ...\n2 Excuting \"$cmdString\" ..." if ($verbose>=2);
	system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	print " Finished!\n" if ($verbose>=2);
    }
    # caculate the accuracy
    my $aug_out="$workDir/training/test/augustus.1.out";
    my $target_1 = accuracy_calculator($aug_out);
    print "1 Ab initio prediction accuracy of AUGUSTUS without optimizing, without UTRs, is: $target_1\n" if ($verbose>=1);

    # optimize metaparameters of AUGUSTUS
    if (!uptodate(["$workDir/training/training.gb.train.test", "$workDir/training/training.gb.onlytrain",
		   "$configDir/$metaName"], ["$workDir/training/optimize.out", "$configDir/$paraName"])){
	chdir "$workDir/training/" or die ("Can not chdir to $workDir/training/.\n");
	$string=find("optimize_augustus.pl");
	if($t_b_o==0){
	    $cmdString="perl $string --rounds=$optrounds --species=$species $workDir/training/training.gb.train.test --metapars=$configDir/$metaName > optimize.out";
	} else{
	    $cmdString="perl $string --rounds=$optrounds --species=$species $workDir/training/training.gb.train.test --onlytrain=$workDir/training/training.gb.onlytrain --metapars=$configDir/$metaName > optimize.out";
	}
	print "1 Optimizing meta parameters of AUGUSTUS\n" if ($verbose>=1);
	print "2 Excuting \"$cmdString\" ..." if ($verbose>=2);
	system("$cmdString")==0 or die("failed to execute: $cmdString!\n"); 
	print " Finished!\n" if ($verbose>=2);
	print "2 You can find all info about optimizing in $workDir/training/optimize.out!\n" if ($verbose>=2);
    } else {
	print "1 Skipping optimization of AUGUSTUS metaparameters.\n" if ($verbose>=1);
    }

    if (!uptodate(["training.gb.train"],["train.withoutCRF.out", "./test/augustus.2.withoutCRF.out"])){
	# train augustus again with optimized parameters
	chdir "$workDir/training" or die ("Could not change directory to $workDir/training\n");
	$cmdString = "etraining --species=$species training.gb.train 1>train.withoutCRF.out 2>train.withoutCRF.err";
	print "2 Running $cmdString ..." if ($verbose>1);
	system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	print " Finished!\n" if ($verbose>1);
	
	# another test on the holdout set, the output file augustus.2.out should have the better results
	print "2 augustus --species=$species training.gb.test > ./test/augustus.2.withoutCRF.out\n" if ($verbose>1);
	system("augustus --species=$species training.gb.test > ./test/augustus.2.withoutCRF.out")==0 or die("failed to execute: $!\n");
    }
    # calculate the accuracy again
    $aug_out="$workDir/training/test/augustus.2.withoutCRF.out";
    my $target_2=accuracy_calculator($aug_out);
    print "1 The accuracy after optimizing without CRF-etraining is $target_2\n" if ($verbose>=1);
    if($target_2-$target_1<0){
	print "\n\n0 WARNING: optimization did not improve accuracy! \n\n";
    }
    
    # copy files to avoid overwriting
    print "2 coping parameter files to $species*.withoutCRF\n" if ($verbose>=2);
    chdir "$configDir";
    for(("$species"."_exon_probs.pbl","$species"."_igenic_probs.pbl", "$species"."_intron_probs.pbl")){
	if (-f $_){ # species_intron_probs.pbl may not exist if 'intronless'
	    print "2 cp $_ $_".'.withoutCRF'."\n" if ($verbose>=3);
	    system("cp $_ $_".'.withoutCRF')==0 or die("failed to execute: cp $_ $_".'.withoutCRF');
	}
    }

    my $target_3 = 0;
    if ($CRFtrain) {
	# train augustus again with CRF
	print "1 Starting training with training as Conditional Random Field (CRF)\n" if ($verbose>=1);
	chdir "$workDir/training" or die ("Could not change directory to $workDir/training\n");
	$cmdString="etraining --species=$species training.gb.train --CRF=1 1>train.CRF.out 2>train.CRF.err";
	print "1 Running $cmdString ..." if ($verbose>0);
	system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	print " Finished!\n" if ($verbose>0);

	# another test on the holdout set, the output file augustus.2.out should have the better results 
	print "3 augustus --species=$species training.gb.test > ./test/augustus.2.CRF.out\n" if ($verbose>2); 
	system("augustus --species=$species training.gb.test > ./test/augustus.2.CRF.out")==0 or die("failed to execute: $!\n");

	# calculate the accuracy again
	$aug_out="$workDir/training/test/augustus.2.CRF.out";
	$target_3 = accuracy_calculator($aug_out);
	print "1 The accuracy with CRF is $target_3\n" if ($verbose>0);
	
	if($target_2>$target_3){
	    print "\n\n1 ####### CRF performance is worse than HMM performance #######\n\n\n" if ($verbose>=1);
	    }
    
	# cp config files
	print "2 copy parameter files to $species*.CRF\n" if ($verbose>=2);
	chdir "$configDir";
	for(("$species"."_exon_probs.pbl","$species"."_igenic_probs.pbl", "$species"."_intron_probs.pbl")){
	    print "2 cp $_ $_".'.CRF'."\n" if ($verbose>1);
	    system("cp $_ $_".'.CRF')==0 or die("failed to execute: $!\n");
	}
    
	# if the accuracy doesn't improve with CRF, overwrite the config files with the HMM parameters from last etraining
	if($target_2>$target_3){
	    chdir "$configDir";
	    for(("$species"."_exon_probs.pbl","$species"."_igenic_probs.pbl", "$species"."_intron_probs.pbl")){
		print "3 rm $_\n" if ($verbose >= 3);
		system("rm $_")==0 or die("failed to execute: $!\n");
		print "3 cp $_".".withoutCRF $_ \n" if ($verbose >= 3);
		system("cp $_".".withoutCRF $_")==0 or die("failed to execute: $!\n");
	    }
	    chdir "$workDir/training";
	}
    } else {
	print "2 Skipping CRF training\n" if ($verbose>1);
    }

    if($verbose>0){
	print "1 WARNING: Optimizing didn't improve the accuracy, maybe need to adjust $species"."_metarpars!\n" if($target_2-$target_1<=0 && $target_3-$target_1<=0);
    }
    
    print "2 Finished autoAugTrain without UTR at " . (scalar localtime()) .
	". All files are stored in $workDir\n" if ($verbose>=2);
}

   ############################################################################################################ 
   #                                        train Augustus with UTR                                           #
   ############################################################################################################

sub trainWithUTR{
    my $estali=shift;
    my $species=shift;
    my $aug=shift;
    my $string;               # temp string for full name of help script
    # other options (???)

    die("Error: Option \"species\" not defined. Please specify a valid species!\n$usage") unless ($species);
    
    print "\n1 ####### Now constructing a training set for Untranslated Regions (UTRs).... ########\n" if ($verbose>=1);

    # check if $estali $aug specified with existent files
    chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
    $aug = checkFile($aug, "augustus", $usage);
    $estali = checkFile($estali, "EST alignment", $usage);

    checkWorkDir();
    
    chdir "$workDir/training" or die ("Error: could not change directory to $workDir/training.\n");
    mkdir "utr" unless(-d "utr");


    # make augustus.gtf
    chdir "$workDir/training/utr/" or die ("Error: could not change directory to $workDir/training/utr/\n");
    if (!uptodate(["$positionWD/autoAugPred_hints/predictions/augustus.gtf"],["genemodels.gtf"])){
	system("cp $positionWD/autoAugPred_hints/predictions/augustus.gtf genemodels.gtf");
    }
    
    # search all start and stop codons from genemodels.gtf and write them to the file stops.and.starts.gff
    if (!uptodate(["genemodels.gtf"],["stops.and.starts.gff"])){
	print "3 Now extracting all stop and start codons from genemodels.gtf to stops.and.starts.gff\n" if ($verbose>=3);
	open(CODON, " > stops.and.starts.gff");
	open(AUGUSTUS, "genemodels.gtf") or die ("Cannot open the file genemodels.gtf, please check if it exists.\n");
        while(defined(my $i=<AUGUSTUS>)){
	    if($i=~/\t(start_codon|stop_codon)\t/) {print CODON $i};
	}
	close(CODON);
	close(AUGUSTUS);
	print "3 File stops.and.starts.gff has been created.\n" if ($verbose>=3);
    }
        
    # make a utr-training set
    if (!uptodate([$genome, $estali, "stops.and.starts.gff"],["utr.gb", "utr.gff"])){
	$string=find("makeUtrTrainingSet.pl");
	print "3 Found script $string.\n" if ($verbose>=3);
	$perlCmdString="perl $string stops.and.starts.gff $genome $estali utr";
	print "2 Running command: $perlCmdString ..." if ($verbose>=2);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString\n");
	print " Finished!\n" if ($verbose>=2);
    } else {
	print "2 Reusing existing UTR training set.\n" if ($verbose>=2);
    }

    # make gbrowse file utr.train.gbrowse
    if (!uptodate(["utr.gff"], ["$workDir/gbrowse/utr.train.gbrowse"])){
	$string = find("utrgff2gbrowse.pl");
	print "3 Running \"cat utr.gff | perl $string\"..." if ($verbose>=3);
	system("cat utr.gff | perl $string > utr.train.gbrowse")==0 or die ("failed to execute: $!\n");
	print " Finished! Made file utr.train.gbrowse\n" if ($verbose>=3);
	system("mv utr.train.gbrowse ../../gbrowse")==0 or die ("failed to execute: $!\n");
	print "3 Moved utr.train.gbrowse to $workDir/gbrowse\n" if ($verbose>=3);
    }

   
    if (!uptodate(["utr.gff", "genemodels.gtf"], ["bothutr.lst", "bothutr.test.gb"])){
	# extract subset of genes, where we have both UTRs
	open(UTR, "utr.gff") or die ("Can not open utr.gff\n");
	open(TRLST, ">tr.lst");
	while(<UTR>){
	    s/.*\t(\S+UTR)\t.*transcript_id \"(\S+)\".*/$2\t$1/;
	    print TRLST;
	}
	close(UTR);
	close(TRLST);
	system("cat tr.lst | sort -u > tr_temp.lst")==0 or die("failed to execute: $!\n");
	system("rm tr.lst")==0 or die("failed to execute: $!\n");
	system("mv tr_temp.lst tr.lst")==0 or die("failed to execute: $!\n");
	open(TR, "tr.lst") or die ("Can not open tr.lst!\n");
	open(BOTH, "> bothutr.lst");
	my $g;
	while(<TR>){
	    split; 
	    print BOTH "$_[0]\n" if ($g eq $_[0]); 
	    $g=$_[0];
	}
	close(TR);
	close(BOTH);
	print "3 The subset of genes, where we have both UTRs is bothutr.lst.\n" if ($verbose>=3);
  
	# sort for gff2gbSmallDNA.pl
	system("cat utr.gff genemodels.gtf > genes.gtf_temp")==0 or die ("failed to execute: $!\n");
	open(GENES, "genes.gtf_temp") or die ("Can not open the file genes.gtf_temp");
	open(WRITEGENES, " > genes.gtf_unsort ");
	while(<GENES>){
	    if(/(CDS|UTR)\t/){print WRITEGENES}
	}
	close(GENES);
	close(WRITEGENES);
	system("cat genes.gtf_unsort | sort -n -k 4,4 | sort -s -k 10,10 | sort -s -k 1,1 > genes.gtf")==0 or die ("failed to execute: $!\n");
	$string=find("gff2gbSmallDNA.pl");
	print "3 Found script $string.\n" if ($verbose>=3);
	$perlCmdString="perl $string genes.gtf ../../../seq/genome.fa $flanking_DNA bothutr.test.gb --good=bothutr.lst 1>gff2gbSmallDNA.stdout 2>gff2gbSmallDNA.stderr";
	print "3 Running \"$perlCmdString\" ..." if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	print " Finished!\n" if ($verbose>=3);
    }
    ##################################################################################################
    # make a training set for optimization including m genes with both UTRs and n genes without UTRs #
    # where m is the smaller one between 150 and the number of blocks in bothutr.test.gb             #
    # n is the  smaller one between 50 and the number of blocks in training.gb.train.test            #
    ##################################################################################################
    if (!uptodate(["bothutr.test.gb", "../training.gb.train.test"] , ["train.gb", "onlytrain.gb"])){
	# evaluate m
	my $m;
	# count the block number in bothutr.test.gb
	my $count=0;
	open(TEMP1, "<bothutr.test.gb") or die {"Can not open the file bothutr.test.gb! \n"};
	while(<TEMP1>){
	    $count++ if /LOCUS/;
	}
	close(TEMP1);
	if($count>=150){$m=150}  else {$m=$count};
	if($count<50){
		die( "ERROR: Number of UTR training examples is smaller than 50. Abort UTR training. If this is the only error message, the AUGUSTUS parameters for your species were optimized ok, but you are lacking UTR parameters. Do not attempt to predict genes with UTRs for this species using the current parameter set!\n");
		exit;
	}
	
	# evaluate n
	my $n;
	# count the block number in training.gb.train.test
	$count=0;
	system("cp ../training.gb.train.test t.gb")==0 or die("failed to execute: $!\n");
	open(TEMP2,"t.gb") or die {"Can not open the file t.gb! \n"};
	while(<TEMP2>){
	    $count++ if /LOCUS/;
	}
	close(TEMP2);
	if($count>=50){$n=50} else {$n=$count};
    
	# now extract traininging set for UTR model
	$string=find("randomSplit.pl");
	print "3 Found script $string.\n" if ($verbose>=3);
	$perlCmdString="perl $string bothutr.test.gb $m";
	print "3 Running \"$perlCmdString\"..." if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	print " Finished!\n" if ($verbose>=3);
	$perlCmdString="perl $string t.gb $n";
	print "3 Running \"$perlCmdString\"..." if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	print " Finished!\n" if ($verbose>=3);
	my $delete;
	open(GB, "t.gb.test") or die ("Can not open file t.gb.test!\n");
	open(NOMRNA, "> t.nomrna.test.gb");
	while(<GB>){
	    $delete=1 if /mRNA/; 
	    $delete=0 if /CDS/; 
	    print NOMRNA if (!$delete);
	}
	close(GB);
	close(NOMRNA);
	system("cat t.nomrna.test.gb bothutr.test.gb.test > train.gb")==0 or die("failed to execute: $!\n");
	
	# count how many genes are contained in train.gb
	my $counter_gen=0;
	open(TS, "train.gb");
	while(<TS>){
	    $counter_gen++ if(/^     CDS             /);
	}
	close TS;
	print "1 Have constructed a training set train.gb for UTRs with $counter_gen genes\n" if ($verbose>=1);
	system("rm t.gb.train t.gb t.gb.test t.nomrna.test.gb")==0 or die("failed to execute: $!\n");
	system("grep LOCUS train.gb | perl -pe 's/^LOCUS\s+(\S+)\s+.*/$1/' > train.gb.lst")==0 or die("failed to execute: $!\n");
	print "3 Made file train.gb.lst under $workDir/training/utr/\n" if ($verbose>=3);
	#why do we need train.gb.lst ???#
    
	########################################################
	# create onlytrain training set only used for training #
	########################################################
       
	system("cp ../training.gb.onlytrain ./")==0 or die("failed to execute: $!\n");
	open(ONLYTRAIN, "training.gb.onlytrain") or die ("Can not open training.gb.onlytrain.");
	open(CDSONLY, "> cdsonly.gb");
	
	# delete the mRNA part up to the next CDS tag
	my $delete=0;
	while(<ONLYTRAIN>){
	    $delete=1 if /mRNA/; 
	    $delete=0 if /CDS/; 
	    print CDSONLY if (!$delete);
	}
	close(ONLYTRAIN);
	close(CDSONLY);
	system("rm training.gb.onlytrain")==0 or die("failed to execute: $!\n");
	
	# construct the disjoint sets: remove training UTR genes from onlytrain UTR gene set (train.utronly.gb)
	open(TRAIN, "train.gb") or die ("Error: could not open the file train.gb");
	open(REMOVE, " > remove.lst");
	$/="\n//\n";
	while(<TRAIN>){
	    if (/LOCUS\s+(\S+)_\d+-\d+ .*gene="(\S+)\.t\d+"/s){print REMOVE "$1_$2\n"}
	}
	close(TRAIN);
	close(REMOVE);
	$/="\n";
	$string=find("filterGenes.pl");
	print "3 Found script $string.\n" if ($verbose>=3);
	$perlCmdString="perl $string remove.lst utr.gb > train.utronly.gb";
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	system("cat cdsonly.gb train.utronly.gb > onlytrain.gb")==0 or die("failed to execute: $!\n");
	print "2 Made onlytrain.gb under $workDir/training/utr/\n" if ($verbose>=2);

	# set UTR value to "on" in the file $species_parameters.cfg
	chdir "$configDir" or die ("Error: could not change directory to $configDir.");
	open(PARA, "$paraName") or die ("Can not open file $paraName.");
	open(PARATEMP, "> cfg_temp") or die ("Can not open file cfg_temp.");
	while(defined(my $i=<PARA>)){
	    if($i=~/(UTR\s+)(off)(\s+#)/){$i=~ s/(UTR\s+)(off)(\s+#)/$1on$3/; print PARATEMP $i;  }
	    else{
		print PARATEMP $i;
	}
	}
	close(PARA);
	close(PARATEMP);
	system("rm $paraName")==0 or die("failed to execute: $!\n");
	system("mv cfg_temp $paraName")==0 or die("failed to execute: $!\n");
	chdir "$workDir/training/utr" or die ("Error: could not chdir to $workDir/training/utr!\n"); 
    } else {
	print "2 Skipping rest of training set construction.\n" if ($verbose>=2);
    }
    # start optimization script
    if (!uptodate(["train.gb", "onlytrain.gb"], ["optimize.utr.out"])){
	$string=find("optimize_augustus.pl");
	print "3 Found script $string.\n" if ($verbose>=3);
	$perlCmdString="perl $string --rounds=$optrounds --species=$species --trainOnlyUtr=1 --onlytrain=onlytrain.gb  --metapars=$configDir/$metaUtrName train.gb --UTR=on > optimize.utr.out";
	print "1 Now optimizing meta parameters of AUGUSTUS for the UTR model." .
	    " This will likely run for a long time..\n" if ($verbose>=1 && $optrounds > 0);
	print "1 Running \"$perlCmdString\"..." if ($verbose>=1);
	system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
	print " finished" if ($verbose>=3);
	system("cat train.gb ../training.gb.train > train.all.gb")==0 or die("failed to execute: $!\n");
    } else {
	print "1 Skipping UTR parameter optimization. Already up to date.\n" if ($verbose>=1);
    }
    # training with final parameters
    $perlCmdString = "etraining --species=$species --UTR=on train.all.gb > etrain.out";
    print "\n2 Running $perlCmdString ..." if ($verbose>=2);
    system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
    if ($CRFtrain){
	# TODO: check what is better at the end
	system("etraining --species=$species --UTR=on train.all.gb --CRF=1 1>etrain.out 2>etrain.err")==0 or die("failed to execute: $!\n");
    }

    print " Finished\n" if ($verbose>=2);
}


########################sub functions###################

# for trainWithUTR,  check if the input working directory contains all necessary files and sub directories

sub checkWorkDir{

    # checking directories
#    foreach my $d ("$positionWD/seq", "$workDir/training", "$positionWD/cdna/alignments", 
#		  "$positionWD/predictions/hints.E"){
#	if (! -d $d){
#	    print STDERR "Directory $d missing.\n";
#	    exit(1);
#	}
#    }
    # checking files
    foreach my $d ("$workDir/training/training.gb.train.test", "$workDir/training/training.gb.onlytrain", 
		   "$positionWD/autoAugPred_hints/predictions/augustus.gtf"){
	if (! -f $d){
	    print STDERR "File $d missing.\n";
	    exit(1);
	}
    }
    
    print "3 File check OK.\n" if ($verbose>=3);   
}

# calculate the result of test-AUGUSTUS
    
sub accuracy_calculator{
    my $aug_out=shift;
    my ($nu_sen, $nu_sp, $ex_sen, $ex_sp, $gen_sen, $gen_sp);
    open(AUGOUT, "$aug_out") or die ("Could not open $aug_out!\n");
    while(<AUGOUT>){
        if(/^nucleotide level\s*\|\s*(\S+)\s*\|\s*(\S+)/){
            $nu_sen=$1;
            $nu_sp=$2;
        }
        if(/^exon level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $ex_sen=$1;
            $ex_sp=$2;
        }
        if(/^gene level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $gen_sen=$1;
            $gen_sp=$2;
        }
    }
    my $target=(3*$nu_sen+2*$nu_sp+4*$ex_sen+3*$ex_sp+2*$gen_sen+1*$gen_sp)/15;
    return $target;
}

#                                                                                                                                                                                            # Check fasta header formatting                                                                                                                                                              #                                                                                                                                                                                             
sub check_fasta_headers{
    my $fastaFile=shift;
    my $someThingWrongWithHeader = 0;
    my $spaces = 0;
    my $orSign = 0;
    my $stdStr = "This may later on cause problems! If the pipeline turns out to crash, please clean up the fasta headers, e.g. by using simplifyFastaHeaders.pl. This message will be suppressed from now on!\n";
    print "1 Checking fasta headers in file $fastaFile...\n" if ($verbose>=1);
    open(FASTA, "<", $fastaFile) or die("Could not open fasta file $fastaFile!\n");
    while(<FASTA>){
        chomp;
        if($_=~m/\s/){
            if($spaces == 0){
                print "1 - WARNING: Detected whitespace in fasta header of file $fastaFile. ".$stdStr;
                $spaces++;
            }
        }
        if($_=~m/\|/){
            if($orSign == 0){
                print "1 - WARNING: Detected "|" in fasta header of file $fastaFile. ".$stdStr;
                $orSign++;
            }
        }
        if($_=!m/[>a-zA-Z0123456789]/){
            if($someThingWrongWithHeader==0){
                print "1 - WARNING: Fasta headers in file $fastaFile seem to contain non-letter and non-number characters. That means they may contain some kind of special character. ".$stdStr;
                $someThingWrongWithHeader++;
            }
        }
    }
    close(FASTA) or die("Could not close fasta file $fastaFile!\n");
}


#
# convert from protein FASTA format to Genbank format using Scipio
#
sub scipio_conversion{
    my $trainingset = shift;
    chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
    $genome=checkFile($genome, "fasta", $usage);
    check_fasta_headers($genome);
    check_fasta_headers($trainingset);
    chdir "$workDir/training/" or die("Error: could not change directory to $workDir/training/!\n");
    $cmdString = "scipio.pl $genome $trainingset > scipio.yaml 2> scipio.err"; 
    print "3 $cmdString ..." if ($verbose>2); 
    system("$cmdString")==0 or die("Program aborted. Possibly \"scipio\" is not installed or not in your PATH");
    print "Finished.\n" if ($verbose>2);
    $cmdString = "cat scipio.yaml | yaml2gff.1.4.pl --filterstatus=\"incomplete\"> scipio.gff 2> yaml2gff.err"; 
    print "3 $cmdString ..." if ($verbose>2);
    system("$cmdString")==0 or die("Command aborted. Possibly \"scipio\" is not installed or not in your PATH\n");
    print "Finished.\n" if ($verbose>2);
    my $pathstr = find("scipiogff2gff.pl");
    $cmdString = "$pathstr --in=scipio.gff --out=traingenes.gff";
    print "3 $cmdString ..." if ($verbose>2);
    system("$cmdString")==0 or die("Command aborted.\n");
    print "Finished.\n" if ($verbose>2);
    # check whether gff file contains any entries                                                                                                                                     
    my $gffLines = 0;
    open(GFFforCounting, "<", "traingenes.gff") or die("Cannot open training gene structure gff file traingenes.gff!\n");
    $gffLines++ while (<GFFforCounting>);
    close(GFFforCounting) or die("Cannot close training gene structure gff file traingenes.gff!\n");
    if($gffLines==0){
	 print STDERR "Number of lines in gff file produced by scipio was zero! This is likely to cause problems because no training gene genbank entries can be created! etraining will crash when the training gene genbank file is empty!\n";
    }
    # finally convert to Genbank format
    $pathstr = find("gff2gbSmallDNA.pl");
    $perlCmdString="perl $pathstr traingenes.gff $genome $flanking_DNA training.gb";
    print "3 $perlCmdString\n" if ($verbose>2);
    system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
    print "3 training.gb with genbank format has been created under $workDir/seq/training/.\n" if ($verbose>2);      
}
