#!/usr/bin/perl
#

##########################################################################################################################
#                                                                                                                        #
# autoAugPred.pl                                                                                                         #
# run AUGUSTUS automatically                                                                                             #
#                                                                                                                        #
# usage:                                                                                                                 #
# autoAugPred.pl [OPTIONS] --genome=genome.fa --species=sname                                                            #
# autoAugPred.pl [OPTIONS] --genome=genome.fa --species=sname --noninteractive                                           #
# autoAugPred.pl [OPTIONS] --genome=genome.fa --species=sname --continue --workingdir=/path/to/wd/                       #
#                                                                                                                        # 
##########################################################################################################################

use Getopt::Long;
use Cwd;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);

BEGIN {
    my $path=rel2abs($0);
    $0=$path;
    our $directory = dirname($path);
}
use lib $directory;
use helpMod qw(find checkFile relToAbs uptodate);
use strict;

my $path=dirname(rel2abs($0));        # path of this script

my $workDir;
my $species;
my $utr;                 # default value: without "utr"
my $hints;               # default value: without "hints"
my $extrinsiccfg;        # extrinsic.cfg file
my $genome;              # fasta file
my $cwd=cwd();           # current directory
my $positionWD=$cwd;     # position where working directory placed, default: current working directory
my $workDir;             # working directory
my $splitDir;            # where splited contigs stored
my $shellDir;            #
my $continue;
my $noninteractive;
my $cname='fe';          # cluster name with default "fe"
my $clusterENVDefs = "export SGE_ROOT=/opt/sge";  # commands to define environment variables on cluster
my $SGEqstatPath = "/opt/sge/bin/lx24-amd64/";    # path to executables on Sun Grid Engine (qsub, qstat)
my $verbose=2;           # verbose level
my $singleCPU=0;         # run sequentially on a single CPU
my $nodeNum=20;          # node number in your cluster
my $user_config_path;    # AUGUSTUS_CONFIG_PATH if specified on the command line
my $AUGUSTUS_CONFIG_PATH;# command line path, or environment variable, if not specified.
my $perlCmdString;       # to store perl commands
my $useexisting=0;       # start with and change existing config, parameter and result files
my $help=0;

# absolute path of this script
my $absPath=dirname(relToAbs($0));

# usage

my $usage =  <<_EOH_;

Name: autoAugPred.pl

Function: predict genes with AUGUSTUS on genomes

Usage:

autoAugPred.pl [OPTIONS] --genome=genome.fa --species=sname
autoAugPred.pl [OPTIONS] --genome=genome.fa --species=sname --hints=hintsfile 

--genome=genome.fa             fasta file with DNA sequences for training

--species=sname                species name as used by AUGUSTUS

--continue                     after cluster jobs are finished, continue to compile prediction sets

options:
--useexisting                  use and change the present config and parameter files if they exist for 'species'
--singleCPU                    run sequentially on a single CPU instead of parallel jobs on a cluster
--noninteractive               for Sun Grid Engine users, who have configurated an openssh key
                               with this option AUGUSTUS is executed automatically on the SGE cluster
--workingdir=/path/to/wd/      in the working directory results and temporary files are stored.
--utr                          switch it on to run AUGUSTUS with untranslated regions. Off by default
--hints=hintsfile              run AUGUSTUS using extrinsic information from hintsfile
--extrinsiccfg=hintcfgfile     configuration file with parameters (boni/mali) for hints.
                               default \$AUGUSTUS_CONFIG_PATH/config/extrinsic/extrinsic.cfg
--verbose                      the verbosity level
--remote=clustername           specify the SGE cluster name for noninteractive, default "fe"
--AUGUSTUS_CONFIG_PATH=path    path to augustus/config directory. default: use environment variable
--help                         print this usage info
_EOH_

if (@ARGV < 1){
    print $usage;
    exit(0);
}

GetOptions('workingdir=s' => \$positionWD,
	   'species=s' => \$species,
	   'genome=s' => \$genome,
 	   'utr!' => \$utr,
	   'hints=s' => \$hints,
	   'extrinsiccfg=s' => \$extrinsiccfg,
	   'continue'=> \$continue,
	   'noninteractive!'=> \$noninteractive,
	   'singleCPU!'=> \$singleCPU,
	   'remote=s'=> \$cname,
	   'verbose+' => \$verbose,
	   'nodeNum=i' => \$nodeNum,
	   'cname=s' => \$cname,
	   'AUGUSTUS_CONFIG_PATH=s' => \$user_config_path,
	   'useexisting!' => \$useexisting,
	   'help!' => \$help
	   );

print "\n";

if ($help){
    print $usage;
    exit(0);
}

if (!$user_config_path && !defined($ENV{'AUGUSTUS_CONFIG_PATH'})){
    die("Neither environment variable AUGUSTUS_CONFIG_PATH defined nor specified on the command line.");
}
if ($user_config_path){
    $AUGUSTUS_CONFIG_PATH = checkFile($user_config_path, "AUGUSTUS_CONFIG_PATH", $usage);
} else {
    $AUGUSTUS_CONFIG_PATH = $ENV{'AUGUSTUS_CONFIG_PATH'};
}

if ($singleCPU && $noninteractive){
    print "Options --singleCPU and --noninteractive cannot be used together. Will use --singleCPU only and run sequentially.\n";
    $noninteractive = 0;
}

# check $positionWD
$positionWD = relToAbs($positionWD);
die("The working directory not found! Please specify a valid one! \n") unless (-d $positionWD);

# build directory structure
$workDir = dirBuilder($positionWD,$utr,$hints);
$shellDir = "$workDir/shells";

my $splitN = prepareScript($genome, $species, $positionWD, $utr, $hints);
my @dopreds = (); # list of job indices that have to be (re)done
for (my $i=1; $i <= $splitN; $i++){
    push @dopreds, $i if (!uptodate(["$workDir/shells/aug$i"], ["$workDir/shells/aug$i.out"]) || !`tail aug$i.out | grep command | wc -l`);
}

if($noninteractive){
    # write $workDir to autoAug.log
    print "3 cd $workDir\n" if ($verbose>=3);
    chdir "$workDir" or die("Error: Could not change directory to $workDir.\n");
    noninteractive($shellDir, $cname, \@dopreds) if (@dopreds);
    continue_aug($shellDir, $utr, $hints);
} else { # interactive
    if (@dopreds){
	if ($singleCPU){
	    print "1 running augustus jobs aug" . join(" aug", @dopreds) . " sequentially now\n" if ($verbose >= 1);
	    chdir "$workDir/shells/";
	    for (my $i=1; $i <= $splitN; $i++){
		print "2 running aug$i\n" if ($verbose >= 2);
		system ("./aug$i");
	    }
	    print "1 done with augustus jobs\n" if ($verbose >= 1);
	    continue_aug($shellDir, $utr, $hints);
	} else {
	    print "\nPlease start the augustus jobs $workDir/shells/aug* manually now.\n";
	    print "Either by running them sequentially on a single PC or by submitting ". 
		"these jobs to a compute cluster." .
		" When the jobs are done, simply rerun your original command with --useexisting\n";
	    exit(0);
	}
    } else {
	print "1 Predictions already there. Reusing them.\n" if ($verbose>=1);
    }
    continue_aug($shellDir, $utr, $hints);
}
print "\n";

##################################
#  prepare  scripts for cluster  #
##################################

sub prepareScript{
    my $genome=shift;
    my $species=shift;
    my $positionWD=shift;
    my $utr=shift;
    my $hints=shift;

    # check fasta file
    chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
    $genome=checkFile($genome, "fasta", $usage);
   
    # check whether species is specified
    die("\nPlease specify the species name!\n\n$usage\n") unless ($species);
    
    # check whether the directory with species name exists under $AUGUSTUS_CONFIG_PATH/species
    die("\nError: the config directory for this species doesn't exist under $AUGUSTUS_CONFIG_PATH/species!\n") unless (-d "$AUGUSTUS_CONFIG_PATH/species/$species");
    
    # check and prepare hint files
    
    if ($hints){
	chdir "$cwd" or die("Error: could not change diectory to $cwd.\n");
	$hints=checkFile($hints, "hints", $usage);      # overwrite $hints with absolute path
	if (!uptodate(["$hints"],["$positionWD/hints/hints.gff"])){
	    system("ln -fs $hints $positionWD/hints/hints.gff")==0 or die ("Failed to create hints link.\n");
	}
	if ($extrinsiccfg){
	    # overwrite $extrinsiccfg with absolute path
	    $extrinsiccfg = checkFile($extrinsiccfg, "extrinsic configuration", $usage); 
	} else {
	    $extrinsiccfg = checkFile("$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.E.cfg",
				      "extrinsic configuration", $usage); ;
		print ("2 Using default hint configuration file extrinsiccfg=$extrinsiccfg\n") if ($verbose>=2);
	}
    }

    
    # in $workDir/seq/ build a link with name genome.fa to $genome
    # not used ?!
    # if (!uptodate([$genome], ["$positionWD/seq/genome.fa"])){
    #	system("ln -fs $genome $$positionWD/seq/genome.fa")==0 or die ("failed to execute: $!\n");
    #	print("3 Created a link to $genome with the name genome.fa under $positionWD/seq/\n") if ($verbose>=3);
    # }

    my $string;          # temp string for perl-scripts, which will be called later
    chdir "$positionWD/seq/split" or die ("Cound not change directory to $positionWD/seq/split.\n");
    my $splitDir=cwd();  # where splited contigs stored
    
    # split fasta into subsets
    opendir(SPLITDIR, $splitDir) or die ("Could not open directory $splitDir.\n");
    my @splitfiles = ();
    while (defined(my $fafile = readdir(SPLITDIR))){
	push @splitfiles, $fafile if ($fafile =~ /\.fa/); # directory '.' is incluced in list
    }
    closedir(SPLITDIR);
    if (@splitfiles == 0 || !uptodate([$genome], \@splitfiles)){
	$string = find("splitMfasta.pl");
	
	# calculate minsize
	my $minsize = `grep -v ">" $genome | wc -c`;
	$minsize = int($minsize/$nodeNum)+1;
	
	print "2 splitting genome sequence into subsets of size >= $minsize bp\n" if ($verbose>=2);
	$perlCmdString="perl $string --minsize=$minsize $genome --outputpath=$splitDir";
	print "3 $perlCmdString ..." if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $!\n");
	print "2 The splited fasta files have been placed under $splitDir\n" if ($verbose>=2);
    } else {
	print "2 Using existing split genome FASTA files.\n" if ($verbose>=2);
    }

    my $splitN = `ls $splitDir/*.fa | wc -l`;
    chomp $splitN;
    chdir "$shellDir" or die ("Error: Could not change directory to $shellDir!\n");
    
    my $i;
    for($i=1; $i<=$splitN; $i++){
	if (!uptodate([],["aug$i"])){# creat not previously existing job scripts
	    my $opt_string = "";  # a string for different situations
	    if ($utr) {
		$opt_string .= "--UTR=on --print_utr=on";
	    } else {
		$opt_string .= "--UTR=off"; 
	    }
	    if ($hints){
		$opt_string .= " --hintsfile=../../hints/hints.gff --extrinsicCfgFile=$extrinsiccfg"; 
	    }
	    
	    my $aug="$AUGUSTUS_CONFIG_PATH/../src/augustus";
	    my $base=basename($genome);
	    $base =~ s/(\.fa|\.fna|\.fasta)$//;     # adjust the function in splitMfasta.pl
	    open(AUG, ">aug$i") or die ("Could not open file aug$i!\n");
	    my $augString="$aug $opt_string --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH --exonnames=on --species=$species $splitDir/$base.split.$i.fa > $workDir/shells/aug$i.out";
	    print AUG "$augString";
	    close AUG;
	    system("chmod +x aug$i")==0 or die("failed to execute: $!\n");
	}
    }
    

    if(!$continue && !$noninteractive && ! -f "shellForAug"){
	open(SH, ">shellForAug");
	print SH 'for (( i=1; i<='."$splitN".'; i++))'."\n";
	print SH 'do'."\n";
        print SH '    echo \'PATH="${PATH}":\'"$AUGUSTUS_CONFIG_PATH/../src/" > aug$i_temp'."\n";
        print SH '    echo "export PATH" >> aug$i_temp'."\n";
        print SH '    cat "aug$i" >> aug$i_temp'."\n";
        print SH '    mv aug$i_temp aug$i'."\n";
        print SH '    chmod +x aug$i'."\n";
        print SH '    qsub -cwd "aug$i"'."\n";  #just for sun grid engine
	print SH 'done'."\n";
	close SH;
    }

    system("chmod +x shellForAug") if (-f "shellForAug");    
    print "2 The shell scripts for cluster have been prepared under $shellDir\n" if ($verbose>=2);

    return $splitN;
}

############
# continue #
############

sub continue_aug{
    my $workDir=shift;
    my $utr=shift;
    my $hints=shift;

    # overwrite $workDir with absolute path
    $workDir=relToAbs($workDir);
    die("Error: The working directory not found! Please specify a valid one! \n") unless (-d $workDir);

    my $gFile;
    my $source;
    if (!$utr && !$hints){
	$gFile="augustus.abinitio.gbrowse";
	$source="AUG-ABINIT";
    }
    if (!$utr &&  $hints){
	$gFile="augustus.E.gbrowse";
	$source="AUG-HINTS";
    }
    if ( $utr &&  $hints){
	$gFile="augustus.UTR.gbrowse";
	$source="AUG-UTR";
    }
    # cat result in aug.out

    print "3 cd $workDir\n" if ($verbose>=3);
    chdir "$workDir" or die ("Could not change directory to $workDir!\n");
   
    my @infiles = ();
    my $i=1;
    while(-f "aug$i.out"){
	push @infiles, "aug$i.out";
	$i++;
    }
    
    my $out='aug.out';
    if (!uptodate(\@infiles, [$out, "../predictions/augustus.gff", "../predictions/augustus.gtf",
			      "$workDir/../gbrowse/$gFile"])){
	if(-f $out){
	    print "3 rm $out\n" if ($verbose>=3);
	    system("rm $out")==0 or die("failed to execute: $!\n");
	}
	
	my $i=1;
	my $notfinished=0;
	while(-f "aug$i.out"){
	    print "3 cat aug$i.out >> $out\n" if ($verbose>=3);
	    system("cat aug$i.out >> $out")==0 or die("failed to execute: $!\n");
	    if (`tail aug$i.out | grep command | wc -l` != 1){
		$notfinished++;
		print STDERR "Job aug$i not properly finished.\n";
	    }
	    $i++;
	}
	die("$notfinished augustus job(s) not properly finished.") if ($notfinished);

	# create file augustus.gff
	my $string=find("join_aug_pred.pl");
	print "3 cat $out | $string > augustus.gff\n" if ($verbose>=3);
	system("cat $out | $string > augustus.gff")==0 or die("failed to execute: $!\n");
	
	# make file augustus.aa
	my $string = find("getAnnoFasta.pl");
	$perlCmdString="perl $string augustus.gff";
	print "3 $perlCmdString\n" if ($verbose>=3);
	system("$perlCmdString")==0 or die ("failed to execute: $!\n"); # make augustus.aa
	if ( -f "augustus.aa" ) {# may not exist if nothing was predicted
	    print "3 mv augustus.aa ../predictions\n" if ($verbose>=3);
	    system("mv augustus.aa ../predictions")==0 or die("failed to execute: mv augustus.aa ../predictions\n");
	}

	# make augustus.gtf
	
	open(GTF, ">augustus.gtf") or die("Could not open file augustus.gtf");  # augustus.gtf
	open(AUG, "augustus.gff");
	while(<AUG>){
	    print GTF if /\tAUGUSTUS\t/;
	}
	close GTF;
	close AUG;
	print "3 mv augustus.gtf ../predictions/\n" if ($verbose>=3);
	system("mv augustus.gtf ../predictions/")==0 or die("failed to execute: $!\n");
	print "3 mv augustus.gff ../predictions/\n" if ($verbose>=3);
	system("mv augustus.gff ../predictions/")==0 or die("failed to execute: $!\n");
	
	# make augustus.abinitio.gbrowse/augustus.E.gbrowse/augustus.UTR.gbrowse 
	
	$string=find("augustus2gbrowse.pl");
	print "3 $workDir/../predictions\n" if ($verbose>=3);
	chdir "$workDir/../predictions";

	my $cmdString="cat augustus.gtf | $string | perl -pe 's/AUGUSTUS\t/"."$source\t/' > $gFile";
	print "3 $cmdString\n" if ($verbose>=3);
	system("$cmdString")==0 or die("failed to execute: $!\n");
	print "3 mv $gFile ../gbrowse\n" if ($verbose>=3);
	system("mv $gFile ../gbrowse")==0 or die("failed to execute: $!\n");
        
    } else {
	print "3 Skipping prediction postprocessing.\n" if ($verbose>=1);
    }
    print "2 Done with \"autoAugPred.pl\"\n" if ($verbose>=2);
}

##############################
#  build directory structure #
##############################
    
sub dirBuilder{
    my $position=shift;         # position to place the main directory
    my $utr=shift;
    my $hints=shift;
    my $prefix;                 # prefix of main directory
    my $mainDir;                # main directory
    
    # show error information and stop the program if the $position could not be found
    # overwrite $position with absolute path
    $position=relToAbs($position);          
    die("The working directory not found! Please specify a valid one! \n") unless (-d $position);
    chdir "$position" or die ("Can not chdir to $position!\n");

    # check the write permission of $positio before building of the main directory
    if(!(-w $position)){
        print "Don\'t have write permission at $position!\n";
        die("Please use command \"chmod\" to change permission or specify another position.\n");
    }

    # build main directory
    $mainDir = "$position/autoAugPred_abinitio" if(!$utr && !$hints);
    $mainDir = "$position/autoAugPred_hints" if(!$utr && $hints);
    $mainDir = "$position/autoAugPred_hints_utr" if($utr && $hints);
    $mainDir = "$position/autoAugPred_utr" if($utr && !$hints);
    if (-d $mainDir && !$useexisting && !$continue){
	print STDERR "$mainDir already exists. Please rename or start with --useexisting.";
	exit(1);
    }
    if (! -d $mainDir) {
	system ("mkdir -p $mainDir");
	print "3 Created directory $mainDir.\n" if ($verbose>=3);
    }

    # build all necessary directories
    chdir "$mainDir" or die ("Could not change directory to $mainDir!\n");
    for(("gbrowse","shells", "predictions")){mkdir "$_" if (! -d $_);}
    system ("mkdir -p $position/seq/split");
    system ("mkdir -p $position/hints/");
    print "3 All necessary directories have been created under $mainDir\n" if ($verbose>=3);
    return $mainDir;
}
    
####################
#  noninteractive  #
####################

sub noninteractive{
    my $workDir=shift;
    my $cname=shift;
    my $dopreds=shift; # array of job indices {1,2,...,N} that aren't done yet (e.g. in case of previous crash)

    if($workDir){
	chdir "$workDir" or die ("Error: Could not change directory to $workDir. Please specify a valid one!\n");
    }
    else{
	die("Error: Missing working direcory!\n$usage");
    }

    if(-f "your_job.tmp"){
	system("rm your_job.tmp")==0 or die ("failed to execute: $!\n");
    }

    open(TEMP, "> your_job.tmp") or die ("Error: could not open a temp file to write your job number!\n");
    foreach my $j (@{$dopreds}){
	my $string_1="ssh $cname ";
	my $string_2="\"cd $workDir; $clusterENVDefs;${SGEqstatPath}qsub -cwd -e aug$j.err ";
	my $string_3="aug$j\"";
	my $cmd="$string_1"."$string_2"."$string_3"." 1>qsub$j.stdout 2>qsub$j.stderr";
	print "3 $cmd\n" if ($verbose>=3);
	system("$cmd")==0 or die ("failed to execute: $!\n");
	$cmd="cat qsub$j.stdout";
	my $output = `$cmd 2>/dev/null`;
	print "3 $output" if ($verbose>=3);
        # exception (???)
	if($output=~/^Your job (\d+)/){
	    print TEMP $1."\n";
	}
	$j++;
	sleep(1);
    }
    close TEMP;
    waitJob($workDir, "your_job.tmp");
    print "3 rm your_job.tmp\n" if ($verbose>=3);
    system("rm your_job.tmp")==0 or die ("failed to execute: $!\n");
}

################################################
# wait until done with all jobs in SGE cluster #
################################################

sub waitJob{
    my $position=shift;
    my $file=shift;
    my $list="$position/$file";       
    my $count=0;
    my $sleep_now = 60;  # check job status after this many seconds at first
    my $sleep_max = 600; # eventuelly check only after this many seconds
    my $getrunning="ssh $cname \"$clusterENVDefs;${SGEqstatPath}qstat| grep -f $list |wc -l 2>/dev/null\"";
    my $QSTAT_STATUS=`$getrunning`;
    $QSTAT_STATUS=~ s/\n//;
    while($QSTAT_STATUS!=0){
        if($count==1){
            print "3 QSTAT: still running $QSTAT_STATUS job(s) at " . (scalar localtime()) . "\n" if ($verbose>=3);
            $count=0;
        }
        sleep($sleep_now);
	$sleep_now += 30 if ($sleep_now < $sleep_max);
        $QSTAT_STATUS=`$getrunning`;
        $QSTAT_STATUS=~ s/\n//;
        $count++;
    }
    print "3 QSTAT: done\n" if ($verbose>=3);
}
