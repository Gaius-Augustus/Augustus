#!/usr/bin/perl

#############################################################
# optimize_augustus
# train augustus and automatically optimize the meta parameters
#
# usage: optimize_augustus.pl --species=myspecies train.gb [parameters]
#
#
#
# Mario Stanke, 23.04.2007
#############################################################

use strict;
use IO::File;

my %cmdpars = ( 'species'              => '',
		'train.gb'             => '',
		'metapars'             => '',
		'rounds'               => '5',
		'onlytrain'            => '',
		'kfold'                => '8',
		'pstep'                => '5',
		'AUGUSTUS_CONFIG_PATH' => '',
		'cpus'		       => '1',
		'opt_trans_matrix'     => '',
		'matrix_constraints'   => '',
		'UTR'                  => '',
		'aug_exec_dir'         => '',
		'trainOnlyUtr'         => '',
		'noTrainPars'          => '',
                'translation_table'    => '',
                'genemodel'            => '',
                '/Constant/min_coding_len' => '',
                'sens_spec_bias'       => '1',

                # command line arguments for optimizing cgp parameters
		'treefile'             => '',
		'alnfile'              => '',
                'dbaccess'             => '',
		'speciesfilenames'     => '',
		'eval_exec_dir'        => '',
                'eval_against'         => '',
		'jg'                   => '',
                'jg_exec_dir'          => '',
		'stopCodonExcludedFromCDS' => '',
                'chunksize'            => '5000000');


my $cgp = 0;
my $joingenes=0;

$SIG{INT} = \&got_interrupt_signal;

sub got_interrupt_signal {
    print STDERR "$0 was interrupted.\n";
    if(!$cgp){
	if ($cmdpars{'opt_trans_matrix'} ne ''){
	    if (-s $cmdpars{'opt_trans_matrix'} . ".curopt") {
		system("cp $cmdpars{'opt_trans_matrix'}.curopt $cmdpars{'opt_trans_matrix'}");
		print STDERR "I replaced the transition matrix in $cmdpars{'opt_trans_matrix'} "
		    . "with the currently optimal matrix.\n";
	    }
	} else {
	    print STDERR "\nPlease retrain augustus with the new parameters using etraining.\n";
	}
    }
    exit(1);

}

my $usage = <<'ENDUSAGE';
optimize_augustus.pl    train augustus and automatically optimize the meta parameters.
                        Beside the conventional training of parameters subject to a single
                        genome (USAGE 1), this script can also be used for training of cgp
                        (comparative gene prediction) parameters (USAGE 2).
                        The cgp parameter optimization requires the installation of the 
                        external program Eval for the calculation of accuracy values.
                        For download and installation instructions, see section 9 of README-cgp.txt

USAGE 1 --- single species parameter training and optimization (main usage for most users)

optimize_augustus.pl --species=myspecies train.gb

      myspecies                prefix of the species name
      train.gb                 genbank file for training with bona fide gene structures

OPTIONS

    --metapars=metapars.cfg    metapars.cfg contains the names and their ranges of the
                               meta parameters that are subject to optimization.
                               (default: generic_metapars.cfg)
    --rounds=r                 r is the number of rounds (default: 5)
    --cpus=n                   n is the number of CPUs to use (default: 1)
    --onlytrain=onlytrain.gb   an optional genbank file that is used in addition to train.gb
                               but only for etrain not for intermediate evaluation of accuracy.
                               These genes may e.g. be incomplete.
    --kfold=k                  Make a k-fold cross validation (default: 8)
    --pstep=p                  For integer and floating parameters start with p tests equidistributed
                               in the allowed range of values (default: 5)
    --AUGUSTUS_CONFIG_PATH=d   Specify the config directory d if not set as environment variable
    --opt_trans_matrix=s       Optimize the transition matrix file s. s must be the transition file used.
                               e.g. ../species/nt/generic/generic_trans_shadow_partial.pbl
    --matrix_constraints=s     A file with try list, normed list and bindings.
    --UTR=on                   Turn untranslated region model on for training and prediction.
    --aug_exec_dir=d           Path to augustus and etraining executable. If not specified
                               it must be in \$PATH environment variable.
    --trainOnlyUtr=1           Use this option, if the exon, intron and intergenic models need not be trained. (default: 0)
    --noTrainPars=1            Use this option, if the parameters to optimize do not affect training. The training step (etraining) is omitted completely. (default: 0)
    --sens_spec_bias=f         increase sensitivity weight by factor f. (default: 1)

USAGE 2 --- optimizing cgp (comparative gene prediction) parameters ---

optimize_augustus.pl --species=myspecies --treefile=tree.nwk --alnfile=aln.file --speciesfilenames=genomes.tbl eval.gtf

   myspecies                    prefix of the model species
   tree.nwk                     a phylogenetic tree of the species in NEWICK format
   aln.maf                      an alignment of the genomes in MAF format
   genomes.tbl                  a text file containing the locations of the genomes. Each line starts with the species identifier
                                followed by the location of the corresponding genome file, e.g.
                                hg19 /path/to/human/genome.fa
                                mm9  /path/to/mouse/genome.fa
                                ...
   eval.gtf                     annotation file in GTF format. Accuracy values are computed by comparing the predictions against this reference set.

OPTIONS

    --stopCodonExcludedFromCDS=1 Use this option, if the stop codons are excluded from the CDS features in 'eval.gtf' (default: 0).
    --eval_exec_dir=d            Directory that contains the executable evaluate_gtf.pl from the eval package.
                                 If not specified it must be in \$PATH environment variable.
    --eval_against=s             s is the species identifier to which 'eval.gtf' belongs to. Caution, if not specified, the
                                 reference species in the alignment (first s-line in Maf block) is assumed.
    --chunksize=n                when more than 1 CPU is used, the alignment is split into multiple smaller alignments each approximately covering a
                                 sequence range of length n (in the reference genome). The prediction step is executed in parallel on these chunks (default: 5000000).
    --dbaccess=db                retrieve genomes either from a MySQL or from an SQLITE database. In the SQLITE case, 'db' is a database file
                                 with extension .db, e.g. --dbaccess=vertebrates.db. In the MySQL case, 'db' is a string that contains the connection
                                 information, e.g. --dbaccess=dbname,host,user,passwd (the parameter --speciesfilenames is not required, here).
    --jg=1                       Use this option, if you want to filter out duplicates from the prediction with the external tool 'joingenes' (default: 0,
                                 however --jg=1 is recommended). The tool 'joingenes' is part of the augustus package and can be found in the 'auxprogs' folder.
    --jg_exec_dir=d              Directory that contains the exectuable 'joingenes' (only required when --jg=1)
    --metapars=metapars.cgp.cfg  see usage 1 above (default: generic_metapars.cgp.cfg)
    --cpus=n                     see usage 1 above
    --pstep=p                    see usage 1 above
    --AUGUSTUS_CONFIG_PATH=d     see usage 1 above
    --aug_exec_dir=d             see usage 1 above
    --sens_spec_bias=f           see usage 1 above

ENDUSAGE

my $be_silent = "--/augustus/verbosity=0 --/ExonModel/verbosity=0 --/IGenicModel/verbosity=0 --/IntronModel/verbosity=0 --/UtrModel/verbosity=0 --/genbank/verbosity=0";

##############################################################
# Check the command line
##############################################################

if ($#ARGV<0) {
    print "$usage";
    exit;
}

foreach (@ARGV) {
    if (/--(\w+)=(.*)/){
	if (!exists($cmdpars{$1})){
	    print "unknown parameter: " . $1 . "\n$usage";
	    exit;
	}
	$cmdpars{$1}=$2;
    } else {
	if ($cmdpars{"train.gb"} eq ''){
	    $cmdpars{"train.gb"}=$_;
	} else {
	    print "There can only be one parameter without an explicit name.\n$usage";
	    exit;
	}
    }
}
 
# check whether training of cgp parameters is turned on
if ($cmdpars{"alnfile"} ne "" && $cmdpars{"treefile"} ne "" && ($cmdpars{"speciesfilenames"} ne "" || $cmdpars{"dbaccess"} ne "")){
    $cgp = 1;
    if ($cmdpars{"eval_against"} eq ""){
	if( qx(grep -m 1 "^s" "$cmdpars{'alnfile'}") =~ /^s\s(\S+)\.(\S+)/ ){
	    print STDERR "no species for evaluation specified. Assuming the reference species " . $1 . " from the maf file.\n";
	    $cmdpars{"eval_against"} = $1;
	} else {
	    print "no species for evaluation specified. Please specify to which species the evaluation file corresponds\n$usage";
	    exit;
	}
    }
    $cmdpars{'noTrainPars'} = 1;
    $cmdpars{"opt_trans_matrix"}='';
    $cmdpars{'onlytrain'}='';
    if ($cmdpars{'eval_exec_dir'} =~ /.[^\/]$/) {
	$cmdpars{'eval_exec_dir'} .= '/';
    }
    # check whether the eval package is properly installed (needed for evaluation of the gene prediction)
    if (qx(which "$cmdpars{'eval_exec_dir'}evaluate_gtf.pl") !~ /evaluate_gtf.pl$/){
	die ("eval is not executable. Please add the directory which contains the executable evaluate_gtf.pl to the PATH environment variable or specify the path with --eval_exec_dir.");
    }
    if (qx("$cmdpars{'eval_exec_dir'}evaluate_gtf.pl" 2>&1) =~ /^Can\'t\slocate\s(\w+\.pm)/ ){
	die ("eval is not executable. The perl libary " . $1 . " cannot be located.\n" . 
             "Please add the directory which contains " . $1 . " to the PERL5LIB environment variable, e.g. add the following line to your .bashrc file:\n\n" . 
	     "export PERL5LIB=\$PERL5LIB:/path/to/" . $1 . "\n\n");
    }
    # join genes
    if ($cmdpars{'jg_exec_dir'} =~ /.[^\/]$/) {
	$cmdpars{'jg_exec_dir'} .= '/';
    }
    if ($cmdpars{'jg'} eq '1'){
	$joingenes=1;
    }
    # check whether joingenes is properly installed
    if ($joingenes && qx(which "$cmdpars{'jg_exec_dir'}joingenes") !~ /joingenes$/){
	die ("joingenes is not executable. Please add the directory which contains the executable joingenes to the PATH environment variable or specify the pa
th with --jg_exec_dir.");
    } 
} elsif (!( $cmdpars{"alnfile"} eq "" && $cmdpars{"treefile"} eq "" && $cmdpars{"speciesfilenames"} eq "" && $cmdpars{"dbaccess"} eq "") ){
    die ("For Optimizing cgp parameters you must specify parameters alnfile, treefile\n" . 
         "and one of the following combinations of parameters\n\n" .
        "--dbaccess (retrieving genomes from a MySQL db)\n" .
         "--speciesfilenames and --dbaccess (retrieving genomes from an SQLite db)\n" .
         "--speciesfilenames (retrieving genomes from flat files)\n\n" .
         "Otherwise, none of these parameters are specifieds.\n");
}

if ($cmdpars{"train.gb"} eq ""){
    if(!$cgp){
	print "training file missing\n$usage";
    } else {
	print "annotation file missing\n$usage";
    }
    exit;
}

if ($cmdpars{"species"} eq ""){
    print "no species specified\n$usage";
    exit;
}

if ($cmdpars{"kfold"}<1 && !$cgp) {
    die ("must be at least two fold cross validation");
}

if ($cmdpars{'aug_exec_dir'} =~ /.[^\/]$/) {
    $cmdpars{'aug_exec_dir'} .= '/';
}

# check whether augustus and etraining are executable
if (qx(which "$cmdpars{'aug_exec_dir'}augustus") !~ /augustus$/){
    die ("augustus is not executable. Please add the directory which contains the executable augustus to the PATH environment variable or specify the path with --aug_exec_dir.");
}

if ($cmdpars{'opt_trans_matrix'} ne ''){
    $cmdpars{'noTrainPars'} = 1;
}

if ($cmdpars{'noTrainPars'} eq '') {
    if (qx(which "$cmdpars{'aug_exec_dir'}etraining") !~ /etraining$/){
	die ("etraining is not executable. Please add the directory which contains the executable etraining to the PATH environment variable or specify the path with --aug_exec_dir.");
    }
}

if($cmdpars{'sens_spec_bias'} < 0){
    die("sens_spec_bias parameter must be positive");
}

my $pars="";
if ($cmdpars{"UTR"} eq "on"){
    $pars="--UTR=on";
}
if ($cmdpars{"translation_table"} > 1){
    $pars = $pars." --translation_table=".$cmdpars{"translation_table"};
}
if (length($cmdpars{"genemodel"}) > 1){
    $pars = $pars." --genemodel=".$cmdpars{"genemodel"};
}
if (not($cmdpars{"/Constant/min_coding_len"} eq "")){
    $pars = $pars." --/Constant/min_coding_len=".$cmdpars{"/Constant/min_coding_len"};
}
if($cmdpars{"treefile"} ne ""){
   $pars = $pars." --treefile=".$cmdpars{"treefile"};
}
if($cmdpars{"dbaccess"} ne ""){
   $pars = $pars." --dbaccess=".$cmdpars{"dbaccess"};
}
if($cmdpars{"speciesfilenames"} ne ""){
   $pars = $pars." --speciesfilenames=".$cmdpars{"speciesfilenames"};
}
if (defined($cmdpars{"stopCodonExcludedFromCDS"}) && $cmdpars{"stopCodonExcludedFromCDS"} ne "" && $cmdpars{"stopCodonExcludedFromCDS"} ne "0"){
   $pars = $pars." --stopCodonExcludedFromCDS=true";
}

my $modelrestrict = "";
if (defined($cmdpars{'trainOnlyUtr'}) && $cmdpars{'trainOnlyUtr'} ne "" && $cmdpars{'trainOnlyUtr'} ne "0"){
    $modelrestrict = "--/EHMMTraining/statecount=2 --/EHMMTraining/state00=intronmodel --/EHMMTraining/state01=utrmodel --/IntronModel/outfile=/dev/null";
    # the intron model needs to be trained because the UTR model depends on it
    # possibly change --/IntronModel/outfile=/dev/null if that is incompatible with some UNIX systems
    $pars="--UTR=on" if ($pars eq "");
}

# check whether this perl module for paralell execution is installed
my $got_ForkManager = 0;
eval { require Parallel::ForkManager };
unless ($@) {
  $got_ForkManager = 1;
  Parallel::ForkManager->import();
}

if (!$got_ForkManager && $cmdpars{'cpus'} > 1){
    print STDERR "The perl module Parallel::ForkManager is required to run optimize_augustus.pl in parallel.\n";
    print STDERR "Install this module first. On Ubuntu linux install with\n";
    print STDERR "sudo apt-get install libparallel-forkmanager-perl\n";
    print STDERR "Will now run sequentially (--cpus=1)...\n";
    $cmdpars{'cpus'} = 1;
}

##############################################################
# Create temporary files and folders
##############################################################


my $optdir = "tmp_opt_$cmdpars{'species'}";
system("rm -rf $optdir;mkdir $optdir");
opendir(TMPDIR, $optdir) or die ("Could not open $optdir");

# only needed in cgp
my @gfflines = ();

# Split train.gb into k_fold equal portions
if(!$cgp){
    print "Splitting training file into $cmdpars{'kfold'} buckets...\n";
    open (TRAINGB, <$cmdpars{"train.gb"}>) or die ("Could not open $cmdpars{'train.gb'}");
    
    my @seqlist = ();
    @seqlist  = <TRAINGB>;
    my @namelines = grep /^LOCUS   +/, @seqlist;
    @seqlist = ();
    my @names=();
    
    if (@namelines < $cmdpars{"kfold"}) {
	print "Number of training sequences is too small\n";
	exit;
    }
    
    foreach (@namelines) {
	/LOCUS +([^ ]+) */;
	push @names, $1;
    }
    
    my $bucket=0;
    my %bucketmap=();
    srand(88);
    while ($#names >= 0) {
	my $rand = rand (@names);
	$bucketmap{$names[$rand]}=$bucket;
	$bucket++;
	if ($bucket == $cmdpars{"kfold"}){
	    $bucket = 0;
	}
	splice @names, $rand, 1;          # delete array element
    }
    my $handle;
    my @fh=();
    for ($bucket=1; $bucket<=$cmdpars{"kfold"}; $bucket++) {
	$handle = IO::File->new(">$optdir/bucket$bucket.gb") or die("Could not open bucket$bucket.gb");
	push @fh, $handle;
    }
    
    $/="\n//\n"; # this causes a huge single chunk when DOS carriage returns are used at line breaks
    seek (TRAINGB, 0, 0);
    my $nloci = 0;
    while(<TRAINGB>) {
	my $gendaten=$_;
	m/^LOCUS +(\S+) .*/;
	my $genname=$1;
	
	$bucket = $bucketmap{$genname};
	my $handle = $fh[$bucket];
	print $handle $gendaten;
	$nloci++;
    }
    
    foreach my $handle (@fh) {
	close $handle;
    }
    
    if ($nloci < @namelines){
	die ("Genbank input file appears to have fewer records than expected.\n" . 
	     "This could be a consequence of using DOS (Windows) carriage return symbols at line breaks.");
    }
    
    # create training sets for cross-validation (parallel version)
    if ($cmdpars{'cpus'} > 1){
	for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
	    # make the temporary training and testing files
	    system("rm -f $optdir/curtrain-$k");
	    for (my $m=1; $m <= $cmdpars{"kfold"}; $m++) {
		if ($m != $k) {
		    system("cat $optdir/bucket$m.gb >> $optdir/curtrain-$k");
		}
	    }
	    if ($cmdpars{'onlytrain'} ne ''){
		system ("cat $cmdpars{'onlytrain'} >> $optdir/curtrain-$k");
	    }
	}	
    }
} else { # prepare training of cgp parameters
    # splitting MAF file for parallel computing
    splitMaf();
    # reading in annotation file and checking if it is in a valid GTF format
    print "reading in $cmdpars{'train.gb'}...\n";
    open (ANNO, <$cmdpars{"train.gb"}>) or die ("Could not open $cmdpars{'train.gb'} for reading: $!");
    while(<ANNO>){
	if (/^\s*\#.*/ || /^\s*$/){ # skip comment lines and empty lines
	    next;
	}
	my @line = split(/\t/,$_);
	if(@line < 9){
	    die ("Not GTF format in the following line:\n$_\n");
	}
	if($line[2] eq "CDS" || $line[2] eq "stop_codon" || $line[2] eq "start_codon"){
	    if ($line[8] !~ /transcript_id\s"?[^";]+"?;/){
		die ("Not GTF format in the following line:\n$_\ntranscript_id not found.\n");
	    }
	    if ($line[8] !~ /gene_id\s"?[^";]+"?;/){
		die ("Not GTF format in the following line:\n$_\ngene_id not found.\n");
	    }
	    push @gfflines, $_;
	}	
    }
    close(ANNO);
    # TODO: if stop codons are excluded from CDS features, add --stopCodonExcludedFromCDS=true to $pars
}

##############################################################
# Read in the meta parameters
##############################################################

my @metastartvalues = ();
my @metaparnames = ();
my @metaparranges = ();
my $configdir;
if ($cmdpars{"AUGUSTUS_CONFIG_PATH"} ne "") {
    $configdir = $cmdpars{"AUGUSTUS_CONFIG_PATH"};
} else {
    exists($ENV{AUGUSTUS_CONFIG_PATH}) or die("Environment variable AUGUSTUS_CONFIG_PATH not set.");
    $configdir = $ENV{AUGUSTUS_CONFIG_PATH};
}
if ($configdir !~ /\/$/){
    $configdir .= "/";
}

###############################################
# open the file with the parameters to optimize
#
# either metapars.cfg or opt_transition_matrix
###############################################
my $metaparsfilename;
my $n=0; # number of states
my @trans; # transition matrix, array of array references

if ($cmdpars{'opt_trans_matrix'} eq '') {# optimize meta parameters
   if ($cmdpars{'metapars'} ne '') {
       $metaparsfilename = $cmdpars{'metapars'};
   } elsif(!$cgp) {
       $metaparsfilename = $configdir . "species/generic/generic_metapars.cfg";
   }
   else{
        $metaparsfilename = $configdir . "species/generic/generic_metapars.cgp.cfg";
   }
   open META, '<' , <${metaparsfilename}> or die ("Could not open $metaparsfilename");
   print "Reading in the meta parameters used for optimization from $metaparsfilename...\n";

   $/="\n";
   while (<META>) {
       my ($parname, $range);
       if (/^\s*\#/){
	   next;
       }
       if (/^\s*(\S+)\s+(.*)\s*/) {
	   $parname = $1;
	   $range=$2;
	   if (idx(\@metaparnames, $1) != -1 ){
	       die ("Meta parameter $1 occurs twice in $metaparsfilename.");
	   } else {
	       if ($range =~ /"([^"]+)"-"([^"]+)"/) {
		   push @metaparnames, $parname;
		   push @metaparranges, ['intrange', $1, $2];
	       } elsif ($range =~ /"([^"]+)"_"([^"]+)"/) {
		   push @metaparnames, $parname;
		   push @metaparranges, ['floatrange', $1, $2];
	       } else {
		   my @tokens = split /\s+/, $2;
		   my @list = ();
		   foreach (@tokens){
		       s/^"(.*)"$/$1/;
		       push @list, $_;
		   }
		   push @metaparnames, $parname;
		   push @metaparranges, ['list', @list];
	       }
	   }
       }
   }
} else {# read in transition matrix for optimization
   open TRANS, '<' , <$cmdpars{'opt_trans_matrix'}> 
      or die ("Could not open transition matrix file $cmdpars{'opt_trans_matrix'}");
   print "Reading in the transition matrix...\n";
   $/="\n";
   while (<TRANS>) {
       my ($from, $to, $prob);
       if (/^\s*\#/){
	   next;
       }
       if ($n==0 && /(\d+)/){
	   $n=$1;
	   print "Transition matrix has dimension ${n}x${n}.\n";
       }
       if (/^\s*(\d+)\s+(\d+)\s*(\S+)/) {
	   $from = $1;
	   $to   = $2;
	   $prob = $3;
	   #print "trans[$from][$to]=$prob\n";
	   if ($from < 0 || $from >= $n || $to < 0 || $to >= $n ){
	       print "State of transition matrix out of bounds ($n) for transition $from->$to:$prob\n";
	   }
	   if ($prob < 0) {
	       print "Error: negative probability in transition $from->$to:$prob\n";
	   }
	   if (!defined ($trans[$from])){
	       $trans[$from]=[];
	   }
	   $trans[$from][$to] = $prob;  
       }
   }
}
#printmetaranges(\@metaparnames, \@metaparranges);


# open species_parameters.cfg
my (@spcfilelines, @transfilelines);
my  $speciesdir = $configdir . "species/$cmdpars{'species'}/";
my  $species_cfg_filename = $speciesdir . $cmdpars{'species'} . "_parameters.cfg";
if ($cmdpars{'opt_trans_matrix'} eq ''){
    # make a copy of the original parameter file
    my $y=1;
    while ($y<20 && sysopen(ORIG, "$species_cfg_filename.orig$y",O_WRONLY|O_EXCL)){
	$y++;
    }
    if ($y<20) {
	close(ORIG);
	system("cp $species_cfg_filename $species_cfg_filename.orig$y"); 
    } else {
	die("Too many $species_cfg_filename.orig copies. Please delete some.");
    }

    if(-e "$species_cfg_filename"){
	open(SPCCFG, "<$species_cfg_filename") or die ("Could not open $species_cfg_filename");
    }else{die "File $species_cfg_filename does not seem to exist!\n";}
    print "Reading in the starting meta parameters from $species_cfg_filename...\n";
    $/="\n";
    @spcfilelines = <SPCCFG>;
    close (SPCCFG);
    foreach (@spcfilelines) {
	my ($parname, $value);
	if (/^\s*\#.*/ || /^\s*$/){ # skip comment lines
	    next;
	}
	if (/^\s*(\S+)\s+(\S*)\s*/) {
	    $parname = $1;
	    $value = $2;
	    my $index = idx(\@metaparnames, $1);
	    if($index != -1){
		$metastartvalues[$index] = $value;
	    }
	}
    }

    for (my $i=0; $i<= $#metaparnames; $i++) {
	if (!defined $metastartvalues[$i]){
	    die ("No start value for parameter $metaparnames[$i] found in file $species_cfg_filename.\n\
Maybe you misspelled this parameter in $metaparsfilename.\n");
	}
    }
} else {
    # make a copy of the original transition matrix file
    my $y=1;
    while ($y<40 && sysopen(ORIG, "$cmdpars{'opt_trans_matrix'}.orig$y",O_WRONLY|O_EXCL)){
	$y++;
    }
    if ($y<40) {
	close(ORIG);
	system("cp $cmdpars{'opt_trans_matrix'} $cmdpars{'opt_trans_matrix'}.orig$y"); 
    } else {
	die("Too many $cmdpars{'opt_trans_matrix'}.orig copies. Please delete some.");
    }

    open(TRANS, "<$cmdpars{'opt_trans_matrix'}") or die ("Could not open $cmdpars{'opt_trans_matrix'}");
    $/="\n";
    @transfilelines = <TRANS>;
    close (TRANS);
}

# check if file $cmdpars{'onlytrain'} exists
if ($cmdpars{'onlytrain'} ne ''){
    $cmdpars{'onlytrain'} =~ s/^~/$ENV{HOME}/;
    open(ONLYTRAIN, <$cmdpars{'onlytrain'}>) or die ("Could not open $cmdpars{'onlytrain'}");
    close(ONLYTRAIN);
}

#######################################################################################
# initialize and first test
#######################################################################################
my @curoptmeta = @metastartvalues;
my ($finished, @testlist, @testlisttargets, $opttarget, $optvalue);
my @snsp = evalsnsp(@curoptmeta);
print "\n";
my $target = sprintf("%.4f", gettarget(@snsp));
$opttarget = $target;
print "starting accuracy: ". join(", ", @snsp) . ", starting target: $target\n";
my $found_improvement;
my @bindings;

#######################################################################################
# optimization loop for meta parameters
#######################################################################################

if ($cmdpars{'opt_trans_matrix'} eq ''){
    my ($a,$b);
    my (@testmeta);
    for (my $r=0; $r<$cmdpars{'rounds'}; $r++) {
	$found_improvement = 0;
	for (my $idx=0; $idx<= $#metaparnames; $idx++) {
	    print "improving parameter $metaparnames[$idx] curently set to $curoptmeta[$idx]\n";
	    @testmeta = @curoptmeta;
	    # set the initial min and max of the range to test
	    if ($metaparranges[$idx][0]  ne 'list') {
		$a = $metaparranges[$idx][1];
		$b = $metaparranges[$idx][2];
		print "$a-$b\n";
	    }
	    $finished = 0;
	    while (!$finished) {
		$finished = 1;
		# generate a list of values to test
		@testlist = ();
		if ($metaparranges[$idx][0] eq 'list') {
		    @testlist = @{$metaparranges[$idx]};
		    shift @testlist;
		} elsif ($metaparranges[$idx][0] eq 'floatrange') {
		    for (my $n=0; $n < $cmdpars{'pstep'}; $n++) {
			push @testlist, $a+$n*($b-$a)/($cmdpars{'pstep'}-1);
		    }
		} else{ # round the values
		    for (my $n=0; $n < $cmdpars{'pstep'}; $n++) {
			my $tv = int($a+$n*($b-$a)/($cmdpars{'pstep'}-1));
			if ($tv ne $testlist[$#testlist]) {
			    push @testlist, $tv;
			}
		    }
		}
		@testlisttargets = ();
		print "$metaparnames[$idx]: checking values " . join("\t", @testlist) . "\n";
		foreach my $testvalue (@testlist){
		    $testmeta[$idx] = $testvalue; # set the parameter to the testvalue
		    @snsp = evalsnsp(@testmeta);
		    print "\n";
		    $target = sprintf("%.4f", gettarget(@snsp));
		    push @testlisttargets , $target;
		    if ($target > $opttarget){ # found improvement
			$optvalue = $testvalue;
			@curoptmeta = @testmeta;
			$opttarget = $target;
			$found_improvement = 1;
			print "found improvement: ". join(", ", @snsp) . ", optimal target: $target\n";
			print "changing $metaparnames[$idx] to $optvalue\n";
			printmetavalues(\@metaparnames, \@testmeta);
			savenewpars(@curoptmeta);
			$finished = 0;
		    }
		}
		print "values  " . join("\t", @testlist) . "\n";
		print "targets " . join("\t", @testlisttargets) . "\n";
		if ($finished == 0){
		    # determine whether further improvements are possible at all
		    # and compute the new range boundaries
		    if ($metaparranges[$idx][0] eq 'list') {
			$finished = 1;
		    } else {
			my ($newa, $newb);
			$newa = $optvalue - ($b-$a)/($cmdpars{'pstep'}-1);
			$newb = $optvalue + ($b-$a)/($cmdpars{'pstep'}-1);
			$newa = ($newa < $a)? $a : $newa;
			$newb = ($newb > $b)? $b : $newb;
			$a = $newa;
			$b = $newb;
			if ($metaparranges[$idx][0] eq 'intrange'){
			    $a = int($a+1);	
			    $b = int($b);
			    if ($b<$a) {
				$finished = 1;
			    }
			}
		    }
		}
	    }
	}
#	if (!$found_improvement && $r<$cmdpars{'rounds'}-1) {
#	    print "Could not further improve. Skipping last ". ($cmdpars{'rounds'}-$r-1) ." rounds\n";
#	    last;
#	}
    }
} else {
#######################################################################################
# optimization loop for transition probabilities
#######################################################################################
    my (@trylist, @normedlist);
    if ($cmdpars{'matrix_constraints'} ne ''){
	@trylist = getStateList("TRY");
	@normedlist = getStateList("NORMED");
    } else {
	@trylist = (0 .. $n);       # optimize all transitions
	@normedlist = (0 .. $n);    # normalize all transitions
    }
    print "Try list: " . (join " ", @trylist) . "\n";
    print "Normed list: " . (join " ", @normedlist) . "\n";
    getBindings();

    print "Optimizing transitions from these states.\n";
    my @curopttrans = @trans;
    save_trans_matrix(\@curopttrans, $cmdpars{'opt_trans_matrix'} . '.curopt');
    my @testtrans;
    for (my $r=0; $r<$cmdpars{'rounds'}; $r++) {
	print "Improvement round/cycle " . ($r+1) . "\n";
	$found_improvement = 0;
	foreach my $idx (@trylist){
	    my $normed=0;
	    if (grep /^$idx$/, @normedlist){
		$normed=1;
	    }
	    my $normsum;
	    if ($normed) {
		$normsum=0;
		for(my $j=0; $j<$n; $j++) {
		    $normsum += $curopttrans[$idx][$j];
		}
	    }		
	    my @tolist;
	    my @transvec;
	    for (my $j=0; $j<$n; $j++) {
		if ($curopttrans[$idx][$j]>0){
		    push @tolist, $j;
		    push @transvec, $curopttrans[$idx][$j];
		}
	    }
	    # skip state if it is normed and just one transition out of this state is possible.
	    next unless (!$normed || @tolist > 1);
	    
	    print "Improving transitions out of state $idx\n";
	    print "Nonzero transitions from state $idx: ". join (" ", grep ($_>0, @transvec)) . "\n";
	    
	    # make a list with all the varied probability vectors to try
	    my @tryvectors = getVariedTransVectors(\@transvec, $normsum, $normed);
	    print "Trying " . scalar(@tryvectors) . " " .
		($normed? "normed":"unnormed") . " variations of transition vector.\n";
	    # change each transition probability
	    # evaluate the accuracy
	    foreach my $varieddist (@tryvectors) {
		print "Try varied distribution " . join (" ", @{$varieddist}) . " ";
		# create varied transition matrix
		copyMatrix(\@testtrans,\@curopttrans, $n);
		for (my $k=0; $k < @tolist; $k++) {
		    $testtrans[$idx][$tolist[$k]] = $varieddist->[$k];
		}
		realizeBindings($idx, \@testtrans, 0);
		# save it to the file
		save_trans_matrix(\@testtrans, $cmdpars{'opt_trans_matrix'});
		# start an evaluation run
		@snsp = evalsnsp();
		$target = sprintf("%.4f", gettarget(@snsp));
		print "\ttarget=$target\n";
		if ($target > $opttarget){ # found improvement
		    $opttarget = $target;
		    $found_improvement = 1;
		    print "*** Found improvement: ". join(", ", @snsp) . ", optimal target: $target\n";
		    print "changing trans. probs out of state $idx from " 
			. (join " ", (grep ($_ > 0, @{$curopttrans[$idx]})))
			. " to " . (join " ", (grep ($_>0, @{$testtrans[$idx]}))) . "\n";
		    copyMatrix(\@curopttrans, \@testtrans, $n);
		    save_trans_matrix(\@curopttrans, $cmdpars{'opt_trans_matrix'} . ".curopt");
		    $finished = 0;
		} else { # no improvement
		    # restore file with currently optimal transition matrix
		    save_trans_matrix(\@curopttrans, $cmdpars{'opt_trans_matrix'});
		}
	    }
	}
	if (!$found_improvement && $r<$cmdpars{'rounds'}-1) {
	    print "Could not further improve. Skipping last ". ($cmdpars{'rounds'}-$r-1) ." rounds\n";
	    last;
	}
    }
}

#######################################################################################
# final training
#######################################################################################

if ($cmdpars{'noTrainPars'} eq ''){
    # delete the temporary *pbl files
    for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
	system ("rm -f $speciesdir/exon-tmp$k.pbl $speciesdir/intron-tmp$k.pbl $speciesdir/igenic-tmp$k.pbl $speciesdir/utr-tmp$k.pbl");
    }
    # delete the training files for cross-validation (can be large)
    if ($cmdpars{'cpus'} > 1){
	for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
	    # delete the temporary training files
	    system("rm -f $optdir/curtrain-$k $optdir/predictions-$k.txt");
	}
    }
    print "Making final training with the optimized parameters.\n";
    # make the joint training file (train.gb and onlytrain.gb)
    system("rm -f $optdir/curtrain $optdir/curtest");
    system("cp $cmdpars{'train.gb'} $optdir/curtrain");
    if ($cmdpars{'onlytrain'} ne ''){
	system ("cat $cmdpars{'onlytrain'} >> $optdir/curtrain");
    }
    my $cmd = "etraining --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $be_silent $optdir/curtrain $pars $modelrestrict";
    print "$cmd\n";
    system($cmd);
    system("rm -f $optdir/curtrain");
}


#######################################################################################
# subroutines
#######################################################################################

################################################
# evalsnsp: determine the values
# base sn, base sp, exon sn, exon sp, gene sn, gene sp, tss medianDiff, tts medianDiff
# sn: sensitivity, sp: specificity
# given a set of metaparameter values
################################################

my %storedsnsp = {}; # hash with the stored sn and sp array references

sub evalsnsp {
    my @values = @_;
    my ($cbsn,$cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd); # accuracy values of current bucket
    my ($gbsn,$gbsp, $gesn, $gesp, $ggsn, $ggsp, $gsmd, $gtmd); # total accuracy values
    $gbsn=$gbsp=$gesn=$gesp=$ggsn=$ggsp=$gsmd=$gtmd=0;
    my $argument = ""; 
    if ($cmdpars{'opt_trans_matrix'} eq '') {
	# make the parameters string for the command line
	for (my $i=0; $i<= $#metaparnames; $i++) {
	    $argument=$argument." --".$metaparnames[$i]."=".$values[$i];
	}
	# print "argument:$argument\n";
	# check if accuracy has already been computed for this parameter combination
	if (exists($storedsnsp{$argument})){
	    print "retrieve accuracy from previous computation";
	    return @{$storedsnsp{$argument}};
	}
    }
    if(!$cgp){
	# Loop over the buckets and chose bucket k as the one for testing.
	# All other buckets are taken for training if appropriate.
	if ($cmdpars{'cpus'} <= 1) {
	    print "bucket ";
	    for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
		# make the temporary training and testing files
		system("rm -f $optdir/curtrain $optdir/curtest");
		system("cp $optdir/bucket$k.gb $optdir/curtest");
		for (my $m=1; $m<=$cmdpars{"kfold"}; $m++) {
		    if ($m != $k) {
			system("cat $optdir/bucket$m.gb >> $optdir/curtrain");
		    }
		}
		if ($cmdpars{'onlytrain'} ne ''){
		    system ("cat $cmdpars{'onlytrain'} >> $optdir/curtrain");
		}
		if ($cmdpars{'noTrainPars'} eq '') {# no need to retrain if the trans matrix is optimized or this option is otherwise explicitly set.
		    system("$cmdpars{'aug_exec_dir'}etraining --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $argument $pars $be_silent $modelrestrict $optdir/curtrain");
		}
		
		system("$cmdpars{'aug_exec_dir'}augustus --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $argument $pars $optdir/curtest > $optdir/predictions.txt");
		
		open (PRED, "<$optdir/predictions.txt");
		while (<PRED>){
		    if(/nucleotide level \| +(\S+) \| +(\S+) \|$/){
			($cbsn, $cbsp) = ($1,$2);
		    }
		    if (/exon level \|.*-- \| +(\S+) \| +(\S+) \|$/){
			($cesn, $cesp) = ($1,$2);
		    }
		    if (/gene level \|.* \| +(\S+) \| +(\S+) \|$/){
			($cgsn, $cgsp) = ($1,$2);
		    }
		    if (/TSS \|.* \|.* \|.* \| +(\S+) \|$/){
			$csmd = $1;
		    }
		    if (/TTS \|.* \|.* \|.* \| +(\S+) \|$/){
			$ctmd = $1;
		    }
		}
		close(PRED);
		
		
		if ($cbsn eq "" || $cbsp eq "" ||$cesn eq "" ||$cesp eq "" ||$cgsn eq "" ||$cgsp eq "") {
		    die ("Could not read the accuracy values out of predictions.txt when processing bucket $k.");
		}
		print "$k ";
		#print "accuracy on bucket$k: $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd\n";
		$gbsn += $cbsn;
		$gbsp += $cbsp;
		$gesn += $cesn;
		$gesp += $cesp;
		$ggsn += $cgsn;
		$ggsp += $cgsp;
		$gsmd += $csmd; 
		$gtmd += $ctmd; 
	    }
	} else { # parallel
	    print "bucket ";
	    my $pm = new Parallel::ForkManager($cmdpars{'cpus'});
	    for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
		# fork and return the pid for the child:
		my $pid = $pm->start and next;
		# this part is done in parallel by the child process
		
		my $pbloutfiles = "--/ExonModel/outfile=exon-tmp$k.pbl --/IntronModel/outfile=intron-tmp$k.pbl --/IGenicModel/outfile=igenic-tmp$k.pbl --/UtrModel/outfile=utr-tmp$k.pbl";
		my $pblinfiles = "--/ExonModel/infile=exon-tmp$k.pbl --/IntronModel/infile=intron-tmp$k.pbl --/IGenicModel/infile=igenic-tmp$k.pbl --/UtrModel/infile=utr-tmp$k.pbl";
		if (defined($cmdpars{'trainOnlyUtr'}) && $cmdpars{'trainOnlyUtr'} ne "" && $cmdpars{'trainOnlyUtr'} ne "0"){
		    $pblinfiles = "--/UtrModel/infile=utr-tmp$k.pbl";
		}
		
		if ($cmdpars{'noTrainPars'} eq '') {# No need to retrain if the trans matrix is optimized or noTrainPars=1 set explicitly.
		    my $cmd = "$cmdpars{'aug_exec_dir'}etraining --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $argument $pars "
			. "$be_silent $modelrestrict $pbloutfiles $optdir/curtrain-$k";
		    system($cmd);
#		unlink $optdir/curtrain-$k;
		} else {
		    $pblinfiles = ""; # training did not take place, so the $pbloutfiles have not beeen created and cannot be used for prediction
		}
		
		system("$cmdpars{'aug_exec_dir'}augustus --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $argument $pars $pblinfiles $optdir/bucket$k.gb > $optdir/predictions-$k.txt");
		print "$k ";
		
		$pm->finish; # terminate the child process
	    }
	    $pm->wait_all_children;
	    # compute accuracy
	    for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
		open (PRED, "<$optdir/predictions-$k.txt");
		while (<PRED>){
		    if(/nucleotide level \| +(\S+) \| +(\S+) \|$/){
			($cbsn, $cbsp) = ($1,$2);
		    }
		    if (/exon level \|.*-- \| +(\S+) \| +(\S+) \|$/){
			($cesn, $cesp) = ($1,$2);
		    }
		    if (/gene level \|.* \| +(\S+) \| +(\S+) \|$/){
			($cgsn, $cgsp) = ($1,$2);
		    }
		    if (/TSS \|.* \|.* \|.* \| +(\S+) \|$/){
			$csmd = $1;
		    }
		    if (/TTS \|.* \|.* \|.* \| +(\S+) \|$/){
			$ctmd = $1;
		    }
		}
		close(PRED);
		
		if ($cbsn eq "" || $cbsp eq "" ||$cesn eq "" ||$cesp eq "" ||$cgsn eq "" ||$cgsp eq "") {
		    die ("Could not read the accuracy values out of predictions.txt when processing bucket $k.");
		}
		#print "accuracy on bucket$k: $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd\n";
		$gbsn += $cbsn;
		$gbsp += $cbsp;
		$gesn += $cesn;
		$gesp += $cesp;
		$ggsn += $cgsn;
		$ggsp += $cgsp;
		$gsmd += $csmd;
		$gtmd += $ctmd;
	    }
	} 	
    #print "\n";
    $gbsn = sprintf("%.4f", $gbsn/$cmdpars{'kfold'});
    $gbsp = sprintf("%.4f", $gbsp/$cmdpars{'kfold'});
    $gesn = sprintf("%.4f", $gesn/$cmdpars{'kfold'});
    $gesp = sprintf("%.4f", $gesp/$cmdpars{'kfold'});
    $ggsn = sprintf("%.4f", $ggsn/$cmdpars{'kfold'});
    $ggsp = sprintf("%.4f", $ggsp/$cmdpars{'kfold'});
    $gsmd = sprintf("%.2f", $gsmd/$cmdpars{'kfold'});
    $gtmd = sprintf("%.2f", $gtmd/$cmdpars{'kfold'});    

    } else { # cgp parameter optimization
	print "maf ";
	my $pm = new Parallel::ForkManager($cmdpars{'cpus'});
	for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {
	    # fork and return the pid for the child:
	    my $pid = $pm->start and next;
	    # this part is done in parallel by the child process
	    system("rm -rf $optdir/pred$k; mkdir $optdir/pred$k");
	    system("$cmdpars{'aug_exec_dir'}augustus --species=$cmdpars{'species'} --AUGUSTUS_CONFIG_PATH=$configdir $argument $pars --alnfile=$optdir/splitMaf/$k.maf --/CompPred/outdir=$optdir/pred$k >$optdir/pred$k/aug.out");
	    print "$k ";
	    $pm->finish; # terminate the child process
	}
	$pm->wait_all_children;
	# calcuate accuracy values with eval package
	($gbsn,$gbsp, $gesn, $gesp, $ggsn, $ggsp, $gsmd, $gtmd) = evalCGP();
    }
        
    my @returnarray =($gbsn, $gbsp, $gesn, $gesp, $ggsn, $ggsp, $gsmd, $gtmd);
    $storedsnsp{$argument}= \@returnarray;
    
    return @returnarray;
}

######################################################################################
# gettarget: get an optimization target value from
# base sn, base sp, exon sn, exon sp, gene sn, gene sp, tss medianDiff, tts medianDiff
# feel free to change the weights
######################################################################################

sub gettarget {
    my ($bsn, $bsp, $esn, $esp, $gsn, $gsp, $smd, $tmd) = @_;
    my $f = $cmdpars{'sens_spec_bias'};
    return (3*$f*$bsn + (3/$f)*$bsp + 4*$f*$esn + (4/$f)*$esp + 2*$f*$gsn + (2/$f)*$gsp + 40/($smd+40) + 40/($tmd+40))/(9*($f+(1/$f))+2);
#   return (3*$bsn + 9*$bsp + 4*$esn + 12*$esp + 2*$gsn + 6*$gsp)/36;
}

################################################
# savenewpars: replace parameters in the file $species_cfg_filename
################################################
sub savenewpars{ 
    open (SPCCFG, ">$species_cfg_filename") or die ("Could not open $species_cfg_filename");
    print "Writing new parameters to $species_cfg_filename...\n";
    my $parname;
    foreach my $line (@spcfilelines) {
	if ($line =~ /^\s*\#.*/ || $line =~ /^\s*$/) {
	    print SPCCFG $line; # print unchanged line
	    next;
	}
	# format:
	# parname   value   # comment
	# or
	# parname   value
	if ($line =~ /^(\s*)(\S+)(\s+)(\S+)(.*)$/) {
	    $parname = $2;
	    my $index = idx(\@metaparnames, $parname);
	    if($index != -1){
		print SPCCFG $1 . $parname . $3 . $_[$index] . $5 . "\n";
	    } else {
		print SPCCFG $line;
	    }
	} else {
	    print SPCCFG $line;
	}
    }
    close (SPCCFG);
}

################################################
# save_trans_matrix: replace existing transition 
# probabilities with the given new ones
################################################
sub save_trans_matrix {
    my $newtransref = shift;
    my $filename = shift;
    open (TRANS, ">$filename") or die ("Could not open $filename for writing.");
    foreach my $line (@transfilelines) {
	if ($line =~ /^(\s*)(\d+)(\s+)(\d+)(\s+)(\S+)/) {
	    my $value = 0;
	    $value = $newtransref->[$2][$4] unless (! defined $newtransref->[$2][$4]);
	    print TRANS $1 . $2 . $3 . $4 . $5 . $value . "\n";
	} else {
	    print TRANS $line;
	}
    }
    close TRANS;
}


sub printmetavalues {
    my ($names, $values) = @_;
    for (my $i=0; $i<= $#$names; $i++){
	print $names->[$i] . "\t" . $values->[$i] . "\n";
    }
}

sub printmetaranges {
    my ($names, $ranges) = @_;
    print "metapar ranges $#$names:\n";
    for (my $i=0; $i <= $#$names; $i++){
	print $names->[$i] . "\t" . join(" ", @{$ranges->[$i]}) . "\n";
    }
}

#
# idx: find the index to an element in an array of strings
#

sub idx {
    my $arrayref = shift;
    my $element = shift;

    for (my $i=0; $i<=$#$arrayref; $i++) {
	if ($arrayref->[$i] eq $element) {
	    return $i;
	}
    }
    return -1;
}

#
# norm a vector to sum up to $normsum
# if $normed is true
sub norm {
    my ($vecref, $normsum, $normed) = @_;
    if ($normed) {
	my $sum=0;
	foreach my $item (@{$vecref}){
	    $sum += $item;
	}
	if ($sum > 0 && $normsum > 0) {
	    foreach my $item (@{$vecref}){
		$item *= $normsum/$sum;
	    }
	}
    }
    # round to 6 places
    foreach my $item (@{$vecref}){
	$item = sprintf("%.6f", $item);
	$item =~ s/(\.\d*[1-9])0+$/$1/; # remove trailing zeros
    }
}

#
# read try list or normed list from matrix_constraints file
#
sub getStateList {
    my $identifyer = shift; # 'TRY' or 'NORMED'
    if (!(-s $cmdpars{'matrix_constraints'})) {
	print "Could not find the file with the transition matrix constraints: $cmdpars{'matrix_constraints'}\n";
	return;
    }
    my @statelist;
    open CONSTR, "<$cmdpars{'matrix_constraints'}" or die;
    my $scanning=0;
    my $all=0;
    while(<CONSTR>){
	if (/\[$identifyer\]/) {
	    $scanning=1;
	} elsif (/\s*\[.*\]/) {
	    $scanning=0;
	}
	if ($scanning && /^\s*(\d+)\s*/) {
	    push @statelist, $1;
	}
	if ($scanning && /^\s*all/) {
	    $all=1;
	}
    }
    close CONSTR;
    if ($all) {
	@statelist = (0..$n);
    }
    return @statelist;
}

#
# read bindings list from matrix_constraints file
# into global variables
# 
#
sub getBindings {
   if (!(-s $cmdpars{'matrix_constraints'})) {
	print "Could not find the file with the transition matrix constraints: $cmdpars{'matrix_constraints'}\n";
	return;
    }
   my @bindingsstr;
   open BINDINGS, "<$cmdpars{'matrix_constraints'}" or die;
   my $scanning=0;
   while(<BINDINGS>){
       if (/\[BINDINGS\]/) {
	   $scanning=1;
       } elsif (/\s*\[.*\]/) {
	   $scanning=0;
       }
       if ($scanning){
	   $_ =~ s/#.*//;
	   if (/^\s*(\(.*\))/ || /^(MC\S*)/) {
	       push @bindingsstr, $1;
	   }
       }
   }
   close BINDINGS;
   print "bindings: " . (join "\n", @bindingsstr) . "\n";
   my $btype;

   # @bindings is a list of references of bindings. 
   # Each binding is a list of
   # - the binding equation string
   # - a type, (either '+' or '/')
   # - leftTrans
   # - rightTrans
   # - allStates, the list of states the transitions originate in
   # leftTrans and rightTrans are the transition lists on the left-hand side and right-hand side of the equation
   # They each are lists of transitions.
   # lt1 + lt2 + ... + ltA = rt1 + rt2 + ... + rtB
   # A and B can be 1.
   # example (0,24)+(0,25)=(0,65)+(0,70)
   # 
   foreach my $bstr (@bindingsstr) {
       if ($bstr =~ '/') {
	   $btype = '/';
       } elsif ($bstr =~ /^MC/) {
	   $btype = 'M';
       } else {
	   $btype = '+';
       }
       my @eqsides = split /=/, $bstr;
       if (@eqsides != 2) {
	   print "Error: Wrong format in $cmdpars{'matrix_constraints'}. Each binding must be an equation.\n";
	   return;
       }
       my @leftTrans = parseTransList($eqsides[0]);  # left hand side of equation
       my @rightTrans = parseTransList($eqsides[1]); # right hand side of equation
       my @allStates=();
       my %seen = ();
       foreach my $trans ((@leftTrans,@rightTrans)){
	   push @allStates, $trans->[0] unless $seen{$trans->[0]}++; # push the originating state
       }
       push @bindings, [$bstr, $btype, \@leftTrans, \@rightTrans, \@allStates];
   }
}

#
# realizeBindings
# 
# returns error message string if there are any
sub realizeBindings {
    my $state = shift;                # line that was just changed
    my $trans = shift;                # transition matrix that may need to be adjusted
    my $doNothingButComplain = shift; # only error msgs
    my $errmsg;
    print "realizeBindings $state\n";
    foreach my $binding (@bindings) {
	my ($bstr, $btype, $leftTrans, $rightTrans, $allStates) = @{$binding};
	if (grep /^$state$/, @{$allStates}) {
	    print "binding applies: $bstr\n";
	    if ($btype eq '+') {
		if (@{$allStates} == 1) {
		    # binding applies just to $state
		    # lt1 + lt2 + ... + ltA   rt1 + rt2 + ... + rtB
		    # ---------------------   ---------------------
		    #       lhs                      rhs
		    my $lhs=0;
		    my $rhs=0;
		    foreach my $t (@{$leftTrans}){
			$lhs += $trans->[$state][$t->[1]];
		    }
		    foreach my $t (@{$rightTrans}){
			$rhs += $trans->[$state][$t->[1]];
		    }
		    #print "actual value lhs=$lhs, rhs=$rhs\n";
		    if ($lhs != $rhs) {
			if ($doNothingButComplain) {
			    $errmsg .= "Binding $bstr not satisfied: $lhs != $rhs\n";
			} else {
			    # rescale the probabilities
			    if ($lhs>0 && $rhs>0){
				foreach my $t (@{$leftTrans}){
				    $trans->[$state][$t->[1]] *= ($lhs+$rhs)/2/$lhs;
				}
				foreach my $t (@{$rightTrans}){
				    $trans->[$state][$t->[1]] *= ($lhs+$rhs)/2/$rhs;
				}
				roundVector(@{$trans->[$state]});
				print "After applying binding: " . (join " " , grep {$_!=0} @{$trans->[$state]}) . "\n";
			    }
			}
		    }
		} else {
		    # binding applies to other states as well
		    # change only the other states and not the transition probabilities from this state
		    #
		    my $restsum; #sum of of all transition probabilities except the ones from $state. Those are fixed.
		    foreach my $t (@{$rightTrans}){
			if ($t->[0]==$state){
			    $restsum += $trans->[$t->[0]][$t->[1]];
			}
		    }
		    foreach my $t (@{$leftTrans}){
			if ($t->[0]==$state){
			    $restsum -= $trans->[$t->[0]][$t->[1]];
			}
		    }
		    # compute the actual restsum
		    my $actualrestsum;
		    foreach my $t (@{$rightTrans}){
			if ($t->[0]!=$state){
			    $actualrestsum -= $trans->[$t->[0]][$t->[1]];
			}
		    }
		    foreach my $t (@{$leftTrans}){
			if ($t->[0]!=$state){
			    $actualrestsum += $trans->[$t->[0]][$t->[1]];
			}
		    }
		    if ($restsum*$actualrestsum>0) {
			# rescale the other states' transition probabilities
			my %added = ();
			foreach my $t (@{$rightTrans}, @{$leftTrans}){
			    if ($t->[0]!=$state){
				$added{$t->[0]} += $trans->[$t->[0]][$t->[1]] * ($restsum/$actualrestsum - 1);
				$trans->[$t->[0]][$t->[1]] *= $restsum/$actualrestsum;
			    }
			}
			# now adjust the other states transition vectors
			foreach my $ostate (@{$allStates}){
			    if ($ostate != $state) {
				print "trans probs of other state $ostate: " . 
				    (join " " , grep {$_!=0} @{$trans->[$ostate]}) . "\n";
				# rescale all other transition probabilities, such 
				# that their sum is decreased by $added{$ostate}
				# print "sum must be decreased by $added{$ostate}\n";
				my %otochanged = ();
				foreach my $t (@{$rightTrans}, @{$leftTrans}){
				    if ($t->[0]==$ostate){
					$otochanged{$t->[1]} = 1;
				    }
				}
				my $osum;
				for (my $j=0; $j<$n; $j++) {
				    if (!$otochanged{$j}) {
					$osum += $trans->[$ostate][$j];
				    }
				}
				# now renorm so that the other trans probs sum up to $osum-$added{$ostate}
				if ($osum>0 && $osum-$added{$ostate}>0){
				    for (my $j=0; $j<$n; $j++) {
					if (!$otochanged{$j}) {
					    if ($trans->[$ostate][$j]>0) {
						print "($ostate,$j) is changed from $trans->[$ostate][$j] to ";
						$trans->[$ostate][$j] *= ($osum-$added{$ostate})/$osum;
						print "$trans->[$ostate][$j]\n";
					    }
					}
				    }
				    
				}
				roundVector(@{$trans->[$ostate]});
				print "fixed trans probs of other state $ostate: " . 
				    (join " " , grep {$_!=0} @{$trans->[$ostate]}) . "\n";
			    }
			}
		    } else {
			# Don't fix distribution matrix. Should not happen.
		    }
		}

	    } elsif ($btype eq '/') {
		# not yet implemented
	    } elsif ($btype eq 'M') {
		# markov chain binding
                # make sure that one 3x3 Markov chain transition matrix is from the time-reversed
                # process of another Markov chain. Assume stationarity.
		
		my $leftm = [[$trans->[$leftTrans->[0]->[0]][$leftTrans->[0]->[1]],
			      $trans->[$leftTrans->[1]->[0]][$leftTrans->[1]->[1]],
			      $trans->[$leftTrans->[2]->[0]][$leftTrans->[2]->[1]]],
			     [$trans->[$leftTrans->[3]->[0]][$leftTrans->[3]->[1]],
			      $trans->[$leftTrans->[4]->[0]][$leftTrans->[4]->[1]],
			      $trans->[$leftTrans->[5]->[0]][$leftTrans->[5]->[1]]],
			     [$trans->[$leftTrans->[6]->[0]][$leftTrans->[6]->[1]],
			      $trans->[$leftTrans->[7]->[0]][$leftTrans->[7]->[1]],
			      $trans->[$leftTrans->[8]->[0]][$leftTrans->[8]->[1]]]];
		my $rightm = [[$trans->[$rightTrans->[0]->[0]][$rightTrans->[0]->[1]],
			       $trans->[$rightTrans->[1]->[0]][$rightTrans->[1]->[1]],
			       $trans->[$rightTrans->[2]->[0]][$rightTrans->[2]->[1]]],
			      [$trans->[$rightTrans->[3]->[0]][$rightTrans->[3]->[1]],
			       $trans->[$rightTrans->[4]->[0]][$rightTrans->[4]->[1]],
			       $trans->[$rightTrans->[5]->[0]][$rightTrans->[5]->[1]]],
			      [$trans->[$rightTrans->[6]->[0]][$rightTrans->[6]->[1]],
			       $trans->[$rightTrans->[7]->[0]][$rightTrans->[7]->[1]],
			       $trans->[$rightTrans->[8]->[0]][$rightTrans->[8]->[1]]]];
		
		my $revLeft=0; # 1 iff the left Markov chain needs to be adjusted
		for (my $i=0; $i<9; $i++) {
		    if ($rightTrans->[$i]->[0] == $state){
			$revLeft=1;
		    }
		}
		my ($m, $rm);
		if ($revLeft) {
		    $m = $rightm;
		    $rm = $leftm;
		} else {
		    $m = $leftm;
		    $rm = $rightm;
		}

		reverseMatrix($m, $rm);
		if ($revLeft) {
		    $trans->[$leftTrans->[0]->[0]][$leftTrans->[0]->[1]] = $rm->[0][0];
		    $trans->[$leftTrans->[1]->[0]][$leftTrans->[1]->[1]] = $rm->[0][1];
		    $trans->[$leftTrans->[2]->[0]][$leftTrans->[2]->[1]] = $rm->[0][2];
		    $trans->[$leftTrans->[3]->[0]][$leftTrans->[3]->[1]] = $rm->[1][0];
		    $trans->[$leftTrans->[4]->[0]][$leftTrans->[4]->[1]] = $rm->[1][1];
		    $trans->[$leftTrans->[5]->[0]][$leftTrans->[5]->[1]] = $rm->[1][2];
		    $trans->[$leftTrans->[6]->[0]][$leftTrans->[6]->[1]] = $rm->[2][0];
		    $trans->[$leftTrans->[7]->[0]][$leftTrans->[7]->[1]] = $rm->[2][1];
		    $trans->[$leftTrans->[8]->[0]][$leftTrans->[8]->[1]] = $rm->[2][2];
		} else {
		    $trans->[$rightTrans->[0]->[0]][$rightTrans->[0]->[1]] = $rm->[0][0];
		    $trans->[$rightTrans->[1]->[0]][$rightTrans->[1]->[1]] = $rm->[0][1];
		    $trans->[$rightTrans->[2]->[0]][$rightTrans->[2]->[1]] = $rm->[0][2];
		    $trans->[$rightTrans->[3]->[0]][$rightTrans->[3]->[1]] = $rm->[1][0];
		    $trans->[$rightTrans->[4]->[0]][$rightTrans->[4]->[1]] = $rm->[1][1];
		    $trans->[$rightTrans->[5]->[0]][$rightTrans->[5]->[1]] = $rm->[1][2];
		    $trans->[$rightTrans->[6]->[0]][$rightTrans->[6]->[1]] = $rm->[2][0];
		    $trans->[$rightTrans->[7]->[0]][$rightTrans->[7]->[1]] = $rm->[2][1];
		    $trans->[$rightTrans->[8]->[0]][$rightTrans->[8]->[1]] = $rm->[2][2];
		}
	    }
	}
    }
}

#
# parseTransList
# parses a string like (0,24)+(0,25) or (26,28)/(26,29)
# into a list of pairs
sub parseTransList{
    my $tstr = shift;
    my @transList = ();
    foreach my $tok (split /\)[\(\)\+\*,]*\(/, $tstr) {
	$tok =~ /(\d+),\s*(\d+)/;
	push @transList, [$1,$2];
    }
    return @transList;
}


#
# getVariedTransVectors
# parameters: \@transvec, $normsum, $normed
# returns a list of references to varied transition probability vectors
#
sub getVariedTransVectors {
    my $transvec = shift;
    my $normsum = shift;
    my $normed = shift;
    my @tryvectors = (); #list of references to transition probability vectors
    my $end; # vary components up to this index
    my $maxfactor = 2; # factor is between 1 and this
    my $variations = 3;
    if (!$normed || @{$transvec} > 2) {
	$end = @{$transvec};
    } else {
	$end=1;
    }
    # simply enlarge and decrease each single transition for now
    for(my $k=0; $k < $end; $k++) {
	for (my $r=0; $r<$variations; $r++){
	    # change element k
	    my @vartransvec = @{$transvec};
	    my $factor = 1+rand($maxfactor-1.0);
	    if (rand(2)<1) {
		$factor = 1/$factor;
	    }
	    $vartransvec[$k] = $transvec->[$k] * $factor;
	    norm(\@vartransvec, $normsum, $normed);
	    push @tryvectors, \@vartransvec;
	}
    }
    return @tryvectors;
}


#
# copyMatrix
# copy a tranistion matrix
#
sub copyMatrix {
    my $new = shift;
    my $old = shift;
    my $size = shift;
    for (my $i=0; $i<$size; $i++) {
	$new->[$i]=[];
	for (my $j=0; $j<$size; $j++) {
	    push @{$new->[$i]}, $old->[$i][$j];
	}
    }
}


#
# roundVector
#
sub roundVector {
    my $sum=0;
    foreach my $v (@_) {
	$sum += $v;
	$v = sprintf("%.6f", $v);
    }
    my $roundedsum = sprintf("%.6f", $sum); # is usually 1
    my $newsum=0;
    my $maxel=0;
    for (my $i=0; $i < @_; $i++) {
	$newsum += $_[$i];
	if ($_[$i]> $_[$maxel]) {
	    $maxel = $i;
	}
    }
    # in case the rounding changed the sum a little, adjust the largest element only
    if ($newsum != $roundedsum){
	@_[$maxel] += $roundedsum-$newsum;
    }
}

#
# matrixSquare
#
sub matrixSquare {
    my $m = shift; # reference to quadratic matrix
    
    my $n = @{$m}; # dimension nxn
    my @ret;
    for (my $i=0; $i<$n; $i++) {
	$ret[$i]=[];
	for (my $j=0; $j<$n; $j++) {
	    for (my $k=0; $k<$n; $k++) {
		$ret[$i][$j] += $m->[$i][$k] * $m->[$k][$j];
	    }
	}
    }
    return \@ret;
}

#
# printMatrix
#

sub printMatrix {
    my $m = shift;
    my $n = @{$m}; # dimension nxn
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$n; $j++) {
	    print $m->[$i][$j] . "\t";
	}
	print "\n";
    }
}

#
# reverseMatrix
#
#
sub reverseMatrix {
    my $m = shift;
    my $rm = shift;
		
    # first scale both matrices so they are real transition matrices
    # and remember the scaling factor
    # 
    my (@scale, @scaler);
    for (my $i=0; $i<@{$m}; $i++) {	
	for (my $j=0; $j<@{$m->[$i]}; $j++){
	    $scale[$i] += $m->[$i][$j];
	}
	for (my $j=0; $j<@{$m->[$i]}; $j++){
	    $m->[$i][$j] /= $scale[$i];
	}
    }
    for (my $i=0; $i<@{$rm}; $i++) {	
	for (my $j=0; $j<@{$rm->[$i]}; $j++){
	    $scaler[$i] += $rm->[$i][$j];
	}
	for (my $j=0; $j<@{$rm->[$i]}; $j++){
	    $rm->[$i][$j] /= $scaler[$i];
	}
    }
    #print "scale vector " . join (" ", @scale) . "\n";
    #print "reverse scale vector " . join (" ", @scaler) . "\n";
    
    # find stationary distribution 
    # approximately by taking a sufficiently large
    # matrix power of the transition matrix

    #print "m=\n";
    #printMatrix($m);
    #print "rm=\n";
    #printMatrix($rm);
    my $pm;
    copyMatrix(\@{$pm}, $m, scalar(@{$m}));
    for (my $i=0; $i<5; $i++) {
	$pm=matrixSquare($pm);
    }
    #print "power matrix\n";
    #printMatrix($pm);
    my @pi = @{$pm->[0]}; # stationary distribution
    #print "stationary distribution: " . join(" ", @pi) . "\n";
    # rtime-everse other matrix rm
    # according to the formula
    # Q_ij = pi_j * P_ji / pi_i  
    # Q: transition matrix of time-reversed chain
    # P: transition matrix of normal chain
    # pi: stationary distribution
    for (my $i=0; $i<@{$rm}; $i++) {
	for (my $j=0; $j<@{$rm->[$i]}; $j++){
	    $rm->[$i][$j] = $pi[$j] * $m->[$j][$i] / $pi[$i];
	}
	roundVector(@{$rm->[$i]});
    }
    print "reversed chain:\n";
    printMatrix($rm);
    # rescale both chains so the sums are like before
    for (my $i=0; $i<@{$rm}; $i++) {	
	for (my $j=0; $j<@{$rm->[$i]}; $j++){
	    $rm->[$i][$j] *= $scaler[$i];
	}
	for (my $j=0; $j<@{$m->[$i]}; $j++){
	    $m->[$i][$j] *= $scale[$i];
	}	
    }
    #print "rescaled reversed chain:=\n";
    #printMatrix($rm);
}

###########################################################################
# splitMaf:
# splits the alignment into smaller chunks for parallel computing.
# Splitting is done with respect to a reference genome
# (first s-line in a MAF block). The chunk size, that is the length of the
# sequence segment in the reference genome that covers the alignment chunk,
# can be chosen and is approximately constant over all alignment chunks. 
###########################################################################

sub splitMaf{

    if ( $cmdpars{"chunksize"} < 1000000 ) {
	print STDERR "Warning: the chunksize must be at least 1000000bps.\n";
	print STDERR "resuming with chunksize=1000000.\n";
	$cmdpars{"chunksize"} = 1000000;
    }

    $cmdpars{"kfold"} = 1; # number of chunks

    # directory of alignment chunks
    my $mafDir = "$optdir/splitMaf";
    system("rm -rf $mafDir; mkdir $mafDir");

    open (MAF, <$cmdpars{"alnfile"}>) or die ("Could not open $cmdpars{'alnfile'} for reading: $!");
    open (MAFCHUNK, ">$mafDir/$cmdpars{'kfold'}.maf") or die ("Could not open $mafDir/$cmdpars{'kfold'}.maf for writing: $!");
    
    if($cmdpars{'cpus'} <= 1){ # in the serial case just make a temporary copy of the whole alignment
	system("cp $cmdpars{'alnfile'} $mafDir/$cmdpars{'kfold'}.maf");
    } else {
	
	my ($chr, $start, $len); # coordinates of the previous maf block in the reference genome
	my $seqlen = 0;
	my $newblock = 1; # a new block starts
		
	while(<MAF>){
	    if(/^\s*$/){ # end of alignment block
		if($seqlen > $cmdpars{"chunksize"}){ # if seqlen exceeds the max chunk size, start new alignment chunk
		    close(MAFCHUNK);
		    $cmdpars{"kfold"}++;
		    open (MAFCHUNK, ">$mafDir/$cmdpars{'kfold'}.maf") or die ("Could not open $mafDir/$cmdpars{'kfold'}.maf");
		    ($chr, $start, $len);
		    $seqlen = 0;
		}
		else{
		    print MAFCHUNK $_;
		}
		$newblock = 1;
	    }
	    if(/^s\s/){
		if($newblock){ # first-row reference species
		    my @alnRow = split(/\s+/,$_);
		    if(defined($chr) && $alnRow[1] eq $chr){
			$seqlen += ( $alnRow[2] - ($start + $len) + $alnRow[3] );
		    } else {		    
			$seqlen += $alnRow[3];
		    }
		    ($chr, $start, $len) = ($alnRow[1], $alnRow[2], $alnRow[3]);
		}
		print MAFCHUNK $_;
		$newblock = 0;
	    }
	    if(/^a\s/){
		print MAFCHUNK $_;
	    }
	}
	close(MAFCHUNK);		
	close(MAF);
	print "splitting Maf file in $cmdpars{'kfold'} pieces for parallel computing.\n";
    }
}

###################################################################################################
# evalCGP:
# evaluates a prediction in GTF format against an annotation
# using the external evaluation package Eval by Evan Keibler and Michael R. Brent¹
# and returns accuracy values (SN and SP on gene, exon and nucleotide level)
# evalCGP only compares gene features on given genomic intervals that
# are parsed from the prediction files.

# ¹(Eval: A software package for analysis of genome annotations. BMC Bioinformatics 4: 50 (2003))
###################################################################################################

sub evalCGP{

    my ($cbsn,$cbsp, $cesn, $cesp, $cgsn, $cgsp); # accuracy values

    # temporary directory that contains the prediction and the annotation split by seq
    my $gffDir = "$optdir/tempGFF";
    system ("rm -rf $gffDir; mkdir $gffDir");
    
    my @intervals=(); # hash of genomic intervals
    my %seqlist=(); # hash of sequences (only keys, no values)
    
    # join predictions from parallel runs and make gene IDs unique
    open (JOINPRED, ">$optdir/pred.gtf") or die("Could not open $optdir/pred.gtf for writing: $!");
    my $geneID = 0;
   
    for (my $k=1; $k<=$cmdpars{"kfold"}; $k++) {    
	open (PRED, "<$optdir/pred$k/$cmdpars{'eval_against'}.cgp.gff") or die ("Could not open $optdir/pred$k/$cmdpars{'eval_against'}.cgp.gff for reading: $!");
	while(<PRED>){
	    if(/prediction on sequence range (\w+):(\d+)\-(\d+)/){
		push @intervals,[$1, $2, $3]; # store genomic intervals on which gene prediction is executed
		print JOINPRED $_;
	    }
	    if(/\tgene\t/){
		$geneID++;
	    }
	    s/g\d+/g$geneID/g; # replace gene ID with new unique gene ID
	    if(/\t(CDS|stop_codon|start_codon)\t/){
		print JOINPRED $_;
	    }	
	}
	close(PRED);
    }
    close(JOINPRED);

    # join overlapping genomic intervals
    # sort intervals by 1. chromosome, 2. start
    @intervals = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @intervals;

    my @joined=(); # array of joined genomic intervals
    my ($chr, $start, $end);
    foreach (@intervals){
	if(defined($chr) && $_->[0] eq $chr && $_->[1] - $end <= 0 ) { # overlap between the last and the current interval
	    if($end < $_->[2]){
		$end = $_->[2];
	    }
	} else {
	    if (defined($chr)){
		push @joined,[$chr, $start, $end];
		if (!$seqlist{$chr}){ # add new sequences to seqlist
		    $seqlist{$chr} = 1;
		}
	    }
	    ($chr, $start, $end) = ($_->[0], $_->[1], $_->[2]);
	}
    }
    push @joined,[$chr, $start, $end];
    if (!$seqlist{$chr}){ # add new sequences to seqlist
	$seqlist{$chr} = 1;
    }
    
    # make a new annotation file that only contains features from the training set that are completely contained in one of the intervals
    # (if necessary, this can be done faster with a single loop over the intervals, requires presorting of @gfflines)
    open(ANNO, '>', "$optdir/anno.gtf") or die ("Could not open $optdir/anno.gtf for writing: $!");
    foreach my $line (@gfflines){
	my @gffline = split(/\t/,$line);
	my ($chr, $start, $end)=($gffline[0], $gffline[3], $gffline[4]);
	foreach my $i (@joined){
	    if($chr eq $i->[0] && !($start > $i->[2]) && !($end < $i->[1]) ){
		print ANNO $line;
		last;
	    }
	}
    }
    close(ANNO);

    # filter prediction for duplicates and merge genes (uses external tool 'joingenes')
    if($joingenes){
	system("mv $optdir/pred.gtf $optdir/pred.unfiltered.gtf");
	system("$cmdpars{'jg_exec_dir'}joingenes -g $optdir/pred.unfiltered.gtf -o $optdir/pred.gtf");
    }
    
    # split annotation and prediction file by seqs and prepare
    # list files that contain the directories of the GTF files being compared (required by eval)
    system ("rm -f $optdir/annotation_list");
    system ("rm -f $optdir/prediction_list");
    foreach my $seq (keys %seqlist) {
	system ("echo '$gffDir/$seq.anno.gtf' >> $optdir/annotation_list");
	system ("echo '$gffDir/$seq.pred.gtf' >> $optdir/prediction_list");
	system ("grep \"^$seq\\b\" $optdir/pred.gtf > $gffDir/$seq.pred.gtf");
	system ("grep \"^$seq\\b\" $optdir/anno.gtf > $gffDir/$seq.anno.gtf");
    }
    # call evaluate_gtf
    my @eval_report = qx(evaluate_gtf.pl $optdir/annotation_list $optdir/prediction_list);
    # parse eval output
    foreach (@eval_report){
	if(/^Gene Sensitivity\s+(\S+)%/){
	    $cgsn = $1 * 0.01;
	}
	if(/^Gene Specificity\s+(\S+)%/){
	    $cgsp = $1 * 0.01;
	}
	if(/^Exon Sensitivity\s+(\S+)%/){
	    $cesn = $1 * 0.01;
	}
	if(/^Exon Specificity\s+(\S+)%/){
	    $cesp = $1 * 0.01;
	}
	if(/^Nucleotide Sensitivity\s+(\S+)%/){
	    $cbsn = $1 * 0.01;
	}
	if(/^Nucleotide Specificity\s+(\S+)%/){
	    $cbsp = $1 * 0.01;
	    last;
	}
    }
    if ($cbsn eq "" || $cbsp eq "" ||$cesn eq "" ||$cesp eq "" ||$cgsn eq "" ||$cgsp eq "") {
	die ("Could not retrieve accuracy values from external call of eval_gtf.pl.");
    }
    return ($cbsn,$cbsp, $cesn, $cesp, $cgsn, $cgsp, -1, -1);
}
