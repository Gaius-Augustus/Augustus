#!/usr/bin/env perl
#
# take a multiple fasta amino acid sequence file
# and output a file that is non-redundant:
#
# the percent identity value between each pair of different sequence is below a threshold
# and each sequence name occurrs only once
#
# This script is used by the braker.pl pipeline.
# Please be extremely careful when changing this script because the braker.pl
# pipeline may fail upon custom modification of this script.
# In case of doubt, contact katharina.hoff@uni-greifswald.de
#
# Mario Stanke & Katharina Hoff, last modification on May 31 2019

use strict;
use Getopt::Long;
use File::Which qw(which where);
use File::Spec::Functions qw(rel2abs);
use Cwd 'abs_path';
use File::Path qw(rmtree make_path);
use File::Basename qw(dirname basename);
use Parallel::ForkManager; # native blastp parallelization keeps most nodes idle most of the time, therefore data parallelize in perl
use Scalar::Util qw(openhandle);
use POSIX;

#
# The maximum percent identiy allowed between any two sequences.
# Here defined as (number of identical aa)/(length of shorter sequence)
#
my $max_percent_id = 0.8;
my $BLAST_PATH;
my $DIAMOND_PATH;
my $blast_path;
my $diamond_path;
my $CPU = 1;
my $diamond;
my $v = 0;
my $help;
my $cmdString;

my $usage = "aa2nonred.pl -- make a protein file non-redundant\n";
$usage .= "Usage: aa2nonred.pl input.fa output.fa\n";
$usage .= "In output.fa the percent identity value between each pair of \n";
$usage .= "When removing redundant sequences, priority is given to the sequence occuring last.\n";
$usage .= "Options:\n";
$usage .= "--maxid=f         maximum percent identity between to sequences\n";
$usage .= "                  (#identical aa) / (length of shorter sequence) default: 0.8\n";
$usage .= "--BLAST_PATH=s    path to blast (only implemented for NCBI BLAST)\n";
$usage .= "--DIAMOND_PATH=s  path to diamond\n"; 
$usage .= "--cores=n         number of cores to be used by NCBI BLAST or DIAMOND\n";
$usage .= "--diamond         use DIAMOND istead of NCBI BLAST\n";
$usage .= "--verbosity=n     verbosity level for information printed to stdout\n";
$usage .= "--help            print this help message\n";

GetOptions(
    'maxid:f'  => \$max_percent_id,
    'BLAST_PATH=s' => \$blast_path,
    'DIAMOND_PATH=s' => \$diamond_path,
    'cores=i'  => \$CPU,
    'diamond!' => \$diamond,
    'verbosity=i' => \$v,
    'help!'    => \$help
);

if ($help) {
    print $usage;
    exit(1);
}

if (scalar(@ARGV) < 2){
    print "Too few input arguments!\n";
    print $usage;
    exit(1);
}

if($diamond){
    set_DIAMOND_PATH();
}else{
    set_BLAST_PATH();
}


my $inputfilename  = $ARGV[0];
my $outputfilename = $ARGV[1];

###########################################################################################
#
# split input file
#
###########################################################################################

my $splitDir; # for parallelization
my @splitFiles;
my $SPLITF;
if ( $CPU > 1 && not($diamond)) {

    my $nFastaEntries = 0;
    # counter number of fasta entries
    open( INPUT, "<$inputfilename" ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $inputfilename!\n");
    while (<INPUT>) {
        if ($_ =~ m/^>/) {
            $nFastaEntries++;
        }
    }
    close ( INPUT ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $inputfilename!\n");

    $splitDir = dirname( abs_path($inputfilename) ) . "/split_blast";
    make_path ($splitDir) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to create directory $splitDir!\n");
    my $maxEntries = ceil($nFastaEntries/$CPU);
    my $fileNr = 1;
    my $nEntries = 0;
    my $filename;
    open( INPUT, "<$inputfilename" ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $inputfilename!\n");
    while (<INPUT>) {
        if ($_ =~ m/^>(\S+)/) {
            if ( defined ( openhandle($SPLITF)) && $nEntries == $maxEntries) {
                close ($SPLITF) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $filename!\n");
                $fileNr++;
                $nEntries = 0;
            }
            if($nEntries==0){
                $filename = "$splitDir/split_$fileNr.fa";
                push @splitFiles, $filename;
                open ( $SPLITF, ">", $filename) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open file $filename!\n");
            }
            $nEntries++;
            print $SPLITF $_;

        }else{
            if ( defined ( openhandle($SPLITF) ) ) {
                print $SPLITF $_;
            }
        }

    }
    if ( defined ( openhandle($SPLITF)) ) {
        close ($SPLITF) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $filename!\n");
    }
    close ( INPUT ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $inputfilename!\n");
    foreach(@splitFiles){
        print "$_\n";
    }
}

###########################################################################################
#
# make a file with unique sequence names
#
###########################################################################################

open( INPUT, "<$inputfilename" ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $inputfilename!\n");
my $tempdbname = "$inputfilename.nonreddb";
open( TEMP, ">$tempdbname " ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $tempdbname!\n");

my %seqnames = ();
my $seqname;
$/ = "\n\>";
while (<INPUT>) {
    $_ =~ s/\>//g;
    $_ =~ m/^(\S+)/;
    $seqname = $1;
    if ( !exists( $seqnames{$seqname} )  && length($seqname) > 1) {
        $seqnames{$seqname} = $';    #'
        print TEMP ">$seqname" . $seqnames{$seqname} . "";
    }
}
close (INPUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $inputfilename!\n");
close (TEMP) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $tempdbname!\n");

###########################################################################################
#
# blast that file against itself
#
###########################################################################################

my $tempoutfile = "$inputfilename.blastout";

if($diamond){
    if($CPU > 1){
        $cmdString = "$DIAMOND_PATH/diamond makedb --in $tempdbname -d $tempdbname --threads $CPU";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
    }else{
        $cmdString = "$DIAMOND_PATH/diamond makedb --in $tempdbname -d $tempdbname --threads 1";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
    }
}else{
    ## NCBI blast
    $cmdString = "$BLAST_PATH/makeblastdb -in $tempdbname -dbtype prot -parse_seqids -out $tempdbname";
    system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
}

if ( $CPU == 1 ) {
    if($diamond){
        $cmdString = "$DIAMOND_PATH/diamond blastp --db $tempdbname --outfmt 0 --query $tempdbname --out $tempoutfile";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");

    }else{
        $cmdString = "$BLAST_PATH/blastp -query $tempdbname -db $tempdbname > $tempoutfile";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
    }
}else{
    if($diamond){
        $cmdString = "$DIAMOND_PATH/diamond blastp --db $tempdbname --outfmt 0 --query $tempdbname --threads $CPU --out $tempoutfile";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
    }else{
        my $pm = new Parallel::ForkManager($CPU);
        foreach ( @splitFiles ) {
            my $pid = $pm->start and next;
            $cmdString = "$BLAST_PATH/blastp -query $_ -db $tempdbname > $_.blastout";
            system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
            $pm->finish;
        }
        $pm->wait_all_children;
        foreach ( @splitFiles ) {
            $cmdString = "cat $_.blastout >> $tempoutfile";
            system($cmdString) == 0 or die("ERROR in file " . __FILE__ . " at line " . __LINE__ ."\nFailed to execute: $cmdString!\n");
        }
    }
}


###########################################################################################
#
# parse the blast/diamond output
#
###########################################################################################

open( BLASTOUT, "<$tempoutfile" ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $tempoutfile!\n");
$/ = "\nQuery= ";
my ( $query, $target, $qlen, $tlen, $numid, $minlen );
while (<BLASTOUT>) {
    if(not($diamond)){
        next unless / producing /;
    }
    $_ =~ m/(\S+)\n+Length=(\d+)/;
    $query = $1;
    $qlen  = $2;
    print STDOUT "query=$query, qlen=$qlen\n" if ($v>0);
    while ( $_ =~ m/>(\S+).*\nLength=(\d+)\n.*\n.*\n Identities = (\d+)/g) {
        $target = $1;
        $tlen   = $2;
        $numid  = $3;
        print STDOUT "target=$target, tlen=$tlen, numid=$numid\n" if ($v>0);
        $minlen = ( $qlen < $tlen ) ? $qlen : $tlen;
        if ( $minlen == 0 ) { $minlen = 0.0000001; }
        if ( $numid / $minlen > $max_percent_id ) {    # conflict: too similar
            if ( $query ne $target ) {
                print STDOUT "$query and $target are similar\n" if ($v>0);
                if (   exists( $seqnames{$query} )
                    && exists( $seqnames{$target} ) )
                {
                    delete( $seqnames{$query} );
                    print STDOUT "delete $query\n" if ($v>0);
                }
            }
        }
    }

}
close (BLASTOUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $tempoutfile!\n");

###########################################################################################
#
# output the nonredundant file
#
###########################################################################################
open( OUT, ">", $outputfilename ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not open $outputfilename!\n");

foreach $seqname ( keys %seqnames ) {
    print OUT ">$seqname" . $seqnames{$seqname} ;
}

close( OUT ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close $outputfilename!\n");
###########################################################################################
#
# Clean up
#
###########################################################################################
unlink ( rel2abs($tempdbname) );
if(not($diamond)){
    unlink ( rel2abs($tempdbname).".phr" );
    unlink ( rel2abs($tempdbname).".pin" );
    unlink ( rel2abs($tempdbname).".pog" );
    unlink ( rel2abs($tempdbname).".psd" );
    unlink ( rel2abs($tempdbname).".psi" );
    unlink ( rel2abs($tempdbname).".psq" );
}else{
    unlink (rel2abs($tempdbname).".dmnd" );
}
unlink ( rel2abs($tempoutfile) );
if ($CPU > 1 && not($diamond)) {
    rmtree( ["$splitDir"] );
}

###########################################################################################
#
# finding blastp and makeblastdb executables
#
###########################################################################################
sub set_BLAST_PATH {
    my $prtStr;
    # try to get path from ENV
    if ( defined( $ENV{'BLAST_PATH'} ) ) {
        if ( -e $ENV{'BLAST_PATH'} ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Found environment variable \$BLAST_PATH.\n";
            print STDOUT $prtStr;
            $BLAST_PATH = $ENV{'BLAST_PATH'};
        }
    }
    else {
        $prtStr
            = "\# "
            . (localtime)
            . ": Did not find environment variable \$BLAST_PATH\n";
        print STDOUT $prtStr;
    }

    # try to get path from command line
    if ( defined($blast_path) ) {
        my $last_char = substr( $blast_path, -1 );
        if ( $last_char eq "\/" ) {
            chop($blast_path);
        }
        if ( -d $blast_path ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Setting \$BLAST_PATH to command line argument ";
            $prtStr .= "--BLAST_PATH value $blast_path.\n";
            print STDOUT $prtStr;
            $BLAST_PATH = $blast_path;
        }
        else {
            $prtStr
                = "\# "
                . (localtime)
                . ": WARNING: Command line argument --BLAST_PATH was ";
            $prtStr
                .= "supplied but value $blast_path is not a directory. Will not set ";
            $prtStr .= "\$BLAST_PATH to $blast_path!\n";
            print STDOUT $prtStr;
        }
    }

    # try to guess
    if ( not( defined($BLAST_PATH) )
        || length($BLAST_PATH) == 0 )
    {
        $prtStr
            = "\# "
            . (localtime)
            . ": Trying to guess \$BLAST_PATH from location of blastp";
        $prtStr .= " executable that is available in your \$PATH.\n";
        print STDOUT $prtStr;
        my $epath = which 'blastp';
        if ( -d dirname($epath) ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Setting \$BLAST_PATH to "
                . dirname($epath) . "\n";
            print STDOUT $prtStr;
            $BLAST_PATH = dirname($epath);
        }
        else {
            $prtStr
                = "\# "
                . (localtime)
                . ": WARNING: Guessing the location of \$BLAST_PATH ";
            $prtStr
                .= "failed. " . dirname($epath) . " is not a directory!\n";
            print STDOUT $prtStr;
        }
    }

    if ( not( defined($BLAST_PATH) ) ) {
        my $blast_err;
        $blast_err .= "There are 3 alternative ways to set this variable for "
                   .  " aa2nonred.pl:\n"
                   .  "   a) provide command-line argument "
                   .  "--BLAST_PATH=/your/path\n"
                   .  "   b) use an existing environment variable "
                   .  "\$BLAST_PATH\n"
                   .  "      for setting the environment variable, run\n"
                   .  "           export BLAST_PATH=/your/path\n"
                   .  "      in your shell. You may append this to your "
                   .  ".bashrc or .profile file in\n"
                   .  "      order to make the variable available to all your "
                   .  "bash sessions.\n"
                   .  "   c) aa2nonred.pl can try guessing the location of "
                   .  "\$BLAST_PATH from the\n"
                   .  "      location of a blastp executable that is "
                   .  "available in your \$PATH variable.\n"
                   .  "      If you try to rely on this option, you can check "
                   .  "by typing\n"
                   .  "           which blastp\n"
                   .  "      in your shell, whether there is a blastp "
                   .  "executable in your \$PATH\n";
        $prtStr = "\# " . (localtime) . " ERROR in file " . __FILE__ ." at line ". __LINE__ ."\n\$BLAST_PATH not set!\n";
        print STDERR $prtStr;
        print STDERR $blast_err;
        exit(1);
    }
    if ( not ( -x "$BLAST_PATH/blastp" ) ) {
        print STDERR "\# " . (localtime) . " ERROR in file " . __FILE__ ." at line ". __LINE__ ."\n$BLAST_PATH/blastp is not an executable file!\n";
        exit(1);
    }elsif( not ( -x "$BLAST_PATH/makeblastdb" ) ){
        print STDERR "\# " . (localtime) . " ERROR in file " . __FILE__ ." at line ". __LINE__ ."\n$BLAST_PATH/makeblastdb is not an executable file!\n";
        exit(1);
    }
}

###########################################################################################
#
# finding diamond executable
#
###########################################################################################
sub set_DIAMOND_PATH {
    my $prtStr;
    # try to get path from ENV
    if ( defined( $ENV{'DIAMOND_PATH'} ) ) {
        if ( -e $ENV{'DIAMOND_PATH'} ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Found environment variable \$DIAMOND_PATH.\n";
            print STDOUT $prtStr;
            $DIAMOND_PATH = $ENV{'DIAMOND_PATH'};
        }
    }
    else {
        $prtStr
            = "\# "
            . (localtime)
            . ": Did not find environment variable \$DIAMOND_PATH\n";
        print STDOUT $prtStr;
    }

    # try to get path from command line
    if ( defined($diamond_path) ) {
        my $last_char = substr( $diamond_path, -1 );
        if ( $last_char eq "\/" ) {
            chop($diamond_path);
        }
        if ( -d $diamond_path ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Setting \$DIAMOND_PATH to command line argument ";
            $prtStr .= "--DIAMOND_PATH value $diamond_path.\n";
            print STDOUT $prtStr;
            $DIAMOND_PATH = $diamond_path;
        }
        else {
            $prtStr
                = "\# "
                . (localtime)
                . ": WARNING: Command line argument --DIAMOND_PATH was ";
            $prtStr
                .= "supplied but value $diamond_path is not a directory. Will not set ";
            $prtStr .= "\$DIAMOND_PATH to $diamond_path!\n";
            print STDOUT $prtStr;
        }
    }

    # try to guess
    if ( not( defined($DIAMOND_PATH) )
        || length($DIAMOND_PATH) == 0 )
    {
        $prtStr
            = "\# "
            . (localtime)
            . ": Trying to guess \$DIAMOND_PATH from location of diamond";
        $prtStr .= " executable that is available in your \$PATH.\n";
        print STDOUT $prtStr;
        my $epath = which 'diamond';
        if ( -d dirname($epath) ) {
            $prtStr
                = "\# "
                . (localtime)
                . ": Setting \$DIAMOND_PATH to "
                . dirname($epath) . "\n";
            print STDOUT $prtStr;
            $DIAMOND_PATH = dirname($epath);
        }
        else {
            $prtStr
                = "\# "
                . (localtime)
                . ": WARNING: Guessing the location of \$DIAMOND_PATH ";
            $prtStr
                .= "failed. " . dirname($epath) . " is not a directory!\n";
            print STDOUT $prtStr;
        }
    }

    if ( not( defined($DIAMOND_PATH) ) ) {
        my $diamond_err;
        $diamond_err .= "There are 3 alternative ways to set this variable for "
                     .  " aa2nonred.pl:\n"
                     .  "   a) provide command-line argument "
                     .  "--DIAMOND_PATH=/your/path\n"
                     .  "   b) use an existing environment variable "
                     .  "\$DIAMOND_PATH\n"
                     .  "      for setting the environment variable, run\n"
                     .  "           export DIAMOND_PATH=/your/path\n"
                     .  "      in your shell. You may append this to your "
                     .  ".bashrc or .profile file in\n"
                     .  "      order to make the variable available to all your "
                     .  "bash sessions.\n"
                     .  "   c) aa2nonred.pl can try guessing the location of "
                     .  "\$DIAMOND_PATH from the\n"
                     .  "      location of a diamond executable that is "
                     .  "available in your \$PATH variable.\n"
                     .  "      If you try to rely on this option, you can check "
                     .  "by typing\n"
                     .  "           which diamond\n"
                     .  "      in your shell, whether there is a diamond "
                     .  "executable in your \$PATH\n";
        $prtStr = "\# " . (localtime) . " ERROR in file " . __FILE__ ." at line ". __LINE__ ."\n\$DIAMOND_PATH not set!\n";
        print STDERR $prtStr;
        print STDERR $diamond_err;
        exit(1);
    }
    if ( not ( -x "$DIAMOND_PATH/diamond" ) ) {
        print STDERR "\# " . (localtime) . " ERROR in file " . __FILE__ ." at line ". __LINE__ ."\n$DIAMOND_PATH/diamond is not an executable file!\n";
        exit(1);
    }
}
