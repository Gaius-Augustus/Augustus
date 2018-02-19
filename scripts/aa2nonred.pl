#!/usr/bin/perl
#
# take a multiple fasta amino acid sequence file
# and output a file that is non-redundant:
#
# the percent identity value between each pair of different sequence is below a threshold
# and each sequence name occurrs only once
#
# Mario Stanke & Katharina Hoff, last modification on Feb 19th 2018

use strict;
use Getopt::Long;
use File::Which qw(which where);
use File::Spec::Functions qw(rel2abs);

#
# The maximum percent identiy allowed between any two sequences.
# Here defined as (number of identical aa)/(length of shorter sequence)
#
my $max_percent_id = 0.8;
my $wublast        = 0;
my $BLAST_PATH;
my $blast_path;

my $usage = "aa2nonred.pl -- make a protein file non-redundant\n";
$usage .= "Usage: aa2nonred.pl input.fa output.fa\n";
$usage .= "In output.fa the percent identity value between each pair of \n";
$usage .= "When removing redundant sequences, priority is given to the sequence occuring last.\n";
$usage .= "Options:\n";
$usage .= "--maxid=f       maximum percent identity between to sequences\n";
$usage .= "                (#identical aa) / (length of shorter sequence) default: 0.8\n";
$usage .= "--wublast       flag to turn on WUBLAST instead of NCBI-BLAST (i.e. setdb, blastp)\n";
$usage .= "--BLAST_PATH=s  path to blast (only implemented for NCBI BLAST)\n";

GetOptions(
    'maxid:f'  => \$max_percent_id,
    'wublast!' => \$wublast,
    'BLAST_PATH=s' => \$blast_path
);

if ( $#ARGV != 1 ) {
    die "Unknown option\n\n$usage";
}

if ( !$wublast ) {
    set_BLAST_PATH();
}

my $inputfilename  = $ARGV[0];
my $outputfilename = $ARGV[1];

###########################################################################################
#
# make a file with unique sequence names
#
###########################################################################################

open( INPUT, "<$inputfilename" ) or die("Could not open $inputfilename");
my $tempdbname = "$inputfilename.nonreddb";
open( TEMP, ">$tempdbname " ) or die("Could not open $tempdbname");

my %seqnames = ();
my $seqname;
$/ = "\n>";
while (<INPUT>) {
    s/\n>//g;
    s/^>//g;
    /^(.*)\n/g;
    $seqname = $1;
    if ( !exists( $seqnames{$seqname} ) ) {
        $seqnames{$seqname} = $';    #'
        print TEMP ">$seqname\n" . $seqnames{$seqname} . "\n";
    }
}

###########################################################################################
#
# blast that file against itself
#
###########################################################################################

my $tempoutfile = "$inputfilename.blastout";

if ($wublast) {

    #my $blastdir = "~/wublast";
    system("setdb $tempdbname");
    system("blastp $tempdbname $tempdbname -E=1e-5 -B=100 > $tempoutfile");
}
else {
    ## NCBI blast
    system("$BLAST_PATH/makeblastdb -in $inputfilename -dbtype prot -parse_seqids -out $tempdbname");
    system("$BLAST_PATH/blastp -query $inputfilename -db $tempdbname > $tempoutfile");
}


###########################################################################################
#
# parse the blast output
#
###########################################################################################

open( BLASTOUT, "<$tempoutfile" ) or die("Could not open $tempoutfile");
$/ = "\nQuery= ";
my ( $query, $target, $qlen, $tlen, $numid, $minlen );
while (<BLASTOUT>) {
    next unless / producing /;
    /\s*(.*)\n.*\((\d+) letters.*\)/;
    $query = $1;
    $qlen  = $2;
    print STDOUT "query=$query, qlen=$qlen\n";
    while (
        (   $wublast
            && />(.*)\n\s+Length = (\d+)\n.*\n.*\n Identities = (\d+)/g
        )
        || ( !$wublast
            && />(.*)\n\s+Length = (\d+)\n.*\n.*\n Identities = (\d+)/g )
        )
    {
        $target = $1;
        $tlen   = $2;
        $numid  = $3;
        print STDOUT "target=$target, tlen=$tlen, numid=$numid\n";
        $minlen = ( $qlen < $tlen ) ? $qlen : $tlen;
        if ( $minlen == 0 ) { $minlen = 0.0000001; }
        if ( $numid / $minlen > $max_percent_id ) {    # conflict: too similar
            if ( $query ne $target ) {
                print STDOUT "$query and $target are similar\n";
                if (   exists( $seqnames{$query} )
                    && exists( $seqnames{$target} ) )
                {
                    delete( $seqnames{$query} );
                    print STDOUT "delete $query\n";
                }
            }
        }
    }
}

###########################################################################################
#
# output the nonredundant file
#
###########################################################################################
open( OUTPUT, ">$outputfilename" ) or die("Could not open $outputfilename");

foreach $seqname ( keys %seqnames ) {
    print OUTPUT ">$seqname\n" . $seqnames{$seqname} . "\n";
}


###########################################################################################
#
# Clean up
#
###########################################################################################
unlink ( rel2abs($tempdbname) );
unlink ( rel2abs($tempoutfile) );


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
        print STDOUT "After the which! $epath\n"
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
                   .  " aa2nonred.pl:\n";
        $blast_err .= "   a) provide command-line argument "
                   .  "--BLAST_PATH=/your/path\n";
        $blast_err .= "   b) use an existing environment variable "
                   .  "\$BLAST_PATH\n";
        $blast_err .= "      for setting the environment variable, run\n";
        $blast_err .= "           export BLAST_PATH=/your/path\n";
        $blast_err .= "      in your shell. You may append this to your "
                   .  ".bashrc or .profile file in\n";
        $blast_err .= "      order to make the variable available to all your "
                   .  "bash sessions.\n";
        $blast_err .= "   c) aa2nonred.pl can try guessing the location of "
                   .  "\$BLAST_PATH from the\n";
        $blast_err .= "      location of a blastp executable that is "
                   .  "available in your \$PATH variable.\n";
        $blast_err .= "      If you try to rely on this option, you can check "
                   .  "by typing\n";
        $blast_err .= "           which blastp\n";
        $blast_err .= "      in your shell, whether there is a blastp "
                   .  "executable in your \$PATH\n";
        $prtStr = "\# " . (localtime) . " ERROR: \$BLAST_PATH not set!\n";
        print STDERR $prtStr;
        print STDERR $blast_err;
        exit(1);
    }
    if ( not ( -x "$BLAST_PATH/blastp" ) ) {
        print STDERR "\# " . (localtime) . " ERROR: $BLAST_PATH/blastp is not an executable file!\n";
        exit(1);
    }elsif( not ( -x "$BLAST_PATH/makeblastdb" ) ){
        print STDERR "\# " . (localtime) . " ERROR: $BLAST_PATH/makeblastdb is not an executable file!\n";
        exit(1);
    }
}