#!/usr/bin/perl
#
# compute the accuracy values of a set of predictions
# against a set of annotations
#
# This script is used by the braker.pl pipeline.
# Please be extremely careful when changing this script because the braker.pl
# pipeline may fail upon custom modification of this script.
# In case of doubt, contact katharina.hoff@uni-greifswald.de
#
# Katharina Hoff & Mario Stanke, March 8th 2018

use strict;
use warnings;
use File::Path 'rmtree';
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use FindBin;
use File::Which qw(which where);    # exports which() and where()

my $usage = "$0 -- compute the prediction accuracy\n";
$usage .= "Usage: $0 seqlist annotation.gtf prediction.gtf\n";
$usage .= "Options:\n";
$usage
    .= "--EVAL_PATH   path to eval package (will try to guess from "
    . "availability of evaluate_gtf.pl, but if that does not work because you "
    . "don't have this script in your path, set with this command line option).\n";

if ( $#ARGV < 2 ) {
    die "Too few options\n\n$usage";
}
my $seqlistfilename = $ARGV[0];
my $annofilename    = $ARGV[1];
my $predfilename    = $ARGV[2];
my $help;
my $evalpath;
GetOptions(
    'EVAL_PATH=s' => \$evalpath,
    'help!'       => \$help
);

if ($help) { print $usage; }

if ($evalpath) {
    my $last_char = substr( $evalpath, -1 );
    if ( $last_char eq "\/" ) {
        chop($evalpath);
    }
}
else {
    my $epath = which 'evaluate_gtf.pl';
    $evalpath = dirname( abs_path($epath) );
}

if ( not( -e "$evalpath/evaluate_gtf.pl" ) ) {
    die("Cannot find $evalpath/evaluate_gtf.pl");
}

open( SEQLIST, "<$seqlistfilename" )
    or die("Could not open $seqlistfilename");

################################################################################
#
# create for each sequence in seqlist
# a gtf file with annotation and a gtf file with prediction
#
################################################################################
# create a temporary directory
my $dirname = "tempgtf";
system("rm -rf $dirname; mkdir $dirname");

# create the two list files of gtffiles
system("rm -f annotation_list");
system("rm -f prediction_list");

my @seqlist = <SEQLIST>;

close(SEQLIST) or die("Could not close $seqlistfilename");

foreach my $seq (@seqlist) {
    chomp $seq;
    my $cmdStr = "grep \"^$seq\\b\" $annofilename | grep -P \"\\t\\+\\t\" > "
               . "$dirname/$seq"
               . "_plus.anno.gtf";
    system($cmdStr);
    $cmdStr = "grep \"^$seq\\b\" $predfilename | grep -P \"\\t\\+\\t\" > "
            . "$dirname/$seq"
            . "_plus.pred.gtf";
    system($cmdStr);
    system( "echo '$dirname/$seq" . "_plus.anno.gtf' >> annotation_list" );
    system( "echo '$dirname/$seq" . "_plus.pred.gtf' >> prediction_list" );
    system(   "grep \"^$seq\\b\" $annofilename | grep -P \"\\t-\\t\" > "
            . "$dirname/$seq"
            . "_minus.anno.gtf" );
    system(   "grep \"^$seq\\b\" $predfilename | grep -P \"\\t-\\t\"> "
            . "$dirname/$seq"
            . "_minus.pred.gtf" );
    system( "echo '$dirname/$seq" . "_minus.anno.gtf' >> annotation_list" );
    system( "echo '$dirname/$seq" . "_minus.pred.gtf' >> prediction_list" );

}

# call evaluate_gtf
#
system( "perl -I $evalpath $evalpath/evaluate_gtf.pl annotation_list prediction_list 2> /dev/null" );

# clean up
unlink("annotation_list");
unlink("prediction_list");
rmtree( ["$dirname"] );
