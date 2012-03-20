#!/usr/bin/perl
#
# wig2hints.pl
#
# convert a file with conservation scores in wiggle format to a file with 
# CDSpart hints for AUGUSTUS
#
# Mario Stanke 4.7.2007


use strict;
use Getopt::Long;

if( grep /^--?h$/, @ARGV ) {
	exec("perldoc $0");
}

my $source = "w2h";
my $strand = ".";
my $type = "CDSpart";
my $src = "X";
my $pri;
my $oldseqname = "";
my $width = 40;
my $margin = 20;
my $radius = 0;
my $minscore = 0.4;
my $minthresh = 0.4;
my $maxchunksize = 1100000; # to save memory
my $minchunksize = 1000000; # 
my $offset = 0;
my $av;
my $prune = 0;
my $step = 1;
my $UCSC = "";
my $rn = 0; # running number to make group field unique for ucsc browser
my ($seqname, $offset, $v);
my @values = ();

GetOptions('width:i'=>\$width,
	   'radius:f'=>\$radius,
	   'margin:i'=>\$margin,
	   'minscore:f'=>\$minscore,
	   'minthresh:f'=>\$minthresh,
	   'type:s'=>\$type,
	   'prune:f'=>\$prune,
	   'pri:i'=>\$pri,
	   'strand:s'=>\$strand,
           'UCSC:s'=>\$UCSC,
	   'src:s'=>\$src);

if ($UCSC ne ""){
    open (TRACK, ">$UCSC") or die ("Could not open file $UCSC for writing.");
}

while (<STDIN>) {
    if (/fixedStep chrom=(\S+)\s+start=(\d+) step=(\d+)/) {
	$seqname = $1;
	$offset = $2;
	$step = $3;
	if ($seqname ne $oldseqname && $oldseqname ne "") {
	    printScores();
	    $offset = 0;
	}
	$oldseqname = $seqname;
    } elsif (/variableStep chrom=(\S+)/){
	$seqname = $1;
	if ($seqname ne $oldseqname && $oldseqname ne "") {
	    printScores();
	    $offset = 0;
	}
	$oldseqname = $seqname;
    } elsif (/^([\.0-9]+)$/){
	$v = $1; 
	for (my $i=0; $i<$step; $i++) {
	    $values[$offset+$i] = $v;
	}
	$offset += $step;
    } elsif (/^(\d+)\s*([\.0-9]+)$/){
	$values[$1] = $2;
	if (($2 < $minthresh && $1 > $offset + $minchunksize) || ($1 > $offset + $maxchunksize)){
	    printScores();
	}
    }
}
printScores();


sub printScores {
    my ($a, $b, $aprime, $bprime) = ($offset,0,0,0);
    my $m;
    while ($a < @values) {
	while ($a<@values && $values[$a] < $minthresh){
	    $a++;
	}
	$b=$a;
	while ($b<@values && $values[$b] >= $minthresh) {
	    $b++;
	}
	$b--;
	# [a,b] is contiguous island above minthresh now
	#
	#       --------------------------------------------------------------------------------------------
	#       a                                                                                          b
	#       | margin |rest|   width   |   width    |   width    |   width    |   width   |rest| margin |
	# hints                     |           |            |            |           |
	# print "a=$a, b=$b\n";
	$aprime = $a + $margin;
	$bprime = $b - $margin;
	if ($aprime>$bprime) {
            $aprime=$bprime=int(($aprime+$bprime)/2);
        }

	#
	# big exon "mountains" often have a "foot" area
	# these two lines are trying to remove it, if $prune is >0
	#
	$av = 0;
        for (my $j=$aprime; $j<=$bprime; $j++){
            $av += $values[$j];
        }
        $av /= ($bprime-$aprime+1);

	while ($prune && $prune * $av > $values[$aprime] && $aprime < $bprime) { $aprime++;}
	while ($prune && $prune * $av > $values[$bprime] && $aprime < $bprime) { $bprime--;}

	my $len = $bprime-$aprime+1;
	my $n = int($len/$width);
	$n=1 if ($n<1);
	for (my $i=0; $i<$n; $i++) {
	    $m = int($aprime+($len - ($n-1)*$width)/2+$i*$width+0.5);
	    $av = 0;
	    for (my $j=0; $j<$width; $j++){
		$av += $values[$m+$j-int($width/2)];
	    }
	    $av = sprintf("%.3f", $av / $width);
	    if ($av > $minscore) {
		my ($from, $to) = (int($m-$radius+0.5), int($m+$radius+0.5));
		print "$oldseqname\t$source\t$type\t$from\t$to\t$av\t$strand\t.\t";
		print "src=$src;";
		print "pri=$pri;" if (defined($pri));	
		print "mult=" . int($av) . ";\n";
		print TRACK "$oldseqname\t$source\t$type\t$from\t$to\t$av\t$strand\t.\t$av $m " . ++$rn . "\n" if ($UCSC ne "");
	    }
	}
	$a = $b+1;
    }
    $offset = $a;
    @values = ();
}


close TRACK if ($UCSC ne "");


__END__

=head1 NAME

wig2hints -- convert a file with conservation scores in wiggle format to a file with 
             CDSpart/exonpart hints for AUGUSTUS

=head1 SYNOPSIS

wig2hints.pl < score.wig > hints.gff
    score.wig    - a file in wiggle format, e.g. 
                   fixedStep chrom=Nt-gDNA_Contig1538203 start=829 step=1
                   0.084
                   0.111
                   0.267
                   0.471
                   ...
                   fixedStep chrom=Nt-gDNA_Contig1538203 start=873 step=1
                   0.048
                   0.102
                   0.128
                   ...

             

=head1 OPTIONS

    --width=n         Default: 40
    --margin=n        Default: 20
    --minthresh=p     Default: 0.4 horizontal cutoff for determining intervals
    --minscore=p      Default: 0.4 consider only intervals with an average score at least this high
    --type=s          Default: print this in the third column (default:"CDSpart")
    --pri=n           print priority n in last column
    --radius=n        length of hint interval is 2*radius+1
    --strand=s        Default: "." print this in the strand column
    --prune=f         Default: 0 (off). Remove boundary areas that have coverage less than f*100% of the average over the island
    --src=s           Default: "X" print src=s in the last column
    --UCSC=s          Filename for track of UCSC genome browser custom track
=head1 DESCRIPTION
    

=cut
