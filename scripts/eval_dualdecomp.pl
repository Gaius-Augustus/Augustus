#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long; # for parameter specification on the command line

if  (eval {require Statistics::R;1;} ne 1) { # for making R plots
    # if module can't load
    print "perl module Statistics::R not installed.\n";
    print "Please install, e.g. via CPAN. On command line, type:\n\n";
    print "perl -MCPAN -e 'install Statistics::R'\n";
} else {
    Statistics::R->import();
}

my $usage = <<'ENDUSAGE';
eval_dualdecomp.pl    evaluate the effectivness of dual decomposition either on a single
                      AUGUSTUS output file or a directory of AUGUSTUS output files.
SYNOPSIS

eval_dualdecomp.pl [ --f=input_file | --d=input_directory ]

    --f=<file>                          intput AUGUSTUS file 
    --d=<dir>                           directory of input AUGUSTUS files (recognized by .out file extension)

OPTIONS

    --help                              output this help message
    --hist_iterations=<out.pdf>         output histogram of iterations to out.pdf
    --hist_errors=<out.pdf>             output histogram of error estimates to out.pdf for all cases, where
                                        no convergence is achieved.
    --err_per_iter=<out.pdf>            plots the average percentage of initial error against the iterations to out.pdf.
                                        If after a certain number of iterations the error no further drops, 
                                        early stopping is recommended, i.e. in the next run, the number of 
                                        DD iterations can be restricted to this number of iterations.
    --t=<foat>                          threshold for percentage of initial error. For all cases with an estimated
                                        error higher than this threshold, the evolution of primal an dual values
                                        are plotted against the iterations. This helps debugging cases with a
                                        high error estimate. The threshold is between [0-100] (default: 5)
    --outdir=<dir>                      put all plots in this output directory


DESCRIPTION
      
  Example:

    eval_dualdecomp.pl --d=augouts --hist_iterations=iterations.pdf --hist_errors=errors.pdf
    eval_dualdecomp.pl --d=augouts --hist_iterations=iterations.pdf --hist_errors=errors.pdf --err_per_iter=error_per_iter.pdf --outdir=out
    eval_dualdecomp.pl --f=aug.out --t=0

ENDUSAGE


my ($hist_iter, $hist_err, $err_per_iter, $DIR, $FILE, $help); # options
my $t = 5;  # threshold for percentage of initial error
my @docs = ();
my $outdir = "";

GetOptions('d=s'=>\$DIR,
	   'f=s'=>\$FILE,
	   'hist_iterations=s' =>\$hist_iter,
	   'hist_errors=s' =>\$hist_err,
	   'err_per_iter=s' =>\$err_per_iter,
	   't=f' =>\$t,
	   'outdir=s' =>\$outdir,
	   'help!'=>\$help);

if (defined($DIR)){
    $DIR =~ s/\/$//g;
    opendir(my $DH, $DIR) or die "cannot open directory $DIR: $!";
    @docs = grep(/\.out$/,readdir($DH));
    @docs = map "$DIR/$_", @docs;
}
if (defined($FILE)){
    push @docs, $FILE;
}

if(!defined($DIR) && !defined($FILE)){
    print "either an input directory (option --d) or input file (option --f) is required.\n$usage";
    exit(0);
}
if ($outdir ne ""){
    system ("rm -rf $outdir; mkdir $outdir");
    $outdir =~ s/([^\/])$/$1\//g;
}


my $id = 0;                 # gene range ID
my @avgErr = ();            # average error per iteration

my @conv_iter = ();         # stores number of iterations of all cases that did converge
my @conv_best_p = ();       # stores the iterations of the best primal values of all cases that did converge

my @not_conv_iter = ();     # stores number of iterations of all cases that did not converge
my @not_conv_best_p = ();   # stores the iterations of the best primal values of all cases that did not converge
my @not_conv_error = ();    # stores the estimated errors of all cases that did not converge

my $alreadyOpt = 0; # counts number of cases that are already optimal without DD

my @rounds = (); # number of cases that converged per round

# hash of gene ranges
#    keys: gene range nr
#    values: hash reference
#            keys: iter               array of iterations
#                  primal             array of primal values
#                  dual               array of dual values
#                  stepsize           array if step sizes
#                  inconsistencies    array of inconsistencies

my %geneRanges = ();


foreach my $file (@docs) {
    open (RES, $file) or die "could not open file $file\n";
    my $contents = 0;
    while(<RES>){
	if(/dual decomposition on gene Range (\d+)/){
	    $id++;
	    $geneRanges{$id} = {"iter"=>[], "primal"=>[], "dual"=>[], "stepsize"=>[]};
	}
	next if(!/^round\titer/ && !$contents);
	if(/^dual decomposition reduced/){
	    my $numIter = scalar @{$geneRanges{$id}{"primal"}};
	    if($numIter <= 1){ # already optimal cases, no DD required
		$alreadyOpt++;
		$contents=0;
		delete $geneRanges{$id};
		next;
	    }	    
	    my ($idx_max, $max) = findMaxValueIndex(\@{$geneRanges{$id}{"primal"}});
	    my ($idx_min, $min) = findMinValueIndex(\@{$geneRanges{$id}{"dual"}});
		
	    # calculate percentage of initial error
	    my $perc_initial_err = ($min - @{$geneRanges{$id}{"primal"}}[0]) > 0 ? ($min - $max)*100 / ($min - $geneRanges{$id}{"primal"}->[0]) : 0;
	    if($perc_initial_err < 0){ # rounding error
		$perc_initial_err = 0;
	    }
	    # plot dual and primal values against iterations for all gene ranges with an
	    # error greater than threshold t
	    if($perc_initial_err > $t){
		plot_dual_vs_primal($perc_initial_err);
	    }
	    if($geneRanges{$id}{"inconsistencies"} == 0 || $min-$max < 1e-8){ # convergence achieved
		push @conv_iter, $numIter;
		push @conv_best_p, $idx_max;
		$rounds[$geneRanges{$id}{"round"}]++;
	    }
	    else{ # not converged, stores errors 
		push @not_conv_error, $perc_initial_err;
		push @not_conv_iter, $numIter;
		push @not_conv_best_p, $idx_max; 
	    }
	    delete $geneRanges{$id};
	    $contents=0;
	    next;
	}
	if($contents){
	    my @l = split(/\s+/,$_);
	    my ($round, $iteration, $stepsize, $primal, $dual, $inconsist) = ($l[0], $l[1], $l[2], $l[3], $l[4], $l[5]);
	    if(!defined($geneRanges{$id}{"total_iter"})){
		$geneRanges{$id}{"total_iter"}=0;
	    }
	    else{
		$geneRanges{$id}{"total_iter"}++;
	    }
	    $geneRanges{$id}{"round"} = $round;
	    push @{$geneRanges{$id}{"primal"}}, $primal;
	    push @{$geneRanges{$id}{"dual"}}, $dual;
	    push @{$geneRanges{$id}{"stepsize"}}, $stepsize;
	    $geneRanges{$id}{"inconsistencies"}=$inconsist;
	    if(!defined($geneRanges{$id}{"best_primal"}) || $geneRanges{$id}{"best_primal"} < $primal){
		$geneRanges{$id}{"best_primal"} = $primal;
	    }
	    if(!defined($geneRanges{$id}{"best_dual"}) || $geneRanges{$id}{"best_dual"} > $dual){
		$geneRanges{$id}{"best_dual"} = $dual;
	    }
	    my $best_d = $geneRanges{$id}{"best_dual"};
	    my $best_p = $geneRanges{$id}{"best_primal"};
	    my $initial_p = $geneRanges{$id}{"primal"}->[0];
	    my $e = 0;
	    if(($best_d - $initial_p) > 0){
		$e = ($best_d - $best_p)*100 / ($best_d - $initial_p);
	    }
	    if($e < 0){
		$e = 0;  # rounding error
	    }
	    push @{$avgErr[$geneRanges{$id}{"total_iter"}]} , $e;
	}
	$contents=1;
    }
}
my @iter = (@conv_iter, @not_conv_iter);
my @best_p = (@conv_best_p, @not_conv_best_p);

my ($idx_error, $max_error) = (0,0);
if(@not_conv_error){
 ($idx_error, $max_error) = findMaxValueIndex(\@not_conv_error);
}
my @total_error = (@not_conv_error, ((0) x (scalar @conv_iter)));

print "\n$alreadyOpt gene Ranges were discarded, because they were already optimal prior to Dual Decomposition\n\n";
printf("+-------------------------------------+----------+-----------+-----------+----------------------+ \n");
printf("| gene Ranges      |        No.       | avg iter | avg error | max error | avg iter best primal |\n");
printf("+-------------------------------------+----------+-----------+-----------+----------------------+\n");
printf("| e-convergence    | %6u (%6.2f%%) |   %4u   |  %5.2f%%   |  %5.2f%%   |        %4u          |\n", scalar @conv_iter, (scalar @conv_iter *100 / scalar @iter),avg(\@conv_iter),0,0,avg(\@conv_best_p));
printf("| no e-convergence | %6u (%6.2f%%) |   %4u   |  %5.2f%%   |  %5.2f%%   |        %4u          |\n", scalar @not_conv_iter, (scalar @not_conv_iter *100 / scalar @iter), avg(\@not_conv_iter), avg(\@not_conv_error), $max_error, avg(\@not_conv_best_p));
printf("+-------------------------------------+----------+-----------+-----------+----------------------+\n");
printf("| total            | %6u (%6.2f%%) |   %4u   |  %5.2f%%   |  %5.2f%%   |        %4u          |\n",  scalar @iter, 100, avg(\@iter), avg(\@total_error), $max_error, avg(\@best_p));
printf("+-------------------------------------+----------+-----------+-----------+----------------------+\n\n");
print "No. of e-convergences per round of Dual Decomposition\n\n";
foreach my $i (0 .. $#rounds) {
    print "round $i - $rounds[$i]\n";
}
print "\nIf after n rounds the No. of e-convergences is only very small, the rounds\n";
print "can be restricted to n in the next run.\n\n";
    
# plot histogram of iterations
if(defined($hist_iter)){
    my $R = Statistics::R->new();
    $R->set('font',"Helvetica");
    $R->set('filename', $outdir . $hist_iter);
    $R->set('data',\@iter);
    $R->run(q`(if(require(extrafont)){library(extrafont); font <- "LM Roman 10"})`);
    $R->run(q`pdf(filename, family=font)`);
    $R->run(q`par(cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1)`);
    $R->run(q`hist(data,breaks=50,col="red",xlab="Iterations", main="")`);
    $R->run(q`dev.off()`);
}

# plot histogram of errors
if(defined($hist_err) && @not_conv_error){
    my $R = Statistics::R->new();
    $R->set('filename', $outdir . $hist_err);
    $R->set('font',"Helvetica");
    $R->set('data',\@not_conv_error);
    $R->run(q`(if(require(extrafont)){library(extrafont); font <- "LM Roman 10"})`);
    $R->run(q`pdf(filename, family=font)`);
    $R->run(q`par(cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1)`);
    $R->run(q`hist(data,breaks=50,col="red",xlab="% of initial error", main="")`);
    $R->run(q`dev.off()`);
}

# plot average percentage of initial error against iterations
if(defined($err_per_iter)){
    my $n = scalar @{$avgErr[0]};
    my @avgs = ();
    for my $i (@avgErr){ # compute average error per iteration
	my $avg = sum($i) / $n;
	push @avgs, $avg;
	if( scalar @avgs > 1  && ($avgs[$#avgs-1] - $avgs[$#avgs]  < 0.00001)){
	    last;
	}
    }
    my @x = (0..$#avgs);

    my $R = Statistics::R->new();
    $R->set('font',"Helvetica");
    $R->set('avgs',\@avgs);
    $R->set('iter',\@x);
    $R->set('filename', $outdir . $err_per_iter);
    $R->run(q`(if(require(extrafont)){library(extrafont); font <- "LM Roman 10"})`);
    $R->run(q`pdf(filename, family=font)`);
    $R->run(q`plot(iter,avgs, type="l", lwd=2, col="blue", ylim=c(min(avgs), max(avgs)), xlim=c(0,max(iter)), xlab=expression(paste("iteration ",italic(t))), ylab="average % of initial error", cex.axis=2, cex.lab=2, cex.main=2)`);
    $R->run(q`dev.off()`);	    
}

# plot dual and primal values against iterations
sub plot_dual_vs_primal{
    my $err = shift;
    my $filename = $outdir . "gr_" . $id . ".pdf";
    my @x = (0..$geneRanges{$id}{"total_iter"});
    my $R = Statistics::R->new();
    $R->set('font',"Helvetica");
    $R->set('iter',\@x);
    $R->set('primal',\@{$geneRanges{$id}{"primal"}});
    $R->set('dual',\@{$geneRanges{$id}{"dual"}});
    
    $R->set('filename', $filename);
    $R->set('error', $err);
    $R->run(q`(if(require(extrafont)){library(extrafont); font <- "LM Roman 10"})`);
    $R->run(q`pdf(filename, family=font)`);
    $R->run(q`plot(iter,primal, type="l", lwd=2, col="blue", ylim=c(min(primal), max(dual)), xlab=expression(paste("iteration ",italic(t))), ylab="value", cex.axis=2, cex.lab=2, main=bquote("Percentage of initial error" ~ italic(epsilon) == .(error) ~ "%"), cex.main=2)`);
    $R->run(q`lines(iter,dual, type="l", lwd=2, col="red")`);
    $R->run(q`legend("bottomright",c(expression(paste("current primal ",italic(p^t))),expression(paste("current dual ",italic(d^t)))),col=c("blue", "red"), lwd=c(2,2), cex=2, bty="n")`);
    $R->run(q`points(iter[which.max(primal)], max(primal), cex = 1.5, pch = 20)`);
    $R->run(q`text(iter[which.max(primal)], max(primal), labels = expression(italic(p)[best]), cex= 1.5, pos=3)`);
    $R->run(q`dev.off()`);	
}

sub findMaxValueIndex{
    my $idx;
    my $max;
    my @array = @{$_[0]};
    for my $i (0 .. $#array){
	if(!defined($max) || $max < $array[$i]){
	    $idx = $i;
	    $max = $array[$i];
	}
    }
    return ($idx, $max);
}

sub findMinValueIndex{
    my $idx;
    my $min;
    my @array = @{$_[0]};
    for my $i (0 .. $#array){
	if(!defined($min) || $min > $array[$i]){
	    $idx = $i;
	    $min = $array[$i];
	}
    }
    return ($idx, $min);
}

sub avg{
    my @array = @{$_[0]};
    if(@array){
	my $sum = sum(\@array);
	$sum /= scalar @array;
	return $sum;
    }
    return 0;
}

sub sum{
    my $sum = 0;
    my @array = @{$_[0]};
    for my $i (@array){
	$sum+=$i;
    }
    return $sum;
}
