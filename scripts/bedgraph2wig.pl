#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


#Method to convert the given hash into WIG format and printing it to a handle. Make sure the handle is either STDOUT or a filehandle
sub convertAndWriteWig {
  my $data = shift;
  my $handle = shift;
  my $span = shift;
  if ($span) {
    for (my $i=0;$i < scalar(@{$data->{'start'}});$i++) {
      if(@{$data->{'coverage'}}[$i] != 0) {
        print($handle "variableStep chrom=$data->{'chromosome'} span=",@{$data->{end}}[$i]-@{$data->{start}}[$i],"\n@{$data->{start}}[$i] @{$data->{'coverage'}}[$i]\n");
      }
    }
  }
  else {
    print($handle qq(variableStep chrom=$data->{'chromosome'}\n));
    for (my $i=0;$i < scalar(@{$data->{'start'}});$i++) {
      for (my $j= @{$data->{'start'}}[$i]; $j < @{$data->{'end'}}[$i]+1;$j++) {
	if(@{$data->{'coverage'}}[$i] != 0) {
	  print($handle qq($j @{$data->{'coverage'}}[$i]\n));
	}
      }
    }
  }
}



sub main {
  my $usage = "USAGE:\n\t";
     $usage.= "bedgraph2wig.pl --bedgraphfile=<filename> [Output is STDOUT]\n\t";
     $usage.= "bedgraph2wig.pl --bedgraphfile=<filename> --outputfile=<filename>\n\n\n";
     $usage.= "Parameters:\n";
     $usage.= "\t--bedgraphfile: The path to the bedgraphffile\n";
     $usage.= "\t--outputfile: The file where the output should be saved. If not supplied, the output is printed to STDOUT\n";
     $usage.= "\t--span: Use span notation (see http://genome.ucsc.edu/goldenPath/help/wiggle.html\n"; 

  my $bedgraphfile;
  my $outputfile;
  my $help;
  my $span;  
  
  #test if the input is ok
  if (scalar(@ARGV)==0) {print "$usage\n"; exit(0);}
  my $found;
  foreach (@ARGV) {
    if ($_ =~ m/--bedgraphfile/) {
      $found=1;
    }
  }
  if (!$found) {
    print "$usage\n"; exit(0);
  }
  
  #fetch params
  GetOptions('bedgraphfile=s' => \$bedgraphfile,
             'outputfile=s' => \$outputfile,
	     'help!' => \$help,
	     'span!' => \$span);
   
  if ($help) { print "$usage\n"; exit(0); } 
 
  if (defined($bedgraphfile)) {
    open(BEDGRAPH, "<$bedgraphfile") || die "Couldn't open $bedgraphfile\n";
  }
  else {
    die "No begraphfile supplied\n";
  }
  
  #Check if we should print to STDOUT or a file. We don't print to a file if it exists.
  my $fh=\*STDOUT;
  if (defined($outputfile)) { 
    if (!(-e $outputfile)) {
      open $fh, '>>', $outputfile or $fh=\*STDOUT; 
    } 
    else {
      print(STDERR  "Outputfile allready exists, printing to STDOUT\n");
    }

  }
  #initialize the data hash with arrays and strings we need
  my %data=('chromosome'=>"",start=>[],end=>[],coverage=>[]); 

  my $lastchromosome;
  #we dont want the first line to be converted, its the track line 
  my $infostring = <BEDGRAPH>;
  chomp($infostring);
  $infostring =~ s/.*name=(".*")/$1/i;
  print($fh "track type=wiggle_0 name=$infostring \n");
  while(<BEDGRAPH>) {
      print "$_";
    #print($fh "$_\n");
    my ($chromosome,$start,$end,$coverage) = split(/\s+/,$_);
    # print($fh "$chromosome  $start  $end  $coverage \n");
    if (defined($lastchromosome)) {
      { 
        if($chromosome ne $lastchromosome) {
          &convertAndWriteWig(\%data,$fh,$span);
	}
      }
    }

    $data{'chromosome'} = $chromosome;
    push(@{$data{'start'}},$start);
    push(@{$data{'end'}},$end);
    push(@{$data{'coverage'}},$coverage);
    $lastchromosome = $chromosome; 
  }    
    

}

&main;
