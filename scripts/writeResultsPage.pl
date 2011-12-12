#!/usr/bin/perl

# This scripts writes an HTML results page and does all the necessary file modifications for making results of the AUGUSTUS training web server application available.

# ID is the random project ID of a grails project
# species-name
# db-file is the grails database flatfile 
# grails-out is the directory where grails places the autoAug output files
# www-out is the directory where results shall be made available to end users (i.e. an apache directory)
	
my $usage = "writeResultsPage ID species-name db-file grails-out www-out AugustusConfigPath AugustusScriptsPath final-flag\n";

if(@ARGV != 8){
        print $usage;
        exit;
}

## automatically configured variables
my $projectID = $ARGV[0];
my $species = $ARGV[1];
my $dbFile = $ARGV[2];
my $grailsOut = $ARGV[3];
my $wwwOut = $ARGV[4];
my $AUGUSTUS_CONFIG_PATH = $ARGV[5];
my $svnScripts = $ARGV[6];
my $final = $ARGV[7];

## retrieve all submission data information that will later be used to create the results page

my $dbLine = `grep $projectID $dbFile`;
my @tdbLine = split(/\[|\]/, $dbLine);
my $submissionDate = $tdbLine[1];

my $projectWebOutDir = $wwwOut."/$projectID";

if($final == 0){
	## create the webserver output directory
	`mkdir $projectWebOutDir`;
	if (not(-d "$projectWebOutDir")){
		print STDERR "Creating the output directory in for apache failed! Check writing permissions!\n";
		exit;
	}
	## create emtpy results page
	my $head = $svnScripts."/webserver-results.head";
	my $body = $svnScripts."/webserver-results.body";
	my $tail = $svnScripts."/webserver-results.tail";
	my $segmentFile1 = $grailsOut."/$projectID/segment1";
	my $segmentFile2 = $grailsOut."/$projectID/segment2";
	open(SEG1, ">", $segmentFile1) or die "Could not open file $segmentFile1!\n";
	print SEG1 "Results for job $projectID\n";
	close(SEG1) or die "Could not close file $segmentFile1!\n";
	open(SEG2, ">", $segmentFile2) or die "Could not open file $segmentFile2!\n";
	print SEG2 "<a href=\"index.html\" class=\"contentpagetitle\">Results for job $projectID</a>\n</td>\n</tr>\n</table>\n";
	print SEG2 "<p>There are no results available, yet!</p>\n";
	close(SEG2) or die "Could not close file $segmentFIle2!\n";
	$cmdStr = "cat $head $segmentFile1 $body $segmentFile2 $tail > $projectWebOutDir/index.html";
	`$cmdStr`;
	
}

if($final == 1){
## create a paramter folder in web-out folder
if($projectID =~ m/^t/){
	print STDOUT "We are in training mode\n";
	my $projectWebOutParams = $wwwOut."/$projectID/$species";
	`mkdir $projectWebOutParams`;
	if (not(-d "$projectWebOutParams")){
		print STDERR "Creating the parameter directory in for apache failed! Check writing permissions!\n";
		exit;
	}

	## copy 'adapted' parameter files to web-out parameter folder

	my $cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_parameters.cfg";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_parameters.cfg";
	system "$cmdStr\n";
	if (not(-e $cfgFile)){
		print STDERR "$projectWebOutParams/$species"."_parameters.cfg was not written! Check writing permissions!\n";
	}
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_metapars.cfg";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_metapars.cfg";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_metapars.utr.cfg";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_metapars.utr.cfg";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_exon_probs.pbl";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_exon_probs.pbl";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_exon_probs.pbl.withoutCRF";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_exon_probs.pbl.withoutCRF";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_igenic_probs.pbl";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_igenic_probs.pbl";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_igenic_probs.pbl.withoutCRF";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_igenic_probs.pbl.withoutCRF";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_intron_probs.pbl";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_intron_probs.pbl";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_intron_probs.pbl.withoutCRF";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_intron_probs.pbl.withoutCRF";
	system "$cmdStr\n";
	$cfgFile = $AUGUSTUS_CONFIG_PATH."/species/$projectID/$projectID"."_weightmatrix.txt";
	if (not(-e $cfgFile)){
		print STDERR "$cfgFile does not exist!\n";
	}
	$cmdStr = "cat $cfgFile | perl -pe 's/$projectID/$species/;' > $projectWebOutParams/$species"."_weightmatrix.txt";
	system "$cmdStr\n";

	## pack parameters
	$cmdStr = "cd $projectWebOutDir; tar -czvf parameters.tar.gz $species &> /dev/null;";
	`$cmdStr`;
	if (not(-e "$projectWebOutDir/parameters.tar.gz")){
		print STDERR "$projectWebOutDir/parameters.tar.gz was not packed!\n";
	}

	## remove original parameter directoy from apache directory
	$cmdStr = "rm -r $projectWebOutDir/$species";
	system "$cmdStr\n";

	## copy and pack training gene file
	my $trainingFile = $grailsOut."/$projectID/autoAug/trainingSet/training/training.gb";
	if (not(-e $trainingFile)){
		print STDERR "$trainingFile does not exist!\n";
	}
	$cmdStr = "cp $trainingFile $projectWebOutDir/training.gb; cd $projectWebOutDir; gzip training.gb &> /dev/null;";
	`$cmdStr`;
	if (not(-e $projectWebOutDir."/training.gb.gz")){
	print STDERR "$projectWebOutDir/training.gb.gz was not packed!\n";
	}


	## copy and pack ab-initio output file
	my $ab_initio_webDir = $projectWebOutDir."/ab_initio";
	$cmdStr = "mkdir $ab_initio_webDir";
	system $cmdStr;
	my $ab_initio_grailsDir = $grailsOut."/$projectID/autoAug/autoAugPred_abinitio";
	if (not(-d "$ab_initio_grailsDir")){
		print STDERR "AutoAug did not produce ab initio predictions!\n";
		exit;
	}
	$cmdStr = "cp  $ab_initio_grailsDir/predictions/* $ab_initio_webDir; cp $ab_initio_grailsDir/gbrowse/* $ab_initio_webDir;";
	`$cmdStr`;
	$cmdStr = "cd $projectWebOutDir; tar -czvf ab_initio.tar.gz ab_initio;";
	`$cmdStr`;
	$cmdStr = "rm -r $ab_initio_webDir;";
	`$cmdStr`;

	## copy and pack hints predictions - if they exist	
	my $hintPredsExistFlag = 0;
	if(-d $grailsOut."/$projectID/autoAug/autoAugPred_hints"){
		$hintsPredsExistFlag = 1;
	}
	if($hintsPredsExistFlag==1){
		my $hintsPred_webDir = $projectWebOutDir."/hints_pred";
		$cmdStr = "mkdir $hintsPred_webDir";
		system $cmdStr;
		my $hintsPred_grailsDir = $grailsOut."/$projectID/autoAug/autoAugPred_hints";
		if (not(-d "$hintsPred_grailsDir")){
			print STDERR "AutoAug did not produce hints predictions!\n";
		}
		$cmdStr = "cp  $hintsPred_grailsDir/predictions/* $hintsPred_webDir; cp $hintsPred_grailsDir/gbrowse/* $hintsPred_webDir;";
		`$cmdStr`;
		$cmdStr = "cd $projectWebOutDir; tar -czvf hints_pred.tar.gz hints_pred;";
		`$cmdStr`;
		$cmdStr = "rm -r $hintsPred_webDir;";
		`$cmdStr`;	
	}

	## copy and pack UTR predictions - if they exist
	my $utrPredsExistFlag = 0;
	if(-d $grailsOut."/$projectID/autoAug/autoAugPred_utr"){
		$utrPredsExistFlag = 1;
	}else{
		print STDOUT "AutoAug did not produce utr predictions!\n";
	}
	if($utrPredsExistFlag==1){
		my $utrPred_webDir = $projectWebOutDir."/utr_pred";
		$cmdStr = "mkdir $utrPred_webDir";
		system $cmdStr;
		my $utrPred_grailsDir = $grailsOut."/$projectID/autoAug/autoAugPred_utr";
		$cmdStr = "cp  $utrPred_grailsDir/predictions/* $utrPred_webDir; cp $utrPred_grailsDir/gbrowse/* $utrPred_webDir;";
		`$cmdStr`;
		$cmdStr = "cd $projectWebOutDir; tar -czvf utr_pred.tar.gz utr_pred;";
		`$cmdStr`;
		$cmdStr = "rm -r $utrPred_webDir;";
		`$cmdStr`;	
	}	

	## copy and pack UTR hint predictions - if they exist
	my $utrHintsPredsExistFlag = 0;
	if(-d $grailsOut."/$projectID/autoAug/autoAugPred_hints_utr"){
		$utrHintsPredsExistFlag = 1;
	}else{
		print STDOUT "AutoAug did not produce utr and hint predictions!\n";
	}
	if($utrHintsPredsExistFlag==1){
		my $utrHintsPred_webDir = $projectWebOutDir."/hints_utr_pred";
		$cmdStr = "mkdir $utrHintsPred_webDir";
		system $cmdStr;
		my $utrHintsPred_grailsDir = $grailsOut."/$projectID/autoAug/autoAugPred_hints_utr";
		$cmdStr = "cp  $utrHintsPred_grailsDir/predictions/* $utrHintsPred_webDir; cp $utrHintsPred_grailsDir/gbrowse/* $utrHintsPred_webDir;";
		`$cmdStr`;
		$cmdStr = "cd $projectWebOutDir; tar -czvf hints_utr_pred.tar.gz hints_utr_pred;";
		`$cmdStr`;
		$cmdStr = "rm -r $utrHintsPred_webDir;";
		`$cmdStr`;	
	}	

	## copy log and error file
	my $errorFile = $grailsOut."/$projectID/AutoAug.err";
	my $logFile = $grailsOut."/$projectID/AutoAug.log";
	$cmdStr = "cp $errorFile $projectWebOutDir; cp $logFile $projectWebOutDir;";
	`$cmdStr`;
}else{
	## pack and copy augustus predictions
	if(-d $grailsOut."/$projectID/augustus"){
		$cmdStr = "cd $grailsOut/$projectID; tar -czvf augustus.tar.gz augustus;";
		`$cmdStr`;
		$cmdStr = "cp $grailsOut/$projectID/augustus.tar.gz $projectWebOutDir/augustus.tar.gz";
		`$cmdStr`;
	}
}

## create index-html page
my $head = $svnScripts."/webserver-results.head";
my $body = $svnScripts."/webserver-results.body";
my $tail = $svnScripts."/webserver-results.tail";
my $segmentFile1 = $grailsOut."/$projectID/segment1";
my $segmentFile2 = $grailsOut."/$projectID/segment2";


if($projectID =~ m/^t/){
	open(SEG1, ">", $segmentFile1) or die "Could not open file $segmentFile1!\n";
	print SEG1 "Training Results for job $projectID\n";
	close(SEG1) or die "Could not close file $segmentFile1!\n";

	open(SEG2, ">", $segmentFile2) or die "Could not open file $segmentFile2!\n";
	print SEG2 "<a href=\"index.html\" class=\"contentpagetitle\">Training results for job $projectID</a>\n</td>\n</tr>\n</table>\n";
	print SEG2 "<div class=\"main\" id=\"main\">\n<p>On this page, you find all relevant results to your AUGUSTUS training run $projectID for species $species, first submitted to our web server application on $submissionDate.</p>\n";
	print SEG2 "<h1>Files for download</h1>\n<table>\n<tr><td><b>Log-file</b></td><td><a href=\"AutoAug.log\">AutoAug.log</a></td></tr>\n<tr><td><b>Error-file</b></td><td><a href=\"AutoAug.err\">AutoAug.err</a></td></tr>\n";
	print SEG2 "<tr>\n<td><b>Species parameter archive</b>&nbsp;&nbsp;</td>\n<td><a href=\"parameters.tar.gz\">parameters.tar.gz</a></td>\n</tr>\n<tr>\n<td><b>Training genes</b>&nbsp;&nbsp;</td>\n";
	print SEG2 "<td><a href=\"training.gb.gz\">training.gb.gz</a></td>\n</tr>\n<tr>\n<td><b>Ab initio predictions</b></td>\n<td><a href=\"ab_initio.tar.gz\">ab_initio.tar.gz</a></td>\n</tr>\n";
	if($hintsPredsExistFlag==1){
	print SEG2 "<tr>\n<td><b>Predictions with hints</b></td>\n<td><a href=\"hints_pred.tar.gz\">hints_pred.tar.gz</a></td>\n</tr>\n";
	}
	if($utrPredsExistFlag==1){
	print SEG2 "<tr>\n<td><b>Predictions with UTR</b></td>\n<td><a href=\"utr_pred.tar.gz\">utr_pred.tar.gz</a></td>\n</tr>\n";
	}
	if($utrHintsPredsExistFlag==1){
	print SEG2 "<tr>\n<td><b>Predictions with hints and UTR &nbsp;&nbsp;</b></td>\n<td><a href=\"hints_utr_pred.tar.gz\">hints_utr_pred.tar.gz</a></td>\n</tr>\n";
	}
	print SEG2 "</table>\n<br><br>\n";
	close(SEG2) or die "Could not close file $segmentFile2!\n";
}else{
	open(SEG1, ">", $segmentFile1) or die "Could not open file $segmentFile1!\n";
	print SEG1 "Prediction Results for job $projectID\n";
	close(SEG1) or die "Could not close file $segmentFile1!\n";
	open(SEG2, ">", $segmentFile2) or die "Could not open file $segmentFile2!\n";
	print SEG2 "<a href=\"index.html\" class=\"contentpagetitle\">Prediction results for job $projectID</a>\n</td>\n</tr>\n</table>\n";
	print SEG2 "<div class=\"main\" id=\"main\">\n<p>On this page, you find all relevant results to your AUGUSTUS prediction run $projectID, first submitted to our web server application on $submissionDate.</p>\n";
	print SEG2 "<h1>Files for download</h1>\n<table>\n<tr><td><b>Prediction archive</b></td><td><a href=\"augustus.tar.gz\">augustus.tar.gz</a></td></tr>\n";
	print SEG2 "</table>\n<br><br>\n";
	close(SEG2) or die "Could not close file $segmentFile2!\n";
}

$cmdStr = "cat $head $segmentFile1 $body $segmentFile2 $tail > $projectWebOutDir/index.html";
`$cmdStr`;
}
