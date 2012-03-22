<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Help</title>         
    </head>
    <body>
<!-- ***** Start: Kopfbereich ********************************************// -->
<p class="unsichtbar">
  <a href="#inhalt" title="Directly to Contents">Directly to Contents</a>
</p>

<div id="navigation_oben">
  <a name="seitenanfang"></a>
  <table width="100%" border="0" cellpadding="0" cellspacing="1">
    <tr>
      <td nowrap="nowrap">
        <a href="http://www.uni-greifswald.de" target="_blank" class="mainleveltop_" >University of Greifswald</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.mnf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Faculty</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.math-inf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Institute</a><span class="mainleveltop_">&nbsp;|&nbsp;</span><a href="http://bioinf.uni-greifswald.de/" target="_blank" class="mainleveltop_">Bioinformatics Group</a>      </td>
    </tr>
  </table>
</div>

<div id="banner">
   <div id="banner_links">
       <a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Institut f&uuml;r Mathematik und Informatik"><img src="images/header.gif" alt="Directly to home" /> </a>
   </div>
   <div id="banner_mitte">
      <div id="bannertitel1">
        Bioinformatics Web Server at University of Greifswald
      </div>
      <div id="bannertitel2">
        Gene Prediction with AUGUSTUS <b><font color="#ffb22a" size=3>beta</font></b>
      </div>
   </div>
   <div id="banner_rechts">
     <a href="http://www.math-inf.uni-greifswald.de/mathe/index.php/geschichte-und-kultur/167" title="Voderberg-Doppelspirale">
     <img src="images/spirale.gif" align="left" />
     </a>
   </div>
</div>
<div id="wegweiser">
  Navigation for: &nbsp; &nbsp;<span class="breadcrumbs pathway">
    Training Tutorial
</span>

  <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Kopfbereich *********************************************// -->

<!-- ***** Start: Koerper ************************************************// -->
<div id="koerper">

  <div id="linke_spalte">
    <ul class="menu">
         <li><div id="linksMenuText">AUGUSTUS Web Server Navigation</div></li>
         <li><a href="index.gsp"><span>Introduction</span></a></li>
         <li><a href="about.gsp"><span>About AUGUSTUS</span></a></li>
         <li><a href="accuracy.gsp"><span>Accuracy</span></a></li>
         <li id="current"><a href="trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
         <li><g:link controller="training" action="create"><span>Submit Training</span></g:link></li>
         <li><a href="predictiontutorial.gsp"><span>Prediction Tutorial</span></a></li>
         <li><g:link controller="prediction" action="create"><span>Submit Prediction</span></g:link></li>
         <li><a href="help.gsp"><span>Help</span></a></li>
         <li><a href="datasets.gsp"><span>Datasets for Download</span></a></li>
         <li><a href="predictions_for_download.gsp"><span>Predictions for Download</span></a></li>
         <li><a href="references.gsp"><span>Links & References</span></a></li>
         <li><a href="impressum.gsp"><span>Impressum</span></a></li>
	  <li>&nbsp;</li>
         <li><div id="linksMenuText">Other AUGUSTUS Resources</div></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus">AUGUSTUS Wiki</a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Forum.Forum">AUGUSTUS Forum</a></li>
         <li><a href="http://bioinf.uni-greifswald.de/augustus/binaries/">Download AUGUSTUS</a></li>
         <li><a href="http://bioinf.uni-greifswald.de/augustus">Old AUGUSTUS gene prediction web server</a></li>
	  <li>&nbsp;</li>
         <li><div id="linksMenuText">Other Links</div></li>
         <li><a href="http://bioinf.uni-greifswald.de"><span>Bioinformatics Group Greifswald</span></a></li>
     </ul>
  </div>
 <div id="mittel_spalte">
<div class="inhalt" id="inhalt">
<!-- ***** Start: Content ************************************************// -->
<table class="contentpaneopen">
      <tr>
	<td class="contentheading" width="100%">
	  <a href="trainingtutorial.gsp" class="contentpagetitle">AUGUSTUS Training Tutorial</a>
        </td>
      </tr>
</table>
<p>This website explains step-by-step how to use the AUGUSTUS training web server application to train AUGUSTUS parameters for you individual species of interest. You find a similar tutorial on how to predict genes with pre-trained AUGUSTUS parameters <a href="predictiontutorial.gsp">here (click)</a>.</p>

<p>Functionalities of the AUGUSTUS training web server application are (with a single run):</p>

<p>
<ul>
<li>The generation of training gene structures (if no training gene structures are available)</li>
<li>Optimization of AUGUSTUS parameters according to the training gene structures</li>
<li>Prediction of genes in a supplied genome file with the newly optimized AUGUSTUS parameters. Genes will be predicted <i>ab initio</i> and with <i>hints</i> (the latter only if a cDNA file is provided).</li>
</ul>
</p>

<p><b>Before submitting a training job</b> for your species of interest, please check whether parameters have already been trained and have been made publicly available for your species at <a href="predictiontutorial.gsp#param_id">our species overview table</a></p>

<hr>
<br>

<div id="contents"><h1><a href="#contents">Contents</a></h1></div>

<p>
<a href="#job_submission">1 - Job Submission in general</a><br>
<a href="#finding_form">1.1 - Finding the training submission form</a><br>
<a href="#general_data">1.2 - Filling in general job data</a><br>
<a href="#email">1.2.1 - E-mail address</a><br>
<a href="#email_purpose">1.2.1.1 - What your e-mail address is used for</a><br>
<a href="#species">1.2.2 - Species name</a><br>
<a href="#species_form">1.2.2.1 - The species name format</a><br>
<a href="#species_purpose">1.2.2.2 - Purpose of entering a species name</a><br>
<a href="#genome_file">1.2.3 - Genome file</a><br>
<a href="#genome_file_format">1.2.3.1 - Genome file format</a></a><br>
<a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a><br>
<a href="#genome_file_purpose">1.2.3.3 - What the genome file is used for</a><br>
<a href="#optional_obligatory">1.3 - Optional obligatory fields</a><br>
<a href="#cDNA">1.3.1 - cDNA file</a><br>
<a href="#cDNA_format">1.3.1.1 - cDNA file format</a><br>
<a href="#cDNA_upload">1.3.1.2 - cDNA file upload options</a><br>
<a href="#cDNA_purpose">1.3.1.3 - What cDNA files are used for</a><br>
<a href="#protein_file">1.3.2 - Protein file</a><br>
<a href="#protein_file_format">1.3.2.1 - Protein file format</a><br>
<a href="#protein_file_upload">1.3.2.2 - Protein file upload options</a><br>
<a href="#protein_file_purpose">1.3.2.3 - What the protein file is used for</a><br>
<a href="#structure_file">1.3.3 - Training gene structure file</a><br>
<a href="#structure_file_format">1.3.3.1 - Training gene structure file format</a><br>
<a href="#structure_file_format_gbk">1.3.3.1.1 - Training gene structure file in genbank format</a><br>
<a href="#structure_file_format_gff">1.3.3.1.2 - Training gene structure file in gff format</a><br>
<a href="#structure_file_purpose">1.3.3.2 - What training gene structure files are used for</a><br>
<a href="#verification">1.4 - Verfification that you are a human</a><br>
<a href="#submitt">1.5 - The submitt button</a><br>
<a href="#exampledata">1.6 - Example data files</a><br><br>
<a href="#job_status">2 - What happens after submission</a><br>
<a href="#duplication">2.1 - Submission duplication</a><br>
<a href="#error">2.2 - Errors during training</a><br><br>
<a href="#results">3 - Training Results</a><br>
<a href="#results_autoAuglog">3.1 - The AutoAug.log file</a><br>
<a href="#results_autoAugerr">3.2 - The AutoAug.err file</a><br>
<a href="#results_parameters">3.3 - The parameters.tar.gz archive</a><br>
<a href="#traininggb">3.4 - The training.gb.gz file</a><br>
<a href="#predarchives">3.5 - Gene prediction archives</a><br><br>
<a href="#references">References</a><br>

</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="job_submission"><h1><a href="#job_submission">1 - Job Submission in general</a></h1></div>

<p>The pipeline invoked by submitting a job to the AUGUSTUS trainng web server application is complex. The following image gives an overview of processes that will possibly be invoked - depending on the input file combination:</p>


<p>
<table border="2" cellspacing="0" cellpadding="2">
<tr><td><img src="images/autoAug.jpg" alt="image of the autoAug.pl pipeline processes"></td></tr>
</table>
</p>

<p>The sharp-edge fields are input values of the submission form, the round edged fields are background processes/applied software tools/file formats/result files.</p>

<p>The input fields of the AUGUSTUS training web server application form are: E-mail, Species name, <b>Genome file</b>, <b>cDNA file</b>, <b>Protein file</b> and <b>Training gene structure file</b>. The actual processes invoked by job submission depend on the combination of input files.

<p>
<ul>
<li><b>{genome file, cDNA file}</b><br>In this case, the cDNA file is used to create training gene strucutures with <a href="http://pasa.sourceforge.net/">PASA</a> [<a href="#ref1">1</a>]. If cDNA end quality is sufficient, also a UTR training set will be created (this is currently the only possibility to train UTR parameters using this web server application). After parameters have been trained using the so created gene structure file, the cDNAs will additionally be used to create <i>hints</i>. Hints are extrinsic evidence for gene structures that are used during gene prediction <i>with hints</i>. The mapping tool <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="#ref2">2</a>]is used during hint generation. <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> is available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server, only! Finally, <a href="http://bioinf.uni-greifswald.de/augustus/">AUGUSTUS</a> [<a href="#ref4">4</a>] is used to predict genes in the genome file <i>ab initio</i> and <i>with hints</i>.</li>

<li><b>{genome file, protein file}</b><br>In this case, the protein file is used to create a training gene set using <a href="http://www.webscipio.org/">Scipio</a>. <a href="http://www.webscipio.org/">Scipio</a> [<a href="#ref3">3</a>] uses <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="#ref2">2</a>]. <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> is available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server, only! After parameter optimization, <a href="http://bioinf.uni-greifswald.de/augustus/">AUGUSTUS</a> [<a href="#ref4">4</a>] is used to predict genes in the genome sequence <i>ab initio</i>.</li>
					<li><b>{genome file, gene structure file}</b><br>In this case, the gene structure file is used as a training gene set. Gene structure files can be provided in two different formats: genbank format and gff format. If a genbank file is submitted, there is no dependency between training gene structure file and genome file. Parameters are then optimized based on the genbank training gene structure file. If a gff file is submitted, the gff must comply with the genome file entries. Training gene sequences are extraced from the genome file prior parameter optimization. Finally, <a href="http://bioinf.uni-greifswald.de/augustus/">AUGUSTUS</a> [<a href="#ref4">4</a>] is used for <i>ab initio</i> gene prediction. The submission of training gene structure file in combination with a genome file is <b><font color="f40b0b">open to all users!</font></b></li>
					<li><b>{genome file, cDNA file, protein file}</b><br>In this case, the protein file will be used to create a training gene set using <a href="http://www.webscipio.org/">Scipio</a> [<a href="#ref3">3</a>]. No UTR training set will be created. cDNA sequences will be used as evidence for gene prediction, only. Since both <a href="http://www.webscipio.org/">Scipio</a> and hint generation employ <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="#ref2">2</a>], this file combination is available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server, only! Finally, <a href="http://bioinf.uni-greifswald.de/augustus/">AUGUSTUS</a> [<a href="#ref4">4</a>] is used to predict genes in the genome file <i>ab initio</i> and <i>with hints</i>.</li>
					<li><b>{genome file, cDNA file, gene structure file}</b><br>In this case, the gene structure file is used as a training gene set. cDNA sequences will be used as evidence for prediction in form of hints that are generated with the help of <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a>. Since hint generation employs <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="#ref2">2</a>], this file combination is available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server, only! Finally, <a href="http://bioinf.uni-greifswald.de/augustus/">AUGUSTUS</a> [<a href="#ref4">4</a>] is used to predict genes in the genome file <i>ab initio</i> and <i>with hints</i>.</li>
				</ul>
</p>

<p>In the following, you find detailed instructions for submitting an AUGUSTUS training job.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="finding_form"><h2><a href="#finding_form">1.1 - Finding the training submission form</a></h2></div>

<p>You find the AUGUSTUS training submission form by clicking on the following link in the left side navigation bar:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/submitt-link.jpg" alt="image of submission link"></td></tr>
</table>
</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="general_data"><h2><a href="#general_data">1.2 - Filling in general job data</a></h2></div>

<p>This section describes all fields that must be filled in for every job submission, i.e. fields that are obligatory.</p> 

<div id="email"><h3><a href="#email">1.2.1 - E-mail address</a></h3></div> 

<p>At first, you have to enter a <b>valid e-mail address</b>:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/email-training.jpg" alt="image of e-mail address field"></td></tr>
</table>
</p>


<p>The e-mail address is an obligatory field. Currently, there is no option for submitting AUGUSTUS training jobs anonymously.</p>

<div id="email_purpose"><h4><a href="#email_purpose">1.2.1.1 - What your e-mail address is used for</a></h4></div>

<p> We save use your e-mail address for the following purposes:</p>

<p>
<ul>
<li>Confirming your job submission</li>
<li>Confirming successful file upload (for large files)</li>
<li>Sending you link to the results of your job (you will never get the link to the results if you do not enter a <b>valid</b> address to which you have access)</li>
<li>Informing you about any problems that might occur during your particular AUGUSTUS training job</li>
</ul>
</p>

<p>We do <b>not</b> use your e-mail address to send you any <i>spam</i>, i.e. about web service updates. We do not share your e-mail addresses with any third parties.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="species"><h3><a href="#species">1.2.2 - Species name</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/species-name.jpg" alt="image of species name field"></td></tr>
</table>
</p>


<div id="species_form"><h3><a href="#species_form">1.2.2.1 - The species name format</a></h3></div>

<p>The species name must be a string of at most 30 characters length. It must be unique in our database and it may <b>not</b> contain spaces.</p>

<div id="species_purpose"><h3><a href="#species_purpose">1.2.2.2 - Purpose of entering a species name</a></h3></div>

<p>The species name is the name of the species for whose genome you want to train AUGUSTUS. The species name is an obligatory parameter that AUGUSTUS needs in order to find the correct parameters that shall be applied for predicting genes in a specified genomic sequence. Considering that AUGUSTUS training is such a time consuming process, our objective is to know the names of species for which AUGUSTUS was trained in order to make the trained parameters available to the public so that others who are interested in the same species as you do not have to rerun the training process.</p>

<p>However, if you do not want to reveal the true species name, you may use any other string shorter than 30 characters as a species name. Species names must be unique on our system, i.e. if the string of your choice is already existing in our system, you will get a message that you have to choose another species name.</p>

<p>We are <b>not</b> redistributing the original sequence data that you submitted to our web server application. But we are redistributing the trained parameters and the species name and any other kind of results that your computation may have produced.</p><p><b>Example:</b></p><p>If person 1 submits a sequence data set for training and names it <i>hypothetical_species</i>, and a second person tries later to train AUGUSTUS with exactly the same sequence files (new data upload) and names the species <i>some_other_name</i>, the second person will be redirected to the results of the original training run of person 1, and the species name <i>hypothetical_species</i> will be publicly readable to person 2.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="genome_file"><h3><a href="#genome_file">1.2.3 - Genome file</a></h3></div>

<p>The genome file is an obligatory file for training AUGUSTUS.</p>

<div id="genome_file_format"><h4><a href="#genome_file_format">1.2.3.1 - Genome file format</a></h4></div>

The genome file must contain the genome sequence in (multiple) fasta format. Every <i>unique</i> header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver application. White spaces in the sequence header might cause problems if the first word after the leading character <b>></b> is identical for several fasta entries. We generally recommend short, unique, non-white-space containing fasta headers.<br><br>
            <b>Correct file format example:</b>

<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1">            <pre class="example">
>Chr.1
CCTCCTCCTGTTTTTCCCTCAATACAACCTCATTGGATTATTCAATTCAC
CATCCTGCCCTTGTTCCTTCCATTATACAGCTGTCTTTGCCCTCTCCTTC
TCTCGCTGGACTGTTCACCAACTCTCAGCCCGCGATCCCAATTTCCAGAC
AACCCATCTTATCAGCTTGGCCACGGCCTCGACCCGAACAGACCGGCGTC
CAGCGAGAAGAGCGTCGCCTCGACGCCTCTGCTTGACCGCACCTTGATGC
TCAAGACTTATCGCGATGCCAAGAAGCGTCTCATCATGTTCGACTACGA
>Chr.2
CGAAACGGGCACCTATACAACGATTGAAACCATTATTCAAGCTCAGCAAG
CGTCTATGCTAGCGGTTATTGCGAGCACTTCAGCGGTTGCTACTACGACT
ACTACTTGATAAATGAAACGGCTATAAAAGAGGCTGGGGCAAAAGTATGT
TAGTTGAAGGGTGACCTGAACGATGAATCGGTCGAATTTTTTATTGGCAG
AGGGAAGGTAGGTTTACTCAATTTAGTTACTTCTAGCCGTTGATTGGAGG
AGCGCAAGCGACGAGGAGGCTCATCGGCCGCCCGCGGAAAGCGTAGTCT
TACACGGAAATCAACGGCGGTGTCATAAGCGAG
>Chr.3
.....
            </pre></font></td></tr>
</table>
</p>

<p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for genome file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>.</p>

<div id="genome_file_upload"><h4><a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a></h4></div>

<p>The AUGUSTUS training web server application offers two possiblities for transferring the genome file to the server: <i>Upload a file</i> and <i>specify a web link to file</i>.</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-genome.jpg" alt="image of genome submission field"></td></tr>
</table>
</p>

             <p>
             <ul>
                <li>For <b>small files</b>, please click on the <i>Choose File</i> or <i>Browse</i>-button and select a file on your harddrive.<br>If you experience a <i>Connection timeout</i> (because your file was too large for this type of upload - the size is browser dependent), please use the option for large files!</li>
                <li><b>Large files</b> can be retrieved from a <b>public</b> web link. Deposit your sequence file at a http or ftp server and specify the valid URL to your sequence file in the training submission form. Our server will fetch the file from the given address upon job submission. (File size limit: currently 1 GB. Please contact us in case you want to upload a bigger genome file.) You will be notified by e-mail when the file upload from web-link is finished (i.e. you can delete the file from the public server after you received that e-mail).</li>
            </ul>
	     </p>
	    <p>
            <b>You cannot do both at the same time!</b> You must <b>either</b> select a file on your harddrive <b>or</b> give a web link!</p>

<div id="genome_file_purpose"><h4><a href="#genome_file_purpose">1.2.3.3 - What the genome file is used for</a></h4></div>

<p>The genome file can be used for two purposes:</p>

<p>
<ul>
<li>Training gene structure generation (if no training gene structure file was submitted)</li>
<li>Target sequence for gene prediction with the new AUGUSTUS parameters (always)</li>
</ul>
</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="optional_obligatory"><h2><a href="#optional_obligatory">1.3 - Optional obligatory fields</a></h2></div>

<p>This section describes a number of fields from which <b>at least one</b>  must be specified for training AUGUSTUS.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="cDNA"><h3><a href="#cDNA">1.3.1 - cDNA file</a></h3></div>

<p>This feature is only for personal, academic, and non-profit use as this is required by the BLAT license.</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-cDNA.jpg" alt="image of cDNA submission field"></td></tr>
</table>
</p>


<div id="cDNA_format"><h3><a href="#cDNA_format">1.3.1.1 - cDNA file format</a></h3></div>

<p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. An example for correct cDNA file format is given at <a href="#genome_file_format">1.2.3.1 - Genome file format</a>.</p>

<p>It is currently possible to submitt assembled RNA-seq transcripts instead of or mixed with ESTs as a cDNA/EST file. However, you should be aware that the success of creating training genes with assembled RNA-seq data depends very much on the assembly quality, and that RNA-seq files are often much bigger than EST or cDNA files, which increases runtime of a training job. In order to keep runtime of your training job as low as possible, you should remove all assembled RNA-seq transcripts from your file that do not map to the submitted genome sequence. (In principle, this holds true for EST and cDNA files, too, but there, the problem is not as pronounced due to a smaller number of sequences.)</p>

<p>It is currently not allowed to upload RNA-seq raw sequences. (We filter for the average length of cDNA fasta entries and may reject the entire training job in case the sequences are on average too short, i.e. shorter than 400 bp.)</p>

<p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for cDNA file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>. The maximal supported file size is 1 GB.</p>

<div id="cDNA_upload"><h3><a href="#cDNA_upload">1.3.1.2 - cDNA file upload options</a></h3></div>

<p>There are two options for cDNA file upload: upload from your local harddrive, or upload from a public http or ftp server. Please see <a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a> for a more detailed description of upload options.</p>

<div id="cDNA_purpose"><h3><a href="#cDNA_purpose">1.3.1.3 - What cDNA files are used for</a></h3></div>

<p>The cDNA file can be used for two purposes by the AUGUSTUS training web server application:</p>

<p>
<ul>
<li>for generating training gene structures (if no protein file and no training gene structure file are specified)</li>
<li>for generating extrinsic evidence for gene structures in the gene prediction process (always)</li>
</ul>
</p>


<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="protein_file"><h3><a href="#protein_file">1.3.2 - Protein file</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-protein.jpg" alt="image of protein submission field"></td></tr>
</table>
</p>

<div id="protein_file_format"><h4><a href="#protein_file_format">1.3.2.1 - Protein file format</a></h4></div>

            <p>The protein file is a multiple fasta file that contains protein sequences, e.g. from a closely related species. Allowed sequence characters: <b>A a R r N n D d C c E e Q q G g H h I i L l K k M m F f P p S s T t W w Y y V v B b Z z J j X x</b>. Empty lines are not allowed but will simply be removed from the file by the webserver application. <br><br></p>
            <p><b>Correct file format example:</b>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1"><pre class="example">
>protein1
maaaafgqlnleepppiwgsrsvdcfekleqigegtygqvymakeiktgeivalkkirmd
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqvlhrdikgsnllidnegnlkl
adfglarsyshdhtgnltnrvitlwyrppelllgatkygp
>protein2
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqv
>protein3
...
            </pre></font></td></tr>
</table>
            </p>

<div id="protein_file_upload"><h4><a href="#protein_file_upload">1.3.2.2 - Protein file upload options</a></h4></div>

<p>There are two options for protein file upload: upload from your local harddrive, or upload from a public http or ftp server. Please see <a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a> for a more detailed description of upload options.</p>

<div id="protein_file_purpose"><h4><a href="#protein_file_purpose">1.3.2.3 - What the protein file is used for</a></h4></div>

<p>The protein file is used for generating training gene structures by mapping the protein sequences against the supplied genome sequence. You may e.g. upload a file with protein sequences from a closely related species in order to obtain training genes and AUGUSTUS parameters for the new sspecies whose genome you uploaded.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="structure_file"><h3><a href="#structure_file">1.3.3 - Training gene structure file</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/structure-file.jpg" alt="image of structure file upload field"></td></tr>
</table>
</p>

<div id="structure_file_format"><h4><a href="#structure_file_format">1.3.3.1 - Training gene structure file format</a></h4></div>

<p>Training gene structure files can be submitted in two different formats: Genbank format or gff format.</p>

<div id="structure_file_format_gbk"><h5><a href="#structure_file_format_gbk">1.3.3.1.1 - Training gene structure file in genbank format</a></h5></div>

<p>Gene structures in genbank format must contain the coding sequence parts and flanking regions. Flanking regions are important because AUGUSTUS is supposed to differentiate between genes and intergenic regions. The length of flanking regions depends on the length of genes in the target genome. In our pipeline, flanking regions are set to the average gene length (exceptionally applying the extreme limits between 1000 and 10000 nt). It is very important to make sure that the flanking regions do not contain any other protein coding gene parts, i.e. we recommend to trim flanking regions in a way that will exclude other CDS parts.</p>

<p><b>Correct file format example (condensed view, the three dots represent further lines of sequences):</b>

<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1">            <pre class="example">
LOCUS       Chr.1_1-159458   159458 bp  DNA
FEATURES             Location/Qualifiers
     source          1..159458
     CDS             complement(join(2421..2655,3858..4005,4080..4235,5569..5857
                     ,10316..10534,155240..155458))
                     /gene="1474336"
BASE COUNT     49195 a   29117 c  28985 g   49950 t   2211 n
ORIGIN
        1 aaaatacatc acaatacatt taattcactt tccatcatcg agattaacga aaattattta
       61 aaatatcgaa gatgaaaata tcctcaagat gatactgaac ggctaagaaa aatacatcac
      121 acaactttaa ttcattttcc atcatcgaga ttaacgaaaa gaaaaaattt taactcccta
...
   159301 atacgccacc aggtatttcg cctgattgtt cctcgaatat cttctctctc tctatatata
   159361 tatatattac ttggcacgat aatcgtcgaa tcgttattta taaattgctt catctatcgc
   159421 gatatttttg caacaactct cgcttttctc tttccatt
//
LOCUS       Chr.1_313992-323129   9138 bp  DNA
FEATURES             Location/Qualifiers
     source          1..9138
     CDS             join(4001..4048,4989..5138)
                     /gene="194551"
BASE COUNT     2829 a   1502 c  1750 g   2948 t   109 n
ORIGIN
        1 ttttccttct ttcttttttt tttatttaca ttaatgagaa ttttcgcaaa tatttcatcg
       61 ctgccatcct tttttttcct cgacgtcaat cacgcgacac atttgttaga gaaatggatt
      121 ttaatcttga aaaaagaaaa atacaaatgc caacgcattt caaatccttt cctattatta
...
     9001 tcaacgaaac aaataattgc ttcacaaaat atcgcacgta acaacaatat agacttcaat
     9061 attcaacaat tcttttcctt tatacacaaa gatacacaaa atataaaagt tttaatactt
     9121 caacttcaac gaaacagg
//
            </pre></font></td></tr>
</table>
</p>


<p>If you want to train UTRs, you have to additionally incorporate mRNA information in your genbank file.</p>

<p><b>Correct file format example (including UTR training):</b>

<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1">            <pre class="example">
LOCUS       scf7180001240730_g20   526 bp  DNA
FEATURES             Location/Qualifiers
     source          1..526
     mRNA            99..125
     CDS             99..99
BASE COUNT     164 a   99 c  68 g   195 t
ORIGIN
        1 gtgacggagc ccaaggacga gcccgtgccc tcagagccca cgtccgacgt gaggcccgcg
       61 ccagcgcccc tcccgccgcc cgtcgcagcc actgcttaga ctttactaat ataaacattg
      121 aaaatatttt gtgttttatt tccaatcatt gaattataat cctattataa tataactaac
      181 attcgtaatt ttacaaaata actatgcaaa ttattttgta ttttcgtttt aaattatact
      241 tttcatataa atttctacaa atcttattca agaccataag tatccgctcg ctctacttcg
      301 ggcatttcct ttatttatat cttatttgac ttattttgat tatttaggct tatgttttcg
      361 atactattga aaacagaaaa taatttcata taattaataa tatattttca attaatatat
      421 ttaacaaata tttgtatagt tcaagcggac aaatccgttc ccatagtatt tatataaatt
      481 ttaatttaga gtaataacag tttgctgtat tgttgtagtc aaatac
//
LOCUS       scf7180001240751_g30   876 bp  DNA
FEATURES             Location/Qualifiers
     source          1..876
     mRNA            complement(401..777)
     CDS             complement(777..777)
BASE COUNT     300 a   136 c  116 g   324 t
ORIGIN
        1 aatgtaggaa aatgaaatat ttatttaaat tgttattatc acttcttcgc tctagtgtct
       61 tggcaaagcg cggcgttgag ttcagcctct cacacgcaat gcctccagaa ttcggcgaaa
      121 tgtgggggac agagtgtatt aacactaagt tccctcagcc acgactggtg aaattatata
      181 ttcagtttgt atactattac tcatgcaaac acttcatcat actttcactc aatcagtaaa
      241 gcataatatt ttatttaata ttgtttatca atactatttc cttgttgtta aatattattt
      301 tatttattat attaaattaa aatgtcaaaa ttaaaagtag gtgatgattt attactatct
      361 tttctatcca agaaaaaaaa gacacactga aacaattgta atttttgtta tgtttttatt
      421 acttaatatt attataaaaa tttgtaaata cgaaataaaa tagatagacg taataatatt
      481 tatttgttag ttaataataa taatgataat tacgaaagat acaagaaata tgcataaatg
      541 agtgttatat tatgtatttt atgagaatat aaatataaaa actgtcattg attatatttt
      601 ctaaatactt tcattttatg gcttgctggc ttttcaattt ccttatgttt cagcttttca
      661 ctcaatagag cgaaaccttc atcgacatgt aagccaatag aacaattaca aactaacttt
      721 attacatcag tcttttcatt tctttaagct tcaggcaaat atcatctaaa tgcctttcaa
      781 ctcgctacta acatcgcgtc gttatataaa tcagtgtata cggaattaaa cctgtcatgt
      841 ctcttgcaag acgtgtctgc tgttgtcacg cacaca
//
            </pre></font></td></tr>
</table>
</p>

<div id="structure_file_format_gff"><h5><a href="#structure_file_format_gff">1.3.3.1.2 - Training gene structure file in gff format</a></h5></div>

<p>Training gene structure in gff format must comply with the fasta entry names of the genome file.</p>

<p>In general, gff format must contain the following columns (The columns are separated by <b>tabulators</b>):</p>

<p>
            <OL TYPE="1">
            <li>The <b>sequence names</b> must be found in the fasta headers of sequences in the genome file. 
	    <li>The <b>source</b> tells with which software/process the gene structure was generated (you can fill in whatever you like).
            <li>The <b>feature</b> may for AUGUSTUS training be CDS, 5'-UTR or 3'-UTR. 
            <li><b>Start</b> is the beginning position of the line's feature, counting the first position of a sequence as position 1.
            <li><b>Stop</b> position, must be at least as large as start position.
            <li>The <b>score</b> must be a number but the number is irrelevant to our web server applications.
	    <li>The <b>strand</b> denotes whether the gene is located on the forward (+) or on the reverse (-) strand.
            <li><b>Frame</b> is the reading frame, can be denoted as '.' if unknown or irrelevant. For exonpart and exon this is as defined as follows: On the forward strand it is the number of bases after (begin position 1) until the next codon boundary comes (0, 1 or 2). On the reverse strand it is the number of bases before (end position + 1) the next codon boundary comes (0, 1 or 2).
            <li><b>Attribute</b> contains a <i>transcript identifier</i>. All gff-entries belonging to one transcript must contain the same transcript identifier in the last column.
            </OL>
</p>

<p><b>Correct file format example (without UTR):</b>

<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1"><pre class="example">
Chr.1	mySource	CDS	1767	1846	1.000	-	0	transcript_id "1597_1"
Chr.1	mySource	CDS	1666	1709	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	1486	1605	1.000	-	2	transcript_id "1597_1"
Chr.1	mySource	CDS	1367	1427	1.000	-	2	transcript_id "1597_1"
Chr.1	mySource	CDS	1266	1319	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	1145	1181	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	847	1047	1.000	-	0	transcript_id "1597_1"
Chr.2	mySource	CDS	9471	9532	1.000	+	0	transcript_id "1399_2"
Chr.2	mySource	CDS	9591	9832	1.000	+	1	transcript_id "1399_2"
Chr.2	mySource	CDS	9885	10307	1.000	+	2	transcript_id "1399_2"
Chr.2	mySource	CDS	10358	10507	1.000	+	2	transcript_id "1399_2"
Chr.2	mySource	CDS	10564	10643	1.000	+	2	transcript_id "1399_2"
</pre></font></td></tr>
</table>
    </p>

<p><b>Correct file format example (with UTR):</b>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1"><pre class="example">
Chr.1	mySource	5'-UTR	277153	277220	45	+	.	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	277221	277238	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278100	278213	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278977	279169	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279630	279648	0.94	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279734	279768	0.94	+	1	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	280307	280344	1	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	3'-UTR	280345	280405	78	+	.	transcript_id "g22472.t1"; gene_id "g22472";
</pre></font></td></tr>
</table>
</p>

<div id="structure_file_purpose"><h4><a href="#structure_file_purpose">1.3.3.2 - What training gene structure files are used for</a></h4></div>

<p>While the cDNA and protein file upload options fulfill the purpose of training gene generation in context with a genome sequence file, the gene structure file upload form provides a flexible interface to users who generated training gene structures externally instead of using our webserver application for this purpose. Supplied training gene structures are directly used for optimizing species-specific AUGUSTUS parameters.</p>

<p>If the training gene structure file is in gff-format, sequence names in the gff file have to match fasta headers in the genome file. If the training gene structure file is in genbank format, no direct dependency between training gene strcuture file and genome file exists, i.e. you could submitt training genes that were generated from a different chromosome than the genome file that you submitt.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="verification"><h2><a href="#verification">1.4 - Verification that you are a human</a></h2></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/verification.jpg" alt="image of verification field"></td></tr>
</table>
</p>

<p>Trying to avoid abuse of our web server application through bots, we implemented a <i>captcha</i>. The captcha is an image that contains a string. You have to type the string from the image into the field next to the image.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="submitt"><h2><a href="#submitt">1.5 - The submitt button</a></h2></div>


<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/train-submitt.jpg" alt="image of submission button"></td></tr>
</table>
</p>

<p>After filling out the appropriate fields in the submission form, you have to click on the button that says "Start Training" at the bottom of the page. It might take a while until you are redirected to the status page of your job. The reason is that we are checking various file formats prior job acceptance, and that the transfer of files from your local harddrive to our server might take a while. Please be patient and wait until you are redirected to the status page! Do not click the button more than once (it won't do any harm but it also doesn't speed up anything).</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="exampledata"><h2><a href="#exampledata">1.6 - Example data files</a></h2></div>

<p>In the following, we provide some correctly formatted, compatible example data files:</p>

<p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/chr1to3.fa">http://bioinf.uni-greifswald.de/trainaugustus/examples/chr1to3.fa</a> - This file may be used as a <b>Genome file</b>. It contains the first three chromosomes of <i>Mus musculus</i> from GenBank (modified headers).</p>

<p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/estsChr1to3.fa">http://bioinf.uni-greifswald.de/trainaugustus/examples/estsChr1to3.fa</a> - This file may be used as a <b>cDNA file</b>. It contains EST sequences of <i>Mus musculus</i> from GenBank (modified headers).</p>

<p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/rattusProteinsChr1to3.fa">http://bioinf.uni-greifswald.de/trainaugustus/examples/rattusProteinsChr1to3.fa</a> - This file may be used as a <b>Protein file</b>. It contains proteins of <i>Rattus norvegicus</i> from GenBank (modified headers) that map to chromosomes 1 to 3.</p>

<p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/traingenes.gb">http://bioinf.uni-greifswald.de/trainaugustus/examples/traingenes.gb</a> - This file contains 221 training genes for <i>Mus musculus</i> in GenBank format (including flanking regions).</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="job_status"><h1><a href="#job_status">2 - What happens after submission</a></h1></div>

<p>After you click the "Start Training" button, the web server application first validates whether the combination of your input fields is generally correct. If you did anything wrong, you will be redirected to the training submission form and an error message will be displayed at the top of the page.</p>

<p>If all fields were filled in correctly, the application is actually initiated. You will receive an e-mail that confirms your job submission and that contains a link to the job status page.</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-job-status.jpg" alt="image of job status page"></td></tr>
</table>
</p>


<p>In the beginning, the status page will display that your job has been <b>submitted</b>. This means, the web server application is currently uploading your files and validating file formats. After a while, the status will change to <b>waiting for execution</b>. This means that all file formats have been confirmed and an actually AUGUSTUS training job has been submitted to our grid engine, but the job is still pending in the queue. Depending on waiting queue length, this status may persist for a while. Please contact us in case you job is pending for more than one month. Later, the job status will change to <b>computing</b>. This means the job is currently computing. When the page displays <b>finished</b>, all computations have been finished and a website with your job's results has been generated.</p>

<p>You will receive an e-mail with the link to the results of your job. This link is publicly accessible but hard to guess.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="duplication"><h2><a href="#duplication">2.1 - Submission duplication</a></h2></div>

<p>Since training AUGUSTUS is a very resource consuming process, we try to avoid data duplication. In case you or somebody else tries to submitt exactly the same input file combination more than once, the duplicated job will be deleted and the submitter of the redundant job will receive an e-mail with the link to the previously submitted job status and with a link to the results page (which may still be empty in case the duplicated job has not finished computing, yet.)</p>

<p>You web browser will be redirected to a page with the following content if data duplication occured:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/job-does-not-exist-training.jpg" alt="image of non existing job status"></td></tr>
</table>
</p>


<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="error"><h2><a href="#error">2.2 - Errors during training</a></h2></div>

<p>You should automatically receive an e-mail in case an error occurs during the AUGUSTUS training process. The admin of this server is also notified by e-mail about errors. We will get in touch with you, again, after we figured out what caused the error.</p>

<p>Since the web server application is currently in beta testing phase, completely unexpected errors might still occur. Therefore we ask for you help on reporting any weird webserver application behaviour to augustus-web@uni-greifswald.de. Please include which actions from your side exactly caused the error, and also copy the Grails exception message into your e-mail in case such a message was displayed.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="results"><h1><a href="#results">3 - Training Results</a></h1></div>

<p>After job computations have finished, you will receive an e-mail with a link to the results of your submission. The linked web page may look similar to this:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-results-example.jpg" alt="image of results example"></td></tr>
</table>
</p>


<p>This page will should contain the files <b>AutoAug.log</b> and <b>AutoAug.err</b>. If no other files are displayed, you will find the reason for this behaviour in one of those files (e.g. if the creation of a sufficient number of training genes was not possible using the data that you supplied).</p>

<p>If training was successful, at least the files <b>parameters.tar.gz</b> and <b>ab_initio.tar.gz</b> should be displayed.</p>

<p>If no training gene structure file was supplied, <b>training.gb.gz</b> should be linked.</p>

<p>If a cDNA file was submitted, <b>hints_pred.tar.gz</b> should be shown.</p>

<p>To save files on your local system, click the file link with the right mouse button and select "Save Link As..." (or similar).</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="results_autoAuglog"><h2><a href="#results_autoAuglog">3.1 - The AutoAug.log file</a></h2></div>

<p>This file contains all logging messages of the AUGUSTUS training and prediction processed invoked by the AUGUSTUS training web server application. It is a plain text file, i.e. you should be able to open it with any text editor or even your web browser. In your own interest, you should check the AutoAug.log file for the number of training genes that were generated from you input files (except if you submitted a training gene structure file), for the number UTR training genes, and for gene prediction accuracy.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="results_autoAugerr"><h2><a href="#results_autoAugerr">3.2 - The AutoAug.err file</a></h2></div>

<p>This file contains all error messages of the AUGUSTUS training and prediction processed invoked by the AUGUSTUS training web server application. It is a plain text file, i.e. you should be able to open it with any text editor or even your web browser. In your own interest, you should check the AutoAug.err file. If this file is not empty, something did go wrong and you probably shouldn't trust the results blindly - in case any were produced at all.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="results_parameters"><h2><a href="#results_parameters">3.3 - The parameters.tar.gz archive</a></h2></div>

<p>This file contains the parameters that were trained for your species of interest based on the submitted input files. The AUGUSTUS parameter archive always contains the following files for your <i>species</i> (the term species will be replaced by the species name that you entered when submitting a training job):</p>

<p>
<ul>
<li><i>species</i>_parameters.cfg
<li><i>species</i>_metapars.cfg
<li><i>species</i>_metapars.utr.cfg
<li><i>species</i>_exon_probs.pbl.withoutCRF
<li><i>species</i>_exon_probs.pbl
<li><i>species</i>_weightmatrix.txt
<li><i>species</i>_intron_probs.pbl
<li><i>species</i>_intron_probs.pbl.withoutCRF
<li><i>species</i>_igenic_probs.pbl
<li><i>species</i>_igenic_probs.pbl.withoutCRF
</ul>
</p>

<p>If you want to use the parameters for predicting genes with a local AUGUSTUS installation, place the species parameter folder into your <tt>$AUGUSTUS_CONFIG_PATH</tt> species folder. Be aware that in general, our web server application does not optimize the CRF parameters. Therefore, do not run AUGUSTUS with CRF-flag with these parameters. UTR parameters have only been optimized for your species if this website contains an archive utr_pred.tar.gz. It does not make sense to locally predict genes with UTR if the UTR parameters were not optimized for your species.</p>

<p>If you want to run further predictions for your species of interest in other genomic sequences or with other hints on our web server application, you do not need to save the parameter archive locally for re-upload. Just remember the ID of your job and specify it on the prediction interface. You find your job's ID in the headline of this web site at "Training results for job ...".</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="traininggb"><h2><a href="#traininggb">3.4 - The training.gb.gz file</a></h2></div>

<p>
This file contains generated training genes that were used to optimize AUGUSTUS parameters in genbank format. The file is compressed with gzip. To unpack it, type <tt>gunzip training.gb.gz</tt> in your shell. You find more information about gzip at <a href="http://www.gzip.org/">gzip.org</a>.
</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="predarchives"><h2><a href="#predarchives">3.5 - Gene prediction archives</a></h2></div>

<p>
The files <b>ab_inito.tar.gz</b> and <b>pred_hints.tar.gz</b> are gene prediction archives. Their content depends on the input file combination.
</p>

<h4>Files that are always contained in gene prediction archives:</h4>

<p>
<ul>
<li>*.gff - gene predictions in gff format</li>
</ul>
</p>

<h4>Files that may optionally be contained in gene prediction archives:</h4>
<p>
<ul>
<li>*.gtf - gene predictions in gtf format
<li>*.aa - gene predictions as protein fasta sequences
<li>*.codingseq - gene predictions as CDS DNA fasta sequences
<li>*.cdsexons - predicted exons in DNA fasta sequences
<li>*.mrna - predicted mRNA sequences (with UTRs) in DNA fasta sequences
<li>*.gbrowse - gene prediction track for the GBrowse genome browser
</ul>
</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="references"><h1><a href="#references">References</a></h1></div>

<div id="ref1"><p>
[1] Haas, B.J., Delcher, A.L., Mount, S.M., Wortman, J.R., Smith Jr, R.K., Jr., Hannick, L.I., Maiti, R., Ronning, C.M., Rusch, D.B., Town, C.D. et al. (2003) Improving the Arabidopsis genome annotation using maximal transcript alignment assemblies. Nucleic Acids Res, 31, 5654-5666.
</p></div>

<div id="ref2"><p>
[2] Kent, W.J. (2002) BLAT—The BLAST-Like Alignment Tool. Genome Res, 12, 656-664.
</p></div>

<div id="ref3"><p>
[3] Keller, O., Odronitz, F., Stanke, M., Kollmar, M., Waack, S. (2008) Scipio: Using protein sequences to determine the precise exon/intron structures of genes and their orthologs in closely related species. BMC Bioinformatics 9, 278.
</p></div>

<div id="ref4"><p>
[4] Stanke, M., Diekhans, M., Baertsch, R., Haussler, D. (2008) Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, doi: 10.1093/bioinformatics/btn013
</p></div>


<!-- ***** End: Content ************************************************// -->
</div>
</div>
 </div>
</div>
<div id="rechte_spalte">
    <div class="linien_div">
      <h5 class="ueberschrift_spezial">CONTACT</h5>
      <strong>Institute for Mathematics und Computer Sciences</strong><br/>
      <strong>Bioinformatics Group</strong><br />
      Walther-Rathenau-Stra&szlig;e 47<br />
      17487 Greifswald<br />
      Germany<br />
      Tel.: +49 (0)3834 86 - 46 24<br/>
      Fax:  +49 (0)3834 86 - 46 40<br /><br />
      <a href="mailto:augustus-web@uni-greifswald.de" title="E-Mail augustus-web@uni-greifswald.de, opens the standard mail program">augustus-web@uni-greifswald.de</a>
    </div>

    <div class="beendeFluss"></div>
</div>
</div>

<!-- ***** Ende: Koerper *************************************************// -->
<!-- ***** Start: Fuss ***************************************************// -->
<div id="fuss">
  <div id="fuss_links"><p class="copyright">&copy; 2011 University of Greifswald</p></div>
  <div id="fuss_mitte">
  <div class="bannergroup">
 </div>
 </div>
 <div id="fuss_rechts" >
  <ul>
   <li>
    <a href="#seitenanfang">
     <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
     Top of page
    </a>
   </li>
  </ul>
 </div>
 <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Fuss ***************************************************// -->

    </body>
</html>
