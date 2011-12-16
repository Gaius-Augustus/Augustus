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
         <li><a href="index.gsp"><span>Introduction</span></a></li>
         <li><a href="trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
         <li><g:link controller="training" action="create"><span>Submit Training</span></g:link></li>
         <li id="current"><a href="predictiontutorial.gsp"><span>Prediction Tutorial</span></a></li>
         <li><g:link controller="prediction" action="create"><span>Submit Prediction</span></g:link></li>
         <li><a href="help.gsp"><span>Help</span></a></li>
         <li><a href="references.gsp"><span>Links & References</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de"><span>Bioinformatics Group</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/impressum.html"><span>Impressum</span></a></li>
     </ul>
  </div>
 <div id="mittel_spalte">
<div class="inhalt" id="inhalt">
<!-- ***** Start: Content ************************************************// -->
<table class="contentpaneopen">
      <tr>
	<td class="contentheading" width="100%">
	  <a href="trainingtutorial.gsp" class="contentpagetitle">AUGUSTUS Prediction Tutorial</a>
        </td>
      </tr>
</table>
<p>This website explains step-by-step how to use the AUGUSTUS prediction web server application to predict genes in a genomic sequence. You find a similar tutorial on how to train AUGUSTUS parameters <a href="trainingtutorial.gsp">here (click)</a>.</p>

<p>Functionalities of the AUGUSTUS prediction web server application are (with a single run):</p>

<p>
<ul>
<li>Generation of hints (if a cDNA file is supplied)</li>
<li>Prediction of genes in a genome sequence file using the supplied parameters. Genes will be predicted <i>ab initio</i> and with <i>hints</i> (the latter only if a cDNA and/or hint file is provided).</li>
</ul>
</p>

<hr>
<br>

<div id="contents"><h1><a href="#contents">Contents</a></h1></div>

<p>
<a href="#job_submission">1 - Job submission in general</a><br>
<a href="#finding_form">1.1 - Finding the prediction submission form</a><br>
<a href="#general_data">1.2 - Filling in general job data</a><br>
<a href="#email">1.2.1 - E-mail address</a><br>
<a href="#params">1.2.2 - AUGUSTUS species parameters</a><br>
<a href="#param_archive">1.2.2.1 - Uploading an archive file</a><br>
<a href="#param_id">1.2.2.2 - Project identifier</a><br>
<a href="#param_id">1.2.2.3 - What the AUGUSTUS species parameters are used for</a><br>
<a href="#genome_file">1.2.3 - Genome file</a><br>
<a href="#genome_file_format">1.2.3.1 - Genome file format</a></a><br>
<a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a><br>
<a href="#genome_file_purpose">1.2.3.3 - What the genome file is used for</a><br>
<a href="#optional">1.3 - Optional fields</a><br>
<a href="#cDNA">1.3.1 - cDNA file</a><br>
<a href="#cDNA_format">1.3.1.1 - cDNA file format</a><br>
<a href="#cDNA_upload">1.3.1.2 - cDNA file upload options</a><br>
<a href="#cDNA_purpose">1.3.1.3 - What cDNA files are used for</a><br>
<a href="#hints">1.3.2 - Hints file</a><br>
<a href="#hints_format">1.3.1.1 - Hints file format</a><br>
<a href="#hints_purpose">1.3.1.2 - What hints files are used for</a><br>
<a href="#utr">1.3.3 - UTR prediction</a><br>
<a href="#strand">1.3.4 - Strand specific prediction</a><br>
<a href="#alternative">1.3.5 - Alternative transcripts</a><br>
<a href="#allowed_structure">1.3.6 - Allowed gene structure</a><br>
<a href="#verification">1.4 - Verfification that you are a human</a><br>
<a href="#submitt">1.5 - The submitt button</a><br><br>
<a href="#job_status">2 - What happens after submission</a><br>
<a href="#duplication">2.1 - Submission duplication</a><br>
<a href="#error">2.2 - Errors during prediction</a><br><br>
<a href="#results">3 - Prediction Results</a><br>

</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="job_submission"><h1><a href="#job_submission">1 - Job submission in general</a></h1></div>

<p>The pipeline invoked by submitting a job to the AUGUSTUS prediction web server application is straight forward. If a cDNA file is supplied, hints are first generated from this cDNA file. If no cDNA file is supplied, AUGUSTUS is immediately called with the specified parameters.</p>

<p>The input fields of the AUGUSTUS prediction web server application form are: E-mail, AUGUSTUS species parameters, <b>Genome file</b>, <b>cDNA file</b>, <b>Hints file</b> and a number of options in form of checkboxes.</p>

<p>Please be aware that the submission of cDNA files will invoke the software <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="trainingtutorial#ref2">2</a>], which is on our server available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b>, only.


<p>In the following, you find detailed instructions for submitting an AUGUSTUS prediction job.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="finding_form"><h2><a href="#finding_form">1.1 - Finding the prediction submission form</a></h2></div>

<p>You find the AUGUSTUS prediction submission form by clicking on the following link in the left side navigation bar:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/submitt-link-pred.jpg" alt="image of submission link"></td></tr>
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
<tr><td><img src="images/email-prediction.jpg" alt="image of e-mail address field"></td></tr>
</table>
</p>


<p>The e-mail address is an obligatory field. Currently, there is no option for submitting AUGUSTUS jobs prediction through this server anonymously. If you prefer anonymous job submission, please check out our <a href="http://bioinf.uni-greifswald.de/augustus/submission">old AUGUSTUS gene prediction web server application</a>. It's functionality is slightly different, e.g. one cannot use one's own parameters, but maybe you'll find there what you are looking for.</p>

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

<div id="params"><h3><a href="#params">1.2.2 - AUGUSTUS species parameters</a></h3></div>

<p>The web server application offers you two options to specify with parameter set you want to use for predicting genes with AUGUSTUS. You can either uploaded a <tt>*.tar.gz</tt> parameter archive from your local harddrive, or you can specify the job ID of a previously finished AUGUSTUS web server application training run.</p>


<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/parameters.jpg" alt="image of parameter upload field"></td></tr>
</table>
</p>

<div id="param_archive"><h3><a href="#param_archive">1.2.2.1 - Uploading an archive file</a</h3></div>

<p>A <tt>*.tar.gz</tt> archive with a folder containing the following files is required for predicting genes in a new genome with pre-trained parameters:<p>

<p>
<ul>
<li><i>species</i>/<i>species</i>_parameters.cfg
<li><i>species</i>/<i>species</i>_metapars.cfg
<li><i>species</i>/<i>species</i>_metapars.utr.cfg
<li><i>species</i>/<i>species</i>_exon_probs.pbl.withoutCRF
<li><i>species</i>/<i>species</i>_exon_probs.pbl
<li><i>species</i>/<i>species</i>_weightmatrix.txt
<li><i>species</i>/<i>species</i>_intron_probs.pbl
<li><i>species</i>/<i>species</i>_intron_probs.pbl.withoutCRF
<li><i>species</i>/<i>species</i>_igenic_probs.pbl
<li><i>species</i>/<i>species</i>_igenic_probs.pbl.withoutCRF
</ul>
</p>
<p>where <i>species</i> is replaced by the name of the species you trained AUGUSTUS for (e.g. carrot would result it <i>carrot</i>/<i>carrot</i>_parameters.cfg). The additional <i>species</i> before the slash means that all those files must reside in a directory that is called <i>species</i> before you tar and gzip it. If you simply tar and gzip the folder that contains parameters of an AUGUSTUS training run, everything should work fine.</p>

<div id="param_id"><h3><a href="#param_id">1.2.2.2 - Project identifier</a></h3></div>

<p>If you trained AUGUSTUS on this webserver, you may instead of downloading and re-uploading a parameter archive, simply specify the project identifier of this training run. You find the project identifier for example in the job confirmation e-mail. It starts either with <i>train</i> or with <i>pred</i> and is followed by 8 digits.</p>

<p>In addition to using parameters that you trained yourself, you may also use pre-trained parameters for the following species:

<p>
<table border="1">
<font size="1">
<tr>
<td><font size="1"><b>Species</b></font></td><td><font size="1"><b>Project identifier</b></font></td><td><font size="1"><b>Courtesy of</b></font></td>
</tr>
<tr>
<td><font size="1"><b><i>Animals</i></b></font></td><td><font size="1"><b></b></font></td><td><font size="1"><b></b></font></td>
</tr>
<tr>
<td><font size="1">Acyrthosiphon pisum</font></td><td><font size="1">pea_aphid</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Aedes aegypti</font></td><td><font size="1">aedes</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Amphimedon queenslandica</font></td><td><font size="1">amphimedon</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Apis mellifera</font></td><td><font size="1">honeybee1</font></td><td><font size="1">Katharina Hoff and Mario Stanke</font></td>
</tr>
<tr>
<td><font size="1">Bombus terrestris</font></td><td><font size="1">bombus_terrestris2</font></td><td><font size="1">Katharina Hoff</font></td>
</tr>
<tr>
<td><font size="1">Brugia malayi</font></td><td><font size="1">brugia</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Caenorhabditis elegans</font></td><td><font size="1">caenorhabditis</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Callorhinchus milli</font></td><td><font size="1">elephant_shark</font></td><td><font size="1">Tereza Manousaki and Shigehiro Kuraku</font></td>
</tr>
<tr>
<td><font size="1">Drosophila melanogaster</font></td><td><font size="1">fly</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Heliconius melpomene</font></td><td><font size="1">heliconius_melpomene1</font></td><td><font size="1">Sebastian Adler and Katharina Hoff</font></td>
</tr>
<tr>
<td><font size="1">Homo sapiens</font></td><td><font size="1">human</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Nasonia vitripennis</font></td><td><font size="1">nasonia</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Petromyzon marinus</font></td><td><font size="1">lamprey</font></td><td><font size="1">Falk Hildebrand and Shigehiro Kuraku</font></td>
</tr>
<tr>
<td><font size="1">Schistosoma mansoni</font></td><td><font size="1">schistosoma</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Tribolium castaneum</font></td><td><font size="1">tribolium</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Trichinella spiralis</font></td><td><font size="1">trichinella</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1"><b><i>Alveolata</i></b></font></td><td><font size="1"><b></b></font></td><td><font size="1"><b></b></font></td>
</tr>
<tr>
<td><font size="1">Tetrahymena thermophila</font></td><td><font size="1">tetrahymena</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Toxoplasma gondii</font></td><td><font size="1">toxoplasma</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1"><b><i>Plants and algae</i></b></font></td><td><font size="1"><b></b></font></td><td><font size="1"><b></b></font></td>
</tr>
<tr>
<td><font size="1">Arabidopsis thaliana</font></td><td><font size="1">arabidopsis</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Chlamydomonas reinhartii</font></td><td><font size="1">chlamy2011</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Galdieria sulphuraria</font></td><td><font size="1">galdieria</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Zea mays</font></td><td><font size="1">maize</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1">Solanum lycopersicum</font></td><td><font size="1">tomato</font></td><td><font size="1"></font></td>
</tr>
<tr>
<td><font size="1"><b><i>Fungi</i></b></font></td><td><font size="1"><b></b></font></td><td><font size="1"><b></b></font></td>
</tr>
<tr>
<td><font size="1">Aspergillus fumigatus</font></td><td><font size="1">aspergillus_fumigatus</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Aspergillus nidulans</font></td><td><font size="1">aspergillus_nidulans</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Aspergillus oryzae</font></td><td><font size="1">aspergillus_oryzae</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Aspergillus terreus</font></td><td><font size="1">aspergillus_terreus</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Botrytis cinerea</font></td><td><font size="1">botrytis_cinerea</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Candida albicans</font></td><td><font size="1">candida_albicans</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Candida guilliermondii</font></td><td><font size="1">candida_guilliermondii</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Candida tropicalis</font></td><td><font size="1">candida_tropicalis</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Chaetomium globosum</font></td><td><font size="1">chaetomium_globosum</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Coccidioides immitis</font></td><td><font size="1">coccidioides_immitis</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Coprinus cinereus</font></td><td><font size="1">coprinus</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Cryptococcus neoformans</font></td><td><font size="1">cryptococcus_neoformans_neoformans_B</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Debaryomyces hansenii</font></td><td><font size="1">debaryomyces_hansenii</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Encephalitozoon cuniculi</font></td><td><font size="1">encephalitozoon_cuniculi_GB</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Eremothecium gossypii</font></td><td><font size="1">eremothecium_gossypii</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Fusarium graminearum</font></td><td><font size="1">fusarium_graminearum</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Histoplasma capsulatum</font></td><td><font size="1">histoplasma_capsulatum</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Kluyveromyces lactis</font></td><td><font size="1">kluyveromyces_lactis</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Laccaria bicolor</font></td><td><font size="1">laccaria_bicolor</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Lodderomyces elongisporus</font></td><td><font size="1">lodderomyces_elongisporus</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Magnaporthe grisea</font></td><td><font size="1">magnaporthe_grisea</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Neurospora crassa</font></td><td><font size="1">neurospora_crassa</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Phanerochaete chrysosporium</font></td><td><font size="1">phanerochaete_chrysosporium</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Pichia stipitis</font></td><td><font size="1">pichia_stipitis</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Phizopus oryzae</font></td><td><font size="1">rhizopus_oryzae</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Saccharomyces cerevisiae</font></td><td><font size="1">saccharomyces_cerevisiae_S288C</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Schizosaccharomyces pombe</font></td><td><font size="1">schizosaccharomyces_pombe</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Ustilago maydis</font></td><td><font size="1">ustilago_maydis</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
<tr>
<td><font size="1">Verticillium longisporum</font></td><td><font size="1">verticillium_longisporum1</font></td><td><font size="1">Katharina Hoff and Mario Stanke</font></td>
</tr>
<tr>
<td><font size="1">Yarrowia lipolytica</font></td><td><font size="1">yarrowia_lipolytica</font></td><td><font size="1">Jason Stajich</font></td>
</tr>
</font>
</table>
</p>

<p>
Please let us know whether you want to have parameters that you trained for a certain species to be included in this public list! If they are included in this list, they will also be distributed with the upcoming AUGUSTUS release.</p>
</p>


<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="genome_file"><h3><a href="#genome_file">1.2.3 - Genome file</a></h3></div>

<p>The genome file is an obligatory file for predicting genes with AUGUSTUS.</p>

<div id="genome_file_format"><h4><a href="#genome_file_format">1.2.3.1 - Genome file format</a></h4></div>

The genome file must contain the genome sequence in (multiple) fasta format. Every <i>unique</i> header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver application. White spaces in the sequence header might cause problems if the first word after the leading character <b>></b> is identical for several fasta entries. We generally recommend short, unique, non-white-space containing fasta headers.<br><br>
            <b>Correct file format example:</b>


<table border="2" cellspacing="0" cellpadding="0">
<tr><td>            <pre class="example"><font size="1">
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
            </font></pre></td></tr>
</table>


</p>

<p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for genome file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>.</p>

<div id="genome_file_upload"><h4><a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a></h4></div>

<p>The AUGUSTUS prediction web server application offers two possiblities for transferring the genome file to the server: <i>Upload a file</i> and <i>specify a web link to file</i>.</p>

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

<p>The genome file is used as a template for gene prediction, it is the sequence in which you want to predict genes.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="optional"><h2><a href="#optional">1.3 - Optional fields</a></h2></div>

<p>This section describes a number of fields that are optional for predicting genes with AUGUSTUS.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="cDNA"><h3><a href="#cDNA">1.3.1 - cDNA file</a></h3></div>

<p>This feature is available only <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> as this is required by the BLAT license.</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/training-cDNA.jpg" alt="image of cDNA submission field"></td></tr>
</table>
</p>


<div id="cDNA_format"><h3><a href="#cDNA_format">1.3.1.1 - cDNA file format</a></h3></div>

<p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. An example for correct cDNA file format is given at <a href="#genome_file_format">1.2.3.1 - Genome file format</a>.</p>

<p>It is currently possible to submitt assembled RNA-seq transcripts instead of or mixed with ESTs as a cDNA/EST file. However, you should be aware that RNA-seq files are often much bigger than EST or cDNA files, which increases runtime of a prediction job. In order to keep runtime of your prediction job as low as possible, you should remove all assembled RNA-seq transcripts from your file that do not map to the submitted genome sequence. (In principle, this holds true for EST and cDNA files, too, but there, the problem is not as pronounced due to a smaller number of sequences.)</p>

<p>It is currently not allowed to upload RNA-seq raw sequences. (We filter for the average length of cDNA fasta entries and may reject the entire training job in case the sequences are on average too short, i.e. shorter than 400 bp.)</p>

<p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for cDNA file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>. The maximal supported file size is 1 GB.</p>

<div id="cDNA_upload"><h3><a href="#cDNA_upload">1.3.1.2 - cDNA file upload options</a></h3></div>

<p>There are two options for cDNA file upload: upload from your local harddrive, or upload from a public http or ftp server. Please see <a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a> for a more detailed description of upload options.</p>

<div id="cDNA_purpose"><h3><a href="#cDNA_purpose">1.3.1.3 - What cDNA files are used for</a></h3></div>

<p>The cDNA file is used for generating extrinsic evidence for gene structures in the gene prediction process, also called <i>hints</i></li>
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


<div id="hints"><h3><a href="#hints">1.3.2 - Hints file</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/hints-field.jpg" alt="image of hints file upload field"></td></tr>
</table>
</p>

<div id="structure_file_format"><h4><a href="#structure_file_format">1.3.3.1 - Hints file format</a></h4></div>

	    <p>It is possible to submit an externally created file that contains extrinsic evidence for gene structures in gff format.</p>
	    <p>In general, gff files must contain the following columns (the columns are separated by <b>tabulators</b>):</p>

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
            <li>For usage as hint, <b>Attribute</b> must contain the string <b>source=M</b> (for manual). Other sources, such EST or protein, are possible, but only in the command line version of AUGUSTUS. Source types other than <b>M</b> are ignored by AUGUSTUS web server applications.
            </OL>
</p>
	    <p>
            <br>
            <b>Correct format example:</b>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><font size="1"><pre class="example">
HS04636 anchor  exonpart        500     506     0       -       .       source=M
HS04636 anchor  exon            966     1017    0       +       0       source=M
HS04636 anchor  start           966     968     0       +       0       source=M
HS04636 anchor  dss             2199    2199    0       +       .       source=M
HS04636 anchor  stop            7631    7633    0       +       0       source=M
HS04636 anchor  intronpart      7631    7633    0       +       0       source=M
            </pre></font></td></tr>
</table>   
            </p>

<div id="hints_purpose"><h4><a href="#hints_purpose">1.3.1.2 - What hints files are used for</a></h4></div>

<p>The hints file is used as extrinsic evidence that supports gene structure prediction. You can generate hints yourself based on any alignment program and information resource (e.g. ESTs, RNA-seq data, peptides, proteins, ...) that appears suitable to you.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="utr"><h3><a href="#utr">1.3.3 - UTR prediction</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/utr-checkbox.jpg" alt="image of utr checkbox"></td></tr>
</table>
</p>

<p>It takes significantly more time to predict UTRs but in addition to reporting UTRs, it usually is also a little more accurate on the coding regions when ESTs are given as extrinsic evidence.</p>

<p>UTR prediction is only possible if UTR parameter files exist for your species. Even if UTR parameter files exist for a species, you should make sure, that they are <i>species specific</i>, i.e. have actually been optimized for your target species. It is a waste of time to predict UTRs with <i>general</i> (template) parameters.</p>

<p>If no UTR parameter files exist for your species but you enables UTR prediction in the form, the web server application will overrule the choice to predict UTRs by simply not predicting any UTRs.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="strand"><h3><a href="#strand">1.3.4 - Strand specific prediction</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/strands.jpg" alt="image of the strand checkboxes"></td></tr>
</table>
</p>

<p>By default, AUGUSTUS predicts genes in both strands but you may alter this behavior by checking another radio button in this field to predict genes in the forward (+) or reverse (-) strand, only.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="alternative"><h3><a href="#alternative">1.3.5 - Alternative transcripts</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/alternative.jpg" alt="image of alternative transcript radio buttons"></td></tr>
</table>
</p>


<p>By default, AUGUSTUS does not predict any alternative transcripts.</p> <p>If you select <b>few</b>, then the following AUGUSTUS parameters are set to result in the prediction of relatively few alternative transcripts:<br><tt>--alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 <br>--maxtracks=2</tt>.</p> <p> If you select <b>medium</b> the AUGUSTUS parameters are set to <br><tt>--alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.4 <br>--maxtracks=3"</tt>.</p> <p> If you select <b>many</b>, AUGUSTUS parameters are set to <br><tt>--alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 <br>--maxtracks=20"</tt>. </p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="allowed_structure"><h3><a href="#allowed_structure">1.3.6 - Allowed gene structure</a></h3></div>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/allowed-structures.jpg" alt="image of allowed gene structure buttons"></td></tr>
</table>
</p>


<p><b>Predict any number of (possibly partial) genes:</b> This option is set by default. AUGUSTUS may predict no gene at all, one or more genes. The genes at the boundaries of the input sequence may be partial. Partial here means that not all of the exons of a gene are contained in the input sequence, but it is assumed that the sequence starts or ends in a non-coding region.<p>

<p><b>Predict only complete genes:</b> AUGUSTUS assumes that the input sequence does not start or end within a gene. Zero or more complete genes are predicted.</p>

<p><b>Predict only complete genes - at least one:</b> As the previous option. But AUGUSTUS predicts at least one gene (if possible).</p>

<p><b>Predict exactly one complete gene:</b> AUGUSTUS assumes that the sequence contains exactly one complete gene. Note: This feature does not work properly in combination with alternative transcripts. </p>

<p><b>Ignore conflicts with other strand:</b> By default AUGUSTUS assumes that no genes - even on opposite strands - overlap. Indeed, this usually is the case but sometimes an intron contains a gene on the opposite strand. In this case, or when AUGUSTUS makes a false prediction on the one strand because it falsely thinks there is a conflicting gene on the other strand, AUGUSTUS should be run with this option set. It then predicts the genes on each strand separately and independently. This may lead to more false positive predictions, though.</p>

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
<tr><td><img src="images/pred-submitt.jpg" alt="image of submission button"></td></tr>
</table>
</p>


<p>After filling out the appropriate fields in the submission form, you have to click on the button that says "Start Predicting" at the bottom of the page. It might take a while until you are redirected to the status page of your job. The reason is that we are checking various file formats prior job acceptance, and that the transfer of files from your local harddrive to our server might take a while. Please be patient and wait until you are redirected to the status page! Do not click the button more than once (it won't do any harm but it also doesn't speed up anything).</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="job_status"><h1><a href="#job_status">2 - What happens after submission</a></h1></div>

<p>After you click the "Start Predicting" button, the web server application first validates whether the combination of your input fields is generally correct. If you did anything wrong, you will be redirected to the training submission form and an error message will be displayed at the top of the page.</p>

<p>If all fields were filled in correctly, the application is actually initiated. You will receive an e-mail that confirms your job submission and that contains a link to the job status page.</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/prediction-job-status.jpg" alt="image of job status page"></td></tr>
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

<p>Since predicting genes wiht AUGUSTUS may under certain circumstances be is a very resource consuming process, we try to avoid data duplication. In case you or somebody else tries to submitt exactly the same input file combination more than once, the duplicated job will be deleted and the submitter of the redundant job will receive an e-mail with the link to the previously submitted job status and with a link to the results page (which may still be empty in case the duplicated job has not finished computing, yet.)</p>

<p>Your web browser will be redirected to a page with the following content if data duplication occured:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/job-does-not-exist-prediction.jpg" alt="image of non existing job status"></td></tr>
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

<div id="error"><h2><a href="#error">2.2 - Errors during prediction</a></h2></div>

<p>You will receive an e-mail in case an error occurs during the AUGUSTUS gene prediction process. The admin of this server is also notified by e-mail about errors. We will get in touch with you, again, after we figured out what caused the error.</p>

<p>Since the web server application is currently in beta testing phase, unexpected errors might still occur. Therefore we ask for you help on reporting any unexpected errors.</p>

<h3>Reasons that should not lead to the assumption of error occurence:</h3>

<p><ul><li>Your job waiting for execution shorter than two months</li>
<li>Your job is computing for shorter than two weeks</li>
</ul></p>

<h3>Reasons to report errors:</h3>

<p><ul><li>Your job waiting for execution longer than two months</li>
<li>Your job is computing longer than two weeks</li>
<li>The job status page shows "finished" but you did not receive an e-mail with the link to results</li>
<li>A Grails execption error is displayed</li>
<li>You results page is empty although you received a confirmation e-mail that says your job finished</li>
<li>...</li>
</ul></p>

<p>Please report unexpected errors to augustus-web@uni-greifswald.de. Please include which actions from your side exactly caused the error, and also copy the Grails exception message into your e-mail in case such a message was displayed.</p>

<p><a href="#seitenanfang">
<img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
Top of page
</a>
<br>
</p>
<hr>
<br>

<div id="results"><h1><a href="#results">3 - Prediction Results</a></h1></div>

<p>After job computations have finished, you will receive an e-mail with a link to the results of your submission. The linked web page may look similar to this:</p>

<p>
<table border="2" cellspacing="0" cellpadding="0">
<tr><td><img src="images/prediction-results-example.jpg" alt="image of results example"></td></tr>
</table>
</p>

<p>This page should contain the file <b>augustus.tar.gz</b>. Please make a "right click" on the link and select "Save As" (or similar) to save the file on your local harddrive.</p>

<p><b>augustus.tar.gz</b> is a gene prediction archive and its content depends on the input file combination. You can unpack the archive by typing <tt>tar -xzvf *.tar.gz</tt> into your shell. (You find more information about the software tar at the <a href="http://www.gnu.org/s/tar/">GNU tar website</a>.)
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
