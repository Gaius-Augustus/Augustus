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
    Help
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
         <li><a href="trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
         <li><g:link controller="training" action="create"><span>Submit Training</span></g:link></li>
         <li><a href="predictiontutorial.gsp"><span>Prediction Tutorial</span></a></li>
         <li><g:link controller="prediction" action="create"><span>Submit Prediction</span></g:link></li>
         <li id="current"><a href="help.gsp"><span>Help</span></a></li>
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
    <a name="inhalt" id="inhalt"></a>
    <table class="contentpaneopen">
      <tr>
	<td class="contentheading" width="100%">
	  <a href="accuracy.gsp" class="contentpagetitle">Help</a>
        </td>
      </tr>
    </table>


            <p>This website contains <i>short instructions</i> and some <i>frequently asked questions</i> concerning <ul><li>the training of AUGUSTUS and <li>predicting genes in a new genomic sequence with pre-trained parameters.</ul> <br>For more detailed instructions, please read <a href="trainingtutorial.gsp">Training Tutorial</a> and <a href="predictiontutorial.gsp">Prediction Tutorial</a>.</p>
            <hr>
            <br>
	    <h2>Contents</h2>
	    <p>
	      <a href="#noResults">Why do I not get any results?</a><br>
	      <a href="#buisy">Why is the server busy?</a><br>
              <a href="#species_name">What is the species name?</a><br>
	      <a href="#email">Why should I give my e-mail address?</a><br>
              <a href="#upload_link">File upload versus web link</a><br>
	      <a href="#most_common_problem">Instructions for fasta headers</a><br>
              <a href="#which_files">Which files must or can I submit for training AUGUSTUS?</a><br>
              <a href="#which_files_pred">Which files are required for predicting genes in a new genome?</a><br>
              <a href="#genome_file">Genome file</a><br>
              <a href="#cDNA">cDNA file</a><br>
              <a href="#protein">Protein file</a><br>
              <a href="#structure">Training gene structure file</a><br>
              <a href="#hints">Hints file</a><br>
              <a href="#archive">Parameter archive</a><br>
              <a href="#project_id">What is the project identifier?</a><br>
              <a href="#job_status">What does my job status mean?</a><br>
              <a href="#utr">UTR prediction: yes or no?</a><br>
              <a href="#allowedGeneStructure">Allowed gene structure<a><br>
              <a href="#list">What about that data duplication?</a><br>
              <a href="#accuracy">Why is the prediction accuracy in the genome of my species not as good as I expected?</a><br>
              <a href="#privacy">What about data privacy and security?</a><br>
	      <a href="#results_pred">Gene prediction results</a><br>
	      <a href="#results_train">Training results</a><br>
              <a href="#commercial">I am not from academia/not non-profit. What can I do?</a><br>
              <a href="#dog">Why do I see a running dog after pressing the submission button?</a><br>
            </p>
            <hr>
	    <br>
	    <h2 id="noResults">Why do I not get any results?</h2>
		<p>
		<ul>
			<li><b>Did an obvious error occur?</b><br>Please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a> if you are not sure what the error message is telling you.<br>
			  <p>One frequently occuring error in the AutoAug.err file is the following:</p>

<p><b>The file with UTR parameters for train****** does not seem to exist.</b> This likely means that the UTR model has not beeen trained yet for train******.</p>

			  <p>This error message tells you that no UTR parameters were trained for your species. If no other error messages are contained above the first UTR error message, the general results of your job are ok, you simply did not get UTR parameters and thus no predictions with UTR.</p>

                          <p><b>Illegal division by zero at /usr/local/augustus/trunks/scripts/autoAugTrain.pl line 241.<br>
failed to execute: No such file or directory<br></b> This error occurs when not training gene structures were generated/available. This may be caused by one of the following circumstances:<br>

<ul>
<li>You supplied a genome and protein file. In this case, Scipio was not able to generate any <i>complete</i> gene structures from the data set. In most cases, some <i>incomplete</i> gene structures were produced, but since they frequently cause crashes in the augustus training routine, we do not use them within the web service.</li>

<li>The files that you supplied had long an complex fasta headers. This causes problems with PASA and Scipio. Take care that the fasta headers in all your files are unique, short, do not contain whitespaces or special characters.</li>
</ul>

</li>
			<li><b>Did you submitt your job a long time ago and it seems to be "stuck" at the status of "computing"?</b><br>Please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a> to inquire whether your job is really still running.</li>
			<li><b>Did your job finish but there are just no parameters or predictions?</b><br>The quality of results depends on the quality and combination of your input data. If the input data did e.g. not provide sufficient information for generating training genes, then no AUGUSTUS parameters will be optimized for your species, and no predictions will be made. In case of the gene prediction web server application, it is also possible that your submitted genome sequence does not contain any protein coding genes.</li>
		</ul>
		</p>
	    <p></p>
<p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <h2 id="buisy">Why is the server busy?</h2>
            <p>Training AUGUSTUS is a very resource and time consuming process. We use a grid engine queuing system with a limited number of waiting slots. If we estimate that the time from job submission to computation start might be very long, our web server might display a message that our server is buisy. The submission of new jobs is then disabled (prediction and training submission will both be disabled). Please wait one or two weeks before you try a new submission. If the problem persists longer than a month, or if your job is urgent, please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a>.</p>
	    <p><a href="#seitenanfang">
	      <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
	      Top of page
	    </a>
	    <br>
	    </p>
            <hr>
	    <br>
            <div id="species_name"><h2>What is the species name?</h2></div>
<p>The species name is the name of the species for whose genome you want to train AUGUSTUS. The species name is an obligatory parameter. Considering that AUGUSTUS training is such a time consuming process, our objective is to know the names of species for which AUGUSTUS was trained in order to make the trained parameters available to the public so that others who are interested in the same species as you do not have to rerun the training process. (We will only explicitely publish your parameter set with the next AUGUSTUS release after confirming via e-mail that you agree to this.)</p>

<p>However, if you do not want to reveal the true species name, you may use any other string shorter than 30 characters as a species name.</p>

<p>The species name is not allowed to contain spaces!</p>


            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="email"><h2>Why should I give my e-mail address?</h2></div>
	    <p>Unlike many other bioinformatics web services, the AUGUSTUS web server application is not an implementation of a fail-safe procedure. Particularly the assembly of a training gene set from extrinsic data (ESTs and protein sequences) and a genome sequence may not always work perfectly. Our pipeline may issue warnings or errors, and sometimes, we need to get some feedback from you via e-mail in order to figure out what is the problem with your particular input data set.</p>

            <p>In addition, training and running AUGUSTUS are rather time consuming processes that may take up to several weeks (depending on the input data). It may be more convenient to receive an e-mail notification about your job having finished, than checking the status page over and over, again.</p>
	    <p> Therefore, we strongly recommend that you enter an e-mail adress.</p>
            <p> If supplied, we use your e-mail address for the following purposes:</p>

<p>
<ul>
<li>Confirming your job submission</li>
<li>Confirming successful file upload (for large files via ftp/http)</li>
<li>Notifying you about your job having finished</li>
<li>Informing you about any problems that might occur during your particular job and asking questions about that job in order to solve those problems</li>
<li>Asking you whether we should include your species parameters into the next AUGUSTUS release (applies to training, only)</li>
</ul>
</p>

<p>We do <b>not</b> use your e-mail address to send you any <i>spam</i>, i.e. about web service updates. We do <b>not</b> share your e-mail address with any third parties.</p>

<p>Job submission without giving an email adress is possible but discouraged.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="upload_link"><h2>File upload versus web link</h2></div>
            <p> The AUGUSTUS training and prediction web server application offers in some cases two possiblities for transferring files to the server: <i>Upload a file</i> and <i>specify a web link to file</i>.</p>
             <p>
             <ul>
                <li>For <b>small files</b>, please click on the Browse-button and select a file on your harddrive.<br>If you experience a <i>Connection timeout</i> (because your file was too large for this type of upload), please use the option for large files!</li>
                <li><b>Large files</b> can be retrieved from a <b>public</b> web link. Specify a valid ftp or http URL to your sequence file. Our server will fetch the file from the given address.</li>
            </ul>
	     </p>
	    <p>
            <b>You cannot do both at the same time!</b> For each file type (e.g. the genome file), you must <b>either</b> select a file on your harddrive <b>or</b> give a web link!</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
	    <div id="most_common_problem"><h2>Instructions for fasta headers</h2></div>
            <p>We observed that most problems with generating training genes for training AUGUSTUS are caused by fasta headers in the sequence files. Some of the tools in our pipeline will truncate fasta headers if they are too long or contain spaces, or contain special characters. This definitely leads to a lot of warning messages in the AutoAug.err file, and it may also lead to non-unique fasta entry names, which will lead to a crash of the pipeline. We therefore strongly recommend that you adhere to the following rules for fasta headers when using our web services:
	      <ul>
		<li>no whitespaces in the headers</li>
		<li>no special characters in the headers (e.g. <tt>!#@&|;</tt>)</li>
		<li>make the headers as short as possible</li>
		<li>let headers not start with a number but with a letter</li>
		<li>let headers contain letters and numbers, only</li>
	      </ul>
	    </p>
	    <p>In the following we give some header examples that will not cause problems:<br><br>
	      <tt>>entry1</tt><br>
	      <tt>>contig1000</tt><br>
	      <tt>>est20</tt><br>
	      <tt>>scaffold239</tt>
	      </p>
	    <p>The following kinds of headers will cause at least warning messages but probably also a pipeline crash:<br><br>
	      <tt>>contig1 length=1000 Arabidopsis thaliana</tt><br>
	      <tt>>gi|123344545|some_protein|some_species</tt><br>
	      <tt>>Drosophila melanogaster scaffold 10000</tt><br>
	      </p>
	    <p>If you have a fasta file with unsuitable headers and you do not know how to modify them automatically, you may use the Perl script <a href="http://bioinf.uni-greifswald.de/bioinf/downloads/simplifyFastaHeaders.pl">simplifyFastaHeaders.pl</a>. After saving it on your local Unix system, first check whether the location of Perl in the first line of the script is correct for your system (<tt>#!/usr/bin/perl</tt>). If Perl is installed in another location, you need to modify that line! Then, execute the script with the following parameters:<br><br>
<tt>perl simplifyFastaHeaders.pl in.fa nameStem out.fa header.map</tt><br><br>
<ul>
<li><tt>in.fa</tt> is the input fasta file, it must already be in valid fasta format</li>
<li><tt>nameStem</tt> is a character descriptor that will be used as a start for all simplified headers, e.g. <tt>est</tt>, or <tt>contig</tt>, or <tt>protein</tt>, etc. Be aware that fasta headers must always be unique, so choose different nameStems for genome and cDNA and protein file!</li>
<li><tt>out.fa</tt></li> is the output fasta file with simplified fasta headers, this file can be processed by our web service.</li>
<li><tt>header.map</tt></li> is a map that contains the simplied header and the original header in a tabular separated format.</li>
</ul>
</p>
<p><b>Why is the simplification of fasta headers not a built in function of the web service?</b> The reason is that we think you should be able to recognize the predictions later on! Gene predictions will be made available in gff format, which contains the sequence name in the first column. Therefore, you should modify the fasta headers yourself, before submitting data to the web service!</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
	    <div id="which_files"><h2>Which files must or can I submit for training AUGUSTUS?</h2></div>
            <p>You need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li><b>at least one out of the following files:</b>  <a  href="#cDNA">cDNA file</a>, <a href="#structure">training gene structure file</a>, and <a href="#protein">protein file</a>.</li>
            </ul>
	    </p>
	    <p>
            Please consider that training AUGUSTUS is a time and resource consuming process. For optimal results, you should specify as much information as possible for a single training run instead of starting the AUGUSTUS training multiple times with different file combinations! If you have a lot of EST data, we recommend that you submitt ESTs instead of protein sequences since ESTs will likely allow the generation of a UTR training set.</p>
	                <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="which_files_pred"><h2>Which files are required for predicting genes in a new genome?</h2></div>
            <p>For predicting genes in a new genome with already trained parameters, you need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li>a <a href="#archive">parameter archive</a>. Instead of uploading the archive, you may also enter a valid <a href="#project_id">project identifier</a> in case you trained AUGUSTUS on this web server and the training has already finished; or you may select pre-trained parameter set from the drop down menu.
            </ul></p>
	    <p>
            You may in addition specify an <a href="#cDNA">EST/cDNA file</a> and/or a <a href="#structure">hints file</a> that will be used as extrinsic evidence for predicting genes.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="genome_file"><h2>Genome file</h2></div>
            <p>The genome file is an obligatory file for training AUGUSTUS and for making predictions with pre-trained parameters in a new genome. It must contain the genome sequence in (multiple) fasta format. Every header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver applications.</p>
<p>Headers must be <b>unique</b> within a file! We recommend that you use <b>short fasta headers</b>. Headers like<br>
<pre>
>gi|382483733|gb|GZ667513.1|GW667513 SSH_BP_47 Some species
Wicked root cDNA library Some species cDNA clone SSH_BP_47 
similar to Putative NADH-cytochrome B5 reductase, mRNA sequence
</pre><br>
are likely to cause a lot of warning messages. An example for a short header created from the too long header above:<br><br>
<pre>
>GZ667513.1
</pre>
<br><br>
            <b>Correct file format example:</b>
            <pre class="example">
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
            </pre>
	    </p>
	                <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="cDNA"><h2>cDNA file</h2></div>
            <p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. See <a href="#genome_file">Genome file</a> for a format example. Upload of a cDNA file to our web server application will invoke the software <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="trainingtutorial.gsp#ref2">2</a>], which is on our webserver application only available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b>.</p>
	    <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="protein"><h2>Protein file</h2></div>
            <p>The protein file is a multiple fasta file that contains protein sequences as supporting evidence for genes. Allowed sequence characters: <b>A a R r N n D d C c E e Q q G g H h I i L l K k M m F f P p S s T t W w Y y V v B b Z z J j X x</b>. Empty lines are not allowed but will simply be removed from the file by the webserver application. <br><br></p>
            <p><b>Correct file format example:</b>
            <pre class="example">
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
            </pre></p>
	   <p>Submitting a protein file to our AUGUSTUS training web server application will invoke <a href="http://www.webscipio.org/">Scipio</a> [<a href="trainingtutorial.gsp#ref3">3</a>], which uses <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="trainingtutorial.gsp#ref2">2</a>]. Therefore, protein file upload is only available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server application.</p>
           <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="structure"><h2>Training gene structure file</h2></div>
	    <p>You can submit your own, externally created training gene structures to the AUGUSTUS training web server application. Regardless of the format, gene structure files are not allowed to contain java metacharacters like "*" or "?".</p>

<p>Training gene structure files can be submitted in two different formats: Genbank format or gff format.</p>

<h3>Training gene structure file in genbank format</h3>

<p>Gene structures in genbank format must contain the coding sequence parts and flanking regions. Flanking regions are important because AUGUSTUS is supposed to differentiate between genes and intergenic regions. The length of flanking regions depends on the length of genes in the target genome. In our pipeline, flanking regions are set to the average gene length (exceptionally applying the extreme limits between 1000 and 10000 nt). It is very important to make sure that the flanking regions do not contain any other protein coding gene parts, i.e. we recommend to trim flanking regions in a way that will exclude other CDS parts.</p>

<p><b>Correct file format example (condensed view, the three dots represent further lines of sequences):</b>
            <pre class="example">
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
            </pre></p>


<p>If you want to train UTRs, you have to additionally incorporate mRNA information in your genbank file.</p>

<p><b>Correct file format example (including UTR training):</b>
            <pre class="example">
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
            </pre></p>

<h3>Training gene structure file in gff format</h3>

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
            <pre class="example">
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
</pre></p>

<p><b>Correct file format example (with UTR):</b>
            <pre class="example">
Chr.1	mySource	5'-UTR	277153	277220	45	+	.	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	277221	277238	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278100	278213	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278977	279169	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279630	279648	0.94	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279734	279768	0.94	+	1	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	280307	280344	1	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	3'-UTR	280345	280405	78	+	.	transcript_id "g22472.t1"; gene_id "g22472";
</pre></p>


           <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="hints"><h2>Hints file</h2></div>

	    <p>For the gene prediction web server application, it is possible to submit an externally created file that contains extrinsic evidence for gene structures in gff format.</p>
	    <p>In general, gff format must contain the following columns (The columns are separated by <b>tabulators</b>):</p>

<p>
            <OL TYPE="1">
            <li>The <b>sequence names</b> must be found in the fasta headers of sequences in the genome file. 
	    <li>The <b>source</b> tells with which software/process the gene structure was generated (you can fill in whatever you like).
            <li>The <b>feature</b> may for AUGUSTUS gene prediction be 
		<ul>
		<li><tt>start</tt> - translation start, specifies an interval that contains the start codon. The interval can be larger than 3 nucleotides, in which case every ATG in the interval gets a bonus.</li>
		<li><tt>stop</tt> - translation end  (stop codon)</li>
		<li><tt>tss</tt> - transcription start site</li>
		<li><tt>tts</tt> - transcription termination site</li>
		<li><tt>ass</tt> - acceptor (3') splice site, the last intron position</li>
		<li><tt>dss</tt> - donor (5') splice site, the first intron position</li>
		<li><tt>exonpart</tt> - part of an exon in the biological sense.</li>
		<li><tt>exon</tt> - complete exon in the biological sense.</li>
		<li><tt>intronpart</tt> - introns both between coding and non-coding exons.</li>
		<li><tt>intron</tt> - complete intron in the biological sense</li>
		<li><tt>CDSpart</tt> - part of the coding part of an exon. (CDS = coding sequence)</li>
		<li><tt>CDS</tt> - coding part of an exon with exact boundaries. For internal exons of a multi exon gene this is identical to the biological boundaries of the exon. For the first and the last coding exon the boundaries are the boundaries of the coding sequence (start, stop).</li>
		<li><tt>UTRpart</tt> - The hint interval must be included in the UTR part of an exon.</li>
		<li><tt>UTR</tt> - exact boundaries of a UTR exon or the untranslated part of a partially coding exon.</li>
		<li><tt>irpart</tt> - intergenic region part. The bonus applies to every base of the intergenic region. If UTR prediction is turned on (--UTR=on) then UTR is considered genic.
		<li><tt>nonexonpart</tt> -  intergenic region or intron.</li>
		<li><tt>genicpart</tt> - everything that is not intergenic region, i.e. intron or exon or UTR if applicable.</li>
		</ul>
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
            <pre class="example">
HS04636 anchor  exonpart        500     506     0       -       .       source=M
HS04636 anchor  exon            966     1017    0       +       0       source=M
HS04636 anchor  start           966     968     0       +       0       source=M
HS04636 anchor  dss             2199    2199    0       +       .       source=M
HS04636 anchor  stop            7631    7633    0       +       0       source=M
HS04636 anchor  intronpart      7631    7633    0       +       0       source=M
            </pre>
            </p>
	    <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <div id="archive"><h2>Parameter archive</h2></div>
	    <p>
            A *.tar.gz archive with a folder containing the following files is required for predicting genes in a new genome with pre-trained parameters:</p>
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
<p>where <i>species</i> is replaced by the name of the species you trained AUGUSTUS for (e.g. carrot would result it <i>carrot</i>/<i>carrot</i>_parameters.cfg). The additional <i>species</i> before the slash means that all those files must reside in a directory that is called <i>species</i>  (or in our example: <i>carrot</i>) before you tar and gzip it. If you simply tar and gzip the folder that contains parameters of an AUGUSTUS training run, everything should work fine.</p>
                       <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="project_id"><h2>What is the project identifier?</h2></div>
            <p>If you trained AUGUSTUS on this webserver, you may instead of uploading a parameter archive, simply specify the project identifier of this training run. You find the project identifier for example in the subject line for your training confirmation e-mail, where it says <i>Your AUGUSTUS training job project_id</i>. Project identitfiers typically consist of the letters <i>pred</i> or <i>train</i>, followed by a random string of 8 digits resulting in for example <i>train345kljD4</i>.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="job_status"><h2>What does my job status mean?</h2></div>
		<p>
In the beginning, the status page will display that your job has been submitted. This means, the web server application is currently uploading your files and validating file formats. After a while, the status will change to waiting for execution. This means that all file formats have been confirmed and an AUGUSTUS training job has been submitted to our grid engine, but the job is still pending. Depending on waiting queue length, this status may persist for a while. Please contact us in case you job is pending for more than one month. Later, the job status will change to computing. This means the job is currently computing. When the page displays finished, all computations have been finished and a website with your job's results has been generated.</p>

<p>You will receive an e-mail with the link to the results of your job when computations are finished if you specified an email adress.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="utr"><h2>UTR prediction: yes or no?</h2></div>
		<p>It takes significantly more time to predict UTRs but in addition to reporting UTRs, it usually is also a little more accurate on the coding regions when ESTs are given as extrinsic evidence.</p>

<p>UTR prediction is only possible if UTR parameter files exist for your species. Even if UTR parameter files exist for a species, you should make sure, that they are species specific, i.e. have actually been optimized for your target species. It is a waste of time to predict UTRs with general (template) parameters.</p>

<p>UTR prediction is only supported in combination with the following two gene structure constraints:</p>
<ul>
<li>predict any number of (possibly partial) genes</li>
<li>only predict complete genes</li>
</ul>
<p>UTR prediction is not possible in combination with the gene structure constraints:</p>
<ul>
<li>only predict complete genes - at least one</li>
<li>predict exactly one gene</li>
</ul>

<p>If no UTR parameter files exist for your species but you enables UTR prediction in the form, the web server application will overrule the choice to predict UTRs by simply not predicting any UTRs.</p>

<h3>Species for which UTR parameters are available:</h3>

<ul>
  <li>Acyrthosiphon pisum (pea_aphid)</li>
  <li>Amphimedon queenslandica (amphimedon)</li>
  <li>Apis mellifera (honeybee1)</li>
  <li>Bombus terrestris (bombus_terrestris2)</li>
  <li>Caenorhabditis elegans (caenorhabditis)</li>
  <li>Drosophila melanogaster (fly)</li>
  <li>Homo sapiens (human)</li>
  <li>Trichinella spiralis (trichinella)</li>
  <li>Toxoplasma gondii (toxoplasma)</li>
  <li>Arabidopsis thaliana (arabidopsis)</li>
  <li>Chlamydomonas reinhartii (chlamy2011)</li>
  <li>Galdieria sulphuraria (galdieria)</li>
  <li>Solanum lycopersicum (tomato)</li>
</ul>

 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="allowedGeneStructure"><h2>Allowed gene structure</h2></div>

<p><b>Predict any number of (possibly partial) genes:</b> This option is set by default. AUGUSTUS may predict no gene at all, one or more genes. The genes at the boundaries of the input sequence may be partial. Partial here means that not all of the exons of a gene are contained in the input sequence, but it is assumed that the sequence starts or ends in a non-coding region.</p>

<p><b>Predict only complete genes:</b> AUGUSTUS assumes that the input sequence does not start or end within a gene. Zero or more complete genes are predicted.</p>

<p><b>Predict only complete genes - at least one:</b> As the previous option. But AUGUSTUS predicts at least one gene (if possible).</p>

<p><b>Predict exactly one complete gene:</b> AUGUSTUS assumes that the sequence contains exactly one complete gene. Note: This feature does not work properly in combination with alternative transcripts.</p>

<p><b>Ignore conflicts with other strand:</b> By default AUGUSTUS assumes that no genes - even on opposite strands - overlap. Indeed, this usually is the case but sometimes an intron contains a gene on the opposite strand. In this case, or when AUGUSTUS makes a false prediction on the one strand because it falsely thinks there is a conflicting gene on the other strand, AUGUSTUS should be run with this option set. It then predicts the genes on each strand separately and independently. This may lead to more false positive predictions, though.</p>

 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="list"><h2>What about that data duplication?</h2></div>
<p>We are trying to avoid data duplication. If you submitted some data that was already submitted before, by you or by somebody else, we will display a link to the previously submitted job.</p>
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="accuracy"><h2>Why is the prediction accuracy in the genome of my species not as good as I expected?</h2></div>
<p>Gene prediction accuracy of AUGUSTUS in the genome of a certain species depends on the quality of training genes that were used for optimizing species specific parameters. The pipeline behind our AUGUSTUS training web server application offers a fully automated way of generating training genes, but it does not replace manual quality checks on the training genes that are often needed for improving the training gene set quality.</p>

<p>In order to improve accuracy, you could manually inspect the generated training genes and select a trustworthy subset and try retraining AUGUSTUS with this subset. It also helps to compare the training gene set to other sources of evidence that are not supported by our web server application, e.g. RNA-seq data.</p>
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="privacy"><h2>What about data privacy and security</h2></div>
<p>The results of your job submission (i.e. in case of the training web server application that means log files, trained parameters, training genes, ab initio gene predictions and gene prediction with hints; or in case of the prediction web server application the augustus prediction archive) are publicly available. The link to your job status contains a long, pseudo-random string (uuid), and one needs to guess the string in order to get access to the results - but this is not particularly secure!</p>
<p>Other users who submit exactly the same input files as have been submitted before, will be redirected to the results page of the previously submitted job. They do not need to guess the link.</p>
<p>Files that you upload to our server, e.g. sequence files, are not directly made available to anyone. However, if you chose to upload a file via http/ftp link, the link to your file is displayed on the job status page.</p>
<p>We are interested in redistributing high quality parameter sets for novel species with the AUGUSTUS release. We will not do so without your explicit permission.</p>
<p>Our server logs e-mail adresses, IP adresses and all job submission details. We store this data for a limited time in order to be able to trace back errors or e.g. contact you about a permission to publish parameter sets. By submitting a job, you agree that we log this data.</p>
<p>Please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a> if your particular job requires a more secure environment, e.g. as part of a collaboration.</p>
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
	    <div id="results_pred"><h2>Prediction results</h2></div>
<p>After job computations have finished, you will receive an e-mail (if you supplied an e-mail adress). The job status web page may at this point in time look similar to this:</p>

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
<br>
            <b>Format example AUGUSTUS prediction gff file:</b>
            <pre class="example">
# This output was generated with AUGUSTUS (version 2.6).
# AUGUSTUS is a gene prediction tool for eukaryotes written by Mario Stanke (mario.stanke@uni-greifswald.de)
# and Oliver Keller (keller@cs.uni-goettingen.de).
# Please cite: Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008),
# Using native and syntenically mapped cDNA alignments to improve de novo gene finding
# Bioinformatics 24: 637-644, doi 10.1093/bioinformatics/btn013
# reading in the file /var/tmp/augustus/AUG-1855139717/hints.gff ...
# Setting 1group1gene for E.
# Sources of extrinsic information: M E 
# Have extrinsic information about 1 sequences (in the specified range). 
# Initialising the parameters ...
# human version. Use default transition matrix.
# Looks like /var/tmp/augustus/AUG-1855139717/input.fa is in fasta format.
# We have hints for 1 sequence and for 1 of the sequences in the input set.
#
# ----- prediction on sequence number 1 (length = 6483, name = HSACKI10) -----
#
# Delete group HintGroup , 5803-5803, mult= 1, priority= -1 1 features
# Forced unstranded hint group to the only possible strand for 3 groups.
# Deleted 1 groups because some hint was not satisfiable.
# Constraints/Hints:
HSACKI10	anchor	start	182	184	0	+	.	src=M
HSACKI10	anchor	stop	3058	3060	0	+	.	src=M
HSACKI10	anchor	dss	4211	4211	0	+	.	src=M
HSACKI10	b2h	ep	1701	2075	0	.	.	grp=154723761;pri=4;src=E
HSACKI10	b2h	ep	1716	2300	0	+	.	grp=13907559;pri=4;src=E
HSACKI10	b2h	ep	1908	2300	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	ep	3592	3593	0	+	.	grp=13907559;pri=4;src=E
HSACKI10	b2h	ep	3836	3940	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	ep	5326	5499	0	+	.	grp=27937842;pri=4;src=E
HSACKI10	b2h	ep	5805	6157	0	+	.	grp=27937842;pri=4;src=E
HSACKI10	b2h	exon	3142	3224	0	+	.	grp=13907559;pri=4;src=E
HSACKI10	b2h	exon	3142	3224	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	exon	3592	3748	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	anchor	intronpart	5000	5100	0	+	.	src=M
HSACKI10	b2h	intron	2301	3141	0	+	.	grp=13907559;pri=4;src=E
HSACKI10	b2h	intron	2301	3141	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	intron	3225	3591	0	+	.	grp=13907559;pri=4;src=E
HSACKI10	b2h	intron	3225	3591	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	intron	3749	3835	0	+	.	grp=154736078;pri=4;src=E
HSACKI10	b2h	intron	5500	5804	0	+	.	grp=27937842;pri=4;src=E
HSACKI10	anchor	CDS	6194	6316	0	-	0	src=M
HSACKI10	anchor	CDSpart	5900	6000	0	+	.	src=M
# Predicted genes for sequence number 1 on both strands
# start gene g1
HSACKI10	AUGUSTUS	gene	182	3060	0.63	+	.	g1
HSACKI10	AUGUSTUS	transcript	182	3060	0.63	+	.	g1.t1
HSACKI10	AUGUSTUS	start_codon	182	184	.	+	0	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	initial	182	225	1	+	0	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	internal	1691	2300	0.86	+	1	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	terminal	3049	3060	0.74	+	0	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	CDS	182	225	1	+	0	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	CDS	1691	2300	0.86	+	1	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	CDS	3049	3060	0.74	+	0	transcript_id "g1.t1"; gene_id "g1";
HSACKI10	AUGUSTUS	stop_codon	3058	3060	.	+	0	transcript_id "g1.t1"; gene_id "g1";
# coding sequence = [atgatgaaaccctgtctctaccaaaaagacaaaaaattagccagctcaagcaagcactactcttcctcccgcagtggag
# gaggaggaggaggaggaggatgtggaggaggaggaggagtgtcatccctaagaatttctagcagcaaaggctcccttggtggaggatttagctcaggg
# gggttcagtggtggctcttttagccgtgggagctctggtgggggatgctttgggggctcatcaggtggctatggaggattaggaggttttggtggagg
# tagctttcatggaagctatggaagtagcagctttggtgggagttatggaggcagctttggagggggcaatttcggaggtggcagctttggtgggggca
# gctttggtggaggcggctttggtggaggcggctttggaggaggctttggtggtggatttggaggagatggtggccttctctctggaaatgaaaaagta
# accatgcagaatctgaatgaccgcctggcttcctacttggacaaagttcgggctctggaagaatcaaactatgagctggaaggcaaaatcaaggagtg
# gtatgaaaagcatggcaactcacatcagggggagcctcgtgactacagcaaatactacaaaaccatcgatgaccttaaaaatcagagaacaacataa]
# protein sequence = [MMKPCLYQKDKKLASSSKHYSSSRSGGGGGGGGCGGGGGVSSLRISSSKGSLGGGFSSGGFSGGSFSRGSSGGGCFGG
# SSGGYGGLGGFGGGSFHGSYGSSSFGGSYGGSFGGGNFGGGSFGGGSFGGGGFGGGGFGGGFGGGFGGDGGLLSGNEKVTMQNLNDRLASYLDKVRAL
# EESNYELEGKIKEWYEKHGNSHQGEPRDYSKYYKTIDDLKNQRTT]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 20
# CDS exons: 1/3
#      E:   1 
# CDS introns: 0/2
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/0
# hint groups fully obeyed: 0
# incompatible hint groups: 5
#      E:   3 (gi|154723761,gi|13907559,gi|154736078)
#      M:   2 
# end gene g1
###         </pre>
<br>
Different kinds of information are printed after the hash signs, e.g. the applied AUGUSTUS version and parameter set, predicted coding sequence and amino acid sequence. Predictions and hints are given in tabulator separated gff format, i.e. the first column contains the target sequence, second column contains the source of the feature, third column contains the feature, forth column contains the feature start, fifth column contains the feature end, sixth column contains a score (if applicable), seventh column contains the strand, eightth column contains the reading frame and nineth column contains either for hints the grouping and source information, or for prediction lines the gene/transcript identifier.
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

<p><a href="http://bioinf.uni-greifswald.de/webaugustus/prediction/show/ff80818136a76dad0136a76fb00b0002">Click here</a> to view a real AUGUSTUS prediction web service output!</p>

    <p>It is important that you check the results of an AUGUSTUS gene prediction run. Do not trust predictions blindly! Prediction accuracy depends on the input sequence quality, on hints quality and on whether a given parameter set fits to the species of the supplied genomic sequence.</p>	    
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
	    </p>
            <hr>
            <br>
            <div id="results_train"><h2>Training results</h2></div>
	    <p>You find a detailed description of training results by <a href="trainingtutorial.gsp#results">clicking here</a>. To view a sample output, <a href="http://bioinf.uni-greifswald.de/webaugustus/training/show/ff80818136a76dad0136a76edf560001">click here</a>!</p>
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>

            <div id="commercial"><h2>I am not from academia/non-profit. What can I do?</h2></div>
		<p>Users who are not from academia or a non-profit organisation, and who are not using our web application for personal purposes, only, have the following options:
<ul>
<li>Run the training web server application with a genome file and an externally created training gene file</li>
<li>Run AUGUSTUS predictions ab initio or with an externally created hint file</li>
<li>Purchase a BLAT license from <a href="http://www.kentinformatics.com/">http://www.kentinformatics.com/</a> and run the autoAug Pipeline locally</li>
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

            <div id="dog"><h2>Why do I see a running dog when pressing the submission button?</h2></div>
<p>As Loriot said (freely translated): <i>Life without a dog is possible, but pointless.</i> ... the animation is simply displayed to make the waiting time during job submission more pleasant ;-)</p>
 <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
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
