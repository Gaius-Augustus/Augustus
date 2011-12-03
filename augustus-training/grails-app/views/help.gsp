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
        Bioinformatics Web Server
      </div>
      <div id="bannertitel2">
        AUGUSTUS Training
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
         <li><a href="index.gsp"><span>Introduction</span></a></li>
         <li><g:link controller="training" action="create"><span>Submitt Training</span></g:link></li>
         <li><g:link controller="prediction" action="create"><span>Submitt Prediction</span></g:link></li>
         <li id="current"><a href="help.gsp"><span>Help</span></a></li>
         <li><a href="references.gsp"><span>Links & References</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de"><span>Bioinformatics Group</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/impressum.html"><span>Impressum</span></a></li>
     </ul>
  </div>
 <div id="mittel_spalte">
<div class="main" id="main">
   <h1><a href="help.gsp">Help</a></h1>
            <div class="main" id="main">
            <p>This website contains instructions and frequently asked questions concerning <ul><li>the training of AUGUSTUS and <li>predicting genes in a new genome with pre-trained parameters.</ul></p>
            <hr>
            <br>
	    <h2>Contents</h2>
	    <p>
	      <a href="#noResults">Why do I not get any results?</a><br>
	      <a href="#buisy">Why is the server buisy?</a><br>
              <a href="#species_name">What is the species name?</a><br>
	      <a href="#email">Why must I give my e-mail address?</a><br>
              <a href="#upload_link">File upload versus web link</a><br>
              <a href="#which_files">Which files must or can I submitt for training AUGUSTUS?</a><br>
              <a href="#which_files_pred">Which files are required for predicting genes in a new genome?</a><br>
              <a href="#genome_file">Genome file</a><br>
              <a href="#cDNA">cDNA file</a><br>
              <a href="#protein">Protein file</a><br>
              <a href="#structure">Training gene structure and hint file</a><br>
              <a href="#archive">Parameter archive</a><br>
              <a href="#project_id">What is the project identifier?</a><br>
              <a href="#job_status">What does my job status mean?</a><br>
              <a href="#utr">Why predicting the UTR?</a><br>
              <a href="#list">Why does my job not exist?</a><br>
            </p>
            <hr>
	    <br>
	    <h2 id="noResults">Why do I not get any results?</h2>
	    <p>The quality of results depends on the quality and combination of you input data. If the input data did e.g. not provide sufficient information for generating training genes, then no AUGUSTUS paramteres will be optimized for your species, and no predictions will be made. In case of the gene prediction web server application, it is also possible that your submitted genome sequence does not contain any protein coding genes.</p>
<p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <h2 id="buisy">Why is the server buisy?</h2>
            <p>Training AUGUSTUS is a very resource and time consuming process. We use a grid engine queuing system with a limited number of waiting slots. If we estimate that the time from job submission to computation start might be very long, our web server might display a message that our server is buisy. The submission of new jobs is then disabled. Please wait one or two weeks before you try a new submission. If the problem persists longer than a month, please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a>.</p>
	    <p><a href="#seitenanfang">
	      <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
	      Top of page
	    </a>
	    <br>
	    </p>
            <hr>
	    <br>
            <div id="species_name"><h2>What is the species name?</h2></div>
                <p>
                The species name is the name of the species for whose genome you want to train AUGUSTUS. If you do not want to reveal the true species name, you may use any other string shorter than 30 characters as a species name. Species names must be unique, i.e. if the string of your choice is already existing in our system, you will get a message that you have to choose another species name.</p>
               <p>We are not redistributing the original sequence data that you submitted to our web server application. However, we are redistributing the trained parameters and the species name.</p><p><b>Example:</b><br> If person 1 submitts a sequence data set for training and names it <i>hypothetical_species</i>, and a second person tries later to train AUGUSTUS with exactly the same sequence files (new data upload) and names the species <i>some_other_name</i>, the second person will be redirected to the results of the original training parameter results of person 1, and the species name <i>hypothetical_species</i> will be publicly readable.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="email"><h2>Why must I give my e-mail address?</h2></div>
            <p>We need your e-mail address in order to send you a message with the link to the results of your computation when your job is finished. This may take a while (up to several weeks). We do not use your e-mail address for any other purposes.</p>
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
            <div id="which_files"><h2>Which files must or can I submitt for training AUGUSTUS?</h2></div>
            <p>You need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li><b>at least one out of the following files:</b>  <a  href="#cDNA">cDNA file</a>, <a href="#structure">training gene structure file</a>, and <a href="#protein">protein file</a>.</li>
            </ul>
	    </p>
	    <p>
            Please consider that training AUGUSTUS is a time and resource consuming process. For optimal results, you should specify as much information as possible for a single training run instead of starting the AUGUSTUS training multiple times with different file combinations!</p>
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
               <li>a <a href="#archive">parameter archive</a>. Instead of uploading the archive, you may also enter a valid <a href="#project_id">project identifier</a> in case you trained AUGUSTUS on this web server and the training has already finished.
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
            <p>The genome file is an obligatory file for training AUGUSTUS and for making predictions with pre-trained parameters in a new genome. It must contain the genome in (multiple) fasta format. Every header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver applications.<br><br>
            <b>Correct file format example:</b>
            <pre class="example">
> Chr. 1
CCTCCTCCTGTTTTTCCCTCAATACAACCTCATTGGATTATTCAATTCAC
CATCCTGCCCTTGTTCCTTCCATTATACAGCTGTCTTTGCCCTCTCCTTC
TCTCGCTGGACTGTTCACCAACTCTCAGCCCGCGATCCCAATTTCCAGAC
AACCCATCTTATCAGCTTGGCCACGGCCTCGACCCGAACAGACCGGCGTC
CAGCGAGAAGAGCGTCGCCTCGACGCCTCTGCTTGACCGCACCTTGATGC
TCAAGACTTATCGCGATGCCAAGAAGCGTCTCATCATGTTCGACTACGA
> Chr 2
CGAAACGGGCACCTATACAACGATTGAAACCATTATTCAAGCTCAGCAAG
CGTCTATGCTAGCGGTTATTGCGAGCACTTCAGCGGTTGCTACTACGACT
ACTACTTGATAAATGAAACGGCTATAAAAGAGGCTGGGGCAAAAGTATGT
TAGTTGAAGGGTGACCTGAACGATGAATCGGTCGAATTTTTTATTGGCAG
AGGGAAGGTAGGTTTACTCAATTTAGTTACTTCTAGCCGTTGATTGGAGG
AGCGCAAGCGACGAGGAGGCTCATCGGCCGCCCGCGGAAAGCGTAGTCT
TACACGGAAATCAACGGCGGTGTCATAAGCGAG
> Chr 3
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
            <p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. See <a href="#genome_file">Genome file</a> for a format example.</p>
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
> protein 1
maaaafgqlnleepppiwgsrsvdcfekleqigegtygqvymakeiktgeivalkkirmd
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqvlhrdikgsnllidnegnlkl
adfglarsyshdhtgnltnrvitlwyrppelllgatkygp
>protein 2
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqv
>protein 3
...
            </pre></p>
           <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="structure"><h2>Training gene structure and hints file</h2></div>
	    <p>You can submitt your own, externally created training gene structures to the AUGUSTUS training web server application.</p>
	    <p>For the gene prediction web server application, it is possible to submitt an externally created file that contains extrinsic evidence for gene structures.</p>
	    <p>Those two files both have to be in gff-format with the following (obligatory) tabulator-separated entries per line:<br><b>seqname source feature start end score strand frame attributes comments</b><br></p>
            <h3>Important details about gff-format for AUGUSTUS training</h3>
	    <p>
            <ul>
            <li>The columns are separated by <b>tabulators</b>. 
            <li>The <b>seqnames</b> in column 1 must be found in the fasta headers of sequences in the genome file. 
            <li>The <b>feature</b> (column 3) may be one of start, stop, exonpart, exon, dss, ass, or intronpart. 
            <li><b>Start</b> is the beginning position of the line's feature, counting the first position of a sequence as position 1.
            <li><b>Stop</b> position, must be at least as large as start position.
            <li>The <b>score</b> must be a number but the number is irrelevant to our web server applications.
            <li><b>Frame</b> is the reading frame, can be denoted as '.' if unknown or irrelevant. For exonpart and exon this is as defined as follows: On the forward strand it is the number of bases after (begin position 1) until the next codon boundary comes (0, 1 or 2). On the reverse strand it is the number of bases before (end position + 1) the next codon boundary comes (0, 1 or 2).
            <li><b>Attribute</b> must contain the string <b>source=M</b> (for manual). Other sources, such EST or protein, are possible, but only in the command line version of AUGUSTUS. Source types other than <b>M</b> are ignored by AUGUSTUS web server applications.
            </ul>
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
               <li><i>species</i>/<i>species</i>_parameters.cfg</li>
               <li><i>species</i>/<i>species</i>_exon_probs.pbl</li>
               <li><i>species</i>/<i>species</i>_igenic_probs.pbl</li>
               <li><i>species</i>/<i>species</i>_intron_probs.pbl</li>
               <li><i>species</i>/<i>species</i>_weightmatrix.txt</li>
            </ul>
	    </p>
            <p>where <i>species</i> is replaced by the name of the species you trained AUGUSTUS for (e.g. <i>carrot</i>). If you simply tar and gzip the folder that contains parameters of an AUGUSTUS training run, everything should work fine.</p>
                       <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="project_id"><h2>What is the project identifier?</h2></div>
            <p>If you trained AUGUSTUS on this webserver, you may instead of uploading a parameter archive, simply specify the project identifier of this training run. You find the project identifier for example in the subject line for your training confirmation e-mail, where it says <i>Your AUGUSTUS training job project_id</i>.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="job_status"><h2>What does my job status mean?</h2></div>
            <p>Your job is in one of five stages:</p>
	    <p>
            <ul>
               <li><b>stage 1</b>: submitted to webserver but not to cluster, yet. This may have one of the following reasons:
               <ul>
                 <li> file upload is still buisy
                 <li> file format validation is still in progress
                 <li>currently not enough computational resources available 
               </ul>
               <li><b>stage 2</b>: submitted to cluster and waiting for execution
               <li><b>stage 3</b>: calculating<br>
               This may take up to 3 weeks.
               <li><b>stage 4</b>: finished 
               <li><b>stage 5</b>: an error occurred
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
            <div id="utr"><h2>Why predicting the UTR?</h2>
            <p>If a model for the untranslated regions (UTRs) is available for the species, they are included in the prediction. It takes significantly more time but in addition to reporting UTRs, it usually is also a little more accurate on the coding regions when ESTs are given as input.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="list"><h2>Why does my job not exist?</h2></div>
            <p>We are trying to avoid data duplication. If you submitted some data that was already submitted before, by you or somebody else, we will delete your job. You receive an e-mail with a link to the job-status and results of training AUGUSTUS on your data. The results link will only be functional in case the computations of that previously submitted job have already finished.</p>

            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
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
