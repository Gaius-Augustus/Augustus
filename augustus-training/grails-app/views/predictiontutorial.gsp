

<!DOCTYPE html>
<html lang="de">
   <head>
      <meta http-equiv="content-type" content="text/html;charset=utf-8" />
      <meta charset="utf-8">
      <meta name="robots" content="INDEX,FOLLOW">
      <meta name="revisit-after" content="7 days">
      <meta name="abstract" content="Bioinformatics Greifswald">
      <meta name="keywords" content="Bioinformatics Greifswald">
      <meta name="description" content="University of Greifswald">
      <meta property="author" content="University of Greifswald">
      <meta name="date" content="2018-07-17">
      <link rel="stylesheet" href="${resource(dir: 'css', file: 'new1.css')}" type="text/css">
      <link rel="stylesheet" href="${resource(dir: 'css', file: 'new2.css')}" type="text/css">
      <title>Bioinformatics Web Server - University of Greifswald</title>
      <meta name="lastModified" content="2018-07-16">
      <meta http-equiv="X-UA-Compatible" content="IE=edge">
      <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">
      <style>.stickyFixed {position: fixed !important;  top:0 !important;} </style>
   </head>
   <body id="page-6289" class="">
      <!-- dark header topbar -->
      <div class="topbar">
      </div>
      <!-- header with uni logo -->
      <header class="header">
         <div class="header__content">
            <div class="header__top-wrapper">
               <!-- left side banner white spacer -->
               <div class="header__submenu"></div>
               <div class="logo">
                  <a href="https://www.uni-greifswald.de/" title="Universität Greifswald" class="logo-main">
                  <img src="http://bioinf.uni-greifswald.de/bioinf/img/uni-greifswald_opt.svg" width="400" height="118"   alt="Universität Greifswald" title="Universität Greifswald" >
                  </a>
               </div>
               <!-- middle part of header -->
               <div class="organization">
                  <a href="http://bioinf.uni-greifswald.de/">
                     <h3>Bioinformatics Web Server</h3>
                  </a>
               </div>
            </div>
            <nav id="nav" class="navigation">
               <ul class="navigation-list navigation-list--table">
                  <li class="navigation-list__item navigation-list__item--level-1 navigation-list__item--active" data-dropdown="true"><a href="http://bioinf.uni-greifswald.de/">Bioinformatics Group</a></li>
                  <li class="navigation-list__item navigation-list__item--level-1" data-dropdown="true"><a href="http://math-inf.uni-greifswald.de/">Mathematics and Computer Science</a></li>
                  <li class="navigation-list__item navigation-list__item--level-1" data-dropdown="true"><a href="https://mnf.uni-greifswald.de/en/faculty/">Faculty of Math and Natural Sciences</a></li>
               </ul>
            </nav>
         </div>
      </header>
      <div class="container">
         <div class="grid">
            <div class="column-1 grid__column grid__column--md-3">
               <ul class="navigation-sub">
                  <li class="navigation-sub__item">
                     <span class="navigation-sub__headline">AUGUSTUS Web Server Navigation</span>
                     <ul class="navigation-sub">
                        <li class="navigation-sub__item"><a href="index.gsp">Introduction</a></li>
                        <li class="navigation-sub__item"><a href="about.gsp">About AUGUSTUS</a></li>
                        <li class="navigation-sub__item"><a href="accuracy.gsp">Accuracy</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="training" action="create">Training Tutorial</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="index.gsp">Submit Training</a></li>
                        <li class="navigation-sub__item"><a href="predictiontutorial.gsp">Prediction Tutorial</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="prediction" action="create">Submit Prediction</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="datasets.gsp">Datasets for Download</a></li>
                        <li class="navigation-sub__item"><a href="references.gsp">Links & References</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/bioinf/impressum.html">Impressum</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/bioinf/datenschutz.html">Data Privacy Protection</a></li>
                     </ul>
                  </li>
                  <li class="navigation-sub__item">
                     <span class="navigation-sub__headline">Other AUGUSTUS Resources</span>
                     <ul class="navigation-sub">
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus">AUGUSTUS Wiki</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/bioinf/forum">AUGUSTUS Forum</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/augustus/downloads/index.php">Download AUGUSTUS</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/augustus">Old AUGUSTUS web server</a></li>
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de/bioinf/braker">BRAKER</a></li>
                     </ul>
                  </li>
                  <li class="navigation-sub__item">
                     <span class="navigation-sub__headline">Other Links</span>
                     <ul class="navigation-sub">
                        <li class="navigation-sub__item"><a href="http://bioinf.uni-greifswald.de">Bioinformatics Greifswald</a></li>
                     </ul>
                  </li>
               </ul>
            </div>
            <div class="column-2 grid__column grid__column--md-9">
               <main class="main-content">
                  <div id="c180465" class="csc-default">
                     <div class="csc-header csc-header-n1">
                        <h1 class="csc-firstHeader">WebAUGUSTUS Prediction Tutorial</h1>
                     </div>
                  </div>
                  <div id="c261665" class="csc-default">
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
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
                     <div id="contents">
                        <h1><a href="#contents">Contents</a></h1>
                     </div>
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
                        <a href="#submitt">1.5 - The submitt button</a><br>
                        <a href="#exampledata">1.6 - Example data files</a><br><br>
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
                     <div id="job_submission">
                        <h1><a href="#job_submission">1 - Job submission in general</a></h1>
                     </div>
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
                     <div id="finding_form">
                        <h2><a href="#finding_form">1.1 - Finding the prediction submission form</a></h2>
                     </div>
                     <p>You find the AUGUSTUS prediction submission form by clicking on the following link in the left side navigation bar:</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/submitt-link-pred.jpg" alt="image of submission link"></td>
                        </tr>
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
                     <div id="general_data">
                        <h2><a href="#general_data">1.2 - Filling in general job data</a></h2>
                     </div>
                     <p>This section describes all fields that should be filled in for every job submission, i.e. fields that are obligatory (except for the email adress, which is optional but strongly recommended).</p>
                     <div id="email">
                        <h3><a href="#email">1.2.1 - E-mail address</a></h3>
                     </div>
                     <p>At first, we recommend that you enter a <b>valid e-mail address</b>:</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/email-prediction.jpg" alt="image of e-mail address field"></td>
                        </tr>
                     </table>
                     </p>
                     <p>It is possible to run AUGUSTUS without giving an e-mail adress but here are some reasons why we recommend supplying an e-mail adress:
                     <ul>
                        <li>Unlike many other bioinformatics web services, the AUGUSTUS web server application is not an implementation of a fail-safe procedure. Our pipeline may issue warnings or errors, and sometimes, we need to get some feedback from you via e-mail in order to figure out what is the problem with your particular input data set.</li>
                        <li>In addition, running AUGUSTUS on large files is rather time consuming processe that may take up to several weeks (depending on the input data). It may be more convenient to receive an e-mail notification about your job having finished, than checking the status page over and over, again</li>
                     </ul>
                     </p>
                     <div id="email_purpose">
                        <h4><a href="#email_purpose">1.2.1.1 - What your e-mail address is used for</a></h4>
                     </div>
                     <p> We use your e-mail address for the following purposes:</p>
                     <p>
                     <ul>
                        <li>Confirming your job submission</li>
                        <li>Confirming successful file upload (for large files via ftp/http link)</li>
                        <li>Notifying you that your job has finished</li>
                        <li>Informing you about any problems that might occur during your particular AUGUSTUS prediction job</li>
                        <li>Asking you for permission to publish parameters with the next AUGUSTUS release</li>
                     </ul>
                     </p>
                     <p>We do <b>not</b> use your e-mail address to send you any <i>spam</i>, i.e. about web service updates. We do not share your e-mail addresses with any third parties.</p>
                     <p>Job submission without giving an email adress is possible but discouraged for large input files.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="params">
                        <h3><a href="#params">1.2.2 - AUGUSTUS species parameters</a></h3>
                     </div>
                     <p>The web server application offers you three options to specify which parameter set you want to use for predicting genes with AUGUSTUS. You can either uploaded a <tt>*.tar.gz</tt> parameter archive from your local harddrive, or you can specify the job ID of a previously finished AUGUSTUS web server application training run, or you can select a pre-trained parameter set through the drop-down menu.</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/parameters.jpg" alt="image of parameter upload field"></td>
                        </tr>
                     </table>
                     </p>
                     <div id="param_archive">
                        <h3>
                           <a href="#param_archive">
                              1.2.2.1 - Uploading an archive file</a
                        </h3>
                     </div>
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
                              <td><font size="1"><b>Species</b></font></td>
                              <td><font size="1"><b>Project identifier</b></font></td>
                              <td><font size="1"><b>Courtesy of</b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Animals</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Acyrthosiphon pisum</font></td>
                              <td><font size="1">pea_aphid</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Aedes aegypti</font></td>
                              <td><font size="1">aedes</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Amphimedon queenslandica</font></td>
                              <td><font size="1">amphimedon</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Apis dorsata</font></td>
                              <td><font size="1">adorsata</font></td>
                              <td><font size="1">Francisco Camara Ferreira</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Apis mellifera</font></td>
                              <td><font size="1">honeybee1</font></td>
                              <td><font size="1">Katharina Hoff and Mario Stanke</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Bombus terrestris</font></td>
                              <td><font size="1">bombus_terrestris2</font></td>
                              <td><font size="1">Katharina Hoff</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Brugia malayi</font></td>
                              <td><font size="1">brugia</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Caenorhabditis elegans</font></td>
                              <td><font size="1">caenorhabditis</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Callorhinchus milli</font></td>
                              <td><font size="1">elephant_shark</font></td>
                              <td><font size="1">Tereza Manousaki and Shigehiro Kuraku</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Camponotus floridanus</font></td>
                              <td><font size="1">camponotus_floridanus</font></td>
                              <td><font size="1">Shishir K Gupta</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Danio rerio</font></td>
                              <td><font size="1">zebrafish</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Drosophila melanogaster</font></td>
                              <td><font size="1">fly</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Heliconius melpomene</font></td>
                              <td><font size="1">heliconius_melpomene1</font></td>
                              <td><font size="1">Sebastian Adler and Katharina Hoff</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Gallus gallus domesticus</font></td>
                              <td><font size="1">chicken</font></td>
                              <td><font size="1">Stefanie Koenig</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Homo sapiens</font></td>
                              <td><font size="1">human</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Nasonia vitripennis</font></td>
                              <td><font size="1">nasonia</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Petromyzon marinus</font></td>
                              <td><font size="1">lamprey</font></td>
                              <td><font size="1">Falk Hildebrand and Shigehiro Kuraku</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Schistosoma mansoni</font></td>
                              <td><font size="1">schistosoma</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Tribolium castaneum</font></td>
                              <td><font size="1">tribolium</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Trichinella spiralis</font></td>
                              <td><font size="1">trichinella</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Alveolata</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Tetrahymena thermophila</font></td>
                              <td><font size="1">tetrahymena</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Toxoplasma gondii</font></td>
                              <td><font size="1">toxoplasma</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Protozoa</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Leishmania tarantolae</font></td>
                              <td><font size="1">leishmania_tarentolae</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Plants and algae</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Arabidopsis thaliana</font></td>
                              <td><font size="1">arabidopsis</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Chlamydomonas reinhartii</font></td>
                              <td><font size="1">chlamy2011</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Galdieria sulphuraria</font></td>
                              <td><font size="1">galdieria</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Triticum aestivum</font></td>
                              <td><font size="1">wheat</font></td>
                              <td><font size="1">Stefanie K&ouml;nig</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Zea mays</font></td>
                              <td><font size="1">maize</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Solanum lycopersicum</font></td>
                              <td><font size="1">tomato</font></td>
                              <td><font size="1"></font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Fungi</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Aspergillus fumigatus</font></td>
                              <td><font size="1">aspergillus_fumigatus</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Aspergillus nidulans</font></td>
                              <td><font size="1">aspergillus_nidulans</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Aspergillus oryzae</font></td>
                              <td><font size="1">aspergillus_oryzae</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Aspergillus terreus</font></td>
                              <td><font size="1">aspergillus_terreus</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Botrytis cinerea</font></td>
                              <td><font size="1">botrytis_cinerea</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Candida albicans</font></td>
                              <td><font size="1">candida_albicans</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Candida guilliermondii</font></td>
                              <td><font size="1">candida_guilliermondii</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Candida tropicalis</font></td>
                              <td><font size="1">candida_tropicalis</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Chaetomium globosum</font></td>
                              <td><font size="1">chaetomium_globosum</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Coccidioides immitis</font></td>
                              <td><font size="1">coccidioides_immitis</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Coprinus cinereus</font></td>
                              <td><font size="1">coprinus</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Cryptococcus neoformans</font></td>
                              <td><font size="1">cryptococcus_neoformans_neoformans_B</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Debaryomyces hansenii</font></td>
                              <td><font size="1">debaryomyces_hansenii</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Encephalitozoon cuniculi</font></td>
                              <td><font size="1">encephalitozoon_cuniculi_GB</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Eremothecium gossypii</font></td>
                              <td><font size="1">eremothecium_gossypii</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Fusarium graminearum</font></td>
                              <td><font size="1">fusarium_graminearum</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Histoplasma capsulatum</font></td>
                              <td><font size="1">histoplasma_capsulatum</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Kluyveromyces lactis</font></td>
                              <td><font size="1">kluyveromyces_lactis</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Laccaria bicolor</font></td>
                              <td><font size="1">laccaria_bicolor</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Lodderomyces elongisporus</font></td>
                              <td><font size="1">lodderomyces_elongisporus</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Magnaporthe grisea</font></td>
                              <td><font size="1">magnaporthe_grisea</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Neurospora crassa</font></td>
                              <td><font size="1">neurospora_crassa</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Phanerochaete chrysosporium</font></td>
                              <td><font size="1">phanerochaete_chrysosporium</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Pichia stipitis</font></td>
                              <td><font size="1">pichia_stipitis</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Phizopus oryzae</font></td>
                              <td><font size="1">rhizopus_oryzae</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Saccharomyces cerevisiae</font></td>
                              <td><font size="1">saccharomyces_cerevisiae_S288C</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Schizosaccharomyces pombe</font></td>
                              <td><font size="1">schizosaccharomyces_pombe</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Ustilago maydis</font></td>
                              <td><font size="1">ustilago_maydis</font></td>
                              <td><font size="1">Jason Stajich</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Verticillium longisporum</font></td>
                              <td><font size="1">verticillium_longisporum1</font></td>
                              <td><font size="1">Katharina Hoff and Mario Stanke</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Yarrowia lipolytica</font></td>
                              <td><font size="1">yarrowia_lipolytica</font></td>
                              <td><font size="1">Jason Stajich, modified by Katharina Hoff</font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Archaea (experimental parameters)</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Sulfolobus solfataricus</font></td>
                              <td><font size="1">sulfolobus_solfataricus</font></td>
                              <td><font size="1">Katharina Hoff</font></td>
                           </tr>
                           <tr>
                              <td><font size="1"><b><i>Bacteria (experimental parameters)</i></b></font></td>
                              <td><font size="1"><b></b></font></td>
                              <td><font size="1"><b></b></font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Escherichia coli</font></td>
                              <td><font size="1">E_coli_K12</font></td>
                              <td><font size="1">Katharina Hoff</font></td>
                           </tr>
                           <tr>
                              <td><font size="1">Thermoanaerobacter tencongensis</font></td>
                              <td><font size="1">thermoanaerobacter_tengcongensis</font></td>
                              <td><font size="1">Katharina Hoff</font></td>
                           </tr>
                        </font>
                     </table>
                     </p>
                     <p>
                        Please let us know whether you want to have parameters that you trained for a certain species to be included in this public list! If they are included in this list, they will also be distributed with the upcoming AUGUSTUS release.
                     </p>
                     </p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="genome_file">
                        <h3><a href="#genome_file">1.2.3 - Genome file</a></h3>
                     </div>
                     <p>The genome file is an obligatory input for predicting genes with AUGUSTUS.</p>
                     <div id="genome_file_format">
                        <h4><a href="#genome_file_format">1.2.3.1 - Genome file format</a></h4>
                     </div>
                     The genome file must contain the genome sequence in (multiple) fasta format. Every <i>unique</i> header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver application. White spaces in the sequence header might cause problems if the first word after the leading character <b>></b> is identical for several fasta entries. We generally recommend short, unique, non-white-space containing fasta headers.<br><br>
                     <b>Correct file format example:</b>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td>
                              <pre class="example"><font size="1">
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
            </font></pre>
                           </td>
                        </tr>
                     </table>
                     </p>
                     <% /* comment to break the static gsp block */ %>
                     <p>The maximal number of scaffolds allowed in a genome file is 250000. If your file contains more scaffolds, please remove all short scaffolds. For training AUGUSTUS, short scaffolds are worthless because no complete training genes can  be generated from them. In terms of prediction, it is possibleto predict genes in short scaffolds. However, those genes will in most cases be incomplete and probably unreliable.</p>
                     <p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for genome file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>.</p>
                     <div id="genome_file_upload">
                        <h4><a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a></h4>
                     </div>
                     <p>The AUGUSTUS prediction web server application offers two possiblities for transferring the genome file to the server: <i>Upload a file</i> and <i>specify a web link to file</i>.</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/training-genome.jpg" alt="image of genome submission field"></td>
                        </tr>
                     </table>
                     </p>
                     <p>
                     <ul>
                        <li>For <b>small files</b>, please click on the <i>Choose File</i> or <i>Browse</i>-button and select a file on your harddrive.<br>If you experience a <i>Connection timeout</i> (because your file was too large for this type of upload - the size is browser dependent), please use the option for large files!</li>
                        <li><b>Large files</b> can be retrieved from a <b>public</b> web link. Deposit your sequence file at a http or ftp server and specify the valid URL to your sequence file in the training submission form. Our server will fetch the file from the given address upon job submission. (File size limit: currently 1 GB. Please contact us in case you want to upload a bigger genome file, <b>links to dropbox are not supported by WebAUGUSTUS</b>.) You will be notified by e-mail when the file upload from web-link is finished (i.e. you can delete the file from the public server after you received that e-mail).</li>
                     </ul>
                     </p>
                     <p>
                        <b>You cannot do both at the same time!</b> You must <b>either</b> select a file on your harddrive <b>or</b> give a web link!
                     </p>
                     <div id="genome_file_purpose">
                        <h4><a href="#genome_file_purpose">1.2.3.3 - What the genome file is used for</a></h4>
                     </div>
                     <p>The genome file is used as a template for gene prediction, it is the sequence in which you want to predict genes.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="optional">
                        <h2><a href="#optional">1.3 - Optional fields</a></h2>
                     </div>
                     <p>This section describes a number of fields that are optional for predicting genes with AUGUSTUS.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="cDNA">
                        <h3><a href="#cDNA">1.3.1 - cDNA file</a></h3>
                     </div>
                     <p>This feature is available only <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> as this is required by the BLAT license.</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/training-cDNA.jpg" alt="image of cDNA submission field"></td>
                        </tr>
                     </table>
                     </p>
                     <div id="cDNA_format">
                        <h3><a href="#cDNA_format">1.3.1.1 - cDNA file format</a></h3>
                     </div>
                     <p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. An example for correct cDNA file format is given at <a href="#genome_file_format">1.2.3.1 - Genome file format</a>.</p>
                     <p>It is currently possible to submitt assembled RNA-seq transcripts instead of or mixed with ESTs as a cDNA/EST file. However, you should be aware that RNA-seq files are often much bigger than EST or cDNA files, which increases runtime of a prediction job. In order to keep runtime of your prediction job as low as possible, you should remove all assembled RNA-seq transcripts from your file that do not map to the submitted genome sequence. (In principle, this holds true for EST and cDNA files, too, but there, the problem is not as pronounced due to a smaller number of sequences.)</p>
                     <p>It is currently not allowed to upload RNA-seq raw sequences. (We filter for the average length of cDNA fasta entries and may reject the entire training job in case the sequences are on average too short, i.e. shorter than 400 bp.)</p>
                     <p>Besides plain fasta format, our server accepts <b>gzipped-fasta</b> format for cDNA file upload. You find more information about gzip at the <a href="http://www.gzip.org/">gzip homepage</a>. Gzipped files have the file ending <tt>*.gz</tt>. The maximal supported file size is 1 GB.</p>
                     <div id="cDNA_upload">
                        <h3><a href="#cDNA_upload">1.3.1.2 - cDNA file upload options</a></h3>
                     </div>
                     <p>There are two options for cDNA file upload: upload from your local harddrive, or upload from a public http or ftp server. Please see <a href="#genome_file_upload">1.2.3.2 - Genome file upload options</a> for a more detailed description of upload options.</p>
                     <div id="cDNA_purpose">
                        <h3><a href="#cDNA_purpose">1.3.1.3 - What cDNA files are used for</a></h3>
                     </div>
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
                     <div id="hints">
                        <h3><a href="#hints">1.3.2 - Hints file</a></h3>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/hints-field.jpg" alt="image of hints file upload field"></td>
                        </tr>
                     </table>
                     </p>
                     <div id="#hints_format">
                        <h4><a href="#structure_file_format">1.3.3.1 - Hints file format</a></h4>
                     </div>
                     <p>It is possible to submit an externally created file that contains extrinsic evidence for gene structures in gff format.</p>
                     <p>In general, gff files must contain the following columns (the columns are separated by <b>tabulators</b>):</p>
                     <p>
                     <OL TYPE="1">
                        <li>The <b>sequence names</b> must be found in the fasta headers of sequences in the genome file. 
                        <li>The <b>source</b> tells with which software/process the gene structure was generated (you can fill in whatever you like).
                        <li>
                           The <b>feature</b> may for AUGUSTUS gene prediction be 
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
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td>
                              <font size="1">
                                 <pre class="example">
HS04636 anchor  exonpart        500     506     0       -       .       source=M
HS04636 anchor  exon            966     1017    0       +       0       source=M
HS04636 anchor  start           966     968     0       +       0       source=M
HS04636 anchor  dss             2199    2199    0       +       .       source=M
HS04636 anchor  stop            7631    7633    0       +       0       source=M
HS04636 anchor  intronpart      7631    7633    0       +       0       source=M
            </pre>
                              </font>
                           </td>
                        </tr>
                     </table>
                     </p>
                     <div id="hints_purpose">
                        <h4><a href="#hints_purpose">1.3.1.2 - What hints files are used for</a></h4>
                     </div>
                     <p>The hints file is used as extrinsic evidence that supports gene structure prediction. You can generate hints yourself based on any alignment program and information resource (e.g. ESTs, RNA-seq data, peptides, proteins, ...) that appears suitable to you.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="utr">
                        <h3><a href="#utr">1.3.3 - UTR prediction</a></h3>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/utr-checkbox.jpg" alt="image of utr checkbox"></td>
                        </tr>
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
                     <div id="strand">
                        <h3><a href="#strand">1.3.4 - Strand specific prediction</a></h3>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/strands.jpg" alt="image of the strand checkboxes"></td>
                        </tr>
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
                     <div id="alternative">
                        <h3><a href="#alternative">1.3.5 - Alternative transcripts</a></h3>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/alternative.jpg" alt="image of alternative transcript radio buttons"></td>
                        </tr>
                     </table>
                     </p>
                     <p>By default, AUGUSTUS does not predict any alternative transcripts.</p>
                     <p>If you select <b>few</b>, then the following AUGUSTUS parameters are set to result in the prediction of relatively few alternative transcripts:<br><tt>--alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 <br>--maxtracks=2</tt>.</p>
                     <p> If you select <b>medium</b> the AUGUSTUS parameters are set to <br><tt>--alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.4 <br>--maxtracks=3"</tt>.</p>
                     <p> If you select <b>many</b>, AUGUSTUS parameters are set to <br><tt>--alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 <br>--maxtracks=20"</tt>. </p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="allowed_structure">
                        <h3><a href="#allowed_structure">1.3.6 - Allowed gene structure</a></h3>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/allowed-structures.jpg" alt="image of allowed gene structure buttons"></td>
                        </tr>
                     </table>
                     </p>
                     <p><b>Predict any number of (possibly partial) genes:</b> This option is set by default. AUGUSTUS may predict no gene at all, one or more genes. The genes at the boundaries of the input sequence may be partial. Partial here means that not all of the exons of a gene are contained in the input sequence, but it is assumed that the sequence starts or ends in a non-coding region.
                     <p>
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
                     <div id="verification">
                        <h2><a href="#verification">1.4 - Verification that you are a human</a></h2>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/verification.jpg" alt="image of verification field"></td>
                        </tr>
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
                     <div id="submitt">
                        <h2><a href="#submitt">1.5 - The submitt button</a></h2>
                     </div>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/pred-submitt.jpg" alt="image of submission button"></td>
                        </tr>
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
                     <div id="exampledata">
                        <h2><a href="#exampledata">1.6 - Example data files</a></h2>
                     </div>
                     <p>In the following, we provide some correctly formatted, compatible example data files:</p>
                     <p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee1.tar.gz">http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee1.tar.gz</a> - This file is an example of a <b>AUGUSTUS species parameter archive file</b>. Please do <b>not</b> upload this archive to our server since the identical parameters are usable through the AUGUSTUS species parameter project identifier <b>honeybee1</b> and a re-upload would simply duplicate this data set. We only provide this file as an example which may help you check your own parameter archive in case incompatibilities with your application might occur. These parameters were optimized for predicting genes in <i>Apis mellifera</i>.</p>
                     <p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/LG16.fa">http://bioinf.uni-greifswald.de/trainaugustus/examples/LG16.fa</a> - This file may be used as a <b>Genome file</b>. It contains linkage group 16 of <i>Apis mellifera</i> from GenBank (modified headers).</p>
                     <p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee-ests.fa">http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee-ests.fa</a> - This file may be used as a <b>cDNA file</b>. It contains 3 ESTs of <i>Apis mellifera</i> from GenBank (modified headers).</p>
                     <p><a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee.hints">http://bioinf.uni-greifswald.de/trainaugustus/examples/honeybee.hints</a> - This file may be used as a <b>Hints file</b>. It contains hints that were generated from <i>Apis mellifera</i> RNA-Seq data for genome file <a href="http://bioinf.uni-greifswald.de/trainaugustus/examples/LG16.fa">LG16.fa</a>.</p>
                     <p>You can insert some of these sample data sets by pressing the "Fill in Sample Data" button:</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/train-sample.jpg" alt="image of sample button"></td>
                        </tr>
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
                     <div id="job_status">
                        <h1><a href="#job_status">2 - What happens after submission</a></h1>
                     </div>
                     <p>After you click the "Start Predicting" button, the web server application first validates whether the combination of your input fields is generally correct. If you submitted an unsupported input combination you will be redirected to the training submission form and an error message will be displayed at the top of the page.</p>
                     <p>If all fields were filled in correctly, the application is actually initiated. You will receive an e-mail that confirms your job submission and that contains a link to the job status page (if you supplied an e-mail adress). You will be redirected to the job status page.</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/prediction-job-status.jpg" alt="image of job status page"></td>
                        </tr>
                     </table>
                     </p>
                     <p>In the beginning, the status page will display that your job has been <b>submitted</b>. This means, the web server application is currently uploading your files and validating file formats. After a while, the status will change to <b>waiting for execution</b>. This means that all file formats have been confirmed and an actually AUGUSTUS training job has been submitted to our grid engine, but the job is still pending in the queue. Depending on waiting queue length, this status may persist for a while. Please contact us in case you job is pending for more than one month. Later, the job status will change to <b>computing</b>. This means the job is currently computing. When the page displays <b>finished</b>, all computations have been finished and a website with your job's results has been generated.</p>
                     <p>You will receive an e-mail when your job has finished (if you supplied an e-mail adress).</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="duplication">
                        <h2><a href="#duplication">2.1 - Submission duplication</a></h2>
                     </div>
                     <p>Since predicting genes wiht AUGUSTUS may under certain circumstances be is a very resource consuming process, we try to avoid data duplication. In case you or somebody else tries to submitt exactly the same input file combination more than once, the duplicated job will be stopped and the submitter of the redundant job will receive information where the status page of the previously submitted job is located.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="error">
                        <h2><a href="#error">2.2 - Errors during prediction</a></h2>
                     </div>
                     <p>You should automatically receive an e-mail in case an error occurs during the AUGUSTUS gene prediction process. The admin of this server is also notified by e-mail about errors. We will get in touch with you, again, after we figured out what caused the error. If you did not supply an e-mail adress, errors are likely to be ignored by the AUGUSTUS webserver development team.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <hr>
                     <br>
                     <div id="results">
                        <h1><a href="#results">3 - Prediction Results</a></h1>
                     </div>
                     <p>After job computations have finished, you will receive an e-mail (if you supplied an e-mail adress). The job status web page may at this point in time look similar to this:</p>
                     <p>
                     <table border="2" cellspacing="0" cellpadding="0">
                        <tr>
                           <td><img src="images/prediction-results-example.jpg" alt="image of results example"></td>
                        </tr>
                     </table>
                     </p>
                     <p>This page should contain the file <b>augustus.tar.gz</b>. Please make a "right click" on the link and select "Save As" (or similar) to save the file on your local harddrive.</p>
                     <p><b>augustus.tar.gz</b> is a gene prediction archive and its content depends on the input file combination. You can unpack the archive by typing <tt>tar -xzvf *.tar.gz</tt> into your shell. (You find more information about the software tar at the <a href="http://www.gnu.org/s/tar/">GNU tar website</a>.)</p>
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
# Initializing the parameters ...
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
HSACKI10        anchor  start   182     184     0       +       .       src=M
HSACKI10        anchor  stop    3058    3060    0       +       .       src=M
HSACKI10        anchor  dss     4211    4211    0       +       .       src=M
HSACKI10        b2h     ep      1701    2075    0       .       .       grp=154723761;pri=4;src=E
HSACKI10        b2h     ep      1716    2300    0       +       .       grp=13907559;pri=4;src=E
HSACKI10        b2h     ep      1908    2300    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     ep      3592    3593    0       +       .       grp=13907559;pri=4;src=E
HSACKI10        b2h     ep      3836    3940    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     ep      5326    5499    0       +       .       grp=27937842;pri=4;src=E
HSACKI10        b2h     ep      5805    6157    0       +       .       grp=27937842;pri=4;src=E
HSACKI10        b2h     exon    3142    3224    0       +       .       grp=13907559;pri=4;src=E
HSACKI10        b2h     exon    3142    3224    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     exon    3592    3748    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        anchor  intronpart      5000    5100    0       +       .       src=M
HSACKI10        b2h     intron  2301    3141    0       +       .       grp=13907559;pri=4;src=E
HSACKI10        b2h     intron  2301    3141    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     intron  3225    3591    0       +       .       grp=13907559;pri=4;src=E
HSACKI10        b2h     intron  3225    3591    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     intron  3749    3835    0       +       .       grp=154736078;pri=4;src=E
HSACKI10        b2h     intron  5500    5804    0       +       .       grp=27937842;pri=4;src=E
HSACKI10        anchor  CDS     6194    6316    0       -       0       src=M
HSACKI10        anchor  CDSpart 5900    6000    0       +       .       src=M
# Predicted genes for sequence number 1 on both strands
# start gene g1
HSACKI10        AUGUSTUS        gene    182     3060    0.63    +       .       g1
HSACKI10        AUGUSTUS        transcript      182     3060    0.63    +       .       g1.t1
HSACKI10        AUGUSTUS        start_codon     182     184     .       +       0       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        initial 182     225     1       +       0       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        internal        1691    2300    0.86    +       1       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        terminal        3049    3060    0.74    +       0       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        CDS     182     225     1       +       0       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        CDS     1691    2300    0.86    +       1       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        CDS     3049    3060    0.74    +       0       transcript_id "g1.t1"; gene_id "g1";
HSACKI10        AUGUSTUS        stop_codon      3058    3060    .       +       0       transcript_id "g1.t1"; gene_id "g1";
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
                     <p>It is important thatyou check the results of an AUGUSTUS gene prediction run. Do not trust predictions blindly! Prediction accuracy depends on the input sequence quality, on hints quality and on whether a given parameter set fits to the species of the supplied genomic sequence.</p>
                     <p><a href="#seitenanfang">
                        <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
                        Top of page
                        </a>
                        <br>
                     </p>
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
                  </div>
               </main>
            </div>
         </div>
      </div>
      <footer class="footer footer--padding-bottom">
         <div class="footer-column footer-column--dark">
            <div class="footer__content-wrapper">
               <div class="footer-bottom">
                  <div class="footer-bottom__copyright">
                     <p>&copy;&nbsp;2018&nbsp; Universität Greifswald</p>
                  </div>
               </div>
            </div>
         </div>
      </footer>
   </body>
</html>

