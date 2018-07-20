

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
      <script type="text/javascript" src="js/md_stylechanger.js"></script> 
      <!-- flooble Expandable Content header start -->
      <script language="javascript">
         // Expandable content script from flooble.com.
         // For more information please visit:
         //   http://www.flooble.com/scripts/expand.php
         // Copyright 2002 Animus Pactum Consulting Inc.
         // Script was customized for this application by author of this HTML document!
         //----------------------------------------------
         var ie4 = false; if(document.all) { ie4 = true; }
         function getObject(id) { if (ie4) { return document.all[id]; } else { return document.getElementById(id); } }
         function toggle(link, divId) { var lText = link.innerHTML; var d = getObject(divId);
         if (lText == 'click to expand') { link.innerHTML = 'click to minimize'; d.style.display = 'block'; }
         else { link.innerHTML = 'click to expand'; d.style.display = 'none'; } }
      </script>
      <!-- flooble Expandable Content header end   -->    
      <script type="text/javascript">
         <!--
             function toggle_visibility(id) {
                var e = document.getElementById(id);
                if(e.style.display == 'block')
                   e.style.display = 'none';
                else
                   e.style.display = 'block';
             }
         //-->
      </script>                    
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
                        <li class="navigation-sub__item"><a href="../index.gsp">Introduction</a></li>
                        <li class="navigation-sub__item"><a href="../about.gsp">About AUGUSTUS</a></li>
                        <li class="navigation-sub__item"><a href="../accuracy.gsp">Accuracy</a></li>
                        <li class="navigation-sub__item"><a href="../trainingtutorial.gsp">Training Tutorial</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="training" action="create">Submit Training</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="../predictiontutorial.gsp">Prediction Tutorial</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="prediction" action="create">Submit Prediction</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="../datasets.gsp">Datasets for Download</a></li>
                        <li class="navigation-sub__item"><a href="../references.gsp">Links & References</a></li>
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
               <g:if test="${flash.message}">
                  <div class="message">${flash.message}</div>
               </g:if>
               <g:hasErrors bean="${predictionInstance}">
                  <div class="errors">
                     <g:renderErrors bean="${predictionInstance}" as="list" />
                  </div>
               </g:hasErrors>
               <g:if test="${flash.error}">
                  <div class="errors">
                     &nbsp; <img src="../images/skin/exclamation.png"> &nbsp; ${flash.error}
                  </div>
               </g:if>
               <main class="main-content">
                  <div id="c180465" class="csc-default">
                     <div class="csc-header csc-header-n1">
                        <h1 class="csc-firstHeader">Data Input for Running AUGUSTUS</h1>
                     </div>
                  </div>
                  <div id="c261665" class="csc-default">
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
                     <noscript>
                        <p><b><span style="color:red">Please enable javascript in your browser in order to display the submission form correctly!</span></b> Form functionality is not affected significantly while javascript is disabled, but it looks less pretty.</p>
                     </noscript>
                     <g:uploadForm action="commit" method="post" name="submissionform">
                        <fieldset>
                        <p>
                        <div class="dialog">
                           <p>Use this form to submit your data for running AUGUSTUS on new genomic data with already available pre-trained parameters.</p>
                           <p>Please read the <a href="../predictiontutorial.gsp">prediction tutorial</a> before submitting a job for the first time. Example data for this form is available <a href="../predictiontutorial.gsp#exampledata">here</a>. You may also use the button below to insert sample data. Please note that you will always need to enter the verification string at the bottom of the page, yourself, in order to submit a job!</p>
                           <p><b>Current problem:</b> Regrettably, our server is currently connected to the internet via a rather unreliable connection. This may cause connection timeouts (caused by server side) when uploading big files. Please use the web link upload option, instead, if you experience such problems. We apologize for the inconvenience!</p>
                           <g:actionSubmit action="fillSample" value="Fill in Sample Data" />
                           <p>We recommend that you specify an <b>E-mail address</b>.</p>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="email_adress">E-mail</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'email_adress','errors')}">
                                       <input type="text" id="email_adress" name="email_adress" value="${fieldValue(bean:predictionInstance,field:'email_adress')}"/> &nbsp;
                                       <g:checkBox name="agree_email" value="${predictionInstance?.agree_email}" />
                                       &nbsp;If I provide an e-mail address, I agree that it will be stored on the server until the computations of my job have finished. I agree to receive e-mails that are related to the particular AUGUSTUS job that I submitted. <a href="../help.gsp#email"><small>Help</small></a>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           You must <b>either</b> upload a *.tar.gz archive with AUGUSTUS species parameters from your computer <b>or</b> specify a project identifier: &nbsp; <a href="../help.gsp#which_files_pred"><small>Help</small></a>
                           <br>
                           <br>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="ArchiveFile"><b>AUGUSTUS species parameters</b> <font color="#FF0000">*</font></label>
                                    </td>
                                    <td valitn="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valitn="top">Upload an archive file  <font size="1">(max. 100 MB)</font>: &nbsp; <a href="../help.gsp#archive"><small>Help</small></a>
                                    </td>
                                    <td valitn="top">
                                       <g:if test="${predictionInstance.has_param_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="ArchiveFile" name="ArchiveFile"/></label>
                                       <g:if test="${predictionInstance.has_param_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="project_id">specify a project identifier:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'project_id','errors')}">
                                       <input type="text" id="project_id" name="project_id" value="${fieldValue(bean:predictionInstance,field:'project_id')}"/> <a href="../help.gsp#project_id"><small>Help</small></a>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="project_id">select an organism:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'species_select','errors')}">
                                       <g:if test="${predictionInstance.has_select == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <g:select name="species_select" from="${[
                                          'Acyrthosiphon pisum (animal)', 
                                          'Aedes aegypti (animal)', 
                                          'Amphimedon queenslandica (animal)', 
                                          'Apis mellifera (animal)', 
                                          'Bombus terrestris (animal)',
                                          'Brugia malayi (animal)', 
                                          'Caenorhabditis elegans (animal)', 
                                          'Callorhinchus milii (animal)', 
                                          'Camponotus floridanus (animal)',
                                          'Danio rerio (animal)',
                                          'Drosophila melanogaster (animal)', 
                                          'Gallus gallus domesticus (animal)',
                                          'Heliconius melpomene (animal)',
                                          'Homo sapiens (animal)', 
                                          'Nasonia vitripennis (animal)', 
                                          'Petromyzon marinus (animal)',
                                          'Rhodnius prolixus (animal)', 
                                          'Schistosoma mansoni (animal)', 
                                          'Tribolium castaneum (animal)', 
                                          'Trichinella spiralis (animal)', 
                                          'Tetrahymena thermophila (alveolata)',
                                          'Toxoplasma gondii (alveolata)',
                                          'Leishmania tarantolae (protozoa)',
                                          'Arabidopsis thaliana (plant)',
                                          'Chlamydomonas reinhardtii (alga)',
                                          'Galdieria sulphuraria (alga)',
                                          'Solaneum lycopersicum (plant)',
                                          'Triticum/wheat (plant)',
                                          'Zea mays (plant)',
                                          'Aspergillus fumigatus (fungus)',
                                          'Aspergillus nidulans (fungus)',
                                          'Aspergillus oryzae (fungus)',
                                          'Aspergillus terreus (fungus)',
                                          'Botrytis cinerea (fungus)',
                                          'Candida albicans (fungus)',
                                          'Candida guilliermondii (fungus)',
                                          'Candida tropicalis (fungus)',
                                          'Chaetomium globosum (fungus)',
                                          'Coccidioides immitis (fungus)',
                                          'Conidiobolus coronatus (fungus)',
                                          'Coprinus cinereus (fungus)',
                                          'Cryptococcus neoformans (fungus)',
                                          'Debarymomyces hansenii (fungus)',
                                          'Encephalitozoon cuniculi (fungus)',
                                          'Eremothecium gossypii (fungus)',
                                          'Fusarium graminearum (fungus)',
                                          'Histoplasma capsulatum (fungus)',
                                          'Kluyveromyces lactis (fungus)',
                                          'Laccaria bicolor (fungus)',
                                          'Lodderomyces elongisporus (fungus)',
                                          'Magnaporthe grisea (fungus)',
                                          'Neurospora crassa (fungus)',
                                          'Phanerochaete chrysosporium (fungus)',
                                          'Pichia stipitis (fungus)',
                                          'Rhizopus oryzae (fungus)',
                                          'Saccharomyces cerevisiae (fungus)',
                                          'Schizosaccharomyces pombe (fungus)',
                                          'Ustilago maydis (fungus)',
                                          'Verticillium longisporum (fungus)',
                                          'Yarrowia lipolytica (fungus)',
                                          'Sulfolobus solfataricus (archaeon)',
                                          'Escherichia coli (bacterium)',
                                          'Thermoanaerobacter tengcongensis (bacterium)'
                                          ]}" 
                                          value="${fieldValue(bean:predictionInstance,field:'species_select')}" noSelection="${['null':'Select One...']}"/>
                                       <a href="../help.gsp#project_id"><small>Help</small></a>
                                       <g:if test="${predictionInstance.has_select == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           You must <b>either</b> upload a genome file from your computer <b>or</b> specify a web link to a genome file: &nbsp; <a href="../help.gsp#upload_link"><small>Help</small></a>
                           <br>
                           <br>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="GenomeFile"><b>Genome file</b> <font color="#FF0000">*</font>&nbsp; (max. 250000 scaffolds) <a href="../help.gsp#genome_file"><small>Help</small></a></label>
                                    </td>
                                    <td valign="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 100 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${predictionInstance.has_genome_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="GenomeFile" name="GenomeFile"/>
                                       <g:if test="${predictionInstance.has_genome_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="genome_ftp_link">specify web link to genome file <font size="1">(max. 1 GB)</font>:</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'genome_ftp_link','errors')}">
                                       <input type="text" id="genome_ftp_link" name="genome_ftp_link" value="${fieldValue(bean:predictionInstance,field:'genome_ftp_link')}"/>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           You may (optionally) also specify one or several of the following files that contain external evidence for protein coding genes: <a href="../help.gsp#which_files_pred"><small>Help</small></a><br><br>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="EstFile"><b>cDNA file</b> &nbsp; <small><b><i>Non-commercial users only</i></b></small> &nbsp;<a href="../help.gsp#cDNA"><small>Help</small></a></label>
                                    </td>
                                    <td valign="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 100 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${predictionInstance.has_est_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="EstFile" name="EstFile"/>
                                       <g:if test="${predictionInstance.has_est_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="est_ftp_link">specify web link to cDNA file <font size="1">(max. 1 GB)</font>:</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'est_ftp_link','errors')}">
                                       <input type="text" id="est_ftp_link" name="est_ftp_link" value="${fieldValue(bean:predictionInstance,field:'est_ftp_link')}"/>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td><br></td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="hint_file"><b>Hints file</b> &nbsp; <a href="../help.gsp#hints"><small>Help</small></a></label>
                                    </td>
                                    <td valign="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 200 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${predictionInstance.has_hint_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="HintFile" name="HintFile"/>
                                       <g:if test="${predictionInstance.has_hint_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           <p>The following checkboxes allow you to modify the gene prediction behavior of AUGUSTUS:</p>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="UtrPrediction"><b>UTR prediction</b> &nbsp; <a href="../help.gsp#utr"><small>Help</small></a></label>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:if test="${predictionInstance.has_utr == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <g:checkBox name="utr" value="${false}" /> predict UTRs (requires species-specific UTR parameters)
                                       <g:if test="${predictionInstance.has_utr == true}"></div></g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <br><label for="StrandPrediction"><b>Report genes on</b></label>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:if test="${predictionInstance.has_strand == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <g:radio name="pred_strand" value="1" checked="true"/> both strands  &nbsp; &nbsp; <g:radio name="pred_strand" value="2"/> forward strand only &nbsp; &nbsp; <g:radio name="pred_strand" value="3"/> reverse strand only
                                       <g:if test="${predictionInstance.has_strand == true}"></div></g:if>
                                    </td>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <br><label for="AlternativeTranscripts"><b>Alternative transcripts:</b></label>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:radio name="alt_transcripts" value="1" checked="true"/>
                                       none  &nbsp; &nbsp; 
                                       <g:radio name="alt_transcripts" value="2"/>
                                       few &nbsp; &nbsp; 
                                       <g:radio name="alt_transcripts" value="3"/>
                                       medium &nbsp; &nbsp; 
                                       <g:radio name="alt_transcripts" value="4"/>
                                       many 
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <br><label for="GeneStructure"><b>Allowed gene structure:</b>&nbsp; <a href="../help.gsp#allowedGeneStructure"><small>Help</small></a></label>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:if test="${predictionInstance.has_structures == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <g:radio name="allowed_structures" value="1" checked="true"/> predict any number of (possibly partial) genes<br>
                                       <g:radio name="allowed_structures" value="2"/> only predict complete genes<br>
                                       <g:radio name="allowed_structures" value="3"/> only predict complete genes - at least one<br>
                                       <g:radio name="allowed_structures" value="4"/>  predict exactly one gene<br><br>
                                       <g:if test="${predictionInstance.has_structures == true}">
                                       </div></g:if>
                                       <g:if test="${predictionInstance.has_conflicts == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <g:checkBox name="ignore_conflicts" value="${false}" /> ignore conflicts with other strand
                                       <g:if test="${predictionInstance.has_conflicts == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           <p>
                              &nbsp; 
                              <g:checkBox name="agree_nonhuman" value="${predictionInstance?.agree_nonhuman}" />
                              &nbsp;<b>I am not submitting personalized human sequence data (mandatory).</b>
                           </p>
                           <p>We use a <b>verification string</b> to figure out whether you are a <b>human</b> person. Please type the text in the image below into the text field next to the image.
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:if test="${predictionInstance.warn == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <img src="${createLink(controller: 'simpleCaptcha', action: 'captcha')}"/> &nbsp; &nbsp; <g:textField name="captcha"/> <font color="#FF0000">*</font>
                                       <g:if test="${predictionInstance.warn == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <p><font color="#FF0000">*</font>) mandatory input arguments</p>
                        </div>
                        <div class="buttons" onclick="toggle_visibility('spinner');">
                           <span class="button"><input class="commit" type="submit" value="Start Predicting" /></span><br>
                        </div>
                     </g:uploadForm>
                     <br>
                     <table>
                        <tr>
                           <td>
                              <div id="spinner" style='display:none;' align="center">
                                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="http://www.animateit.net/data/media/june2010/compete.gif" border="0" alt="Spinner" /><br>
                                 <p>We are processing your request... please do not close this window and do not click on the submission button again, until this message disappears!</p>
                              </div>
                           </td>
                           <td><img src="../images/spacer.jpg" alt="spacer image"></td>
                        </tr>
                     </table>
                     <p>&nbsp;</p>
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

