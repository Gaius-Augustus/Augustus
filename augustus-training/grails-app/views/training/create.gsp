

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
               <g:hasErrors bean="${trainingInstance}">
                  <div class="errors">
                     <g:renderErrors bean="${trainingInstance}" as="list" />
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
                        <h1 class="csc-firstHeader">Data Input for Training AUGUSTUS</h1>
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
                     <g:uploadForm action="commit" method="post" >
                        <fieldset>
                        <p>
                        <div class="dialog">
                           <p>Use this form to submit data for training AUGUSTUS parameters for novel species/new genomic data.</p>
                           <p><b>Before submitting a training job</b> for your species of interest, please check whether parameters have already been trained and have been made publicly available for your species at <a href="../predictiontutorial.gsp#param_id">our species overview table</a></p>
                           <p>Please read the <a href="../trainingtutorial.gsp">training tutorial</a> before submitting a job for the first time. Example data for this form is available <a href="../trainingtutorial.gsp#exampledata">here</a>. You may also use the button below to insert sample data. Please note that you will always need to enter the verification string at the bottom of the page, yourself, in order to submit a job!</p>
                           <g:actionSubmit action="fillSample" value="Fill in Sample Data" />
                           <p>We strongly recommend that you specify an <b>E-mail address</b>! Please read the <a href="../help.gsp#email"><small>Help</small></a> page before submitting a job without e-mail address! You have to give a <b>species name</b>, and a <b>genome file</b>!</p>
                           <p><b>Current problem:</b> Regrettably, our server is currently connected to the internet via a rather unreliable connection. This may cause connection timeouts (caused by server side) when uploading big files. Please use the web link upload option, instead, if you experience such problems. We apologize for the inconvenience!</p>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="email_adress">E-mail</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'email_adress','errors')}">
                                       <input type="text" id="email_adress" name="email_adress" value="${fieldValue(bean:trainingInstance,field:'email_adress')}"/> 
                                       &nbsp; <a href="../help.gsp#email"><small>Help</small></a><br>
                                       <!--<g:if test="${trainingInstance.agree_email == true}"><div class="prop_warn"></g:if>-->
                                       &nbsp; 
                                       <g:checkBox name="agree_email" value="${trainingInstance?.agree_email}" />
                                       &nbsp;If I provide an e-mail address, I agree that it will be stored on the server until the computations of my job have finished. I agree to receive e-mails that are related to the particular AUGUSTUS job that I submitted.
                                       <!--<g:if test="${trainingInstance.agree_email == true}"></div></g:if>-->
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="project_name">Species name <font color="#FF0000">*</font></label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'project_name','errors')}">
                                       <input type="text" maxlength="30" id="project_name" name="project_name" value="${fieldValue(bean:trainingInstance,field:'project_name')}"/> &nbsp; <a href="../help.gsp#species_name"><small>Help</small></a>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <p>There are two options for sequence file (fasta format) transfer:<br>You may <b>either</b> upload data files from your computer <b>or</b> specify web links. &nbsp; <a href="../help.gsp#upload_link"><small>Help</small></a><br><br><font color="#FF0000">Please read our <a href="../help.gsp#most_common_problem">instructions about fasta headers</a> before using this web service!</font> Most problems with this web service are caused by a wrong fasta header format!</p>
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
                                       <g:if test="${trainingInstance.has_genome_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="GenomeFile" name="GenomeFile"/>
                                       <g:if test="${trainingInstance.has_genome_file == true}">
                                          </div>
                                       </g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td><td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="genome_ftp_link">specify web link to genome file <font size="1">(max. 1 GB)</font>:</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'genome_ftp_link','errors')}">
                                       <input type="text" id="genome_ftp_link" name="genome_ftp_link" value="${fieldValue(bean:trainingInstance,field:'genome_ftp_link')}"/> 
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           You need to specify <b>at least one</b> of the following files: <font color="#FF0000">*</font> <a href="../help.gsp#which_fiiles"><small>Help</small></a><br><br>
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="EstFile"><b>cDNA file</b> &nbsp; <small><b><i>Non-commercial users only</i></b></small> &nbsp; <a href="../help.gsp#cDNA"><small>Help</small></a></label>
                                    </td>
                                    <td valign="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 100 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${trainingInstance.has_est_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="EstFile" name="EstFile"/>
                                       <g:if test="${trainingInstance.has_est_file == true}"></div></g:if>
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
                                    <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'est_ftp_link','errors')}">
                                       <input type="text" id="est_ftp_link" name="est_ftp_link" value="${fieldValue(bean:trainingInstance,field:'est_ftp_link')}"/>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td><br></td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="ProteinFile"><b>Protein file</b> &nbsp; <small><b><i>Non-commercial users only</i></b></small> &nbsp; <a href="../help.gsp#protein"><small>Help</small></a></label>
                                    </td>
                                    <td valign="top">
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 100 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${trainingInstance.has_protein_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="ProteinFile" name="ProteinFile"/>
                                       <g:if test="${trainingInstance.has_protein_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td>&nbsp;<b>or</b>&nbsp;</td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="protein_ftp_link">specify web link to protein file <font size="1">(max. 1 GB)</font>:</label>
                                    </td>
                                    <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'protein_ftp_link','errors')}">
                                       <input type="text" id="protein_ftp_link" name="protein_ftp_link" value="${fieldValue(bean:trainingInstance,field:'protein_ftp_link')}"/>
                                    </td>
                                 </tr>
                                 <tr class="prop">
                                    <td><br></td>
                                    <td></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <label for="struct_file"><b>Training gene structure file</b> &nbsp; <a href="../help.gsp#structure"><small>Help</small></a> <font color="#FF0000">(gff or gb format, no gzip!)</font></label>
                                    </td>
                                    <td valign="top"></td>
                                 </tr>
                                 <tr class="prop">
                                    <td valign="top">Upload a file <font size="1">(max. 200 MB)</font>:</td>
                                    <td valign="top">
                                       <g:if test="${trainingInstance.has_struct_file == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <input type="file" id="StructFile" name="StructFile"/> 
                                       <g:if test="${trainingInstance.has_struct_file == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <br>
                           <!-- show some content upon click -->
                           <h2>Possible file combinations [<a title="show/hide" id="exp_file_options_link" href="javascript: void(0);" onclick="toggle(this, 'exp_file_options');"  style="text-decoration: none; color: #006699; ">click to minimize</a>]</h2>
                           <div id="exp_file_options" style="padding: 3px;">
                              <p>
                              <ul>
                                 <li><b>{genome file, cDNA file}</b><br>In this case, the cDNA file is used to create a training gene set. If cDNA quality is sufficient, also a UTR training set will be created.</li>
                                 <li><b>{genome file, protein file}</b><br>In this case, the protein file is used to create a training gene set. No UTR training set can be created.</li>
                                 <li><b>{genome file, gene structure file}</b><br>In this case, the gene structure file is used as a training gene set. If the gene structure file contains UTR elements, also a UTR training set will be created.</li>
                                 <li><b>{genome file, cDNA file, protein file}</b><br>In this case, the protein file will be used to create a training gene set. No UTR training set will be created. cDNA sequences will be used as evidence for prediction, only.</li>
                                 <li><b>{genome file, cDNA file, gene structure file}</b><br>In this case, the gene structure file is used as a training gene set. If the gene structure file contains UTR elements, also a UTR training set will be created. cDNA sequences will be used as evidence for prediction, only.</li>
                              </ul>
                              </p>
                              <h2>File combinations that are currently not supported</h2>
                              <p>
                              <ul>
                                 <li><b>{genome file, cDNA file, protein file, gene structure file}</b></li>
                                 <li><b>{genome file, protein file, gene structure file}</b></li>
                              </ul>
                              </p>
                           </div>
                           <script language="javascript">toggle(getObject('exp_file_options_link'), 'exp_file_options');</script>
                           <!-- end of javascript content on click -->
                           <br>
                           <p>
                              <!--<g:if test="${trainingInstance.agree_nonhuman == true}"><div class="prop_warn"></g:if>-->
                              &nbsp; 
                              <g:checkBox name="agree_nonhuman" value="${trainingInstance?.agree_nonhuman}" />
                              &nbsp;<b>I am not submitting personalized human sequence data (mandatory).</b> <a href="../help.gsp#nonhuman"><small>Help</small></a>
                           </p>
                           <!-- <g:if test="${trainingInstance.agree_nonhuman == true}"></div></g:if></p>-->
                           <p>We use a <b>verification string</b> to figure out whether you are a <b>human</b>. Please type the text in the image below into the text field next to the image.
                           <table>
                              <tbody>
                                 <tr class="prop">
                                    <td valign="top" class="name">
                                       <g:if test="${trainingInstance.warn == true}">
                                          <div class="prop_warn">
                                       </g:if>
                                       <img src="${createLink(controller: 'simpleCaptcha', action: 'captcha')}"/> &nbsp; &nbsp; <g:textField name="captcha"/><font color="#FF0000">*</font>
                                       <g:if test="${trainingInstance.warn == true}"></div></g:if>
                                    </td>
                                 </tr>
                              </tbody>
                           </table>
                           <p><font color="#FF0000">*</font>) mandatory input arguments</p>
                        </div>
                        <div class="buttons" onclick="toggle_visibility('spinner');">
                           <span class="button"><input class="commit" type="submit" value="Start Training" /></span>
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

