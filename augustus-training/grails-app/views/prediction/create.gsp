

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Create Prediction</title>
    <script type="text/javascript" src="js/md_stylechanger.js"></script>
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
        <a href="http://www.uni-greifswald.de" target="_blank" class="mainleveltop_" >University of Greifswald</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.mnf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Faculty</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.math-inf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Institute</a><span class="mainleveltop_">&nbsp;|&nbsp;</span><a href="http://bioinf.uni-greifswald.de/" target="_blank" class="mainleveltop_">Bioinformatics Group</a>
      </td>
    </tr>
  </table>
</div>
<div id="banner">
   <div id="banner_links">
       <a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Institut f&uuml;r Mathematik und Informatik"><img src="../images/header.gif" alt="Directly to home" /> </a>
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
     <img src="../images/spirale.gif" align="left" />
     </a>
   </div>
</div>

<div id="wegweiser">
  Navigation for: &nbsp; &nbsp;<span class="breadcrumbs pathway">
    Submit Prediction</span>

  <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Kopfbereich *********************************************// -->
<!-- ***** Start: Koerper ************************************************// -->
<div id="koerper">

  <div id="linke_spalte">
     <ul class="menu">
         <li><div id="linksMenuText">AUGUSTUS Web Server Navigation</div></li>
         <li><a href="../index.gsp"><span>Introduction</span></a></li>
         <li><a href="../about.gsp"><span>About AUGUSTUS</span></a></li>
         <li><a href="../accuracy.gsp"><span>Accuracy</span></a></li>
         <li><a href="../trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
         <li><g:link controller="training" action="create"><span>Submit Training</span></g:link></li>
         <li><a href="../predictiontutorial.gsp"><span>Prediction Tutorial</span></a></li>
         <li id="current"><g:link controller="prediction" action="create"><span>Submit Prediction</span></g:link></li>
         <li><a href="../help.gsp"><span>Help</span></a></li>
         <li><a href="../datasets.gsp"><span>Datasets for Download</span></a></li>
         <li><a href="../predictions_for_download.gsp"><span>Predictions for Download</span></a></li>
         <li><a href="../references.gsp"><span>Links & References</span></a></li>
         <li><a href="../impressum.gsp"><span>Impressum</span></a></li>
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
            <div class="main" id="main">
			<noscript><p><b><span style="color:red">Please enable javascript in your browser in order to display the submission form correctly!</span></b> Form functionality is not affected significantly while javascript is disabled, but it looks less pretty.</p>
			</noscript>
            <g:uploadForm action="commit" method="post" name="submissionform">
            <fieldset><legend>
                <table class="contentpaneopen">
                  <tr>
                    <td class="contentheading" width="100%">
                      
                        <g:link controller="prediction" action="create">Data Input for Running AUGUSTUS</g:link>
                    </td>
                  </tr>
                </table>
              </legend><p>
                <div class="dialog">
		    <p>Use this form to submit your data for running AUGUSTUS on new genomic data with already available pre-trained parameters.</p>

<p>Please read the <a href="../predictiontutorial.gsp">prediction tutorial</a> before submitting a job for the first time. Example data for this form is available <a href="../predictiontutorial.gsp#exampledata">here</a>. You may also use the button below to insert sample data. Please note that you will always need to enter the verification string at the bottom of the page, yourself, in order to submit a job!</p>

           <g:actionSubmit action="fillSample" value="Fill in Sample Data" />



		    <p>We strongly recommend that you specify an <b>E-mail address</b>! Please read the <a href="../help.gsp#email"><small>Help</small></a> page before submitting a job without e-mail address!</p>
                    <table>
                        <tbody>
                                <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="email_adress">E-mail</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'email_adress','errors')}">
                                    <input type="text" id="email_adress" name="email_adress" value="${fieldValue(bean:predictionInstance,field:'email_adress')}"/> &nbsp; <a href="../help.gsp#email"><small>Help</small></a>
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
                                    <input type="file" id="ArchiveFile" name="ArchiveFile"/></label>
                                </td>
                            </tr>
                            <tr class="prop">
                                <td>&nbsp;<b>or</b>&nbsp;</td><td></td>
                            </tr>
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="project_id">specify a project identifier:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'project_id','errors')}">
                                    <input type="text" id="project_id" name="project_id" value="${fieldValue(bean:predictionInstance,field:'project_id')}"/> <a href="../help.gsp#project_id"><small>Help</small></a>

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
                                    <label for="GenomeFile"><b>Genome file</b> <font color="#FF0000">*</font>&nbsp; <a href="../help.gsp#genome_file"><small>Help</small></a></label>
                                </td>
                                <td valign="top">
                                </td>
                            </tr>
                            <tr class="prop">
                                <td valitn="top">Upload a file <font size="1">(max. 100 MB)</font>:</td>
                                <td valitn="top">
                                    <input type="file" id="GenomeFile" name="GenomeFile"/>
                                </td>
                            </tr>
                            <tr class="prop">
                                <td>&nbsp;<b>or</b>&nbsp;</td><td></td>
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
                                            <input type="file" id="EstFile" name="EstFile"/>
                              </td>
                            </tr>
                            <tr class="prop">
                              <td>&nbsp;<b>or</b>&nbsp;</td><td></td>
                            </tr>
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="est_ftp_link">specify web link to cDNA file <font size="1">(max. 1 GB)</font>:</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:predictionInstance,field:'est_ftp_link','errors')}">
                                    <input type="text" id="est_ftp_link" name="est_ftp_link" value="${fieldValue(bean:predictionInstance,field:'est_ftp_link')}"/>
                                </td>
                            </tr>
                            <tr class="prop"><td><br></td><td></td></tr>
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
                                    <input type="file" id="HintFile" name="HintFile"/>
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
							<g:checkBox name="utr" value="${false}" /> predict UTRs (requires species-specific UTR parameters)
						</td>
					</tr>
					<tr class="prop">
						<td valign="top" class="name">
							<br><label for="StrandPrediction"><b>Report genes on</b></label>
						</td>
					</tr>
					<tr class="prop">
						<td valign="top" class="name">
							<g:radio name="pred_strand" value="1" checked="true"/> both strands  &nbsp; &nbsp; <g:radio name="pred_strand" value="2"/> forward strand only &nbsp; &nbsp; <g:radio name="pred_strand" value="3"/> reverse strand only
						</td>
					<tr class="prop">
						<td valign="top" class="name">
							<br><label for="AlternativeTranscripts"><b>Alternative transcripts:</b></label>
						</td>
					</tr>
					<tr class="prop">
						<td valign="top" class="name">
							<g:radio name="alt_transcripts" value="1" checked="true"/> none  &nbsp; &nbsp; <g:radio name="alt_transcripts" value="2"/> few &nbsp; &nbsp; <g:radio name="alt_transcripts" value="3"/> medium &nbsp; &nbsp; <g:radio name="alt_transcripts" value="4"/> many 
						</td>
					</tr>
					<tr class="prop">
						<td valign="top" class="name">
							<br><label for="GeneStructure"><b>Allowed gene structure:</b>&nbsp; <a href="../help.gsp#allowedGeneStructure"><small>Help</small></a></label>
						</td>
					</tr>
					<tr class="prop">
						<td valign="top" class="name">
							<g:radio name="allowed_structures" value="1" checked="true"/> predict any number of (possibly partial) genes<br>
							<g:radio name="allowed_structures" value="2"/> only predict complete genes<br>
							<g:radio name="allowed_structures" value="3"/> only predict complete genes - at least one<br>
							<g:radio name="allowed_structures" value="4"/>  predict exactly one gene<br><br>
							<g:checkBox name="ignore_conflicts" value="${false}" /> ignore conflicts with other strand
						</td>
					</tr>
				</tbody>
			</table>
                    <br>
		    <p>We use a <b>verification string</b> to figure out whether you are a <b>human</b> person. Please type the text in the image below into the text field next to the image.
			<table>
				<tbody>
					<tr class="prop">
						<td valign="top" class="name">
							<img src="${createLink(controller: 'simpleCaptcha', action: 'captcha')}"/> &nbsp; &nbsp; <g:textField name="captcha"/> <font color="#FF0000">*</font>
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
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="http://www.animateit.net/data/media/june2010/compete.gif" border="0" alt="Spinner" /><br><p>We are processing your request... please do not close this window and do not click on the submission button again, until this message disappears!</p>
					</div>
				</td>
				<td><img src="../images/spacer.jpg" alt="spacer image"></td>
			</tr>
		</table>
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
      <a href="mailto:augustus-web@uni-greifswald.de" title="E-Mail augustus-web@uni-greifswald.de, opens the  standard mail program">augustus-web@uni-greifswald.de</a>
    </div>
    </div>

    <div class="beendeFluss"></div>
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
     <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="../images/top.gif" />
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
