<accuhtml xmlns="http://www.w3.org/1999/xhtml" xml:lang="de-de" lang="de-de" dir="ltr" >
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
	<meta name="layout" content="main" />
        <title>Create Training</title>
	<script type="text/javascript" src="js/md_stylechanger.js"></script>
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
       <a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Institut f&uuml;r Mathematik und Informatik"><img src="images/header.gif" alt="Directly to home" /> </a>
   </div>
   <div id="banner_mitte">
      <div id="bannertitel1">
        Bioinformatics Web Server at University of Greifswald
      </div>
      <div id="bannertitel2">
        Gene Prediction with AUGUSTUS
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
    Introduction</span>

  <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Kopfbereich *********************************************// -->

<!-- ***** Start: Koerper ************************************************// -->
<div id="koerper">

  <div id="linke_spalte">
     <ul class="menu">
         <li><div id="linksMenuText">AUGUSTUS Web Server Navigation</div></li>
         <li id="current"><a href="index.gsp"><span>Introduction</span></a></li>
         <li><a href="about.gsp"><span>About AUGUSTUS</span></a></li>
         <li><a href="accuracy.gsp"><span>Accuracy</span></a></li>
         <li><a href="trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
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
	 <li><a href="http://bioinf.uni-greifswald.de/bioinf/forum">AUGUSTUS Forum</a></li>
         <li><a href="http://bioinf.uni-greifswald.de/augustus/downloads/index.php">Download AUGUSTUS</a></li>
	 <li><a href="http://bioinf.uni-greifswald.de/augustus">Old AUGUSTUS gene prediction web server</a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/braker">BRAKER</a></li>
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
	  <a href="index.gsp" class="contentpagetitle">Welcome to the WebAUGUSTUS</a>
        </td>
      </tr>
    </table>
    <p> AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences. This web server provides an interface for training AUGUSTUS for predicting genes in genomes of novel species. It also enables you to predict genes in a genome sequence with already trained parameters.</p>
    <p>AUGUSTUS usually belongs to the most accurate programs for the species it is trained for. Often it is the most accurate ab initio program. For example, at the independent gene finder assessment (EGASP) on the human ENCODE regions AUGUSTUS was the most accurate gene finder among the tested ab initio programs. At the more recent nGASP (worm), it was among the best in the ab initio and transcript-based categories. See <a href="accuracy.gsp">accuracy statistics</a> for further details.</p>
    <p>Please be aware that gene prediction accuracy of AUGUSTUS always depends on the quality of the training gene set that was used for training species specific parameters. You should not expect the greatest accuracy from fully automated training gene generation as provided by this web server application. Instead, you should manually inspect (and maybe interatively improve) the training gene set.</p>
<p>AUGUSTUS is already trained for a number of genomes and you find the according parameter sets at <a href="predictiontutorial.gsp#param_id">the prediction tutorial</a>. <b>Please check whether AUGUSTUS was already trained for your species before submitting a new training job.</b></p>
    <p><a href="http://bioinf.uni-greifswald.de/augustus/">The Old AUGUSTUS web server</a> offers similar gene prediction services but no parameter training service.</p>
<br><br>
<h2>OK, I got it! Take me straight to...</h2>
<p>
<ul>
  <li><g:link controller="training" action="create"><span>AUGUSTUS training submission</span></g:link></li>
  <li><g:link controller="prediction" action="create"><span>AUGUSTUS prediction submission</span></g:link></li>
</ul>
</p>
  </div>

  <div id="rechte_spalte">
    <div id="schriftgroesse">
      <script type="text/javascript">
	//<![CDATA[
                  document.write('<a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Increase size" onclick="changeFontSize(2); return false;" class="larger">bigger</a><span class="unseen">&nbsp;</span>');
                  document.write('<a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Decrease size" onclick="changeFontSize(-2); return false;" class="smaller">smaller</a><span class="unseen">&nbsp;</span>');
                  document.write('<a href="http://www.math-inf.uni-greifswald.de/mathe/index.php" title="Revert styles to default" onclick="revertStyles(); return false;" class="reset">reset</a></p>');
        //]]>
      </script>
    </div>
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
