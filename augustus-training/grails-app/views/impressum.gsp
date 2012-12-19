<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="de-de" lang="de-de" dir="ltr" >
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
    Impressum</span>

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
         <li><a href="help.gsp"><span>Help</span></a></li>
         <li><a href="datasets.gsp"><span>Datasets for Download</span></a></li>
         <li><a href="predictions_for_download.gsp"><span>Predictions for Download</span></a></li>
         <li><a href="references.gsp"><span>Links & References</span></a></li>
         <li id="current"><a href="impressum.gsp"><span>Impressum</span></a></li>
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
	  <a href="impressum.gsp" class="contentpagetitle">Impressum</a>
        </td>
      </tr>
    </table>

	<h1>Publisher</h1>
	<p>University of Greifswald</p>
	<p>Institute for Mathematics and Computer Science</p>
	<p>Bioinformatics Group</p>
	<p>Represented by</p>
	<p>Prof. Dr. Mario Stanke</p>
	<p>Walther-Rathenau-Stra&szlig;e 47</p>
	<p>17487 Greifswald</p>
	<p><b>Tel.:</b> +49 (0)3834 86- 46 24</p>
	<p><b>Fax:</b> +49 (0)3834 86- 46 40</p>
	<p><b>E-Mail:</b> <a href="mailto:mario.stanke@uni-greifswald.de">mario.stanke@uni-greifswald.de</a></p>

	<h2>Editorial Responsibility according to &sect; 6 MDStV</h2>
	<p>The group leader is responsible for contents of this web page.</p>
	
	<h2>Concept, Design, Technical Implementation</h2>
	
	<p>Dr. Katharina Hoff, Marco Frehse and Sebastian Matuszewski</p>

	<h2>Images</h2>

	<p>University archive, Marco Frehse, Sebastian Matuszewski, Dr. Katharina Hoff, Prof. Dr. Mario Stanke</p>
	<h2>Webmaster</h2>
	<p>Dr. Katharina Hoff</p>
	<h2>Content</h2>
	<p>A web server application for training and running the eukaryotic gene prediction software AUGUSTUS</p>
	<h2>First online</h2>
	<p>December 2011</p>

	<h3>1. Ambit</h3>
	<p>This impressum applies to pages with the URL http://bioinf.uni-greifswald.de/trainaugustus and http://bioinf.uni-greifswald.de/webaugustus. </p>

	<h3>2. Content of these web pages</h3>
	<p>The University of Greifswald does not give any warranty for timeliness, correctness, completeness, runnability or quality of the provided information and services. Liability claims against the university that refer to damages of material or ideal value and that resulted from the usage or non-usage of error-containing and incomplete information, are generally excluded if no verifiable deliberate or coarsely careless behaviour by the university exists.</p>
	<p>All services are free and without committment. The university explicitely reserves the right to change, complete or delete parts of web pages or to stop the entire web service, for intervals of time or ultimately, without separate announcement.</p>
	<p>Please send hints for corrections to <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a></p>

	<h3>3. References and Links</h3>
	<p>The following section is specific to German law and is therefore given in German. It basically expresses that we are not responsible for the content of web pages that we link to from our web pages.</p>

<p><i>Bei direkten oder indirekten Verweisen auf fremde Webseiten ("Hyperlinks"), die au&szlig;erhalb des Verantwortungsbereiches der Universit&auml;t liegen, w&uuml;rde eine Haftungsverpflichtung ausschlie&szlig;lich in dem Fall in Kraft treten, in dem der/die Autor/in von den Inhalten Kenntnis hat und es ihm/ihr technisch m&ouml;glich und zumutbar w&auml;re, die Nutzung im Falle rechtswidriger Inhalte zu verhindern.</i></p>

<p><i>Der Universit&auml;t erkl&auml;rt hiermit ausdr&uuml;cklich, dass zum Zeitpunkt der Linksetzung keine illegalen Inhalte auf den zu verlinkenden Seiten erkennbar waren. Auf die aktuelle und zuk&uuml;nftige Gestaltung, die Inhalte oder die Urheberschaft der verlinkten/verkn&uuml;pften Seiten hat die Universit&auml;t keinerlei Einfluss. Deshalb distanziert sie sich hiermit ausdr&uuml;cklich von allen Inhalten aller verlinkten/verkn&uuml;pften Seiten, die nach der Linksetzung ver&auml;ndert wurden. Diese Feststellung gilt f&uuml;r alle innerhalb des eigenen Internetangebotes gesetzten Links und Verweise sowie f&uuml;r Fremdeintr&auml;ge in von der Universit&auml;t eingerichteten G&auml;steb&uuml;chern, Diskussionsforen, Linkverzeichnissen, Mailinglisten und in allen anderen Formen von Datenbanken, auf deren Inhalt externe Schreibzugriffe m&ouml;glich sind. F&uuml;r illegale, fehlerhafte oder unvollst&auml;ndige Inhalte und insbesondere f&uuml;r Sch&auml;den, die aus der Nutzung oder Nichtnutzung solcherart dargebotener Informationen entstehen, haftet allein der Anbieter der Seite, auf welche verwiesen wurde, nicht derjenige, der &uuml;ber Links auf die jeweilige Ver&ouml;ffentlichung lediglich verweist.</i></p>

<h3>4. Urheber- und Kennzeichenrecht</h3>
<p>The university tries to respect the copyright of all used graphics, fotos, audo files, video sequences, and texts, and to to use self-created graphics, fotos, audio files, video sequences and texts where possible.</p>

<p>All brand- and trademarks that are mentioned within the web service that is labeled by this impressum and possibly brand- and trademarks that are proprietary by thirs parties unlimitedly underlie the respectively valid marking law and the right of owernship of the respectively registered owner. Once cannot conclude that a brand- or trademark is unprotected by the rights of third parties only because the brand- or trademark is mentioned!</p>

<p>The copyright for online services by the university that are labeled with the above mentioned impressum solely stays with the university. Duplication or the usage of graphics, fotos, audio files, video sequences or texts from these web pages in electronic or printed publications is not permitted without explicite approval by the university.</p>

<p>This does not affect the "copying" of files onto personal computers for viewing www-pages in a browser. It also does not affect press releases of the university and software that is under free software licenses. Those are free for usage and modification without special permission. Please notify the webmaster of  improper usage.</p>

<h1>Law Binding Impressum</h1>

<p>This impressum web site is only a rough translation. You find the law binding impressum in German language at <a href="http://www.math-inf.uni-greifswald.de/mathe/index.php/impressum">http://www.math-inf.uni-greifswald.de/mathe/index.php/impressum</a>.</p>
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
