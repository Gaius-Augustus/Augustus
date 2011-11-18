<html>
    <head>
        <META HTTP-EQUIV="Refresh" CONTENT="20">
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Current Prediction Jobs</title>
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
        Bioinformatics Web Server
      </div>
      <div id="bannertitel2">
        AUGUSTUS Training
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
    Current Training Jobs
</span>

  <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Kopfbereich *********************************************// -->

<!-- ***** Start: Koerper ************************************************// -->
<div id="koerper">

  <div id="linke_spalte">
     <ul class="menu">
         <li><a href="../index.gsp"><span>Introduction</span></a></li>
         <li><g:link controller="training" action="create"><span>Submitt Training</span></g:link></li>
         <li><g:link controller="prediction" action="create"><span>Submitt Prediction</span></g:link></li>
         <li><g:link controller="help" action="list"><span>Help</span></g:link></li>
         <li><a href="../references.gsp"><span>Links & References</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de"><span>Bioinformatics Group</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/impressum.html"><span>Impressum</span></a></li>
     </ul>
  </div>

 <div id="mittel_spalte">
<div class="main" id="main">
   <h1><g:link controller="prediction" action="list">Current AUGUSTUS Prediction Jobs</g:link></h1>
            <g:if test="${flash.message}">
            <div class="message">${flash.message}</div>
            </g:if>
            <div class="list">
                <table frame="border" rules="all" cellpadding="5">
                    <thead>
                        <tr>
                                <g:sortableColumn property="id" title="&nbsp;ID&nbsp;" />
                                <g:sortableColumn
                                property="accession_id"
                                title="&nbsp;Project ID&nbsp;" />

                                <g:sortableColumn
                                property="dateCreated"
                                title="&nbsp;Date&nbsp;" />

                                <g:sortableColumn
                                property="job_status"
                                title="&nbsp;Job Status&nbsp;" />

                        </tr>
                    </thead>
                    <tbody>
                    <g:each in="${predictionInstanceList}" status="i" var="predictionInstance">
                        <tr class="${(i % 2) == 0 ? 'odd' : 'even'}">

                            <td><g:link action="show" id="${trainingInstance.id}">${fieldValue(bean:predictionInstance, field:'id')}</g:link></td>

                            <td>${fieldValue(bean:predictionInstance, field:'accession_id')}</td>

                            <td>${fieldValue(bean:predictionInstance, field:'dateCreated')}</td>

                            <td>${fieldValue(bean:predictionInstance, field:'job_status')}</td>

                        </tr>
                    </g:each>
                    </tbody>
                </table>
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
      <a href="mailto:bioinformatik.greifswald@gmail.com" title="E-Mail bioinformatik.greifswald@gmail.com, opens the standard mail program">bioinformatik.greifswald@gmail.com</a>
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
