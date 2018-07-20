

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
                        <li class="navigation-sub__item"><a href="../../index.gsp">Introduction</a></li>
                        <li class="navigation-sub__item"><a href="../../about.gsp">About AUGUSTUS</a></li>
                        <li class="navigation-sub__item"><a href="../../accuracy.gsp">Accuracy</a></li>
                        <li class="navigation-sub__item"><a href="../../trainingtutorial.gsp">Training Tutorial</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="training" action="create">Submit Training</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="../../predictiontutorial.gsp">Prediction Tutorial</a></li>
                        <li class="navigation-sub__item">
                           <g:link controller="prediction" action="create">Submit Prediction</g:link>
                        </li>
                        <li class="navigation-sub__item"><a href="../../datasets.gsp">Datasets for Download</a></li>
                        <li class="navigation-sub__item"><a href="../../references.gsp">Links & References</a></li>
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
                        <h1 class="csc-firstHeader">Training AUGUSTUS<br>Job ${trainingInstance.accession_id}</h1>
                     </div>
                  </div>
                  <div id="c261665" class="csc-default">
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
                  </div>
                  <g:if test="${flash.message}">
                     <div class="message">${flash.message}</div>
                  </g:if>
                  <p>
                     <font color="#FF0000"><b>Please bookmark this page!</b></font> Training AUGUSTUS may take up to several weeks depending on the input data properties. Bookmarking this page ensures that you will be able to return to this page in order to find the results of your job, later.
                  </p>
                  <hr>
                  <h2><font color="#006699">Job Status</font></h2>
                  <p>
                     <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '0' || fieldValue(bean:trainingInstance, field:'job_status') == '1' || fieldValue(bean:trainingInstance, field:'job_status') == '2' || fieldValue(bean:trainingInstance, field:'job_status') == '3'}">
                        <g:if test = "${trainingInstance.old_url == null}">
                  <div style="width:600px;height:30px;border:1px solid #d2d2dc">
                  <p>
                  <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '0'|| fieldValue(bean:trainingInstance, field:'job_status') == '1'}">
                  <b><font color="#006699" size=2>Job submitted</font> <font color="#d2d2dc" size=2>&rarr; waiting for execution &rarr; computing &rarr; finished!</font></b><br>
                  </g:if>
                  <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '2'}">
                  <b><font color="#d2d2dc" size=2>Job submitted</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#006699" size=2>waiting for execution</font> <font color="#d2d2dc" size=2>&rarr; computing &rarr; finished!</font></b><br>
                  </g:if>
                  <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '3'}">
                  <b><font color="#d2d2dc" size=2>Job submitted</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#d2d2dc" size=2>waiting for execution</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#006699" size=2>computing</font> <font color="#d2d2dc" size=2>&rarr; finished!</font></b><br>
                  </g:if>
                  <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '4'}">
                  <b><font color="#d2d2dc" size=2>Job submitted</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#d2d2dc" size=2>waiting for execution</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#d2d2dc" size=2>computing</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#006699" size=2>finished!</font></b><br>
                  </g:if>
                  </p>
                  </div>
                  </g:if>
                  </g:if>
                  </p>
                  <g:if test ="${fieldValue(bean:trainingInstance, field:'job_status') == '4'}">
                     <p>
                        <b><font color="#d2d2dc" size=2>Job submitted</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#d2d2dc" size=2>waiting for execution</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#d2d2dc" size=2>computing</font> <font color="#ffb22a" size=2>&rarr;</font> <font color="#006699" size=2>finished!</font></b>
                     </p>
                  </g:if>
                  <g:if test ="${fieldValue(bean:trainingInstance, field:'job_error') == '5' || fieldValue(bean:trainingInstance, field:'job_status') == '5'}">
                     <p>
                        <b><font color="#f40707" size=2>An error occured when executing this job!</font></b>
                     </p>
                     <g:if test = "${trainingInstance.old_url != null}">
                        <p><b><font color="#FF0000">Data duplication!</font></b> A job with identical data was submitted before. You find the old job at <a href="${trainingInstance.old_url}">${trainingInstance.old_url}</a>.
                        </p>
                     </g:if>
                  </g:if>
                  <g:if test="${trainingInstance.job_status >= '4' && trainingInstance.results_urls != null}">
                     <hr>
                     <h2><font color="#006699">Results</font></h2>
                     ${trainingInstance.results_urls}
                     <p><b>Instructions</b></p>
                     <p>Please download the files listed above by clicking on the links.</p>
                     <p>All files and folders (except for the *.log and *.err file) are compressed. To unpack <tt>*.tar.gz</tt> archives, e.g. on linux type<br><br>
                        <tt>tar -xzvf *.tar.gz</tt><br><br>
                        For unpacking <tt>*.gz</tt> files, e.g. on linux type<br><br>
                        <tt>gunzip *.gz.</tt>
                     </p>
                     <p>Further instructions about results contents are given at the <a href="../../trainingtutorial.gsp">Training Tutorial</a> and the <a href="../../predictiontutorial.gsp">Prediction Tutorial</a>. Typical error messages from the log files are explained at <a href="../../help.gsp#noResults">Help</a>.
                     </p>
                  </g:if>
                  <h2><font color="#006699">Messages</font></h2>
                  <p>
                  <pre>${trainingInstance.message}</pre>
                  </p>
                  <hr>
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

