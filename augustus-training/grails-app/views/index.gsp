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
                        <li class="navigation-sub__item"><a href="trainingtutorial.gsp">Training Tutorial</a></li>
                        <li class="navigation-sub__item"><g:link controller="training" action="create">Submit Training</g:link></li>
                        <li class="navigation-sub__item"><a href="predictiontutorial.gsp">Prediction Tutorial</a></li>
                        <li class="navigation-sub__item"><g:link controller="prediction" action="create">Submit Prediction</g:link></li>
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
                        <h1 class="csc-firstHeader">Welcome to the WebAUGUSTUS Service</h1>
                     </div>
                  </div>
                  <div id="c261665" class="csc-default">
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
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

