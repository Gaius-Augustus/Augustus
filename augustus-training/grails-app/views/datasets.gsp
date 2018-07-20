

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
                        <h1 class="csc-firstHeader">Datasets for Download</h1>
                     </div>
                  </div>
                  <div id="c261665" class="csc-default">
                     <div class="csc-default">
                        <div class="divider">
                           <hr>
                        </div>
                     </div>
                     <p>The following sequence files were used to <em>train</em> AUGUSTUS
                        or to <em>test</em> its accuracy. Some of the datasets are described in the
                        paper &ldquo;Gene Prediction with a Hidden Markov Model and a new Intron Submodel&rdquo;,
                        which was presented at the European Conference on Computational Biology
                        in September 2003 and appeared in the proceedings.
                     </p>
                     <h2>Test sets:</h2>
                     <br>
                     <h3>human:</h3>
                     <br>
                     <h4>h178:</h4>
                     <p>178 single-gene short human sequences <br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/h178.gb.gz">h178.gb.gz</a> (gzipped genbank format)
                     </p>
                     <h4>sag178:</h4>
                     <p>
                        semi artificial genomic sequences from Guigo et al.:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/sag178.gb.gz">sag178.gb.gz</a> (gzipped genbank format)<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/sag178.fa.gz">sag178.fa.gz</a> (gzipped fasta format)<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/sag178.gff">sag178.gff</a> (annotation in gff format)
                     </p>
                     <h3>fly:</h3>
                     <br>
                     <h4>fly100:</h4>
                     <p>100 single gene sequences from FlyBase:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/fly100.gb.gz">fly100.gb.gz</a> (gzipped Genbank format)
                     </p>
                     <h4>adh122:</h4>
                     <p>
                        A 2.9 Mb long sequence from the Drosophila adh region (copied from the
                        <a href="http://www.fruitfly.org/GASP1/data/standard.html">GASP dataset page</a>)<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/adh.fa.gz">adh.fa.gz</a> (gzipped fasta format)<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/adh.std1.gff_corrected">adh.std1.gff_corrected</a> (gff format)<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/adh.std1+3.gff">adh.std1+3.gff</a> (gff format)
                     </p>
                     <h3>Arabidopsis thaliana:</h3>
                     <p>Araset. 74 sequences with 168 genes.<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/araset.gb.gz">araset.gb.gz</a> (gzipped genbank format)
                     </p>
                     <h2>Training sets:</h2>
                     <br>
                     <h3>human:</h3>
                     <p>
                        single gene sequences from genbank (1284 genes):<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/human.train.gb.gz">human.train.gb.gz</a> (gzipped genbank format)
                     <p>11739 human splice sites, originally from Guig&oacute; et al., but filtered for similarities to h178, sag178:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/splicesites.gz">splicesites.gz</a> (gzipped flat file)
                     </p>
                     <h3>fly:</h3>
                     <p>
                        320 single gene sequences from FlyBase, disjoint with fly100:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/fly.train.gb.gz">fly.train.gb.gz</a> (gzipped genbank format)
                     </p>
                     <p>
                        400 single gene sequences from FlyBase, disjoint with adh122:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/adh.train.gb.gz">adh.train.gb.gz</a> (gzipped genbank format)
                     </p>
                     <h3>Arabidopsis:</h3>
                     <p>249 single gene sequences obtained by deleting the sequences from the Araball set which overlap with the sequences from Araset:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/araball.train.gb.gz">araball.train.gb</a> (gzipped Genbank format)
                     </p>
                     <h3>Coprinus cinereus (a fungus):</h3>
                     <p>
                        851 single gene sequences predicted by genewise and compiled by Jason Stajich. 261 genes are complete, 590 genes are incomplete at the 3' end.
                        Genes redundand with those in the Genbank annotations were deleted:<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/cop.genomewise.gb.gz">cop.genomewise.gb.gz</a> (gzipped Genbank format)
                     </p>
                     <p>
                        91 sequences containing 93 genes from Genbank. Genes in Genbank with nothing else than the coding sequence were omitted. Identical or extremely 
                        similar genes in genbank were used only once. This set has first been used as a test set for above training set. The Coprinus version
                        here used :<br>
                        <a href="http://bioinf.uni-greifswald.de/augustus/datasets/cop.gb.clean.gb.gz">cop.gb.clean.gb.gz</a> (gzipped Genbank format)
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

