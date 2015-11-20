

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Create Training</title>         
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
        <a href="http://www.uni-greifswald.de" target="_blank" class="mainleveltop_" >University of Greifswald</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.mnf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Faculty</a><span class="mainleveltop_">&nbsp;|&nbsp; </span><a href="http://www.math-inf.uni-greifswald.de/" target="_blank" class="mainleveltop_" >Institute</a><span class="mainleveltop_">&nbsp;|&nbsp;</span><a href="http://bioinf.uni-greifswald.de/" target="_blank" class="mainleveltop_">Bioinformatics Group</a>      </td>
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
    References
</span>

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
         <li id="current"><a href="references.gsp"><span>Links & References</span></a></li>
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
    <table class="contentpaneopen">
      <tr>
	<td class="contentheading" width="100%">
	  <font color="#006699">References & Links</font>
        </td>
      </tr>
    </table>
            <g:if test="${flash.message}">
            <div class="message">${flash.message}</div>
            </g:if>
            <g:hasErrors bean="${trainingInstance}">
            <div class="errors">
                <g:renderErrors bean="${trainingInstance}" as="list" />
            </div>
            </g:hasErrors>
            <h2>Link to the old AUGUSTUS web server</h2>
            <p>The old AUGUSTUS web server for predictions with pre-trained models is available at <a href="http://augustus.gobics.de">http://augustus.gobics.de</a>.</p>
            <h2>References</h2>
            <p>Please cite the following references when using the AUGUSTUS training web server results in your publication:</p>
	    <ul>
	      <li>K. J. Hoff and M. Stanke (2013)<br>
		<a href="http://nar.oxfordjournals.org/content/early/2013/05/21/nar.gkt418.full?keytype=ref&ijkey=Alnjr4zdAka9FGo">WebAUGUSTUS - a web service for training AUGUSTUS and predicting genes in eukaryotes</a><br>
		Nucleic Acids Res, doi:10.1093/nar/gkt418</li>
	      <li>K. J. Hoff and M. Stanke (2012)<br>
            <a href="http://bioinf.uni-greifswald.de/trainaugustus/references/PAG2012.pdf">TrainAUGUSTUS - a Webserver Application for Parameter Training and Gene Prediction in Eukaryoties</a><br>
            International Plant & Animal XX Conference 2012, U.S.A. (Poster)</li>
	      </ul>

	    <p>Please cite the following references when using results obtained by AUGUSTUS in your publication:</p>
            <ul>
            <li>D. Sommerfeld, T. Lingner, M. Stanke, B. Morgenstern, H. Richter (2009)<br>
<a href="http://www.elsevier.com/wps/find/journaldescription.cws_home/505611/description">AUGUSTUS at MediGRID: Adaption of a Bioinformatics Application to Grid Computing for Efficient Genome Analysis</a><br>
Future Generation Computer Systems 25, 337 - 345. 
            </li>
            <li> Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008)<br>
      <a href="http://bioinformatics.oxfordjournals.org/cgi/reprint/btn013?ijkey=AqOiFZBiTC5VhDS&keytype=ref">Using native and syntenically mapped cDNA alignments to improve de novo gene finding</a><br>
      Bioinformatics, doi: 10.1093/bioinformatics/btn013</li>
            <li>Mario Stanke, Ana Tzvetkova, Burkhard Morgenstern (2006)<br>
      <a href="http://genomebiology.com/2006/7/S1/S11">"AUGUSTUS at EGASP: using EST, protein and genomic alignments for improved gene prediction in the human genome"</a><br>
      BMC Genome Biology, 7(Suppl 1):S11.</li>
            <li>M. Stanke , O. Schöffmann , B. Morgenstern, S. Waack (2006)<br>
      <a href="http://www.biomedcentral.com/1471-2105/7/62/abstract">Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources</a><br>
      BMC Bioinformatics 7, 62.</li>
            <li>Mario Stanke and Burkhard Morgenstern (2005)<br>
      <a href="http://nar.oxfordjournals.org/cgi/reprint/33/suppl_2/W465?ijkey=MZ7JCYGHlXzWLcX&keytype=ref">"AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints"</a>,<br>
      Nucleic Acids Research, 33, W465-W467</li>
    <li>Mario Stanke, Rasmus Steinkamp, Stephan Waack and Burkhard Morgenstern (2004)<br>
      <a href="http://nar.oxfordjournals.org/cgi/content/full/32/suppl_2/W309">"AUGUSTUS: a web server for gene finding in eukaryotes"</a><br>
      Nucleic Acids Research, Vol. 32, W309-W312</li>
    <li>Mario Stanke and Stephan Waack (2003)<br>
      <a href="http://bioinformatics.oxfordjournals.org/cgi/content/abstract/19/suppl_2/ii215">Gene Prediction with a Hidden-Markov Model and a new Intron Submodel.<a/><br>
      Bioinformatics, Vol. 19, Suppl. 2, pages ii215-ii225</li>
    <li>Mario Stanke (2003)<br>
      <a href="http://webdoc.sub.gwdg.de/diss/2004/stanke/index.html">Gene Prediction with a Hidden-Markov Model</a>.<br>
      Ph.D. thesis, Universität Göttingen</li>
           </ul>

      <h2>Genome papers, where AUGUSTUS was used in collaboration</h2>
          <ul>
	    <li>
K. K. Dasmahapatra, J. R. Walters, A. D. Briscoe, J. W. Davey, A. Whibley, N. J. Nadeau, A. V. Zimin, D. S. T. Hughes, L. C. Ferguson, S. H. Martin, C. Salazar, J. J. Lewis, <b>S. Adler</b>, S.-J. Ahn, D. A. Baker, S. W. Baxter, N. L. Chamberlain, R. Chauhan, B. A. Counterman, T. Dalmay, L. E. Gilbert, K. Gordon, D. G. Heckel, H. M. Hines, <b>K. J. Hoff</b>, P. W. H. Holland, E. Jacquin-Joly, F. M. Jiggins, R. T. Jones, D. D. Kapan, P. Kersey, G. Lamas, D. Lawson, D. Mapleson, L. S. Maroja, A. Martin, S. Moxon, W. J. Palmer, R. Papa, A. Papanicolaou, Y. Pauchet, D. A. Ray, N. Rosser, S. L. Salzberg, M. A. Supple, A. Surridge, A. Tenger-Trolander, H. Vogel, P. A. Wilkinson, D. Wilson, J. A. Yorke, F. Yuan, A. L. Balmuth, C. Eland, K. Gharbi, M. Thomson, R. A. Gibbs, Y. Han, J. C. Jayaseelan, C. Kovar, T. Mathew, D. M. Muzny, F. Ongeri, L.-L. Pu, J. Qu, R. L. Thornton, K. C. Worley, Y.-Q. Wu, M. Linares, M. L. Blaxter, R. H. ffrench-Constant, M. Joron, M. R. Kronforst, S. P. Mullen, R. D. Reed, S. E. Scherer, S. Richards, J. Mallet, W. Owen McMillan and C. D. Jiggins (2012)<br>
<a href="http://www.nature.com/nature/journal/vaop/ncurrent/full/nature11041.html">Butterfly genome reveals promiscuous exchange of mimicry adaptations among species</a><br>
Nature 2012, doi:10.1038/nature11041
	      </li>
    <li>M. Berriman, B.J. Haas, P.T. LoVerde, R.A. Wilson, G.P. Dillon, G.C. Cerqueira, S.T. Mashiyama, B. Al-Lazikani, L.F. Andrade, P.D. Ashton, M.A. Aslett, D.C. Bartholomeu, G. Blandin, C.R. Caffrey, A. Coghlan, R. Coulson, T.A. Day, A. Delcher, R. DeMarco, A. Djikeng, T. Eyre, J.A. Gamble, E. Ghedin, Y. Gu, C. Hertz-Fowler, H. Hirai, Y. Hirai, R. Houston, A. Ivens, D.A. Johnston, D. Lacerda, C.D. Macedo, P. McVeigh, Z. Ning, G. Oliveira, J.P. Overington, J. Parkhill, M. Pertea, R.J. Pierce, A.V. Protasio, M.A. Quail, M.-A. Rajandream, J. Roger, M. Saji, S.L. Salzberg, <b>M. Stanke</b>, A.R. Tivey, O. Whit, D.L. William, J. Wortman, W. Wu, M. Zamanian, A. Zerlotini, C.M. Fraser-Liggett, B.G. Barrell, N.M. El-Sayed (2009)<br>
<a href="http://www.nature.com/nature/journal/v460/n7253/abs/nature08160.html">The Genome of the blood fluke Schistosoma mansoni</a><br>
Nature 460, 352-358 
</li>
    <li>Tribolium Genome Sequencing Consortium,<br> <a href="http://www.nature.com/nature/journal/v452/n7190/abs/nature06784.html">The genome of the model beetle and pest Tribolium castaneum</a>.<br> Nature, March 2008.</li>
    <li> Nene, Vishvanath, Wortman, Jennifer R., Lawson, Daniel, Haas, Brian, Kodira, Chinnappa, Tu, Zhijian Jake, Loftus, Brendan, Xi, Zhijong, Megy, Karyn, Grabherr, Manfred, Ren, Quinghu, Zdobnov, Evgeny M., Lobo, Neil F., Campbell, Kathryn S., Brown, Susan E., Bonaldo, Maria F., Zhu, Jingsong, Sinkins, Steven P., Hogenkamp, David G., Amedo, Paulo, Arsenburger, Peter, Atkinson, Peter W., Bidwell, Shelby, Biedler, Jim, Birney, Ewan, Bruggner, Robert V., Costas, Javier, Coy, Monique R., Crabtree, Jonathan, Crawford, Matt, deBruyn, Becky, DeCaprio, David, Eiglmeier, Karin, Eisenstadt, Eric, El-Dorry, Hamza, Gelbart, William M., Gomes, Suely L., Hammond, Martin, Hannick, Linda I., Hogan, James R., Holmes, Michael H., Jaffe, David, Johnston, Spencer J., Kennedy, Ryan C., Koo, Hean, Kravitz, Saul, Kriventseva, Evgenia V., Kulp, David, LaButti, Kurt, Lee, Edward, Li, Song, Lovin, Diane D., Mao, Chunhong, Mauceli, Evan, Menck, Carlos F. M., Miller, Jason R., Montgomery, Philip, Mori, Akio, Nascimento, Ana L., Naveira, Horacio F., Nusbaum, O'Leary, Sinead B., Orvis, Joshua, Pertea, Mihaela, Quesneville, Hadi, Reidenbach, Kyanne R., Rogers, Yu-Hui, Roth, Charles W., Schneider, Jennifer R., Schatz, Michael, Shumway, Martin, <b>Stanke, Mario</b>, Stinson, Eric O., Tubio, Jose M. C., VanZee, Janice P., Verjovski-Almeida, Sergio, Werner, Doreen, White, Owen, Wyder, Stefan, Zeng, Qi, Zhao, Qi, Zhao, Yongmei, Hill, Catherine A., Raikhel, Alexander S., Soares, Marcelo B., Knudson, Dennis L., Lee, Norman H., Galagan, James, Salzberg, Steven L., Paulsen, Ian T., Dimopoulos, George, Collins, Frank H., Bruce, Birren, Fraser-Liggett, Claire M., Severson, David W. (2007)<br>
      <a href="http://www.sciencemag.org/cgi/content/abstract/1138878v1">Genome Sequence of Aedes aegypti, a Major Arbovirus Vector</a><br>
      Science, doi: 10.1126/science.1138878</li>
    <li>Elodie Ghedin, Shiliang Wang, David Spiro, Elisabet Caler, Qi Zhao, Jonathan Crabtree, Jonathan E. Allen, Arthur L. Delcher, David B. Guiliano, Diego Miranda-Saavedra, Samuel V. Angiuoli, Todd Creasy, Paolo Amedeo, Brian Haas, Najib M. El-Sayed, Jennifer R. Wortman, Tamara Feldblyum, Luke Tallon, Michael Schatz, Martin Shumway, Hean Koo, Steven L. Salzberg, Seth Schobel, Mihaela Pertea, Mihai Pop, Owen White, Geoffrey J. Barton, Clotilde K. S. Carlow, Michael J. Crawford, Jennifer Daub, Matthew W. Dimmic, Chris F. Estes, Jeremy M. Foster, Mehul Ganatra, William F. Gregory, Nicholas M. Johnson, Jinming Jin, Richard Komuniecki, Ian Korf, Sanjay Kumar, Sandra Laney, Ben-Wen Li, Wen Li, Tim H. Lindblom, Sara Lustigman, Dong Ma, Claude V. Maina, David M. A. Martin, James P. McCarter, Larry McReynolds, Makedonka Mitreva, Thomas B. Nutman, John Parkinson, Jose M. Peregri-Alvarez, Catherine Poole, Qinghu Ren, Lori Saunders, Ann E. Sluder, Katherine Smith, <b>Mario Stanke</b>, Thomas R. Unnasch, Jenna Ware, Aguan D. Wei, Gary Weil, Deryck J. Williams, Yinhua Zhang, Steven A. Williams, Claire Fraser-Liggett, Barton Slatko, Mark L. Blaxter, and Alan L. Scott (21 September 2007)<br>
      <a href="http://www.sciencemag.org/cgi/content/full/317/5845/1756">Draft Genome of the Filarial Nematode Parasite Brugia malayi</a><br>
      Science, 317 (5845), 1756, doi: 10.1126/science.1145406</li>

</ul>

<div id="contrib"><h2>User contributed species parameters:</h2></div>

<ul>
<li>
The training annotation files of the following species are a courtesy of Jason Stajich:<br>
<i>Aspergillus fumigatus, Aspergillus nidulans, Aspergillus oryzae, Aspergillus terreus, Botrytis cinerea, Candida albicans,
Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Cryptococcus neoformans gattii,
Cryptococcus neoformans neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum,
Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa,
Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Ustilago maydis, Yarrowia lipolytica</i>.
</li>
<li>
The training for lamprey (<i>Petromyzon marinus</i>) was performed by Falk Hildebrand and Shigehiro Kuraku, based on the
genome assembly (PMAR3.0) provided by the Genome Sequencing Center at Washington University School of Medicine (WUGSC)
in St. Louis.
</li>
<li>
The training for elephant shark (<i>Callorhinchus milii</i>) was performed by Tereza Manousaki and Shigehiro Kuraku, based on the genome assembly
(made up of 1.4x whole genome shotgun reads) available at <a href="http://esharkgenome.imcb.a-star.edu.sg/resources.html">http://esharkgenome.imcb.a-star.edu.sg/resources.html</a>.
</li>
<li>
The training for <i>Pneumocystis jirovecii</i> was performed by Marco Pagni, Philippe Hauser et al as described in Hauser PM, Burdet FX, Cisse OH,
Keller L, Taffe P, Sanglard D, Pagni M., Comparative Genomics Suggests that the Fungal Pathogen Pneumocystis Is an Obligate Parasite
Scavenging Amino Acids from Its Host's Lungs. PLoS One. 2010, Dec 20;5(12):e15152. PubMed PMID: 21188143; PubMed Central PMCID: PMC3004796.
</li>
</ul>

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
      <a href="mailto:augustus-web@uni-greifswald.de" title="E-Mail augustus-web@uni-greifswald.de, opens the standard mail program">augustus-web@uni-greifswald.de</a>
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
