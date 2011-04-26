

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Create Training</title>         
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="help.gsp" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[References & Links]</span></h1>
         </div>
         <div id="appnav">
             <ul>
                 <li><a href="index.gsp" >Introduction</a></li>
                 <li><a href="training/create.gsp">Training Submission</a></li>
                 <li><a href="prediction/create.gsp">Prediction Submission</a></li>
                 <li><g:link controller="help" action="list">Help</g:link></li>
                 <li><a href="references.gsp">Links & References</a></li>
                 <li><a href="http://gobics.de/department/" title="Our department's homepage">Department</a></li>
             </ul>
         </div>
        <div class="body">
            <g:if test="${flash.message}">
            <div class="message">${flash.message}</div>
            </g:if>
            <g:hasErrors bean="${trainingInstance}">
            <div class="errors">
                <g:renderErrors bean="${trainingInstance}" as="list" />
            </div>
            </g:hasErrors>
            <div class="main" id="main">
            <h2>Link to the old AUGUSTUS web server</h2>
            The old AUGUSTUS web server for predictions with pre-trained models is available at <a href="http://augustus.gobics.de">http://augustus.gobics.de</a>.
            <h2>References</h2>
            Please cite the following references when using the AUGUSTUS training web server results in your publication:
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
            </div>
    </body>
</html>
