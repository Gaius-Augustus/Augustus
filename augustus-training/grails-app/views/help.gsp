<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Help</title>         
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
    Help
</span>

  <div class="beendeFluss"></div>
</div>
<!-- ***** Ende: Kopfbereich *********************************************// -->

<!-- ***** Start: Koerper ************************************************// -->
<div id="koerper">

  <div id="linke_spalte">
    <ul class="menu">
         <li><a href="index.gsp"><span>Introduction</span></a></li>
         <li><a href="trainingtutorial.gsp"><span>Training Tutorial</span></a></li>
         <li><g:link controller="training" action="create"><span>Submit Training</span></g:link></li>
         <li><a href="predictiontutorial.gsp"><span>Prediction Tutorial</span></a></li>
         <li><g:link controller="prediction" action="create"><span>Submit Prediction</span></g:link></li>
         <li id="current"><a href="help.gsp"><span>Help</span></a></li>
         <li><a href="references.gsp"><span>Links & References</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de"><span>Bioinformatics Group</span></a></li>
         <li><a href="http://bioinf.uni-greifswald.de/bioinf/impressum.html"><span>Impressum</span></a></li>
     </ul>
  </div>
 <div id="mittel_spalte">
<div class="main" id="main">
   <h1><a href="help.gsp">Help</a></h1>
            <div class="main" id="main">
            <p>This website contains <i>short instructions</i> and some <i>frequently asked questions</i> concerning <ul><li>the training of AUGUSTUS and <li>predicting genes in a new genomic sequence with pre-trained parameters.</ul> <br>For more detailed instructions, please read <a href="trainingtutorial.gsp">Training Tutorial</a> and <a href="predictiontutorial.gsp">Prediction Tutorial</a>.</p>
            <hr>
            <br>
	    <h2>Contents</h2>
	    <p>
	      <a href="#noResults">Why do I not get any results?</a><br>
	      <a href="#buisy">Why is the server busy?</a><br>
              <a href="#species_name">What is the species name?</a><br>
	      <a href="#email">Why must I give my e-mail address?</a><br>
              <a href="#upload_link">File upload versus web link</a><br>
              <a href="#which_files">Which files must or can I submit for training AUGUSTUS?</a><br>
              <a href="#which_files_pred">Which files are required for predicting genes in a new genome?</a><br>
              <a href="#genome_file">Genome file</a><br>
              <a href="#cDNA">cDNA file</a><br>
              <a href="#protein">Protein file</a><br>
              <a href="#structure">Training gene structure and hint file</a><br>
              <a href="#hints">Hints file</a><br>
              <a href="#archive">Parameter archive</a><br>
              <a href="#project_id">What is the project identifier?</a><br>
              <a href="#job_status">What does my job status mean?</a><br>
              <a href="#utr">UTR prediction: yes or no?</a><br>
              <a href="#allowedGeneStructure">Allowed gene structure<a><br>
              <a href="#list">Why does my job not exist?</a><br>
            </p>
            <hr>
	    <br>
	    <h2 id="noResults">Why do I not get any results?</h2>
	    <p>The quality of results depends on the quality and combination of you input data. If the input data did e.g. not provide sufficient information for generating training genes, then no AUGUSTUS paramteres will be optimized for your species, and no predictions will be made. In case of the gene prediction web server application, it is also possible that your submitted genome sequence does not contain any protein coding genes.</p>
<p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <h2 id="buisy">Why is the server busy?</h2>
            <p>Training AUGUSTUS is a very resource and time consuming process. We use a grid engine queuing system with a limited number of waiting slots. If we estimate that the time from job submission to computation start might be very long, our web server might display a message that our server is buisy. The submission of new jobs is then disabled. Please wait one or two weeks before you try a new submission. If the problem persists longer than a month, please contact <a href="mailto:augustus-web@uni-greifswald.de">augustus-web@uni-greifswald.de</a>.</p>
	    <p><a href="#seitenanfang">
	      <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
	      Top of page
	    </a>
	    <br>
	    </p>
            <hr>
	    <br>
            <div id="species_name"><h2>What is the species name?</h2></div>
<p>The species name is the name of the species for whose genome you want to train AUGUSTUS. The species name is an obligatory parameter that AUGUSTUS needs in order to find the correct parameters that shall be applied for predicting genes in a specified genomic sequence. Considering that AUGUSTUS training is such a time consuming process, our objective is to know the names of species for which AUGUSTUS was trained in order to make the trained parameters available to the public so that others who are interested in the same species as you do not have to rerun the training process.</p>

<p>However, if you do not want to reveal the true species name, you may use any other string shorter than 30 characters as a species name. Species names must be unique on our system, i.e. if the string of your choice is already existing in our system, you will get a message that you have to choose another species name.</p>

<p>We are <b>not</b> redistributing the original sequence data that you submitted to our web server application. But we are redistributing the trained parameters and the species name and any other kind of results that your computation may have produced.</p><p><b>Example:</b></p><p>If person 1 submits a sequence data set for training and names it <i>hypothetical_species</i>, and a second person tries later to train AUGUSTUS with exactly the same sequence files (new data upload) and names the species <i>some_other_name</i>, the second person will be redirected to the results of the original training run of person 1, and the species name <i>hypothetical_species</i> will be publicly readable to person 2.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="email"><h2>Why must I give my e-mail address?</h2></div>
            <p> We save use your e-mail address for the following purposes:</p>

<p>
<ul>
<li>Confirming your job submission</li>
<li>Confirming successful file upload (for large files)</li>
<li>Sending you link to the results of your job (you will never get the link to the results if you do not enter a <b>valid</b> address to which you have access)</li>
<li>Informing you about any problems that might occur during your particular job</li>
</ul>
</p>

<p>We do <b>not</b> use your e-mail address to send you any <i>spam</i>, i.e. about web service updates. We do<b>not</b> share your e-mail addresses with any third parties.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="upload_link"><h2>File upload versus web link</h2></div>
            <p> The AUGUSTUS training and prediction web server application offers in some cases two possiblities for transferring files to the server: <i>Upload a file</i> and <i>specify a web link to file</i>.</p>
             <p>
             <ul>
                <li>For <b>small files</b>, please click on the Browse-button and select a file on your harddrive.<br>If you experience a <i>Connection timeout</i> (because your file was too large for this type of upload), please use the option for large files!</li>
                <li><b>Large files</b> can be retrieved from a <b>public</b> web link. Specify a valid ftp or http URL to your sequence file. Our server will fetch the file from the given address.</li>
            </ul>
	     </p>
	    <p>
            <b>You cannot do both at the same time!</b> For each file type (e.g. the genome file), you must <b>either</b> select a file on your harddrive <b>or</b> give a web link!</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="which_files"><h2>Which files must or can I submit for training AUGUSTUS?</h2></div>
            <p>You need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li><b>at least one out of the following files:</b>  <a  href="#cDNA">cDNA file</a>, <a href="#structure">training gene structure file</a>, and <a href="#protein">protein file</a>.</li>
            </ul>
	    </p>
	    <p>
            Please consider that training AUGUSTUS is a time and resource consuming process. For optimal results, you should specify as much information as possible for a single training run instead of starting the AUGUSTUS training multiple times with different file combinations!</p>
	                <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="which_files_pred"><h2>Which files are required for predicting genes in a new genome?</h2></div>
            <p>For predicting genes in a new genome with already trained parameters, you need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li>a <a href="#archive">parameter archive</a>. Instead of uploading the archive, you may also enter a valid <a href="#project_id">project identifier</a> in case you trained AUGUSTUS on this web server and the training has already finished.
            </ul></p>
	    <p>
            You may in addition specify an <a href="#cDNA">EST/cDNA file</a> and/or a <a href="#structure">hints file</a> that will be used as extrinsic evidence for predicting genes.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="genome_file"><h2>Genome file</h2></div>
            <p>The genome file is an obligatory file for training AUGUSTUS and for making predictions with pre-trained parameters in a new genome. It must contain the genome in (multiple) fasta format. Every header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n</b>. (Internally, AUGUSTUS will interpret everyting that is not <b>A a T t C c G g</b> as an <b>N</b>!) Empty lines are not allowed. If they occur, they will automatically be removed by the webserver applications.<br><br>
            <b>Correct file format example:</b>
            <pre class="example">
>Chr.1
CCTCCTCCTGTTTTTCCCTCAATACAACCTCATTGGATTATTCAATTCAC
CATCCTGCCCTTGTTCCTTCCATTATACAGCTGTCTTTGCCCTCTCCTTC
TCTCGCTGGACTGTTCACCAACTCTCAGCCCGCGATCCCAATTTCCAGAC
AACCCATCTTATCAGCTTGGCCACGGCCTCGACCCGAACAGACCGGCGTC
CAGCGAGAAGAGCGTCGCCTCGACGCCTCTGCTTGACCGCACCTTGATGC
TCAAGACTTATCGCGATGCCAAGAAGCGTCTCATCATGTTCGACTACGA
>Chr.2
CGAAACGGGCACCTATACAACGATTGAAACCATTATTCAAGCTCAGCAAG
CGTCTATGCTAGCGGTTATTGCGAGCACTTCAGCGGTTGCTACTACGACT
ACTACTTGATAAATGAAACGGCTATAAAAGAGGCTGGGGCAAAAGTATGT
TAGTTGAAGGGTGACCTGAACGATGAATCGGTCGAATTTTTTATTGGCAG
AGGGAAGGTAGGTTTACTCAATTTAGTTACTTCTAGCCGTTGATTGGAGG
AGCGCAAGCGACGAGGAGGCTCATCGGCCGCCCGCGGAAAGCGTAGTCT
TACACGGAAATCAACGGCGGTGTCATAAGCGAG
>Chr.3
.....
            </pre>
	    </p>
	                <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="cDNA"><h2>cDNA file</h2></div>
            <p>The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>A a T t G g C c H h X x R r Y y W w S s M m K k B b V v D d N n U u</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. See <a href="#genome_file">Genome file</a> for a format example. Upload of a cDNA file to our web server applications will invoke the software <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="trainingtutorial.gsp#ref2">2</a>], which is on our webserver application only available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b>.</p>
	    <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="protein"><h2>Protein file</h2></div>
            <p>The protein file is a multiple fasta file that contains protein sequences as supporting evidence for genes. Allowed sequence characters: <b>A a R r N n D d C c E e Q q G g H h I i L l K k M m F f P p S s T t W w Y y V v B b Z z J j X x</b>. Empty lines are not allowed but will simply be removed from the file by the webserver application. <br><br></p>
            <p><b>Correct file format example:</b>
            <pre class="example">
>protein1
maaaafgqlnleepppiwgsrsvdcfekleqigegtygqvymakeiktgeivalkkirmd
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqvlhrdikgsnllidnegnlkl
adfglarsyshdhtgnltnrvitlwyrppelllgatkygp
>protein2
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqv
>protein3
...
            </pre></p>
	   <p>Submitting a protein file to our AUGUSTUS training web server application will invoke <a href="http://www.webscipio.org/">Scipio</a> [<a href="trainingtutorial.gsp#ref3">3</a>], which uses <a href="http://genome.ucsc.edu/cgi-bin/hgBlat?command=start">BLAT</a> [<a href="trainingtutorial.gsp#ref2">2</a>]. Therefore, protein file upload is only available <b><font color="f40b0b">for academic, personal and  non-profit use</font></b> on our web server application.</p>
           <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="structure"><h2>Training gene structure file</h2></div>
	    <p>You can submit your own, externally created training gene structures to the AUGUSTUS training web server application.</p>

<p>Training gene structure files can be submitted in two different formats: Genbank format or gff format.</p>

<h3>Training gene structure file in genbank format</h3>

<p>Gene structures in genbank format must contain the coding sequence parts and flanking regions. Flanking regions are important because AUGUSTUS is supposed to differentiate between genes and intergenic regions. The length of flanking regions depends on the length of genes in the target genome. In our pipeline, flanking regions are set to the average gene length (exceptionally applying the extreme limits between 1000 and 10000 nt). It is very important to make sure that the flanking regions do not contain any other protein coding gene parts, i.e. we recommend to trim flanking regions in a way that will exclude other CDS parts.</p>

<p><b>Correct file format example (condensed view, the three dots represent further lines of sequences):</b>
            <pre class="example">
LOCUS       Chr.1_1-159458   159458 bp  DNA
FEATURES             Location/Qualifiers
     source          1..159458
     CDS             complement(join(2421..2655,3858..4005,4080..4235,5569..5857
                     ,10316..10534,155240..155458))
                     /gene="1474336"
BASE COUNT     49195 a   29117 c  28985 g   49950 t   2211 n
ORIGIN
        1 aaaatacatc acaatacatt taattcactt tccatcatcg agattaacga aaattattta
       61 aaatatcgaa gatgaaaata tcctcaagat gatactgaac ggctaagaaa aatacatcac
      121 acaactttaa ttcattttcc atcatcgaga ttaacgaaaa gaaaaaattt taactcccta
...
   159301 atacgccacc aggtatttcg cctgattgtt cctcgaatat cttctctctc tctatatata
   159361 tatatattac ttggcacgat aatcgtcgaa tcgttattta taaattgctt catctatcgc
   159421 gatatttttg caacaactct cgcttttctc tttccatt
//
LOCUS       Chr.1_313992-323129   9138 bp  DNA
FEATURES             Location/Qualifiers
     source          1..9138
     CDS             join(4001..4048,4989..5138)
                     /gene="194551"
BASE COUNT     2829 a   1502 c  1750 g   2948 t   109 n
ORIGIN
        1 ttttccttct ttcttttttt tttatttaca ttaatgagaa ttttcgcaaa tatttcatcg
       61 ctgccatcct tttttttcct cgacgtcaat cacgcgacac atttgttaga gaaatggatt
      121 ttaatcttga aaaaagaaaa atacaaatgc caacgcattt caaatccttt cctattatta
...
     9001 tcaacgaaac aaataattgc ttcacaaaat atcgcacgta acaacaatat agacttcaat
     9061 attcaacaat tcttttcctt tatacacaaa gatacacaaa atataaaagt tttaatactt
     9121 caacttcaac gaaacagg
//
            </pre></p>


<p>If you want to train UTRs, you have to additionally incorporate mRNA information in your genbank file.</p>

<p><b>Correct file format example (including UTR training):</b>
            <pre class="example">
LOCUS       scf7180001240730_g20   526 bp  DNA
FEATURES             Location/Qualifiers
     source          1..526
     mRNA            99..125
     CDS             99..99
BASE COUNT     164 a   99 c  68 g   195 t
ORIGIN
        1 gtgacggagc ccaaggacga gcccgtgccc tcagagccca cgtccgacgt gaggcccgcg
       61 ccagcgcccc tcccgccgcc cgtcgcagcc actgcttaga ctttactaat ataaacattg
      121 aaaatatttt gtgttttatt tccaatcatt gaattataat cctattataa tataactaac
      181 attcgtaatt ttacaaaata actatgcaaa ttattttgta ttttcgtttt aaattatact
      241 tttcatataa atttctacaa atcttattca agaccataag tatccgctcg ctctacttcg
      301 ggcatttcct ttatttatat cttatttgac ttattttgat tatttaggct tatgttttcg
      361 atactattga aaacagaaaa taatttcata taattaataa tatattttca attaatatat
      421 ttaacaaata tttgtatagt tcaagcggac aaatccgttc ccatagtatt tatataaatt
      481 ttaatttaga gtaataacag tttgctgtat tgttgtagtc aaatac
//
LOCUS       scf7180001240751_g30   876 bp  DNA
FEATURES             Location/Qualifiers
     source          1..876
     mRNA            complement(401..777)
     CDS             complement(777..777)
BASE COUNT     300 a   136 c  116 g   324 t
ORIGIN
        1 aatgtaggaa aatgaaatat ttatttaaat tgttattatc acttcttcgc tctagtgtct
       61 tggcaaagcg cggcgttgag ttcagcctct cacacgcaat gcctccagaa ttcggcgaaa
      121 tgtgggggac agagtgtatt aacactaagt tccctcagcc acgactggtg aaattatata
      181 ttcagtttgt atactattac tcatgcaaac acttcatcat actttcactc aatcagtaaa
      241 gcataatatt ttatttaata ttgtttatca atactatttc cttgttgtta aatattattt
      301 tatttattat attaaattaa aatgtcaaaa ttaaaagtag gtgatgattt attactatct
      361 tttctatcca agaaaaaaaa gacacactga aacaattgta atttttgtta tgtttttatt
      421 acttaatatt attataaaaa tttgtaaata cgaaataaaa tagatagacg taataatatt
      481 tatttgttag ttaataataa taatgataat tacgaaagat acaagaaata tgcataaatg
      541 agtgttatat tatgtatttt atgagaatat aaatataaaa actgtcattg attatatttt
      601 ctaaatactt tcattttatg gcttgctggc ttttcaattt ccttatgttt cagcttttca
      661 ctcaatagag cgaaaccttc atcgacatgt aagccaatag aacaattaca aactaacttt
      721 attacatcag tcttttcatt tctttaagct tcaggcaaat atcatctaaa tgcctttcaa
      781 ctcgctacta acatcgcgtc gttatataaa tcagtgtata cggaattaaa cctgtcatgt
      841 ctcttgcaag acgtgtctgc tgttgtcacg cacaca
//
            </pre></p>

<h3>Training gene structure file in gff format</h3>

<p>Training gene structure in gff format must comply with the fasta entry names of the genome file.</p>

<p>In general, gff format must contain the following columns (The columns are separated by <b>tabulators</b>):</p>

<p>
            <OL TYPE="1">
            <li>The <b>sequence names</b> must be found in the fasta headers of sequences in the genome file. 
	    <li>The <b>source</b> tells with which software/process the gene structure was generated (you can fill in whatever you like).
            <li>The <b>feature</b> may for AUGUSTUS training be CDS, 5'-UTR or 3'-UTR. 
            <li><b>Start</b> is the beginning position of the line's feature, counting the first position of a sequence as position 1.
            <li><b>Stop</b> position, must be at least as large as start position.
            <li>The <b>score</b> must be a number but the number is irrelevant to our web server applications.
	    <li>The <b>strand</b> denotes whether the gene is located on the forward (+) or on the reverse (-) strand.
            <li><b>Frame</b> is the reading frame, can be denoted as '.' if unknown or irrelevant. For exonpart and exon this is as defined as follows: On the forward strand it is the number of bases after (begin position 1) until the next codon boundary comes (0, 1 or 2). On the reverse strand it is the number of bases before (end position + 1) the next codon boundary comes (0, 1 or 2).
            <li><b>Attribute</b> contains a <i>transcript identifier</i>. All gff-entries belonging to one transcript must contain the same transcript identifier in the last column.
            </OL>
</p>

<p><b>Correct file format example (without UTR):</b>
            <pre class="example">
Chr.1	mySource	CDS	1767	1846	1.000	-	0	transcript_id "1597_1"
Chr.1	mySource	CDS	1666	1709	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	1486	1605	1.000	-	2	transcript_id "1597_1"
Chr.1	mySource	CDS	1367	1427	1.000	-	2	transcript_id "1597_1"
Chr.1	mySource	CDS	1266	1319	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	1145	1181	1.000	-	1	transcript_id "1597_1"
Chr.1	mySource	CDS	847	1047	1.000	-	0	transcript_id "1597_1"
Chr.2	mySource	CDS	9471	9532	1.000	+	0	transcript_id "1399_2"
Chr.2	mySource	CDS	9591	9832	1.000	+	1	transcript_id "1399_2"
Chr.2	mySource	CDS	9885	10307	1.000	+	2	transcript_id "1399_2"
Chr.2	mySource	CDS	10358	10507	1.000	+	2	transcript_id "1399_2"
Chr.2	mySource	CDS	10564	10643	1.000	+	2	transcript_id "1399_2"
</pre></p>

<p><b>Correct file format example (with UTR):</b>
            <pre class="example">
Chr.1	mySource	5'-UTR	277153	277220	45	+	.	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	277221	277238	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278100	278213	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	278977	279169	1	+	0	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279630	279648	0.94	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	279734	279768	0.94	+	1	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	CDS	280307	280344	1	+	2	transcript_id "g22472.t1"; gene_id "g22472";
Chr.1	mySource	3'-UTR	280345	280405	78	+	.	transcript_id "g22472.t1"; gene_id "g22472";
</pre></p>


           <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
	    <br>
            <div id="hints"><h2>Hints file</h2></div>

	    <p>For the gene prediction web server application, it is possible to submit an externally created file that contains extrinsic evidence for gene structures in gff format.</p>
	    <p>In general, gff format must contain the following columns (The columns are separated by <b>tabulators</b>):</p>

<p>
            <OL TYPE="1">
            <li>The <b>sequence names</b> must be found in the fasta headers of sequences in the genome file. 
	    <li>The <b>source</b> tells with which software/process the gene structure was generated (you can fill in whatever you like).
            <li>The <b>feature</b> may for AUGUSTUS training be CDS, 5'-UTR or 3'-UTR. 
            <li><b>Start</b> is the beginning position of the line's feature, counting the first position of a sequence as position 1.
            <li><b>Stop</b> position, must be at least as large as start position.
            <li>The <b>score</b> must be a number but the number is irrelevant to our web server applications.
	    <li>The <b>strand</b> denotes whether the gene is located on the forward (+) or on the reverse (-) strand.
            <li><b>Frame</b> is the reading frame, can be denoted as '.' if unknown or irrelevant. For exonpart and exon this is as defined as follows: On the forward strand it is the number of bases after (begin position 1) until the next codon boundary comes (0, 1 or 2). On the reverse strand it is the number of bases before (end position + 1) the next codon boundary comes (0, 1 or 2).
            <li>For usage as hint, <b>Attribute</b> must contain the string <b>source=M</b> (for manual). Other sources, such EST or protein, are possible, but only in the command line version of AUGUSTUS. Source types other than <b>M</b> are ignored by AUGUSTUS web server applications.
            </OL>
</p>
	    <p>
            <br>
            <b>Correct format example:</b>
            <pre class="example">
HS04636 anchor  exonpart        500     506     0       -       .       source=M
HS04636 anchor  exon            966     1017    0       +       0       source=M
HS04636 anchor  start           966     968     0       +       0       source=M
HS04636 anchor  dss             2199    2199    0       +       .       source=M
HS04636 anchor  stop            7631    7633    0       +       0       source=M
HS04636 anchor  intronpart      7631    7633    0       +       0       source=M
            </pre>
            </p>
	    <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <div id="archive"><h2>Parameter archive</h2></div>
	    <p>
            A *.tar.gz archive with a folder containing the following files is required for predicting genes in a new genome with pre-trained parameters:</p>
	    <p>
<ul>
<li><i>species</i>/<i>species</i>_parameters.cfg
<li><i>species</i>/<i>species</i>_metapars.cfg
<li><i>species</i>/<i>species</i>_metapars.utr.cfg
<li><i>species</i>/<i>species</i>_exon_probs.pbl.withoutCRF
<li><i>species</i>/<i>species</i>_exon_probs.pbl
<li><i>species</i>/<i>species</i>_weightmatrix.txt
<li><i>species</i>/<i>species</i>_intron_probs.pbl
<li><i>species</i>/<i>species</i>_intron_probs.pbl.withoutCRF
<li><i>species</i>/<i>species</i>_igenic_probs.pbl
<li><i>species</i>/<i>species</i>_igenic_probs.pbl.withoutCRF
</ul>
</p>
<p>where <i>species</i> is replaced by the name of the species you trained AUGUSTUS for (e.g. carrot would result it <i>carrot</i>/<i>carrot</i>_parameters.cfg). The additional <i>species</i> before the slash means that all those files must reside in a directory that is called <i>species</i> before you tar and gzip it. If you simply tar and gzip the folder that contains parameters of an AUGUSTUS training run, everything should work fine.</p>
                       <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="project_id"><h2>What is the project identifier?</h2></div>
            <p>If you trained AUGUSTUS on this webserver, you may instead of uploading a parameter archive, simply specify the project identifier of this training run. You find the project identifier for example in the subject line for your training confirmation e-mail, where it says <i>Your AUGUSTUS training job project_id</i>.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="job_status"><h2>What does my job status mean?</h2></div>
<p>In the beginning, the status page will display that your job has been <b>submitted</b>. This means, the web server application is currently uploading your files and validating file formats. After a while, the status will change to <b>waiting for execution</b>. This means that all file formats have been confirmed and an actually AUGUSTUS training job has been submitted to our grid engine, but the job is still pending in the queue. Depending on waiting queue length, this status may persist for a while. Please contact us in case you job is pending for more than one month. Later, the job status will change to <b>computing</b>. This means the job is currently computing. When the page displays <b>finished</b>, all computations have been finished and a website with your job's results has been generated.</p>

<p>You will receive an e-mail with the link to the results of your job when computations are finished.</p>
            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="utr"><h2>UTR prediction: yes or no?</h2>
            <p>It takes significantly more time to predict UTRs but in addition to reporting UTRs, it usually is also a little more accurate on the coding regions when ESTs are given as extrinsic evidence.</p>

<p>UTR prediction is only possible if UTR parameter files exist for your species. Even if UTR parameter files exist for a species, you should make sure, that they are <i>species specific</i>, i.e. have actually been optimized for your target species. It is a waste of time to predict UTRs with <i>general</i> (template) parameters.</p>

<p>If no UTR parameter files exist for your species but you enables UTR prediction in the form, the web server application will overrule the choice to predict UTRs by simply not predicting any UTRs.</p>

            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="allowedGeneStructure"><h2>Allowed gene structure</h2>

<p><b>Predict any number of (possibly partial) genes:</b> This option is set by default. AUGUSTUS may predict no gene at all, one or more genes. The genes at the boundaries of the input sequence may be partial. Partial here means that not all of the exons of a gene are contained in the input sequence, but it is assumed that the sequence starts or ends in a non-coding region.<p>

<p><b>Predict only complete genes:</b> AUGUSTUS assumes that the input sequence does not start or end within a gene. Zero or more complete genes are predicted.</p>

<p><b>Predict only complete genes - at least one:</b> As the previous option. But AUGUSTUS predicts at least one gene (if possible).</p>

<p><b>Predict exactly one complete gene:</b> AUGUSTUS assumes that the sequence contains exactly one complete gene. Note: This feature does not work properly in combination with alternative transcripts. </p>

<p><b>Ignore conflicts with other strand:</b> By default AUGUSTUS assumes that no genes - even on opposite strands - overlap. Indeed, this usually is the case but sometimes an intron contains a gene on the opposite strand. In this case, or when AUGUSTUS makes a false prediction on the one strand because it falsely thinks there is a conflicting gene on the other strand, AUGUSTUS should be run with this option set. It then predicts the genes on each strand separately and independently. This may lead to more false positive predictions, though.</p>


            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            <br>
            <div id="list"><h2>Why does my job not exist?</h2></div>
            <p>We are trying to avoid data duplication. If you submitted some data that was already submitted before, by you or somebody else, we will delete your job. You receive an e-mail with a link to the job-status and results of training AUGUSTUS on your data. The results link will only be functional in case the computations of that previously submitted job have already finished.</p>

            <p><a href="#seitenanfang">
              <img hspace="5" height="4" border="0" width="7" alt="Seitenanfang" src="images/top.gif" />
              Top of page
            </a>
            <br>
            </p>
            <hr>
            </div></div>
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
