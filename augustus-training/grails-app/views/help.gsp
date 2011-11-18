


<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Help</title>         
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="show.gsp:" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[Help]</span></h1>
         </div>
         <div id="appnav">
             <ul>
                 <li><a href="../index.gsp" >Introduction</a></li>
                 <li><a href="../training/create.gsp">Training Submission</a></li>
                 <li><a href="../prediction/create.gsp">Prediction Submission</a></li>
                 <li><a href="show.gsp">Help</a></li>
                 <li><a href="../references.gsp">Links & References</a></li>
                 <li><a href="http://gobics.de/department/" title="Our department's homepage">Department</a></li>
             </ul>
         </div>
        <div class="body">
            <div class="main" id="main">
            This website contains instructions and frequently asked questions concerning <ul><li>the training of AUGUSTUS and <li>predicting genes in a new genome with pre-trained parameters.</ul>
            <hr>
            <h2>Why is the server buisy?</h2>
            Training AUGUSTUS is a very resource and time consuming process. Therefore, the server might be buisy. Wait one or two weeks and then try a new submission. If the problem persists, please contact <a href="mailto:augustus-training@gobics.de">augustus-training@gobics.de</a>.
            <hr>
            <div id="species_name"><h2>What is the species name?</h2></div>
                The species name is the name of the species for whose genome you want to train AUGUSTUS. If you do not want to reveal what kind of data you are submitting, you may use any string shorter than 30 characters as a species name.
            <hr>
            <div id="email"><h2>Why must I give my e-mail address?</h2></div>
            We need your e-mail address in order to send you a message when your job is finished. This may take a while (up to several weeks).
            <hr>
            <div id="upload_link"><h2>File upload versus web link</h2></div>
             The AUGUSTUS training and prediction web server applications offers in some cases two possiblities for transferring files to the server: 
             <ul>
                <li><b>Small files:</b> click on the Browse-button and select a file on your harddrive. If you experience a <i>Connection timeout</i>, please use the option for large files!</li>
                <li><b>Large files</b> can be retrieved from a <b>public</b> web link. Specify a valid ftp or http URL.</li>
            </ul>
            <b>You cannot do both!</b> For each file type (e.g. the genome file), you must <b>either</b> select a file on your harddrive <b>or</b> give a web link!
            <hr>
            <div id="which_files"><h2>Which files must or can I submitt for training AUGUSTUS?</h2></div>
            You need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li><b>at least one out of the</b>  <a  href="#cDNA">cDNA file</a>, <a href="#structure">training gene structure file</a>, and <a href="#protein">protein file</a>.</li>
            </ul>
            Please consider that training AUGUSTUS is a time and resource consuming process. Therefore, you should specify as much information as possible for a single training run instead of starting the AUGUSTUS training multiple times with different file combinations!
            <hr>
            <div id="which_files_pred"><h2>Which files are required for predicting genes in a new genome?</h2></div>
            For predicting genes in a new genome with already trained parameters, you need to specify
            <ul>
               <li>a <a href="#genome_file">genome file</a> and</li>
               <li>a <a href="#archive">parameter archive</a>. Instead of uploading the archive, you may also enter a valid <a href="#project_id">project identifier</a> in case you trained AUGUSTUS on this web server and the training has already finished.
            </ul>
            You may in addition specify an <a href="#cDNA">EST/cDNA file</a> and/or a <a href="#structure">hint file</a>.
            <hr>
            <div id="genome_file"><h2>Genome file</h2></div>
            The genome file is an obligatory file for training AUGUSTUS and for making predictions with pre-trained parameters in a new genome. It must contain the genome in (multiple) fasta format. Every header begins with a <b>></b>. The sequence must be DNA. Allowed sequence characters: <b>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn</b>. (Internally, AUGUSTUS will interpret everyting that is not AaTtCcGg as an N!) Empy lines are not allowed. If they occur, they will automatically be removed by the webserver applications.<br><br>
            <b>Example:</b>
            <pre class="example">
> Chr. 1
CCTCCTCCTGTTTTTCCCTCAATACAACCTCATTGGATTATTCAATTCAC
CATCCTGCCCTTGTTCCTTCCATTATACAGCTGTCTTTGCCCTCTCCTTC
TCTCGCTGGACTGTTCACCAACTCTCAGCCCGCGATCCCAATTTCCAGAC
AACCCATCTTATCAGCTTGGCCACGGCCTCGACCCGAACAGACCGGCGTC
CAGCGAGAAGAGCGTCGCCTCGACGCCTCTGCTTGACCGCACCTTGATGC
TCAAGACTTATCGCGATGCCAAGAAGCGTCTCATCATGTTCGACTACGA
> Chr 2
CGAAACGGGCACCTATACAACGATTGAAACCATTATTCAAGCTCAGCAAG
CGTCTATGCTAGCGGTTATTGCGAGCACTTCAGCGGTTGCTACTACGACT
ACTACTTGATAAATGAAACGGCTATAAAAGAGGCTGGGGCAAAAGTATGT
TAGTTGAAGGGTGACCTGAACGATGAATCGGTCGAATTTTTTATTGGCAG
AGGGAAGGTAGGTTTACTCAATTTAGTTACTTCTAGCCGTTGATTGGAGG
AGCGCAAGCGACGAGGAGGCTCATCGGCCGCCCGCGGAAAGCGTAGTCT
TACACGGAAATCAACGGCGGTGTCATAAGCGAG
> Chr 3
.....
            </pre>
            <hr>
            <div id="cDNA"><h2>cDNA file</h2></div>
            The cDNA file is a multiple fasta DNA file that contains e.g. ESTs or full-length cDNA sequences. Allowed sequence characters: <b>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNnUu</b>. Empty lines are not allowed and will be removed from the submitted file by the webserver application. See <a href="#genome_file">Genome file</a> for a format example.
            <hr>
            <div id="protein"><h2>Protein file</h2></div>
            The protein file is a multiple fasta file that contains protein sequences as supporting evidence for genes. Allowed sequence characters: <b>AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx</b>. Empty lines are not allowed but will simply be removed from the file by the webserver application. <br><br>
            <b>Example:</b>
            <pre class="example">
> protein 1
maaaafgqlnleepppiwgsrsvdcfekleqigegtygqvymakeiktgeivalkkirmd
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqvlhrdikgsnllidnegnlkl
adfglarsyshdhtgnltnrvitlwyrppelllgatkygp
>protein 2
neregfpitaireikilkklhhenvihlkeivtspgrdrddqgkpdnnkykggiymvfey
mdhdltgladrpglrftvpqikcymkqlltglhychvnqv
>protein 3
...
            </pre>
            <hr>
            <div id="structure"><h2>Training gene structure and hint file</h2></div>
            The training gene structure file or hint file contains given gene structures in gff-format with the following (obligatory)entries per line:<br><b>seqname source feature start end score strand frame attributes comments</b><br>
            <h3>Important details about gff-format for AUGUSTUS training</h3>
            <ul>
            <li>The columns are separated by <b>tabulators</b>. 
            <li>The <b>seqnames</b> in column 1 must be found in the fasta headers of sequences in the genome file. 
            <li>The <b>feature</b> (column 3) may be one of start, stop, exonpart, exon, dss, ass, or intronpart. 
            <li><b>Start</b> is the begin position, counting the first position of a sequence as position 1.
            <li><b>Stop</b> position, must be at least as large as start position.
            <li>The <b>score</b> must be a number but is irrelevant, here.
            <li><b>Frame</b> is the reading frame, can be denoted as '.' if unknown or irrelevant. For exonpart and exon this is as defined in the GFF format. On the forward strand it is the number of bases after (begin position - 1) until the next codon boundary comes (0, 1 or 2). On the reverse strand it is the number of bases before (end position + 1) until the next codon boundary comes (0, 1 or 2).
            <li><b>Attribute</b> must contain the string <b>source=M</b> (for manual). Other sources such EST or protein homology are possible but only in the command line version of AUGUSTUS. Then the hints may be ignored.
            </ul>
            <br>
	    <b>Example:</b>
            <pre class="example">
HS04636	anchor	exonpart	500	506	0	-	.	source=M
HS04636	anchor	exon	        966	1017	0	+	0	source=M
HS04636	anchor	start	        966	968	0	+	0	source=M
HS04636	anchor	dss	        2199	2199	0	+	.	source=M
HS04636	anchor	stop	        7631	7633	0	+	0	source=M
HS04636	anchor	intronpart	7631	7633	0	+	0	source=M
            </pre>
            </div>
            <hr>
            <div id="archive"><h2>Parameter archive</h2></div>
            A *.tar.gz archive with the following files is required for predicting genes in a new genome with pre-trained parameters:
            <ul>
               <li><i>species</i>_parameters.cfg</li>
               <li><i>species</i>_exon_probs.pbl</li>
               <li><i>species</i>_igenic_probs.pbl</li>
               <li><i>species</i>_intron_probs.pbl</li>
               <li><i>species</i>_weightmatrix.txt</li>
            </ul>
            where <i>species</i> is replaced by the name of the species you trained AUGUSTUS for (e.g. <i>carrot</i>). If you simply tar and gzip the results of an AUGUSTUS training run, everything should work fine.
            <hr>
            <div id="project_id"><h2>What is the project identifier?</h2></div>
            If you trained AUGUSTUS on this webserver, you may instead of uploading a parameter archive, simply specify the project identifier of this training run. You find the project identifier for example in the subject line for your training confirmation e-mail, where it says <i>Your AUGUSTUS training job project_id</i>. But you can also find the project identifier on the results webpage.
            <hr>
            <div id="job_status"><h2>What does my job status mean?</h2></div>
            Your job is in one of five stages:
            <ul>
               <li><b>stage 1</b>: submitted to webserver but not to cluster, yet. This may have one of the following reasons:
               <ul>
                 <li> file upload is still buisy
                 <li> file format validation is still in progress
                 <li>currently not enough computational resources available 
               </ul>
               <li><b>stage 2</b>: submitted to cluster and waiting for execution
               <li><b>stage 3</b>: calculating<br>
               This may take up to 3 weeks.
               <li><b>stage 4</b>: finished 
               <li><b>stage 5</b>: an error occurred
            </ul>
            <hr>
            <div id="utr"><h2>Why predicting the UTR?</h2>
            If a model for the untranslated regions (UTRs) is available for the species, they are included in the prediction. It takes significantly more time but in addition to reporting UTRs, it usually is also a little more accurate on the coding regions when ESTs are given as input.
            <hr>
            <div id="list"><h2>Why does my results page redirect me to a list of jobs?!</h2></div>
            We are trying to avoid data duplication. If you submitted some data that was already submitted before, by you or somebody else, we will delete your job. You receive an e-mail with a link to the results of training AUGUSTUS on your data (but that results page was - as stated in the beginning - created before you even submitted your job). <br>
            In the meantime your brower is redirected to a list of all available projects.
            <hr>
         </div>
<p style="text-align:right;">
            <small>Please direct your questions and comments to <a href="mailto:augustus-training@gobics.de">augustus-training@gobics.de</a></small>
            </p>
        </div>
        </div>
    </body>
</html>
