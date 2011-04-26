

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Create Training</title>         
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="index.gsp" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[Introduction]</span></h1>
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
            <div class="main" id="main">
            <h1>Welcome to the AUGUSTUS training web server</h1>

<p>AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences. This web server provides an interface for training AUGUSTUS on new genomes. It also enables you to predict genes in a genome sequence with already trained parameters.</p>

<p>AUGUSTUS usually belongs to the most accurate programs for the species it is trained for. Often it is the most accurate ab initio program. For example, at the independent gene finder assessment (EGASP) on the human ENCODE regions AUGUSTUS was the most accurate gene finder among the tested ab initio programs. At the more recent nGASP (worm), it was among the best in the ab initio and transcript-based categories. See <a href="http://augustus.gobics.de/accuracy">accuracy statistics</a> for further details.</p>

<p>For more information about AUGUSTUS, have a look at <a href="http://augustus.gobics.de/">the old AUGUSTUS web server</a>. There, you also find the <a href="http://augustus.gobics.de/binaries/">stand alone tool</a> for download. AUGUSTUS is already trained for a number of genomes and you find the according parameter sets at <a href="http://augustus.gobics.de/">the old web server</a>. Please check whether AUGUSTUS was already trained for your species before submitting a new training job.</p>
            </div>
        </div>
<p style="text-align:right;">
            <small>Please direct your questions and comments to <a href="mailto:augustus-training@gobics.de">augustus-training@gobics.de</a></small>
            </p>
    </body>
</html>
