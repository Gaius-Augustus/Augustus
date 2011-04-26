

<html>
    <head>
        <META HTTP-EQUIV="Refresh" CONTENT="5">
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Show Training</title>
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="create.gsp:" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[Training Results]</span></h1>
         </div>
         <div id="appnav">
             <ul>
                 <li><a href="../../index.gsp" >Introduction</a></li>
                 <li><g:link controller="training" action="create">Training Submission</g:link></li>
                 <li><a href="../../prediction/create.gsp">Prediction Submission</a></li>
                 <li><g:link controller="help" action="list">Help</g:link></li>
                 <li><a href="../../references.gsp">Links & References</a></li>
                 <li><a href="http://gobics.de/department/" title="Our department's homepage">Department</a></li>
             </ul>
         </div>
        <div class="body">
            <g:if test="${flash.message}">
            <div class="message">${flash.message}</div>
            </g:if>

            <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '0' || fieldValue(bean:trainingInstance, field:'job_status') == '1' || fieldValue(bean:trainingInstance, field:'job_status') == '2' || fieldValue(bean:trainingInstance, field:'job_status') == '3'}">
            	<h1>Status of Job ${fieldValue(bean:trainingInstance, field:'accession_id')}</h1>
            
            Your job is in stage <g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '0'|| fieldValue(bean:trainingInstance, field:'job_status') == '1'}">1</g:if><g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '2'}">2</g:if><g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '3'}">3</g:if><g:if test = "${fieldValue(bean:trainingInstance, field:'job_status') == '4'}">4</g:if>.<br><br>Explannation:
            <ul>
                <li>stage 1: submitted to webserver but not to cluster, yet
                <li>stage 2: submitted to cluster and waiting for execution
                <li>stage 3: calculating
                <li>stage 4: finished
                <li>stage 5: error
            </ul>
            For more details, see <g:link controller="help" action="list" fragment="job_status">Help</g:link>
            </g:if>
            

            <g:if test ="${fieldValue(bean:trainingInstance, field:'job_status') == '4'}">
            <h1>Results for Job ${fieldValue(bean:trainingInstance, field:'accession_id')}</h1>
               <ul>
                  <li>Your job is finished.
               </ul>
            </g:if>

        </div>
    </body>
</html>
