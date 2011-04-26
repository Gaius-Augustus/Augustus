

<html>
    <head>
        <META HTTP-EQUIV="Refresh" CONTENT="20">
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>List of training jobs</title>
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="create.gsp:" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[List of Training Jobs]</span></h1>
         </div>
         <div id="appnav">
             <ul>
                 <li><a href="../index.gsp" >Introduction</a></li>
                 <li><g:link controller="training" action="create">Training Submission</g:link></li>
                 <li><g:link controller="prediction" action="create">Prediction Submission</g:link></li>
                 <li><g:link controller="help" action="list">Help</g:link></li>
                 <li><a href="../references.gsp">Links & References</a></li>
                 <li><a href="http://gobics.de/department/" title="Our department's homepage">Department</a></li>
             </ul>
         </div>
        <div class="body">
            <g:if test="${flash.message}">
            <div class="message">${flash.message}</div>
            </g:if>
            <div class="list">
                <table frame="border" rules="all" cellpadding="5">
                    <thead>
                        <tr>
                   	        <g:sortableColumn property="id" title="Id" />

                                <g:sortableColumn
                                property="accession_id"
                                title="Project ID" />

                                <g:sortableColumn
                                property="dateCreated"
                                title="Date" />

                                <g:sortableColumn
                                property="job_status"
                                title="Job Status" />
                        
                        </tr>
                    </thead>
                    <tbody>
                    <g:each in="${trainingInstanceList}" status="i" var="trainingInstance">
                        <tr class="${(i % 2) == 0 ? 'odd' : 'even'}">
                        
                            <td><g:link action="show" id="${trainingInstance.id}">${fieldValue(bean:trainingInstance, field:'id')}</g:link></td>
                        
                            <td>${fieldValue(bean:trainingInstance, field:'accession_id')}</td>
                        
                            <td>${fieldValue(bean:trainingInstance, field:'dateCreated')}</td>
                        
                            <td>${fieldValue(bean:trainingInstance, field:'job_status')}</td>

                        </tr>
                    </g:each>
                    </tbody>
                </table>
            </div>
        </div>
    </body>
</html>
