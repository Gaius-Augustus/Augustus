

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <meta name="layout" content="main" />
        <title>Create Training</title>         
    </head>
    <body>
        <div class="headline" id="headline">
                <h1 class="title" id="title"><a href="create.gsp:" id="applink" name="applink">AUGUSTUS</a> <span class="subtitle" id="subtitle">[Training Submission]</span></h1>
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
            <g:hasErrors bean="${trainingInstance}">
            <div class="errors">
                <g:renderErrors bean="${trainingInstance}" as="list" />
            </div>
            </g:hasErrors>
            <g:if test="${flash.error}">
                 <div class="errors">
                    &nbsp; <img src="../images/skin/exclamation.png"> &nbsp; ${flash.error}
                 </div>
            </g:if>
            <div class="main" id="main">
            <g:uploadForm action="commit" method="post" >
            <fieldset><legend><b>Data Input for Training AUGUSTUS</b></legend><p>
                <div class="dialog">
                    <table>
                        <tbody>
                        
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="email_adress">E-mail</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'email_adress','errors')}">
                                    <input type="text" id="email_adress" name="email_adress" value="${fieldValue(bean:trainingInstance,field:'email_adress')}"/> &nbsp; <g:link controller="help" action="list" fragment="email"><small>Help</small></g:link>
                                </td> 
                            </tr> 
                        
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="project_name">Species name</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'project_name','errors')}">
                                    <input type="text" maxlength="30" id="project_name" name="project_name" value="${fieldValue(bean:trainingInstance,field:'project_name')}"/> &nbsp; <g:link controller="help" action="list" fragment="species_name"><small>Help</small></g:link>
                                </td>
                            </tr>
                          </tbody> 
                      </table>
                      <br>
                      You may <b>either</b> upload data files from your computer <b>or</b> specify web links. &nbsp; <g:link controller="help" action="list" fragment="upload_link"><small>Help</small></g:link>
                      <br>
                      <br>
                      <table>
                         <tbody>
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="GenomeFile"><b>Genome file</b></label>
                                </td>
                                <td valitn="top">
                                    <input type="file" id="GenomeFile" name="GenomeFile"/>
                                </td>
                                <td>&nbsp;<b>or</b>&nbsp;</td>
                                <td valign="top" class="name">
                                    <label for="genome_ftp_link">web link to genome file</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'genome_ftp_link','errors')}">
                                    <input type="text" id="genome_ftp_link" name="genome_ftp_link" value="${fieldValue(bean:trainingInstance,field:'genome_ftp_link')}"/> &nbsp; <g:link controller="help" action="list" fragment="genome_file"><small>Help</small></g:link>
                                </td>
                            </tr> 
                          </tbody>
                        </table>
                        <br>
                        You need to specify <b>at least one</b> of these files: <g:link controller="help" action="list" fragment="which_files"><small>Help</small></g:link><br><br>
                        <table>
                          <tbody>
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="EstFile">cDNA file</label>
                                </td>
                                <td valign="top">
                                    <input type="file" id="EstFile" name="EstFile"/>
                                </td>
                            <td>&nbsp;<b>or</b>&nbsp;</td>
                                <td valign="top" class="name">
                                    <label for="est_ftp_link">web link to cDNA file</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'est_ftp_link','errors')}">
                                    <input type="text" id="est_ftp_link" name="est_ftp_link" value="${fieldValue(bean:trainingInstance,field:'est_ftp_link')}"/> &nbsp; <g:link controller="help" action="list" fragment="cDNA"><small>Help</small></g:link>
                                </td>
                            </tr> 
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="protein_file">Protein file</label>
                                </td>
                                <td valign="top">
                                    <input type="file" id="ProteinFile" name="ProteinFile"/>
                                </td>
                            <td>&nbsp;<b>or</b>&nbsp;</td>
                                <td valign="top" class="name">
                                    <label for="protein_ftp_link">web link to protein file</label>
                                </td>
                                <td valign="top" class="value ${hasErrors(bean:trainingInstance,field:'protein_ftp_link','errors')}">
                                    <input type="text" id="protein_ftp_link" name="protein_ftp_link" value="${fieldValue(bean:trainingInstance,field:'protein_ftp_link')}"/> &nbsp; <g:link controller="help" action="list" fragment="protein"><small>Help</small></g:link>
                                </td>
                            </tr> 
                            <tr class="prop">
                                <td valign="top" class="name">
                                    <label for="struct_file">Training gene structure file</label>
                                </td>
                                <td valign="top">
                                    <input type="file" id="StructFile" name="StructFile"/> 
                                </td>
                            <td><g:link controller="help" action="list" fragment="structure"><small>Help</small></g:link></td>
                                <td valign="top" class="name">
                                   
                                </td>
                                <td valign="top">
                                    
                                </td>
                            </tr> 
                        </tbody>
                    </table>
                    <br><br>
                </div>
                <div class="buttons">
                    <span class="button"><input class="commit" type="submit" value="Start Training" /></span>
                </div>
            </g:uploadForm>
            </div>
            <p>&nbsp;</p>
            <p style="text-align:right;">
            <small>Please direct your questions and comments to <a href="mailto:augustus-training@gobics.de">augustus-training@gobics.de</a></small>
            </p>
        </div>
    </body>
</html>
