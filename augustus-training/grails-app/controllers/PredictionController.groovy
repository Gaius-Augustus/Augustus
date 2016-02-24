// The class PredictionController controls everything that is related to submitting a job for predicting genes with pre-trained parameters on a novel genome
//    - it handles the file upload (or wget)
//    - format check
//    - SGE job submission and status checks
//    - rendering of results/job status page
//    - sending E-Mails concerning the job status (submission, errors, finished)

import java.io.BufferedReader
import java.io.File
import java.io.FileInputStream
import java.io.FileReader
import java.io.IOException
import java.io.InputStream
import java.io.InputStreamReader
import java.io.OutputStream
import java.io.Reader
import java.io.StreamTokenizer
import java.io.StringWriter
import java.io.Writer
import java.net.URL
import java.util.zip.GZIPInputStream
import java.net.UnknownHostException

class PredictionController {
	// need to adjust the output dir to whatever working dir! This is where uploaded files and results will be saved.
	def output_dir = "/data/www/augpred/webdata" // should be something in home of webserver user and augustus frontend user.
//	def output_dir = "/data/www/test"
	// this log File contains the "process log", what was happening with which job when.
	def logFile = new File("${output_dir}/pred.log")
	// this log File contains the "database" (not identical with the grails database and simply for logging purpose)
	def dbFile = new File("${output_dir}/augustus-pred-database.log")
	// oldID is a parameter that is used for show redirects (see bottom)
	def oldID
	def oldAccID
	// web-output, root directory to the results that are shown to end users
	def web_output_dir = "/var/www/trainaugustus/prediction-results" // must be writable to webserver application
//	def web_output_dir = "/data/www/test/out"
	def web_output_url = "http://bioinf.uni-greifswald.de/trainaugustus/prediction-results/"
	def war_url = "http://bioinf.uni-greifswald.de/webaugustus/"
	def footer = "\n\n------------------------------------------------------------------------------------\nThis is an automatically generated message.\n\nhttp://bioinf.uni-greifswald.de/webaugustus" // footer of e-mail
	// AUGUSTUS_CONFIG_PATH
	def AUGUSTUS_CONFIG_PATH = "/usr/local/augustus/trunks/config"
	def AUGUSTUS_SCRIPTS_PATH = "/usr/local/augustus/trunks/scripts"
	def BLAT_PATH = "/usr/local/blat/blat"
	def scaffold = Prediction
	// Admin mail for errors
	def admin_email = "katharina.hoff@gmail.com"
	// sgeLen length of SGE queue, when is reached "the server is buisy" will be displayed
	def sgeLen = 20;
	// max button filesize
	def long maxButtonFileSize = 104857600 // 100 MB = 13107200 bytes = 104857600 bit, getFile etc. gives size in bit
	def long preUploadSize
	// max ftp/http filesize
	def long maxFileSizeByWget = 1073741824 // 1 GB = 1073741824 bytes, curl gives size in bytes
	// EST sequence properties (length)
	def int estMinLen = 250
	def int estMaxLen = 20000
	// logging verbosity-level
	def verb = 2 // 1 only basic log messages, 2 all issued commands, 3 also script content
	def cmd2Script
	def cmdStr
	def logDate
        def maxNSeqs = 250000 // maximal number of scaffolds allowed in genome file
        // other variables
	def accession_id
	def prokaryotic = false // flag to determine whether augustus should be run in prokaryotic mode
	// human verification:
	def simpleCaptchaService

	// check whether the server is buisy
	def beforeInterceptor = {
		def String prefixChars ="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_"
		def rnd = new Random()
		def qstatFilePrefix = (1..10).sum{prefixChars[ rnd.nextInt(prefixChars.length()) ]} 
		def qstatFile = new File("${output_dir}/${qstatFilePrefix}.qstatScript")
		cmd2Script = "qstat -u \"*\" | grep qw | wc -l > ${output_dir}/${qstatFilePrefix}.qstatResult 2> /dev/null"
		qstatFile << "${cmd2Script}"
		if(verb > 2){
			logDate = new Date()
			logFile << "${logDate} SGE          v3 - qstatFile << \"${cmd2Script}\"\n"
		}
		if(!qstatFile.exists()){
			logDate = new Date()
			logFile << "SEVERE ${logDate} SGE          v1 - ${qstatFile} does not exist!\n"
		}
		cmdStr = "bash ${output_dir}/${qstatFilePrefix}.qstatScript"
		def qstatStatus = "${cmdStr}".execute()
		if(verb > 2){
			logDate = new Date()
			logFile << "${logDate} SGE          v3 - \"${cmdStr}\"\n"
		}
		qstatStatus.waitFor()
		def qstatStatusResult = new File("${output_dir}/${qstatFilePrefix}.qstatResult").text
		def qstatStatus_array = qstatStatusResult =~ /(\d*)/
		def qstatStatusNumber 
		(1..qstatStatus_array.groupCount()).each{qstatStatusNumber = "${qstatStatus_array[0][it]}"}
		cmdStr = "rm -r ${output_dir}/${qstatFilePrefix}.qstatScript &> /dev/null"
		def delProc = "${cmdStr}".execute()
		if(verb > 2){
			logDate = new Date()
			logFile << "${logDate} SGE          v3 - \"${cmdStr}\"\n"
		}
		delProc.waitFor()
		if(qstatFile.exists()){
			//logDate = new Date()
			//logFile << "SEVERE ${logDate} SGE          v1 - ${qstatFile} was not deleted!\n"
		}
		cmdStr = "rm -r ${output_dir}/${qstatFilePrefix}.qstatResult &> /dev/null"
		delProc = "${cmdStr}".execute()
		if(verb > 2){
			logDate = new Date()
			logFile << "${logDate} SGE          v3 - \"${cmdStr}\"\n"
		}
		delProc.waitFor()
		if(qstatStatusNumber > sgeLen){
			// get date
			def todayTried = new Date()
			// get IP-address
			String userIPTried = request.remoteAddr
			logDate = new Date()
			logFile <<  "${logDate} SGE          v1 - On ${todayTried} somebody with IP ${userIPTried} tried to invoke the Prediction webserver but the SGE queue was longer than ${sgeLen} and the user was informed that submission is currently not possible\n"
			render "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\"/><meta name=\"layout\" content=\"main\" /><title>Submitt Prediction</title><script type=\"text/javascript\" src=\"js/md_stylechanger.js\"></script></head><body><!-- Start: Kopfbereich --><p class=\"unsichtbar\"><a href=\"#inhalt\" title=\"Directly to Contents\">Directly to Contents</a></p><div id=\"navigation_oben\"><a name=\"seitenanfang\"></a><table width=\"100%\" border=\"0\" cellpadding=\"0\" cellspacing=\"1\"><tr><td nowrap=\"nowrap\"><a href=\"http://www.uni-greifswald.de\" target=\"_blank\" class=\"mainleveltop_\" >University of Greifswald</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.mnf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Faculty</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.math-inf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Institute</a><span class=\"mainleveltop_\">&nbsp;|&nbsp;</span><a href=\"http://bioinf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\">Bioinformatics Group</a></td></tr></table></div><div id=\"banner\"><div id=\"banner_links\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php\" title=\"Institut f&uuml;r Mathematik und Informatik\"><img src=\"../images/header.gif\" alt=\"Directly to home\" /> </a></div><div id=\"banner_mitte\"><div id=\"bannertitel1\">Bioinformatics Web Server at University of Greifswald</div><div id=\"bannertitel2\">Gene Prediction with AUGUSTUS</div></div><div id=\"banner_rechts\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php/geschichte-und-kultur/167\" title=\"Voderberg-Doppelspirale\"><img src=\"../images/spirale.gif\" align=\"left\" /></a></div></div><div id=\"wegweiser\">Navigation for: &nbsp; &nbsp;<span class=\"breadcrumbs pathway\">Submitt Prediction</span><div class=\"beendeFluss\"></div></div><!-- Ende: Kopfbereich --><!-- Start: Koerper --><div id=\"koerper\"><div id=\"linke_spalte\"><ul class=\"menu\"><li><div id=\"linksMenuText\">AUGUSTUS Web Server Navigation</div></li><li><a href=\"../index.gsp\"><span>Introduction</span></a></li><li><a href=\"../about.gsp\"><span>About AUGUSTUS</span></a></li><li><a href=\"../accuracy.gsp\"><span>Accuracy</span></a></li><li><a href=\"../trainingtutorial.gsp\"><span>Training Tutorial</span></a></li><li><a href=\"/webaugustus/training/create\"><span>Submit Training</span></a></li><li><a href=\"../predictiontutorial.gsp\"><span>Prediction Tutorial</span></a></li><li><a href=\"/webaugustus/prediction/create\"><span>Submit Prediction</span></a></li><li><a href=\"../help.gsp\"><span>Help</span></a></li><li><a href=\"../datasets.gsp\"><span>Datasets for Download</span></a></li><li><a href=\"../predictions_for_download.gsp\"><span>Predictions for Download</span></a></li><li><a href=\"../references.gsp\"><span>Links & References</span></a></li><li><a href=\"http://bioinf.uni-greifswald.de/bioinf/impressum.html\"><span>Impressum</span></a></li><li>&nbsp;</li><li><div id=\"linksMenuText\">Other AUGUSTUS Resources</div></li><li><a href=\"http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus\">AUGUSTUS Wiki</a></li><li><a href=\"http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Forum.Forum\">AUGUSTUS Forum</a></li><li><a href=\"http://bioinf.uni-greifswald.de/augustus/binaries/\">Download AUGUSTUS</a></li><li><a href=\"http://bioinf.uni-greifswald.de/augustus\">Old AUGUSTUS gene prediction web server</a></li><li>&nbsp;</li><li><div id=\"linksMenuText\">Other Links</div></li><li><a href=\"http://bioinf.uni-greifswald.de\"><span>Bioinformatics Group Greifswald</span></a></li></ul></div><div id=\"mittel_spalte\"><div class=\"main\" id=\"main\"><h1><font color=\"#006699\">The Server is Busy</font></h1><p>You tried to access the AUGUSTUS prediction job submission page.</p><p>Predicting genes with AUGUSTUS is a process that takes a lot of computation time. We estimate that one prediction process requires at most approximately 7 days. Our web server is able to process a certain number of jobs in parallel, and we established a waiting queue. The waiting queue has a limited length, though. Currently, all slots for computation and for waiting are occupied.</p><p>We apologize for the inconvenience! Please try to submitt your job later.</p><p>Feel free to contact us in case your job is particularly urgent.</p></div><p>&nbsp;</p>           </div><div id=\"rechte_spalte\"><div class=\"linien_div\"><h5 class=\"ueberschrift_spezial\">CONTACT</h5><strong>Institute for Mathematics und Computer Sciences</strong><br/><strong>Bioinformatics Group</strong><br />Walther-Rathenau-Stra&szlig;e 47<br />17487 Greifswald<br />Germany<br />Tel.: +49 (0)3834 86 - 46 24<br/>Fax:  +49 (0)3834 86 - 46 40<br /><br /><a href=\"mailto:augustus-web@uni-greifswald.de\" title=\"E-Mail augustus-web@uni-greifswald.de, opens the standard mail program\">augustus-web@uni-greifswald.de</a></div></div><div class=\"beendeFluss\"></div></div><!-- Ende: Koerper --><!-- Start: Fuss --><div id=\"fuss\"><div id=\"fuss_links\"><p class=\"copyright\">&copy; 2011 University of Greifswald</p></div><div id=\"fuss_mitte\"><div class=\"bannergroup\"></div></div><div id=\"fuss_rechts\" ><ul><li><a href=\"#seitenanfang\"><img hspace=\"5\" height=\"4\" border=\"0\" width=\"7\" alt=\"Seitenanfang\" src=\"../images/top.gif\" />Top of page</a></li></ul></div><div class=\"beendeFluss\"></div></div><!-- Ende: Fuss --></body></html>"
			return
		}		
	}

	// the method commit is started if the "Submit Job" button on the website is hit. It is the main method of Prediction Controller and contains a Thread method that will continue running as a background process after the user is redirected to the job status page.

	def fillSample = {
		redirect(action:create, params:[genome_ftp_link:"http://bioinf.uni-greifswald.de/trainaugustus/examples/LG16.fa",project_id:"honeybee1"])
	}



	def commit = {
		def predictionInstance = new Prediction(params)
		if(!(predictionInstance.id == null)){
			flash.error = "Internal error 2. Please contact augustus-web@uni-greifswald.de if the problem persists!"
			redirect(action:create)
			return
		}else{
			// retrieve parameters of form for early save()
			def uploadedGenomeFile = request.getFile('GenomeFile')
			def uploadedParamArch = request.getFile('ArchiveFile')
	      		def uploadedEstFile = request.getFile('EstFile')
			def uploadedStructFile = request.getFile('HintFile')
			if(!(uploadedGenomeFile.empty)){
	        		predictionInstance.genome_file = uploadedGenomeFile.originalFilename
			}
			if(!(uploadedParamArch.empty)){
				predictionInstance.archive_file = uploadedParamArch.originalFilename
			}
			if(!(uploadedEstFile.empty)){
				predictionInstance.est_file = uploadedEstFile.originalFilename
			}
			if(!(uploadedStructFile.empty)){
				predictionInstance.hint_file = uploadedStructFile.originalFilename
			}
			predictionInstance.save()
			// info string for confirmation E-Mail
			def confirmationString
			def mailStr
			confirmationString = "Prediction job ID: ${predictionInstance.accession_id}\n"
			predictionInstance.job_id = 0
			// define flags for file format check, file removal in case of failure
			def archiveExistsFlag = 0
			def speciesNameExistsFlag = 0
			def genomeFastaFlag = 0
			def estFastaFlag = 0
			def estExistsFlag = 0
			def hintGffFlag = 0
			def hintExistsFlag = 0
			def overRideUtrFlag = 0
			def metacharacterFlag = 0
			// species name for AUGUSTUS
			def species
			// delProc is needed at many places
			def delProc
			def st
			def content
			def int error_code
			def urlExistsScript
			def sgeErrSize = 10
			def writeResultsErrSize = 10
			def msgStr
			// get date
			def today = new Date()
			logFile << "${today} ${predictionInstance.accession_id} v1 - AUGUSTUS prediction webserver starting on ${today}\n"
      			// get IP-address
      			String userIP = request.remoteAddr
			logDate = new Date()
      			logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - user IP: ${userIP}\n"

			// flag for redirect to submission form, display warning in appropriate places
			predictionInstance.warn = true
			// parameters for redirecting
			def redirParams=[:] 
			if(predictionInstance.email_adress != null){
				redirParams["email_adress"]="${predictionInstance.email_adress}"
			}
			if(predictionInstance.genome_ftp_link != null){
				redirParams["genome_ftp_link"]="${predictionInstance.genome_ftp_link}"
			}
			if(predictionInstance.est_ftp_link != null){
				redirParams["est_ftp_link"]="${predictionInstance.est_ftp_link}"
			}
			if(predictionInstance.project_id != null){
				redirParams["project_id"]="${predictionInstance.project_id}"
			}
			if(predictionInstance.genome_file != null){
				redirParams["has_genome_file"]="${predictionInstance.warn}"
			}
			if(predictionInstance.est_file != null){
				redirParams["has_est_file"]="${predictionInstance.warn}"
			}
			if(predictionInstance.archive_file != null){
				redirParams["has_param_file"]="${predictionInstance.warn}"
			}
			if(predictionInstance.hint_file != null){
				redirParams["has_hint_file"]="${predictionInstance.warn}"
			}
			if(predictionInstance.species_select != "null"){
				redirParams["has_select"]="${predictionInstance.warn}"
			}
			if(predictionInstance.utr == true){
				redirParams["has_utr"]="${predictionInstance.warn}"	
			}
			if(predictionInstance.pred_strand != 1){
				redirParams["has_strand"]="${predictionInstance.warn}"
			}
			if(predictionInstance.alt_transcripts != 1){
				redirParams["has_transcripts"]="${predictionInstance.warn}"
			}
			if(predictionInstance.allowed_structures != 1){
				redirParams["has_structures"]="${predictionInstance.warn}"
			}
			if(predictionInstance.ignore_conflicts == true){
				redirParams["has_conflicts"]="${predictionInstance.warn}"
			}
			redirParams["warn"]="${predictionInstance.warn}"
			// put redirect procedure into a function
			def cleanRedirect = {
				logDate = new Date()
				if(predictionInstance.email_adress == null){
           				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${predictionInstance.accession_id} by anonymous user with IP ${userIP} is aborted!\n"
				}else{
           				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${predictionInstance.accession_id} by user ${predictionInstance.email_adress} is aborted!\n"
				}
				flash.message = "Info: Please check all fields marked in blue for completeness before starting the prediction job!"
            			redirect(action:create, params:redirParams)
			}
			// clean up directory (delete) function
			def String dirName = "${output_dir}/${predictionInstance.accession_id}"
			def projectDir = new File(dirName)
			def deleteDir = {
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - Project directory is deleted\n"
				cmdStr = "rm -r ${projectDir} &> /dev/null"
				delProc = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"								}
            			delProc.waitFor()
			}
			// log abort function
			def logAbort = {
				logDate = new Date()
					if(predictionInstance.email_adress == null){
           					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${predictionInstance.accession_id} by anonymous user with IP ${userIP} is aborted!\n"
					}else{
           					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${predictionInstance.accession_id} by user ${predictionInstance.email_adress} is aborted!\n"
					}
			}
			//verify that the submitter is a person
			boolean captchaValid = simpleCaptchaService.validateCaptcha(params.captcha)
			if(captchaValid == false){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - The user is probably not a human person.\n"
				flash.error = "The verification string at the bottom of the page was not entered correctly!"
            			cleanRedirect()
           			return
			}
			// utr checkbox
			if(predictionInstance.utr == true){
				overRideUtrFlag = 1
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - User enabled UTR prediction.\n"
			}else{
				overRideUtrFlag = 0
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - User did not enable UTR prediction.\n"
			}
			// get parameter archive file (if available)
			//def uploadedParamArch = request.getFile('ArchiveFile')
			if(!uploadedParamArch.empty){
				// check file size
				preUploadSize = uploadedParamArch.getSize()
				if(preUploadSize <= maxButtonFileSize){
					// actually upload the file
					projectDir.mkdirs()
         				uploadedParamArch.transferTo( new File (projectDir, "parameters.tar.gz"))
					//predictionInstance.archive_file = uploadedParamArch.originalFilename
					confirmationString = "${confirmationString}Parameter archive: ${predictionInstance.archive_file}\n"
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - uploaded parameter archive ${predictionInstance.archive_file} was renamed to parameters.tar.gz and moved to ${projectDir}\n"
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The selected parameter archive file was bigger than ${maxButtonFileSize}.\n"
					flash.error = "Parameter archive file is bigger than ${maxButtonFileSize} bytes, which is our maximal size for file upload from local harddrives via web browser. Please select a smaller file or use the ftp/http web link file upload option."
            				cleanRedirect()
					return
				}
				// get cksum and file size for database
				def archCksumScript = new File("${projectDir}/archive_cksum.sh")
         			def archCksumFile = "${projectDir}/arch.cksum"
				cmd2Script = "cksum ${projectDir}/parameters.tar.gz > ${archCksumFile} 2> /dev/null"
         			archCksumScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v3 - archCksumScript << \"${cmd2Script}\"\n"
				}
				cmdStr = "bash ${projectDir}/archive_cksum.sh"
         			def archCksumProcess = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
         			archCksumProcess.waitFor()
         			def archCksumContent = new File("${archCksumFile}").text
         			def archCksum_array = archCksumContent =~/(\d*) \d* /
         			def archCksum
         			(1..archCksum_array.groupCount()).each{archCksum = "${archCksum_array[0][it]}"}
         			predictionInstance.archive_cksum = "${archCksum}"
         			predictionInstance.archive_size = uploadedParamArch.size
				logDate = new Date()
         			logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - parameters.tar.gz is ${predictionInstance.archive_size} big and has a cksum of ${archCksum}.\n"
				cmdStr = "rm ${projectDir}/arch.cksum"
         			def delProcCksumarch = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
         			delProcCksumarch.waitFor()
				cmdStr = "rm ${projectDir}/archive_cksum.sh &> /dev/null"
         			def delProcCkSharch = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
         			delProcCkSharch.waitFor()
				// check whether the archive contains all relevant files
				def String paramDirName = "${projectDir}/params"
				def paramDir = new File(paramDirName)
				paramDir.mkdirs()
				def checkParamArch = new File("${projectDir}/ckArch.sh")
				cmd2Script = "${AUGUSTUS_SCRIPTS_PATH}/checkParamArchive.pl ${projectDir}/parameters.tar.gz ${paramDirName} > ${projectDir}/archCheck.log 2> ${projectDir}/archCheck.err"
				checkParamArch << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v3 - checkParamArch << \"${cmd2Script}\"\n"
				}
				cmdStr = "bash ${checkParamArch}"
				def checkParamArchRunning = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				checkParamArchRunning.waitFor()
				def archCheckLog = new File("${projectDir}/archCheck.log")
				def archCheckErr = new File("${projectDir}/archCheck.err")
				def archCheckLogSize = archCheckLog.text.size()
				def archCheckErrSize = archCheckErr.text.size()
				// if essential file are missing, redirect to input interface and inform user that the archive was not compatible
				if(archCheckErrSize > 0){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The parameter archive was not compatible.\n"
					deleteDir()
           				flash.error = "Parameter archive ${uploadedParamArch.originalFilename} is not compatible with the AUGUSTUS prediction web server application."
            				cleanRedirect()
           				return
				// if only UTR params are missing, set flag to override any user-defined UTR settings
				}else if(archCheckLogSize > 0){
					overRideUtrFlag = 0 // UTR predictions are now permanently disabled
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - UTR predictions have been disabled because UTR parameters are missing!\n"
				}
				archiveExistsFlag = 1
			}else{predictionInstance.archive_file = "empty"}
			// check whether parameters are available for project_id (previous prediction run)
			logDate = new Date()
			logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The given parameter ID is ${predictionInstance.project_id}\n"
			if(!(predictionInstance.project_id == null)){
				def spec_conf_dir = new File("${AUGUSTUS_CONFIG_PATH}/species/${predictionInstance.project_id}")
				if(!spec_conf_dir.exists()){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The given parameter-string does not exist on our system.\n"
					deleteDir()
           				flash.error = "The specified parameter ID ${predictionInstance.project_id} does not exist on our system."
            				cleanRedirect()
           				return
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Requested ${spec_conf_dir} exists on our system.\n"
					speciesNameExistsFlag = 1
					species = predictionInstance.project_id
				}
				confirmationString = "${confirmationString}AUGUSTUS parameter project identifier: ${predictionInstance.project_id}\n"
			}
			// check whether parameters were supplied in double or triple
			if(predictionInstance.archive_file == "empty" && predictionInstance.project_id == null && predictionInstance.species_select == "null"){
           			flash.error = "You need to specify a parameter archive for upload OR enter a project identifier OR select an organism!"
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - user specified no AUGUSTUS parameters\n"
				deleteDir()
            			cleanRedirect()
           			return
			}else if(predictionInstance.archive_file != "empty" && predictionInstance.project_id != null && predictionInstance.species_select != "null"){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - user specified more than one option for AUGUSTUS parameters\n"
				flash.error = "You specified parameters in three different ways. Please decide for on way! You need to specify a parameter archive for upload OR enter a project identifier OR select an organism!"
				deleteDir()
            			cleanRedirect()
           			return
			}else if(predictionInstance.archive_file != "empty" && predictionInstance.project_id != null){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - user specified more than one option for AUGUSTUS parameters\n"
				flash.error = "You specified parameters as archive file and as project ID. Please decide for on way! You need to specify a parameter archive for upload OR enter a project identifier OR select an organism!"
				deleteDir()
            			cleanRedirect()
				return
			}else if(predictionInstance.project_id != null && predictionInstance.species_select != "null"){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - user specified more than one option for AUGUSTUS parameters\n"
				flash.error = "You specified parameters as project ID and by selecting an organism from the dropdown menu. Please decide for on way! You need to specify a parameter archive for upload OR enter a project identifier OR select an organism!"
				deleteDir()
            			cleanRedirect()
				return
			}else if(predictionInstance.archive_file != "empty" && predictionInstance.species_select != "null"){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - user specified more than one option for AUGUSTUS parameters\n"
				flash.error = "You specified parameters as parameter archive and by selecting an organism from the dropdown menu. Please decide for on way! You need to specify a parameter archive for upload OR enter a project identifier OR select an organism!"
				deleteDir()
            			cleanRedirect()
				return
			}
			// assign parameter set from dropdown menu
			if(predictionInstance.species_select == "Acyrthosiphon pisum (animal)"){
				predictionInstance.project_id = "pea_aphid"				
			}else if(predictionInstance.species_select == "Aedes aegypti (animal)"){
				predictionInstance.project_id = "aedes"
			}else if(predictionInstance.species_select == "Amphimedon queenslandica (animal)"){
				predictionInstance.project_id = "amphimedon"
			}else if(predictionInstance.species_select == "Apis mellifera (animal)"){
				predictionInstance.project_id = "honeybee1"
			}else if(predictionInstance.species_select == "Brugia malayi (animal)"){
				predictionInstance.project_id = "brugia"
			}else if(predictionInstance.species_select == "Caenorhabditis elegans (animal)"){
				predictionInstance.project_id = "caenorhabditis"
			}else if(predictionInstance.species_select == "Callorhinchus milii (animal)"){
				predictionInstance.project_id = "elephant_shark"
			}else if(predictionInstance.species_select == "Drosophila melanogaster (animal)"){
				predictionInstance.project_id = "fly"
			}else if(predictionInstance.species_select == "Gallus gallus domesticus (animal)"){
			      	predictionInstance.project_id = "chicken"
			}else if(predictionInstance.species_select == "Homo sapiens (animal)"){
				predictionInstance.project_id = "human"
			}else if(predictionInstance.species_select == "Petromyzon marinus (animal)"){
				predictionInstance.project_id = "lamprey"
			}else if(predictionInstance.species_select == "Nasonia vitripennis (animal)"){
				predictionInstance.project_id = "nasonia"
			}else if(predictionInstance.species_select == "Schistosoma mansoni (animal)"){
				predictionInstance.project_id = "schistosoma"
			}else if(predictionInstance.species_select == "Tribolium castaneum (animal)"){
				predictionInstance.project_id = "tribolium"
			}else if(predictionInstance.species_select == "Trichinella spiralis (animal)"){
				predictionInstance.project_id = "trichinella"
			}else if(predictionInstance.species_select == "Tetrahymena thermophila (alveolata)"){
				predictionInstance.project_id = "tetrahymena"
			}else if(predictionInstance.species_select == "Toxoplasma gondii (alveolata)"){
				predictionInstance.project_id = "toxoplasma"
			}else if(predictionInstance.species_select == "Leishmania tarantolae (protozoa)"){
				predictionInstance.project_id = "leishmania_tarentolae"
			}else if(predictionInstance.species_select == "Arabidopsis thaliana (plant)"){
				predictionInstance.project_id = "arabidopsis"
			}else if(predictionInstance.species_select == "Chlamydomonas reinhardtii (alga)"){
				predictionInstance.project_id = "chlamy2011"
			}else if(predictionInstance.species_select == "Galdieria sulphuraria (alga)"){
				predictionInstance.project_id = "galdieria"
			}else if(predictionInstance.species_select == "Solaneum lycopersicum (plant)"){
				predictionInstance.project_id = "tomato"
			}else if(predictionInstance.species_select == "Triticum/wheat (plant)"){
			      	predictionInstance.project_id = "wheat"
			}else if(predictionInstance.species_select == "Zea mays (plant)"){
				predictionInstance.project_id = "maize"
			}else if(predictionInstance.species_select == "Aspergillus fumigatus (fungus)"){
				predictionInstance.project_id = "aspergillus_fumigatus"
			}else if(predictionInstance.species_select == "Aspergillus nidulans (fungus)"){
				predictionInstance.project_id = "aspergillus_nidulans"
			}else if(predictionInstance.species_select == "Aspergillus oryzae (fungus)"){
				predictionInstance.project_id = "aspergillus_oryzae"
			}else if(predictionInstance.species_select == "Aspergillus terreus (fungus)"){
				predictionInstance.project_id = "aspergillus_terreus"
			}else if(predictionInstance.species_select == "Botrytis cinerea (fungus)"){
				predictionInstance.project_id = "botrytis_cinerea"
			}else if(predictionInstance.species_select == "Candida albicans (fungus)"){
				predictionInstance.project_id = "candida_albicans"
			}else if(predictionInstance.species_select == "Candida guilliermondii (fungus)"){
				predictionInstance.project_id = "candida_guilliermondii"
			}else if(predictionInstance.species_select == "Candida tropicalis (fungus)"){
				predictionInstance.project_id = "candida_tropicalis"
			}else if(predictionInstance.species_select == "Chaetomium globosum (fungus)"){
				predictionInstance.project_id = "chaetomium_globosum"
			}else if(predictionInstance.species_select == "Coccidioides immitis (fungus)"){
				predictionInstance.project_id = "coccidioides_immitis"
			}else if(predictionInstance.species_select == "Coprinus cinereus (fungus)"){
				predictionInstance.project_id = "coprinus"
			}else if(predictionInstance.species_select == "Cryptococcus neoformans (fungus)"){
				predictionInstance.project_id = "cryptococcus_neoformans_neoformans_B"
			}else if(predictionInstance.species_select == "Debarymomyces hansenii (fungus)"){
				predictionInstance.project_id = "debaryomyces_hansenii"
			}else if(predictionInstance.species_select == "Encephalitozoon cuniculi (fungus)"){
				predictionInstance.project_id = "encephalitozoon_cuniculi_GB"
			}else if(predictionInstance.species_select == "Eremothecium gossypii (fungus)"){
				predictionInstance.project_id = "eremothecium_gossypii"
			}else if(predictionInstance.species_select == "Fusarium graminearum (fungus)"){
				predictionInstance.project_id = "fusarium_graminearum"
			}else if(predictionInstance.species_select == "Histoplasma capsulatum (fungus)"){
				predictionInstance.project_id = "histoplasma_capsulatum"
			}else if(predictionInstance.species_select == "Kluyveromyces lactis (fungus)"){
				predictionInstance.project_id = "kluyveromyces_lactis"
			}else if(predictionInstance.species_select == "Laccaria bicolor (fungus)"){
				predictionInstance.project_id = "laccaria_bicolor"
			}else if(predictionInstance.species_select == "Lodderomyces elongisporus (fungus)"){
				predictionInstance.project_id = "lodderomyces_elongisporus"
			}else if(predictionInstance.species_select == "Magnaporthe grisea (fungus)"){
				predictionInstance.project_id = "magnaporthe_grisea"
			}else if(predictionInstance.species_select == "Neurospora crassa (fungus)"){
				predictionInstance.project_id = "neurospora_crassa"
			}else if(predictionInstance.species_select == "Phanerochaete chrysosporium (fungus)"){
				predictionInstance.project_id = "phanerochaete_chrysosporium"
			}else if(predictionInstance.species_select == "Pichia stipitis (fungus)"){
				predictionInstance.project_id = "pichia_stipitis"
			}else if(predictionInstance.species_select == "Rhizopus oryzae (fungus)"){
				predictionInstance.project_id = "rhizopus_oryzae"
			}else if(predictionInstance.species_select == "Saccharomyces cerevisiae (fungus)"){
				predictionInstance.project_id = "saccharomyces_cerevisiae_S288C"
			}else if(predictionInstance.species_select == "Camponotus floridanus (animal)"){
			      	predictionInstance.project_id = "camponotus_floridanus"
                        }else if(predictionInstance.species_select == "Danio rerio (animal)"){
			        predictionInstance.project_id = "zebrafish"
			}else if(predictionInstance.species_select == "Schizosaccharomyces pombe (fungus)"){
				predictionInstance.project_id = "schizosaccharomyces_pombe"
			}else if(predictionInstance.species_select == "Ustilago maydis (fungus)"){
				predictionInstance.project_id = "ustilago_maydis"
			}else if(predictionInstance.species_select == "Verticillium longisporum (fungus)"){
				predictionInstance.project_id = "verticillium_longisporum1"
			}else if(predictionInstance.species_select == "Yarrowia lipolytica (fungus)"){
				predictionInstance.project_id = "yarrowia_lipolytica"
			}else if(predictionInstance.species_select == "Heliconius melpomene (animal)"){
				predictionInstance.project_id = "heliconius_melpomene1"
			}else if(predictionInstance.species_select == "Bombus terrestris (animal)"){
				predictionInstance.project_id = "bombus_terrestris2"
			}else if(predictionInstance.species_select == "Rhodnius prolixus (animal)"){
				predictionInstance.project_id = "rhodnium"
			}else if(predictionInstance.species_select == "Conidiobolus coronatus (fungus)"){
				predictionInstance.project_id = "Conidiobolus_coronatus"
			}else if(predictionInstance.species_select == "Sulfolobus solfataricus (archaeon)"){
			        predictionInstance.project_id = "sulfolobus_solfataricus"
				prokaryotic = true
			}else if(predictionInstance.species_select == "Escherichia coli (bacterium)"){
			        predictionInstance.project_id = "E_coli_K12"
				prokaryotic = true
			}else if(predictionInstance.species_select == "Thermoanaerobacter tengcongensis (bacterium)"){
                                predictionInstance.project_id = "thermoanaerobacter_tengcongensis"
				prokaryotic = true
                        }
			if(predictionInstance.project_id != null && predictionInstance.species_select != "null"){
				species = predictionInstance.project_id
				confirmationString = "${confirmationString}AUGUSTUS parameter project identifier: ${predictionInstance.project_id}\n"
			}
			logDate = new Date()
			logFile << "${logDate} ${predictionInstance.accession_id} v1 - Parameter set ${predictionInstance.project_id} was assigned through dropdown selection ${predictionInstance.species_select}\n"
			if(predictionInstance.project_id == null && predictionInstance.archive_file == "empty"){
				logDate = new Date()
				logFile << "${logDate} ${predictionInstance.accession_id} v1 - project_id is empty.\n"
				deleteDir()
				flash.error = "No parameters given!"
            			cleanRedirect()
				return

			}

			// upload of genome file
			//def uploadedGenomeFile
			//uploadedGenomeFile = request.getFile('GenomeFile')
     			def seqNames = []
			if(!uploadedGenomeFile.empty){
				// check file size
				preUploadSize = uploadedGenomeFile.getSize()
         			projectDir.mkdirs()
				if(preUploadSize <= maxButtonFileSize){
         				uploadedGenomeFile.transferTo( new File (projectDir, "genome.fa"))
        				//predictionInstance.genome_file = uploadedGenomeFile.originalFilename
					confirmationString = "${confirmationString}Genome file: ${predictionInstance.genome_file}\n"
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The selected genome file was bigger than ${maxButtonFileSize}. Submission rejected.\n"
					flash.error = "Genome file is bigger than ${maxButtonFileSize} bytes, which is our maximal size for file upload from local harddrives via web browser. Please select a smaller file or use the ftp/http web link file upload option."
					return
				}
				if("${uploadedGenomeFile.originalFilename}" =~ /\.gz/){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Genome file is gzipped.\n"
					def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
					cmd2Script = "cd ${projectDir}; mv genome.fa genome.fa.gz &> /dev/null; gunzip genome.fa.gz 2> /dev/null"
					gunzipGenomeScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v3 - gunzipGenomeScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${gunzipGenomeScript}"
					def gunzipGenome = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					gunzipGenome.waitFor()
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Unpacked genome file.\n"
					delProc = "rm ${gunzipGenomeScript} &> /dev/null".execute()
					delProc.waitFor()
				}
				logDate = new Date()
         			logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - uploaded genome file ${uploadedGenomeFile.originalFilename} was renamed to genome.fa and moved to ${projectDir}\n"
				// check number of scaffolds
                                def nSeqFile = new File("${projectDir}/genome_nSeq.sh")
                                cmd2Script = 'grep -c ">" ${projectDir}/genome.fa > ${projectDir}/genome.nSeq'
                                nSeqFile << "${cmd2Script}"
                                def nSeqStatus = "${cmdStr}".execute()
                                nSeqStatus.waitFor()
                                def nSeqResult = new File("${projectDir}/genome.nSeq").text
                                def nSeq_array = nSeqResult =~ /(\d*)/
                                def nSeqNumber
                                (1..nSeq_array.groupCount()).each{nSeqNumber = "${nSeq_array[0][it]}"}
                                cmdStr = "rm ${nSeqFile} ${projectDir}/genome.nSeq &> /dev/null"
                                delProc = "${cmdStr}".execute();
                                delProc.waitFor()
                                if(nSeqNumber > maxNSeqs){
                                       logDate = new Date()
                                       logFile << "${logDate} ${predictionInstance.accession_id} v1 - genome file contains more than ${maxNSeqs} scaffolds. Aborting job."
                                       flash.error = "Genome file contains more than ${maxNSeqs} scaffolds, which is the maximal number of scaffolds that we permit for submission with WebAUGUSTUS. Please remove all short scaffolds from your genome file."
                                       cleanRedirect()
                                       return
                                }

        			// check for fasta format & extract fasta headers for gff validation:
         			new File("${projectDir}/genome.fa").eachLine{line -> 
            			if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }
					if(line =~ /\*/ || line =~ /\?/){
						metacharacterFlag = 1
					}else{
            					if(line =~ /^>/){
               						def len = line.length()
               						seqNames << line[1..(len-1)]
            					}
					}	
         			}
				if(metacharacterFlag == 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The genome file contains metacharacters (e.g. * or ?).\n";
					deleteDir()
          				flash.error = "Genome file contains metacharacters (*, ?, ...). This is not allowed."
            				cleanRedirect()
					return
				}	
         			if(genomeFastaFlag == 1) {
					logDate = new Date()
            				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The genome file was not fasta.\n"
					deleteDir()
            				flash.error = "Genome file ${uploadedGenomeFile.originalFilename} is not in DNA fasta format."
            				cleanRedirect()
					return
	         		} else {
	            			def genomeCksumScript = new File("${projectDir}/genome_cksum.sh")
	            			def genomeCksumFile = "${projectDir}/genome.cksum"
					cmd2Script = "cksum ${projectDir}/genome.fa > ${genomeCksumFile} 2> /dev/null"
	            			genomeCksumScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v3 - genomeCksumScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${projectDir}/genome_cksum.sh"
	            			def genomeCksumProcess = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
	            			genomeCksumProcess.waitFor()
	            			def genomeCksumContent = new File("${genomeCksumFile}").text
	            			def genomeCksum_array = genomeCksumContent =~/(\d*) \d* /
	            			def genomeCksum
	            			(1..genomeCksum_array.groupCount()).each{genomeCksum = "${genomeCksum_array[0][it]}"}
	            			predictionInstance.genome_cksum = "${genomeCksum}"
	            			predictionInstance.genome_size = uploadedGenomeFile.size
					logDate = new Date()
	            			logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - genome.fa is ${predictionInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
					cmdStr = "rm ${projectDir}/genome.cksum &> /dev/null"
	            			def delProcCksumGenome = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
	            			delProcCksumGenome.waitFor()
					cmdStr = "rm ${projectDir}/genome_cksum.sh &> /dev/null"
	            			def delProcCkShGenome = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
	            			delProcCkShGenome.waitFor()
	         		}
			}
			// retrieve beginning of genome file for format check
	      		if(!(predictionInstance.genome_ftp_link == null)){
				confirmationString = "${confirmationString}Genome file: ${predictionInstance.genome_ftp_link}\n"
				logDate = new Date()
	         		logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - genome web-link is ${predictionInstance.genome_ftp_link}\n"
	         		projectDir.mkdirs()
	         		def URL url = new URL("${predictionInstance.genome_ftp_link}");
				// check whether URL exists
				urlExistsScript = new File("${projectDir}/genomeExists.sh")
				cmd2Script = "curl -o /dev/null --silent --head --write-out '%{http_code}\n' \"${predictionInstance.genome_ftp_link}\" > ${projectDir}/genomeExists"
				urlExistsScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v3 - urlExistsScript << \"${cmd2Script}\"\n"
				}
				cmdStr = "bash ${urlExistsScript}"
				def genomeUrlExists = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				genomeUrlExists.waitFor()
				cmdStr = "rm ${urlExistsScript} &> /dev/null"			
				delProc = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				delProc.waitFor()
				content = new File("${projectDir}/genomeExists").text
				st = new Scanner(content)//works for exactly one number in a file
				error_code = st.nextInt();
				if(!(error_code == 200)){
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v1 - The genome URL is not accessible. Response code: ${error_code}.\n"
					deleteDir()
					flash.error = "Cannot retrieve genome file from HTTP/FTP link ${predictionInstance.genome_ftp_link}."
            				cleanRedirect()
					return
				}else{
					logDate = new Date()
					logFile << "${logDate} ${predictionInstance.accession_id} v1 - The genome URL is accessible. Response code: ${error_code}.\n"
				}
	         		// checking web file for DNA fasta format: 
	         		def URLConnection uc = url .openConnection()
				if(!("${predictionInstance.genome_ftp_link}" =~ /\.gz/)){
	         			def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
					try{
	         				def String inputLine=null
                                                def char inputChar=null
                                                def charCounter = 1
                                                while ( ((inputChar = br.read()) != null) && (charCounter <= 1000)) {
                                                        if(inputChar =~ />/){
                                                                     inputLine = br.readLine();
                                                        }else if(!(inputChar =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(inputChar =~ /^$/)){
                                                        genomeFastaFlag = 1
                                                        }
                                                        charCounter = charCounter + 1
                                                }
					}finally{
	         				br.close()
					}
	         			if(genomeFastaFlag == 1) {
						logDate = new Date()
	            				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The first 20 lines in genome file are not fasta.\n"
						deleteDir()
	            				flash.error = "Genome file ${predictionInstance.genome_ftp_link} is not in DNA fasta format."
            					cleanRedirect()
	            				return
	         			}
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The linked genome file is gzipped. Format will be checked later after extraction.\n"
				}
	      		}
	      		// upload of est file
	      		// def uploadedEstFile = request.getFile('EstFile')
	      		if(!uploadedEstFile.empty){
				// check file size
				preUploadSize = uploadedEstFile.getSize()
				if(preUploadSize <= maxButtonFileSize){
	         			projectDir.mkdirs()
	         			uploadedEstFile.transferTo( new File (projectDir, "est.fa"))
	         			//predictionInstance.est_file = uploadedEstFile.originalFilename
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The selected cDNA file was bigger than ${maxButtonFileSize}.\n"
					flash.error = "cDNA file is bigger than ${maxButtonFileSize} bytes, which is our maximal size for file upload from local harddrives via web browser. Please select a smaller file or use the ftp/http web link file upload option."
					deleteDir()
            				cleanRedirect()
					return
				}
				confirmationString = "${confirmationString}cDNA file: ${predictionInstance.est_file}\n"
				if("${uploadedEstFile.originalFilename}" =~ /\.gz/){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - EST file is gzipped.\n"
					def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
					cmd2Script = "cd ${projectDir}; mv est.fa est.fa.gz &> /dev/null; gunzip est.fa.gz 2> /dev/null"
					gunzipEstScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - gunzipEstScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${gunzipEstScript}"
					def gunzipEst = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					gunzipEst.waitFor()
					logFile <<  "${predictionInstance.accession_id} v1 - Unpacked EST file.\n"
					cmdStr = "rm ${gunzipEstScript} &> /dev/null"
					delProc = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProc.waitFor()
				}
	         		logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Uploaded EST file ${uploadedEstFile.originalFilename} was renamed to est.fa and moved to ${projectDir}\n"
	         		// check fasta format
	         		new File("${projectDir}/est.fa").eachLine{line -> 
					if(line =~ /\*/ || line =~ /\?/){
						metacharacterFlag = 1
					}else{
	            				if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNnUu]/) && !(line =~ /^$/)){ 
							estFastaFlag = 1
						}
					}
	         		}
				if(metacharacterFlag == 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The cDNA file contains metacharacters (e.g. * or ?).\n"
					deleteDir()
            				flash.error = "cDNA file contains metacharacters (*, ?, ...). This is not allowed."
            				cleanRedirect()
					return
				}	
	         		if(estFastaFlag == 1) {
            				logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The cDNA file was not fasta.\n"
					deleteDir()
            				flash.error = "cDNA file ${uploadedEstFile.originalFilename} is not in DNA fasta format."
            				cleanRedirect()
            				return
         			} else { estExistsFlag = 1 }
         				def estCksumScript = new File("${projectDir}/est_cksum.sh")
         				def estCksumFile = "${projectDir}/est.cksum"
					cmd2Script = "cksum ${projectDir}/est.fa > ${estCksumFile} 2> /dev/null"
         				estCksumScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - estCksumScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${projectDir}/est_cksum.sh"
         				def estCksumProcess = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
         				estCksumProcess.waitFor()
         				def estCksumContent = new File("${estCksumFile}").text
         				def estCksum_array = estCksumContent =~/(\d*) \d* /
         				def estCksum
         				(1..estCksum_array.groupCount()).each{estCksum = "${estCksum_array[0][it]}"}
         				predictionInstance.est_cksum = "${estCksum}"
         				predictionInstance.est_size = uploadedEstFile.size
         				logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - est.fa is ${predictionInstance.est_size} big and has a cksum of ${estCksum}.\n"
					cmdStr = "rm ${projectDir}/est.cksum &> /dev/null"
         				def delProcCksumEst = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
         				delProcCksumEst.waitFor()
					cmdStr = "rm ${projectDir}/est_cksum.sh &> /dev/null"
         				def delProcCkShEst = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
         				delProcCkShEst.waitFor()
      			}
      			// retrieve beginning of est file for format check
      			if(!(predictionInstance.est_ftp_link == null)){
				confirmationString = "${confirmationString}cDNA file: ${predictionInstance.est_ftp_link}\n"
         			logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - est web-link is ${predictionInstance.est_ftp_link}\n"
         			projectDir.mkdirs()
				estExistsFlag = 1
				if(!("${predictionInstance.est_ftp_link}" =~ /\.gz/)){
         				def URL url = new URL("${predictionInstance.est_ftp_link}");
					// check whether URL exists
					urlExistsScript = new File("${projectDir}/estExists.sh")
					cmd2Script = "curl -o /dev/null --silent --head --write-out '%{http_code}\n' \"${predictionInstance.est_ftp_link}\" > ${projectDir}/estExists"
					urlExistsScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - urlExistsScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${urlExistsScript}"
					def genomeUrlExists = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					genomeUrlExists.waitFor()
					cmdStr = "rm ${urlExistsScript} &> /dev/null"			
					delProc = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProc.waitFor()
					content = new File("${projectDir}/estExists").text
					st = new Scanner(content)//works for exactly one number in a file
					error_code = st.nextInt();
					if(!(error_code == 200)){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The EST URL is not accessible. Response code: ${error_code}.\n"
						deleteDir()
						flash.error = "Cannot retrieve cDNA file from HTTP/FTP link ${predictionInstance.est_ftp_link}."
            					cleanRedirect()
						return
					}else{
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The EST URL is accessible. Response code: ${error_code}.\n"
					}
         				// checking web file for DNA fasta format: 
         				def URLConnection uc = url .openConnection()
         				def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
					try{
         					def String inputLine=null
                                                def char inputChar=null
                                                def charCounter = 1
                                                while ( ((inputChar = br.read()) != null) && (charCounter <= 1000)) {
                                                        if(inputChar =~ />/){
                                                                     inputLine = br.readLine();
                                                        }else if(!(inputChar =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(inputChar =~ /^$/)){
                                                        estFastaFlag = 1
                                                        }
                                                        charCounter = charCounter + 1
                                                }
					}finally{
         					br.close()
					}
         				if(estFastaFlag == 1) {
           					logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The cDNA file was not fasta.\n"
						deleteDir()
            					flash.error = "cDNA file ${predictionInstance.est_ftp_link} is not in DNA fasta format."
            					cleanRedirect()
            					return
         				}
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The linked EST file is gzipped. Format will be checked later after extraction.\n"
				}
      			}
			// get hints file, format check
			// def uploadedStructFile = request.getFile('HintFile')
			if(!uploadedStructFile.empty){
				// check file size
				preUploadSize = uploadedStructFile.getSize()
				if(preUploadSize <= maxButtonFileSize){
					projectDir.mkdirs()
					uploadedStructFile.transferTo( new File (projectDir, "hints.gff"))
					//predictionInstance.hint_file = uploadedStructFile.originalFilename
				}else{
					def long allowedHintsSize = maxButtonFileSize * 2
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The selected Hints file was bigger than ${allowedHintsSize}.\n"
					flash.error = "Hints file is bigger than ${allowedHintsSize} bytes, which is our maximal size for file upload from local harddrives via web browser. Please select a smaller file or use the ftp/http web link file upload option."
            				cleanRedirect()
					return
				}
				confirmationString = "${confirmationString}Hints file: ${predictionInstance.hint_file}\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Uploaded hints file ${uploadedStructFile.originalFilename} was renamed to hints.gff and moved to ${projectDir}\n"
				def gffColErrorFlag = 0
				def gffNameErrorFlag = 0
				def gffSourceErrorFlag = 0
				def gffFeatureErrorFlag = 0
				if(!uploadedGenomeFile.empty){ // if seqNames already exists
					// gff format validation: number of columns 9, + or - in column 7, column 1 muss member von seqNames sein
					def gffArray
					def isElement
					new File("${projectDir}/hints.gff").eachLine{line -> 
						if(line =~ /\*/ || line =~ /\?/){
							metacharacterFlag = 1
						}else{
							gffArray = line.split("\t")
							if(!(gffArray.size() == 9)){ 
								gffColErrorFlag = 1
							}else{
								isElement = 0
								seqNames.each{ seq ->
									if(seq =~ /${gffArray[0]}/){ isElement = 1 }
									if(isElement == 0){ gffNameErrorFlag = 1 }
									if(!("${gffArray[8]}" =~ /source=M/)){gffSourceErrorFlag = 1}
									if(!("${gffArray[2]}" =~ /start$/) && !("${gffArray[2]}" =~ /stop$/) && !("${gffArray[2]}" =~ /tss$/) && !("${gffArray[2]}" =~ /tts$/) && !("${gffArray[2]}" =~ /ass$/) && !("${gffArray[2]}" =~ /dss$/) && !("${gffArray[2]}" =~ /exonpart$/) && !("${gffArray[2]}" =~ /exon$/) && !("${gffArray[2]}" =~ /exon$/) && !("${gffArray[2]}" =~ /intronpart$/) && !("${gffArray[2]}" =~ /intron$/) && !("${gffArray[2]}" =~ /CDSpart$/) && !("${gffArray[2]}" =~ /CDS$/) && !("${gffArray[2]}" =~ /UTRpart$/) && !("${gffArray[2]}" =~ /UTR$/) && !("${gffArray[2]}" =~ /irpart$/) && !("${gffArray[2]}" =~ /nonexonpart$/) && !("${gffArray[2]}" =~ /genicpart$/)){
										gffFeatureErrorFlag = 1
									}
								}
							}
						}
					}
					if(metacharacterFlag == 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The hints file contains metacharacters (e.g. * or ?).\n"
						deleteDir()
            					flash.error = "Hints file contains metacharacters (*, ?, ...). This is not allowed."
            					cleanRedirect()
						return
					}
					if(gffSourceErrorFlag == 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hint file's last column is not in correct format\n"
						flash.error = "Hints file  ${predictionInstance.hint_file} is not in a compatible gff format (the last column does not contain source=M). Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffColErrorFlag == 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hint file does not always contain 9 columns.\n"
						flash.error = "Hints file  ${predictionInstance.hint_file} is not in a compatible gff format (has not 9 columns). Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffNameErrorFlag == 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hint file contains entries that do not comply with genome sequence names.\n"
						flash.error = "Entries in the hints file  ${predictionInstance.hint_file} do not match the sequence names of the genome file. Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffFeatureErrorFlag == 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hint file contains unsupported features.\n"
						flash.error = "Entries in the hints file  ${predictionInstance.hint_file} contain unsupported features. Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if((gffColErrorFlag == 1 || gffNameErrorFlag == 1 || gffSourceErrorFlag == 1 || gffFeatureErrorFlag == 1)){
						deleteDir()
            					cleanRedirect()
						return
					}
				}
				hintExistsFlag = 1
				def structCksumScript = new File("${projectDir}/struct_cksum.sh")
				def structCksumFile = "${projectDir}/struct.cksum"
				cmd2Script = "cksum ${projectDir}/hints.gff > ${structCksumFile} 2> /dev/null"
				structCksumScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - structCksumScript << \"${cmd2Script}\"\n"
				}
				cmdStr = "bash ${projectDir}/struct_cksum.sh"
				def structCksumProcess = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				structCksumProcess.waitFor()
				def structCksumContent = new File("${structCksumFile}").text
				def structCksum_array = structCksumContent =~/(\d*) \d* /
				def structCksum
				(1..structCksum_array.groupCount()).each{structCksum = "${structCksum_array[0][it]}"}
				predictionInstance.hint_cksum = "${structCksum}"
				predictionInstance.hint_size = uploadedStructFile.size
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - hints.gff is ${predictionInstance.hint_size} big and has a cksum of ${structCksum}.\n"
				cmdStr = "rm ${projectDir}/struct.cksum &> /dev/null"
				def delProcCksumStruct = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				delProcCksumStruct.waitFor()
				cmdStr = "rm ${projectDir}/struct_cksum.sh &> /dev/null"
				def delProcCkShStruct = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				delProcCkShStruct.waitFor()
			}
			def radioParameterString
			confirmationString = "${confirmationString}User set UTR prediction: ${predictionInstance.utr}\n"
			// utr
			// check whether utr parameters actually exist:
			def utrParamContent = new File("${AUGUSTUS_CONFIG_PATH}/species/${species}/${species}_utr_probs.pbl")
			if(utrParamContent.exists() == false){
			        overRideUtrFlag = 0;
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - UTR prediction was disabled because UTR parameters do not exist for this species!\n";
			}
			// enable or disable utr prediction in AUGUSTUS command
			if(overRideUtrFlag==1){
				if(predictionInstance.allowed_structures == 1 || predictionInstance.allowed_structures == 2){
					radioParameterString = " --UTR=on"
				}else{
					radioParameterString = " --UTR=off"
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - UTR prediction was disabled due to incompatibility with at least one or exactly one gene predcition\n";
					overRideUtrFlag = 0;
				}
			}else if(overRideUtrFlag==0 && predictionInstance.utr == true){
				confirmationString = "${confirmationString}Server set UTR prediction: false [UTR parameters missing or conflict with allowed gene structure!]\n"
				radioParameterString = " --UTR=off"
			}else{
				radioParameterString = " --UTR=off"
			}
			// strand prediction radio buttons
			if(predictionInstance.pred_strand == 1){
				radioParameterString = "${radioParameterString} --strand=both"
				confirmationString = "${confirmationString}Report genes on: both strands\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction on both strands.\n"
			}else if(predictionInstance.pred_strand == 2){
				confirmationString = "${confirmationString}Report genes on: forward strand only\n"
				radioParameterString = "${radioParameterString} --strand=forward"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction on forward strand, only.\n"
			}else{
				confirmationString = "${confirmationString}Report genes on: reverse strand only\n"
				radioParameterString = "${radioParameterString} --strand=backward"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction on reverse strand, only.\n"
			}
			// alternative transcript radio buttons
			if(predictionInstance.alt_transcripts == 1){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=false"
				confirmationString = "${confirmationString}Alternative transcripts: none\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User disabled prediction of alternative transcripts.\n"
			}else if(predictionInstance.alt_transcripts == 2){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 --maxtracks=2"
				confirmationString = "${confirmationString}Alternative transcripts: few\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction of few alternative transcripts.\n"
			}else if(predictionInstance.alt_transcripts == 3){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.4 --maxtracks=3"
				confirmationString = "${confirmationString}Alternative transcripts: medium\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction of medium alternative transcripts.\n"
			}else{
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=20"
				confirmationString = "${confirmationString}Alternative transcripts: many\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled prediction of many alternative transcripts.\n"
			}
			// gene structure radio buttons
			if(predictionInstance.allowed_structures == 1){
				radioParameterString = "${radioParameterString} --genemodel=partial"
				confirmationString = "${confirmationString}Allowed gene structure: predict any number of (possibly partial) genes\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled the prediction of any number of genes.\n"
			}else if(predictionInstance.allowed_structures == 2){
				radioParameterString = "${radioParameterString} --genemodel=complete"
				confirmationString = "${confirmationString}Allowed gene structure: only predict complete genes\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User disabled the prediction of incomplete genes.\n"
			}else if(predictionInstance.allowed_structures == 3){
				radioParameterString = "${radioParameterString} --genemodel=atleastone"
				confirmationString = "${confirmationString}Allowed gene structure: only predict complete genes - at least one\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User disabled the prediction of incomplete genes and insists on at least one predicted gene.\n"
			}else{
				radioParameterString = "${radioParameterString} --genemodel=exactlyone"
				confirmationString = "${confirmationString}Allowed gene structure: predict exactly one gene\n"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled the prediction of exactly one gene.\n"
			}
			// ignore gene structure conflicts with other strand checkbox
			if(predictionInstance.ignore_conflicts == false){
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User did not enable to ignore strand conflicts.\n"
			}else{
				radioParameterString = "${radioParameterString} --strand=both"
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User enabled to ignore strand conflicts.\n"
			}
			confirmationString = "${confirmationString}Ignore conflictes with other strand: ${predictionInstance.ignore_conflicts}\n"
			// prokaryotic predictions (log information only)
			if(prokaryotic == false){
				logDate = new Date()
                                logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User selected a eukaryotic parameter set.\n"
			}else{
                                logDate = new Date()
                                logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - User selected an experimental prokaryotic parameter set.\n";
			}
			// send confirmation email and redirect
			if(!predictionInstance.hasErrors() && predictionInstance.save()){
				// save new variables in database
				predictionInstance.message = ""
				predictionInstance.save()
				// generate empty results page
				def emptyPageScript = new File("${projectDir}/emptyPage.sh")
				cmd2Script = "${AUGUSTUS_SCRIPTS_PATH}/writeResultsPage.pl ${predictionInstance.accession_id} null ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} 0 &> /dev/null"
				emptyPageScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - emptyPageScript << \"${cmd2Script}\"\n"
				}
				cmdStr = "bash ${projectDir}/emptyPage.sh"
				def emptyPageExecution = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				emptyPageExecution.waitFor()
				predictionInstance.job_status = 0
				logDate = new Date()
				mailStr = "Details of your job:\n\n${confirmationString}\n"
				predictionInstance.message = "----------------------------------------\n${logDate} - Message:\n"
				predictionInstance.message = "${predictionInstance.message}----------------------------------------\n\n${mailStr}"
				if(predictionInstance.email_adress != null){
					msgStr = "Hello!\n\n"
					msgStr = "${msgStr}Thank you for submitting the AUGUSTUS gene prediction "
					msgStr = "${msgStr}job ${predictionInstance.accession_id}.\n\n"
					msgStr = "${msgStr}${mailStr}The status/results page of your job is "
					msgStr = "${msgStr}${war_url}prediction/show/${predictionInstance.id}.\n\n"
					msgStr = "${msgStr}You will be notified via email when the job has finished.\n\nBest regards,\n\n"
					msgStr = "${msgStr}the AUGUSTUS web server team"
					sendMail {
						to "${predictionInstance.email_adress}"
						subject "AUGUSTUS prediction job ${predictionInstance.accession_id}"
						body """${msgStr}${footer}"""
					}
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Confirmation e-mail sent.\n" 
				}else{
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Did not send confirmation e-mail because user stays anonymous, but everything is ok.\n"
				}
				redirect(action:show,id:predictionInstance.id)
			} else {
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - An error occurred in the predictionInstance (e.g. E-Mail missing, see domain restrictions).\n"
				deleteDir()
				logAbort()
				render(view:'create', model:[predictionInstance:predictionInstance])
				return
			}

			//---------------------  BACKGROUND PROCESS ----------------------------
			Thread.start{
				// retrieve genome file
				if(!(predictionInstance.genome_ftp_link == null)){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Checking genome file size with curl prior upload\n"
					projectDir.mkdirs()
					// check whether the genome file is small enough for upload
					def fileSizeScript = new File("${projectDir}/filzeSize.sh")
					cmd2Script = "curl -sI ${predictionInstance.genome_ftp_link} | grep Content-Length | cut -d ' ' -f 2 > ${projectDir}/genomeFileSize 2> /dev/null"
					fileSizeScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - fileSizeScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${fileSizeScript}"
					def retrieveFileSize = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					retrieveFileSize.waitFor()	
					cmdStr = "rm ${fileSizeScript} &> /dev/null"		
					def delSzCrProc = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delSzCrProc.waitFor()
					content = new File("${projectDir}/genomeFileSize").text
					st = new Scanner(content)//works for exactly one number in a file
					def long genome_size;
					genome_size = st.nextLong();
					if(genome_size < maxFileSizeByWget){//1 GB
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Retrieving genome file ${predictionInstance.genome_ftp_link}\n"
						def getGenomeScript = new File("${projectDir}/getGenome.sh")
						cmd2Script = "wget -O ${projectDir}/genome.fa ${predictionInstance.genome_ftp_link} &> /dev/null"
						getGenomeScript << "${cmd2Script}"
						if(verb > 2){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - getGenomeScript << \"${cmd2Script}\"\n"
						}
						cmdStr = "bash ${projectDir}/getGenome.sh"
						def wgetGenome = "${cmdStr}".execute()
						if(verb > 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
						}
						wgetGenome.waitFor()
						cmdStr = "rm ${projectDir}/getGenome.sh &> /dev/null"
						delProc = "${cmdStr}".execute()
						if(verb > 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
						}
						delProc.waitFor()
						if("${predictionInstance.genome_ftp_link}" =~ /\.gz/){
							def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
							cmd2Script = "cd ${projectDir}; mv genome.fa genome.fa.gz &> /dev/null; gunzip genome.fa.gz 2> /dev/null"
							gunzipGenomeScript << "${cmd2Script}"
							if(verb > 2){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - gunzipGenomeScript << \"${cmd2Script}\"\n"
							}
							cmdStr = "bash ${gunzipGenomeScript}"
							def gunzipGenome = "${cmdStr}".execute()
							if(verb > 1){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
							}
							gunzipGenome.waitFor()
							cmdStr = "rm ${gunzipGenomeScript} &> /dev/null"
							delProc = "${cmdStr}".execute()
							if(verb > 1){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
							}
							delProc.waitFor()
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Unpacked genome file.\n"
						}
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - genome file upload finished, file stored as genome.fa at ${projectDir}\n"
                                                // check number of scaffolds (to avoid Java heapspace error in the next step)
                                                def nSeqFile = new File("${projectDir}/genome_nSeq.sh")
                                                cmd2Script = 'grep -c ">" ${projectDir}/genome.fa > ${projectDir}/genome.nSeq'
                                                nSeqFile << "${cmd2Script}"
                                                def nSeqStatus = "${cmdStr}".execute()
                                                nSeqStatus.waitFor()
                                                def nSeqResult = new File("${projectDir}/genome.nSeq").text
                                                def nSeq_array = nSeqResult =~ /(\d*)/
                                                def nSeqNumber
                                                (1..nSeq_array.groupCount()).each{nSeqNumber = "${nSeq_array[0][it]}"}
                                                cmdStr = "rm ${nSeqFile} ${projectDir}/genome.nSeq &> /dev/null"
                                                delProc = "${cmdStr}".execute();
                                                delProc.waitFor()
                                                if(nSeqNumber > maxNSeqs){
                                                        logDate = new Date()
                                                        logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The genome file contains more than ${maxNSeqs} scaffolds. Aborting job.\n";
                                                        deleteDir()
                                                        logAbort()
                                                        mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} for species\n${predictionInstance.project_name} was aborted\nbecause the provided genome file\n${predictionInstance.genome_ftp_link}\ncontains more than ${maxNSeqs} scaffolds. This is not allowed.\n\n"
                                                        logDate = new Date()
                                                        predictionInstance.message = "${predictionInstance.message}-----------------------------"
                                                        predictionInstance.message = "${predictionInstance.message}-----------------\n${logDate}"
                                                        predictionInstance.message = "${predictionInstance.message} - Error Message:\n-----------"
                                                        predictionInstance.message = "${predictionInstance.message}-----------------------------"
                                                        predictionInstance.message = "------\n\n${mailStr}"
                                                        preidctionInstance = predictionInstance.merge()
                                                        predictionInstance.save()
                                                        if(predictionInstance.email_adress != null){
                                                                msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS webserver team"
                                                                sendMail {
                                                                        to "${predictionInstance.email_adress}"
                                                                        subject "Your AUGUSTUS training job ${predictionInstance.accession_id} was aborted"
                                                                        body """${msgStr}${footer}"""
                                                                }
                                                        }
                                                        predictionInstance.results_urls = null
                                                        predictionInstance.job_status = 5
                                                        predictionInstance = predictionInstance.merge()
                                                        predictionInstance.save()
                                                        return
                                                }

						// check for fasta format & get seq names for gff validation:
						new File("${projectDir}/genome.fa").eachLine{line -> 
							if(line =~ /\*/ || line =~ /\?/){
								metacharacterFlag = 1
							}else{
								if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }	
								if(line =~ /^>/){
									def len = line.length()
									seqNames << line[1..(len-1)]
								}
							}
						}
						if(metacharacterFlag == 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The genome file contains metacharacters (e.g. * or ?).\n";
							deleteDir()
							logDate = new Date()
							logAbort()
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided genome file\n${predictionInstance.genome_ftp_link}\ncontains metacharacters (e.g. * or ?). This is not allowed.\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------"
							predictionInstance.message = "${predictionInstance.message}------------------\n${logDate}"
							predictionInstance.message = "${predictionInstance.message} - Error Message:\n----------"
							predictionInstance.message = "${predictionInstance.message}------------------------------"
							predictionInstance.message = "------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
							// delete database entry
							//predictionInstance.delete()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
						if(genomeFastaFlag == 1) {
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The genome file was not fasta.\n"
							deleteDir()
							logAbort()
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided genome file ${predictionInstance.genome_ftp_link} was not in DNA fasta format.\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress == null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
								to "${predictionInstance.email_adress}"
								subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
								body """${msgStr}${footer}"""
								}
							}
							// delete database entry
							//predictionInstance.delete()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
					}else{// actions if remote file was bigger than allowed
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Genome file size exceeds permitted ${maxFileSizeByWget} bytes.\n"
						logAbort()
						mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the genome file size was\nwith ${genome_size} bigger than 1 GB. Please submitt a smaller genome size!\n\n"
						def errorStrMsg = "Hello!\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
						logDate = new Date()
						predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						if(predictionInstance.email_adress != null){
							sendMail {
							to "${predictionInstance.email_adress}"
							subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
							body """${errorStrMsg}${footer}"""
							}
						}
						deleteDir()
						predictionInstance.results_urls = null
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						return

					}
					// check gff format
					def gffColErrorFlag = 0
					def gffNameErrorFlag = 0
					def gffSourceErrorFlag = 0
					if((!uploadedStructFile.empty) &&(!(predictionInstance.genome_ftp_link == null))){ // if seqNames already exists
						// gff format validation: number of columns 9, + or - in column 7, column 1 muss member von seqNames sein
						def gffArray
						def isElement
						new File("${projectDir}/hints.gff").eachLine{line -> 
							if(line =~ /\*/ || line =~ /\?/){
								metacharacterFlag = 1
							}else{



								gffArray = line.split("\t")
								if(!(gffArray.size() == 9)){ 
									gffColErrorFlag = 1 
								}else{
									isElement = 0
									seqNames.each{ seq ->
										if(seq =~ /${gffArray[0]}/){ isElement = 1 }
										if(isElement == 0){ gffNameErrorFlag = 1 }
										if(!("${gffArray[8]}" =~ /source=M/)){gffSourceErrorFlag = 1}
									}
								}
							}
						}
						if(metacharacterFlag == 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The hints file contains metacharacters (e.g. * or ?).\n";
							deleteDir()
							logDate = new Date()
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided hints file\ncontains metacharacters (e.g. * or ?). This is not allowed.\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
							// delete database entry
							//predictionInstance.delete()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
						if(gffColErrorFlag == 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hints file does not always contain 9 columns.\n"
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided hints file\n${predictionInstance.hint_file}\ndid not contain 9 columns in each line. Please make sure the gff-format complies\nwith the instructions in our 'Help' section before submitting another job!\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
						}
						if(gffNameErrorFlag == 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hints file contains entries that do not comply with genome sequence names.\n"
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the sequence names in\nthe provided hints file\n${predictionInstance.hint_file}\ndid not comply with the sequence names in the supplied genome file\n${predictionInstance.genome_ftp_link}.\nPlease make sure the gff-format complies with the instructions in our 'Help' section\nbefore submitting another job!\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
						}
						if(gffSourceErrorFlag ==1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Hints file contains entries that do not have source=M in the last column.\n"
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the last column of your\nhints file\n${predictionInstance.hint_file}\ndoes not contain the content source=M. Please make sure the gff-format complies with\nthe instructions in our 'Help' section before submitting another job!\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}						
						}
						if((gffColErrorFlag == 1 || gffNameErrorFlag == 1 || gffSourceErrorFlag ==1)){
							deleteDir()
							logAbort()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
					}
					def genomeCksumScript = new File("${projectDir}/genome_cksum.sh")
					def genomeCksumFile = "${projectDir}/genome.cksum"
					cmd2Script = "cksum ${projectDir}/genome.fa > ${genomeCksumFile} 2> /dev/null"
					genomeCksumScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${projectDir}/genome_cksum.sh"
					def genomeCksumProcess = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					genomeCksumProcess.waitFor()
					def genomeCksumContent = new File("${genomeCksumFile}").text
					def genomeCksum_array = genomeCksumContent =~/(\d*) \d* /
					def genomeCksum
					(1..genomeCksum_array.groupCount()).each{genomeCksum = "${genomeCksum_array[0][it]}"}
					predictionInstance.genome_cksum = "${genomeCksum}"
					genomeCksum_array = genomeCksumContent =~/\d* (\d*) /
					predictionInstance.genome_size
					(1..genomeCksum_array.groupCount()).each{predictionInstance.genome_size = "${genomeCksum_array[0][it]}"}
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - genome.fa is ${predictionInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
					cmdStr = "rm ${projectDir}/genome.cksum &> /dev/null"
					def delProcCksumGenome = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProcCksumGenome.waitFor()
					cmdStr = "rm ${projectDir}/genome_cksum.sh &> /dev/null"
					def delProcCkShGenome = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProcCkShGenome.waitFor()
				} // end of if(!(predictionInstance.genome_ftp_link == null))				
				

				// retrieve EST file
				if(!(predictionInstance.est_ftp_link == null)){
					// check whether the EST file is small enough for upload
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Checking cDNA file size with curl prior upload\n"
					def fileSizeScript = new File("${projectDir}/filzeSize.sh")
					cmd2Script = "curl -sI ${predictionInstance.est_ftp_link} | grep Content-Length | cut -d ' ' -f 2 > ${projectDir}/estFileSize 2> /dev/null"
					fileSizeScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - fileSizeScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${fileSizeScript}"
					def retrieveFileSize = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					retrieveFileSize.waitFor()	
					cmdStr = "rm ${fileSizeScript} &> /dev/null"		
					def delSzCrProc = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delSzCrProc.waitFor()
					content = new File("${projectDir}/estFileSize").text
					st = new Scanner(content)//works for exactly one number in a file
					def long est_size;
					est_size = st.nextLong();
					if(est_size < maxFileSizeByWget){//1 GB
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Retrieving EST/cDNA file ${predictionInstance.est_ftp_link}\n"
						def getEstScript = new File("${projectDir}/getEst.sh")
						cmd2Script = "wget -O ${projectDir}/est.fa ${predictionInstance.est_ftp_link} &> /dev/null"
						getEstScript << "${cmd2Script}"
						if(verb > 2){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - getEstScript << \"${cmd2Script}\"\n"
						}
						cmdStr = "bash ${projectDir}/getEst.sh"
						def wgetEst = "${cmdStr}".execute()
						if(verb > 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
						}
						wgetEst.waitFor()
						cmdStr = "rm ${projectDir}/getEst.sh &> /dev/null"
						delProc = "${cmdStr}".execute()
						if(verb > 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
						}
						delProc.waitFor()
						if("${predictionInstance.est_ftp_link}" =~ /\.gz/){
							def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
							cmd2Script = "cd ${projectDir}; mv est.fa est.fa.gz &> /dev/null; gunzip est.fa.gz 2> /dev/null"
							gunzipEstScript << "${cmd2Script}"
							if(verb > 2){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - gunzipEstScript << \"${cmd2Script}\"\n"
							}
							cmdStr = "bash ${gunzipEstScript}"
							def gunzipEst = "${cmdStr}".execute()
							if(verb > 1){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
							}
							gunzipEst.waitFor()	
							cmdStr = "rm ${gunzipEstScript} &> /dev/null"		
							delProc = "${cmdStr}".execute()
							if(verb > 1){
								logDate = new Date()
								logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
							}
							delProc.waitFor()
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Unpacked EST file.\n"
						}
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - EST/cDNA file upload finished, file stored as est.fa at ${projectDir}\n"
						// check for fasta format:
						new File("${projectDir}/est.fa").eachLine{line -> 
							if(line =~ /\*/ || line =~ /\?/){
								metacharacterFlag = 1
							}else{
								if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ 
									estFastaFlag = 1 
								}
							}
						}
						if(metacharacterFlag == 1){
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The cDNA file contains metacharacters (e.g. * or ?).\n";
							deleteDir()
							logAbort()
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided cDNA file\n${predictionInstance.est_ftp_link}\ncontains metacharacters (e.g. * or ?). This is not allowed.\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
							// delete database entry
							//predictionInstance.delete()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
						if(estFastaFlag == 1) {
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The EST/cDNA file was not fasta.\n"
							deleteDir()
							logAbort()
							mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided cDNA file\n${predictionInstance.est_ftp_link}\nwas not in DNA fasta format.\n\n"
							logDate = new Date()
							predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							if(predictionInstance.email_adress != null){
								msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
								sendMail {
									to "${predictionInstance.email_adress}"
									subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${msgStr}${footer}"""
								}
							}
							// delete database entry
							//predictionInstance.delete()
							predictionInstance.results_urls = null
							predictionInstance.job_status = 5
							predictionInstance = predictionInstance.merge()
							predictionInstance.save()
							return
						}
					}else{// actions if remote file was bigger than allowed
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - EST file size exceeds permitted ${maxFileSizeByWget} bytes. Abort job.\n"
						mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the cDNA file size was\nwith ${est_size} bigger than 1 GB. Please submitt a smaller cDNA size!\n\n"
						def errorStrMsg = "Hello!\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
						logDate = new Date()
						predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						if(predictionInstance.email_adress != null){
							sendMail {
									to "${predictionInstance.email_adress}"
									subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${errorStrMsg}${footer}"""
							}
						}
						deleteDir()
						predictionInstance.results_urls = null
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						return

					}
					def estCksumScript = new File("${projectDir}/est_cksum.sh")
					def estCksumFile = "${projectDir}/est.cksum"
					cmd2Script = "cksum ${projectDir}/est.fa > ${estCksumFile} 2> /dev/null"
					estCksumScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - estCksumScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${projectDir}/est_cksum.sh"
					def estCksumProcess = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					estCksumProcess.waitFor()
					def estCksumContent = new File("${estCksumFile}").text
					def estCksum_array = estCksumContent =~/(\d*) \d* /
					def estCksum
					(1..estCksum_array.groupCount()).each{estCksum = "${estCksum_array[0][it]}"}
					predictionInstance.est_cksum = "${estCksum}"
					estCksum_array = estCksumContent =~/\d* (\d*) /
					predictionInstance.est_size
					(1..estCksum_array.groupCount()).each{predictionInstance.est_size = "${estCksum_array[0][it]}"}
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - est.fa is ${predictionInstance.est_size} big and has a cksum of ${estCksum}.\n"
					cmdStr = "rm ${projectDir}/est.cksum &> /dev/null"
					def delProcCksumEst = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProcCksumEst.waitFor()
					cmdStr = "rm ${projectDir}/est_cksum.sh &> /dev/null"
					def delProcCkShEst = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					delProcCkShEst.waitFor()
				} // end of if(!(predictionInstance.est_ftp_link == null))

				// check whether EST file is NOT RNAseq, i.e. does not contain on average very short entries
				def int nEntries = 0
				def int totalLen = 0
				if(estExistsFlag == 1){
					new File("${projectDir}/est.fa").eachLine{line -> 
						if(line =~ /^>/){
							nEntries = nEntries + 1
						}else{
							totalLen = totalLen + line.size()
						}
					}
					def avEstLen = totalLen/nEntries
					if(avEstLen < estMinLen){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - EST sequences are on average shorter than ${estMinLen}, suspect RNAseq raw data.\n"
						logAbort()
						mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the sequences in your\ncDNA file have an average length of ${avEstLen}. We suspect that sequences files\nwith an average sequence length shorter than ${estMinLen} might contain RNAseq\nraw sequences. Currently, our web server application does not support the integration\nof RNAseq raw sequences. Please either assemble your sequences into longer contigs,\nor remove short sequences from your current file, or submitt a new job without\nspecifying a cDNA file.\n\n"
						def errorStrMsg = "Hello!\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
						logDate = new Date()
						predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						if(predictionInstance.email_adress != null){

							sendMail {
									to "${predictionInstance.email_adress}"
									subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
									body """${errorStrMsg}${footer}"""
							}
						}
						deleteDir()
						predictionInstance.results_urls = null
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						return
					}else if(avEstLen > estMaxLen){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - EST sequences are on average longer than ${estMaxLen}, suspect non EST/cDNA data.\n"
						logAbort()
						mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the sequences in your\ncDNA file have an average length of ${avEstLen}. We suspect that sequence\nfiles with an average sequence length longer than ${estMaxLen} might not contain\nESTs or cDNAs. Please either remove long sequences from your current file, or\nsubmitt a new job without specifying a cDNA file.\n\n"
						def errorStrMsg = "Hello!\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
						logDate = new Date()
						predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						if(predictionInstance.email_adress != null){
							sendMail {
								to "${predictionInstance.email_adress}"
								subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
								body """${errorStrMsg}${footer}"""
							}
						}
						deleteDir()
						// delete database entry
						//predictionInstance.delete()
						predictionInstance.results_urls = null
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						return
					}
				}

				// confirm file upload via e-mail
				if((!(predictionInstance.genome_ftp_link == null)) || (!(predictionInstance.est_ftp_link == null))){
					mailStr = "We have retrieved all files that you specified, successfully. You may delete them\nfrom the public server, now, without affecting the AUGUSTUS prediction job.\n\n"
					logDate = new Date()
					predictionInstance.message = "${predictionInstance.message}----------------------------------------\n${logDate} - Message:\n----------------------------------------\n\n${mailStr}"
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
					if(predictionInstance.email_adress != null){
						msgStr = "Hello!\n\n${mailStr}Best regards,\n\nthe AUGUSTUS web server team"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "File upload has been completed for AUGUSTUS prediction job ${predictionInstance.accession_id}"
							body """${msgStr}${footer}"""
						}
					}
				}

				// File formats appear to be ok. 
				// check whether this job was submitted before:
				def grepScript = new File("${projectDir}/grepScript.sh")
				def grepResult = "${projectDir}/grep.result"
				cmd2Script = "grep \"\\(Genome-Cksum: \\[${predictionInstance.genome_cksum}\\] Genome-Filesize: \\[${predictionInstance.genome_size}\\]\\)\" ${dbFile} | grep \"\\(EST-Cksum: \\[${predictionInstance.est_cksum}\\] EST-Filesize: \\[${predictionInstance.est_size}\\]\\)\" | grep \"\\(Hint-Cksum: \\[${predictionInstance.hint_cksum}\\] Hint-Filesize: \\[${predictionInstance.hint_size}\\] Parameter-String: \\[${predictionInstance.project_id}\\]\\)\" | grep \"\\(Parameter-Cksum: \\[${predictionInstance.archive_cksum}\\] Parameter-Size: \\[${predictionInstance.archive_size}\\] Server-Set-UTR-Flag: \\[${overRideUtrFlag}\\]\\)\" | grep \"\\(Report-Genes: \\[${predictionInstance.pred_strand}\\] Alternative-Transcripts: \\[${predictionInstance.alt_transcripts}\\] Gene-Structures: \\[${predictionInstance.allowed_structures}\\] Ignore-Conflicts: \\[${predictionInstance.ignore_conflicts}\\]\\)\"  > ${grepResult} 2> /dev/null"
				cmdStr = "bash ${projectDir}/grepScript.sh"
				grepScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - grepScript << \"${cmd2Script}\"\n"
				}
				def grepJob = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				grepJob.waitFor()
				def grepContent = new File("${grepResult}").text
				if(grepContent =~ /Genome-Cksum/){
					//job was submitted before. Send E-Mail to user with a link to the results.
					def id_array = grepContent =~ /Grails-ID: \[(\w*)\] /
					oldID
					(0..id_array.groupCount()).each{oldID = "${id_array[0][it]}"}
					def oldAccScript = new File("${projectDir}/oldAcc.sh")
					def oldAccResult = "${projectDir}/oldAcc.result"
					cmd2Script = "grep \"Grails-ID: \\[${oldID}\\]\" ${dbFile} | perl -ne \"@t = split(/\\[/); @t2 = split(/\\]/, \\\$t[4]); print \\\$t2[0];\" > ${oldAccResult} 2> /dev/null"
					oldAccScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - oldAccScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${projectDir}/oldAcc.sh"
					def oldAccScriptProc = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					oldAccScriptProc.waitFor()
					def oldAccContent = new File("${oldAccResult}").text
					mailStr = "You submitted job ${predictionInstance.accession_id}.\nThe job was aborted because the files that you submitted were submitted, before.\n\n"
					predictionInstance.old_url = "${war_url}prediction/show/${oldID}"
					logDate = new Date()
					predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}"
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
					if(predictionInstance.email_adress != null){
						msgStr = "Hello!\n\n${mailStr}The old job with identical input files and identical parameters"
						msgStr = "${msgStr} is available at\n${war_url}prediction/show/${oldID}.\n\nBest regards,\n\n"
						msgStr = "${msgStr}the AUGUSTUS web server team"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "AUGUSTUS prediction job ${predictionInstance.accession_id} was submitted before as job ${oldAccContent}"
							body """${msgStr}${footer}"""
						}
					}
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Data are identical to old job ${oldAccContent} with Accession-ID ${oldAccContent}.\n"
					deleteDir()
					logAbort()
					predictionInstance.results_urls = null
					predictionInstance.job_status = 5
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
					return
				} // end of job was submitted before check

				//Write DB file: 
				dbFile << "Date: [${today}] User-IP: [${userIP}] Grails-ID: [${predictionInstance.id}] Accession-ID: [${predictionInstance.accession_id}] Genome-File: [${predictionInstance.genome_file}] Genome-FTP-Link: [${predictionInstance.genome_ftp_link}] Genome-Cksum: [${predictionInstance.genome_cksum}] Genome-Filesize: [${predictionInstance.genome_size}] EST-File: [${predictionInstance.est_file}] EST-FTP-Link: [${predictionInstance.est_ftp_link}] EST-Cksum: [${predictionInstance.est_cksum}] EST-Filesize: [${predictionInstance.est_size}] Hint-File: [${predictionInstance.hint_file}] Hint-Cksum: [${predictionInstance.hint_cksum}] Hint-Filesize: [${predictionInstance.hint_size}] Parameter-String: [${predictionInstance.project_id}] Parameter-File: [${predictionInstance.archive_file}] Parameter-Cksum: [${predictionInstance.archive_cksum}] Parameter-Size: [${predictionInstance.archive_size}] Server-Set-UTR-Flag: [${overRideUtrFlag}] User-Set-UTR-Flag: [${predictionInstance.utr}] Report-Genes: [${predictionInstance.pred_strand}] Alternative-Transcripts: [${predictionInstance.alt_transcripts}] Gene-Structures: [${predictionInstance.allowed_structures}] Ignore-Conflicts: [${predictionInstance.ignore_conflicts}]\n"

				//rename and move parameters
				if(!uploadedParamArch.empty){
					def mvParamsScript = new File("${projectDir}/mvParams.sh")
					cmd2Script = "${AUGUSTUS_SCRIPTS_PATH}/moveParameters.pl ${projectDir}/params ${predictionInstance.accession_id} ${AUGUSTUS_CONFIG_PATH}/species 2> /dev/null\n"
					mvParamsScript << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - mvParamsScript << \"${cmd2Script}\"\n"
					}
					cmdStr = "bash ${mvParamsScript}"
					def mvParamsRunning = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					mvParamsRunning.waitFor()
					species = "${predictionInstance.accession_id}"
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Moved uploaded parameters and renamed species to ${predictionInstance.accession_id}\n"
				}
				//Create sge script:
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Writing SGE submission script.\n"
				def sgeFile = new File("${projectDir}/aug-pred.sh")
				// write command in script (according to uploaded files)
				sgeFile << "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n\n"
				def cmdStr = "mkdir ${projectDir}/augustus\n"
				if(estExistsFlag == 1){
					cmdStr = "${cmdStr}${BLAT_PATH} -noHead ${projectDir}/genome.fa ${projectDir}/est.fa ${projectDir}/est.psl\n"
					cmdStr = "${cmdStr}cat ${projectDir}/est.psl | sort -n -k 16,16 | sort -s -k 14,14 > ${projectDir}/est.s.psl\n"
					cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/blat2hints.pl --in=${projectDir}/est.s.psl --out=${projectDir}/est.hints --source=E\n"
					cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/blat2gbrowse.pl ${projectDir}/est.s.psl ${projectDir}/est.gbrowse\n"
				}
				if(hintExistsFlag == 1){
					cmdStr = "${cmdStr}cat ${projectDir}/hints.gff >> ${projectDir}/est.hints\n"
				}
				if((hintExistsFlag == 1) || (estExistsFlag == 1)){
					radioParameterString = "${radioParameterString} --hintsfile=${projectDir}/est.hints --extrinsicCfgFile=${AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.ME.cfg"
				}
				cmdStr = "${cmdStr}cd ${projectDir}/augustus\naugustus --species=${species} ${radioParameterString} ${projectDir}/genome.fa --codingseq=on --exonnames=on > ${projectDir}/augustus/augustus.gff\n"
				cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/getAnnoFasta.pl --seqfile=${projectDir}/genome.fa ${projectDir}/augustus/augustus.gff\n"
				cmdStr = "${cmdStr}cat ${projectDir}/augustus/augustus.gff | perl -ne 'if(m/\\tAUGUSTUS\\t/){print;}' > ${projectDir}/augustus/augustus.gtf\n"
				cmdStr = "${cmdStr}cat ${projectDir}/augustus/augustus.gff | ${AUGUSTUS_SCRIPTS_PATH}/augustus2gbrowse.pl > ${projectDir}/augustus/augustus.gbrowse\n"
				cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/writeResultsPage.pl ${predictionInstance.accession_id} null ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} 1 2> ${projectDir}/writeResults.err"
				sgeFile << "${cmdStr}"
				
				// write submission script
				def submissionScript = new File("${projectDir}/submitt.sh")
				def fileID = "${projectDir}/jobID"
				cmd2Script = "cd ${projectDir}; qsub aug-pred.sh > ${fileID} 2> /dev/null"
				submissionScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - submissionScript << \"${cmd2Script}\"\n"
				}
				// submitt job
				cmdStr = "bash ${projectDir}/submitt.sh"
				def jobSubmission = "${cmdStr}".execute()
				if(verb > 1){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
				}
				jobSubmission.waitFor()
				// get job ID
				content = new File("${fileID}").text
				def jobID_array = content =~/Your job (\d*)/
				def jobID
				(1..jobID_array.groupCount()).each{jobID = "${jobID_array[0][it]}"}
				predictionInstance.job_id = jobID
				logDate = new Date()
				logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${jobID} submitted.\n"
				// check for job status
				predictionInstance.job_status = 1 // submitted
				predictionInstance = predictionInstance.merge()
				predictionInstance.save()
				def statusScript = new File("${projectDir}/status.sh")
				def statusFile = "${projectDir}/job.status"
				cmd2Script = "cd ${projectDir}; qstat -u \"*\" | grep aug-pred | grep ${jobID} > ${statusFile} 2> /dev/null"
				statusScript << "${cmd2Script}"
				if(verb > 2){
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - statusScript <<\"${cmd2Script}\"\n"
				}
				def statusContent
				def statusCheck 
				def qstat = 1
				def runFlag = 0;

				while(qstat == 1){
					sleep(300000) // 300000 = 5 minutes
					cmdStr = "bash ${projectDir}/status.sh"
					statusCheck = "${cmdStr}".execute()
					statusCheck.waitFor()
					sleep(100)
					statusContent = new File("${statusFile}").text
					if(statusContent =~ /qw/){ 
						predictionInstance.job_status = 2 
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
					}else if( statusContent =~ /  r  / ){
						predictionInstance.job_status = 3
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						if(runFlag == 0){
							today = new Date() 
							logDate = new Date()
							logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${jobID} begins running at ${today}.\n"
						}
						runFlag = 1
					}else if(!statusContent.empty){
						predictionInstance.job_status = 3
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						logFile << "${logDate} ${predictionInstance.accession_id} v1 - Job ${jobID} is neither in qw nor in r status but is still on the grid!\n"
					}else{
						predictionInstance.job_status = 4
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
						qstat = 0
						today = new Date()
 						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job ${jobID} left SGE at ${today}.\n"
					}
			   	}
				// collect results link information
				if(new File("${web_output_dir}/${predictionInstance.accession_id}/predictions.tar.gz").exists()){
					predictionInstance.results_urls = "<p><b>Prediction archive</b>&nbsp;&nbsp;<a href=\"${web_output_url}${predictionInstance.accession_id}/predictions.tar.gz\">predictions.tar.gz</a><br></p>"
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
				}
			   	// check whether errors occured by log-file-sizes
				if(new File("${projectDir}/aug-pred.sh.e${jobID}").exists()){
					sgeErrSize = new File("${projectDir}/aug-pred.sh.e${jobID}").size()
				}else{
					sgeErrSize = 10
 					logDate = new Date()
					logFile <<  "SEVERE ${logDate} ${predictionInstance.accession_id} v1 - segErrFile was not created. Setting size to default value 10.\n"
				}
				if(new File("${projectDir}/writeResults.err").exists()){
					writeResultsErrSize = new File("${projectDir}/writeResults.err").size()
				}else{
					writeResultsErrSize = 10
 					logDate = new Date()
					logFile <<  "SEVERE ${logDate} ${predictionInstance.accession_id} v1 - writeResultsErr was not created. Setting size to default value 10.\n"
				}
				if(sgeErrSize==0 && writeResultsErrSize==0){
					mailStr = "Your AUGUSTUS prediction job ${predictionInstance.accession_id} finished.\n\n"
					logDate = new Date()
					predictionInstance.message = "${predictionInstance.message}----------------------------------------\n${logDate} - Message:\n----------------------------------------\n\n${mailStr}"
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
					if(predictionInstance.email_adress == null){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Computation was successful. Did not send e-mail to user because no e-mail adress was supplied.\n"
					}
					if(predictionInstance.email_adress != null){
						msgStr = "Hello!\n\n${mailStr}You find the results at "
						msgStr = "${msgStr}${war_url}prediction/show/${predictionInstance.id}.\n\nBest regards,\n\n"
						msgStr = "${msgStr}the AUGUSTUS web server team"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "AUGUSTUS prediction job ${predictionInstance.accession_id} is complete"
							body """${msgStr}${footer}"""
						}
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Sent confirmation Mail that job computation was successful.\n"
					}
					// unpack with 7z x XA2Y5VMJ.tar.7z
					// tar xvf XA2Y5VMJ.tar
					def packResults = new File("${output_dir}/pack${predictionInstance.accession_id}.sh")
					cmd2Script = "cd ${output_dir}; tar -czvf ${predictionInstance.accession_id}.tar.gz ${predictionInstance.accession_id} &> /dev/null"
					packResults << "${cmd2Script}"
					if(verb > 2){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v3 - packResults << \"${cmd2Script}\"\n"
					}
					//packResults << "cd ${output_dir}; tar cf - ${predictionInstance.accession_id} | 7z a -si ${predictionInstance.accession_id}.tar.7z; rm -r ${predictionInstance.accession_id};"
					cmdStr = "bash ${output_dir}/pack${predictionInstance.accession_id}.sh"
					def cleanUp = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					cleanUp.waitFor()
					cmdStr = "rm ${output_dir}/pack${predictionInstance.accession_id}.sh &> /dev/null"
					cleanUp = "${cmdStr}".execute()
					if(verb > 1){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v2 - \"${cmdStr}\"\n"
					}
					deleteDir()
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - job directory was packed with tar/gz.\n"
					logDate = new Date()
					logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Job completed. Result: ok.\n"
				}else{ 
					if(sgeErrSize > 0){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - a SGE error occured!\n";
						msgStr = "Hi ${admin_email}!\n\nJob: ${predictionInstance.accession_id}\n"
						msgStr = "${msgStr}IP: ${userIP}\n"
						msgStr = "${msgStr}E-Mail: ${predictionInstance.email_adress}\n"
						msgStr = "${msgStr}Link: ${war_url}prediction/show/${predictionInstance.id}\n\n"
						msgStr = "${msgStr}An SGE error occured. Please check manually what's wrong. "
						if(predictionInstance.email_adress == null){
							msgStr = "${msgStr}The user has not been informed."
							sendMail {
							to "${admin_email}"
							subject "Error in AUGUSTUS prediction job ${predictionInstance.accession_id}"
							body """${msgStr}${footer}"""	
							}
						}else{
							msgStr = "${msgStr}The user has been informed."
							sendMail {
							to "${admin_email}"
							subject "Error in AUGUSTUS prediction job ${predictionInstance.accession_id}"
							body """${msgStr}${footer}"""	
							}
						}
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
					}else{
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - an error occured during writing results!\n";
						msgStr = "Hi ${admin_email}!\n\nJob: ${predictionInstance.accession_id}\n"
						msgStr = "${msgStr}IP: ${userIP}\n"
						msgStr = "${msgStr}E-Mail: ${predictionInstance.email_adress}\n"
						msgStr = "${msgStr}Link: ${war_url}prediction/show/${predictionInstance.id}\n\n"
						msgStr = "${msgStr}An error occured during writing results.. Please check manually what's wrong. "
						if(predictionInstance.email_adress == null){
							msgStr = "${msgStr} The user has not been informed."
							sendMail {
								to "${admin_email}"
								subject "Error in AUGUSTUS prediction job ${predictionInstance.accession_id}"
								body """${msgStr}${footer}"""
							}
						}else{
							msgStr = "${msgStr} The user has been informed."
							sendMail {
								to "${admin_email}"
								subject "Error in AUGUSTUS prediction job ${predictionInstance.accession_id}"
								body """${msgStr}${footer}"""
							}
						}
						predictionInstance.job_status = 5
						predictionInstance = predictionInstance.merge()
						predictionInstance.save()
					}
					mailStr = "An error occured while running the AUGUSTUS prediction job ${predictionInstance.accession_id}.\n\n"
					logDate = new Date()
					predictionInstance.message = "${predictionInstance.message}----------------------------------------------\n${logDate} - Error Message:\n----------------------------------------------\n\n${mailStr}Please contact augustus-web@uni-greifswald.de if you want to find out what went wrong.\n\n"
					predictionInstance = predictionInstance.merge()
					predictionInstance.save()
					if(predictionInstance.email_adress == null){
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - The job is in an error state. Cound not send e-mail to anonymous user because no email adress was supplied.\n"
					}else{
						msgStr = "Hello!\n\n${mailStr}The administrator of the AUGUSTUS web server has been informed and"
						msgStr = "${msgStr} will get back to you as soon as the problem is solved.\n\nBest regards,\n\n"
						msgStr = "${msgStr}the AUGUSTUS web server team"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "An error occured while executing AUGUSTUS prediction job ${predictionInstance.accession_id}"
							body """${msgStr}${footer}"""
						}
						logDate = new Date()
						logFile <<  "${logDate} ${predictionInstance.accession_id} v1 - Sent confirmation Mail, the job is in an error state.\n"
					}
				}
			}
			//------------ END BACKGROUND PROCESS ----------------------------------
		} // end of (!(predictionInstance.id == null))
	}// end of commit
} // end of Controller
