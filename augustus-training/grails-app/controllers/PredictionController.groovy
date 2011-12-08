// The class PredictionController controls everything that is related to submitting a job for predicting genes with pre-trained parameters on a novel genome
//    - it handles the file upload (or wget)
//    - format check
//    - SGE job submission and status checks
//    - rendering of results/job status page
//    - sending E-Mails concerning the job status (submission, errors, finished)

class PredictionController {
	// need to adjust the output dir to whatever working dir! This is where uploaded files and results will be saved.
	def output_dir = "/data/www/augpred/webdata" // should be something in home of webserver user and augustus frontend user.
	// this log File contains the "process log", what was happening with which job when.
	def logFile = new File("${output_dir}/augustus-prediction.log")
	// this log File contains the "database" (not identical with the grails database and simply for logging purpose)
	def dbFile = new File("${output_dir}/augustus-pred-database.log")
	// oldID is a parameter that is used for show redirects (see bottom)
	def oldID
	def oldAccID
	// web-output, root directory to the results that are shown to end users
	def web_output_dir = "/var/www/trainaugustus/prediction-results" // must be writable to webserver application
	// AUGUSTUS_CONFIG_PATH
	def AUGUSTUS_CONFIG_PATH = "/usr/local/augustus/trunks/config"
	def AUGUSTUS_SCRIPTS_PATH = "/usr/local/augustus/trunks/scripts"
	def BLAT_PATH = "/usr/local/blat/blat"
	def scaffold = Prediction
	// Admin mail for errors
	def admin_email = "katharina.hoff@gmail.com"
	// sgeLen length of SGE queue, when is reached "the server is buisy" will be displayed
	def sgeLen = 8;

	// human verification:
	def simpleCaptchaService

	// check whether the server is buisy
	def beforeInterceptor = {
		def String prefixChars ="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_"
		def rnd = new Random()
		def qstatFilePrefix = (1..10).sum{prefixChars[ rnd.nextInt(prefixChars.length()) ]} 
		def qstatFile = new File("${output_dir}/${qstatFilePrefix}.qstatScript")
		qstatFile << "qstat | grep qw | wc -l > ${output_dir}/${qstatFilePrefix}.qstatResult"
		def qstatStatus = "bash ${output_dir}/${qstatFilePrefix}.qstatScript".execute()
		qstatStatus.waitFor()
		def qstatStatusResult = new File("${output_dir}/${qstatFilePrefix}.qstatResult").text
		def qstatStatus_array = qstatStatusResult =~ /(\d*)/
		def qstatStatusNumber 
		(1..qstatStatus_array.groupCount()).each{qstatStatusNumber = "${qstatStatus_array[0][it]}"}
		def delProc = "rm -r ${output_dir}/${qstatFilePrefix}.qstatScript".execute()
		delProc.waitFor()
		delProc = "rm -r ${output_dir}/${qstatFilePrefix}.qstatResult".execute()
		delProc.waitFor()

		if(qstatStatusNumber > sgeLen){
			// get date
			def todayTried = new Date()
			// get IP-address
			String userIPTried = request.remoteAddr
			logFile <<  "SGE      On ${todayTried} somebody with IP ${userIPTried} tried to invoke the Training webserver but the SGE queue was longer than ${sgeLen} and the user was informed that submission is currently not possible\n"
			render "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\"/><meta name=\"layout\" content=\"main\" /><title>Submitt Prediction</title><script type=\"text/javascript\" src=\"js/md_stylechanger.js\"></script></head><body><!-- ***** Start: Kopfbereich ********************************************// --><p class=\"unsichtbar\"><a href=\"#inhalt\" title=\"Directly to Contents\">Directly to Contents</a></p><div id=\"navigation_oben\"><a name=\"seitenanfang\"></a><table width=\"100%\" border=\"0\" cellpadding=\"0\" cellspacing=\"1\"><tr><td nowrap=\"nowrap\"><a href=\"http://www.uni-greifswald.de\" target=\"_blank\" class=\"mainleveltop_\" >University of Greifswald</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.mnf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Faculty</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.math-inf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Institute</a><span class=\"mainleveltop_\">&nbsp;|&nbsp;</span><a href=\"http://bioinf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\">Bioinformatics Group</a></td></tr></table></div><div id=\"banner\"><div id=\"banner_links\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php\" title=\"Institut f&uuml;r Mathematik und Informatik\"><img src=\"../images/header.gif\" alt=\"Directly to home\" /> </a></div><div id=\"banner_mitte\"><div id=\"bannertitel1\">Bioinformatics Web Server at University of Greifswald</div><div id=\"bannertitel2\">Gene Prediction with AUGUSTUS</div></div><div id=\"banner_rechts\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php/geschichte-und-kultur/167\" title=\"Voderberg-Doppelspirale\"><img src=\"../images/spirale.gif\" align=\"left\" /></a></div></div><div id=\"wegweiser\">Navigation for: &nbsp; &nbsp;<span class=\"breadcrumbs pathway\">Submitt Prediction</span><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Kopfbereich *********************************************// --><!-- ***** Start: Koerper ************************************************// --><div id=\"koerper\"><div id=\"linke_spalte\"><ul class=\"menu\"><li><a href=\"../index.gsp\"><span>Introduction</span></a></li><li><a href=\"/augustus-training/training/create\"><span>Submitt Training</span></a></li><li><a href=\"/augustus-training/prediction/create\"><span>Submitt Prediction</span></a></li><li><a href=\"../help.gsp\"><span>Help</span></a></li><li><a href=\"../references.gsp\"><span>Links & References</span></a></li><li><a href=\"http://bioinf.uni-greifswald.de\"><span>Bioinformatics Group</span></a></li><li><a href=\"http://bioinf.uni-greifswald.de/bioinf/impressum.html\"><span>Impressum</span></a></li></ul></div><div id=\"mittel_spalte\"><div class=\"main\" id=\"main\"><h1><font color=\"#006699\">The Server is Busy</font></h1><p>You tried to access the AUGUSTUS prediction job submission page.</p><p>Predicting genes with AUGUSTUS is a process that takes a lot of computation time. We estimate that one prediction process requires approximately 7 days. Our web server is able to process a certain number of jobs in parallel, and we established a waiting queue. The waiting queue has a limited length, though. Currently, all slots for computation and for waiting are occupied.</p><p>We apologize for the inconvenience! Please try to submitt your job in a couple of weeks, again.</p><p>Feel free to contact us in case your job is particularly urgent.</p></div><p>&nbsp;</p>           </div><div id=\"rechte_spalte\"><div class=\"linien_div\"><h5 class=\"ueberschrift_spezial\">CONTACT</h5><strong>Institute for Mathematics und Computer Sciences</strong><br/><strong>Bioinformatics Group</strong><br />Walther-Rathenau-Stra&szlig;e 47<br />17487 Greifswald<br />Germany<br />Tel.: +49 (0)3834 86 - 46 24<br/>Fax:  +49 (0)3834 86 - 46 40<br /><br /><a href=\"mailto:augustus-web@uni-greifswald.de\" title=\"E-Mail augustus-web@uni-greifswald.de, opens the standard mail program\">augustus-web@uni-greifswald.de</a></div></div><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Koerper *************************************************// --><!-- ***** Start: Fuss ***************************************************// --><div id=\"fuss\"><div id=\"fuss_links\"><p class=\"copyright\">&copy; 2011 University of Greifswald</p></div><div id=\"fuss_mitte\"><div class=\"bannergroup\"></div></div><div id=\"fuss_rechts\" ><ul><li><a href=\"#seitenanfang\"><img hspace=\"5\" height=\"4\" border=\"0\" width=\"7\" alt=\"Seitenanfang\" src=\"../images/top.gif\" />Top of page</a></li></ul></div><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Fuss ***************************************************// --></body></html>"
			return
		}		
	} 

	// the method commit is started if the "Submit Job" button on the website is hit. It is the main method of Prediction Controller and contains a Thread method that will continue running as a background process after the user is redirected to the job status page.

	def commit = {
		def predictionInstance = new Prediction(params)
		if(!(predictionInstance.id == null)){
			redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
			return
		}else{
			// info string for confirmation E-Mail
			def confirmationString
			// this_project_id that is used internally by pipeline as species name.
			def this_project_id = "web" + predictionInstance.accession_id
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
			// species name for AUGUSTUS
			def species
			// get date
			def today = new Date()
			logFile << "${predictionInstance.accession_id} AUGUSTUS prediction webserver starting on ${today}\n"
      			// get IP-address
      			String userIP = request.remoteAddr
      			logFile <<  "${predictionInstance.accession_id} user IP: ${userIP}\n"
			//verify that the submitter is a person
			boolean captchaValid = simpleCaptchaService.validateCaptcha(params.captcha)
			if(captchaValid == false){
				logFile << "${predictionInstance.accession_id} The user is probably not a human person. Job aborted.\n"
				flash.error = "The verification string at the bottom of the page was not entered correctly!"
            			redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
           			return
			}
			// utr checkbox
			if(predictionInstance.utr == true){
				overRideUtrFlag = 1
				logFile << "${predictionInstance.accession_id} User enabled UTR prediction.\n"
			}else{
				overRideUtrFlag = 0
				logFile << "${predictionInstance.accession_id} User did not enable UTR prediction.\n"
			}
			// get parameter archive file (if available)
			def uploadedParamArch = request.getFile('ArchiveFile')
			def String dirName = "${output_dir}/${predictionInstance.accession_id}"
			def projectDir = new File(dirName)
			if(!uploadedParamArch.empty){
				// actually upload the file
				projectDir.mkdirs()
         			uploadedParamArch.transferTo( new File (projectDir, "parameters.tar.gz"))
				predictionInstance.archive_file = uploadedParamArch.originalFilename
				confirmationString = "${confirmationString}Parameter archive: ${predictionInstance.archive_file}\n"
				logFile <<  "${predictionInstance.accession_id} uploaded parameter archive ${predictionInstance.archive_file} was renamed to parameters.tar.gz and moved to ${projectDir}\n"
				// get cksum and file size for database
				def archCksumScript = new File("${projectDir}/archive_cksum.sh")
         			def archCksumFile = "${projectDir}/arch.cksum"
         			archCksumScript << "cksum ${projectDir}/parameters.tar.gz > ${archCksumFile}"
         			def archCksumProcess = "bash ${projectDir}/archive_cksum.sh".execute()
         			archCksumProcess.waitFor()
         			def archCksumContent = new File("${archCksumFile}").text
         			def archCksum_array = archCksumContent =~/(\d*) \d* /
         			def archCksum
         			(1..archCksum_array.groupCount()).each{archCksum = "${archCksum_array[0][it]}"}
         			predictionInstance.archive_cksum = "${archCksum}"
         			predictionInstance.archive_size = uploadedParamArch.size
         			logFile <<  "${predictionInstance.accession_id} parameters.tar.gz is ${predictionInstance.archive_size} big and has a cksum of ${archCksum}.\n"
         			def delProcCksumarch = "rm ${projectDir}/arch.cksum".execute()
         			delProcCksumarch.waitFor()
         			def delProcCkSharch = "rm ${projectDir}/archive_cksum.sh".execute()
         			delProcCkSharch.waitFor()
				// check whether the archive contains all relevant files
				def String paramDirName = "${projectDir}/params"
				def MakeParamDir = "mkdir ${paramDirName}".execute()
				MakeParamDir.waitFor()
				def checkParamArch = new File("${projectDir}/ckArch.sh")
				checkParamArch << "checkParamArchive.pl ${projectDir}/parameters.tar.gz ${paramDirName} > ${projectDir}/archCheck.log 2> ${projectDir}/archCheck.err"
				def checkParamArchRunning = "bash ${checkParamArch}".execute()
				checkParamArchRunning.waitFor()
				def archCheckLog = new File("${projectDir}/archCheck.log")
				def archCheckErr = new File("${projectDir}/archCheck.err")
				def archCheckLogSize = archCheckLog.text.size()
				def archCheckErrSize = archCheckErr.text.size()
				// if essential file are missing, redirect to input interface and inform user that the archive was not compatible
				if(archCheckErrSize > 0){
					logFile <<  "${predictionInstance.accession_id} The parameter archive was not compatible. Project directory ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
            				delProc.waitFor()
           				logFile <<  "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
           				flash.error = "Parameter archive ${uploadedParamArch.originalFilename} is not compatible with the AUGUSTUS prediction web server application."
            				redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
           				return
				// if only UTR params are missing, set flag to override any user-defined UTR settings
				}else if(archCheckLogSize > 0){
					overRideUtrFlag = 0 // UTR predictions are now permanently disabled
					logFile <<  "${predictionInstance.accession_id} UTR predictions have been disabled because UTR parameters are missing!\n"
				}
				archiveExistsFlag = 1
			}else{predictionInstance.archive_file = "empty"}
	
			// check whether parameters are available for project_id (previous prediction run)
			logFile <<  "${predictionInstance.accession_id} The given parameter ID is ${predictionInstance.project_id}\n"
			if(!(predictionInstance.project_id == null)){
				def spec_conf_dir = new File("${AUGUSTUS_CONFIG_PATH}/species/${predictionInstance.project_id}")
				if(!spec_conf_dir.exists()){
					logFile <<  "${predictionInstance.accession_id} The given parameter-string does not exist on our system. Project directory ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
            				delProc.waitFor()
           				logFile <<  "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
           				flash.error = "The specified parameter ID ${predictionInstance.project_id} does not exist on our system."
            				redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
           				return
				}else{
					logFile <<  "${predictionInstance.accession_id} Requested ${spec_conf_dir} exists on our system.\n"
					speciesNameExistsFlag = 1
					species = predictionInstance.project_id
				}
				confirmationString = "${confirmationString}AUGUSTUS parameter project identifier: ${predictionInstance.project_id}\n"
			}
	
			// upload of genome file
			def uploadedGenomeFile = request.getFile('GenomeFile')
     			def seqNames = []
			if(!uploadedGenomeFile.empty){
         			projectDir.mkdirs()
         			uploadedGenomeFile.transferTo( new File (projectDir, "genome.fa"))
        			predictionInstance.genome_file = uploadedGenomeFile.originalFilename
				confirmationString = "${confirmationString}Genome file: ${predictionInstance.genome_file}\n"
				if("${uploadedGenomeFile.originalFilename}" =~ /\.gz/){
					logFile <<  "${predictionInstance.accession_id} Genome file is gzipped.\n"
					def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
					gunzipGenomeScript << "cd ${projectDir}; mv genome.fa genome.fa.gz; gunzip genome.fa.gz"
					def gunzipGenome = "bash ${gunzipGenomeScript}".execute()
					gunzipGenome.waitFor()
					logFile <<  "${predictionInstance.accession_id} Unpacked genome file.\n"
					def delProc = "rm ${gunzipGenomeScript}".execute()
					delProc.waitFor()
				}
         			logFile <<  "${predictionInstance.accession_id} uploaded genome file ${uploadedGenomeFile.originalFilename} was renamed to genome.fa and moved to ${projectDir}\n"
        			// check for fasta format & extract fasta headers for gff validation:
         			new File("${projectDir}/genome.fa").eachLine{line -> 
            			if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }
            				if(line =~ /^>/){
               					def len = line.length()
               					seqNames << line[1..(len-1)]
            				}
         			}
         			if(genomeFastaFlag == 1) {
            				logFile <<  "${predictionInstance.accession_id} The genome file was not fasta. Project directory ${projectDir} is deleted (rm -r).\n"
            				def delProc = "rm -r ${projectDir}".execute()
            				delProc.waitFor()
            				logFile <<  "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
            				flash.error = "Genome file ${uploadedGenomeFile.originalFilename} is not in DNA fasta format."
            				redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
					return
	         		} else {
	            			def genomeCksumScript = new File("${projectDir}/genome_cksum.sh")
	            			def genomeCksumFile = "${projectDir}/genome.cksum"
	            			genomeCksumScript << "cksum ${projectDir}/genome.fa > ${genomeCksumFile}"
	            			def genomeCksumProcess = "bash ${projectDir}/genome_cksum.sh".execute()
	            			genomeCksumProcess.waitFor()
	            			def genomeCksumContent = new File("${genomeCksumFile}").text
	            			def genomeCksum_array = genomeCksumContent =~/(\d*) \d* /
	            			def genomeCksum
	            			(1..genomeCksum_array.groupCount()).each{genomeCksum = "${genomeCksum_array[0][it]}"}
	            			predictionInstance.genome_cksum = "${genomeCksum}"
	            			predictionInstance.genome_size = uploadedGenomeFile.size
	            			logFile <<  "${predictionInstance.accession_id} genome.fa is ${predictionInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
	            			def delProcCksumGenome = "rm ${projectDir}/genome.cksum".execute()
	            			delProcCksumGenome.waitFor()
	            			def delProcCkShGenome = "rm ${projectDir}/genome_cksum.sh".execute()
	            			delProcCkShGenome.waitFor()
	         		}
			}
	
			// retrieve beginning of genome file for format check
	      		if(!(predictionInstance.genome_ftp_link == null)){
				confirmationString = "${confirmationString}Genome file: ${predictionInstance.genome_ftp_link}\n"
	         		logFile <<  "${predictionInstance.accession_id} genome web-link is ${predictionInstance.genome_ftp_link}\n"
	         		projectDir.mkdirs()
	         		// checking web file for DNA fasta format: 
	         		def URL url = new URL("${predictionInstance.genome_ftp_link}");
	         		def URLConnection uc = url .openConnection()
				if(!("${predictionInstance.genome_ftp_link}" =~ /\.gz/)){
	         			def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
	         			def String inputLine=null
	         			def lineCounter = 1;
	         			while ( ((inputLine = br.readLine()) != null) && (lineCounter <= 20)) {
	            				if(!(inputLine =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(inputLine =~ /^$/)){ genomeFastaFlag = 1 }
	            				lineCounter = lineCounter + 1
	         			}
	         			br.close()
	         			if(genomeFastaFlag == 1) {
	            				logFile <<  "${predictionInstance.accession_id} The first 20 lines in genome file are not fasta.\n"
	            				def delProc = "rm -r ${projectDir}".execute()
	            				delProc.waitFor()
	            				logFile << "${predictionInstance.accession_id} Project directory ${projectDir} is deleted.\n${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"	
	            				flash.error = "Genome file ${predictionInstance.genome_ftp_link} is not in DNA fasta format."
	            				redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
	            				return
	         			}
				}else{
					logFile <<  "${predictionInstance.accession_id} The linked genome file is gzipped. Format will be checked later after extraction.\n"
				}
	      		}
	
	      		// upload of est file
	      		def uploadedEstFile = request.getFile('EstFile')
	      		if(!uploadedEstFile.empty){
	         		projectDir.mkdirs()
	         		uploadedEstFile.transferTo( new File (projectDir, "est.fa"))
	         		predictionInstance.est_file = uploadedEstFile.originalFilename
				confirmationString = "${confirmationString}cDNA file: ${predictionInstance.est_file}\n"
				if("${uploadedEstFile.originalFilename}" =~ /\.gz/){
					logFile <<  "${predictionInstance.accession_id} EST file is gzipped.\n"
					def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
					gunzipEstScript << "cd ${projectDir}; mv est.fa est.fa.gz; gunzip est.fa.gz"
					def gunzipEst = "bash ${gunzipEstScript}".execute()
					gunzipEst.waitFor()
					logFile <<  "${predictionInstance.accession_id} Unpacked EST file.\n"
					def delProc = "rm ${gunzipEstScript}".execute()
					delProc.waitFor()
				}
	         		logFile << "${predictionInstance.accession_id} Uploaded EST file ${uploadedEstFile.originalFilename} was renamed to est.fa and moved to ${projectDir}\n"
	         		// check fasta format
	         		new File("${projectDir}/est.fa").eachLine{line -> 
	            			if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNnUu]/) && !(line =~ /^$/)){ estFastaFlag = 1 }
	         		}	
	         		if(estFastaFlag == 1) {
            				logFile << "${predictionInstance.accession_id} The cDNA file was not fasta. ${projectDir} (rm -r) is deleted.\n"
            				def delProc = "rm -r ${projectDir}".execute()
            				delProc.waitFor()
            				logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
            				flash.error = "cDNA file ${uploadedEstFile.originalFilename} is not in DNA fasta format."
            				redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
            				return
         			} else { estExistsFlag = 1 }
         				def estCksumScript = new File("${projectDir}/est_cksum.sh")
         				def estCksumFile = "${projectDir}/est.cksum"
         				estCksumScript << "cksum ${projectDir}/est.fa > ${estCksumFile}"
         				def estCksumProcess = "bash ${projectDir}/est_cksum.sh".execute()
         				estCksumProcess.waitFor()
         				def estCksumContent = new File("${estCksumFile}").text
         				def estCksum_array = estCksumContent =~/(\d*) \d* /
         				def estCksum
         				(1..estCksum_array.groupCount()).each{estCksum = "${estCksum_array[0][it]}"}
         				predictionInstance.est_cksum = "${estCksum}"
         				predictionInstance.est_size = uploadedEstFile.size
         				logFile <<  "${predictionInstance.accession_id} est.fa is ${predictionInstance.est_size} big and has a cksum of ${estCksum}.\n"
         				def delProcCksumEst = "rm ${projectDir}/est.cksum".execute()
         				delProcCksumEst.waitFor()
         				def delProcCkShEst = "rm ${projectDir}/est_cksum.sh".execute()
         				delProcCkShEst.waitFor()
      			}

      			// retrieve beginning of est file for format check
      			if(!(predictionInstance.est_ftp_link == null)){
				confirmationString = "${confirmationString}cDNA file: ${predictionInstance.est_ftp_link}\n"
         			logFile << "${predictionInstance.accession_id} est web-link is ${predictionInstance.est_ftp_link}\n"
         			projectDir.mkdirs()
				estExistsFlag = 1
				if(!("${predictionInstance.est_ftp_link}" =~ /\.gz/)){
         				// checking web file for DNA fasta format: 
         				def URL url = new URL("${predictionInstance.est_ftp_link}");
         				def URLConnection uc = url .openConnection()
         				def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
         				def String inputLine=null
         				def lineCounter = 1
         				while ( ((inputLine = br.readLine()) != null) && (lineCounter <= 20)) {
            					if(!(inputLine =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNnUu]/) && !(inputLine =~ /^$/)){ estFastaFlag = 1 }
            					lineCounter = lineCounter + 1
         				}
         				br.close()
         				if(estFastaFlag == 1) {
           					logFile << "${predictionInstance.accession_id} The cDNA file was not fasta. ${projectDir} is deleted (rm -r).\n"
            					def delProc = "rm -r ${projectDir}".execute()
         					delProc.waitFor()
            					logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
            					flash.error = "cDNA file ${predictionInstance.est_ftp_link} is not in DNA fasta format."
            					redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
            					return
         				}
				}else{
					logFile <<  "${predictionInstance.accession_id} The linked EST file is gzipped. Format will be checked later after extraction.\n"
				}
      			}

			// get hints file, format check
			def uploadedStructFile = request.getFile('HintFile')
			if(!uploadedStructFile.empty){
				projectDir.mkdirs()
				uploadedStructFile.transferTo( new File (projectDir, "hints.gff"))
				predictionInstance.hint_file = uploadedStructFile.originalFilename
				confirmationString = "${confirmationString}Hints file: ${predictionInstance.hint_file}\n"
				logFile << "${predictionInstance.accession_id} Uploaded hints file ${uploadedStructFile.originalFilename} was renamed to hints.gff and moved to ${projectDir}\n"
				def gffColErrorFlag = 0
				def gffNameErrorFlag = 0
				def gffSourceErrorFlag = 0
				if(!uploadedGenomeFile.empty){ // if seqNames already exists
					// gff format validation: number of columns 9, + or - in column 7, column 1 muss member von seqNames sein
					def gffArray
					def isElement
					new File("${projectDir}/hints.gff").eachLine{line -> 
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
					if(gffSourceErrorFlag == 1){
						logFile << "${predictionInstance.accession_id} Hint file's last column is not in correct format\n"
						flash.error = "Hints file  ${predictionInstance.hint_file} is not in a compatible gff format (the last column does not contain source=M). Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffColErrorFlag == 1){
						logFile << "${predictionInstance.accession_id} Hint file does not always contain 9 columns.\n"
						flash.error = "Hints file  ${predictionInstance.hint_file} is not in a compatible gff format (has not 9 columns). Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffNameErrorFlag == 1){
						logFile << "${predictionInstance.accession_id} Hint file contains entries that do not comply with genome sequence names.\n"
						flash.error = "Entries in the hints file  ${predictionInstance.hint_file} do not match the sequence names of the genome file. Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if((gffColErrorFlag == 1 || gffNameErrorFlag == 1 || gffSourceErrorFlag == 1)){
						logFile << "${predictionInstance.accession_id} ${projectDir} (rm -r) is deleted.\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
						redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
						return
					}
				}
				hintExistsFlag = 1
				def structCksumScript = new File("${projectDir}/struct_cksum.sh")
				def structCksumFile = "${projectDir}/struct.cksum"
				structCksumScript << "cksum ${projectDir}/hints.gff > ${structCksumFile}"
				def structCksumProcess = "bash ${projectDir}/struct_cksum.sh".execute()
				structCksumProcess.waitFor()
				def structCksumContent = new File("${structCksumFile}").text
				def structCksum_array = structCksumContent =~/(\d*) \d* /
				def structCksum
				(1..structCksum_array.groupCount()).each{structCksum = "${structCksum_array[0][it]}"}
				predictionInstance.hint_cksum = "${structCksum}"
				predictionInstance.hint_size = uploadedStructFile.size
				logFile <<  "${predictionInstance.accession_id} struct.fa is ${predictionInstance.hint_size} big and has a cksum of ${structCksum}.\n"
				def delProcCksumStruct = "rm ${projectDir}/struct.cksum".execute()
				delProcCksumStruct.waitFor()
				def delProcCkShStruct = "rm ${projectDir}/struct_cksum.sh".execute()
				delProcCkShStruct.waitFor()
			}
			def radioParameterString
			confirmationString = "${confirmationString}User set UTR prediction: ${predictionInstance.utr}\n"
			// utr
			if(overRideUtrFlag==1){
				radioParameterString = " --UTR=on"
			}else{
				confirmationString = "${confirmationString}Server set UTR prediction: false [UTR parameters missing!]\n"
				radioParameterString = " --UTR=off"
			}
			// strand prediction radio buttons
			if(predictionInstance.pred_strand == 1){
				radioParameterString = "${radioParameterString} --strand=both"
				confirmationString = "${confirmationString}Report genes on: both strands\n"
				logFile << "${predictionInstance.accession_id} User enabled prediction on both strands.\n"
			}else if(predictionInstance.pred_strand == 2){
				confirmationString = "${confirmationString}Report genes on: forward strand only\n"
				radioParameterString = "${radioParameterString} --strand=forward"
				logFile << "${predictionInstance.accession_id} User enabled prediction on forward strand, only.\n"
			}else{
				confirmationString = "${confirmationString}Report genes on: reverse strand only\n"
				radioParameterString = "${radioParameterString} --strand=backward"
				logFile << "${predictionInstance.accession_id} User enabled prediction on reverse strand, only.\n"
			}
			// alternative transcript radio buttons
			if(predictionInstance.alt_transcripts == 1){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=false"
				confirmationString = "${confirmationString}Alternative transcripts: none\n"
				logFile << "${predictionInstance.accession_id} User disabled prediction of alternative transcripts.\n"
			}else if(predictionInstance.alt_transcripts == 2){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 --maxtracks=2"
				confirmationString = "${confirmationString}Alternative transcripts: few\n"
				logFile << "${predictionInstance.accession_id} User enabled prediction of few alternative transcripts.\n"
			}else if(predictionInstance.alt_transcripts == 3){
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.4 --maxtracks=3"
				confirmationString = "${confirmationString}Alternative transcripts: medium\n"
				logFile << "${predictionInstance.accession_id} User enabled prediction of medium alternative transcripts.\n"
			}else{
				radioParameterString = "${radioParameterString} --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=20"
				confirmationString = "${confirmationString}Alternative transcripts: many\n"
				logFile << "${predictionInstance.accession_id} User enabled prediction of many alternative transcripts.\n"
			}
			// gene structure radio buttons
			if(predictionInstance.allowed_structures == 1){
				radioParameterString = "${radioParameterString} --genemodel=partial"
				confirmationString = "${confirmationString}Allowed gene structure: predict any number of (possibly partial) genes\n"
				logFile << "${predictionInstance.accession_id} User enabled the prediction of any number of genes.\n"
			}else if(predictionInstance.allowed_structures == 2){
				radioParameterString = "${radioParameterString} --genemodel=complete"
				confirmationString = "${confirmationString}Allowed gene structure: only predict complete genes\n"
				logFile << "${predictionInstance.accession_id} User disabled the prediction of incomplete genes.\n"
			}else if(predictionInstance.allowed_structures == 3){
				radioParameterString = "${radioParameterString} --genemodel=atleastone"
				confirmationString = "${confirmationString}Allowed gene structure: only predict complete genes - at least one\n"
				logFile << "${predictionInstance.accession_id} User disabled the prediction of incomplete genes and insists on at least one predicted gene.\n"
			}else{
				radioParameterString = "${radioParameterString} --genemodel=exactlyone"
				confirmationString = "${confirmationString}Allowed gene structure: predict exactly one gene\n"
				logFile << "${predictionInstance.accession_id} User enabled the prediction of exactly one gene.\n"
			}
			// ignore gene structure conflicts with other strand checkbox
			if(predictionInstance.ignore_conflicts == false){
				logFile << "${predictionInstance.accession_id} User did not enable to ignore strand conflicts.\n"
			}else{
				radioParameterString = "${radioParameterString} --strand=both"
				logFile << "${predictionInstance.accession_id} User enabled to ignore strand conflicts.\n"
			}
			confirmationString = "${confirmationString}Ignore conflictes with other strand: ${predictionInstance.ignore_conflicts}\n"
			// send confirmation email and redirect
			if(!predictionInstance.hasErrors() && predictionInstance.save()){
				predictionInstance.job_status = 0
				sendMail {
					to "${predictionInstance.email_adress}"
					subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id}"
					body """Hello!

Thank you for submitting the AUGUSTUS gene prediction job ${predictionInstance.accession_id}.

Details of your job:

${confirmationString}
The job status is available at http://bioinf.uni-greifswald.de/trainaugustus/prediction/show/${predictionInstance.id}

You will receive a link to the results via email when the job has finished.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
				}


				logFile << "${predictionInstance.accession_id} Confirmation e-mail sent.\n"          
				redirect(action:show,id:predictionInstance.id)
			} else {
				logFile << "${predictionInstance.accession_id} An error occurred in the predictionInstance (e.g. E-Mail missing, see domain restrictions).\n"
				logFile << "${predictionInstance.accession_id} ${projectDir} is deleted (rm -r).\n"
				def delProc = "rm -r ${projectDir}".execute()
				delProc.waitFor()
				logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
				render(view:'create', model:[predictionInstance:predictionInstance])
				return
			}

			//---------------------  BACKGROUND PROCESS ----------------------------
			Thread.start{
				// retrieve genome file
				if(!(predictionInstance.genome_ftp_link == null)){
					logFile <<  "${predictionInstance.accession_id} Retrieving genome file ${predictionInstance.genome_ftp_link}\n"
					projectDir.mkdirs()
					def wgetGenome = "wget -O ${projectDir}/genome.fa ${predictionInstance.genome_ftp_link}".execute()
					wgetGenome.waitFor()
					if("${predictionInstance.genome_ftp_link}" =~ /\.gz/){
						def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
						gunzipGenomeScript << "cd ${projectDir}; mv genome.fa genome.fa.gz; gunzip genome.fa.gz"
						def gunzipGenome = "bash ${gunzipGenomeScript}".execute()
						gunzipGenome.waitFor()			
						def delProc = "rm ${gunzipGenomeScript}".execute()
						delProc.waitFor()
						logFile <<  "${predictionInstance.accession_id} Unpacked genome file.\n"
					}
					logFile <<  "${predictionInstance.accession_id} genome file upload finished, file stored as genome.fa at ${projectDir}\n"
					// check for fasta format & get seq names for gff validation:
					new File("${projectDir}/genome.fa").eachLine{line -> 
						if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }
						if(line =~ /^>/){
							def len = line.length()
							seqNames << line[1..(len-1)]
						}
					}
					if(genomeFastaFlag == 1) {
						logFile <<  "${predictionInstance.accession_id} The genome file was not fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile <<  "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided genome file ${predictionInstance.genome_ftp_link} was not in DNA fasta format.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""	
						}
						// delete database entry
						predictionInstance.delete()
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
						if(gffColErrorFlag == 1){
							logFile << "${predictionInstance.accession_id} Hints file does not always contain 9 columns.\n"
							sendMail {
								to "${predictionInstance.email_adress}"
								subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
								body """Hello!

Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided hints file ${predictionInstance.hint_file} did not contain 9 columns in each line. Please make sure the gff-format complies with the instructions in our 'Help' section before submitting another job!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus/
"""
							}
						}
						if(gffNameErrorFlag == 1){
							logFile << "${predictionInstance.accession_id} Hints file contains entries that do not comply with genome sequence names.\n"
							sendMail {
								to "${predictionInstance.email_adress}"
								subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
								body """Hello!

Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the sequence names in the provided hints file ${predictionInstance.hint_file} did not comply with the sequence names in the supplied genome file ${predictionInstance.genome_ftp_link}. Please make sure the gff-format complies with the instructions in our 'Help' section before submitting another job!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
							}
						}
						if(gffSourceErrorFlag ==1){
							logFile << "${predictionInstance.accession_id} Hints file contains entries that do not have source=M in the last column.\n"
							sendMail {
								to "${predictionInstance.email_adress}"
								subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
								body """Hello!

Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the last column of your hints file ${predictionInstance.hint_file} does not contain the content source=M. Please make sure the gff-format complies with the instructions in our 'Help' section before submitting another job!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
							}
						
						}
						if((gffColErrorFlag == 1 || gffNameErrorFlag == 1 || gffSourceErrorFlag ==1)){
							logFile << "${predictionInstance.accession_id} ${projectDir} is deleted (rm -r).\n"
							def delProc = "rm -r ${projectDir}".execute()
							delProc.waitFor()
							logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
							// delete database entry
							predictionInstance.delete()
							return
						}
					}
					def genomeCksumScript = new File("${projectDir}/genome_cksum.sh")
					def genomeCksumFile = "${projectDir}/genome.cksum"
					genomeCksumScript << "cksum ${projectDir}/genome.fa > ${genomeCksumFile}"
					def genomeCksumProcess = "bash ${projectDir}/genome_cksum.sh".execute()
					genomeCksumProcess.waitFor()
					def genomeCksumContent = new File("${genomeCksumFile}").text
					def genomeCksum_array = genomeCksumContent =~/(\d*) \d* /
					def genomeCksum
					(1..genomeCksum_array.groupCount()).each{genomeCksum = "${genomeCksum_array[0][it]}"}
					predictionInstance.genome_cksum = "${genomeCksum}"
					genomeCksum_array = genomeCksumContent =~/\d* (\d*) /
					predictionInstance.genome_size
					(1..genomeCksum_array.groupCount()).each{predictionInstance.genome_size = "${genomeCksum_array[0][it]}"}
					logFile <<  "${predictionInstance.accession_id} genome.fa is ${predictionInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
					def delProcCksumGenome = "rm ${projectDir}/genome.cksum".execute()
					delProcCksumGenome.waitFor()
					def delProcCkShGenome = "rm ${projectDir}/genome_cksum.sh".execute()
					delProcCkShGenome.waitFor()
				} // end of if(!(predictionInstance.genome_ftp_link == null))				
				

				// retrieve EST file
				if(!(predictionInstance.est_ftp_link == null)){
					logFile <<  "${predictionInstance.accession_id} Retrieving EST/cDNA file ${predictionInstance.est_ftp_link}\n"
					def wgetEst = "wget -O ${projectDir}/est.fa ${predictionInstance.est_ftp_link}".execute()
					wgetEst.waitFor()
					if("${predictionInstance.est_ftp_link}" =~ /\.gz/){
						def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
						gunzipEstScript << "cd ${projectDir}; mv est.fa est.fa.gz; gunzip est.fa.gz"
						def gunzipEst = "bash ${gunzipEstScript}".execute()
						gunzipEst.waitFor()			
						def delProc = "rm ${gunzipEstScript}".execute()
						delProc.waitFor()
						logFile <<  "${predictionInstance.accession_id} Unpacked EST file.\n"
					}
					logFile <<  "${predictionInstance.accession_id} EST/cDNA file upload finished, file stored as est.fa at ${projectDir}\n"
					// check for fasta format:
					new File("${projectDir}/est.fa").eachLine{line -> if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ estFastaFlag = 1 }}
					if(estFastaFlag == 1) {
						logFile <<  "${predictionInstance.accession_id} The EST/cDNA file was not fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile <<  "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted!\n"
						sendMail {
							to "${predictionInstance.email_adress}"
							subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS prediction job ${predictionInstance.accession_id} was aborted because the provided cDNA file ${predictionInstance.est_ftp_link} was not in DNA fasta format.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
						}
						// delete database entry
						predictionInstance.delete()
						return
					}
					def estCksumScript = new File("${projectDir}/est_cksum.sh")
					def estCksumFile = "${projectDir}/est.cksum"
					estCksumScript << "cksum ${projectDir}/est.fa > ${estCksumFile}"
					def estCksumProcess = "bash ${projectDir}/est_cksum.sh".execute()
					estCksumProcess.waitFor()
					def estCksumContent = new File("${estCksumFile}").text
					def estCksum_array = estCksumContent =~/(\d*) \d* /
					def estCksum
					(1..estCksum_array.groupCount()).each{estCksum = "${estCksum_array[0][it]}"}
					predictionInstance.est_cksum = "${estCksum}"
					estCksum_array = estCksumContent =~/\d* (\d*) /
					predictionInstance.est_size
					(1..estCksum_array.groupCount()).each{predictionInstance.est_size = "${estCksum_array[0][it]}"}
					logFile <<  "${predictionInstance.accession_id} est.fa is ${predictionInstance.est_size} big and has a cksum of ${estCksum}.\n"
					def delProcCksumEst = "rm ${projectDir}/est.cksum".execute()
					delProcCksumEst.waitFor()
					def delProcCkShEst = "rm ${projectDir}/est_cksum.sh".execute()
					delProcCkShEst.waitFor()
				} // end of if(!(predictionInstance.est_ftp_link == null))

				// File formats appear to be ok. 
				// check whether this job was submitted before:
				def grepScript = new File("${projectDir}/grepScript.sh")
				def grepResult = "${projectDir}/grep.result"
				grepScript << "grep \"\\(Genome-Cksum: \\[${predictionInstance.genome_cksum}\\] Genome-Filesize: \\[${predictionInstance.genome_size}\\]\\).*\\(EST-Cksum: \\[${predictionInstance.est_cksum}\\] EST-Filesize: \\[${predictionInstance.est_size}\\]\\).*\\(Hint-Cksum: \\[${predictionInstance.hint_cksum}\\] Hint-Filesize: \\[${predictionInstance.hint_size}\\] Parameter-String: \\[${predictionInstance.project_id}\\]\\).*\\(Parameter-Cksum: \\[${predictionInstance.archive_cksum}\\] Parameter-Size: \\[${predictionInstance.archive_size}\\] Server-Set-UTR-Flag: \\[${overRideUtrFlag}\\]\\).*\\(Report-Genes: \\[${predictionInstance.pred_strand}\\] Alternative-Transcripts: \\[${predictionInstance.alt_transcripts}\\] Gene-Structures: \\[${predictionInstance.allowed_structures}\\] Ignore-Conflicts: \\[${predictionInstance.ignore_conflicts}\\]\\)\" ${dbFile} > ${grepResult}\n"
				def grepJob = "bash ${projectDir}/grepScript.sh".execute()
				grepJob.waitFor()
				def grepContent = new File("${grepResult}").text
				if(grepContent =~ /Genome-Cksum/){
					//job was submitted before. Send E-Mail to user with a link to the results.
					def id_array = grepContent =~ /Grails-ID: \[(\d*)\] /
					oldID
					(0..id_array.groupCount()).each{oldID = "${id_array[0][it]}"}
					def oldAccScript = new File("${projectDir}/oldAcc.sh")
					def oldAccResult = "${projectDir}/oldAcc.result"
					oldAccScript << "grep \"Grails-ID: \\[${oldID}\\]\" ${dbFile} | perl -ne \"@t = split(/\\[/); @t2 = split(/\\]/, \\\$t[4]); print \\\$t2[0];\" > ${oldAccResult}"
					def oldAccScriptProc = "bash ${projectDir}/oldAcc.sh".execute()
					oldAccScriptProc.waitFor()
					def oldAccContent = new File("${oldAccResult}").text	      
					sendMail {
						to "${predictionInstance.email_adress}"
						subject "Your AUGUSTUS prediction job ${predictionInstance.accession_id} was submitted before as job ${oldAccContent}"
						body """Hello!

You submitted job ${predictionInstance.accession_id}. The job was aborted because the files that you submitted were submitted, before. 

The job status of the previously submitted job is available at http://bioinf.uni-greifswald.de/trainaugustus/prediction/show/${oldID}

The results are available at http://bioinf.uni-greifswald.de/trainaugustus/prediction-results/${oldAccContent}/index.html (Results are only available in case the previously submitted job's computations have finished, already.)

Thank you for using AUGUSTUS!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
					}
					logFile << "${predictionInstance.accession_id} Data are identical to old job ${oldAccContent} with Accession-ID ${oldAccContent}. ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile << "${predictionInstance.accession_id} Job ${this_project_id} by user ${predictionInstance.email_adress} is aborted, the user is informed!\n"
					predictionInstance.delete()
					return
				} // end of job was submitted before check

				//Write DB file: 
				dbFile << "Date: [${today}] User-IP: [${userIP}] Grails-ID: [${predictionInstance.id}] Accession-ID: [${predictionInstance.accession_id}] Genome-File: [${predictionInstance.genome_file}] Genome-FTP-Link: [${predictionInstance.genome_ftp_link}] Genome-Cksum: [${predictionInstance.genome_cksum}] Genome-Filesize: [${predictionInstance.genome_size}] EST-File: [${predictionInstance.est_file}] EST-FTP-Link: [${predictionInstance.est_ftp_link}] EST-Cksum: [${predictionInstance.est_cksum}] EST-Filesize: [${predictionInstance.est_size}] Hint-File: [${predictionInstance.hint_file}] Hint-Cksum: [${predictionInstance.hint_cksum}] Hint-Filesize: [${predictionInstance.hint_size}] Parameter-String: [${predictionInstance.project_id}] Parameter-File: [${predictionInstance.archive_file}] Parameter-Cksum: [${predictionInstance.archive_cksum}] Parameter-Size: [${predictionInstance.archive_size}] Server-Set-UTR-Flag: [${overRideUtrFlag}] User-Set-UTR-Flag: [${predictionInstance.utr}] Report-Genes: [${predictionInstance.pred_strand}] Alternative-Transcripts: [${predictionInstance.alt_transcripts}] Gene-Structures: [${predictionInstance.allowed_structures}] Ignore-Conflicts: [${predictionInstance.ignore_conflicts}]\n"

				//rename and move parameters
				if(!uploadedParamArch.empty){
					def mvParamsScript = new File("${projectDir}/mvParams.sh")
					mvParamsScript << "${AUGUSTUS_SCRIPTS_PATH}/moveParameters.pl ${projectDir}/params ${predictionInstance.accession_id} ${AUGUSTUS_CONFIG_PATH}/species 2> ${projectDir}/mvParams.err"
					def mvParamsRunning = "bash ${mvParamsScript}".execute()
					mvParamsRunning.waitFor()
					species = "${predictionInstance.accession_id}"
					logfile << "${predictionInstance.accession_id} Moved uploaded parameters and renamed species to ${predictionInstance.accession_id}\n"
				}
				//Create sge script:
				logFile << "${predictionInstance.accession_id} Writing SGE submission script.\n"
				def sgeFile = new File("${projectDir}/web-aug.sh")
				// write command in script (according to uploaded files)
				sgeFile << "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n\n"
				def cmdStr = "mkdir ${projectDir}/augustus\n"
				if(estExistsFlag == 1){
					cmdStr = "${cmdStr}${BLAT_PATH} -noHead ${projectDir}/genome.fa ${projectDir}/est.fa ${projectDir}/est.psl\n"
					cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/blat2hints.pl --in=${projectDir}/est.psl --out=${projectDir}/est.hints --source=E\n"
					cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/blat2gbrowse.pl ${projectDir}/est.psl ${projectDir}/est.gbrowse\n"
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
				cmdStr = "${cmdStr}${AUGUSTUS_SCRIPTS_PATH}/writeResultsPage.pl ${predictionInstance.accession_id} null ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH}"
				sgeFile << "${cmdStr}"
				

			}
			//------------ END BACKGROUND PROCESS ----------------------------------
		} // end of (!(predictionInstance.id == null))
	}// end of commit
} // end of Controller
