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
   def AUGUSTUS_CONFIG_PATH = "/usr/local/augustus/trunks/config";
   def AUGUSTUS_SCRIPTS_PATH = "/usr/local/augustus/trunks/scripts";
   def scaffold = Prediction
   // the method commit is started if the "Submit Job" button on the website is hit. It is the main method of Prediction Controller and contains a Thread method that will continue running as a background process after the user is redirected to the job status page.

   def commit = {
      def predictionInstance = new Prediction(params)
      if(!(predictionInstance.id == null)){
         redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
         return
      }else{
		// this_project_id that is used internally by pipeline as species name.
		def this_project_id = "web" + predictionInstance.accession_id
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
		// get date
		def today = new Date()
		logFile << "${predictionInstance.accession_id} AUGUSTUS training webserver starting on ${today}\n"
      		// get IP-address
      		String userIP = request.remoteAddr
      		logFile <<  "${predictionInstance.accession_id} user IP: ${userIP}\n"

		// get parameter archive file (if available)
		def uploadedParamArch = request.getFile('ArchiveFile')
		def String dirName = "${output_dir}/${predictionInstance.accession_id}"
		def projectDir = new File(dirName)
		if(!uploadedParamArch.empty){
			// actually upload the file
			projectDir.mkdirs()
         		uploadedParamArch.transferTo( new File (projectDir, "parameters.tar.gz"))
			predictionInstance.archive_file = uploadedArchFile.originalFilename
			logFile <<  "${predictionInstance.accession_id} uploaded parameter archive ${predictionInstance.archive_file} was renamed to parameters.tar.gz and moved to ${projectDir}\n"
			// check whether the archive contains all relevant files
			def String paramDirName = "${projectDir}/params"			
			def paramDir = new File(paramDirName)
			checkParamArch = "checkParamArchive.pl ${projectDir}/parameters.tar.gz paramDirName > ${projectDir}/archCheck.log 2> ${projectDir}/archCheck.err".execute()
			checkParamArch.waitFor()
			def archCheckLog = new File("${projectDir}/archCheck.log")
			def archCheckErr = new File("${projectDir}/archCheck.err")
			def archCheckLogSize = archCheckLog.text.size()
			def archCheckErrSize = archCheckLog.text.size()
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
				overRideUtrFlag = 1 // UTR predictions are now permanently disabled
				logFile <<  "${predictionInstance.accession_id} UTR predictions have been disabled because UTR parameters are missing!\n"
			}			
		}

		// check whether parameters are available for project_id (previous training run)
		if(!(predictionInstance.project_id == null)){

		}


      }
   }

}
