// The class TrainingController controls everything that is related to submitting a job for training AUGUSTUS  through the webserver:
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


class TrainingController {
	// need to adjust the output dir to whatever working dir! This is where uploaded files and results will be saved.
	def output_dir = "/data/www/augtrain/webdata" // should be something in home of webserver user and augustus frontend user.
	// this log File contains the "process log", what was happening with which job when.
	def logFile = new File("${output_dir}/augustus-training.log")
	// this log File contains the "database" (not identical with the grails database and simply for logging purpose)
	def dbFile = new File("${output_dir}/augustus-database.log")
	// oldID is a parameter that is used for show redirects (see bottom)
	def oldID
	def oldAccID
	// web-output, root directory to the results that are shown to end users
	def web_output_dir = "/var/www/trainaugustus/training-results" // must be writable to webserver application
	// AUGUSTUS_CONFIG_PATH
	def AUGUSTUS_CONFIG_PATH = "/usr/local/augustus/trunks/config";
	def AUGUSTUS_SCRIPTS_PATH = "/usr/local/augustus/trunks/scripts";
	def scaffold = Training
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
			render "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\"/><meta name=\"layout\" content=\"main\" /><title>Submitt Training</title><script type=\"text/javascript\" src=\"js/md_stylechanger.js\"></script></head><body><!-- ***** Start: Kopfbereich ********************************************// --><p class=\"unsichtbar\"><a href=\"#inhalt\" title=\"Directly to Contents\">Directly to Contents</a></p><div id=\"navigation_oben\"><a name=\"seitenanfang\"></a><table width=\"100%\" border=\"0\" cellpadding=\"0\" cellspacing=\"1\"><tr><td nowrap=\"nowrap\"><a href=\"http://www.uni-greifswald.de\" target=\"_blank\" class=\"mainleveltop_\" >University of Greifswald</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.mnf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Faculty</a><span class=\"mainleveltop_\">&nbsp;|&nbsp; </span><a href=\"http://www.math-inf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\" >Institute</a><span class=\"mainleveltop_\">&nbsp;|&nbsp;</span><a href=\"http://bioinf.uni-greifswald.de/\" target=\"_blank\" class=\"mainleveltop_\">Bioinformatics Group</a></td></tr></table></div><div id=\"banner\"><div id=\"banner_links\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php\" title=\"Institut f&uuml;r Mathematik und Informatik\"><img src=\"../images/header.gif\" alt=\"Directly to home\" /> </a></div><div id=\"banner_mitte\"><div id=\"bannertitel1\">Bioinformatics Web Server at University of Greifswald</div><div id=\"bannertitel2\">Gene Prediction with AUGUSTUS</div></div><div id=\"banner_rechts\"><a href=\"http://www.math-inf.uni-greifswald.de/mathe/index.php/geschichte-und-kultur/167\" title=\"Voderberg-Doppelspirale\"><img src=\"../images/spirale.gif\" align=\"left\" /></a></div></div><div id=\"wegweiser\">Navigation for: &nbsp; &nbsp;<span class=\"breadcrumbs pathway\">Submitt Training</span><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Kopfbereich *********************************************// --><!-- ***** Start: Koerper ************************************************// --><div id=\"koerper\"><div id=\"linke_spalte\"><ul class=\"menu\"><li><a href=\"../index.gsp\"><span>Introduction</span></a></li><li><a href=\"/augustus-training/training/create\"><span>Submitt Training</span></a></li><li><a href=\"/augustus-training/prediction/create\"><span>Submitt Prediction</span></a></li><li><a href=\"../help.gsp\"><span>Help</span></a></li><li><a href=\"../references.gsp\"><span>Links & References</span></a></li><li><a href=\"http://bioinf.uni-greifswald.de\"><span>Bioinformatics Group</span></a></li><li><a href=\"http://bioinf.uni-greifswald.de/bioinf/impressum.html\"><span>Impressum</span></a></li></ul></div><div id=\"mittel_spalte\"><div class=\"main\" id=\"main\"><h1><font color=\"#006699\">The Server is Busy</font></h1><p>You tried to access the AUGUSTUS training job submission page.</p><p>Training parameters for gene prediction is a process that takes a lot of computation time. We estimate that one training process requires approximately 10 days. Our web server is able to process a certain number of jobs in parallel, and we established a waiting queue. The waiting queue has a limited length, though. Currently, all slots for computation and for waiting are occupied.</p><p>We apologize for the inconvenience! Please try to submitt your job in a couple of weeks, again.</p><p>Feel free to contact us in case your job is particularly urgent.</p></div><p>&nbsp;</p>           </div><div id=\"rechte_spalte\"><div class=\"linien_div\"><h5 class=\"ueberschrift_spezial\">CONTACT</h5><strong>Institute for Mathematics und Computer Sciences</strong><br/><strong>Bioinformatics Group</strong><br />Walther-Rathenau-Stra&szlig;e 47<br />17487 Greifswald<br />Germany<br />Tel.: +49 (0)3834 86 - 46 24<br/>Fax:  +49 (0)3834 86 - 46 40<br /><br /><a href=\"mailto:augustus-web@uni-greifswald.de\" title=\"E-Mail augustus-web@uni-greifswald.de, opens the standard mail program\">augustus-web@uni-greifswald.de</a></div></div><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Koerper *************************************************// --><!-- ***** Start: Fuss ***************************************************// --><div id=\"fuss\"><div id=\"fuss_links\"><p class=\"copyright\">&copy; 2011 University of Greifswald</p></div><div id=\"fuss_mitte\"><div class=\"bannergroup\"></div></div><div id=\"fuss_rechts\" ><ul><li><a href=\"#seitenanfang\"><img hspace=\"5\" height=\"4\" border=\"0\" width=\"7\" alt=\"Seitenanfang\" src=\"../images/top.gif\" />Top of page</a></li></ul></div><div class=\"beendeFluss\"></div></div><!-- ***** Ende: Fuss ***************************************************// --></body></html>"
			return
		}		
	} 

	// the method commit is started if the "Submit Job" button on the website is hit. It is the main method of Training Controller and contains a Thread method that will continue running as a background process after the user is redirected to the job status page.
	def commit = {
		def trainingInstance = new Training(params)
		if(!(trainingInstance.id == null)){
         		redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
			return
		}else{
			// project_id that is used internally by pipeline as species name.
			def project_id = "web" + trainingInstance.accession_id
			trainingInstance.job_id = 0
			// define flags for file format check, file removal in case of failure
			def genomeFastaFlag = 0
			def estFastaFlag = 0
			def estExistsFlag = 0
			def structureGffFlag = 0
			def structureGbkFlag = 0
			def structureExistsFlag = 0
			def proteinFastaFlag = 0
			def proteinExistsFlag = 0

			// get date
			def today = new Date()
			logFile << "${trainingInstance.accession_id} AUGUSTUS training webserver starting on ${today}\n"
			// get IP-address
			String userIP = request.remoteAddr
			logFile <<  "${trainingInstance.accession_id} user IP: ${userIP}\n"

			//verify that the submitter is a person
			boolean captchaValid = simpleCaptchaService.validateCaptcha(params.captcha)
			if(captchaValid == false){
				logFile << "${trainingInstance.accession_id} The user is probably not a human person. Job aborted.\n"
				flash.error = "The verification string at the bottom of the page was not entered correctly!"
            			redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}"])
           			return
			}

			// upload of genome file
			def uploadedGenomeFile = request.getFile('GenomeFile')
			def seqNames = []
			def String dirName = "${output_dir}/${trainingInstance.accession_id}"
			def projectDir = new File(dirName)
			if(!uploadedGenomeFile.empty){
				projectDir.mkdirs()
				uploadedGenomeFile.transferTo( new File (projectDir, "genome.fa"))
				trainingInstance.genome_file = uploadedGenomeFile.originalFilename
				if("${uploadedGenomeFile.originalFilename}" =~ /\.gz/){
					logFile <<  "${trainingInstance.accession_id} Genome file is gzipped.\n"
					def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
					gunzipGenomeScript << "cd ${projectDir}; mv genome.fa genome.fa.gz; gunzip genome.fa.gz"
					def gunzipGenome = "bash ${gunzipGenomeScript}".execute()
					gunzipGenome.waitFor()
					logFile <<  "${trainingInstance.accession_id} Unpacked genome file.\n"
					def delProc = "rm ${gunzipGenomeScript}".execute()
					delProc.waitFor()
				}
				logFile <<  "${trainingInstance.accession_id} uploaded genome file ${uploadedGenomeFile.originalFilename} was renamed to genome.fa and moved to ${projectDir}\n"
				// check for fasta format & extract fasta headers for gff validation:
				new File("${projectDir}/genome.fa").eachLine{line -> 
					if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }
					if(line =~ /^>/){
						def len = line.length()
						seqNames << line[1..(len-1)]
					}
				}
				if(genomeFastaFlag == 1) {
					logFile <<  "${trainingInstance.accession_id} The genome file was not fasta. Project directory ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile <<  "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
					flash.error = "Genome file ${uploadedGenomeFile.originalFilename} is not in DNA fasta format."
					redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
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
					trainingInstance.genome_cksum = "${genomeCksum}"
					trainingInstance.genome_size = uploadedGenomeFile.size
					logFile <<  "${trainingInstance.accession_id} genome.fa is ${trainingInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
					def delProcCksumGenome = "rm ${projectDir}/genome.cksum".execute()
					delProcCksumGenome.waitFor()
					def delProcCkShGenome = "rm ${projectDir}/genome_cksum.sh".execute()
					delProcCkShGenome.waitFor()
				}	
			}

			// retrieve beginning of genome file for format check
			if(!(trainingInstance.genome_ftp_link == null)){
				logFile <<  "${trainingInstance.accession_id} genome web-link is ${trainingInstance.genome_ftp_link}\n"
				projectDir.mkdirs()
				// checking web file for DNA fasta format: 
				def URL url = new URL("${trainingInstance.genome_ftp_link}");
				def URLConnection uc = url .openConnection()
				if(!("${trainingInstance.genome_ftp_link}" =~ /\.gz/)){
					def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
					def String inputLine=null
					def lineCounter = 1;
					while ( ((inputLine = br.readLine()) != null) && (lineCounter <= 20)) {
						if(!(inputLine =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(inputLine =~ /^$/)){ genomeFastaFlag = 1 }
							lineCounter = lineCounter + 1
					}
					br.close()
					if(genomeFastaFlag == 1) {
						logFile <<  "${trainingInstance.accession_id} The first 20 lines in genome file are not fasta.\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Project directory ${projectDir} is deleted.\n${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						flash.error = "Genome file ${trainingInstance.genome_ftp_link} is not in DNA fasta format."
						redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
						return
					}
				}else{
					logFile <<  "${trainingInstance.accession_id} The linked genome file is gzipped. Format will be checked later after extraction.\n"
				}
			}
 
			// upload of est file
			def uploadedEstFile = request.getFile('EstFile')
			if(!uploadedEstFile.empty){
				projectDir.mkdirs()
				uploadedEstFile.transferTo( new File (projectDir, "est.fa"))
				trainingInstance.est_file = uploadedEstFile.originalFilename
				if("${uploadedEstFile.originalFilename}" =~ /\.gz/){
					logFile <<  "${trainingInstance.accession_id} EST file is gzipped.\n"
					def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
					gunzipEstScript << "cd ${projectDir}; mv est.fa est.fa.gz; gunzip est.fa.gz"
					def gunzipEst = "bash ${gunzipEstScript}".execute()
					gunzipEst.waitFor()
					logFile <<  "${trainingInstance.accession_id} Unpacked EST file.\n"
					def delProc = "rm ${gunzipEstScript}".execute()
					delProc.waitFor()
				}
				logFile << "${trainingInstance.accession_id} Uploaded EST file ${uploadedEstFile.originalFilename} was renamed to est.fa and moved to ${projectDir}\n"
				// check fasta format
				new File("${projectDir}/est.fa").eachLine{line -> 
					if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNnUu]/) && !(line =~ /^$/)){ estFastaFlag = 1 }
				}
				if(estFastaFlag == 1) {
					logFile << "${trainingInstance.accession_id} The cDNA file was not fasta. ${projectDir} (rm -r) is deleted.\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
					flash.error = "cDNA file ${uploadedEstFile.originalFilename} is not in DNA fasta format."
					redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
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
				trainingInstance.est_cksum = "${estCksum}"
				trainingInstance.est_size = uploadedEstFile.size
				logFile <<  "${trainingInstance.accession_id} est.fa is ${trainingInstance.est_size} big and has a cksum of ${estCksum}.\n"
				def delProcCksumEst = "rm ${projectDir}/est.cksum".execute()
				delProcCksumEst.waitFor()
				def delProcCkShEst = "rm ${projectDir}/est_cksum.sh".execute()
				delProcCkShEst.waitFor()
			}

			// retrieve beginning of est file for format check
			if(!(trainingInstance.est_ftp_link == null)){
				logFile << "${trainingInstance.accession_id} est web-link is ${trainingInstance.est_ftp_link}\n"
				projectDir.mkdirs()
				if(!("${trainingInstance.est_ftp_link}" =~ /\.gz/)){
					// checking web file for DNA fasta format: 
					def URL url = new URL("${trainingInstance.est_ftp_link}");
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
						logFile << "${trainingInstance.accession_id} The cDNA file was not fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						flash.error = "cDNA file ${trainingInstance.est_ftp_link} is not in DNA fasta format."
						redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
						return
					}
				}else{
					logFile <<  "${trainingInstance.accession_id} The linked EST file is gzipped. Format will be checked later after extraction.\n"
				}
			}

			// upload of structure file
			def uploadedStructFile = request.getFile('StructFile')
			if(!uploadedStructFile.empty){
				projectDir.mkdirs()
				uploadedStructFile.transferTo( new File (projectDir, "training-gene-structure.gff"))
				trainingInstance.struct_file = uploadedStructFile.originalFilename
				logFile << "${trainingInstance.accession_id} Uploaded training gene structure file ${uploadedStructFile.originalFilename} was renamed to training-gene-structure.gff and moved to ${projectDir}\n"
				def gffColErrorFlag = 0
				def gffNameErrorFlag = 0
				if(!uploadedGenomeFile.empty){ // if seqNames already exists
					// gff format validation: number of columns 9, + or - in column 7, column 1 muss member von seqNames sein
					def gffArray
					def isElement
					new File("${projectDir}/training-gene-structure.gff").eachLine{line -> 
						if(line =~ /^LOCUS/){
							structureGbkFlag = 1 
						}
						if(structureGbkFlag == 0){
							gffArray = line.split("\t")
							if(!(gffArray.size() == 9)){ gffColErrorFlag = 1 }
							isElement = 0
							seqNames.each{ seq ->
								if(seq =~ /${gffArray[0]}/){ isElement = 1 }
								if(isElement == 0){ gffNameErrorFlag = 1 }
							}
						}
					}
					if(gffColErrorFlag == 1 && structureGbkFlag == 0){
						logFile << "${trainingInstance.accession_id} Training gene structure file does not always contain 9 columns.\n"
						flash.error = "Training gene structure file  ${trainingInstance.struct_file} is not in a compatible gff format (has not 9 columns). Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if(gffNameErrorFlag == 1 && structureGbkFlag == 0){
						logFile << "${trainingInstance.accession_id} Training gene structure file contains entries that do not comply with genome sequence names.\n"
						flash.error = "Entries in the training gene structure file  ${trainingInstance.struct_file} do not match the sequence names of the genome file. Please make sure the gff-format complies with the instructions in our 'Help' section!"
					}
					if((gffColErrorFlag == 1 || gffNameErrorFlag == 1) && structureGbkFlag == 0){
						logFile << "${trainingInstance.accession_id} ${projectDir} (rm -r) is deleted.\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
						return
					}
				}
				structureExistsFlag = 1
				def structCksumScript = new File("${projectDir}/struct_cksum.sh")
				def structCksumFile = "${projectDir}/struct.cksum"
				structCksumScript << "cksum ${projectDir}/training-gene-structure.gff > ${structCksumFile}"
				def structCksumProcess = "bash ${projectDir}/struct_cksum.sh".execute()
				structCksumProcess.waitFor()
				def structCksumContent = new File("${structCksumFile}").text
				def structCksum_array = structCksumContent =~/(\d*) \d* /
				def structCksum
				(1..structCksum_array.groupCount()).each{structCksum = "${structCksum_array[0][it]}"}
				trainingInstance.struct_cksum = "${structCksum}"
				trainingInstance.struct_size = uploadedStructFile.size
				logFile <<  "${trainingInstance.accession_id} struct.fa is ${trainingInstance.struct_size} big and has a cksum of ${structCksum}.\n"
				def delProcCksumStruct = "rm ${projectDir}/struct.cksum".execute()
				delProcCksumStruct.waitFor()
				def delProcCkShStruct = "rm ${projectDir}/struct_cksum.sh".execute()
				delProcCkShStruct.waitFor()
			}

			// upload of protein file
			def cRatio = 0
			def uploadedProteinFile = request.getFile('ProteinFile')
			if(!uploadedProteinFile.empty){
				projectDir.mkdirs()
				uploadedProteinFile.transferTo( new File (projectDir, "protein.fa"))
				trainingInstance.protein_file = uploadedProteinFile.originalFilename
				if("${uploadedProteinFile.originalFilename}" =~ /\.gz/){
					logFile <<  "${trainingInstance.accession_id} Protein file is gzipped.\n"
					def gunzipProteinScript = new File("${projectDir}/gunzipProtein.sh")
					gunzipProteinScript << "cd ${projectDir}; mv protein.fa protein.fa.gz; gunzip protein.fa.gz"
					def gunzipProtein = "bash ${gunzipProteinScript}".execute()
					gunzipProtein.waitFor()
					logFile <<  "${trainingInstance.accession_id} Unpacked Protein file.\n"
					def delProc = "rm ${gunzipProteinScript}".execute()
					delProc.waitFor()
				}
				logFile << "${trainingInstance.accession_id} Uploaded protein file ${uploadedProteinFile.originalFilename} was renamed to protein.fa and moved to ${projectDir}\n"
				// check fasta format
				// check that file contains protein sequence, here defined as not more than 5 percent C or c
				def cytosinCounter = 0 // C is cysteine in amino acids, and cytosine in DNA.
				def allAminoAcidsCounter = 0
				new File("${projectDir}/protein.fa").eachLine{line -> 
					if(!(line =~ /^[>AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx ]/) && !(line =~ /^$/)){ proteinFastaFlag = 1 }
					if(!(line =~ /^>/)){
						line.eachMatch(/[AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/){allAminoAcidsCounter = allAminoAcidsCounter + 1}
						line.eachMatch(/[Cc]/){cytosinCounter = cytosinCounter + 1}
					}
				}
				cRatio = cytosinCounter/allAminoAcidsCounter
				if (cRatio >= 0.05){
					logFile << "${trainingInstance.accession_id} The protein file was with cysteine ratio ${cRatio} not recognized as protein file (probably DNA sequence). ${projectDir} is deleted.\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
					flash.error = "Your protein file was not recognized as a protein file. It may be DNA file. The training job was not started. Please contact augustus@uni-greifswald.de if you are completely sure this file is a protein fasta file."
					redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
					return
				}
				if(proteinFastaFlag == 1) {
					logFile << "${trainingInstance.accession_id} The protein file was not protein fasta. ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
					flash.error = "Protein file ${uploadedProteinFile.originalFilename} is not in protein fasta format."
					redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
					return
				}
				proteinExistsFlag = 1
				def proteinCksumScript = new File("${projectDir}/protein_cksum.sh")
				def proteinCksumFile = "${projectDir}/protein.cksum"
				proteinCksumScript << "cksum ${projectDir}/protein.fa > ${proteinCksumFile}"
				def proteinCksumProcess = "bash ${projectDir}/protein_cksum.sh".execute()
				proteinCksumProcess.waitFor()
				def proteinCksumContent = new File("${proteinCksumFile}").text
				def proteinCksum_array = proteinCksumContent =~/(\d*) \d* /
				def proteinCksum
				(1..proteinCksum_array.groupCount()).each{proteinCksum = "${proteinCksum_array[0][it]}"}
				trainingInstance.protein_cksum = "${proteinCksum}"
				trainingInstance.protein_size = uploadedProteinFile.size
				logFile <<  "${trainingInstance.accession_id} protein.fa is ${trainingInstance.protein_size} big and has a cksum of ${proteinCksum}.\n"
				def delProcCksumProtein = "rm ${projectDir}/protein.cksum".execute()
				delProcCksumProtein.waitFor()
				def delProcCkShProtein = "rm ${projectDir}/protein_cksum.sh".execute()
				delProcCkShProtein.waitFor()
			}

			// retrieve beginning of protein file for format check 
			if(!(trainingInstance.protein_ftp_link == null)){
				logFile << "${trainingInstance.accession_id} protein web-link is ${trainingInstance.protein_ftp_link}\n"
				projectDir.mkdirs()
				if(!("${trainingInstance.protein_ftp_link}" =~ /\.gz/)){
					// checking web file for protein fasta format: 
					def URL url = new URL("${trainingInstance.protein_ftp_link}");
					def URLConnection uc = url .openConnection()
					def BufferedReader br = new BufferedReader(new InputStreamReader(uc.getInputStream()))
					def String inputLine=null
					def lineCounter = 1;
					def cytosinCounter = 0 // C is cysteine in amino acids, and cytosine in DNA.
					def allAminoAcidsCounter = 0
					while ( ((inputLine = br.readLine()) != null) && (lineCounter <= 50)) {
						if(!(inputLine =~ /^[>AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/) && !(inputLine =~ /^$/)){ proteinFastaFlag = 1 }
						if(!(inputLine =~ /^>/)){
							inputLine.eachMatch(/[AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/){allAminoAcidsCounter = allAminoAcidsCounter + 1}
							inputLine.eachMatch(/[Cc]/){cytosinCounter = cytosinCounter + 1}
						}
					}
					br.close()
					cRatio = cytosinCounter/allAminoAcidsCounter
					if (cRatio >= 0.05){
						logFile << "${trainingInstance.accession_id} The protein file was with cysteine ratio ${cRatio} not recognized as protein file (probably DNA sequence).\n ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "\n${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						flash.error = "Protein file ${trainingInstance.protein_ftp_link} does not contain protein sequences."
						redirect(action:create)
						return
					}
					if(proteinFastaFlag == 1) {
						logFile << "${trainingInstance.accession_id} The protein file was not protein fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						flash.message = "Protein file ${trainingInstance.protein_ftp_link} is not in protein fasta format."
						redirect(action:create, params:[email_adress:"${trainingInstance.email_adress}", project_name:"${trainingInstance.project_name}"])
						return
					}
				}else{
					logFile <<  "${trainingInstance.accession_id} The linked Protein file is gzipped. Format will be checked later after extraction.\n"
				}
			}

			// send confirmation email and redirect
			if(!trainingInstance.hasErrors() && trainingInstance.save()){
				trainingInstance.job_status = 0
				sendMail {
					to "${trainingInstance.email_adress}"
					subject "Your AUGUSTUS training job ${trainingInstance.accession_id}"
					body """Hello!

Thank you for submitting a job to train AUGUSTUS for the species ${trainingInstance.project_name}. The job status is available at http://bioinf.uni-greifswald.de/trainaugustus/training/show/${trainingInstance.id}

You will be notified by e-mail when the job is finished.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
				}
				logFile << "${trainingInstance.accession_id} Confirmation e-mail sent.\n"          
				redirect(action:show,id:trainingInstance.id)
			} else {
				logFile << "${trainingInstance.accession_id} An error occurred in the trainingInstance (e.g. E-Mail missing, see domain restrictions).\n"
				logFile << "${trainingInstance.accession_id} ${projectDir} is deleted (rm -r).\n"
				def delProc = "rm -r ${projectDir}".execute()
				delProc.waitFor()
				logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
				render(view:'create', model:[trainingInstance:trainingInstance])
				return
			}

			//---------------------  BACKGROUND PROCESS ----------------------------
			Thread.start{
				// retrieve genome file
				if(!(trainingInstance.genome_ftp_link == null)){
					logFile <<  "${trainingInstance.accession_id} Retrieving genome file ${trainingInstance.genome_ftp_link}\n"
					projectDir.mkdirs()
					def wgetGenome = "wget -O ${projectDir}/genome.fa ${trainingInstance.genome_ftp_link}".execute()
					wgetGenome.waitFor()
					if("${trainingInstance.genome_ftp_link}" =~ /\.gz/){
						def gunzipGenomeScript = new File("${projectDir}/gunzipGenome.sh")
						gunzipGenomeScript << "cd ${projectDir}; mv genome.fa genome.fa.gz; gunzip genome.fa.gz"
						def gunzipGenome = "bash ${gunzipGenomeScript}".execute()
						gunzipGenome.waitFor()			
						def delProc = "rm ${gunzipGenomeScript}".execute()
						delProc.waitFor()
						logFile <<  "${trainingInstance.accession_id} Unpacked genome file.\n"
					}
					logFile <<  "${trainingInstance.accession_id} genome file upload finished, file stored as genome.fa at ${projectDir}\n"
					// check for fasta format & get seq names for gff validation:
					new File("${projectDir}/genome.fa").eachLine{line -> 
						if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ genomeFastaFlag = 1 }
						if(line =~ /^>/){
							def len = line.length()
							seqNames << line[1..(len-1)]
						}
					}
					if(genomeFastaFlag == 1) {
						logFile <<  "${trainingInstance.accession_id} The genome file was not fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile <<  "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						sendMail {
							to "${trainingInstance.email_adress}"
							subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the provided genome file ${trainingInstance.genome_ftp_link} was not in DNA fasta format.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""	
						}
						// delete database entry
						trainingInstance.delete()
						return
					}
					// check gff format
					def gffColErrorFlag = 0
					def gffNameErrorFlag = 0
					if((!uploadedStructFile.empty) &&(!(trainingInstance.genome_ftp_link == null))){ // if seqNames already exists
						// gff format validation: number of columns 9, + or - in column 7, column 1 muss member von seqNames sein
						def gffArray
						def isElement
						new File("${projectDir}/training-gene-structure.gff").eachLine{line -> 
							if(line =~ /^LOCUS/){
								structureGbkFlag = 1 
							}
							if(structureGbkFlag == 0){
								gffArray = line.split("\t")
								if(!(gffArray.size() == 9)){ gffColErrorFlag = 1 }
								isElement = 0
								seqNames.each{ seq ->
									if(seq =~ /${gffArray[0]}/){ isElement = 1 }
									if(isElement == 0){ gffNameErrorFlag = 1 }
								}
							}
						}
						if(gffColErrorFlag == 1 && structureGbkFlag == 0){
							logFile << "${trainingInstance.accession_id} Training gene structure file does not always contain 9 columns.\n"
							sendMail {
								to "${trainingInstance.email_adress}"
								subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
								body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the provided training gene structure file ${trainingInstance.struct_file} did not contain 9 columns in each line. Please make sure the gff-format complies with the instructions in our 'Help' section before submitting another job!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus/
"""
							}
						}
						if(gffNameErrorFlag == 1 && structureGbkFlag == 0){
							logFile << "${trainingInstance.accession_id} Training gene structure file contains entries that do not comply with genome sequence names.\n"
							sendMail {
								to "${trainingInstance.email_adress}"
								subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
								body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the sequence names in the provided training gene structure file ${trainingInstance.struct_file} did not comply with the sequence names in the supplied genome file ${trainingInstance.genome_ftp_link}. Please make sure the gff-format complies with the instructions in our 'Help' section before submitting another job!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
							}
						}
						if((gffColErrorFlag == 1 || gffNameErrorFlag == 1) && structureGbkFlag == 0){
							logFile << "${trainingInstance.accession_id} ${projectDir} is deleted (rm -r).\n"
							def delProc = "rm -r ${projectDir}".execute()
							delProc.waitFor()
							logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
							// delete database entry
							trainingInstance.delete()
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
					trainingInstance.genome_cksum = "${genomeCksum}"
					genomeCksum_array = genomeCksumContent =~/\d* (\d*) /
					trainingInstance.genome_size
					(1..genomeCksum_array.groupCount()).each{trainingInstance.genome_size = "${genomeCksum_array[0][it]}"}
					logFile <<  "${trainingInstance.accession_id} genome.fa is ${trainingInstance.genome_size} big and has a cksum of ${genomeCksum}.\n"
					def delProcCksumGenome = "rm ${projectDir}/genome.cksum".execute()
					delProcCksumGenome.waitFor()
					def delProcCkShGenome = "rm ${projectDir}/genome_cksum.sh".execute()
					delProcCkShGenome.waitFor()
				} // end of if(!(trainingInstance.genome_ftp_link == null))

				// retrieve EST file
				if(!(trainingInstance.est_ftp_link == null)){
					logFile <<  "${trainingInstance.accession_id} Retrieving EST/cDNA file ${trainingInstance.est_ftp_link}\n"
					def wgetEst = "wget -O ${projectDir}/est.fa ${trainingInstance.est_ftp_link}".execute()
					wgetEst.waitFor()
					if("${trainingInstance.est_ftp_link}" =~ /\.gz/){
						def gunzipEstScript = new File("${projectDir}/gunzipEst.sh")
						gunzipEstScript << "cd ${projectDir}; mv est.fa est.fa.gz; gunzip est.fa.gz"
						def gunzipEst = "bash ${gunzipEstScript}".execute()
						gunzipEst.waitFor()			
						def delProc = "rm ${gunzipEstScript}".execute()
						delProc.waitFor()
						logFile <<  "${trainingInstance.accession_id} Unpacked EST file.\n"
					}					
					logFile <<  "${trainingInstance.accession_id} EST/cDNA file upload finished, file stored as est.fa at ${projectDir}\n"
					// check for fasta format:
					new File("${projectDir}/est.fa").eachLine{line -> if(!(line =~ /^[>AaTtGgCcHhXxRrYyWwSsMmKkBbVvDdNn]/) && !(line =~ /^$/)){ estFastaFlag = 1 }}
					if(estFastaFlag == 1) {
						logFile <<  "${trainingInstance.accession_id} The EST/cDNA file was not fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile <<  "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						sendMail {
							to "${trainingInstance.email_adress}"
							subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the provided cDNA file ${trainingInstance.est_ftp_link} was not in DNA fasta format.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
						}
						// delete database entry
						trainingInstance.delete()
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
					trainingInstance.est_cksum = "${estCksum}"
					estCksum_array = estCksumContent =~/\d* (\d*) /
					trainingInstance.est_size
					(1..estCksum_array.groupCount()).each{trainingInstance.est_size = "${estCksum_array[0][it]}"}
					logFile <<  "${trainingInstance.accession_id} est.fa is ${trainingInstance.est_size} big and has a cksum of ${estCksum}.\n"
					def delProcCksumEst = "rm ${projectDir}/est.cksum".execute()
					delProcCksumEst.waitFor()
					def delProcCkShEst = "rm ${projectDir}/est_cksum.sh".execute()
					delProcCkShEst.waitFor()
				} // end of if(!(trainingInstance.est_ftp_link == null))

				// retrieve protein file
				if(!(trainingInstance.protein_ftp_link == null)){
					logFile <<  "${trainingInstance.accession_id} Retrieving protein file ${trainingInstance.protein_ftp_link}\n"
					def wgetProtein = "wget -O ${projectDir}/protein.fa ${trainingInstance.protein_ftp_link}".execute()
					wgetProtein.waitFor()
					if("${trainingInstance.protein_ftp_link}" =~ /\.gz/){
						def gunzipProteinScript = new File("${projectDir}/gunzipProtein.sh")
						gunzipProteinScript << "cd ${projectDir}; mv protein.fa protein.fa.gz; gunzip protein.fa.gz"
						def gunzipProtein = "bash ${gunzipProteinScript}".execute()
						gunzipProtein.waitFor()			
						def delProc = "rm ${gunzipProteinScript}".execute()
						delProc.waitFor()
						logFile <<  "${trainingInstance.accession_id} Unpacked protein file.\n"
					}
					logFile <<  "${trainingInstance.accession_id} Protein file upload finished, file stored as protein.fa at ${projectDir}\n"
					// check for fasta protein format:
					def cytosinCounter = 0 // C is cysteine in amino acids, and cytosine in DNA.
					def allAminoAcidsCounter = 0
					new File("${projectDir}/protein.fa").eachLine{line -> 
						if(!(line =~ /^[>AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx ]/) && !(line =~ /^$/)){ proteinFastaFlag = 1 }
						if(!(line =~ /^>/)){
							line.eachMatch(/[AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/){ allAminoAcidsCounter = allAminoAcidsCounter + 1 }
							line.eachMatch(/[Cc]/){ cytosinCounter = cytosinCounter + 1 }
						}
					}
					cRatio = cytosinCounter/allAminoAcidsCounter
					if (cRatio >= 0.05){
						logFile << "${trainingInstance.accession_id} The protein file was with cysteine ratio ${cRatio} not recognized as protein file (probably DNA sequence). ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						sendMail {
							to "${trainingInstance.email_adress}"
							subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the provided protein file ${trainingInstance.protein_ftp_link} is suspected to contain DNA instead of protein sequences.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""	
						}
						// delete database entry
						trainingInstance.delete()
						return
					}
					if(proteinFastaFlag == 1) {
						logFile << "${trainingInstance.accession_id} The protein file was not protein fasta. ${projectDir} is deleted (rm -r).\n"
						def delProc = "rm -r ${projectDir}".execute()
						delProc.waitFor()
						logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted!\n"
						sendMail {
							to "${trainingInstance.email_adress}"
							subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was aborted"
							body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} was aborted because the provided protein file ${trainingInstance.protein_ftp_link} is not in fasta format.

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
						}
						// delete database entry
						trainingInstance.delete()
						return
					}
					def proteinCksumScript = new File("${projectDir}/protein_cksum.sh")
					def proteinCksumFile = "${projectDir}/protein.cksum"
					proteinCksumScript << "cksum ${projectDir}/protein.fa > ${proteinCksumFile}"
					def proteinCksumProcess = "bash ${projectDir}/protein_cksum.sh".execute()
					proteinCksumProcess.waitFor()
					def proteinCksumContent = new File("${proteinCksumFile}").text
					def proteinCksum_array = proteinCksumContent =~/(\d*) \d* /
					def proteinCksum
					(1..proteinCksum_array.groupCount()).each{proteinCksum = "${proteinCksum_array[0][it]}"}
					trainingInstance.protein_cksum = "${proteinCksum}"
					proteinCksum_array = proteinCksumContent =~/\d* (\d*) /
					trainingInstance.protein_size
					(1..proteinCksum_array.groupCount()).each{trainingInstance.protein_size = "${proteinCksum_array[0][it]}"}
					logFile <<  "${trainingInstance.accession_id} protein.fa is ${trainingInstance.protein_size} big and has a cksum of ${proteinCksum}.\n"
					def delProcCksumProtein = "rm ${projectDir}/protein.cksum".execute()
					delProcCksumProtein.waitFor()
					def delProcCkShProtein = "rm ${projectDir}/protein_cksum.sh".execute()
					delProcCkShProtein.waitFor()
				} // end of (!(trainingInstance.protein_ftp_link == null))

				// File formats appear to be ok. 
				// check whether this job was submitted before:
				def grepScript = new File("${projectDir}/grepScript.sh")
				def grepResult = "${projectDir}/grep.result"
				grepScript << "grep \"\\(Genome-Cksum: \\[${trainingInstance.genome_cksum}\\] Genome-Filesize: \\[${trainingInstance.genome_size}\\]\\).*\\(EST-Cksum: \\[${trainingInstance.est_cksum}\\] EST-Filesize: \\[${trainingInstance.est_size}\\]\\).*\\(Protein-Cksum: \\[${trainingInstance.protein_cksum}\\] Protein-Filesize: \\[${trainingInstance.protein_size}\\]\\)*\\(Struct-Cksum: \\[${trainingInstance.struct_cksum}\\]\\)\" ${dbFile} > ${grepResult}\n"
				def grepJob = "bash ${projectDir}/grepScript.sh".execute()
				grepJob.waitFor()
				def grepContent = new File("${grepResult}").text
				if(grepContent =~ /Struct-Cksum/){
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
						to "${trainingInstance.email_adress}"
						subject "Your AUGUSTUS training job ${trainingInstance.accession_id} was submitted before as job ${oldAccContent}"
						body """Hello!

You submitted job ${trainingInstance.accession_id} for species ${trainingInstance.project_name}. The job was aborted because the files that you submitted were submitted, before. 

The job status of the previously submitted job is available at http://bioinf.uni-greifswald.de/trainaugustus/training/show/${oldID}

The results are available at http://bioinf.uni-greifswald.de/trainaugustus/training-results/${oldAccContent}/index.html (Results are only available in case the previously submitted job's computations have finished, already.)

Thank you for using AUGUSTUS!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
					}
					logFile << "${trainingInstance.accession_id} Data are identical to old job ${oldAccContent} with Accession-ID ${oldAccContent}. ${projectDir} is deleted (rm -r).\n"
					def delProc = "rm -r ${projectDir}".execute()
					delProc.waitFor()
					logFile << "${trainingInstance.accession_id} Job ${project_id} by user ${trainingInstance.email_adress} is aborted, the user is informed!\n"
					trainingInstance.delete()
					return
				} // end of job was submitted before check
				//Write DB file: 
				dbFile << "Date: [${today}] User-IP: [${userIP}] Grails-ID: [${trainingInstance.id}] Accession-ID: [${trainingInstance.accession_id}] Genome-File: [${trainingInstance.genome_file}] Genome-FTP-Link: [${trainingInstance.genome_ftp_link}] Genome-Cksum: [${trainingInstance.genome_cksum}] Genome-Filesize: [${trainingInstance.genome_size}] EST-File: [${trainingInstance.est_file}] EST-FTP-Link: [${trainingInstance.est_ftp_link}] EST-Cksum: [${trainingInstance.est_cksum}] EST-Filesize: [${trainingInstance.est_size}] Protein-File: [${trainingInstance.protein_file}] Protein-FTP-Link: [${trainingInstance.protein_ftp_link}] Protein-Cksum: [${trainingInstance.protein_cksum}] Protein-Filesize: [${trainingInstance.protein_size}] Training-Structure-File: [${trainingInstance.struct_file}] Struct-Cksum: [${trainingInstance.struct_cksum}] Struct-Filesize: [${trainingInstance.struct_size}]\n"
				//Create a test sge script:
				logFile << "${trainingInstance.accession_id} Writing SGE submission script.\n"
				def sgeFile = new File("${projectDir}/web-aug.sh")
				// write command in script (according to uploaded files)
				sgeFile << "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n\n"
				// this has been checked, works.
				if( estExistsFlag ==1 && proteinExistsFlag == 0 && structureExistsFlag == 0){
					sgeFile << "autoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --cdna=${projectDir}/est.fa --pasa -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err\n"
					// this is currently tested
				}else if(estExistsFlag == 0 && proteinExistsFlag == 0 && structureExistsFlag == 1){
					sgeFile << "autoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --trainingset=${projectDir}/training-gene-structure.gff -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err"
					// this is currently tested
				}else if(estExistsFlag == 0 && proteinExistsFlag == 1 && structureExistsFlag == 0){
					sgeFile << "autoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --trainingset=${projectDir}/protein.fa -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err\n"
					// all following commands still need testing
				}else if(estExistsFlag == 1 && proteinExistsFlag == 1 && structureExistsFlag == 0){
					sgeFile << "autoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --cdna=${projectDir}/est.fa --trainingset=${projectDir}/protein.fa -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err\n\n"
				}else if(estExistsFlag == 1 && proteinExistsFlag == 0 && structureExistsFlag == 1){
					sgeFile << "autoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --cdna=${projectDir}/est.fa --trainingset=${projectDir}/training-gene-structure.gff -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err"
				}else if(estExistsFlag == 0 && proteinExistsFlag == 1 && structureExistsFlag == 1){
					sgeFile << "echo 'Simultaneous protein and structure file support are currently not implemented. Using the structure file, only.'\n\nautoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --trainingset=${projectDir}/training-gene-structure.gff -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err"
				}else if(estExistsFlag == 1 && proteinExistsFlag == 1 && structureExistsFlag == 1){
					sgeFile << "echo Simultaneous protein and structure file support are currently not implemented.\n\nUsing the structure file, only.'\n\nautoAug.pl --genome=${projectDir}/genome.fa --species=${trainingInstance.accession_id} --trainingset=${projectDir}/training-gene-structure.gff -v --singleCPU --workingdir=${projectDir} > ${projectDir}/AutoAug.log 2> ${projectDir}/AutoAug.err\n\nwriteResultsPage.pl ${trainingInstance.accession_id} ${trainingInstance.project_name} ${dbFile} ${output_dir} ${web_output_dir} ${AUGUSTUS_CONFIG_PATH} ${AUGUSTUS_SCRIPTS_PATH} > ${projectDir}/writeResults.log 2> ${projectDir}/writeResults.err"
				}
				// write submission script
				def submissionScript = new File("${projectDir}/submitt.sh")
				def fileID = "${projectDir}/jobID"
				submissionScript << "cd ${projectDir}; qsub web-aug.sh > ${fileID}"
				// submitt job
				def jobSubmission = "bash ${projectDir}/submitt.sh".execute()
				jobSubmission.waitFor()
				// get job ID
				def content = new File("${fileID}").text
				def jobID_array = content =~/Your job (\d*)/
				def jobID
				(1..jobID_array.groupCount()).each{jobID = "${jobID_array[0][it]}"}
				trainingInstance.job_id = jobID
				logFile << "${trainingInstance.accession_id} Job ${jobID} submitted.\n"
				// check for job status
				trainingInstance.job_status = 1 // submitted
				trainingInstance = trainingInstance.merge()
				trainingInstance.save()
				def statusScript = new File("${projectDir}/status.sh")
				def statusFile = "${projectDir}/job.status"
				statusScript << "cd ${projectDir}; qstat|grep web-aug.sh |grep ${jobID} > ${statusFile}"
				def statusContent
				def statusCheck 
				def qstat = 1
				def runFlag = 0;

				while(qstat == 1){
					statusCheck = "bash ${projectDir}/status.sh".execute()
					statusCheck.waitFor()
					statusContent = new File("${statusFile}").text
					if(statusContent =~ /qw/){ 
						trainingInstance.job_status = 2 
						trainingInstance = trainingInstance.merge()
						trainingInstance.save()
					}else if( statusContent =~ /  r  / ){
						trainingInstance.job_status = 3
						trainingInstance = trainingInstance.merge()
						trainingInstance.save()
						if(runFlag == 0){
							today = new Date() 
							logFile << "${trainingInstance.accession_id} Job ${jobID} begins running at ${today}.\n"
						}
						runFlag = 1
					}else{
						trainingInstance.job_status = 4
						trainingInstance = trainingInstance.merge()
						trainingInstance.save()
						qstat = 0
						today = new Date()
 						logFile << "${trainingInstance.accession_id} Job ${jobID} left SGE at ${today}.\n"
					}
					sleep 5000
			   	}

			   	// check whether errors occured by log-file-sizes
			   	def autoAugErrFile = new File("${projectDir}/AutoAug.err")
				def sgeErrFile = new File("${projectDir}/web-aug.sh.e${jobID}")
				def writeResultsErrFile = new File("${projectDir}/writeResults.err")
				def autoAugErrSize = autoAugErrFile.text.size()
				def sgeErrSize = sgeErrFile.text.size()
				def writeResultsErrSize = writeResultsErrFile.text.size()
				if(autoAugErrSize==0 && sgeErrSize==0 && writeResultsErrSize==0){
					sendMail {
						to "${trainingInstance.email_adress}"
						subject "Your AUGUSTUS training job ${trainingInstance.accession_id} is complete"
						body """Hello!

Your AUGUSTUS training job ${trainingInstance.accession_id} for species ${trainingInstance.project_name} is complete. You find the results at http://bioinf.uni-greifswald.de/trainaugustus/training-results/${trainingInstance.accession_id}/index.html .

Thank you for using AUGUSTUS!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
					}
					logFile << "${trainingInstance.accession_id} Sent confirmation Mail that job computation was successful.\n"
					// unpack with 7z x XA2Y5VMJ.tar.7z
					// tar xvf XA2Y5VMJ.tar
					def packResults = new File("${output_dir}/pack${trainingInstance.accession_id}.sh")
					packResults << "cd ${output_dir}; tar cf - ${trainingInstance.accession_id} | 7z a -si ${trainingInstance.accession_id}.tar.7z; rm -r ${trainingInstance.accession_id};"
					def cleanUp = "bash ${output_dir}/pack${trainingInstance.accession_id}.sh".execute()
					cleanUp.waitFor()
					def cleanUp2 = "rm ${output_dir}/pack${trainingInstance.accession_id}.sh".execute()
					cleanUp2.waitFor()
					logFile << "${trainingInstance.accession_id} autoAug directory was packed with tar/7z.\n"
					logFile << "${trainingInstance.accession_id} Job completed. Result: ok.\n"
				}else{
					if(autoAugErrSize > 0){
						logFile << "${trainingInstance.accession_id} an error occured when autoAug.pl was executed!\n"; 
						sendMail {
						to "${admin_email}"
						subject "Error in AUGUSTUS training job ${trainingInstance.accession_id}"
						body """Hi ${admin_email}!

Job: ${trainingInstance.accession_id}
E-Mail: ${trainingInstance.email_adress}
Species: ${trainingInstance.project_name}
Link: http://bioinf.uni-greifswald.de/trainaugustus/training-results/${trainingInstance.accession_id}/index.html

An error occured in the autoAug pipeline. Please check manually what's wrong. The user has been informed.
"""
						}
						trainingInstance.job_status = 5
					}else if(sgeErrFile > 0){
						logFile << "${trainingInstance.accession_id} a SGE error occured!\n";
						sendMail {
						to "${admin_email}"
						subject "Error in AUGUSTUS training job ${trainingInstance.accession_id}"
						body """Hi ${admin_email}!

Job: ${trainingInstance.accession_id}
E-Mail: ${trainingInstance.email_adress}
Species: ${trainingInstance.project_name}
Link: http://bioinf.uni-greifswald.de/trainaugustus/training-results/${trainingInstance.accession_id}/index.html

An SGE error occured. Please check manually what's wrong. The user has been informed.
"""
						}
						trainingInstance.job_status = 5
					}else{
						logFile << "${trainingInstance.accession_id} an error occured during writing results!\n";
						sendMail {
							to "${admin_email}"
							subject "Error in AUGUSTUS training job ${trainingInstance.accession_id}"
							body """Hi ${admin_email}!

Job: ${trainingInstance.accession_id}
E-Mail: ${trainingInstance.email_adress}
Species: ${trainingInstance.project_name}
Link: http://bioinf.uni-greifswald.de/trainaugustus/training-results/${trainingInstance.accession_id}/index.html

An error occured during writing results. Please check manually what's wrong. The user has been informed.
"""
						}
					}
					sendMail {
						to "${trainingInstance.email_adress}"
						subject "Your AUGUSTUS training job ${trainingInstance.accession_id} is in an error state"
						body """Hello!

An error occured while training AUGUSTUS for species ${trainingInstance.project_name} (Job ${trainingInstance.accession_id}).

The administrator of the AUGUSTUS training web server has been informed and will get back to you as soon as the problem is solved.

Thank you for using AUGUSTUS!

Best regards,

the AUGUSTUS training web server team

http://bioinf.uni-greifswald.de/trainaugustus
"""
					}

					logFile << "${trainingInstance.accession_id} Sent confirmation Mail, the job is in an error state.\n"
				}
			}
			//------------ END BACKGROUND PROCESS ----------------------------------
		}
	}
}
