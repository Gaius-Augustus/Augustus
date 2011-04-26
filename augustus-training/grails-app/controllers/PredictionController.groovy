// The class PredictionController controls everything that is related to submitting a job for predicting genes with pre-trained parameters on a novel genome
//    - it handles the file upload (or wget)
//    - format check
//    - SGE job submission and status checks
//    - rendering of results/job status page
//    - sending E-Mails concerning the job status (submission, errors, finished)

class PredictionController {
   // need to adjust the output dir to whatever working dir! This is where uploaded files and results will be saved.
   def output_dir = "/gobics/home/katharina/test"
   // this log File contains the "process log", what was happening with which job when.
   def logFile = new File("${output_dir}/augustus-prediction.log")
   // this log File contains the "database" (not identical with the grails database and simply for logging purpose)
   def dbFile = new File("${output_dir}/augustus-pred-database.log")
   // this log File contains the "database" of the training webserver (must be identical with the file in TrainingController.groovy!)
   def dbTrainFile = new File("${output_dir}/augustus-database.log")

   def scaffold = Prediction
   // the method commit is started if the "Submit Job" button on the website is hit. It is the main method of Prediction Controller and contains a Thread method that will continue running as a background process after the user is redirected to the job status page.
   def commit = {
      def predictionInstance = new Prediction(params)
      if(!(predictionInstance.id == null)){
         redirect(action:create, params:[email_adress:"${predictionInstance.email_adress}"])
         return
      }else{
         // send confirmation email and redirect
         if(!predictionInstance.hasErrors() && predictionInstance.save()){
         predictionInstance.job_status = 0
         sendMail {
            to "${trainingInstance.email_adress}"
            subject "Your AUGUSTUS prediction job ${predictionInstance.project_id}"
            body """Hello!

Thank you for submitting a job to predict genes in a new genome with AUGUSTUS. The job status is available at http://localhost:8080/augustus-training/prediction/show/${predictionInstance.id}

You will be notified by e-mail when the job is finished.

Best regards,

the AUGUSTUS prediction web server team

http://localhost:8080/augustus-training
"""
         }
         redirect(action:show,id:predictionInstance.id)
         }else{
            render(view:'create', model:[predictionInstance:predictionInstance])
            return
         }
      }
   }

}