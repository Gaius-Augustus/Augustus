// Domain file belonging to ../controllers/PredictionController.groovy
// This file contains the variables of all database columns.

class Prediction {
   static searchable = true
   String id 
   static mapping ={
	id generator:'uuid'
   }
   String email_adress
   String project_id // the species name according to a formerly computed webserver training run (random 8 character string)
   String genome_file
   String genome_ftp_link
   String genome_cksum = 0
   String genome_size = 0
   String archive_file
   String archive_cksum = 0
   String archive_size = 0
   String est_file
   String est_ftp_link
   String est_cksum = 0
   String est_size = 0
   String hint_file
   String hint_cksum = 0
   String hint_size = 0
   String job_id // SGE Job ID will be determined by controller
   String job_status // SGE job status will be determined by
   //generate a random (and unique) string for links to results here
  private static String validChars ="ABCDEFGHJKLMNPQRSTUVWXYZ123456789_abcdefghijkmnpqrstuvqxyz"
   private int IDlength=8
   int maxIndex = validChars.length()
   def rnd = new Random()
   String bef_accession_id = (1..IDlength).sum{ 
      validChars[ rnd.nextInt(maxIndex) ] 
   } 
   String accession_id = "pred${bef_accession_id}"
   Date dateCreated
   Boolean utr = false
   Integer pred_strand = 1
   Integer alt_transcripts = 1
   Integer allowed_structures = 1
   String species_select
   String results_urls
   String message
   Boolean ignore_conflicts = false
   String old_url
   // some flags for empty fields in form
   Boolean warn
   Boolean has_genome_file
   Boolean has_est_file
   Boolean has_param_file
   Boolean has_hint_file
   Boolean has_select
   Boolean has_utr
   Boolean has_strand
   Boolean has_transcripts
   Boolean has_structures
   Boolean has_conflicts
   static constraints = {
      //accession_id(unique:true) // may (unlikely) cause problems if the grails database ever gets lost.
      email_adress(email:true,blank:true,nullable:true)
      genome_file(nullable:true, blank:true, validator: { val, obj ->
        if (obj.genome_file == null && obj.genome_ftp_link == null) {
           return 'training.genome_file.no_genome_file'
        } else if (!(obj.genome_ftp_link == null) && !(obj.genome_file == null)) {
           return 'training.genome_file.not_both'
        }else if (obj.genome_ftp_link =~ /dropbox/) {
           return 'prediction.genome_ftp_link.no_dropbox'
        }
	 //else if ((obj.project_id == null) && (obj.archive_file == null)) {
        //   return 'prediction.genome_file.archive_or_id'
        //}
      })
      genome_ftp_link(nullable:true, blank:true, url:true)
      old_url(nullable:true)
      project_id(nullable:true, blank:true, size:3..50)
      species_select(nullable:true)
      est_file(nullable:true, blank:true, validator: { val, obj ->
         if (!(obj.est_file == null) && !(obj.est_ftp_link == null)) {
            return 'training.est_file.not_both'
         } else if (obj.est_ftp_link =~ /dropbox/) {
                return 'prediction.genome_ftp_link.no_dropbox'
            }

      })
      est_ftp_link(nullable:true, blank:true, url:true)
      hint_file(nullable:true, blank:true)
      genome_cksum(nullable:true)
      genome_size(nullable:true)
      est_cksum(nullable:true)
      est_size(nullable:true)
      hint_cksum(nullable:true)
      hint_size(nullable:true)
      archive_cksum(nullable:true)
      archive_size(nullable:true)
      archive_file(nullable:true)
      job_id(nullable:true)
      job_status(nullable:true)
      results_urls(nullable:true)
      message(maxSize:1000000000, nullable:true)
      utr()
      dateCreated()
      has_genome_file(nullable:true)
      has_est_file(nullable:true)
      warn(nullable:true)
      has_param_file(nullable:true)
      has_hint_file(nullable:true)
      has_select(nullable:true)
      has_utr(nullable:true)
      has_strand(nullable:true)
      has_transcripts(nullable:true)
      has_structures(nullable:true)
      has_conflicts(nullable:true)
   }
}
