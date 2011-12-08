// Domain file belonging to ../controllers/TrainingController.groovy
// This file contains the variables of all database columns.

class Training {
   static searchable = true
   
   String email_adress
   String project_name
   String genome_file
   String genome_ftp_link
   String est_file
   String est_ftp_link
   String job_id // SGE Job ID will be determined by controller
   String job_status // SGE job status will be determined by controller
   String protein_file
   String protein_ftp_link
   String struct_file
   //generate a random (and unique) string for links to results here
   private static String validChars ="ABCDEFGHJKLMNPQRSTUVWXYZ123456789_abcdefghijkmnpqrstuvqxyz"
   private int IDlength=8
   int maxIndex = validChars.length()
   def rnd = new Random()
   String bef_accession_id = (1..IDlength).sum{ 
      validChars[ rnd.nextInt(maxIndex) ] 
   } 
   String accession_id = "train${bef_accession_id}"

   // checksum columns in database, values will be determined by controller
   String genome_cksum = 0
   String est_cksum = 0
   String protein_cksum = 0
   String struct_cksum = 0
   String genome_size = 0
   String est_size = 0
   String protein_size = 0
   String struct_size = 0
   Date dateCreated

   static constraints ={
       accession_id(unique:true) // may (unlikely) cause problems if the grails database ever gets lost.
       email_adress(email:true,blank:false,nullable:false)
       project_name(blank:false, unique:true, maxSize:30)
       genome_file(nullable:true, blank:true, validator: { val, obj ->
            if (obj.genome_file == null && obj.genome_ftp_link == null) {
                return 'training.genome_file.no_genome_file'
            } else if (!(obj.genome_ftp_link == null) && !(obj.genome_file == null)) {
                return 'training.genome_file.not_both'
            } else if ((obj.est_file == null) && (obj.est_ftp_link == null) && (obj.struct_file == null) && (obj.protein_file == null) && (obj.protein_ftp_link == null)) {
                return 'training.genome_file.at_least_one'
            } else if (!(obj.protein_file == null) && !(obj.struct_file == null)) {
                return 'training.protein_file.not_struct'
            } else if (!(obj.protein_ftp_link == null) && !(obj.struct_file == null)) {
                return 'training.protein_file.not_struct'
            }
       })
       genome_ftp_link(nullable:true, blank:true, url:true)
       est_file(nullable:true, blank:true, validator: { val, obj ->
           if (!(obj.est_file == null) && !(obj.est_ftp_link == null)) {
                return 'training.est_file.not_both'
           }
       })
       est_ftp_link(nullable:true, blank:true, url:true)
       protein_file(nullable:true, blank:true, validator: { val, obj ->
           if (!(obj.protein_file == null) && !(obj.protein_ftp_link == null)) {
                return 'training.protein_file.not_both'
           }
       })
       struct_file(nullable:true, blank:true)
       protein_ftp_link(nullable:true, blank:true, url:true)
       genome_cksum(nullable:true)
       est_cksum(nullable:true)
       protein_cksum(nullable:true)
       struct_cksum(nullable:true)
       genome_size(nullable:true)
       est_size(nullable:true)
       protein_size(nullable:true)
       struct_size(nullable:true)
       job_id(nullable:true)
       job_status(nullable:true)
       dateCreated()
    }
}
