#include "jg_transcript.h"
#include "jg_ios.h"

using namespace std;

int temp_gene_number = 0;			// needed for the adjusted case (to test); renames the genes from g0 to gn; if n+1 is the number of input transcripts

void load_error(string const &error)
{
    cerr << "Load error: " << error << endl;
    cerr << "The file is not in a correct gff or gtf format." << endl;
    exit( EXIT_FAILURE );
}

void load_warning(string const &warning)
{
    cerr << "Load warning: " << warning << endl;
    cerr << "This warning may affect the result." << endl;
}

void load(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename, int &priority)
{
    // loads gtf file in a special data structure
    // mainly there is a list where all transcripts are saved; later this programm work on pointers to this list elements
    // every transcript gets a pointer to his gene in the gene hash map and backwards
    // this function will not distinguish between first column notation "chr1" and "1" or another 
    // this function doesnt tests whether there are transcripts with equal name; transcripts with equal name in a row will be read as one transcript; a divided transcript may be read as different transcripts
    // if their is a notation for prediction ranges BEFORE a transcript in the form "#Xprediction on sequence rangeX:NYNYZ" it will save this for the transcript where "X" can be every char exept ":", "Y" have to be a char out of " ", "-", "(", ")", "b", "p", and "Z" can be every char and the number of X,Y and Z is arbitrary at every position and "N" are the prediction range integers with the lower one first
    ifstream infile(filename.c_str());

    if(!infile){load_error("Die Datei "+filename+" exisitiert nicht");}

    char buff[1024]; 
    char copybuff[1024];
    char *temp;
    char *temp_inside;
    pair<int,int> pred_range;
    unordered_map<string,Transcript> transcript_hash;

    infile.getline(buff, 1024);
    list<string> unknownFeatures;
    while (infile){
	// if line dont starts with '#'
	if (buff[0]!='#'){
	    if ((strstr(buff, "gene_id")!=NULL) && (strstr(buff, "transcript_id")!=NULL)){
		Exon exon;
		Transcript transcript;
		transcript.priority = priority;

		Gene gene;
		strncpy(copybuff, buff, 1014);
		if (strstr(buff, "\t")==NULL) {
		    load_error("Line not tab separated.");
		}
		temp = strtok(buff, "\t");
		if (temp){
		    exon.chr = temp;
		}else 
		    load_error("Can not read sequence name.");
		temp = strtok(NULL, "\t");
		if (temp)
		    transcript.source = temp;
		else
		    load_error("Can not read second column.");
		temp = strtok(NULL, "\t");
		if (temp)
		    exon.feature = temp;
		else
		    load_error("Can not read feature.");
		if (exon.feature != "CDS" && exon.feature != "start_codon" && exon.feature != "stop_codon" && exon.feature != "exon" && exon.feature != "UTR"){
		    if (exon.feature != "gene" && exon.feature != "transcript" && exon.feature != "intron"){
			list<string>::iterator fit = find(unknownFeatures.begin(),unknownFeatures.end(),exon.feature);
			if (fit == unknownFeatures.end()){
			    unknownFeatures.push_back(exon.feature);
			}
		    }
		    continue;
		}
		temp = strtok(NULL, "\t");
		if (temp)
		    exon.from = atoi(temp);
		else 
		    load_error("Can not read start position.");
		temp = strtok(NULL, "\t");
		if (temp)
		    exon.to = atoi(temp);
		else
		    load_error("Can not read end position.");
		temp = strtok(NULL, "\t");
		if (temp) 
		    exon.score = atof(temp);
		else 
		    load_error("Can not read score.");
		temp = strtok(NULL, "\t");
		if (!temp)
		    load_error("Can not read strand.");
		if (strcmp(temp, "+") == 0)
		    transcript.strand = '+';
		else if (strcmp(temp, "-") == 0)
		    transcript.strand = '-';
		else {
		    transcript.strand = '.';
		    load_error("The strand of this transcript is unknown.");
		}
		temp = strtok(NULL, "\t");
		if (!temp)
		    load_error("Can not read frame.");
		if (strcmp(temp, "0") == 0)
		    exon.frame = 0;
		else if (strcmp(temp, "1") == 0)
		    exon.frame = 1;
		else if (strcmp(temp, "2") == 0)
		    exon.frame = 2;
		else 
		    exon.frame = -1;
		temp = strtok(NULL, "\t");
		if (temp){
		    char attribute[300];
		    strncpy(attribute, temp, 300);
		    string gene_id;
		    string transcript_id;
		    temp_inside = strtok(attribute, "\"");
		    while (temp_inside && (gene_id.empty() || transcript_id.empty())){
			if (strstr(temp_inside, "transcript_id")!=NULL){
			    temp_inside = strtok(NULL, "\"");
			    transcript_id = temp_inside;
			    if (transcript_hash.find(transcript_id) == transcript_hash.end()){
				transcript.t_id = transcript_id;
				transcript_hash[transcript_id] = transcript;
			    }else if (transcript.strand != transcript_hash[transcript_id].strand){
				load_warning("One transcript on different strands.");
			    }
			    if (exon.feature == "start_codon"){
				if (transcript_hash[transcript_id].tl_complete.first)
				    transcript_hash[transcript_id].separated_codon.first = true;
				if (transcript_hash[transcript_id].strand == '+'){
				    if (!transcript_hash[transcript_id].tl_complete.first || (transcript_hash[transcript_id].tl_complete.first && transcript_hash[transcript_id].tis > exon.from))
					transcript_hash[transcript_id].tis = exon.from;
				}
				if (transcript_hash[transcript_id].strand == '-'){
				    if (!transcript_hash[transcript_id].tl_complete.first || (transcript_hash[transcript_id].tl_complete.first && transcript_hash[transcript_id].tis < exon.to))
					transcript_hash[transcript_id].tis = exon.to;
				}
				transcript_hash[transcript_id].tl_complete.first = true;
			    }else if (exon.feature == "stop_codon"){
				if (transcript_hash[transcript_id].tl_complete.second){
				    transcript_hash[transcript_id].separated_codon.second = true;}
				if (transcript_hash[transcript_id].strand == '+'){
				    if (!transcript_hash[transcript_id].tl_complete.second || (transcript_hash[transcript_id].tl_complete.second && transcript_hash[transcript_id].tes < exon.to))
					transcript_hash[transcript_id].tes = exon.to;
				}
				if (transcript_hash[transcript_id].strand == '-'){
				    if (!transcript_hash[transcript_id].tl_complete.second || (transcript_hash[transcript_id].tl_complete.second && transcript_hash[transcript_id].tes > exon.from))
					transcript_hash[transcript_id].tes = exon.from;
				}
				transcript_hash[transcript_id].tl_complete.second = true;
			    }else{
				transcript_hash[transcript_id].exon_list.push_front(exon);
			    }
			    if (pred_range.first && pred_range.second){
				transcript_hash[transcript_id].pred_range = pred_range;
			    }
			}
			if (strstr(temp_inside, "gene_id")!=NULL){
			    temp_inside = strtok(NULL, "\"");
			    gene_id = temp_inside;
			    if (gene_map.count(gene_id) == 0)
				{
				    gene.g_id = gene_id;
				    gene_map[gene_id] = gene;
				}
			}
			temp_inside = strtok(NULL, "\"");
		    }
		    transcript_hash[transcript_id].parent = &gene_map[gene_id];
		}
		else 
		    load_error("Can not read last column.");
	    }
	    else{}	// cerr << "A line without gene_id and/or transcript_id." << endl;
	}else {			// if line starts with '#'
	    strncpy(copybuff, buff, 1014);
	    if (strstr(buff, "sequence range")!=NULL) {
		temp = strtok(buff, ":");
		temp = strtok(NULL, ":- ()bp");
		pred_range.first = atoi(temp);
		temp = strtok(NULL, ":- ()bp");
		pred_range.second = atoi(temp);
		if ((pred_range.second - pred_range.first) <= 0){
		    load_warning("There is a non positiv prediction range.");
		}
	    } // at this position other optional options could be read from gff-#-lines with else if
	}
	infile.getline(buff, 1024);
    }
    for(auto pointer = transcript_hash.begin(); pointer != transcript_hash.end(); pointer++)
	{
	    transcript_list.push_front((*pointer).second);
	    gene_map[(*pointer).second.parent->g_id].children.push_front(&transcript_list.front());
	}
    if (!unknownFeatures.empty()){
	for (list<string>::iterator it = unknownFeatures.begin(); it != unknownFeatures.end(); it++){
	    load_warning("There is the unexpected feature \"" + *it + "\" that is not equal to \"CDS\", \"UTR\", \"exon\", \"intron\", \"gene\", \"transcript\", \"start_codon\" and \"stop_codon\". This feature is going to be ignored.");
	}
    }
}

void saveOverlap(list<Transcript*> &overlap, string outFileName, Properties &properties)
{
    // outputs overlap at the end of an existing file in gff format
    // every first and last outfile-line is adjusted and might be change back (to comments above)
    if (overlap.size() == 0) {return;}
    fstream outfile;
    outfile.open(outFileName, ios::out | ios::app);
    outfile << "# overlap start --------------------------------------------------------------------------------" << endl;
    outfile << "# this overlap has " << overlap.size() << " different transcripts" << endl;
    // write by transcripts:
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if (find(properties.suppList.begin(),properties.suppList.end(),(*it)->priority) != properties.suppList.end()){continue;}
	outfile << "# " << (*it)->t_id << " is supported by " << (*it)->supporter.size() << " other predictions" << endl;
	outfile << "# core-transcrpit has priority " << (*it)->priority << endl;
	if ((*it)->strand == '+' && (*it)->tl_complete.first){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    //outfile << "chr" << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "start_codon" << "\t";
	    outfile << (*it)->tis << "\t";
	    outfile << ((*it)->tis+2) << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    // outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    outfile << "transcript_id \"" << "jg" << temp_gene_number << ".t1" << "\"; gene_id \"" << "jg" << temp_gene_number << "\";" << endl;
	}
	else if ((*it)->strand == '-' && (*it)->tl_complete.second){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    //outfile << "chr" << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "stop_codon" << "\t";
	    outfile << (*it)->tes << "\t";
	    outfile << ((*it)->tes+2) << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    // outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    outfile << "transcript_id \"" << "jg" << temp_gene_number << ".t1" << "\"; gene_id \"" << "jg" << temp_gene_number << "\";" << endl;
	}
	for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
	    //cout << (*it_inside).feature << " " <<  (*it_inside).frame << " " <<  (*it_inside).to << " " << (*it_inside).from << " " << temp_gene_number << endl;
	    outfile << (*it)->exon_list.front().chr << "\t";
	    //outfile << "chr" << (*it_inside).chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << (*it_inside).feature << "\t";
	    outfile << (*it_inside).from << "\t";
	    outfile << (*it_inside).to << "\t";
	    outfile << (*it_inside).score << "\t";
	    outfile << (*it)->strand << "\t";
	    if ((*it_inside).frame != -1)
		outfile << (*it_inside).frame << "\t";
	    else
		outfile << "." << "\t";
	    // outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    outfile << "transcript_id \"" << "jg" << temp_gene_number << ".t1" << "\"; gene_id \"" << "jg" << temp_gene_number << "\";" << endl;
	}
	if ((*it)->strand == '-' && (*it)->tl_complete.first){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    //outfile << "chr" << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "start_codon" << "\t";
	    outfile << ((*it)->tis-2) << "\t";
	    outfile << (*it)->tis << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    // outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    outfile << "transcript_id \"" << "jg" << temp_gene_number << ".t1" << "\"; gene_id \"" << "jg" << temp_gene_number << "\";" << endl;
	}
	else if ((*it)->strand == '+' && (*it)->tl_complete.second){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    //outfile << "chr" << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "stop_codon" << "\t";
	    outfile << ((*it)->tes-2) << "\t";
	    outfile << (*it)->tes << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    // outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    outfile << "transcript_id \"" << "jg" << temp_gene_number << ".t1" << "\"; gene_id \"" << "jg" << temp_gene_number << "\";" << endl;
	}
	temp_gene_number++;
    }
    outfile.close();
}

void save(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, Properties &properties)		// pointer to gene_map is only for "write by genes" necessary
{
    // this is and old save version, which saves by transcript list or by gene hash map but only saves the exon features (no start/stop codon)
    fstream outfile;
    outfile.open(properties.outFileName, ios::out);		// delete what is allready in the file (ios::out | ios::app - dont delete what is allready in the file)
	
    // write by transcripts:
    for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
	for (list<Exon>::iterator it_inside = (*it).exon_list.begin(); it_inside != (*it).exon_list.end(); it_inside++){
	    outfile << (*it_inside).chr << "\t";
	    outfile << (*it).source << "\t";
	    outfile << (*it_inside).feature << "\t";
	    outfile << (*it_inside).from << "\t";
	    outfile << (*it_inside).to << "\t";
	    outfile << (*it_inside).score << "\t";
	    outfile << (*it).strand << "\t";
	    if ((*it_inside).frame != -1)
		outfile << (*it_inside).frame << "\t";
	    else
		outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it).t_id << "\"; gene_id \"" << (*it).parent->g_id << "\";" << endl;
	}
    }
    // write by genes (unordered):
    /*
      for (auto pointer = gene_map.begin(); pointer != gene_map.end(); pointer++){
      for (list<Transcript*>::iterator it = pointer->second.children.begin(); it != pointer->second.children.end(); it++){
      for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
      outfile << (*it_inside).chr << "\t";
      outfile << (*it)->source << "\t";
      outfile << (*it_inside).feature << "\t";
      outfile << (*it_inside).from << "\t";
      outfile << (*it_inside).to << "\t";
      outfile << (*it_inside).score << "\t";
      outfile << (*it).strand << "\t";
      if ((*it_inside).frame != -1)
      outfile << (*it_inside).frame << "\t";
      else
      outfile << "." << "\t";
      outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << pointer->second.g_id << "\";" << endl;
      }
      }
      }
    */
    outfile.close();
}
