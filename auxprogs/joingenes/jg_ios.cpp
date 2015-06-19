#include "jg_transcript.h"
#include "jg_ios.h"

using namespace std;

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

void load(unordered_map<string,Gene*> &geneMap, string &filename, int &priority, unordered_map<string,bool> &taxaMap, Properties &properties)
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
    unordered_map<string,Transcript*> thisFileTranscriptMap;
    unordered_map<string,Gene*> thisFileGeneMap;

    unordered_map<string,bool> taxaMapTemp;
    infile.getline(buff, 1024);
    list<string> unknownFeatures;
    while (infile){
	// if line dont starts with '#'
	if (buff[0]!='#'){
	    //if ((strstr(buff, "gene_id")!=NULL) && (strstr(buff, "transcript_id")!=NULL)){
		Exon exon;
		Transcript* transcript = new Transcript;
		transcript->inputFile = filename;
		(*transcript).priority = priority;
		strncpy(copybuff, buff, 1014);
		if (strstr(buff, "\t")==NULL) {
		    cerr << "ErrorLine: " << buff << endl;
		    load_error("Line not tab separated.");
		}
		temp = strtok(buff, "\t");
		if (temp){
		    exon.chr = temp;
		}else 
		    load_error("Can not read sequence name.");
		temp = strtok(NULL, "\t");
		if (temp)
		    (*transcript).source = temp;
		else
		    load_error("Can not read second column.");
		temp = strtok(NULL, "\t");
		if (temp)
		    exon.feature = temp;
		else
		    load_error("Can not read feature.");
		if (exon.feature != "CDS" && exon.feature != "start_codon" && exon.feature != "stop_codon" && exon.feature != "exon" && exon.feature != "UTR" && exon.feature != "3'-UTR" && exon.feature != "5'-UTR" && exon.feature != "tss" && exon.feature != "tts"  && exon.feature != "intron"){
		    if (exon.feature != "gene" && exon.feature != "transcript"){
			list<string>::iterator fit = find(unknownFeatures.begin(),unknownFeatures.end(),exon.feature);
			if (fit == unknownFeatures.end()){
			    unknownFeatures.push_back(exon.feature);
			}
		    }
		    infile.getline(buff, 1024);
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
		    (*transcript).strand = '+';
		else if (strcmp(temp, "-") == 0)
		    (*transcript).strand = '-';
		else {
		    (*transcript).strand = '.';
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
		    Gene* gene = new Gene;

		    // get the gene_id and transcript_id of the gtf format... If there is non, create a new one.
		    while (/*temp_inside && */(gene_id.empty() || transcript_id.empty())){
			if (((strstr(temp, "gene_id")==NULL) && (strstr(temp, "transcript_id")==NULL)) || !temp_inside){
			    gene_id = nextFreeGeneID(properties, &thisFileGeneMap);
			    (*gene).g_id = gene_id;
			    transcript_id = nextFreeTxID(gene, properties, &thisFileTranscriptMap);
			    transcript->originalId = "none";
			}
			if (strstr(temp_inside, "transcript_id")!=NULL){
			    temp_inside = strtok(NULL, "\"");
			    if (temp_inside){
				transcript_id = temp_inside;
				transcript->originalId = temp_inside;
			    }else{
				load_error("Missing id behind the flag transcript_id.");
			    }
			}
			if (strstr(temp_inside, "gene_id")!=NULL){
			    temp_inside = strtok(NULL, "\"");
			    if (temp_inside){
				gene_id = temp_inside;
			    }else{
				load_error("Missing id behind the flag gene_id in file "+filename+" at "+exon.chr+" from "+to_string((long long int) exon.from)+" to "+to_string((long long int) exon.to)+".");
			    }
			}
			temp_inside = strtok(NULL, "\"");
		    }

		    if (thisFileTranscriptMap.find(transcript_id) == thisFileTranscriptMap.end()){
			(*transcript).t_id = transcript_id;
			thisFileTranscriptMap[transcript_id] = transcript;
		    }else if ((*transcript).strand != thisFileTranscriptMap[transcript_id]->strand){
			load_warning("One transcript on different strands.");
		    }




		    if (pred_range.first && pred_range.second){
			exon.predRange = pred_range;
		    }


		    if (exon.feature == "start_codon"){
			if (thisFileTranscriptMap[transcript_id]->tl_complete.first)
			    thisFileTranscriptMap[transcript_id]->separated_codon.first = true;
			if (thisFileTranscriptMap[transcript_id]->strand == '+'){
			    if (!thisFileTranscriptMap[transcript_id]->tl_complete.first || (thisFileTranscriptMap[transcript_id]->tl_complete.first && thisFileTranscriptMap[transcript_id]->tis > exon.from))
				thisFileTranscriptMap[transcript_id]->tis = exon.from;
			}
			if (thisFileTranscriptMap[transcript_id]->strand == '-'){
			    if (!thisFileTranscriptMap[transcript_id]->tl_complete.first || (thisFileTranscriptMap[transcript_id]->tl_complete.first && thisFileTranscriptMap[transcript_id]->tis < exon.to))
				thisFileTranscriptMap[transcript_id]->tis = exon.to;
			}
			thisFileTranscriptMap[transcript_id]->tl_complete.first = true;
		    }else if (exon.feature == "stop_codon"){
			if (thisFileTranscriptMap[transcript_id]->tl_complete.second){
			    thisFileTranscriptMap[transcript_id]->separated_codon.second = true;}
			if (thisFileTranscriptMap[transcript_id]->strand == '+'){
			    if (!thisFileTranscriptMap[transcript_id]->tl_complete.second || (thisFileTranscriptMap[transcript_id]->tl_complete.second && thisFileTranscriptMap[transcript_id]->tes < exon.to)){
				//thisFileTranscriptMap[transcript_id]->tes = exon.to;
			    }
			}
			if (thisFileTranscriptMap[transcript_id]->strand == '-'){
			    if (!thisFileTranscriptMap[transcript_id]->tl_complete.second || (thisFileTranscriptMap[transcript_id]->tl_complete.second && thisFileTranscriptMap[transcript_id]->tes > exon.from)){
				//thisFileTranscriptMap[transcript_id]->tes = exon.from;
			    }
			}
			thisFileTranscriptMap[transcript_id]->tl_complete.second = true;
			thisFileTranscriptMap[transcript_id]->stop_list.push_back(exon);		// TESTLARS
		    }else if (exon.feature == "tss"){
			if (exon.from != exon.to){load_error("\"tss\" can only take place on one position.");}
			thisFileTranscriptMap[transcript_id]->tss = exon.from;
			thisFileTranscriptMap[transcript_id]->tx_complete.first = true;
		    }else if (exon.feature == "tts"){
			if (exon.from != exon.to){load_error("\"tts\" can only take place on one position.");}
			thisFileTranscriptMap[transcript_id]->tts = exon.to;
			thisFileTranscriptMap[transcript_id]->tx_complete.second = true;
		    }else if (exon.feature == "intron"){
			thisFileTranscriptMap[transcript_id]->intron_list.push_front(exon);
		    }else{
			thisFileTranscriptMap[transcript_id]->exon_list.push_front(exon);
		    }
		    if (pred_range.first && pred_range.second){
			thisFileTranscriptMap[transcript_id]->pred_range = pred_range;
		    }

		    (*gene).g_id = gene_id;
		    if (thisFileGeneMap.count(gene_id) == 0){
			thisFileGeneMap[gene_id] = gene;
		    }

		    thisFileTranscriptMap[transcript_id]->parent = thisFileGeneMap[gene_id];
		}else 
		    load_error("Can not read last column.");
	    //}else{}	// cerr << "A line without gene_id and/or transcript_id." << endl;
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
    for(auto pointer = thisFileGeneMap.begin(); pointer != thisFileGeneMap.end(); pointer++)
    {
	if ((*properties.geneMap).find( pointer->second->g_id ) != (*properties.geneMap).end()){
	    pointer->second->g_id = nextFreeGeneID(properties, &thisFileGeneMap);
	    (*properties.geneMap)[pointer->second->g_id] = (*pointer).second;
	    // maybe fill thisFileGeneMap with new gene, but is not necessary
	}else{
	    (*properties.geneMap)[pointer->second->g_id] = (*pointer).second;
	}
    }
    for(auto pointer = thisFileTranscriptMap.begin(); pointer != thisFileTranscriptMap.end(); pointer++)
    {
        /*if (thisFileGeneMap.find(pointer->second->parent->g_id) == thisFileGeneMap.end()){
		load_error("Something went wrong in the load function, because one transcript is part of an unknown gene.");
	}*/
	(*properties.geneMap)[pointer->second->parent->g_id]->children.push_front(pointer->second);
        if ((*properties.transcriptMap).find(pointer->second->t_id) != (*properties.transcriptMap).end()){
	    pointer->second->t_id = nextFreeTxID(pointer->second->parent, properties, &thisFileTranscriptMap);
	    (*properties.transcriptMap)[pointer->second->t_id] = pointer->second;
	}else{
	    (*properties.transcriptMap)[pointer->second->t_id] = pointer->second;
	}
    }

    if (!unknownFeatures.empty()){
	for (list<string>::iterator it = unknownFeatures.begin(); it != unknownFeatures.end(); it++){
	    load_warning("Did not expect feature \"" + *it + "\". Known features are \"CDS\", \"UTR\", \"3'-UTR\", \"5'-UTR\", \"exon\", \"intron\", \"gene\", \"transcript\", \"tss\", \"tts\", \"start_codon\" and \"stop_codon\". This feature is going to be ignored.");
	}
    }
    taxaMap.insert(taxaMapTemp.begin(),taxaMapTemp.end());
}

/*bool taxonInList(list<Transcript*> List, string Taxon){
    for(list<Transcript*>::iterator it = List.begin(); it != List.end(); it++){
	if ((*it)->t_id == Taxon){
	    return true;
	}
    }
    return false;
}

void mergeTranscripts(Properties &properties, Transcript* acc, Transcript* don){
    if (acc->parent->g_id != don->parent->g_id){
	load_error("Only Transcripts from same gene can be merged.");
    }
    if (acc->t_id != don->t_id){
	load_error("Something went wrong in the mergeTranscript function, because only transcripts with same name should be merged.");
    }
    if (acc->strand != don->strand){
	load_error("Transcripts on different strands are not mergeable.");
    }
    if (acc->priority != don->priority){
	load_error("Transcripts with different priorities are not mergeable.");
    }
    if (){}


//    int tis;
//    int tes;
//    pair<bool,bool> tx_complete;
//    pair<bool,bool> tl_complete;

// merge exon list:
//    list<Exon> exon_list;
}

void mergeGenes(Properties &properties, Gene* acc, Gene* don){
    if (acc->g_id != don->g_id){
	load_error("Something went wrong in the mergeGenes function, because only genes with same name should be joined.");
    }
    if (acc->nrOfTx < don->nrOfTx){
	acc->nrOfTx = don->nrOfTx;
    }
    for(list<Transcript*>::iterator it = don->children.begin(); it != don->children.end(); it++){
	if (taxonInList(acc->children, (*it)->t_id)){
	    mergeTranscripts();f
	}else{
	    if ((*properties.transcriptMap).find((*it)->t_id) == (*properties.transcriptMap).end()){
		acc->children.push_front(*it);
		
	    }else{
		load_error("Something went wrong in the load function, because one transcript is part of two different gene. (maybe enable option newTaxa)");
	    }

	}
    }
}*/

void renameTaxa(list<Transcript*> &overlap, Properties &properties){
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->parent->nrOfPrintedTx == 0){
	    properties.nrOfPrintedGenes++;
	    (*it)->parent->g_id = "jg" + to_string((long long int) properties.nrOfPrintedGenes);
	}
	(*it)->parent->nrOfPrintedTx++;
	(*it)->t_id = (*it)->parent->g_id + ".t" + to_string((long long int) (*it)->parent->nrOfPrintedTx);
    }
}

void saveOverlap(list<Transcript*> &overlap, string outFileName, Properties &properties)
{
    renameTaxa(overlap, properties);
    // outputs overlap at the end of an existing file in gff format
    // every first and last outfile-line is adjusted and might be change back (to comments above)
    if (overlap.size() == 0) {return;}
    fstream outfile;
    outfile.open(outFileName, ios::out | ios::app);
    outfile << "# overlap start --------------------------------------------------------------------------------" << endl;
    outfile << "# this overlap has " << overlap.size() << " different transcripts" << endl;
    // write by transcripts:
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){

	if (find(properties.supprList.begin(),properties.supprList.end(),(*it)->priority) != properties.supprList.end()){
	    outfile << "# " << (*it)->t_id << " is suppressed." << endl;
	    continue;
	}
	outfile << "# This transcript "<< (*it)->t_id << " is derived from " << (*it)->originalId << " from the input file " << (*it)->inputFile << endl;
	outfile << "# It is supported by " << (*it)->supporter.size() << " other predicted genes" << endl;
	outfile << "# the core of this joined transcript has priority " << (*it)->priority << endl;
	if (!(*it)->joinpartner.first.empty())
	    outfile << "# transcrpit has been joined at 5'-side with " << (*it)->joinpartner.first << endl;
	if (!(*it)->joinpartner.second.empty())
	    outfile << "# transcrpit has been joined at 3'-side with " << (*it)->joinpartner.second << endl;

	/*if (!(*it)->parent->printed){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "gene" << "\t";
	    int minStart = (*it)->getTxStart();
	    int maxStop = (*it)->getTxEnd();
	    for (list<Transcript*>::iterator it_inside = (*it)->parent->children.begin(); it_inside != (*it)->parent->children.end(); it_inside++){
		if (minStart > (*it_inside)->getTxStart()){minStart=(*it_inside)->getTxStart();}
		if (maxStop < (*it_inside)->getTxEnd()){maxStop=(*it_inside)->getTxEnd();}
	    }
	    outfile << minStart << "\t";
	    outfile << maxStop << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << (*it)->parent->g_id << endl;
	    (*it)->parent->printed = true;
	}*/

	if (!(*it)->exon_list.empty()){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "transcript" << "\t";
	    outfile << (*it)->getTxStart() << "\t";
	    outfile << (*it)->getTxEnd() << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << (*it)->t_id << endl;
	}

	if ((*it)->strand == '+' && (*it)->tss!=-1){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "tss" << "\t";
	    outfile << (*it)->tss << "\t";
	    outfile << (*it)->tss << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	if ((*it)->strand == '-' && (*it)->tts!=-1){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "tts" << "\t";
	    outfile << (*it)->tts << "\t";
	    outfile << (*it)->tts << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
	    if ((*it_inside).feature == "CDS"){break;}
	    outfile << (*it_inside).chr << "\t";
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
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	if ((*it)->strand == '+' && (*it)->tl_complete.first){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "start_codon" << "\t";
	    outfile << (*it)->tis << "\t";
	    outfile << ((*it)->tis+2) << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;

	}
	else if ((*it)->strand == '-' && (*it)->tl_complete.second){
	/*    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "stop_codon" << "\t";
	    outfile << (*it)->tes << "\t";
	    outfile << ((*it)->tes+2) << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
*/

	  for (list<Exon>::iterator it_inside = (*it)->stop_list.begin(); it_inside != (*it)->stop_list.end(); it_inside++){
	    if ((*it_inside).feature != "stop_codon"){continue;}
	    outfile << (*it_inside).chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << (*it_inside).feature << "\t";
	    outfile << (*it_inside).from << "\t";
	    outfile << (*it_inside).to << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    if ((*it_inside).frame != -1)
		outfile << (*it_inside).frame << "\t";
	    else
		outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	  }
	}

	if ((*it)->intron_list.size() > 0 && !properties.join){
	    (*it)->intron_list.sort();
	    for (list<Exon>::iterator it_inside = (*it)->intron_list.begin(); it_inside != (*it)->intron_list.end(); it_inside++){
		outfile << (*it_inside).chr << "\t";
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
		outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	    }
	}

	for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
	    if ((*it_inside).feature != "CDS"){continue;}
	    outfile << (*it_inside).chr << "\t";
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
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;

	}
	if ((*it)->strand == '-' && (*it)->tl_complete.first){
	    outfile << (*it)->exon_list.back().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "start_codon" << "\t";
	    outfile << ((*it)->tis-2) << "\t";
	    outfile << (*it)->tis << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	else if ((*it)->strand == '+' && (*it)->tl_complete.second){
/*	    outfile << (*it)->exon_list.back().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "stop_codon" << "\t";
	    outfile << ((*it)->tes-2) << "\t";
	    outfile << (*it)->tes << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << '0' << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
*/
	  for (list<Exon>::iterator it_inside = (*it)->stop_list.begin(); it_inside != (*it)->stop_list.end(); it_inside++){
	    if ((*it_inside).feature != "stop_codon"){continue;}
	    outfile << (*it_inside).chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << (*it_inside).feature << "\t";
	    outfile << (*it_inside).from << "\t";
	    outfile << (*it_inside).to << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    if ((*it_inside).frame != -1)
		outfile << (*it_inside).frame << "\t";
	    else
		outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	  }
	}
	bool backUTR = false;
	for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
	    if ((*it_inside).feature != "CDS" && !backUTR){continue;}
	    if ((*it_inside).feature == "CDS"){backUTR = true; continue;}
	    outfile << (*it_inside).chr << "\t";
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
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	if ((*it)->strand == '+' && (*it)->tts!=-1){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "tts" << "\t";
	    outfile << (*it)->tts << "\t";
	    outfile << (*it)->tts << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
	if ((*it)->strand == '-' && (*it)->tss!=-1){
	    outfile << (*it)->exon_list.front().chr << "\t";
	    outfile << (*it)->source << "\t";
	    outfile << "tss" << "\t";
	    outfile << (*it)->tss << "\t";
	    outfile << (*it)->tss << "\t";
	    outfile << '.' << "\t";
	    outfile << (*it)->strand << "\t";
	    outfile << "." << "\t";
	    outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
	}
    }
    outfile.close();
}

string nextFreeGeneID(Properties &properties, unordered_map<string,Gene*>* addGeneMap){
    string geneID;
    while (geneID == "" || (*properties.geneMap).find(geneID) != (*properties.geneMap).end() || (addGeneMap != NULL && (*addGeneMap).find(geneID) != (*addGeneMap).end())){
	geneID = "g" + to_string((long long int) properties.unknownCount);
	properties.unknownCount++;
    }
    return geneID;
}

string nextFreeTxID(Gene* gene, Properties &properties, unordered_map<string,Transcript*>* addTranscriptMap){
    string txID;
    while (txID == "" || (*properties.transcriptMap).find(txID) != (*properties.transcriptMap).end() || (addTranscriptMap != NULL && (*addTranscriptMap).find(txID) != (*addTranscriptMap).end())){
	txID = gene->g_id + ".t" + to_string((long long int) gene->nrOfTx);
	gene->nrOfTx++;
    }
    return txID;
}

