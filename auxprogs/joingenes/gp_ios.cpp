#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

void load(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename, int &priority)
{
	ifstream infile(filename.c_str());

	char buff[1024]; 
	char copybuff[1024];
	char *temp;
	char *temp_inside;
	
	//unsigned int k = 0;
	//unsigned int grenze = 100000;

	infile.getline(buff, 1024);
	while (infile){	// && k < grenze
	//k++;
	// if NOT line starts with '#'	
		if (buff[0]!='#'){
			if ((strstr(buff, "gene_id")!=NULL) && (strstr(buff, "transcript_id")!=NULL)){
				Exon exon;
				Transcript transcript;
				transcript.priority = priority;
				Gene gene;
				strncpy(copybuff, buff, 1014);
				if (strstr(buff, "\t")==NULL) {
					cerr << "Line not tab separated." << endl;
				}
				temp = strtok(buff, "\t");
				if (temp)
					exon.chr = temp;
				else 
					cerr << "Could not read sequence name." << endl;
				temp = strtok(NULL, "\t");
				if (temp)
					transcript.source = temp;
				else
					cerr << "Could not read second column." << endl;
				temp = strtok(NULL, "\t");
				if (temp)
					exon.feature = temp;
				else
					cerr << "Could not read frame." << endl;
				temp = strtok(NULL, "\t");
				if (temp)
					exon.from = atoi(temp);
				else 
					cerr << "Could not read start position." << endl;
				temp = strtok(NULL, "\t");
				if (temp)
					exon.to = atoi(temp);
				else
					cerr << "Could not read end position." << endl;
				temp = strtok(NULL, "\t");
				if (temp) 
					exon.score = atof(temp);
				else 
					cerr << "Could not read score." << endl;
				temp = strtok(NULL, "\t");
				if (!temp)
					cerr << "Could not read strand." << endl;
				if (strcmp(temp, "+") == 0)
					transcript.strand = '+';
				else if (strcmp(temp, "-") == 0)
					transcript.strand = '-';
				else 
					transcript.strand = '.';
				temp = strtok(NULL, "\t"); // frame
				if (!temp)
					cerr << "Could not read frame." << endl;
				if (strcmp(temp, "0") == 0)
					exon.frame = 0;
				else if (strcmp(temp, "1") == 0)
					exon.frame = 1;
				else if (strcmp(temp, "2") == 0)
					exon.frame = 2;
				else 
					exon.frame = -1;
				temp = strtok(NULL, "\t");
				// if file is gff:
				if (temp){
					char attribute[300];
					strncpy(attribute, temp, 300);
					string gene_id;
					string transcript_id;
					temp_inside = strtok(attribute, "\"");
					while (temp_inside && (gene_id.empty() || transcript_id.empty())){
						if (strstr(temp_inside, "transcript_id")!=NULL){
							temp_inside = strtok(NULL, "\"");
							//cout << "transcript_id: " << temp_inside << endl;
							transcript_id = temp_inside;
							if (transcript_list.empty()){
								transcript.t_id = transcript_id;
								transcript_list.push_front(transcript);
							}else if (transcript_list.front().t_id != transcript_id){		// it is necessary, that features from the same transcript are in a row
								transcript.t_id = transcript_id;
								transcript_list.push_front(transcript);
							}else if (transcript.strand != transcript_list.front().strand){
								cout << "care, one transcript on different strands" << endl;
							}
							if (exon.feature == "start_codon"){
								transcript_list.front().start_codon = exon.from;
							}else if (exon.feature == "stop_codon"){
								transcript_list.front().stop_codon = exon.from;
							}else{
								transcript_list.front().exon_list.push_front(exon);
							}
							//transcript_list.front().frame = (exon.from + exon.frame) % 3;
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
					transcript_list.front().parent = &gene_map[gene_id];

					if (gene_map[gene_id].children.empty()){
						gene_map[gene_id].children.push_front(&transcript_list.front());
					}else if (gene_map[gene_id].children.front()->t_id != transcript_id){
						gene_map[gene_id].children.push_front(&transcript_list.front());
					}
				}
				else 
					cerr << "Could not read last column." << endl;

			}
			else{}	//cerr << "A line without gene_id and/or transcript_id." << endl;
		}
		infile.getline(buff, 1024);
	}
	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
		(*it).start = (*it).getstart();
		(*it).stop = (*it).getstop();
	} 
}

void save_overlap(list<Transcript*> &overlap, string &filename)		// pointer to gene_map is only for "write by genes" necessary
{
	fstream outfile;
	outfile.open(filename, ios::out | ios::app);
	outfile << "# overlap start --------------------------------------------------------------------------------" << endl;
	outfile << "# this overlap has " << overlap.size() << " different transcripts" << endl;
// write by transcripts:
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		outfile << "# " << (*it)->t_id << " is supported by " << (*it)->supporter.size() << " other predictions" << endl;
		outfile << "# it has ";
		if ((*it)->start_codon)
			outfile << "a start_codon";
		else
			outfile << "no start_codon";
		if ((*it)->stop_codon)
			outfile << " and a stop_codon";
		else
			outfile << " and no stop_codon";
		outfile << endl;
		for (list<Exon>::iterator it_inside = (*it)->exon_list.begin(); it_inside != (*it)->exon_list.end(); it_inside++){
			outfile << (*it_inside).chr << "\t";
			outfile << (*it)->source << "\t";
			outfile << (*it_inside).feature << "\t";
			outfile << (*it_inside).from << "\t";
			outfile << (*it_inside).to << "\t";
			outfile << (*it_inside).score << "\t";
			outfile << (*it)->strand << "\t";
			outfile << (*it_inside).frame << "\t";
			outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << (*it)->parent->g_id << "\";" << endl;
		}
	}
	outfile.close();
}

void save(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename)		// pointer to gene_map is only for "write by genes" necessary
{
	fstream outfile;
	outfile.open(filename, ios::out);		// delete what is allready in the file (ios::out | ios::app - dont delete what is allready in the file)
	
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
			outfile << (*it_inside).frame << "\t";
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
				outfile << (*it_inside).frame << "\t";
				outfile << "transcript_id \"" << (*it)->t_id << "\"; gene_id \"" << pointer->second.g_id << "\";" << endl;
			}
		}
	}
*/
	outfile.close();
}
