#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

int overallnumber=0;

void seek_overlaps(list<Transcript> &transcript_list){			// needs a sorted list of transcripts
	list<Transcript*> overlap;
	int max_base = transcript_list.front().stop;
	
	string filename = "overlap_out";
	fstream outfile;
	outfile.open(filename, ios::out);		// delete content of filename
	outfile.close();

	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
		if ((*it).start < max_base){
			overlap.push_back(&*it);
			if (max_base < (*it).stop){
				max_base = (*it).stop;
			}
		}
		else{
			if (overlap.size() > 1)
			{
				//cout << overlap.size() << " changed to ";
				decide(overlap);
				//cout << overlap.size() << endl;
				//save_overlap(overlap, filename);
			}
			overlap.clear();
			max_base = (*it).stop;
			overlap.push_front(&*it);
		}
	}
}

void decide(list<Transcript*> &overlap)
{
	//cout << overlap.size() << " ";
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		list<Transcript*>::iterator it_temp = it;
		it_temp++;
		for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
			if (compare_transcripts(*it, *it_inside)){
				//cout << (*it)->t_id << " and " << (*it_inside)->t_id << " are equal" << endl;
				(*it)->supporter.push_front(*it_inside);
				(*it_inside)->supporter.push_front(*it);
				if ((*it)->priority < (*it_inside)->priority){
					it = overlap.erase(it);
					it_inside = it_temp;
				}else{
					it_inside = overlap.erase(it_inside);
					it_inside--;
				}
			}
		}
	}

	/*cout << "overlap: ";
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		cout << (*it)->t_id << " ";
	}
	cout << endl;*/

	//cout << overlap.size() << " ";
	// test für is_part_of():
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		list<Transcript*>::iterator it_temp = it;
		it_temp++;
		for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
			pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
			if (who_is_part.first == true){
				if (who_is_part.second == true)
					cout << "SOMETHING WENT WRONG! Because this shouldnt happen after erasing equal transcripts" << endl;	// this case could be taken instead of "compare_transcripts"
				else{
					it = overlap.erase(it);
					it_inside = it;
					//cout << "t1 ist part of t2" << endl;
					/*cout << "t2--------------- " << (*it_inside)->supporter.size() << " supporter " << (*it_inside)->t_id << endl;
					for (list<Exon>::iterator it_iii = (*it_inside)->exon_list.begin(); it_iii != (*it_inside)->exon_list.end(); it_iii++){	
						cout << (*it_iii).from << "\t";
						cout << (*it_iii).to << "\t";
						cout << endl;
					}
					cout << "t1--------------- " << (*it)->supporter.size() << " supporter " << (*it)->t_id << endl;
					for (list<Exon>::iterator it_iii = (*it)->exon_list.begin(); it_iii != (*it)->exon_list.end(); it_iii++){	
						cout << (*it_iii).from << "\t";
						cout << (*it_iii).to << "\t";
						cout << endl;
					}
					cout << "-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-" << endl;*/
				}
			}else{
				if (who_is_part.second == true){
					it_inside = overlap.erase(it_inside);
					it_inside--;
					//cout << "t2 ist part of t1" << endl;
					/*cout << "t1--------------- " << (*it)->supporter.size() << " supporter " << (*it)->t_id << endl;
					for (list<Exon>::iterator it_iii = (*it)->exon_list.begin(); it_iii != (*it)->exon_list.end(); it_iii++){	
						cout << (*it_iii).from << "\t";
						cout << (*it_iii).to << "\t";
						cout << endl;
					}
					cout << "t2--------------- " << (*it_inside)->supporter.size() << " supporter " << (*it_inside)->t_id << endl;
					for (list<Exon>::iterator it_iii = (*it_inside)->exon_list.begin(); it_iii != (*it_inside)->exon_list.end(); it_iii++){	
						cout << (*it_iii).from << "\t";
						cout << (*it_iii).to << "\t";
						cout << endl;
					}
					cout << "-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-" << endl;*/
				}
			}
		}
	}

	//cout << overlap.size() << " ";
	//cout << "---------------------------------------------------------" << endl;
	//join:
	overlap.sort(compare_priority);
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		//cout << (*it)->start << " " << (*it)->stop << " " << (*it)->priority << endl;
		list<Transcript*>::iterator it_temp = it;
		it_temp++;
		for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
			join(*it, *it_inside);
		}
	}
}

bool compare_transcripts(Transcript* t1, Transcript* t2)
{
	bool equal = false;
	if (t1->strand == t2->strand && t1->frame == t2->frame && t1->start == t2->start && t1->stop == t2->stop && t1->exon_list.size() == t2->exon_list.size()){
		equal = true;
		list<Exon>::iterator it2 = t2->exon_list.begin();
		for (list<Exon>::iterator it1 = t1->exon_list.begin(); it1 != t1->exon_list.end(); it1++){
			if ((*it1).from == (*it2).from && (*it1).to == (*it2).to){
			}else{
				equal = false;
				break;
			}
			it2++;
		}
	}
	return equal;
}

pair<bool,bool> is_part_of(Transcript* t1, Transcript* t2)
{
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello" << endl;}
	bool t1_is_part = false;
	bool t2_is_part = false;
	if (t1->strand == t2->strand && t1->frame == t2->frame){
		t1_is_part = true;
		t2_is_part = true;
		if (t1->start_codon && t2->start_codon && t1->start_codon != t2->start_codon){
			t1_is_part = false;
			t2_is_part = false;
		}
		if (t1->stop_codon && t2->stop_codon && t1->stop_codon != t2->stop_codon){
			t1_is_part = false;
			t2_is_part = false;
		}
		list<Exon>::iterator it1 = t1->exon_list.begin();
		list<Exon>::iterator it2 = t2->exon_list.begin();
		while (t1_is_part == true || t2_is_part == true){
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello1" << endl;}
			if ((*it1).from == (*it2).from){
				if((*it1).to == (*it2).to){
					it1++;
					it2++;
				}
				else{
					t1_is_part = false;
					t2_is_part = false;
					break;
				}
			}else if ((*it1).from > (*it2).from){
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello2" << endl;}
				t2_is_part = false;
				it2++;
			}else{
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello3" << endl;}
				t1_is_part = false;
				it1++;
			}
			if (it1 == t1->exon_list.end() && it2 == t2->exon_list.end()){
				break;
			}else{
				if (it1 == t1->exon_list.end() && !(it2 == t2->exon_list.end()))
				{
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello4" << endl;}
					t2_is_part = false;
					break;
				}
				if (it2 == t2->exon_list.end() && !(it1 == t1->exon_list.end()))
				{
// if ((t1->t_id == "g7213.t1" || t1->t_id == "g7215.t1") && (t2->t_id == "g7213.t1" || t2->t_id == "g7215.t1")){cout << "hello5" << endl;}
					t1_is_part = false;
					break;
				}
			}
		}
	}
	return make_pair(t1_is_part, t2_is_part);
}

bool compare_priority(Transcript* lhs, Transcript* rhs){
	return ( lhs->priority > rhs->priority );
}

void join(Transcript* t1, Transcript* t2){	// (in Arbeit)
		if (three_compatible(t1,t2)){
			/*overallnumber++;
			cout << overallnumber << endl;*/
			/*if (overallnumber < 4){
				cout << "-----------------------------" << endl;
				cout << t1->t_id << ", Strand: " << t1->strand << ", Priority: " << t1->priority << endl;
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					cout << (*it).from << " " << (*it).to << endl;
				}
				cout << t2->t_id << ", Strand: " << t2->strand << ", Priority: " << t2->priority << endl;
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					cout << (*it).from << " " << (*it).to << endl;
				}
				//pair<bool,bool> ttt = is_part_of(t1,t2);
				//cout << ttt.first << " " << ttt.second << endl;
			}*/
			
			/*if (t1->strand == '+'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
				}
			}else if (t1->strand == '-'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					while (t1->exon_list.front().from != (*it).from){
						// welches ist welches??? muss vorher abgefragt werden
					}
				}
			}else{}		// no '+' or '-' strand*/
		}
		if (five_compatible(t1,t2)){
			overallnumber++;
			//cout << overallnumber << endl;
				if (overallnumber < 4){
				cout << "-----------------------------" << endl;
				cout << t1->t_id << ", Strand: " << t1->strand << ", Priority: " << t1->priority << endl;
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					cout << (*it).from << " " << (*it).to << endl;
				}
				cout << t2->t_id << ", Strand: " << t2->strand << ", Priority: " << t2->priority << endl;
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					cout << (*it).from << " " << (*it).to << endl;
				}
				//pair<bool,bool> ttt = is_part_of(t1,t2);
				//cout << ttt.first << " " << ttt.second << endl;
			}
		}
}

bool three_compatible(Transcript* t1, Transcript* t2){		// stop_codon-Seite	// Translation 5'->3'
	if (t1->strand == t2->strand && t1->frame == t2->frame){
		if (!(t1->stop_codon) && t2->stop_codon)
		{
			if (t1->strand == '+'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					if (t1->exon_list.back().to == (*it).to){
						return true;
					}
				}
			}else if (t1->strand == '-'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					if (t1->exon_list.front().from == (*it).from){
						return true;
					}
				}
			}else{}		// no '+' or '-' strand
		}else if(t1->stop_codon && !(t2->stop_codon)){
			if (t1->strand == '+'){
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					if (t2->exon_list.back().to == (*it).to){
						return true;
					}
				}
			}else if (t1->strand == '-'){
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					if (t2->exon_list.front().from == (*it).from){
						return true;
					}
				}
			}else{}		// no '+' or '-' strand*/ return false;
		}else{
			return false;
		}
	}else{
		return false;
	}
	return false;
}


bool five_compatible(Transcript* t1, Transcript* t2){	// (in Arbeit)		// start_codon-Seite
	if (t1->strand == t2->strand && t1->frame == t2->frame){
		if (!(t1->start_codon) && t2->start_codon)
		{
			if (t1->strand == '-'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					if (t1->exon_list.back().to == (*it).to){
						return true;
					}
				}
			}else if (t1->strand == '+'){
				for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
					if (t1->exon_list.front().from == (*it).from){
						return true;
					}
				}
			}else{}		// no '+' or '-' strand
		}else if(t1->start_codon && !(t2->start_codon)){
			if (t1->strand == '-'){
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					if (t2->exon_list.back().to == (*it).to){
						return true;
					}
				}
			}else if (t1->strand == '+'){
				for (list<Exon>::iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
					if (t2->exon_list.front().from == (*it).from){
						return true;
					}
				}
			}else{}		// no '+' or '-' strand*/ return false;
		}else{
			return false;
		}
	}else{
		return false;
	}
	return false;
}

void check_frame_annotation(list<Transcript> &transcript_list){	// (in Arbeit)
	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
		bool in_frame = true;
		for (list<Exon>::iterator it1 = (*it).exon_list.begin(); it1 != (*it).exon_list.end(); it1++){
			if ((*it).frame != ((*it1).from + (*it1).frame) % 3){
				cerr << "ERROR in Frame" << endl;
				in_frame = false;
			}
			else
				cerr << "Correct Frame" << endl;
		}
		cout << in_frame << endl;
	}
}

// check strand:
/*	Transcript* t1,t2;
	if (t1->start_codon && t1->stop_codon){
		if (t1->start_codon > t1->stop_codon){
			if ('-' != t1->strand)
				cerr << "Stranderror" << endl;
		}
		else{
			if ('+' != t1->strand)
				cerr << "Stranderror" << endl;
		}
	}
*/
