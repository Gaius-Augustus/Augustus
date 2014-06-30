#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

// int overallnumber=0;		// only for semantic tests
// int overallnumber2=0;	// only for semantic tests

void seek_overlaps(list<Transcript> &transcript_list, string &outfilename){
// devide a transcript list in overlaps and start works at these overlaps
	transcript_list.sort();
	list<Transcript> new_transcripts;
	list<Transcript*> overlap;
	int max_base = transcript_list.front().stop;
	fstream outfile;
	outfile.open(outfilename, ios::out);		// delete content of file filename
	outfile.close();

	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
/*if ((*it).pred_range.first && (*it).pred_range.second){
cout << (*it).strand << " ----------------------------" << endl;
cout << "Transcript range: " << (*it).start << " " << (*it).stop << " (" << (*it).stop - (*it).start << ")" << endl;
cout << "Prediction range: " << (*it).pred_range.first << " " << (*it).pred_range.second << " (" << (*it).pred_range.second - (*it).pred_range.first << ")" << endl;
cout << "Distances: " << (*it).start - (*it).pred_range.first << " " << (*it).pred_range.second - (*it).stop << endl;
}*/
		if ((*it).start < max_base){
			overlap.push_back(&*it);
			if (max_base < (*it).stop){
				max_base = (*it).stop;
			}
		}
		else{
			{
				work_at_overlap(overlap, new_transcripts);
				save_overlap(overlap, outfilename);
			}
			overlap.clear();
			max_base = (*it).stop;
			overlap.push_front(&*it);
		}
	}
}

void work_at_overlap(list<Transcript*> &overlap, list<Transcript> &new_transcripts)
{
// calls methods for joining transcripts with the target that most of transcripts are complete (have start and stop codon) and delete doublings and other unwanted transcripts from the overlap
	search_n_destroy_doublings(overlap);

	list<Transcript*> new_overlap_part_stop;
	join_stop(overlap, new_transcripts, new_overlap_part_stop);
	overlap.merge(new_overlap_part_stop);

	list<Transcript*> new_overlap_part_start;
	join_start(overlap, new_transcripts, new_overlap_part_start);
	overlap.merge(new_overlap_part_start);

	search_n_destroy_doublings(overlap);

	overlap.sort(compare_priority);
	int highest_complete_priority = 0;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if ((*it)->start_codon && (*it)->stop_codon){
			highest_complete_priority = (*it)->priority;
			break;
		}
	}
	if (highest_complete_priority){
		for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
			if (highest_complete_priority > (*it)->priority){
				it = overlap.erase(it);
				it--;
			}
			else if ((!((*it)->start_codon) || !((*it)->stop_codon)) && (highest_complete_priority == (*it)->priority)){
				it = overlap.erase(it);
				it--;
			}
		}
	}
}

void search_n_destroy_doublings(list<Transcript*> &overlap){
// delete all transcripts that are completly part of another transcript (in particular all exons are also in the other transcript); in case of equality the one with the lesser priority will be deleted
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		list<Transcript*>::iterator it_temp = it;
		it_temp++;
		for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
			if (compare_transcripts(*it, *it_inside)){
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
				}
			}else{
				if (who_is_part.second == true){
					it_inside = overlap.erase(it_inside);
					it_inside--;
				}
			}
		}
	}
}



void join_stop(list<Transcript*> &overlap, list<Transcript> &new_transcripts, list<Transcript*> &new_overlap_part){
// devides an overlap in a list of stop codon donors and in an other list of stop codon acceptors and joins every pair (one of each list) if they are combinable 
	list<Transcript*> donor_stop;
	list<Transcript*> acceptor_stop;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if ((*it)->stop_codon){
			donor_stop.push_front(*it);
		}else{
			acceptor_stop.push_front(*it);
		}
	}
	for (list<Transcript*>::iterator it = acceptor_stop.begin(); it != acceptor_stop.end(); it++){
		for (list<Transcript*>::iterator it_donor = donor_stop.begin(); it_donor != donor_stop.end(); it_donor++){
			if ((*it)->strand == '+' && (*it_donor)->strand == '+'){
				if (is_combinable(*it, *it_donor, (*it)->strand, '3')){
					Transcript tx_new;
					joining(*it, *it_donor, (*it)->strand, '3', tx_new);
					new_transcripts.push_back(tx_new);
					new_overlap_part.push_back(&new_transcripts.back());
				}
			}
			else if ((*it)->strand == '-' && (*it_donor)->strand == '-'){
				if (is_combinable(*it, *it_donor, (*it)->strand, '5')){
					Transcript tx_new;
					joining(*it, *it_donor, (*it)->strand, '5', tx_new);
					new_transcripts.push_back(tx_new);
					new_overlap_part.push_back(&new_transcripts.back());
				}
			}
		}
	}
}

void join_start(list<Transcript*> &overlap, list<Transcript> &new_transcripts, list<Transcript*> &new_overlap_part){
// devides an overlap in a list of start codon donors and in an other list of start codon acceptors and joins every pair (one of each list) if they are combinable 
	list<Transcript*> donor_start;
	list<Transcript*> acceptor_start;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if ((*it)->start_codon){
			donor_start.push_front(*it);
		}else{
			acceptor_start.push_front(*it);
		}
	}
	for (list<Transcript*>::iterator it = acceptor_start.begin(); it != acceptor_start.end(); it++){
		for (list<Transcript*>::iterator it_donor = donor_start.begin(); it_donor != donor_start.end(); it_donor++){
			if ((*it)->strand == '+' && (*it_donor)->strand == '+'){
				if (is_combinable(*it, *it_donor, (*it)->strand, '5')){
					Transcript tx_new;
					joining(*it, *it_donor, (*it)->strand, '5', tx_new);
					new_transcripts.push_back(tx_new);
					new_overlap_part.push_back(&new_transcripts.back());
				}
			}
			else if ((*it)->strand == '-' && (*it_donor)->strand == '-'){
				if (is_combinable(*it, *it_donor, (*it)->strand, '3')){
					Transcript tx_new;
					joining(*it, *it_donor, (*it)->strand, '3', tx_new);
					new_transcripts.push_back(tx_new);
					new_overlap_part.push_back(&new_transcripts.back());
				}
			}
		}
	}
}

void joining(Transcript* t1, Transcript* t2, char strand, char side, Transcript &tx_new){
// joins transcripts in one direction so that every suitable exon will be transferred and returns a new "joined" transcript without deleting the old ones
	int minimum_intron_length = 20;								// also defined in is_combinable
	tx_new = *t1;
	if ((strand == '+' && side == '3') || (strand == '-' && side == '5')){
		tx_new.stop_codon = t2->stop_codon;
		if (strand == '+'){tx_new.stop_codon = t2->stop_codon;}else{tx_new.start_codon = t2->start_codon;}
		tx_new.stop = t2->stop;
		bool are_at_add_part = false;
		for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
				are_at_add_part = true;
			}
			if (are_at_add_part){
				tx_new.exon_list.push_back(*it);
			}
		}
	}else if ((strand == '+' && side == '5') || (strand == '-' && side == '3')){
		if (strand == '+'){tx_new.start_codon = t2->start_codon;}else{tx_new.stop_codon = t2->stop_codon;}
		tx_new.start = t2->start;
		list<Exon> temp_exon_list;
		bool are_at_add_part = true;
		for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
				tx_new.exon_list.merge(temp_exon_list);
				are_at_add_part = false;
				break;
			}
			if (are_at_add_part){
				temp_exon_list.push_back(*it);
			}
		}
	}
}

bool is_combinable(Transcript const* t1, Transcript const* t2, char strand, char side){
// is true,	if the first exon which is minimal the minimum_intron_length away from the last exon of the other transcript in the appropriate direction
// 			&& the exons at these positions are frame-compatible
// 			&& the transcripts are overlapping					// maybe we can improve something here, that combinable non-overlaping transcripts gets true (but be carefull)
	int minimum_intron_length = 20;								// also defined in joining
	if ((t1->stop < t2->start) || (t1->start > t2->stop)){
		return false;
	}
	if ((strand == '+' && side == '3') || (strand == '-' && side == '5')){
		for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
				if ((strand == '+') && ((*it).frame == (3 - ( ((*t1).exon_list.back().to - (*t1).exon_list.back().from + 1) - (*t1).exon_list.back().frame) % 3) % 3)){
					return true;
				}else if ((strand == '-') && ((*t1).exon_list.back().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
					return true;
				}else{
					return false;
				}
			}
		}
	}else if ((strand == '+' && side == '5') || (strand == '-' && side == '3')){
		for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
				if (it != t2->exon_list.begin()){
					it--;
				}else{
					return false;
				}
				if ((strand == '-') && ((*it).frame == (3 - ( ((*t1).exon_list.back().to - (*t1).exon_list.back().from + 1) - (*t1).exon_list.back().frame) % 3) % 3)){
					return true;
				}else if ((strand == '+') && ((*t1).exon_list.back().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
					return true;
				}else{
					return false;
				}
			}
		}
	}else{									// no '+' or '-' strand
		return false;						// here we have an unexpected case and so we dont wanna combine
	}
	return false;
}

bool compare_transcripts(Transcript const* t1, Transcript const* t2)
{
// is true, if both transcripts are equal in exons, strand and start/stop
	if (t1->strand == t2->strand && t1->start == t2->start && t1->stop == t2->stop && t1->exon_list.size() == t2->exon_list.size()){
		list<Exon>::const_iterator it2 = t2->exon_list.begin();
		for (list<Exon>::const_iterator it1 = t1->exon_list.begin(); it1 != t1->exon_list.end(); it1++){
			if ((*it1).from == (*it2).from && (*it1).to == (*it2).to){
			}else{
				return false;
			}
			it2++;
		}
	}else{
		return false;
	}
	return true;
}

pair<bool,bool> is_part_of(Transcript const* t1, Transcript const* t2)
{
// is true,false or false,true if one transcript contains the other completely
// is true,true if the transcripts are equal in exons			// this case could completly replace compare_transcripts
	bool t1_is_part = false;
	bool t2_is_part = false;
	if (t1->strand == t2->strand){
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
		list<Exon>::const_iterator it1 = t1->exon_list.begin();
		list<Exon>::const_iterator it2 = t2->exon_list.begin();
		while (t1_is_part == true || t2_is_part == true){
			if ((*it1).from == (*it2).from){
				if((*it1).to == (*it2).to && (*it1).frame == (*it2).frame){
					it1++;
					it2++;
				}
				else{
					t1_is_part = false;
					t2_is_part = false;
					break;
				}
			}else if ((*it1).from > (*it2).from){
				t2_is_part = false;
				it2++;
			}else{
				t1_is_part = false;
				it1++;
			}
			if (it1 == t1->exon_list.end() && it2 == t2->exon_list.end()){
				break;
			}else{
				if (it1 == t1->exon_list.end() && !(it2 == t2->exon_list.end()))
				{
					t2_is_part = false;
					break;
				}
				if (it2 == t2->exon_list.end() && !(it1 == t1->exon_list.end()))
				{
					t1_is_part = false;
					break;
				}
			}
		}
	}
	return make_pair(t1_is_part, t2_is_part);
}

bool compare_priority(Transcript const* lhs, Transcript const* rhs){
// is true, if the priority of the first element is higher; used to sort elements by priority
	return ( lhs->priority > rhs->priority );
}

bool check_frame_annotation(Transcript const &transcript){
// returns true, if the frame annotation of the transcript is correct
	list<Exon>::const_iterator it2;
	if (transcript.strand == '+'){
		for (list<Exon>::const_iterator it1 = transcript.exon_list.begin(); it1 != transcript.exon_list.end(); it1++){
			it2 = it1;
			it2++;
			if (it2 == transcript.exon_list.end())
				break;
			if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
				return false;
			}
		}
	}else if (transcript.strand == '-'){
		for (list<Exon>::const_iterator it1 = transcript.exon_list.end(); it1 != transcript.exon_list.begin(); it1--){
			if (it1 == transcript.exon_list.end()){
				it1--;
				if (it1 == transcript.exon_list.begin())
					break;
			}
			it2 = it1;
			it2--;
			if (it2 == transcript.exon_list.begin())
				break;
			if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
				return false;
			}else
			if (it2 == transcript.exon_list.begin())
				break;
		}
	}
	return true;
}
