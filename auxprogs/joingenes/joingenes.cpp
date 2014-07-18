#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

int overallnumber=0;		// only for semantic tests
int overallnumber2=0;		// only for semantic tests
list<Point> points;			// only for semantic tests
int count1 = 0;				// only for semantic tests
int value = 0;				// only for semantic tests

void seek_overlaps(list<Transcript> &transcript_list, string &outfilename){
// devide a transcript list in overlaps and start works at these overlaps
	transcript_list.sort();
	list<Transcript> new_transcripts;
	list<Transcript*> overlap;

	fstream outfile;
	outfile.open(outfilename, ios::out);		// delete content of file filename
	outfile.close();

	string filename2 = "/home/lars/lars/test_data/eval_test.txt";
	fstream outfile2;
	outfile2.open(filename2, ios::out);
	outfile2.close();

	int max_base = (*transcript_list.front().hbt);

	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
		if ((*(*it).lbt) < max_base){
			overlap.push_back(&*it);
			if (max_base < (*(*it).hbt)){
				max_base = (*(*it).hbt);
			}
		}
		else{
			//eval_gtf(overlap);

			work_at_overlap(overlap, new_transcripts);
			save_overlap(overlap, outfilename);
			overlap.clear();
			max_base = (*(*it).hbt);
			overlap.push_front(&*it);
		}
	}
	// scatter_plot(points);			// depends to eval_gtf(overlap) and is only for tests about generangesites
	// cout << "Average error of single species prediction: " << (float)value/count1 << endl;		// only for prediction evaluation
}

void work_at_overlap(list<Transcript*> &overlap, list<Transcript> &new_transcripts)
{
// calls methods for joining transcripts with the target that most of transcripts are complete (have start and stop codon) and delete doublings and other unwanted transcripts from the overlap

	search_n_destroy_doublings(overlap);

	list<Transcript*> new_overlap_part_stop;
	join(overlap, new_transcripts, new_overlap_part_stop, '3');
	overlap.merge(new_overlap_part_stop);

	search_n_destroy_doublings(overlap);

	list<Transcript*> new_overlap_part_start;
	join(overlap, new_transcripts, new_overlap_part_start, '5');
	overlap.merge(new_overlap_part_start);

	search_n_destroy_doublings(overlap);

	overlap.sort(compare_priority);
	int highest_complete_priority = 0;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if ((*it)->tl_complete.first && (*it)->tl_complete.second){
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
			else if ((!(*it)->tl_complete.first || !(*it)->tl_complete.second) && (highest_complete_priority == (*it)->priority)){
				it = overlap.erase(it);
				it--;
			}
		}
	}

	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if (!check_frame_annotation(*(*it))){
			cout << "After joining: " << endl;
			cerr << (*it)->t_id << " has wrong frame annotation." << endl;
			//display_error("Frames are not correct.");
		}
	}
}

void search_n_destroy_doublings(list<Transcript*> &overlap){
// delete all transcripts that are completly part of another transcript (in particular all exons are also in the other transcript); in case of equality the one with the lesser priority will be deleted
	if (overlap.size() > 1){
		for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
			list<Transcript*>::iterator it_temp = it;
			it_temp++;
			for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
				pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
				if (who_is_part.first == true){
					if (who_is_part.second == true){
						(*it)->supporter.push_front(*it_inside);
						(*it_inside)->supporter.push_front(*it);
						if ((*it)->priority < (*it_inside)->priority){
							it = overlap.erase(it);
							it_inside = it;
						}else{
							it_inside = overlap.erase(it_inside);
							it_inside--;
						}
					}else{
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
}

void join(list<Transcript*> &overlap, list<Transcript> &new_transcripts, list<Transcript*> &new_overlap_part, char side){
// devides an overlap in a list of start/stop codon donors and in an other list of start/stop codon acceptors and joins every pair (one of each list) if they are combinable
	list<Transcript*> donor;
	list<Transcript*> acceptor;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if (side == '3'){
			if ((*it)->tl_complete.second){
				donor.push_front(*it);
			}else{
				acceptor.push_front(*it);
			}
		}else if (side == '5'){
			if ((*it)->tl_complete.first){
				donor.push_front(*it);
			}else{
				acceptor.push_front(*it);
			}
		}
	}
	for (list<Transcript*>::iterator it = acceptor.begin(); it != acceptor.end(); it++){
		for (list<Transcript*>::iterator it_donor = donor.begin(); it_donor != donor.end(); it_donor++){
			if ((*it)->strand == (*it_donor)->strand){
				int fitting_case = is_combinable(*it, *it_donor, (*it)->strand, side);
				if (fitting_case){
					Transcript tx_new;
					joining(*it, *it_donor, (*it)->strand, side, tx_new, fitting_case);
					new_transcripts.push_back(tx_new);
					new_overlap_part.push_back(&new_transcripts.back());
				}
			}
		}
	}
}

void joining(Transcript* t1, Transcript* t2, char strand, char side, Transcript &tx_new, int fitting_case){
// joins transcripts in one direction so that every suitable exon will be transferred and returns a new "joined" transcript without deleting the old ones
	int minimum_intron_length = 20;								// also defined in is_combinable
	tx_new = *t1;
//	cout << "----------- " << fitting_case << endl;
//	cout << tx_new.exon_list.size() << endl;
	if (strand == '+'){
		tx_new.lbt = &tx_new.tis;
		tx_new.hbt = &tx_new.tes;
	}else{
		tx_new.lbt = &tx_new.tes;
		tx_new.hbt = &tx_new.tis;
	}
	bool are_at_add_part;
	list<Exon> temp_exon_list;
	temp_exon_list.clear();
	//cout << temp_exon_list.size() << endl;
	switch (fitting_case){
		case 0:
			cerr << "WARNING: shouldnt happen (in joining())!" << endl;
			break;
		case 1:
			if (strand == '+'){
				tx_new.tes = t2->tes;
				tx_new.tl_complete.second = t2->tl_complete.second;
			}else{
				tx_new.tis = t2->tis;
				tx_new.tl_complete.first = t2->tl_complete.first;
			}
			are_at_add_part = false;
			for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
				if (are_at_add_part){
					tx_new.exon_list.push_back(*it);
				}
				if ((*t1).exon_list.back().to >= (*it).from && (*t1).exon_list.back().to <= (*it).to){
					tx_new.exon_list.back().to = (*it).to;
					if (strand == '-'){tx_new.exon_list.back().frame = (*it).frame;}
					are_at_add_part = true;
				}
			}
			break;
		case 2:
			if (strand == '+'){
				tx_new.tes = t2->tes;
				tx_new.tl_complete.second = t2->tl_complete.second;
			}else{
				tx_new.tis = t2->tis;
				tx_new.tl_complete.first = t2->tl_complete.first;
			}
			are_at_add_part = false;
			for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
				if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
					are_at_add_part = true;
				}
				if (are_at_add_part){
					tx_new.exon_list.push_back(*it);
				}
			}
			break;
		case 3:
			if (strand == '+'){
				tx_new.tis = t2->tis;
				if (t2->tl_complete.first){tx_new.tl_complete.first = true;}
			}else{
				tx_new.tes = t2->tes;
				if (t2->tl_complete.second){tx_new.tl_complete.second = true;}
			}
			//are_at_add_part = true;
			for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
				if ((*t1).exon_list.front().from >= (*it).from && (*t1).exon_list.front().from <= (*it).to){
					tx_new.exon_list.front().from = (*it).from;
					if (strand == '+'){tx_new.exon_list.front().frame = (*it).frame;}
					tx_new.exon_list.merge(temp_exon_list);
					//are_at_add_part = false;
					break;
				}
				//if (are_at_add_part){
					temp_exon_list.push_back(*it);
				//}
			}
			break;
		case 4:
			if (strand == '+'){
				tx_new.tis = t2->tis;
				if (t2->tl_complete.first){tx_new.tl_complete.first = true;}
			}else{
				tx_new.tes = t2->tes;
				if (t2->tl_complete.second){tx_new.tl_complete.second = true;}
			}
			//are_at_add_part = true;
			for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
				if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
					tx_new.exon_list.merge(temp_exon_list);
					//are_at_add_part = false;
					break;
				}
				//if (are_at_add_part){
					temp_exon_list.push_back(*it);
				//}
			}
			break;
		default:
			cerr << "WARNING: unexpected case (in joining())!" << endl;
	}
	//cout << tx_new.exon_list.size() << endl;
}

int is_combinable(Transcript const* t1, Transcript const* t2, char strand, char side){
// is true,	if the first exon which is minimal the minimum_intron_length away from the last exon of the other transcript in the appropriate direction
// 			&& the exons at these positions are frame-compatible
// 			&& the transcripts are overlapping					// maybe we can improve something here, that combinable non-overlaping transcripts gets true (but be carefull)
	int minimum_intron_length = 20;								// also defined in joining
	if ((*(*t1).hbt < *(*t2).lbt) || (*(*t1).lbt > *(*t2).hbt)){
		return 0;
	}
	if ((strand == '+' && side == '3') || (strand == '-' && side == '5')){
		for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if (((*t1).exon_list.back().to >= (*it).from) && ((*t1).exon_list.back().to <= (*it).to)){
				if ((*t1).exon_list.back().frame == -1 && (*it).frame == -1){return true;}
				else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return false;}
				if ((strand == '+') && ((3 - (((*t1).exon_list.back().to - (*t1).exon_list.back().from) - (*t1).exon_list.back().frame) % 3) % 3) == ((3 - (((*t1).exon_list.back().to - (*it).from) - (*it).frame) % 3) % 3) ){
					return 1;
				}else if ((strand == '-') && (*t1).exon_list.back().frame == ((3 - (((*it).to - (*t1).exon_list.back().to) - (*it).frame) % 3) % 3) ){
					return 1;
				}else{
					return 0;
				}
			}else{
				if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
					if ((strand == '+') && ((*it).frame == (3 - ( ((*t1).exon_list.back().to - (*t1).exon_list.back().from + 1) - (*t1).exon_list.back().frame) % 3) % 3)){
						return 2;
					}else if ((strand == '-') && ((*t1).exon_list.back().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
						return 2;
					}else{
						return 0;
					}
				}
			}
		}
	}else if ((strand == '+' && side == '5') || (strand == '-' && side == '3')){
		for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
			if (((*t1).exon_list.front().from >= (*it).from) && ((*t1).exon_list.front().from <= (*it).to)){
				if ((*t1).exon_list.front().frame == -1 && (*it).frame == -1){return true;}
				else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return false;}
				if ( (strand == '-') && ((3 - (((*t1).exon_list.front().to - (*t1).exon_list.front().from) - (*t1).exon_list.front().frame) % 3) % 3) == ((3 - (((*it).to - (*t1).exon_list.front().from) - (*it).frame) % 3) % 3) ){
					return 3;
				}else if ( (strand == '+') && (*t1).exon_list.front().frame == ((3 - (((*t1).exon_list.front().from - (*it).from) - (*it).frame) % 3) % 3) ){
					return 3;
				}else{
					return 0;
				}
			}else{
				if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
					if (it != t2->exon_list.begin()){
						it--;
					}else{
						return 0;
					}
					if ((strand == '-') && ((*it).frame == (3 - ( ((*t1).exon_list.front().to - (*t1).exon_list.front().from + 1) - (*t1).exon_list.front().frame) % 3) % 3)){
						return 4;
					}else if ((strand == '+') && ((*t1).exon_list.front().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
						return 4;
					}else{
						return 0;
					}
				}
			}
		}
	}else{									// no '+' or '-' strand
		return 0;						// here we have an unexpected case and so we dont wanna combine
	}
	return 0;
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
		if (t1->tl_complete.first && t2->tl_complete.first && t1->tis != t2->tis){
			t1_is_part = false;
			t2_is_part = false;
		}
		if (t1->tl_complete.second && t2->tl_complete.second && t1->tes != t2->tes){
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
// returns true, if the priority of the first element is higher; used to sort transcripts by priority
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
			if ((*it1).frame != -1 && (*it2).frame != -1){
				if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
					return false;
				}
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
			if ((*it1).frame != -1 && (*it2).frame != -1){
				if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
					return false;
				}else
				if (it2 == transcript.exon_list.begin())
					break;
			}
		}
	}
	return true;
}

void eval_gtf(list<Transcript*> &overlap){
	string filename = "/home/lars/lars/test_data/eval_test.txt";
	fstream outfile;
	outfile.open(filename, ios::out | ios::app);

	list<Transcript*> annotation;
	list<Transcript*> prediction;
	list<Transcript*> single;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if ((*it)->priority == 3){
			annotation.push_back(*it);
		}
		search_n_destroy_doublings(annotation);
		if ((*it)->pred_range.first && (*it)->pred_range.second){
			prediction.push_back(*it);
		}
		search_n_destroy_doublings(prediction);
		if ((*it)->priority == 1){
			single.push_back(*it);
		}
	}
	
	if (annotation.empty()){
		for (list<Transcript*>::iterator itp = prediction.begin(); itp != prediction.end(); itp++){
			for (list<Exon>::iterator itpe = (*itp)->exon_list.begin(); itpe != (*itp)->exon_list.end(); itpe++){
				(*itpe).toless = 0;
				(*itpe).tomany = (*itpe).to - (*itpe).from + 1;

				(*itpe).penalty = (*itpe).toless + (*itpe).tomany;
				(*itpe).distance = min((*itpe).from - (*itp)->pred_range.first, (*itp)->pred_range.second - (*itpe).to);
				outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
			}
		}
		for (list<Transcript*>::iterator its = single.begin(); its != single.end(); its++){
			for (list<Exon>::iterator itpe = (*its)->exon_list.begin(); itpe != (*its)->exon_list.end(); itpe++){
				(*itpe).toless = 0;
				(*itpe).tomany = (*itpe).to - (*itpe).from + 1;

				(*itpe).penalty = (*itpe).toless + (*itpe).tomany;
				count1++;
				value += (*itpe).penalty;
			}
		}
	}else{
		for (list<Transcript*>::iterator itp = prediction.begin(); itp != prediction.end(); itp++){
			for (list<Exon>::iterator itpe = (*itp)->exon_list.begin(); itpe != (*itp)->exon_list.end(); itpe++){
				int min_pen = -1;
				for (list<Transcript*>::iterator ita = annotation.begin(); ita != annotation.end(); ita++){
					for (list<Exon>::iterator itae = (*ita)->exon_list.begin(); itae != (*ita)->exon_list.end(); itae++){
						if ((*itpe).to >= (*itae).from && (*itpe).from <= (*itae).to && (*itp)->strand == (*ita)->strand && strandeq(*itpe,*itae,(*itp)->strand)){
							int more = 0;
							int lesser = 0;
							if ((*itpe).to > (*itae).to){
								more += (*itpe).to - (*itae).to;
							}else if ((*itpe).to < (*itae).to){
								lesser += (*itae).to - (*itpe).to;
							}
							if ((*itpe).from > (*itae).from){
								lesser += (*itpe).from - (*itae).from;
							}else if ((*itpe).from < (*itae).from){
								more += (*itae).from - (*itpe).from;
							}
							if (min_pen == -1 || min_pen > more + lesser){
								(*itpe).toless = lesser;
								(*itpe).tomany = more;
								min_pen = (*itpe).toless + (*itpe).tomany;
							}
						}else{
							if (min_pen == -1 || min_pen > (*itpe).to - (*itpe).from + 1){
								(*itpe).toless = 0;
								(*itpe).tomany = (*itpe).to - (*itpe).from + 1;
								min_pen = (*itpe).toless + (*itpe).tomany;
							}
						}
					}
				}
				(*itpe).penalty = min_pen;
				(*itpe).distance = min((*itpe).from - (*itp)->pred_range.first, (*itp)->pred_range.second - (*itpe).to);
				outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
			}
		}
		for (list<Transcript*>::iterator its = single.begin(); its != single.end(); its++){
			for (list<Exon>::iterator itpe = (*its)->exon_list.begin(); itpe != (*its)->exon_list.end(); itpe++){
				int min_pen = -1;
				for (list<Transcript*>::iterator ita = annotation.begin(); ita != annotation.end(); ita++){
					for (list<Exon>::iterator itae = (*ita)->exon_list.begin(); itae != (*ita)->exon_list.end(); itae++){
						if ((*itpe).to >= (*itae).from && (*itpe).from <= (*itae).to && (*its)->strand == (*ita)->strand && strandeq(*itpe,*itae,(*its)->strand)){
							int more = 0;
							int lesser = 0;
							if ((*itpe).to > (*itae).to){
								more += (*itpe).to - (*itae).to;
							}else if ((*itpe).to < (*itae).to){
								lesser += (*itae).to - (*itpe).to;
							}
							if ((*itpe).from > (*itae).from){
								lesser += (*itpe).from - (*itae).from;
							}else if ((*itpe).from < (*itae).from){
								more += (*itae).from - (*itpe).from;
							}
							if (min_pen == -1 || min_pen > more + lesser){
								(*itpe).toless = lesser;
								(*itpe).tomany = more;
								min_pen = (*itpe).toless + (*itpe).tomany;
							}
						}else{
							if (min_pen == -1 || min_pen > (*itpe).to - (*itpe).from + 1){
								(*itpe).toless = 0;
								(*itpe).tomany = (*itpe).to - (*itpe).from + 1;
								min_pen = (*itpe).toless + (*itpe).tomany;
							}
						}
					}
				}
				(*itpe).penalty = min_pen;
				count1++;
				value += (*itpe).penalty;
			}
		}
	}
	outfile.close();
}

bool strandeq(Exon ex1, Exon ex2, char strand){
	if (strand == '+'){
		if (ex1.frame == -1 && ex2.frame == -1){return true;}
		if (ex1.from >= ex2.from)
			if (((ex1.from - ex2.from) + ex1.frame) % 3 == ex2.frame){return true;}else{return false;}
		else
			if (((ex2.from - ex1.from) + ex2.frame) % 3 == ex1.frame){return true;}else{return false;}
	}else{
		if (ex1.frame == -1 && ex2.frame == -1){return true;}
		if (ex1.to >= ex2.to)
			if (((ex1.to - ex2.to) + ex2.frame) % 3 == ex1.frame){return true;}else{return false;}
		else
			if (((ex2.to - ex1.to) + ex1.frame) % 3 == ex2.frame){return true;}else{return false;}
	}
	return false;
}

void output_exon_list(Transcript const* tx){
	cout << ">> Exon_list of " << tx->t_id << ": " << tx->tl_complete.first << " " << tx->tl_complete.second << " " << tx->strand << endl;
	for (list<Exon>::const_iterator it = tx->exon_list.begin(); it != tx->exon_list.end(); it++){
		cout << (*it).from << "\t" << (*it).to << "\t" << (*it).frame << "\t" << (*it).feature << "\t" << tx->tis << "\t" << tx->tes << endl;
	}
}

void scatter_plot(list<Point> input){
	int terminal_high = 62;
	int terminal_length = 180;
	double xmax = 0;
	double ymax = 0;
	int xpos = 0;
	int ypos = terminal_high;
	for (list<Point>::iterator it = input.begin(); it != input.end(); it++){
	//	cout << (*it).i_x << " " << (*it).i_y << endl;
	}

	for (list<Point>::iterator it = input.begin(); it != input.end(); it++){
		if ((*it).i_x > 100)
			(*it).i_x = 100;
		if ((*it).i_y > 100)
			(*it).i_y = 100;
		if (xmax < (*it).i_x)
			xmax = (*it).i_x;
		if (ymax < (*it).i_y)
			ymax = (*it).i_y;
	}
	for (list<Point>::iterator it = input.begin(); it != input.end(); it++){
		(*it).x = (int)((terminal_length * (*it).i_x)/xmax);
		(*it).y = (int)((terminal_high * (*it).i_y)/ymax);
	}
	input.sort();
	for (list<Point>::iterator it = input.begin(); it != input.end(); it++){
	//	cout << (*it).i_x << " " << (*it).i_y << endl;
	}
	cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	for (list<Point>::iterator it = input.begin(); it != input.end(); it++){
		if (ypos > (*it).y){
			for (unsigned int i = 0; i < (ypos - (*it).y); i++){cout << endl;}
			ypos = (*it).y;
			xpos = 0;
		}
		if (xpos < (*it).x){
			for (unsigned int i = 0; i < ((*it).x-xpos); i++){cout << " ";}
			xpos = (*it).x;
			cout << "X";
			xpos++;
		}else if (xpos == (*it).x){
			cout << "X";
			xpos++;
		}
	}
	cout << endl;
}
