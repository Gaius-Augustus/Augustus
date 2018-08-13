/*
 * \file UTRs.cpp
 */

#include "Genomic_Data.hpp"
#include "Supporting_Methods.hpp"
#include "flex_vec.hpp"
#include "Coord_Transform.hpp"
#include "Compute_UTRs.hpp"
#include "UTRs.hpp"
#include "Global.hpp"

#include <vector>
#include <string>
#include <map>

#include <iostream>
#include <fstream>

#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <iterator>

#include <boost/tuple/tuple.hpp>
#include <boost/assign/std/vector.hpp>

std::vector<UTRs::UTR> UTRs::s_all_UTRs;
UTRs::Parameters UTRs::s_parameters;

using namespace std;
using namespace boost;
using namespace boost::assign;



/**
 * @brief Sorts introns in order from smallest start to greatest start position.
 */
bool sorting_introns_s(Genomic_Data::Intron i, Genomic_Data::Intron j) { return i.start < j.start; }

/**
 * @brief Sorts introns in order from smallest end to greatest end position.
 */
bool sorting_introns_e(Genomic_Data::Intron i, Genomic_Data::Intron j) { return i.end < j.end; }

/**
 * @brief Sorts introns depending on the strand. First minus '-' strand and than plus '+' strand.
 */
bool sort_strand_intron(Genomic_Data::Intron i) { return i.strand == "-"; }

/**
 * @brief Sorts introns in order from most possible intron to least possible intron.
 *
 * Sorts introns in order from most possible intron (highest multiplicity value) to least possible intron
 * (lowest multiplicity value).
 */
bool sorting_intron_mult(Genomic_Data::Intron i, Genomic_Data::Intron j) { return i.mult > j.mult; }

/**
 * @brief Sorts coding region boundaries in order from smallest start to greatest start position.
 */
bool sorting_coding_region(Genomic_Data::CRB i, Genomic_Data::CRB j) {
	return i.codon_pos < j.codon_pos; }

/**
 * @brief Sorts coding region boundaries depending on the strand. First minus '-' strand and second plus '+' strand.
 */
bool sort_strand_coding_region(Genomic_Data::CRB i) { return i.strand == "-"; }


/**
 * @brief Sorts repeats in order from smallest start to greatest start position.
 */
bool sorting_repeat(Genomic_Data::Repeat i, Genomic_Data::Repeat j) { return i.start < j.start; }



/**
 * @brief Finds first start or stop codon with codon position greater than j.
 */
template <class Coding_Region_Boundary> struct finding_coding_region_greater :
	binary_function<Genomic_Data::CRB,unsigned,bool> {
		bool operator() (const Genomic_Data::CRB& i, const unsigned j) const
		{ return i.codon_pos > j; }
};

/**
 * @brief Finds first start or stop codon with codon position smaller than j.
 */
template <class Coding_Region_Boundary> struct finding_coding_region_lesser :
	binary_function<Genomic_Data::CRB,unsigned,bool> {
		bool operator() (const Genomic_Data::CRB& i, const unsigned j) const
		{ return i.codon_pos < j; }
};

/**
 * @brief Finds first intron with a start position greater than j.
 */
template <class Intron> struct finding_intron_greater_s : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.start > j; }
};

/**
 * @brief Finds first intron with a start position smaller than j.
 */
template <class Intron> struct finding_intron_lesser_s : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.start < j; }
};

/**
 * @brief Finds first intron with end position greater than j.
 */
template <class Intron> struct finding_intron_greater_e : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.end > j; }
};

/**
 * @brief Finds first intron with end position smaller than j.
 */
template <class Intron> struct finding_intron_lesser_e : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{
	  return i.end < j; }
};

/**
 * @brief Finds first repeat with a start position greater than j.
 */
template <class Repeat> struct finding_repeat : binary_function<Genomic_Data::Repeat,unsigned,bool> {
	bool operator() (const Genomic_Data::Repeat& i, const unsigned j) const
	{ return i.start > j; }
};



void UTRs::clear() {
	s_all_UTRs.clear();
}



void UTRs::remove_overlapping_introns(vector<Genomic_Data::Intron>* introns) {
	//cannot be empty (because otherwise this function will not be called)
	sort (introns->begin(), introns->end(), sorting_intron_mult);
	vector<int> to_delete(introns->size(),0);
	for (unsigned i = 0; i<introns->size(); ++i){
	  for(unsigned j = i+1; j<introns->size(); ++j){
	    if(to_delete[i] != 1 && to_delete[j] != 1){
	      bool remove = Genomic_Data::Intron::is_overlapping((*introns)[i], (*introns)[j]);
	      if (remove) {
		if((*introns)[i].mult > (*introns)[j].mult){
		  to_delete[j] = 1;
		}else{
		  to_delete[i] = 1;
		}
	      }
	    }
	  }
	}
	int erase_offset = 0;
	for (unsigned i = 0; i< to_delete.size(); ++i){
	  if(to_delete[i] == 1){
	    introns->erase (introns->begin()+i-erase_offset);
	    erase_offset++;
	  }
	}
	
}



unsigned UTRs::get_max_pos(unsigned pos, int dir, unsigned scaffold_size, unsigned read_length,
		vector<Genomic_Data::CRB>* orfs) {
	vector<Genomic_Data::CRB>::iterator it;
	it =
		(dir == 1) ? find_if(orfs->begin(), orfs->end(),
							bind2nd(finding_coding_region_greater<Genomic_Data::CRB>(), pos)) :
					 find_if(orfs->begin(), orfs->end(),
							 bind2nd(finding_coding_region_lesser<Genomic_Data::CRB>(), pos));

	if (it == orfs->end()) {
		return (dir == 1) ? scaffold_size - read_length : read_length + 1; //END
	}

	return (unsigned)((*it).codon_pos - dir); //cannot be negative, because gene found after POS
}



bool UTRs::repeat_in_UTR(unsigned start, unsigned end, vector<Genomic_Data::Repeat>* repeats) {
	vector<Genomic_Data::Repeat>::iterator it, it_pre;
	it = find_if(repeats->begin(), repeats->end(), bind2nd(finding_repeat<Genomic_Data::Repeat>(), start));
	it_pre = it - 1 ; //predecessor of it

	if (repeats->empty()) {
		return false; //if no repeats, then no repeat in interval
	}

	if ( it == repeats->end()) { 	//no repeat with start greater than START found,
		                            //but the repeat before could be in the interval
		 return (*it_pre).end >= start; //it_pre exists, because repeats not empty an it == repeats->end
	}									//end of it_pre is greater than START -> repeat in interval

	if ( ((*it).start <= end)) { //start of repeat it is greater than START and end is smaller than END
		                         //-> repeat in interval
		return true;
	}

	if ( it != repeats->begin()) { //if it == repeats->begin(), than it_pre is not defined
		return (*it_pre).end >= start; //end of it_pre is greater than START -> repeat in interval
	}

	return false;

}



vector<int> UTRs::find_introns_in_range(unsigned start, unsigned end, int dir, vector<Genomic_Data::Intron>* introns) {
	vector<int> indices;

	if (introns->empty()) {
		return indices;
	}

	vector<Genomic_Data::Intron>::iterator it;
	if (dir == 1 ) { //left -> right
		if ((*introns->begin()).strand == "+") {
			it = find_if(introns->begin(), introns->end(), bind2nd(finding_intron_greater_s<Genomic_Data::Intron>(), start));
			while ( it != introns->end() && (*it).start < end) {
				indices += it - introns->begin();
				++it;
			}
		}
		else {
			it = find_if(introns->begin(), introns->end(), bind2nd(finding_intron_greater_e<Genomic_Data::Intron>(), start));
			while ( it != introns->end() && (*it).end < end) {
				indices += it - introns->begin();
				++it;
			}
		}
	}
	else { //right -> left
		vector<Genomic_Data::Intron> reverse_introns = *introns;
		reverse(reverse_introns.begin(), reverse_introns.end()); //done this in COORD Transform

		if ((*reverse_introns.begin()).strand == "+") {
			it = find_if(reverse_introns.begin(), reverse_introns.end(), bind2nd(finding_intron_lesser_e<Genomic_Data::Intron>(), start));
			while ( it != reverse_introns.end() && (*it).end > end) {
				indices += it - reverse_introns.begin();
				++it;
			}
		}
		else {
			it = find_if(reverse_introns.begin(), reverse_introns.end(), bind2nd(finding_intron_lesser_s<Genomic_Data::Intron>(), start));
			while ( it != reverse_introns.end() && (*it).start > end) {
				indices += it - reverse_introns.begin();
				++it;
			}
		}
		// reverse indices (because they refer to reverse_introns, but later, the indices are used to access unreversed introns
		vector<int> reversed_indices;
		for(unsigned w=0; w<indices.size(); ++w){
		  reversed_indices += reverse_introns.size()-indices[w]-1;
		}
		indices = reversed_indices;
	}
	return indices;
}



void UTRs::compute_UTRs(const Genomic_Data::Scaff_plus_Gen_Data& curr_scaffold, const map<unsigned,double>& wiggle_data) {
	Genomic_Data::Scaff_plus_Gen_Data curr_genome = curr_scaffold;

	sort(curr_genome.introns.begin(), curr_genome.introns.end(), sorting_introns_s);
	vector<Genomic_Data::Intron>::iterator intron_split_idx =
			stable_partition (curr_genome.introns.begin(), curr_genome.introns.end(), sort_strand_intron);
	const vector<Genomic_Data::Intron> INTRON_MINUS (curr_genome.introns.begin(), intron_split_idx);
	const vector<Genomic_Data::Intron> INTRON_PLUS (intron_split_idx, curr_genome.introns.end());

	sort(curr_genome.crbs.begin(), curr_genome.crbs.end(), sorting_coding_region);
	vector<Genomic_Data::CRB>::iterator orf_split_idx =
			stable_partition (curr_genome.crbs.begin(), curr_genome.crbs.end(), sort_strand_coding_region);
	vector<Genomic_Data::CRB> orf_minus (curr_genome.crbs.begin(), orf_split_idx);
	vector<Genomic_Data::CRB> orf_plus (orf_split_idx, curr_genome.crbs.end());

	sort(curr_genome.repeats.begin(), curr_genome.repeats.end(), sorting_repeat);

	const map<unsigned,double>& wd = wiggle_data; //wd[i] = c[i]
	const unsigned GENOME_SIZE = curr_genome.sequence.size();
	for (unsigned i = 0; i < curr_genome.crbs.size(); ++i) {
		//setting direction of computation, left -> right 1 or right -> left -1
		const int DIR =
			( (curr_genome.crbs[i].feature == "stop") && (curr_genome.crbs[i].strand == "+") ) ||
			( (curr_genome.crbs[i].feature == "start") && (curr_genome.crbs[i].strand == "-") ) ? 1 : -1;
		vector<Genomic_Data::Intron> curr_introns =
			(curr_genome.crbs[i].strand == "-") ? INTRON_MINUS : INTRON_PLUS;
		vector<Genomic_Data::CRB> curr_orfs =
			(curr_genome.crbs[i].strand == "-") ? orf_minus : orf_plus;
		if (DIR == -1) {
			reverse(curr_orfs.begin(), curr_orfs.end());
		}

		const unsigned POS = curr_genome.crbs[i].codon_pos;
		//starting window size in the gene so the computed values are more accurate
		//======== W out of bounds ==========
		//W is to big, POS - W*DIR outside of genome
		if ( (int)(POS -s_parameters.window_size*DIR) <= 0 || (int)(POS - s_parameters.window_size*DIR) >= (int)GENOME_SIZE) {
			continue;
		}

		//======= RESETTING LIMIT= ==========	    
		const unsigned MAX_POS = get_max_pos(POS, DIR, GENOME_SIZE, s_parameters.read_length, &curr_orfs);

		//=========== REMOVING OVERLAPPING INTRONS ===========
		//		cout << "Before filtering, we have the following introns in the list:" << endl;
		//for(unsigned v=0; v< curr_introns.size();++v){
		//  cout << "start " << curr_introns[v].start << " end " << curr_introns[v].end << " strand " << curr_introns[v].strand << endl;
		//}
		vector<int> poss_intron_idx = find_introns_in_range(POS, MAX_POS, DIR, &curr_introns);
		//cout << POS << " " << MAX_POS << " " << DIR << endl;
		if (poss_intron_idx.empty()) {
			curr_introns.clear(); //no introns in the computation range -> no introns overlap in this range
		}
		else {
			vector<Genomic_Data::Intron> new_introns;
			for (unsigned i = 0; i < poss_intron_idx.size(); ++i ) {
				new_introns += curr_introns[poss_intron_idx[i]];
			}
			curr_introns = new_introns;
			if (curr_introns.size() > 1) { //with only one intron, no overlapping introns possible
			  //for(unsigned q=0; q<curr_introns.size();++q){
			  //  cout << curr_introns[q].start << " " << curr_introns[q].end << endl;
			  //}
				remove_overlapping_introns(&curr_introns);
				sort(curr_introns.begin(), curr_introns.end(), sorting_introns_s);
			}
			//			cout << "after filtering we have the following introns in the list:" << endl;
			//for(unsigned v=0; v< curr_introns.size();++v){
			//  cout << "start " << curr_introns[v].start << " end " << curr_introns[v].end << " strand " << curr_introns[v].strand << endl;
			//}

		}

		//========= COORD TRANSFORMATION =============
		Coord_Transform Coord_Trans (POS, DIR, s_parameters.limit, MAX_POS, s_parameters.window_size, GENOME_SIZE, &curr_introns);

		//================================ COMPUTATION =========================================

		//==================== WINDOW ==============================
		if (Coord_Trans.get_last_vec_pos() <= (int)s_parameters.window_size) {
			//cout << "UTR cannot be long enough, because computation range is not at least W long!" << endl;
			continue;
		}
		Flex_Vec<double> c(-s_parameters.window_size-1, s_parameters.limit-1);
		c[-s_parameters.window_size-1] = 0;
		for (int j = -s_parameters.window_size; j <= (int)Coord_Trans.get_last_vec_pos(); j++) {
			if (wd.count( Coord_Trans.pos_transform(j) ) == 1) { //position exits
				c[j] = wd.find( Coord_Trans.pos_transform(j) )->second + c[j-1]; //position j in vector, position pos+k*dir in genome
			}
			else {
				if (wd.count( Coord_Trans.pos_transform(j) ) > 1) {
					ERR_STRM << "Error in wiggle file: Duplicate Wiggle Position!" << endl;
					cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
					abort();
				}
				c[j] = c[j-1]; //position j in vector, position pos+j*dir in genome
			}
		}

		//==================== ZERO COV ========================
		int j; //transcript end

		if (!s_parameters.zero_cov) {
			vector<double> d;
			for (unsigned k = 0; k <= Coord_Trans.get_last_vec_pos()-s_parameters.window_size; ++k) {
				d += (c[k+s_parameters.window_size] + c[k-s_parameters.window_size] - 2*c[k]) / (s_parameters.window_size * s_parameters.window_size);   // pow(w, 2)
			}

			//point of steepest decline in coverage
			//most likely transcript end
			vector<double>::iterator jt = min_element(d.begin(), d.end());
			j = jt - d.begin();
		}
		else {
			//==================== ZERO COV ========================
			bool zero_cov_found = false;
			j = 0; //With Coord Transform, no Introns
			while (!zero_cov_found && (j <= (int)Coord_Trans.get_last_vec_pos()) ) {
				if (wd.count( Coord_Trans.pos_transform(j) ) == 1) { //position exits -> not Zero
					++j;
				}
				else {
					zero_cov_found = true;
				}
			}

			if (!zero_cov_found) {
			//	cout << "No zero coverage found for scaffold " << curr_genome.name
			//		 << " and start or stop codon at position " << POS << "." << endl;
				continue;
			}
		}

		//============================================= EVALUATION =======================================

		//========= UTR MIN LENGTH TEST ===============
		if ( abs( (int)(Coord_Trans.pos_transform(j) - POS) ) < (int)s_parameters.min_length ) {
			//cout << "Discard UTR, because UTR is not at least MIN_LENGTH long!" << endl;
			continue;
		}

		//========= UTR LENGTH TEST =========
		//shortening UTR, because UTR+V is longer than computation range
		if ((int)(j + s_parameters.drop_window_size) >  Coord_Trans.get_last_vec_pos()) {
			if ( (int)(Coord_Trans.get_last_vec_pos()-s_parameters.drop_window_size) > 0) {
				j = Coord_Trans.get_last_vec_pos()-s_parameters.drop_window_size;
			}
			else {
				//cout << "Discard UTR because UTR is not at least V long!" << endl;
				continue;
			}
		}

		double awc; //average window coverage
		double aec; //average exon UTR coverage (dependent on coordination transform? Maybe this doesn't work with introns in
		            //an UTR, because of the window size (not tested, only tested without introns)

		//average coverage in window j+1...j+V is
		awc =  (c[j + s_parameters.drop_window_size]-c[j])/s_parameters.drop_window_size;
		//compute average coverage in window 0..j (excluding introns)
		aec = (c[j]-c[0])/j;

		//========= MINIMAL AVERAGE COVERAGE UTR ================================0
		if (aec < (int)s_parameters.min_average_cov) {
			//Discard UTR, because average coverage of UTR is to low!
			continue;
		}

		//======== AVERAGE WINDOW COVERAGE TEST ================
		if (aec*s_parameters.p_win < awc) {
			//cout << "Discarding UTR, because average coverage in window j+1...j+W is not at most P percent of the coverage in 0..i!" << endl;
			continue;
		}

		//======== REPEAT TEST ===============
		bool repeat_ex = repeat_in_UTR(Coord_Trans.pos_transform(0)-s_parameters.window_size, Coord_Trans.pos_transform(j)+s_parameters.window_size, &curr_genome.repeats);
		if (repeat_ex) {
			//cout << "Discarding UTR, because a repeat exists in range 0...j!" << endl;
			continue;
		}

		//========= INTRON TEST ===============
		vector<int> intron_indices = find_introns_in_range(POS, Coord_Trans.pos_transform(j), DIR, &curr_introns);
		//could be less than in the interval (POS, POS+TRUE_LIMIT)
		for (unsigned l = 0; l < intron_indices.size(); ++l) {
			unsigned start =
				(curr_genome.crbs[i].strand == "+") ? curr_introns[intron_indices[l]].start :
				         											curr_introns[intron_indices[l]].end;
			unsigned end =
				(curr_genome.crbs[i].strand == "+") ? curr_introns[intron_indices[l]].end :
					        										curr_introns[intron_indices[l]].start;

			double sum_int_cov = 0;
			for (unsigned m = start; m <= end; ++m) {
				if (wd.count(m) == 1) { //position exits
					sum_int_cov += wd.find(m)->second;
				}
				else {
					if (wd.count(m) > 1) {
						ERR_STRM << "Error in wiggle file: Duplicate Wiggle Position!" << endl;
						cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
						abort();
					} //else 0
				}
			}

			double aic = sum_int_cov / (end-start+1);//average intron coverage
			//average intron <= P_INT*average exon coverage
			if (aec*s_parameters.p_int >= aic) {

				//mult >= P_MULT*average exon coverage
				if (aec*s_parameters.p_mult > curr_introns[intron_indices[l]].mult) {
					//cout << "Discarding UTR, because multiplicity of intron is to low!" << endl;
					continue;
				}
			}
			else {
				//cout << "Discarding UTR, because average coverage of intron is to high!" << endl;
				continue;
			}
		}

		//============= TEST N-Region =================
		char char_at_j = curr_genome.sequence[Coord_Trans.pos_transform(j)-1];
		if (char_at_j == 'N') {
			//cout << "Discarding UTR, because at j is a N-region!" << endl;
			continue;
		}

		//============ UTR prepared for output ==================
		//UTR not discarded
		UTRs::UTR this_utr;

		this_utr.name = curr_genome.name;
		this_utr.start =
				(DIR == 1) ? POS + 1 : Coord_Trans.pos_transform(j);
		this_utr.end =
				(DIR == 1) ? Coord_Trans.pos_transform(j) : POS - 1;
		this_utr.feature =
			(curr_genome.crbs[i].feature == "start") ? "5'-UTR" : "3'-UTR";
		this_utr.strand = curr_genome.crbs[i].strand;
		this_utr.group =  curr_genome.crbs[i].grp_att;
		this_utr.all_introns.clear();
		for (unsigned n = 0; n < intron_indices.size(); ++n) {
			this_utr.all_introns += curr_introns[intron_indices[n]];
		}
		//========== TEST whether UTR begins or ends with an intron  ================
		// It is - in theory - possible that a UTR begins with an intron but
		// for simplicity, we discard those introns, too (ending in an intron is
		// impossible).
		//cout << "The potential end of this UTR is : "	<< this_utr.end	<< endl;
		bool ends_in_intron = false;
		bool begins_in_intron = false;
		  for (unsigned m = 0; m < this_utr.all_introns.size(); ++m){
		  //cout << "An intron from " << this_utr.all_introns[m].start << " to " << this_utr.all_introns[m].end << endl;
		  if(this_utr.end == this_utr.all_introns[m].end){
		    ends_in_intron = true;
		  }else if(this_utr.start >= this_utr.all_introns[m].start){
		    begins_in_intron = true;
		  }
		}
		if((ends_in_intron == true) || (begins_in_intron == true)){
		  continue;
		}
		//UTR not discarded
		UTRs::s_all_UTRs += this_utr;

	}
}



unsigned UTRs::size() {
	return s_all_UTRs.size();
}



void UTRs::output() {
	ofstream output(s_parameters.output_fname.c_str());

	for (unsigned i = 0; i < s_all_UTRs.size(); ++i) {

		//no introns in UTR, only one output line
		if (s_all_UTRs[i].all_introns.empty()) {
			output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << s_all_UTRs[i].feature << "\t" << s_all_UTRs[i].start << "\t"
				   << s_all_UTRs[i].end << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;
		}
		else {
			unsigned start; //start position in gff file

			//distinction between the + and - strand, because of the different orientation of introns start and end position
			//output for + strand and existing introns
			if (s_all_UTRs[i].strand == "+") {
				sort(s_all_UTRs[i].all_introns.begin(), s_all_UTRs[i].all_introns.end(), sorting_introns_s);  //sorting introns into right order
				start = s_all_UTRs[i].start;

				//for each intron two output lines, til intron and then the intron itself
				for (unsigned j = 0; j < s_all_UTRs[i].all_introns.size(); ++j) {
					output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << s_all_UTRs[i].feature << "\t" <<  start << "\t"
						   << s_all_UTRs[i].all_introns[j].start-1 << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;
					output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << "intron" << "\t" <<  s_all_UTRs[i].all_introns[j].start << "\t"
						   << s_all_UTRs[i].all_introns[j].end << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;

					start = s_all_UTRs[i].all_introns[j].end+1;
				}
			}
			else {//output for - strand and existing introns
				sort(s_all_UTRs[i].all_introns.begin(), s_all_UTRs[i].all_introns.end(), sorting_introns_e); //sorting introns into right order
				start = s_all_UTRs[i].start;

				for (unsigned j = 0; j < s_all_UTRs[i].all_introns.size(); ++j) {
					output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << s_all_UTRs[i].feature << "\t" <<  start << "\t"
						   << s_all_UTRs[i].all_introns[j].end-1 << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;
					output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << "intron" << "\t" <<  s_all_UTRs[i].all_introns[j].end << "\t"
						   << s_all_UTRs[i].all_introns[j].start << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;

					start = s_all_UTRs[i].all_introns[j].start+1;
				}
			}
			//rest of UTR after the last intron
			output << s_all_UTRs[i].name << "\t" << "makeUTR" << "\t" << s_all_UTRs[i].feature << "\t" <<  start << "\t"
				   << s_all_UTRs[i].end << "\t" << "." << "\t" << s_all_UTRs[i].strand << "\t" << "." << "\t" << s_all_UTRs[i].group << endl;
		}
	}

	cout << "Finished writing output file!" << endl;
}
