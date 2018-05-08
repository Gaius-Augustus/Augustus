/**
 * \file Coord_Transform.cpp
 */

#include "flex_vec.hpp"
#include "Coord_Transform.hpp"

#include <string>
#include <vector>

#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;

/**
 * @brief Sort introns in order from smallest starting to greatest starting point.
 */
bool sorting_intron(Genomic_Data::Intron i, Genomic_Data::Intron j) { return i.start < j.start; }

/**
 * @brief Find intron with start position greater than j.
 */
template <class Intron> struct finding_intron_greater_s : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.start >= j; }
};

/**
 * @brief Find intron with start position lesser than j.
 */
template <class Intron> struct finding_intron_lesser_s : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.start < j; }
};

/**
 * @brief Find intron with end position greater than j.
 */
template <class Intron> struct finding_intron_greater_e : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.end > j; }
};

/**
 * @brief Find intron with end position lesser than j.
 */
template <class Intron> struct finding_intron_lesser_e : binary_function<Genomic_Data::Intron,unsigned,bool> {
	bool operator() (const Genomic_Data::Intron& i, const unsigned j) const
	{ return i.end < j; }
};



Coord_Transform::Coord_Transform(unsigned comp_start, int dir, unsigned limit, unsigned max_pos, unsigned W,
		unsigned scaffold_size,	vector<Genomic_Data::Intron>* curr_introns) : m_transform_vec(-W, limit-1) {

	m_win_size = W;
	m_last_vec_pos = -W;

	int i = comp_start - W*dir; //positions in the scaffold
	int j = -W; //indices of the vector m_transform_vec
	int skip_start; //normally next intron start position - W; if no introns no skipping
	int skip_end; //normally end position + W for this next intron
	vector<Genomic_Data::Intron>::iterator it; //next intron after comp_start

	if (curr_introns->empty()) {
		skip_start =
				(dir == 1) ? (scaffold_size + 1) : 0; //no skipping due to intron
	}
	else {
		//find first intron after comp_start
		if (dir == 1) { //left -> right
			if ((*curr_introns->begin()).strand == "+") {
				it = find_if(curr_introns->begin(), curr_introns->end(),
						bind2nd(finding_intron_greater_s<Genomic_Data::Intron>(), i));
			}
			else {
				it = find_if(curr_introns->begin(), curr_introns->end(),
						bind2nd(finding_intron_greater_e<Genomic_Data::Intron>(), i));
			}
		}
		else { //right -> left
			reverse(curr_introns->begin(),curr_introns->end());

			if ((*curr_introns->begin()).strand == "+") {
				it = find_if( curr_introns->begin(), curr_introns->end(),
						bind2nd(finding_intron_lesser_e<Genomic_Data::Intron>(), i));
			}
			else {
				it = find_if( curr_introns->begin(), curr_introns->end(),
						bind2nd(finding_intron_lesser_s<Genomic_Data::Intron>(), i));
			}
		}


		if (it == curr_introns->end()) {
			//no intron after comp_start, I think this case will never occur,
			//because an intron in range exists, because vector is not empty
			skip_start =
					(dir == 1) ? (scaffold_size + 1) : 0;
		}
		else {
			//skip_start and skip_end depending on strand and direction
			if ((*curr_introns->begin()).strand == "+") {
				skip_start =
					( dir == 1) ? (*it).start - W : (*it).end + W;
				skip_end =
					( dir == 1) ? (*it).end + W : (*it).start - W;
			}
			else {
				skip_start =
					( dir == 1) ? (*it).end - W : (*it).start + W;
				skip_end =
					( dir == 1) ? (*it).start + W : (*it).end - W;
			}
		}
	}

	while ( (j <= (int)(limit-1)) ) {
		//i in intron range
		if (i*dir >= skip_start*dir) {
			bool intron_found = false; //next intron not found

			while (!intron_found) {
				++it;
				if ( it == curr_introns->end() ) {
					skip_start =
						(dir == 1) ? scaffold_size + 1 : 0; //no more skipping
					i = skip_end + dir;
					intron_found = true; //next intron found, or in this case no intron,
					                     //after the last one, exists
				}
				else {
					//skip_start depending on strand and direction
					if ((*curr_introns->begin()).strand == "+") {
						skip_start =
							( dir == 1) ? (*it).start - W : (*it).end + W;
					}
					else {
						skip_start =
							( dir == 1) ? (*it).end - W : (*it).start + W;
					}

					//because W could be so big, that before skip_end, which is intron end + W*dir,
					//another intron starts
					if (skip_start*dir > skip_end*dir) {
						i = skip_end + dir;
						intron_found = true;
					}

					//skip_end depending on strand and direction
					if ( (*curr_introns->begin()).strand == "+" ) {
						skip_end =
							( dir == 1) ? (*it).end + W : (*it).start - W;
					}
					else {
						skip_end =
							( dir == 1) ? (*it).start + W : (*it).end - W;
					}
				}
			}
		}
		if ( i*dir <= (int)(max_pos)*dir ) { //otherwise out of bounds
			m_transform_vec[j] = i;
			i = i + dir; //next scaffold position
			m_last_vec_pos = j; //new last_vec_pos
			++j;
		}
		else {
			j = limit; //no more computing, break would work too
		}
	}
}



unsigned Coord_Transform::pos_transform (int pos) {

	if ( (pos <= (int)m_last_vec_pos) && (pos >= (int)(-m_win_size)) ) {
		return m_transform_vec[pos];
	}
	else {
		cerr << "Error (in Coord_Transform.cpp): Position in relative coordinate system does not exist!"
			 <<	"Please report this bug to the responsible administrator!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		//Should never happen when the tool works as it should, only possible with bugs in the tool
		abort(); 	//because this should not happen and it never happened, I excluded this error from the
					//test cases. Otherwise make a test case with the used data and use throw Error() instead
					//of abort (or better: fix the bug.)
	}
}

