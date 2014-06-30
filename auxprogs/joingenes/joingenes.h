#ifndef _JOIN_GENES_HPP
#define _JOIN_GENES_HPP

#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <iostream>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <fstream>

using namespace std;

class Gene;
class Transcript;

class Exon{
	public:
	int from, to;
	string chr;
	string feature;
	float score;
	int frame;

	bool operator<(Exon const& rhs) const {
			return (from < rhs.from);
	}
};

class Transcript{
	public:
	list<Transcript*> supporter;		// at the moment these list only consist the pointer to transcripts, who are totaly equal to this transcript (doesnt matter whether their creation is related)
	list<Exon> exon_list;				// at the moment these list is only used for CDSs and thats the way it is used!
	Gene* parent;
	string source;
	string t_id;
	char strand;
	int start;							// transcript start (lowest base position in transcript) 
	int stop;							// transcript stop (highest base position in transcript)
	int start_codon;					// tis (translation initation site)
	int stop_codon;						// tes (translation end site)
	int priority;
	pair<int,int> pred_range;

	void sort_exon(){
		exon_list.sort();
	}
	int getstart(){
		return exon_list.front().from;
	}
	int getstop(){
		return exon_list.back().to;
	}
	bool operator<(Transcript const& rhs) const {
		if (start != rhs.start)
			return (start < rhs.start);
		else
			return (stop < rhs.stop);
	}
};

class Gene{
	public:
	string g_id;
	list<Transcript*> children;
};

void seek_overlaps(list<Transcript> &transcript_list, string &outfilename);
void work_at_overlap(list<Transcript*> &overlap, list<Transcript> &new_transcripts);
void search_n_destroy_doublings(list<Transcript*> &overlap);
bool compare_transcripts(Transcript const* t1, Transcript const* t2);
pair<bool,bool> is_part_of(Transcript const* t1, Transcript const* t2);
bool compare_priority(Transcript const* lhs, Transcript const* rhs);
void join_stop(list<Transcript*> &overlap, list<Transcript> &new_transcripts, list<Transcript*> &new_overlap_part);
void join_start(list<Transcript*> &overlap, list<Transcript> &new_transcripts, list<Transcript*> &new_overlap_part);
void joining(Transcript* t1, Transcript* t2, char strand, char side, Transcript &tx_new);
bool is_combinable(Transcript const* t1, Transcript const* t2, char strand, char side);
bool check_frame_annotation(Transcript const &transcript);

#endif
