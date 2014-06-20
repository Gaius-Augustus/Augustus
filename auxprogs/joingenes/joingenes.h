#ifndef _JOIN_GENES_HPP
#define _JOIN_GENES_HPP

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include <unordered_map>
#include <sstream>
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
		if (from != rhs.from)
			return (from < rhs.from);
		else
			return (to < rhs.to);
	}
};

class Transcript{
	public:
	list<Transcript*> supporter;
	list<Exon> exon_list;
	Gene* parent;
	string source;
	string t_id;
	char strand;
	string getchr();	// unused	
	bool complete;
	int start;
	int stop;
	int start_codon;
	int stop_codon;
	int priority;
	int frame;

	int getstart(){
		/*	// also gets start but dont sort list, what is needed to compare two transcripts
		int start_temp = exon_list.front().from;
		for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
			if (start_temp > (*it).from)
				start_temp = (*it).from;
		}
		return start_temp;
		*/
		exon_list.sort();			// takes longer but sorts exon_list
		return exon_list.front().from;
	}
	int getstop(){
		int stop_temp = exon_list.front().to;
		for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
			if (stop_temp < (*it).to)
				stop_temp = (*it).to;
		}
		return stop_temp;
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


void seek_overlaps(list<Transcript> &transcript_list);
void decide(list<Transcript*> &overlap);
bool compare_transcripts(Transcript* t1, Transcript* t2);
pair<bool,bool> is_part_of(Transcript* t1, Transcript* t2);
bool compare_priority(Transcript* lhs, Transcript* rhs);
void join(Transcript* t1, Transcript* t2);
bool three_compatible(Transcript* t1, Transcript* t2);
bool five_compatible(Transcript* t1, Transcript* t2);
void check_frame_annotation(list<Transcript> &transcript_list);

/*
	string seqname;		// in Exon
	string source;		// in Transcript
	string feature;		// in Exon
	unsigned int begin;	// in Exon
	unsigned int end;	// in Exon
	double score;		// in Exon
	char strand;		// in Transcript
	char frame;		// in Exon
	string attribute;	// Transcript + Gene
*/

/*definitions out off ensemble.org:
    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
*/

#endif
