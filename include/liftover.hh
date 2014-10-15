/*****************************************************************************\
 * Filename : liftover.hh
 * Authors  : Stefanie Koenig
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 03.10.2014 | Stefanie Koenig       | Creation of the file
\******************************************************************************/

#ifndef _LIFTOVER_HH
#define _LIFTOVER_HH

// project includes
#include "alignment.hh"

/*
 * SeqIntKey encodes the start position and length of a sequence interval in a 64 bit integer
 *   
 * bits           42                 15          7 
 *      |--------------------|---------------|-------|
 *         start position           length     
 *
 * 7 extra bits for additional information (optional)
 * f.e. if sequence interval is an exon, 5 extra bits for type and 2 extra bits for lenMod3
 */
class SeqIntKey{

public:
    SeqIntKey(int_fast64_t start, int_fast64_t len, int_fast64_t addInfo){
	key = (start << 22)  // 42 bits
	    + (len << 7 )    // 15 bits
	    + (addInfo);     //  7 bits
    }
    SeqIntKey(int_fast64_t k) : key(k) {}
    ~SeqIntKey() {}
    inline int_fast64_t getStart() const {return key >> 22;}
    inline int_fast64_t getLen() const {return ((key>>7) & 0x7FFF);}  // 0x7FFF bit mask for retrieving bits 1-15
    inline int_fast64_t getAddInfo() const {return key & 0x7F;}       // 0x7F bit mask for retrieving bits 1-7
    inline int_fast64_t getKey() const {return key;}

private:
    int_fast64_t key;
};

/*
 * Mapping of sequence intervals from positions in the genome
 * to positions in the alignment and vice versa.
 * A sequence interval can be any object T for which functions
 *
 * int_fast64_ T::getKey()                        // returns a SeqIntKey encoding start position and length of the SI
 * T* create(int_fast64_t key , const char* dna)  // creates a new SI from a SeqIntKey after verification (e.g. if SI is an
 *                                                // intron, check whether splice sites are present in the sequence)
 *
 * exist. Sequence positions can be encoded as zero-length SIs.
 * Purpose: lifting over of exons, introns and non-canonical SS to other species
 * (is it possible to generalize this class to lifting over whole gene structures?)
 */
class LiftOver {

private:
    Alignment* alignment;
    vector<int> &offsets;

public:
    LiftOver(Alignment* a, vector<int> &o) : alignment(a), offsets(o) {}
    ~LiftOver() {}

    /*
     * mapping sequence intervals from positions in the genome to positions in the alignment
     * input:  a vector of hashes (keys encode start and length of the SIs in the genome, value are pointers to the SIs)
     *         ( i-th position in the vector points to the hash of SIs of species i)
     * output: a hash that stores for each sequence interval T in alignment space (key encodes start and length of the SI in the alignment)
     *         a list of (speciesIdx, T*), e.g.
     *         alignedSIs["100:200"] = {(0, si0), (3, si3)}
     */
    template <typename T> void projectToAli(vector< map<int_fast64_t,T*> > &seqints, // input
					    map<int_fast64_t, list<pair<int,T*> > > &alignedSIs) { // output

	int k = alignment->rows.size();
	int_fast64_t aliStart, aliEnd;
	int chrStart, chrEnd;
	
	typename map<int_fast64_t, list<pair<int,T*> > >::iterator asi;
	/* 
	 * map all SIs from genome coordinates to alignment coordinates where possible
	 * this search in LINEAR in the length of all sequence intervals
	 * + the number of all alignment fragments
	 */
	for (size_t s=0; s<k; s++){
	    if (alignment->rows[s] == NULL)
		continue;
	    int offset = offsets[s];
	    AlignmentRow *row = alignment->rows[s];
	    vector<fragment>::const_iterator from = row->frags.begin();
	    for(typename map<int_fast64_t,T*>::iterator siit = seqints[s].begin(); siit != seqints[s].end(); ++siit){

		// go the the first fragment that may contain the SI start
		while (from != row->frags.end() && from->chrPos + from->len - 1 < siit->second->getStart() + offset)
		    ++from;
		if (from == row->frags.end())
		    break; // have searched beyond the last alignment fragment => finished 
		SeqIntKey chrKey(siit->first); // key encodes SI start position and length in the genome
		chrStart = chrKey.getStart() + offset;
		aliStart = row->getAliPos(chrStart, from);
		if (aliStart >= 0){ // left sequence boundary mappable
		    chrEnd = chrKey.getLen() + chrStart;
		    aliEnd = row->getAliPos(chrEnd, from);
		    if (aliEnd >= 0){
			// both sequence boundaries were mappable
			// store the SI in the hash
			SeqIntKey aliKey(aliStart,aliEnd-aliStart,chrKey.getAddInfo()); 
			asi = alignedSIs.find(aliKey.getKey()); // key encodes SI start position and length in the alignment
			if (asi == alignedSIs.end()){ // insert new list
			    list<pair<int,T*> > si;
			    si.push_back(pair<int,T*> (s, siit->second));
			    alignedSIs.insert(pair<int_fast64_t,list<pair<int,T*> > >(aliKey.getKey(), si));
			} else {// append new entry to existing list
			    asi->second.push_back(pair<int,T*> (s, siit->second));
			}
		    }
		}
	    }
	}
    }
    /*
     * mapping sequence intervals from positions in the alignment to positions in the genomes
     * input:  a hash that stores for each sequence interval T in alignment space a list of (speciesIdx, T*)
     *         a vector of AnnoSeqs (necessary to verify whether an si is mappable to the genome)
     *         a vector of hashes of SIs. SIs that are mappable to other species are inserted into the
     *         corresponding hash.
     *         
     */
    template <typename T> void projectToGenome(map<int_fast64_t, list<pair<int,T*> > > &alignedSIs, vector<AnnoSequence*> const &seqRanges,
					       vector< map<int_fast64_t,T*> > &seqints) {

	int k = alignment->rows.size();
	vector<vector<fragment>::const_iterator > fragsit(k);

	// move fragment iterators to the start
	for (size_t s=0; s < k; s++){
	    if (alignment->rows[s])
		fragsit[s] = alignment->rows[s]->frags.begin();
	}

	typename map<int_fast64_t, list<pair<int,T*> > >::iterator asi;
	/* 
	 * map all SIs from alignment coordinates to genome coordinates where possible
	 * this search in LINEAR in the length of all sequence intervals
	 * + the number of all alignment fragments
	 */
	for(asi = alignedSIs.begin(); asi != alignedSIs.end(); asi++){
	    SeqIntKey aliKey(asi->first); // key encodes SI start position and length in the alignment
	    int_fast64_t aliStart = aliKey.getStart();
	    int_fast64_t aliEnd = aliKey.getLen() + aliStart;

	    typename list<pair<int,T*> >::iterator siit = asi->second.begin();
	           
	    for(size_t s=0; s<k; s++){
		if(siit != asi->second.end() && siit->first == s){ // SI is already in the hash
		    ++siit;
		    continue;
		}
		if (alignment->rows[s] == NULL)
		    continue;
		AlignmentRow *row = alignment->rows[s];
		// go the the first fragment that may contain aliStart
		while (fragsit[s] != row->frags.end() && fragsit[s]->aliPos + fragsit[s]->len - 1 < aliStart){
		    ++fragsit[s];
		}
		if (fragsit[s] == row->frags.end()){
		    continue; // searched beyond the last fragment
		}
		int_fast64_t chrStart = row->getChrPos(aliStart,fragsit[s]) - offsets[s];
		if (chrStart >= 0){ // left sequence boundary mappable
		    int_fast64_t chrEnd = row->getChrPos(aliEnd, fragsit[s]) - offsets[s];
		    if (chrEnd >= 0){ // both sequence boundaries mappable
			SeqIntKey chrKey(chrStart,chrEnd-chrStart,aliKey.getAddInfo()); // key encodes SI start position and length in the genome
			T* si = create(chrKey.getKey(), seqRanges[s]->sequence, seqRanges[s]->length);
			if(si){
			    asi->second.insert(siit,pair<int,T*> (s, si));
			    // insert mapped SIs into hash of s-th species
			    seqints[s].insert(pair<int_fast64_t,T*>(chrKey.getKey(), si));
			}   
		    }
		}
	    }		
	}
    }
};

#endif   //  _LIFTOVER_HH
