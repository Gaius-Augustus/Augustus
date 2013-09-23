/**********************************************************************
 * file:    genomicMSA.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Alexander Gebauer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 04.04.12| Alexander Gebauer  | creation of the file
 * 01.07.13| Mario Stanke       | overhaul of readAlignment
 * 27.08.13| Mario Stanke       | summarization of aligmnments into gene ranges
 **********************************************************************/

#include "exoncand.hh"
#include "genomicMSA.hh"
#include "orthoexon.hh"
#include "phylotree.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <cfloat>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/exception.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>

using namespace boost;


void eraseListRange(list<int> L, list<int>::reverse_iterator from, list<int>::reverse_iterator to){
    L.erase(to.base(), from.base());
}

void eraseListRange(list<int> L, list<int>::iterator from, list<int>::iterator to){
    L.erase(from, to);
}


int GenomicMSA::maxIntronLen = 100000;
int GenomicMSA::minGeneLen = 2000; // should be at least 1, otherwise empty alignments are not filtered out
    
ostream& operator<< (ostream& strm, const MsaSignature &sig){
    strm << sig.numAli << "\t" << sig.sumAliLen << "\t" << sig.depth << setw(15);
    if (sig.color < NUMCOLNAMES)
	cout << colornames[NUMCOLNAMES - 1 - sig.color];
    else 
	cout << sig.color;
    cout << setw(MsaSignature::maxSigStrLen+1) << sig.sigstr();// << "\t";
    //    for (list<int>::const_iterator it = sig.nodes.begin(); it != sig.nodes.end(); ++it)
    //	cout << *it << " ";
    return strm;
}

bool cmpSigPtr(const MsaSignature *s1, const MsaSignature *s2){
    return *s1 < *s2;
}

/*  *.maf file is read and saved into this->alignment
 *  only the species with names in the given list are considered, the rest is ignored
 */
void GenomicMSA::readAlignment(string alignFilename) {
    int index = 0;
    string rowseq, buffer;
    string completeName;
    char strandChar;
    AlignmentRow *row;
    Alignment *alignBlock;
    map<string, size_t> notExistingSpecies;
    int lenOfChr;
    string speciesName;
    string seqID;
    int chrStart;
    int seqLen; 
    Strand strand;

    numSpecies = rsa->getNumSpecies();
    ifstream Alignmentfile;
    Alignmentfile.open(alignFilename.c_str(), ifstream::in);
    if (!Alignmentfile) {
        cerr << "Could not find the alignment file " << alignFilename << "." << endl;
        throw PropertiesError( "GenomicMSA::readAlignment: Could not open this file!" );
    }

    while (!Alignmentfile.eof()) {
        Alignmentfile >> buffer;
        int numSpeciesFound = 0;
        if (buffer == "s") {
	    alignBlock = new Alignment(numSpecies); // create new empty alignment block
            // loop over the lines of the alignment block, don't search any further if enough species were found
            while (numSpeciesFound < numSpecies && buffer == "s") {
		// reads the name of the species and the seqID from the first column
		// TODO: may have to be checked/adjusted for general .maf files
		Alignmentfile >> completeName >> chrStart >> seqLen >> strandChar >> lenOfChr >> rowseq >> buffer;
		
		// split species name and sequence ID
		for (int i=0; i<completeName.length(); i++) {
		    // seperator is the point '.' for example hs19.chr21, has to be changed
		    if ((completeName[i] == '-') || (completeName[i] == '.')) { 
			speciesName = completeName.substr(0,i);
			seqID = completeName.substr(i+1, string::npos);
			// some input file have a suffix "(..)" that needs to be stripped
			string::size_type p = seqID.find_first_of("(");
			if (p != std::string::npos)
			    seqID = seqID.erase(p); 
			break;
		    }
		    if (i == completeName.length()-1) {
			speciesName = completeName;
			seqID = "unknown";
		    }
		}
		  
		if (strandChar == '+') {
		    strand = plusstrand;
		} else if (strandChar == '-') {
		    strand = minusstrand;
		} else {
		    strand = STRAND_UNKNOWN;
		}
		if (!alignBlock->aliLen)
		    alignBlock->aliLen = rowseq.length();
		else if (alignBlock->aliLen != rowseq.length()) {
		    throw ProjectError("Error in MAF in sequence " + seqID + " at position " + itoa(chrStart) 
				       + ". Alignment row does not agree in length.");
		}
		
		row = new AlignmentRow (seqID, chrStart, strand, rowseq);
		if (seqLen != row->getSeqLen())
		    cerr << "Inconsistenty in .maf file: Sequence length different from number of non-gap characters in row:" 
			 << endl << "speciesName" << "." << seqID << "\t" << chrStart << "\t" << seqLen << "\t" << rowseq << endl;
		
		index = rsa->getIdx(speciesName);
		if (index >=0) { // species name in the white list
		    alignBlock->rows[index] = row; // place at the right position
		    // store chrLen and check whether consistent with previous chrLen
		    try{
			rsa->setLength(index, row->seqID, lenOfChr);
		    } catch (ProjectError e){
			cerr << e.getMessage() << endl << "MAF file inconsistent." << endl;
			throw e;
		    }
		    numSpeciesFound++;
		} else {
		    notExistingSpecies.insert(pair<string, size_t>(speciesName, 1));
		}
	    }
	    alignment.push_back(alignBlock);
	}
    }
    // clean up
    Alignmentfile.close();
    for (map<string,size_t>::iterator it = notExistingSpecies.begin(); it != notExistingSpecies.end(); ++it)
	cerr << "Warning: Species " << it->first << " is not included in the target list of species. These alignment lines are ingored." << endl;
}

/* printMAF
 * print alignment, to stdout if outFname is empty string
 */
void GenomicMSA::printAlignment(string outFname){
    Alignment *aliblock;
    for (list<Alignment*>::iterator alit = alignment.begin(); alit != alignment.end(); alit++) {
	aliblock = *alit;
	cout << *aliblock << endl << endl;
	/*for(vector<AlignmentRow*>::iterator aseqit = aliblock->rows.begin(); aseqit != aliblock->rows.end(); aseqit++){
	    AlignmentRow* row = *aseqit;
	    if (row) {
		cout << setw(15) << row->seqID
		     << "\tstart=" << row->start
		     << "\tseqLen=" << row->seqLen
		     << "\t" << ((row->strand == plusstrand)? "+":"_");
		cout << "\tstarts:";
		for (vector<int *>::iterator startsit = row->cmpStarts.begin(); startsit != row->cmpStarts.end(); ++startsit){
		    cout << **startsit;
		    if (startsit+1 != row->cmpStarts.end())
			cout << ",";
		}
		cout << "\tsequence:";
		for (list<block>::iterator bit = row->sequence.begin(); bit != row->sequence.end();){
		    cout << "(" << bit->begin << "," << bit->length << "," << bit->previousGaps << ")";
		    if (++bit != row->sequence.end())
			cout << ";";
		}
	    } else {
		cout << "no alignment";
	    }
	    cout << endl;
	}
	cout << endl;*/
    }
}

/**
 * Merges pairs of alignments in trivial cases in order to reduce the alignment number without doing any potentially false mergers.
 * The priority is on avoiding mergers, where there is more than one possibility to merge.
 * Not all trivial mergers are guaranteed to be realized, though.
 * Trivial mergers are found by sorting with respect to a number of (unsystematically chosen) different species.
 */
void GenomicMSA::compactify(){
    int numSortRefs = 3; // Sort wrt to at most this many species, a larger value requires more time now
                         // but may save time later by reducing alignments. For genomic alignments, where the first species is 
                         // present in ALL alignments (e.g. UCSC vertebrate), this only needs to be 1
    int numAlis = alignment.size();
    if (numSortRefs > numSpecies)
	numSortRefs = numSpecies;
    for (int i=0; i< numSortRefs; i++){
	// use the first species (probably often important to users)
	// and otherwise evenly spaced indices to the list of species
	size_t s = i * numSpecies / numSortRefs;
	alignment.sort(SortCriterion(s)); // sort by species number s
	list<Alignment*>::iterator ait, bit;
	for (ait = alignment.begin(); ait != alignment.end();) {
	    bit = ait;
	    ++bit;
	    if (bit == alignment.end())
		break;
	    // ait and bit point to neighbors in this order
	    // all rows present in both alignments must be very close neighbors (at most 3 codons distance)
	    if (mergeable(*ait, *bit, 9, 1.0)){ 
		(*ait)->merge(*bit);
		alignment.erase(bit);
	    } else {
		++ait;
	    }
	}
    }
    cout << "MAF compatification reduced the number of aligments from " << numAlis << " to " <<  alignment.size() << " (to "
	 << setprecision(3) <<  (100.0 * alignment.size() / numAlis) << "%)." << endl;
}

/**
 * 
 */
void GenomicMSA::findGeneRanges(){
    AlignmentGraph::out_edge_iterator ei, ei_end, uoi, uoi_end, voi, voi_end;
    AlignmentGraph::in_edge_iterator vi, vi_end, ui, ui_end;
    AlignmentGraph::edge_iterator fi, fi_end;
    vertex_descriptor u, v, w;
    edge_descriptor e;
    Alignment *ua, *va;
    int uid, vid, wid;
    //int maxGeneLen = 1000000 - 2*GeneMSA::padding;
    int numAlis = alignment.size();
    alignment.sort(SortCriterion(0)); // sort by first species
    vector<Alignment*> nodesA(alignment.begin(), alignment.end()); // make vector copy of alignment list for random access
    list<int> nodesI; // indices to nodesA
    AlignmentGraph aliG(numAlis);
    for (int i=0; i < numAlis; i++){
	nodesI.push_back(i); // starting order as in nodesA: sorted by first species
	v = vertex(i, aliG);
	aliG[v].a = alignment.front();
	alignment.pop_front();
    }
    property_map<AlignmentGraph, vertex_index_t>::type index = get(vertex_index, aliG);

    int numSortRefs = numSpecies; // sort wrt to at most this many species

    for (int i=0; i< numSortRefs; i++){
	size_t s = i * numSpecies / numSortRefs;
	if (i>0) // no need to sort again by 1st species
	    nodesI.sort(IdxSortCriterion(nodesA, s)); // sort indices by species number s
	list<int>::iterator ait, bit;
	for (ait = nodesI.begin(); ait != nodesI.end(); ++ait) {
	    bit = ait;
	    ++bit;
	    if (bit == nodesI.end())
		break;
	    uid = *ait;
	    vid = *bit;
	    ua = aliG[vertex(uid, aliG)].a;
	    va = aliG[vertex(vid, aliG)].a;
	    // u and v index neighbors in this order
	    if (ua && ua->rows[s] && va && va->rows[s] && mergeable(ua, va, maxIntronLen, 0.6)){
		// add edge if not already present
		add_edge(uid, vid, aliG);
	    }
	}
    }


    /* a__           a      
     *  \ \           \        implode any edges u->v, where 
     * c-u--v-w  ==>  c-u+v-w  1) alignments u and v are mergeable with mergeableFrac 1.0
     *    \  \            \       and any path through u or v may as well use edge u->v:
     *     ---d            d   2) predecessors(v) \setminus {u} \subset predecessors(u)
     *                         3) successors(u) \setminus {v} \subset successors(v)
     */
    //    AlignmentGraph::edge_iterator e;
    int numNodes = numAlis, numNodesOld;
    int itnr = 0;
    // iterate over all nodes u
    do {
	numNodesOld = numNodes;
	cout << "number of nodes: " << numNodes << endl;
	writeDot(aliG, "aliGraph." + itoa(itnr++) + ".dot");

	for (uid = 0; uid < num_vertices(aliG); uid++){
	    u = vertex(uid, aliG);
	    ua = aliG[u].a;
	    if (!ua)
		continue;
	    bool merged = true;
	    while (out_degree(u, aliG) >0 && merged){
		//cout << "checking " << uid << endl;
		merged = false;
		for (tie(ei, ei_end) = out_edges(u, aliG); ei != ei_end && !merged; ++ei){
		    v = target(*ei, aliG);
		    vid = index[v];
		    va = aliG[v].a;
		    //cout << " v= " << vid << " ";
		    if (va && mergeable(ua, va, maxIntronLen, 1.0, false)){
			// cout << "mergeable" << endl;
			// check whether predecessors(v) \setminus {u} \subset predecessors(u)
			tie(vi, vi_end) = in_edges(v, aliG);
			tie(ui, ui_end) = in_edges(u, aliG);
			bool predIncluded = true;
			while (vi != vi_end && ui != ui_end && predIncluded){
			    if (source(*ui, aliG) < source(*vi, aliG))
				++ui;
			    else if (source(*vi, aliG) == u)
				++ui;
			    else if (source(*ui, aliG) == source(*vi, aliG)){
				++ui;
				++vi;
			    } else
				predIncluded = false;
			}
			// at most the predecessor u of v may be left in the remaining predecessors(v) list
			predIncluded &= (vi == vi_end || (source(*vi, aliG) == u && ++vi == vi_end));
			//cout << "predIncluded = " << predIncluded;
			if (predIncluded){
			    // check whether successors(u) \setminus {v} \subset successors(v)
			    tie(voi, voi_end) = out_edges(v, aliG);
			    tie(uoi, uoi_end) = out_edges(u, aliG);
			    bool succIncluded = true;
			    while (voi != voi_end && uoi != uoi_end && succIncluded){
				if (target(*voi, aliG) < target(*uoi, aliG))
				    ++voi;
				else if (target(*uoi, aliG) == v)
				    ++uoi;
				else if (target(*uoi, aliG) == target(*voi, aliG)){
				    ++uoi;
				    ++voi;
				} else
				    succIncluded = false;
			    }
			    // at most the successor v of u may be left in the remaining successors(u) list
			    predIncluded &= (uoi == uoi_end || (target(*uoi, aliG) == v && ++uoi == uoi_end));
			    //cout << " succIncluded = " << succIncluded;
			    if (succIncluded){
				// append alignment v to u
				ua->merge(va);
				// add edges u->w for all edges v->w
				for (tie(voi, voi_end) = out_edges(v, aliG); voi != voi_end; ++voi){
				    w = target(*voi, aliG);
				    wid =  index[w];
				    add_edge(uid, wid, aliG);
				}
				// delete v from graph
				clear_vertex(v, aliG);
				aliG[v].a = NULL; // deactivate node
				//remove_vertex(v, aliG);
				merged = true;
			    }
			}
		    }
		    //cout << endl;
		}
	    }
	}
	numNodes = 0;
	for (uid = 0; uid <  num_vertices(aliG); uid++)
	    if (aliG[vertex(uid, aliG)].a)
		numNodes++;
    } while (numNodes < numNodesOld);

    // take all singleton alignments that are now long enough
    alignment.clear();
    int numNodes2 = 0;
    for (uid = 0; uid < num_vertices(aliG); uid++){
	Alignment* a = aliG[vertex(uid, aliG)].a;
	if (a){
	    if (out_degree(uid, aliG) == 0 && in_degree(uid, aliG) == 0){
		if (a->aliLen >= minGeneLen) // discard alignments that are too short to hold at least a short gene
		    alignment.push_back(a);
		aliG[vertex(uid, aliG)].a = NULL; // deactivate node
	    } else
		numNodes2++;
	}
    }
    
    // make a new graph, only of the active nodes with alignment (consequence of vecS as node container class)
    // store all signatures
    AlignmentGraph aliG2(numNodes2+2); // node 0 is source, node 2 is sink
    add_edge(0, 1, aliG2); // edge from source to sink directly to allow finding the empty path
    map<string, MsaSignature>::iterator sit;
    aliG2[0].a = aliG2[1].a = NULL; // no alignment for source and sink
    aliG2[0].covered = aliG2[1].covered = true; // no need to cover source and sink
    numNodes2 = 2; // node numers 0 and 1 are reserved
    bool esucc;
    for (uid = 0; uid < num_vertices(aliG); uid++){
	Alignment* a = aliG[vertex(uid, aliG)].a;
	if (a){
	    v = vertex(uid, aliG);
	    aliG[v].id = numNodes2; // id node attributes in old graph aliG are now the node indices in the new graph aliG2 and vice versa
	    aliG2[numNodes2].a = a;
	    aliG2[numNodes2].id = uid;
	    aliG2[numNodes2].covered = false;
	    // add signature
	    string sigstr = a->getSignature();
	    MsaSignature &sig = signatures[sigstr];
	    if (sig.numAli == 0){
		sig.depth = a->numFilledRows();
		sig.sigrows.resize(a->numRows(), "");
		for (size_t s = 0; s < sig.sigrows.size(); ++s)
		    if (a->rows[s])
			sig.sigrows[s] = a->rows[s]->seqID + strandChar(a->rows[s]->strand);
		if (sigstr.length() > MsaSignature::maxSigStrLen)
		    MsaSignature::maxSigStrLen = sigstr.length();
	    }
	    sig.numAli += 1;
	    sig.sumAliLen += a->aliLen;
	    //sig.nodes.push_back(uid);
	    tie(e, esucc) = add_edge(0, numNodes2, aliG2); // edge from source to any other node
	    aliG2[e].weight = 0;
	    tie(e, esucc) = add_edge(numNodes2, 1, aliG2); // edge from any node to sink
	    aliG2[e].weight = 0;
	    numNodes2++;
	}
    }
    // add edges to new graph aliG2
    for (tie(fi, fi_end) = edges(aliG); fi != fi_end; ++fi){
	int u = index[source(*fi, aliG)];
	int v = index[target(*fi, aliG)];
	add_edge(aliG[u].id, aliG[v].id, aliG2);
    }

    // sort topologically to detect cycles
    aliG2[graph_bundle].topo.resize(numNodes2);
    try {
	vector<AlignmentGraph::vertex_descriptor > topo;
	topological_sort(aliG2, std::back_inserter(topo));
	int i=0;
	for ( vector<AlignmentGraph::vertex_descriptor >::reverse_iterator ii=topo.rbegin(); ii!=topo.rend(); ++ii){
	    aliG2[graph_bundle].topo[i++] = index[*ii];
	}
    } catch (bad_graph &e) {
	// TODO, remove back edges in DFS and retry until it is a DAG
	cerr << e.what() << endl;
    }
    cout << "A topological ordering: ";
    for (int i=0; i<numNodes2; i++)
	cout << aliG2[graph_bundle].topo[i] << " ";
    cout << endl;

    cout << "number of signatures=" << signatures.size() << endl;
    list<MsaSignature*> siglist;
    for (sit = signatures.begin(); sit != signatures.end(); ++sit)
	siglist.push_back(&sit->second);
    
    siglist.sort(cmpSigPtr);
    int color = 0;
    for (list<MsaSignature*>::iterator it = siglist.begin(); it != siglist.end(); ++it){
	(*it)->color = color++;
	cout << **it << endl;
    }
    vector<AliPath> allPaths;
    writeDot(aliG2, "aliGraph." + itoa(itnr++) + ".dot");
    for (list<MsaSignature*>::iterator sit = siglist.begin(); sit != siglist.end(); ++sit){
	project(aliG2, *sit); // set weights wrt to signature
	cout << (*sit)->sigstr() << " writing aliGraph." + itoa(itnr) + ".dot" << endl;
	writeDot(aliG2, "aliGraph." + itoa(itnr++) + ".dot", *sit);
	int numNewCovered = 1;
	while (numNewCovered > 0){
	    AliPath path = getBestConsensus(aliG2, *sit, numNewCovered);
	    // determine additional value (e.g. weight of newly covered nodes
	    if (numNewCovered > 0){ // and if additional value large enough
		allPaths.push_back(path);
	    }
	}
    }
    cout << "found " << allPaths.size() << " paths" << endl;
    /*
    for (int i=0; i<allPaths.size(); i++)
	cout << allPaths[i] << endl;
    */

    // repeat pruning until nothing changes (usually changes only in the first iteration)
    int totalNumNodes;
    while (prunePaths(allPaths, aliG2)){
	//cout << "allPaths after pruning:" << endl;
	totalNumNodes = 0;
	for (int i=0; i < allPaths.size(); i++){
	    totalNumNodes += allPaths[i].path.size();
	    //cout << allPaths[i] << endl;
	}
	cout << "aliG2 nodes " << num_vertices(aliG2) << " allPaths nodes " << totalNumNodes << " ratio: " << (float) totalNumNodes/num_vertices(aliG2) << endl;
    }

    // for each path, make a single alignment and add to alignment list
    for (int i=0; i < allPaths.size(); i++){
	// cout << "alignment from path " << allPaths[i] << endl;
	list<Alignment* > plist;
	list<int> &p = allPaths[i].path;
	for (list<int>::iterator it = p.begin(); it != p.end(); ++it)
	    plist.push_back(aliG2[*it].a);
	va = mergeAliList(plist, allPaths[i].sig);
	if (va && va->aliLen >= minGeneLen) // discard alignments that are too short to hold at least a short gene
	    alignment.push_back(va);
    }
    alignment.sort(SortCriterion(0)); // sort (arbitrarily) by first species
    cout << "findGeneRanges reduced the number of aligments from " << numAlis << " to " <<  alignment.size() << " (to "
	 << setprecision(3) <<  (100.0 * alignment.size() / numAlis) << "%)." << endl;
    // delete all alignments from the graph nodes
    for (uid = 0; uid < num_vertices(aliG2); uid++){
	Alignment* a = aliG2[uid].a;
	if (a)
	    delete a;
    }
}

bool GenomicMSA::prunePaths(vector<AliPath> &allPaths, AlignmentGraph &g){
    bool changed = false;
    int i, j, m = allPaths.size();

    // this could be more efficient, instead of a running time quadratic in the number of paths
    // one could use a hash to get from each node to a list of paths that contain the node
    for (i=0; i<m; i++){
	for (j=i+1; j<m; j++){
	    changed |= prunePathWrt2Other(allPaths[j], allPaths[j].path.begin(), allPaths[j].path.end(),
					  allPaths[i], allPaths[i].path.begin(), allPaths[i].path.end(), g);
	    changed |= prunePathWrt2Other(allPaths[j], allPaths[j].path.rbegin(), allPaths[j].path.rend(),
					  allPaths[i], allPaths[i].path.rbegin(), allPaths[i].path.rend(), g);
	    changed |= deletePathWrt2Other(allPaths[j], allPaths[i], g);
	    changed |= prunePathWrt2Other(allPaths[i], allPaths[i].path.begin(), allPaths[i].path.end(),
					  allPaths[j], allPaths[j].path.begin(), allPaths[j].path.end(), g);
	    changed |= prunePathWrt2Other(allPaths[i], allPaths[i].path.rbegin(), allPaths[i].path.rend(),
					  allPaths[j], allPaths[j].path.rbegin(), allPaths[j].path.rend(), g);
	    changed |= deletePathWrt2Other(allPaths[i], allPaths[j], g);
	}
    }
    return changed;
}

/*
 * Remove from the left end of p the longest contiguous sequence of alignments that is contiguously included in 'other'
 * if node weight wrt to the respective signature is at least as large in 'other'.
 * p:            a4 a2 a9 a3 a5 a6 a0
 * other:        a1 a4 a2 a9 a7 a8 a0
 * afterwards p: a3 a5 a6 a0
 * TODO: remove p competely if the additionally aligned sequence ranges are not negligible wrt to the alignment length
 * p may become the empty path. Iterators can be reverse_iterators in which case the fragment is removed from the right end.
 */
template< class Iterator >
bool GenomicMSA::prunePathWrt2Other(AliPath &p, Iterator pstart, Iterator pend, 
				    AliPath &other, Iterator ostart, Iterator oend,
				    AlignmentGraph &g){
    //cout << "prunePathWrt2Other(" << p << endl << other << ")" << endl;
    
    // match fron the left end of p
    Iterator pa, oa;
    pa = pstart;
    oa = ostart;
    // find start of match of first node of p in other
    while (oa != oend && (*pa != *oa))
	++oa;
    if (oa != oend){ // match found
	while (oa != oend && pa != pend
	       && *pa == *oa 
	       && weight(g[*oa].a, other.sig) >=  weight(g[*pa].a, p.sig)){
	    ++oa;
	    ++pa;
	}
	// prune path until pa
	// a (future) alternative is introduce two iterators as members in AliPath:
	// the interval [first, last) that is to be kept. Then a piece of the neighboring
	// alignments can still be used to generate a better "overhang".
	if (pstart != pa){
	    eraseListRange(p.path, pstart, pa);
	    //	    cout << "pruned to " << p << endl;
	    return true;
	}
    }
    return false;
}

/*
 * Remove path p competely if the additionally aligned sequence ranges are not negligible wrt to the alignment length
 * p may become the empty path. Iterators can be reverse_iterators in which case the fragment is removed from the right end.
 */
bool GenomicMSA::deletePathWrt2Other(AliPath &p, AliPath &other, AlignmentGraph &g){
    bool superfluous = false;
    set<string> ranges[3]; // 0:p 1:other 2:intersection
    AliPath *ap[2];
    ap[0] = &p;
    ap[1] = &other;
    const MsaSignature *sigs[2];
    sigs[0] = p.sig;
    sigs[1] = other.sig;
    Alignment *a;
    for (int i=0; i<2; i++){
	for (list<int>::iterator it = ap[i]->path.begin(); it != ap[i]->path.end(); ++it){
	    for (int s=0; s<numSpecies; s++){
		a = g[*it].a;
		if (a->rows[s] && sigs[i]->fits(*a, s)){
		    string key = itoa(s) + ":" + a->rows[s]->seqID + itoa(a->rows[s]->strand) + " "
			+ itoa(a->rows[s]->chrStart()) + "-" + itoa(a->rows[s]->chrEnd());
		    ranges[i].insert(key);
		}
	    }
	}
    }
    set_intersection(ranges[0].begin(), ranges[0].end(), ranges[1].begin(), ranges[1].end(), std::inserter( ranges[2], ranges[2].begin()));
    //    for (int i=0; i<3; i++)
    //	cout << "set " << i << " has size " << ranges[i].size() << endl;
    // check whether p has too few additional aligned regions that are not shared by 'other'
    if (ranges[0].size() < 1.1 * ranges[2].size()){
	superfluous = true;
	// cout << p << endl << "is superfluous because of" << endl << other << endl;
    }
    if (superfluous)
	p.path.clear(); // delete all alignments from path, will be discarded later
    
    return superfluous;
}

    


AliPath GenomicMSA::getBestConsensus(AlignmentGraph &g, const MsaSignature *sig, int &numNewCovered){
    AliPath p;
    Alignment *a;
    int N = num_vertices(g);
    vector<int> pred(N); // predecessor for each node on shortest path
    int vid;
    int maxWeight = findBestPath(g);
    numNewCovered = 0;
    vid = g[1].pred; // last node on best path
    string newlyids = "";
    while(vid != 0){ // until source node is reached
	p.path.push_front(vid);
	a = g[vid].a;
	if (!g[vid].covered && a->numFitting(sig) == a->numFilledRows()){
	    g[vid].covered = true; // path covers all alignments with the exact signature or a subset thereof
	    g[vid].weight = weight(g[vid].a, sig);
	    numNewCovered++;
	    newlyids = itoa(vid) + " " + newlyids;
	}
	vid = g[vid].pred;
    }
    if (numNewCovered>0 && false)
	cout << "newly covered " << newlyids << "\t" << numNewCovered << " new nodes. With exact signature: " << maxWeight / g[graph_bundle].maxWexceptCov << endl;
    p.sig = sig;
    return p;
}

int GenomicMSA::findBestPath(AlignmentGraph &g){
    // relax all nodes in topologial order
    int N = num_vertices(g);
    int i, vid, wid;	
    AlignmentGraph::out_edge_iterator ei, ei_end;
    property_map<AlignmentGraph, vertex_index_t>::type index = get(vertex_index, g);
    vertex_descriptor v, w;
    vector<int> maxCumW(N, -INT_MAX); // currently maximal weight of any path from source to each node
    if (g[graph_bundle].topo.size() != N)
	throw ProjectError("Topological order required for findBestPath.");
    // initialize source node
    maxCumW[0] = 0;
    g[0].pred = 0; // source node has itself as predecessor
    for (i=0; i<N; i++){
	vid = g[graph_bundle].topo[i];
	v = vertex(vid, g);
	// relax all outgoing edges from v
	for (tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei){
	    w = target(*ei, g);
	    wid = index[w];
	    if (g[*ei].weight > -INT_MAX && g[wid].weight > -INT_MAX &&
		maxCumW[vid] + g[*ei].weight + g[wid].weight > maxCumW[wid]){
		//cout << "relaxed " << vid << " -> " << wid << endl;
		maxCumW[wid] = maxCumW[vid] + g[*ei].weight + g[wid].weight;
		g[wid].pred = v;
	    }
	}
    }
    //    for (vid=0; vid<N; vid++)
    //	cout << vid << "\t" << g[vid].pred << " " << maxCumW[vid] << endl;

    return  maxCumW[1];
}

void GenomicMSA::project(AlignmentGraph &g, const MsaSignature *sig){
    Alignment *a, *b;
    vertex_descriptor u, v;
    int maxWexceptCov = 0;
    // set node weights
    for (int vid = 0; vid < num_vertices(g); vid++){
	v = vertex(vid, g);
	a = g[v].a;
	g[v].weight = weight(a, sig);
	if (g[v].weight > 0)
	    maxWexceptCov += g[v].weight;
    } 
    AlignmentGraph::edge_iterator ei, ei_end;
    // set edge weights
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
	u = source(*ei, g);
	v = target(*ei, g);
	a = g[u].a;
	b = g[v].a;
	g[*ei].weight = weight(a, b, sig);
	if (g[*ei].weight > 0)
	    maxWexceptCov += g[*ei].weight;
    }
    g[graph_bundle].maxWexceptCov = maxWexceptCov;
    //cout << "maxWexceptCov = " << maxWexceptCov << endl;
    // increase node weights for not yet covered nodes of the given signature, so they are first priority
    for (int vid = 0; vid < num_vertices(g); vid++){
	v = vertex(vid, g);
	if (!g[v].covered && sig->sigstr() == g[vid].a->getSignature())
	    g[v].weight += maxWexceptCov + 1;
    }
}

/*
 * Node weight, when projecting a to sig:
 * -Infty, if less than 2 rows agree with signature, otherwise it is the cumulative Fragment length of the fitting rows
 */
int GenomicMSA::weight(const Alignment *a, const MsaSignature *sig){
    int cumFragLen = 0;
    if (!a)
	return 0;
    if (a->numFitting(sig) < 2)
	return -INT_MAX;
    for (size_t s = 0; s < a->numRows(); s++){
	if (sig->fits(*a, s))
	    cumFragLen += a->rows[s]->getCumFragLen();
    }
    return cumFragLen;
}

/*
 * Edge weight, when projecting a to sig:
 * -Infty, if the projected alignments are not mergeable.
 */
int GenomicMSA::weight(const Alignment *a, const Alignment *b, const MsaSignature *sig){
    int numFitting = 0;
    if (!a || !b)
	return 0;
    for (size_t s = 0; s < a->numRows(); s++){
	if (sig->fits(*a, s) && sig->fits(*b, s)){
	    int dist = b->rows[s]->chrStart() - a->rows[s]->chrEnd() - 1;
	    if (dist >=0 && dist <= maxIntronLen)
		numFitting++;
	}
    }
    if (numFitting < 2)
	return -INT_MAX;
    return 0;
}

void GenomicMSA::writeDot(AlignmentGraph const &g, string fname, MsaSignature const *superSig){
    property_map<AlignmentGraph, vertex_index_t>::type index = get(vertex_index, g);
    ofstream dot(fname);
    dot << "digraph G {" << endl;
    dot << "node[shape=box];" << endl;
    // legend box with species names
    dot << "legend1[label=<<TABLE BORDER=\"0\"><TR><TD COLSPAN=\"2\" BGCOLOR=\"deepskyblue1\"><B>Species</B></TD></TR>" << endl;
    for (size_t s=0; s<numSpecies; s++){
	dot << "<TR><TD align=\"left\"";
	if (s < NUMCOLNAMES)
	    dot << " BGCOLOR=\"" << colornames[s] << "\"";
	dot << ">" << s << "</TD><TD align=\"right\">" << rsa->getSname(s) << "</TD></TR>" << endl;
    }
    dot << "</TABLE>>];" << endl;
    // legend box for alignments
    dot << "legend2[label=<<TABLE BORDER=\"0\"><TR><TD COLSPAN=\"5\" BGCOLOR=\"deepskyblue1\"><B>Alignment</B></TD></TR>";
    dot << "<TR><TD COLSPAN=\"2\">id</TD><TD COLSPAN=\"3\">alignment length</TD></TR>";
    dot << "<TR><TD align=\"left\" BGCOLOR=\"" << colornames[0] << "\">0</TD><TD align=\"right\">chrN</TD>"
	<< "<TD>strand</TD><TD>genomeStartPos</TD><TD>genomeEndPos+1</TD></TR>";
    if (superSig)
	dot << "<TR><TD COLSPAN=\"3\">weight</TD><TD COLSPAN=\"2\">covered</TD></TR>";
    dot << "</TABLE>>];" << endl;

    dot << "rankdir=LR;" << endl;
    
    graph_traits<AlignmentGraph>::vertex_iterator vi, vi_end;
    AlignmentRow *row;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
	int i =  index[*vi]; // index of vertex
	Alignment *ia = g[*vi].a;
	if (ia){
	    MsaSignature sig;
	    string sigbgcol = "";
	    if (signatures.size()>0){
		sig = signatures[ia->getSignature()]; // signature of this node
		if (sig.color < NUMCOLNAMES)
		    sigbgcol = " BGCOLOR=\"" + colornames[NUMCOLNAMES - 1 - sig.color] + "\"";
	    }
	    dot << i <<  "[label=<<TABLE BORDER=\"0\">";
	    dot << "<TR><TD COLSPAN=\"2\">" << i << "=" << g[*vi].id << "</TD><TD COLSPAN=\"3\">" << ia->aliLen << "</TD></TR>";
	    for (size_t s = 0; s < ia->rows.size(); s++){
		dot << "<TR>";
		row = ia->rows[s];
		if (row) {
		    dot << "<TD align=\"right\"";
		    if (s < NUMCOLNAMES){
			dot << " BGCOLOR=\"" << colornames[s] << "\"";
		    }
		    dot << ">" << s << "</TD>";
		    dot << "<TD align=\"left\"" << ((!superSig || superSig->fits(*ia, s))? sigbgcol : "" ) <<  ">" << row->seqID << "</TD>"; 
		    dot << "<TD" << ((!superSig || superSig->fits(*ia, s))? sigbgcol : "" ) << ">" << row->strand << "</TD>";
		    dot << "<TD align=\"right\">" << row->chrStart() << "</TD>";
		    dot << "<TD align=\"right\">" << row->chrEnd() << "</TD>";
		}
		dot << "<TD></TD><TD></TD><TD></TD><TD></TD><TD></TD></TR>";
	    }
	    if (superSig) {
		dot << "<TR><TD COLSPAN=\"3\">";
		if (g[*vi].weight > -INT_MAX)
		    dot << g[*vi].weight;
		else 
		    dot << "-infty";
		dot << "</TD><TD COLSPAN=\"2\">" << g[*vi].covered << "</TD></TR>";
	    }
	    dot << "</TABLE>>";
	    if (superSig){
		if (g[*vi].weight == -INT_MAX)
		    dot << ",style=dotted";
		else {
		    dot << ",color=red";
		    if (sig.sigstr() == superSig->sigstr())
			dot << ",peripheries=2,style=filled";
		}
	    }
	    dot << "];" << endl;
	}
    }

    graph_traits<AlignmentGraph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
	int u = index[source(*ei, g)];
	int v = index[target(*ei, g)];
	Alignment *ua = g[vertex(u, g)].a;
	Alignment *va = g[vertex(v, g)].a;
	if (ua && va){ // do not display edges from source and to sink node
	    // compute average gap length between the two aligments
	    int k=0;
	    int avGapLen = 0;
	    for (size_t s=0; s < numSpecies; s++)
		if (ua->rows[s] && va->rows[s] && ua->rows[s]->seqID == va->rows[s]->seqID
		    && ua->rows[s]->strand == va->rows[s]->strand &&
		    (!superSig || (superSig->fits(*ua, s) && superSig->fits(*va, s)))){
		    int gapLen = va->rows[s]->chrStart() - ua->rows[s]->chrEnd() - 1;
		    if (gapLen >= 0){
			avGapLen += gapLen;
			k++;
		    }
		}
	    if (k > 0)
		avGapLen /= k;
	    dot << u << "->" << v << "[label=" << avGapLen;
	    if (superSig){
		if (g[*ei].weight == -INT_MAX)
		    dot << ",style=dotted";
		else
		    dot << ",color=red";
	    }
	    dot << "];" << endl;
	}
    }
    dot << "}" << endl;
    dot.close();
}



// pops the first alignment from list
GeneMSA* GenomicMSA::getNextGene() {
    if (alignment.empty())
	return NULL;
    GeneMSA *geneRange = new GeneMSA(rsa, alignment.front());
    alignment.pop_front();
    return geneRange;
}

ostream& operator<< (ostream& strm, const AliPath &p){
    for(list<int>::const_iterator pit = p.path.begin(); pit != p.path.end(); ++pit)
	strm << *pit << " ";
    strm << "\t" << p.sig->sigstr();
    return strm;
}
