/***********************************************************************
 * file:    vitmatrix.cc
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  data structures used inside the viterbi matrix
 * author:  Oliver Keller (keller@cs.uni-goettingen.de)
 *
 **********************************************************************/

#include "vitmatrix.hh"

#ifdef DEBUG
#include "pp_scoring.hh"
#endif

#include <iostream>
#include <climits>
#include <cstdlib>

/* --- ViterbiSubmapType methods ----------------------------------- */

SubstateId ViterbiSubmapType::popPredSubstate(SubstateId substate) const {
    if (!theMap)
	return SubstateId();
    SubstateIdMap*& predStateMap = theMap->predSubstates;
    if (!predStateMap)
	return substate;
    SubstateIdMap::iterator it = predStateMap->find(substate);
    if (it == predStateMap->end())
	return substate;
    SubstateId result=it->second;
    predStateMap->erase(it);
    if (predStateMap->empty()) {
	delete predStateMap; 
	predStateMap = 0;
    }
    return result;
}


/*
 * ViterbiSubmapType::merge(other, mult):
 * merge two maps and take the maximum entries, after multiplying
 * entries of other map with _mult_
 * note: dont use this on linked maps since it would merge the other
 *       map into all submaps sharing the link
 *
 */
void ViterbiSubmapType::merge(ViterbiSubmapType& other, Double mult) {
    if (other.empty())
	return;
    mult *= other.factor;
#ifdef DEBUG
    if (is_linked()) 
	// this should not happen
	throw ProjectError("tried to merge into a shared map");
    if (inactive()) 
	// this should not happen
	throw ProjectError("tried to merge into an inactive map");
#endif
    mult /=  factor;
    iterator other_it = other.begin();
    iterator it = theMap->lower_bound(other_it->first);
    for (; other_it != other.end(); ++other_it) {
	SubstateId substate=other_it->first;
	Double p = other_it->second.p * mult;
#ifdef DEBUG
	if (p * factor >1)
	    throw ProjectError("ViterbiSubmapType::merge: have prob > 1 !!! substate="+itoa(substate));
#endif
	while (true) {
	    if (it == end() || it->first > substate) {
		insert(it, substate, p, &other); 
		break;
	    }
	    if (it->first == substate) {
		if (it->second.p < p)
		    it->second.set(p, &other);
		++it; break;
	    }
	    ++it;
	}
    }
}

/*
 * eraseUnneededSubstates:
 *
 * decrement the successor counts for the given substate
 *
 *
 */
int ViterbiSubmapType::eraseUnneededSubstate(SubstateId substate) {
#ifdef DEBUG
    if (succCount(substate)==0)
	throw ProjectError("eraseUnneededSubstate: We cannot delete substate."
			   " Successor count already 0!");
#endif
    ViterbiSubmapEntry& ent = theMap->find(substate)->second;
    if (--ent.succCount > 0)
	return 0;
    ViterbiSubmapType* predecessorMap = ent.predMap;
    theMap->erase(substate);
    substate = popPredSubstate(substate);
    int result = 0;
    if (empty()) {
	// result = linkcount() -1;
	// deactivate();
    }
    if (predecessorMap)
	result += predecessorMap->eraseUnneededSubstate(substate);
    return result;
}

/*
 * registerPredecessors
 *
 * i) delete from bestLeftId, if present, all entries with keys
 *    that do not exist in the substate map
 * ii) increment the successor counts of the predecessor submap entries,
 *     for entries that have a predecessor substate
 */

void ViterbiSubmapType::registerPredecessors()  {
    SubstateIdMap* bestLeftId = theMap->predSubstates;
    if (bestLeftId == 0 || bestLeftId->empty())  {
	// if there is no predecessor map, use identical substates
	if (bestLeftId)
	    clearEmptySubstateMap();
	for (iterator it = begin(); it != end(); ++it)
	    it->second.predMap->succCount(it->first)++;
	return;
    }

    SubstateIdMap::iterator idit = bestLeftId->begin();
    for (iterator it = begin(); it != end(); ++it) {
	while (idit != bestLeftId->end() && idit->first < it->first)
	    bestLeftId->erase(idit++);
	ViterbiSubmapType* pred = it->second.predMap;
	if (pred) {
	    SubstateId& predSubstate = idit->second;
#ifdef DEBUG
	    if (idit == bestLeftId->end() || idit->first > it->first)
		throw ProjectError("registerPredecessor: Missing predecessor");
	    if (pred->get(predSubstate) <= 0) 
		throw ProjectError("registerPredecessors: Predecessor substate does not exist");
#endif
	    pred->succCount(predSubstate)++;
	}
	++idit;
    }
    bestLeftId->erase(idit, bestLeftId->end());
    
    clearEmptySubstateMap();
}

void ViterbiSubmapType::useMaximum(Double& mainProb) {
    iterator it, result=end();
    for (it = begin(); it != end(); ++it)
	if (it->second.p > mainProb) {
	    mainProb = it->second.p;
	    result = it;
	} 
    erase(begin(), result);
    if (result != end())
	erase(++result, end());
}
    
/*
 * eraseAllUnneededSubstates:
 *
 * For each entry, check if the successor counts are positive, and if not,
 * delete the state and call eraseUnnededSubstate to decrease the 
 * successor counts of the predecessor. 
 *
 */
int ViterbiSubmapType::eraseAllUnneededSubstates() const {
    if (empty())
	return 0;
    int result = 0;
    iterator it = begin(); 
    while (it != end()) {
	if (succCount(it) == 0) {
	    ViterbiSubmapType* predecessorMap = it->second.predMap;
	    if (predecessorMap) 
		result += predecessorMap->eraseUnneededSubstate(popPredSubstate(it->first));
#ifdef DEBUG
	    else if ((state>=5 && state<=8) ||
		     (state>=25 && state<=28))
		cerr << (int)state << ",";
#endif
	    theMap->erase(it++);
	    continue;
	} 
#ifdef DEBUG
	if (succCount(it) > 10000 || succCount(it) < 0) 
	    throw ProjectError("Have successor count " + itoa(succCount(it)) + ", something is wrong");
#endif
 	++it;
    }  
#ifdef DEBUG
    if (begin(MAX_BLOCKCOUNT) != end())
	throw ProjectError("Have copy-ids left: " + itoa(begin(MAX_BLOCKCOUNT)->first.slot) + "/" +  itoa(begin(MAX_BLOCKCOUNT)->first.value));
#endif
     // if (empty())
    // 	owner->eraseSubstate(state);
    if (empty()) {
	// result += linkcount() -1;
	// deactivate();
    }
   return result;
}

/*
 * addSubstates: add a new ViterbiSubmap
 *
 * The new map will be deleted if a new map is added for
 * the same state, or if the ViterbiColumn is deleted.
 *
 */
ViterbiSubmapType& ViterbiColumnType::addSubstates(const ViterbiSubmapType& newmap) {
    int state = newmap.state;
    if (has(state)) {
	signed char& subidx = elem(state).substate_idx;
	if (subidx >= 0) { // substate map exists
	    subProbs[subidx] = newmap;
	    // subProbs[subidx].owner = this;
	    return subProbs[subidx];
	}
	subidx = subProbs.size();
    } else // create state 
	addState(state, subProbs.size());
    // create substate map
    subProbs.push_back(newmap);
    // subProbs.back().owner = this;
    return subProbs.back();
}

 
/* --- ViterbiColumnType methods ----------------------------------- */

void ViterbiColumnType::reset(int statecount) {
    delete[] idx;
    maxsize = statecount;
    idx = new signed char[statecount];
    subProbs.clear();
    for (int state=0; state<statecount; state++)
	idx[state]=-1;
    clear();
}

void ViterbiColumnType::eraseSubstates(int state) {
    if (has(state)) {
	signed char& subidx = elem(state).substate_idx;
	eraseSubProbs(subidx);
	subidx = -1;
    }
}

void ViterbiColumnType::erase(int state) {
    if (has(state)) {
	eraseSubProbs(elem(state).substate_idx);
	// if state is not last in the vector,
	// move the last in the vector to its place
	// so we can perform a pop_back()
	// idx[state] is the correspondig slot
	if (idx[state] < size()-1) {
	    // set the new index for the previously last entry
	    idx[back().state] = idx[state];
	    // copy it to its new location 
	    elem(state) = back();
	}
	pop_back();
	idx[state]=-1;
    }
}

void ViterbiColumnType::addState(int state, signed char sub_idx) {
    idx[state]=size();
    int cap = capacity();
    if (cap <= size()) {
	cap = (cap*3)/2; 
	if (cap > maxsize) 
	    cap=maxsize;
	else if (cap == 0)
	    cap = VIT_MIN_CAPACITY;
	reserve(cap);
    }
    push_back(ViterbiEntryType(state, 0, sub_idx));
}


/* --- OptionsList methods ----------------------------------------- */

OptionListItem OptionsList::sample() {	if (size() == 0)
	throw ProjectError("Tried to sample from empty list.");
    if (!(cumprob > 0.0))
	throw ProjectError("Tried to sample from list of impossible options.");
    OptionListItem oli;
    Double z = Double((double) rand() / RAND_MAX) * cumprob * 0.99999; // last factor because of computer arithmetic errors
    Double cumsum = 0.0;
    list<OptionListItem>::iterator it = options.begin();
    bool found = false;
    while (!found && it != options.end()) {
	cumsum += it->probability;
	found = (z < cumsum);
	if (!found)
	    ++it;
    }
    if (!found){
	// check later on this error
	// possible solution: decrease z by a tiny percentage.
	cerr << "Sampling Error which should not happen " << cumsum << " < " << z << endl;
	// throw ProjectError("Sampling Error which should not happen");
	it = options.begin(); // TODO: throw an error again
    }
    oli = *it;
    oli.probability /= cumprob;
    return oli;
}


/* --- DEBUGGING methods -------------------------------------------- */

#ifdef DEBUG
void ViterbiSubmapType::dieOnZero() const {
    if (factor==0)
	throw ProjectError("ViterbiSubmap has zero entry. (state=" + itoa(state) + ", all substates.)");
    if (theMap == 0) return;
    for (iterator it=begin(); it!= end(); it++) 
	if (it->second.p == 0)
	    throw ProjectError("ViterbiSubmap has zero entry. (state=" + itoa(state) + ", substate=" + itoa(it->first) + ")");
}

void ViterbiSubmapType::checkAfterWindow() const {
    if (inactive())
	return;
    for (iterator it = begin();
	 it != end();
	 ++it) 
	if (succCount(it) < 0 || succCount(it) > 10000)
	    throw ProjectError("succCounts[" + itoa(succCount(it)) + "]== " +
			       itoa(succCount(it)));
	else if (it->first.slot >= MAX_BLOCKCOUNT)
	    throw ProjectError("Copy-Index not deleted!");
}

SubstateId ViterbiSubmapType::getPredecessorSubstate(SubstateId substate) const {
    if (theMap && theMap->predSubstates) {
	SubstateIdMap::const_iterator 
	    it = theMap->predSubstates->find(substate);
	if (it != theMap->predSubstates->end())
	    return (it->second);
    }
    return substate;
}

void ViterbiSubmapType::checkCounts() const {
    if (theMap)
	for (iterator it = begin(); it != end(); ++it) {
	    SubstateId predSubstate = getPredecessorSubstate(it->first);
	    if (!predSubstate.undef() && it->second.predMap && 
		it->second.predMap->succCount(predSubstate) <= 0)
		throw ProjectError("checkCounts() failed");
	}
}

void ViterbiColumnType::testIntegrity() const {
    for (int i=0; i < subProbs.size(); i++) 
	if (i != getSubstateIdx(subProbs[i].state))
	    throw ProjectError("ViterbiColumnType: internal error: integrity check failed (subProbs[" + itoa(i) + "] exists but its state should have another index)."); // enable later!!
	// else if (subProbs[i].owner != this)
	//     throw ProjectError("ViterbiColumnType: Submap owner mismatch");
	else
	    subProbs[i].dieOnZero();
    for (int i=0; i < maxsize; i++) {
	int ind = getSubstateIdx(i);
	if (ind != -1 && subProbs[ind].state != i)
	    throw ProjectError("ViterbiColumnType: internal error: integrity check failed (state=" + itoa(i) + " exists but the subProbs at index=" + itoa(ind) + " refer to another state).");
    }
    
}

void ViterbiColumnType::showSubstateInfo(ostream& strm) const {
    testIntegrity();
    
    for (int i=0; i < subProbs.size(); i++) 
	strm << "state " << (int)subProbs[i].state << ": " << subProbs[i].size() << " substates" << "\n";
    map<int,Double> sizes = substateCounts();
    pair<int,Double> argmax = make_pair(-1,0);
    for (map<int,Double>::iterator it = sizes.begin();
	 it!=sizes.end();
	 ++it) 
	if (it->second > argmax.second) 
	    argmax = (*it);
    try {
	const ViterbiSubmapType& mp = getSubstates(argmax.first);
	Double p = get(argmax.first);
	strm << "Biggest submap from state " << argmax.first << ":\n";
	strm << "has " << mp.size() << " entries.\n"
	     << "(main probability): " << p << ":\n";
	// strm << "(entries as factors to main prob:)\n";
	// for (ViterbiSubmapType::const_iterator it = mp.begin();
	//      it != mp.end();
	//      ++it) {
	//     if (it->first < 0)
	// 	strm << "[copy](" << it->first ;
	//     else {
	// 	int b = it->first >> 16;
	// 	int i =  ( it->first | ((-1) << 16) ) + ( 1 << 15);
	// 	strm << "(" << b << "/" << i;
	//     }
	//     strm << "): " <<  it->second * mp.factor / p << endl;
	// }
    } catch (NoSubmapFoundError) {
	int index = -1;
	for (int i=0; i < subProbs.size(); i++) 
	    if (subProbs[i].state == argmax.first)
		index = i;
	cerr << "index=";
	if (index >= 0)
	    cerr << index << "\n";
	else
	    cerr << "(not found)\n";
	cerr << "(should be " << (int)getSubstateIdx(argmax.first) << ")\n";
	cerr << "state=" << argmax.first << ".\n";
	throw;
    }
}
 

#define NODESIZE 20 // 20 bytes are allocated for a rb tree node

template<class S, class T>
inline int mapentrysize(const map<S,T>& mp) {
    return mp.size() * (sizeof(pair<S,T>)+NODESIZE);
}
template<class S, class T>
inline int mapsize(const map<S,T>& mp) {
    return sizeof(map<S,T>) + mapentrysize(mp);
}

GenericDataType ViterbiSubmapType::memoryUsage() {
    GenericDataType result;
    
    // heap usage of a ViterbiSubmap
    if (inactive())
	return result;

    result[0] = result[8] = 1;
    result[1] = theMap->size();
    result[2] = sizeof(ViterbiSubmapBasetype);
    result[3] = mapentrysize(*theMap);
    if (theMap->predSubstates) {
	result[4] = theMap->predSubstates->size();
	result[5] = mapsize(*theMap->predSubstates);
    }
    for (int i=1; i<MEMUSEDBLCOUNT; i++) 
	result[i] /= linkcount();
    return result;
}

GenericDataType ViterbiColumnType::getSubstateMemoryUsage() {
    GenericDataType result;
    for (int i=0; i<subProbs.size(); i++) 
	result += subProbs[i].memoryUsage();
    return result;
}
  
void ViterbiColumnType::countSubmaps(map<ViterbiSubmapBasetype*,int>& linkcounts) {
    for (int j=0; j<subProbs.size(); j++) {
	ViterbiSubmapBasetype* theMap = subProbs[j].theMap;
	if (theMap) {
	    int& newCount = linkcounts[theMap];
	    if (++newCount == theMap->linkcount)
		linkcounts.erase(theMap);
	}
    }
}
 
void ViterbiMatrixType::showSubstateMemoryUsage(int from, int to) {
    if (from < 0) from = 0;
    if (to < from) to = count;
    cerr << "Calculating memory usage for range [" << from 
	 << ".." << to << "]...";
    GenericDataType result;
    double& netsize = result[8];
    double& fullsize = result[0];
    int capmem = 0;
    for (int i=from; i<to; i++) {
	GenericDataType newresult = data[i].getSubstateMemoryUsage();
	result += newresult;
	capmem += data[i].getSubmapCapacity();
    }
   double mem 
	= result[2] // ViterbiSubmapBasetype
	+ result[3] // entries
	+ result[5] // predSubstates
	+ result[7];// successors
    cerr << "done." << endl;
    cerr << "Number of submap roots:            " << netsize << endl
	 << "Links per submap root:             " << fullsize/netsize << endl
	 << "Submap entries per base:           " << result[1]/count << endl
	 << "Capacity for submaps:              " << (capmem/1024.0) << "KB.\n"
	 <<  endl
	 << "  basetype:              " << (result[2]/1024.0) << "KB.\n"
	 << "  entries:               " << (result[3]/1024.0) << "KB.\n"
	 << "  predecessor substates: " << (result[5]/1024.0) << "KB.\n"
	 << "  successor counts:      " << (result[7]/1024.0) << "KB.\n"
	 << "total variable submap memory:      " << (mem/1024.0) << "KB.\n"
	 << "capacity + total variable:         " << ((mem + capmem)/1024.0) << "KB.\n";

    cerr << "(per submap root)" << endl
	 << "  basetype:              " << result[2] / netsize << " Bytes.\n"
	 << "  entries:               " << result[3] / netsize << " Bytes.\n"
	 << "  predecessor substates  " << result[5] / netsize << " Bytes.\n"
	 << "  successor counts       " << result[7] / netsize << " Bytes.\n"
	 << "total (per submap root):           " << mem/netsize << " Bytes.\n"
	 << "Average substate counts: " << endl
	 << "  entries:               " << (result[1] / netsize) << endl
	 << "  predecessor substates: " << (result[4] / netsize) << endl
	 << "  successor count:       " << (result[6] / netsize) << endl;
    
    map<ViterbiSubmapBasetype*, int> checkedLinkCounts;
    for (int i=0; i<count; i++)
	data[i].countSubmaps(checkedLinkCounts);
    if (!checkedLinkCounts.empty()) {
	ViterbiSubmapBasetype* problem = checkedLinkCounts.begin()->first;
	cerr << "Problem @ " << problem;
	throw ProjectError("linkcount to large");
    }
}
#endif

