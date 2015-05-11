/***********************************************************************
 * file:    vitmatrix.hh
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  data structures used inside the viterbi matrix
 * author:  Oliver Keller (keller@cs.uni-goettingen.de)
 *
 **********************************************************************/

#ifndef _VITMATRIX_HH
#define _VITMATRIX_HH

// project includes
#include "types.hh"  // for ProjectError

// standard C/C++ includes
#include <vector>
#include <map>
#include <list>
#include <limits>
#include <climits>

using namespace std;


// constants for Viterbi*Types
#define VIT_MIN_CAPACITY 7
#define MAX_LINKCOUNT 300
#define SUBSTATEVALUE_BITS (8*sizeof(short))
	
struct SubstateId {
    explicit SubstateId(signed char s=-1, short v=0) : slot(s), value(v) {}

    bool operator< (const SubstateId& other) const {
	return 
	    slot < other.slot || (slot == other.slot &&
				  value < other.value);
    }
    bool operator== (const SubstateId& other) const {
	return slot == other.slot && value == other.value;
    }
    bool operator> (const SubstateId& other) const {
	return (other < *this);
    }
    bool operator!= (const SubstateId& other) const {
	return !operator==(other);
    }
    bool operator<= (const SubstateId& other) const {
	return !operator>(other);
    }

    bool undef() {
	return slot==-1 && value==0;
    }
    
    

    int fullId() const {
	// implements a monotonous mapping
	// with SubstateId() -> 0
	return ((slot+1) << SUBSTATEVALUE_BITS) + value;
    }
    void set(int fromId) {
	// reverse mapping int->SubstateId()
	value = fromId;
	slot = ((fromId-value) >> SUBSTATEVALUE_BITS) -1;
    }
    static SubstateId min(int slot) {
	return SubstateId(slot, numeric_limits<short>::min());
    }
    static SubstateId max(int slot) {
	return SubstateId(slot, numeric_limits<short>::max());
    }

    signed char slot;
    short value;
};

#ifdef DEBUG
inline string itoa(SubstateId d) {
    return "(" + itoa(d.slot) + "/" + itoa(d.value) + ")";
}
#endif

typedef map<SubstateId, SubstateId> SubstateIdMap;

#ifdef DEBUG
#define MEMUSEDBLCOUNT 10
struct GenericDataType {
    GenericDataType() {
	for (int i=0; i<MEMUSEDBLCOUNT; i++)
	    vals[i]=0;
    }
    void operator+=(const GenericDataType& other) {
	for (int i=0; i<MEMUSEDBLCOUNT; i++)
	    vals[i]+=other.vals[i];
    }
    double& operator[] (int i) { return vals[i]; }
    double vals[MEMUSEDBLCOUNT];
};
#endif

// exceptions thrown by Viterbi
struct NoSubmapFoundError : public ProjectError {
    NoSubmapFoundError() : 
	ProjectError("Internal error: ViterbiSubmapType has no map!") {}
};
struct SubmapEmptyError : public ProjectError {
    SubmapEmptyError() : 
	ProjectError("Internal error: ViterbiSubmapType is empty!") {}
};


/*
 * AlgorithmVariant: 
 * this type is used to determine the variant implemented in the
 * viterbiForwardAndSampling methods
 *
 */
enum AlgorithmVariant { doViterbiOnly, doViterbiAndForward, doSampling, doBacktracking };
/* any changes made to the order above needs to be reflected in the following 
   functions */
inline bool doViterbi(AlgorithmVariant algovar) {
    return algovar <= doViterbiAndForward;
}
inline bool needForwardTable(AlgorithmVariant algovar) {
    return algovar == doViterbiAndForward || algovar == doSampling;
}
inline AlgorithmVariant doViterbi(bool needForwardTable) {
    return needForwardTable ? doViterbiAndForward : doViterbiOnly;
}


/*
 * ViterbiSubmapType:
 * This class stores the probabilites for substates for a single
 * pair (base position, state)
 * 
 * the actual values are the entries of the map, multiplied by a common factor
 */
class ViterbiSubmapType;

struct ViterbiSubmapEntry {
    ViterbiSubmapEntry(const Double& dbl=0, ViterbiSubmapType* pred=0) :
	p(dbl), predMap(pred), succCount(0) {}
    bool operator<(const ViterbiSubmapEntry& other) const {
	return p<other.p;
    }
    bool operator<=(const ViterbiSubmapEntry& other) const {
	return !(other < *this);
    }
    bool operator>=(const ViterbiSubmapEntry& other) const {
	return !operator<(other);
    }
    bool operator>(const ViterbiSubmapEntry& other) const {
	return other < *this;
    }
    void set(const Double& dbl, ViterbiSubmapType* pred) {
	p =dbl;
	predMap = pred;
    }

    Double p;
    ViterbiSubmapType* predMap;
    short succCount;
};

struct ViterbiSubmapBasetype : public map<SubstateId, ViterbiSubmapEntry > {
    ViterbiSubmapBasetype(int ln) : 
	linkcount(ln), 
	predSubstates(0) {}
    ViterbiSubmapBasetype(const ViterbiSubmapBasetype& other, int ln) :
	map<SubstateId,ViterbiSubmapEntry>(other),
	linkcount(ln),
	predSubstates(0) {}
    ~ViterbiSubmapBasetype() {
	delete predSubstates;
    }

    // multiple instances of ViterbiSubmapType can share the same values;
    // this value tells how many of them exist, so ViterbiSubmapType
    // knows when to delete the pointed to 
    int linkcount;

    // this map stores the predecessor substates, if needed
    SubstateIdMap* predSubstates;
}; // ViterbiSubmapBasetype

class ViterbiSubmapType {
    friend class ViterbiColumnType;
    // ``const'' in this class refers not to the contents of the map
    // but just to the pointer (we want to be able to change
    // substates of otherwise ``constant'' Viterbi columns)
 
public:
    // these are aliases for a map iterator
    typedef ViterbiSubmapBasetype::iterator iterator;
    typedef ViterbiSubmapBasetype::const_iterator const_iterator;

    // constructors
    explicit ViterbiSubmapType(signed char st) : 
	factor(1), state(st), theMap(0) {}
    // copy constructor: does not copy the map pointed to but
    // just adds a link to it
    ViterbiSubmapType(const ViterbiSubmapType& other) :
	factor(other.factor), state(other.state)
    {
	link(other);
    }
    // destructor
    ~ViterbiSubmapType() {
	unlink();
    }
    // assignment: the same behaviour as the copy constructor
    void operator= (const ViterbiSubmapType& other) {
	if (&other == this)
	    return;
	factor = other.factor;
	state = other.state;
	unlink();
	link(other);
    }
    // share another submap's map
    void link (const ViterbiSubmapType& other, Double mult) {
	factor = other.factor * mult;
	unlink();
	link(other);
    }
    // remove the map
    void deactivate() {
	unlink();
	theMap=0;
    }
    int getMainState() const {
	return state;
    }


    // the following methods can be used also when the map is NULL
    bool inactive() const {  // is the map NULL?
	return (theMap==0);
    }
    bool empty() const {
     	return inactive() || theMap->empty();
    }
    bool emptypart(int part) {
    	return inactive() || begin(part) == end(part);
    }
    int size() const {
    	return inactive() ? 0 : theMap->size();
    }
    void activate() {
	if (!theMap)
	    theMap = new ViterbiSubmapBasetype(1);
    }
    int linkcount() const {
	return inactive() ? 0 : theMap->linkcount;
    }
    bool is_linked() {
	return theMap && theMap->linkcount > 1;
    }
    // access predecessor substate map
    SubstateIdMap* getPredSubstateMap() const; 
    SubstateId popPredSubstate(SubstateId substate) const;
    void clearEmptySubstateMap() const;

    void swapMaps(ViterbiSubmapType&);       // exchange values with another instance
 
    void merge(ViterbiSubmapType&, Double);  // merge two maps and take the maximum entries
    int eraseAllUnneededSubstates() const;   // erase substates with succCount==0

    // the following methods do not check if theMap is non-NULL
    iterator begin() const {
    	return theMap->begin();
    }
    iterator begin(int slot) const {
	return theMap->lower_bound(SubstateId::min(slot));
    }
    iterator end() const {
	return theMap->end();
    }
    iterator end(int slot) const {
	return theMap->upper_bound(SubstateId::max(slot));
    }
    iterator bound(SubstateId substate) const {
     	return theMap->lower_bound(substate);
    }
    iterator bound(signed char slot, short value) const {
    	return theMap->lower_bound(SubstateId(slot,value));
    }
    Double get(SubstateId substate) const {
	const_iterator it = theMap->find(substate);
	return it==end() ? 0 : it->second.p * factor;
    }
    iterator insert(iterator where, pair<SubstateId, ViterbiSubmapEntry> what) const {
	// insert the value as is (without dividing it by factor)
    	return theMap->insert(where, what);
    }
    iterator insert(iterator where, SubstateId substate, const Double& value, 
		    ViterbiSubmapType* pred) const {
	// insert the value as is (without dividing it by factor)
    	return insert(where, make_pair(substate, ViterbiSubmapEntry(value,pred)));
    }
    void push_back(SubstateId substate, const Double& value, ViterbiSubmapType* pred) const {
	// insert at end of map
    	insert(end(), substate, value, pred);
    }
    void erase(iterator start, iterator end) const {
	theMap->erase(start, end);
    }
    void erase(SubstateId substate) const {
	theMap->erase(substate);
    }
    short& succCount(SubstateId substate) {
	return (*theMap)[substate].succCount;
    }
    short& succCount(iterator it) {
	return it->second.succCount;
    }
    const short& succCount(iterator it) const {
#ifdef DEBUG
	if (it == end())
	    throw SubmapEmptyError();
#endif
	return it->second.succCount;
    }
    short& succCount() {
    	return succCount(begin());
    }
    void eraseAndDecCounts(iterator it) const; 
    bool setMax(SubstateId substate, Double value, ViterbiSubmapType* pred) const;
    void cloneMap(ViterbiSubmapType* newPredMap, int maxlinkcount=1);
    void cloneMapIfLinkcountExceeded(ViterbiSubmapType* newPredMap) {
	cloneMap(newPredMap, MAX_LINKCOUNT);
    }
    void getMaxSubstate(SubstateId& substate, Double& value) const;
    int eraseUnneededSubstate(SubstateId substate);
    void registerPredecessors();
    void useMaximum(Double& mainProb);

#ifdef DEBUG
    // functions used for debugging only
    GenericDataType memoryUsage();
    void dieOnZero() const;
    void checkAfterWindow() const;
    SubstateId getPredecessorSubstate(SubstateId substate) const;
    void checkCounts() const;
#endif

    Double factor; // multiply this to the entries of the map to get the actual value

private:
    // common to copy constructor and assignment operator:
    void link(const ViterbiSubmapType& other) {
	theMap = other.theMap;
	if (theMap) theMap->linkcount++;
    }
    // common to destructor and assignment operator
    void unlink() {
	if (theMap && --(theMap->linkcount)==0)
	    delete theMap;
    }

    signed char state;             // the main state
    ViterbiSubmapBasetype* theMap; // the entries (probabilities/factor) for each substate
}; // ViterbiSubmapType
    
// inline definitions of ViterbiSubmapType member functions

// returns the predecessor substate map, creates one if none exists
inline SubstateIdMap* ViterbiSubmapType::getPredSubstateMap() const {
    if (!theMap)
	return 0;
    SubstateIdMap*& result = theMap->predSubstates;
    if (!result)
	result = new SubstateIdMap;
    return result;
}

inline void ViterbiSubmapType::clearEmptySubstateMap() const {
    if (theMap && theMap->predSubstates &&
	theMap->predSubstates->empty()) {
	delete theMap->predSubstates;
	theMap->predSubstates = 0;
    }
}

inline void ViterbiSubmapType::swapMaps(ViterbiSubmapType& other) {
    ViterbiSubmapBasetype* temp = theMap;
    theMap = other.theMap;
    other.theMap = temp;
}

inline void ViterbiSubmapType::eraseAndDecCounts(iterator it) const {
    ViterbiSubmapType* pred = it->second.predMap;
    if (pred) {
	SubstateId substate = popPredSubstate(it->first);
	pred->succCount(substate)--;
#ifdef DEBUG
	if (pred->succCount(substate) < 0)
	    throw ProjectError("Have negative succCounts");
#endif
    }
    theMap->erase(it);
}

inline bool ViterbiSubmapType::setMax(SubstateId substate, Double value, ViterbiSubmapType* pred)  const {
    ViterbiSubmapEntry& target = (*theMap)[substate];
    if (target.p * factor < value) {
	target.set(value / factor, pred);
	return true;
    }
    return false;
}

inline void ViterbiSubmapType::cloneMap(ViterbiSubmapType* newPredMap, int maxlinkcount) {
#ifdef DEBUG
    if (newPredMap == 0)
	throw ProjectError("SubstateModel::cloneMap: predMap==0!");
#endif
   if (theMap->linkcount > maxlinkcount) {
	theMap->linkcount--;
	theMap = new ViterbiSubmapBasetype(*theMap,1);
	for (iterator it=begin(); it != end(); ++it) {
	    it->second.predMap = newPredMap;
	    it->second.succCount = 0;
	}
    }
}

inline void ViterbiSubmapType::getMaxSubstate(SubstateId& substate, Double& value) const {
    value = 0; substate.set(0);
    if (empty()) 
	return;
    for (iterator it = begin(); it != end(); ++it)
	if (it->second > value) {
	    substate = it->first;
	    value = it->second.p;
	}
    value *= factor;
}


/*
 * ViterbiEntryType:
 * one entry in the viterbi matrix, including the index for the substate map
 *
 */
struct ViterbiEntryType {
    ViterbiEntryType(signed char st, Double v, signed char i=-1) :
	state(st), value(v), substate_idx(i) {}
    signed char state;
    Double value;
    signed char substate_idx;
};


/*
 * ViterbiColumnType:
 * The entries of the Viterbi matrix belonging to the same base position
 * Only nonzero values are actually stored
 *
 */
typedef vector<ViterbiEntryType> ViterbiColumnBasetype;

class ViterbiColumnType : public ViterbiColumnBasetype {
public:
    
    // constructor and destructor
    ViterbiColumnType() : 
	ViterbiColumnBasetype(), 
	subProbs(), idx(0), maxsize(0) {}
    ~ViterbiColumnType()  {
	delete[] idx;
    }

    void reset(int n);     // delete the column and reserve space for n states
    void erase(int state); // erase the value for a particular state (including substates)

    // return a value without autovivification
    Double get(int state) const {
	return has(state) ? elem(state).value : 0;
    }
    Double get(int state, SubstateId substate) const {
	if (substate.undef()) return get(state);
	try {
	    return getSubstates(state).get(substate);
	} catch (NoSubmapFoundError) {
	    return 0;
	}
    }
    // index operator; uses autovivication in non-const contexts
    Double operator[](int state) const { // no autovivification
	return get(state);
    }
    Double& operator[](int state) { // autovivification if needed
	if (has(state))
	    return elem(state).value;
	addState(state);
	return back().value;
    }

    bool has(int state) const {
#ifdef DEBUG
	if (state < 0 || state >= maxsize)
	    throw ProjectError("ViterbiColumnType: state out of range");
#endif	
	return idx[state] >= 0;
    }
    int size() { // number of elements !=0
	return ViterbiColumnBasetype::size();
    }
    int max_capacity() { // this equals statecount
	return maxsize;
    }

    // access substates
    bool hasSubstates(int state) const {
	int idx = getSubstateIdx(state);
	return idx >= 0 && idx < subProbs.size() && !subProbs[idx].inactive();
    }
    const ViterbiSubmapType& getSubstates(int state) const;
    ViterbiSubmapType& getSubstates(int state);
    void getMaxSubstate(int state, SubstateId& substate, Double& value) const;
    ViterbiSubmapType& addSubstates(const ViterbiSubmapType& newmap);

    // calls a member function of class Tp with a parameter of
    // type ViterbiSubmapType for each submap in the column
    template <class Tp, class FncPtr>
    void foreachSubstate(Tp p, FncPtr f) {
	for (int i=0; i<subProbs.size(); i++) 
	    (p->*f)(subProbs[i]);
    }
    void eraseSubstates(int state);  // erase all the substates for a particular state
    int eraseUnneededSubstates();    // calls eraseUnneededSubstates for each submap
    int removeEmptySubmaps();        

#ifdef DEBUG
    // functions used for debugging only
    int substateCount(int state) const {
	try {
	    return getSubstates(state).size();
	} catch (...) {
	    return 0;
	}
    }
    int totalSubstateCount() const {
	int result=0;
	for (int i=0; i<subProbs.size(); i++) result += subProbs[i].size();
	return result;
    }
    map<int,Double> substateCounts() const {
	map<int,Double> result;
	for (int i=0; i<subProbs.size(); i++) 
	    result[subProbs[i].state] = subProbs[i].size();
	return result;
    }
    void testIntegrity() const;
    void checkAfterWindow() const {
	for (int i=0; i < subProbs.size(); i++) 
	    subProbs[i].checkAfterWindow();
    }
    void showSubstateInfo(ostream& strm) const;
    void write(ostream& strm, const vector<string>&);
    void write(const vector<string>&); 

    int getSubmapCapacity() {
	return subProbs.capacity() * sizeof(ViterbiSubmapType);
    }
    GenericDataType getSubstateMemoryUsage();
    void countSubmaps(map<ViterbiSubmapBasetype*,int>& linkcounts);
#endif
    	    
private:
    const ViterbiEntryType& elem(int state) const {
	return ViterbiColumnBasetype::operator[](idx[state]);
    }
    ViterbiEntryType& elem(int state) {
	return ViterbiColumnBasetype::operator[](idx[state]);
    }
    void addState(int state, signed char sub_idx=-1);
    signed char getSubstateIdx(int state) const {
	return has(state) ? elem(state).substate_idx : -1;
    }
    void eraseSubProbs(int subidx);

    vector<ViterbiSubmapType> subProbs; // the substate probabilities
    signed char* idx;                   // for each state, the vector index
                                        // or -1 if entry is zero
    int maxsize;                        // number of states
};  // ViterbiColumnType

// inline definitions of ViterbiColumnType member functions
inline const ViterbiSubmapType& ViterbiColumnType::getSubstates(int state) const {
    int idx = getSubstateIdx(state);
    if (idx >= 0 && idx < subProbs.size()) {
	const ViterbiSubmapType& result = subProbs[idx];
	if (!result.inactive())
	    return result;
    }
    throw NoSubmapFoundError();
}

inline ViterbiSubmapType& ViterbiColumnType::getSubstates(int state) {
    int idx = getSubstateIdx(state);
    if (idx >= 0 && idx < subProbs.size()) {
	ViterbiSubmapType& result = subProbs[idx];
	if (!result.inactive())
	    return result;
    }
    throw NoSubmapFoundError();
}

inline void ViterbiColumnType::getMaxSubstate(int state, SubstateId& substate, Double& value) const {
    Double mainval = get(state);
    try {
	const ViterbiSubmapType& submap = getSubstates(state);
	submap.getMaxSubstate(substate, value);
	if (value > mainval)
	    return;
    } catch (NoSubmapFoundError) {}
    substate.set(0); value = mainval;
}

inline int ViterbiColumnType::eraseUnneededSubstates() {
    int result = 0;
    for (int i=0; i<subProbs.size(); i++) 
	result += subProbs[i].eraseAllUnneededSubstates();
    return result;
}

inline int ViterbiColumnType::removeEmptySubmaps() {
    int result = 0;
    for (int i=0; i<subProbs.size(); i++)
	if (subProbs[i].empty() && !subProbs[i].inactive()) {
	    subProbs[i].deactivate();
	    result++;
	}
    return result;
}

inline void ViterbiColumnType::eraseSubProbs(int subidx) {
    if (subidx >= 0) {
	// overwrite entry at position subidx
	// with the (previously) last entry
	if (subidx < subProbs.size()-1) {
	    elem(subProbs.back().state).substate_idx = subidx;
	    subProbs[subidx] = subProbs.back();
	}
	// delete last entry in vector
	subProbs.pop_back();
    }
}


/*
 * ViterbiMatrixType:
 * An array of Viterbi columns
 */
class ViterbiMatrixType {
public:
    ViterbiMatrixType () : data(0), count(0) {}
    ~ViterbiMatrixType() {
	delete[] data;
    }
    void assign(int colcount, int colsize) {
	count = colcount;
	delete[] data;
	data = new ViterbiColumnType[count];
	for (int i=0; i<count; i++) {
	    data[i].reset(colsize);
	}
    }
    ViterbiColumnType& operator[] (int i) {
	return data[i];
    }
    const ViterbiColumnType& operator[] (int i) const {
	return data[i];
    }
    int size() const {
	return count;
    }

#ifdef DEBUG
    // functions used for debugging only
    double load() {
	double result=0;
	for (int i=0; i<count; i++) 
	    result += data[i].capacity();
	return result / count;
    }
    double load(int i) {
	return data[i].capacity();
    }
    double used_load() {
	double result=0;
	for (int i=0; i<count; i++) 
	    result += data[i].size();
	return result / count;
    }
    void write(ostream& strm, const vector<string>& stateNames = vector<string>()) {
	for (int i=0; i<count; i++) 
	    data[i].write(strm, stateNames);
    }
    int compare(istream&, const vector<string>& stateNames = vector<string>());
    void showSubstateMemoryUsage(int from = 0, int to = -1);
#endif
		
private:
    ViterbiColumnType* data;
    int count;
}; // ViterbiMatrixType


/*
 * Options lists are used for sampling; items also in backtracking
 */
class OptionListItem {
public:
    OptionListItem(){reset();};
    OptionListItem(int st, int bs, Double p){
	state = st;
	base = predEnd = bs;
	probability = p;
    }
    OptionListItem(int st, int bs, int pe, Double p){
	state = st;
	base = bs;
	predEnd = pe;
	probability = p;
    }
    void reset(){state = TYPE_UNKNOWN; base = predEnd = -INT_MAX;}
    void resetPredEnd() {predEnd = -INT_MAX;}
    int state;
    int base;
    int predEnd; // base at which predecessor state ends, usually identical to base unless states (genes) overlap
    Double probability;
};

inline bool operator< (const OptionListItem& first, const OptionListItem& second) {
    // the largest probability comes first
    return (first.probability > second.probability);
}

class OptionsList{
public:
    OptionsList(){
	cumprob = 0.0;
    }
    ~OptionsList(){
	while (!options.empty())
	    options.pop_front();
    }
    void add(int state, int base, Double probability, int predEnd = -INT_MAX) {
	options.push_back(OptionListItem(state, base, 
					 (predEnd == -INT_MAX)? base : predEnd,
					 probability));
	cumprob += probability;
    }
    void prepareSampling() {
	options.sort();
    }
    int size(){
	return options.size();
    }
    void print(ostream& out) {
	out << "options: " << cumprob << " | "; 
	for (list<OptionListItem>::iterator it = options.begin(); it != options.end(); ++it){
	    out << it->state << " " << it->base << " " << (it->probability/cumprob) << ", ";
	}
	out << endl;
    }
    OptionListItem sample();

private:
    list<OptionListItem> options;
    Double cumprob;
};


#endif    // _VITMATRIX_HH
