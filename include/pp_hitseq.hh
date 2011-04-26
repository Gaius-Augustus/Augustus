/*****************************************************************************\
 * Filename : pp_hitseq.hh
 * Author   : Oliver Keller
 * Project  : Gene Prediction with Protein Family Patterns
 *
 * Description: Data structures for storing ProteinProfile hit sequences
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 20.10.08   | Oliver Keller         | creation of the file
 * 12.11.08   | Oliver Keller         | first version in SVN (rev. 137)
 \******************************************************************************/

#ifndef __PP_HITSEQ_HH
#define __PP_HITSEQ_HH

#include "pp_profile.hh"

#include <algorithm>

namespace PP {
    /*
     * HitSequence:
     * consecutive sequence of block hits in the target sequence
     * factor is the quotient of probabilities (block model vs. exon model)
     * firstBlockIndex is the number of the first block of the sequence,
     * posns contains positions from for block numbers first() to last()
     */
    struct HitSequence {
	HitSequence(int b, int from, const Double& p) :
	    factor(p),
	    firstBlockIndex(b),
	    posns(1, from) 
	{}
	int first() const {
	    return firstBlockIndex;
	}
	int last() const {
	    return firstBlockIndex + posns.size() - 1;
	}
	int firstValue() const {
	    return posns.front();
	}
	int lastValue() const {
	    return posns.back();
	}

	// append a single hit to the sequence
	void add(int from, Double p) {
	    posns.push_back(from);
	    factor *= p;
	}
	// append all hits from another sequence
	void add(const HitSequence& from) {
	    copy(from.posns.begin(), from.posns.end(), back_inserter(posns));
	    factor *= from.factor;
	}

	int operator[] (int b) const {
	    return posns[b - firstBlockIndex];
	}
        Double getFactor() const {
	    return factor;
	}

	// position in the profile for a base coordinate specified by offset; returns
	// (b,i) where b is the last block and i is the distance of lastValue() to the
	// offset. i equals (block length + distance of block end to offset). Usually
	// offset is the dna coordinate of the exon end.
	Position targetPosAt(int offset) const {
	    return Position(last(), (offset +1 - lastValue())/3);
	}
 
	Double factor;
	int firstBlockIndex;
	vector<int> posns;
    }; // PP::HitSequence


    /*
     * HitSequenceNode: data structure for navigation through a 
     * collection of hit sequences
     */
    class HitSequenceNode {
    protected:
	/*
	 * iterators for HitSequenceNode's
	 */
	// ordinary iterator (ordering by first block)
	class iterator {
	public:
	    HitSequence& operator*() {
#ifdef DEBUG
		if (node->item == 0)
		    throw ProjectError("invalid iterator end() dereferenced");
#endif
		return *(node->item);
	    }
	    HitSequence* operator->() {
#ifdef DEBUG
		if (node->item == 0)
		    throw ProjectError("invalid iterator end() dereferenced");
#endif
		return node->item;
	    }
	    bool operator==(iterator other) {
		return node==other.node;
	    }
	    iterator& operator=(iterator other) {
		node=other.node;
		return *this;
	    }
	    bool operator!=(iterator other) {
		return node!=other.node;
	    }
	    iterator& operator++() {
		node = node->next_by_first;
		return *this;
	    }
	    iterator& operator--() {
		node = node->previous_by_first;
		return *this;
	    }
	    void erase() {
#ifdef DEBUG
		if (node->item == 0)
		    throw ProjectError("tried to delete root");
#endif
		node = node->next_by_first;
		delete node->previous_by_first;
	    }
	    friend iterator createNode(const HitSequence&, iterator, iterator);
	    friend iterator createNode(int, int, const Double&, iterator, iterator);
	    iterator(HitSequenceNode* cur) : node(cur) {}
	protected:
	    HitSequenceNode* node; // pointer to node containing current hit sequence
	}; // PP::HitSequenceNode::iterator

	// reverse iterator (ordering by first block)
	struct reverse_iterator : iterator {
	    reverse_iterator(const iterator& it) : iterator(it) {}
	    HitSequence& operator*() {
#ifdef DEBUG
		if (node->previous_by_first->item == 0)
		    throw ProjectError("invalid iterator rend() dereferenced");
#endif
		return *(node->previous_by_first->item);
	    }
	    HitSequence* operator->() {
#ifdef DEBUG
		if (node->previous_by_first->item == 0)
		    throw ProjectError("invalid iterator rend() dereferenced");
#endif
		return node->previous_by_first->item;
	    }
	    reverse_iterator& operator++() {
		node = node->previous_by_first;
		return *this;
	    }
	    reverse_iterator operator++(int) {
		node = node->previous_by_first;
		return reverse_iterator(node->next_by_first);
	    }
	    reverse_iterator& operator--() {
		node = node->next_by_first;
		return *this;
	    }
	    reverse_iterator operator--(int) {
		node = node->next_by_first;
		return iterator(node->previous_by_first);
	    }
	    void erase() {
		delete node->previous_by_first;
	    }
	}; // PP::HitSequenceNode::reverse_iterator

	// iterator for ordering by last block
	struct last_iterator : iterator {
	    last_iterator(const iterator& it) : iterator(it) {}
	    last_iterator& operator++() {
		node = node->next_by_last;
		return *this;
	    }
  	    last_iterator operator++(int) {
 		node = node->next_by_last;
 		return iterator(node->previous_by_last);
 	    }
	    last_iterator& operator--() {
		node = node->previous_by_last;
		return *this;
	    }
	    void erase() {
#ifdef DEBUG
		if (node->item == 0)
		    throw ProjectError("tried to delete root");
#endif
		node = node->next_by_last;
		delete node->previous_by_last;
	    }
	}; // PP::HitSequenceNode::last_iterator

    protected:
	// constructors for HitSequenceNode
	// can only be created by HitSequenceList
	// or by createNode
	HitSequenceNode() : 
	    item(0),
	    previous_by_first(this),
	    next_by_first(this),
	    previous_by_last(this),
	    next_by_last(this) {
#ifdef DEBUG
	    if (previous_by_first != this || next_by_first != this || previous_by_last != this)
		throw ProjectError("This shouldn't happen.");
#endif
	}
	HitSequenceNode(const HitSequence& seq, 
			HitSequenceNode* by_first, HitSequenceNode* by_last) :
	    item(new HitSequence(seq)) 
	{
	    insert_this(by_first, by_last);
	}	    
	HitSequenceNode(int b, int i, const Double& p, 
			HitSequenceNode* by_first,
			HitSequenceNode* by_last) : 
	    item(new HitSequence(b,i,p))
	{
	    insert_this(by_first, by_last);
	}
	// destructor
	~HitSequenceNode() {
	    remove_this();
	    delete item;
	}

 	void insert_this(HitSequenceNode* by_first,
			 HitSequenceNode* by_last) {
	    previous_by_first = by_first->previous_by_first;
	    next_by_first = by_first;

	    previous_by_last = by_last->previous_by_last; 
	    next_by_last = by_last;

	    by_first->previous_by_first
		= previous_by_first->next_by_first
		= by_last->previous_by_last
		= previous_by_last->next_by_last 
		= this;
	}
	void remove_this() {
	    previous_by_first->next_by_first = next_by_first;
	    next_by_first->previous_by_first = previous_by_first;
	    previous_by_last->next_by_last = next_by_last;
	    next_by_last->previous_by_last = previous_by_last;
	}
#ifdef DEBUG
	bool consistent() {
	    return 
		this == previous_by_first->next_by_first && 
		this == next_by_first->previous_by_first && 
		this == previous_by_last->next_by_last && 
		this == next_by_last->previous_by_last;
	}
#endif

	
	friend iterator createNode(const HitSequence& hs, iterator it, iterator lit) {
	    return new HitSequenceNode(hs, it.node, lit.node);
	}
	friend iterator createNode(int b, int i, const Double& p, iterator it, iterator lit) {
	    return new HitSequenceNode(b, i, p, it.node, lit.node);
	}


	HitSequence* item;                   // hit sequence of the node
	HitSequenceNode* previous_by_first;  // hit sequences ordered by first block
	HitSequenceNode* next_by_first;
	HitSequenceNode* previous_by_last;   // hit sequences ordered by last block
	HitSequenceNode* next_by_last;
    };  // PP::HitSequenceNode


    /*
     * HitSequenceList:
     * this class is an entry point to a hit sequence collection
     * the base node of the list is without HitSequence
     */
    class HitSequenceList : private HitSequenceNode {
//	friend class HitSequenceCollection;
      public:
	typedef HitSequenceNode::iterator iterator;
	typedef HitSequenceNode::reverse_iterator reverse_iterator;
	typedef HitSequenceNode::last_iterator last_iterator;

	// constructor and destructor
	HitSequenceList() : 
	    HitSequenceNode(),
	    activeMax(this) {
#ifdef DEBUG
	    if (previous_by_first != this)
		throw ProjectError("This shouldn't happen!");
#endif
	}
	~HitSequenceList() {
	    while (next_by_first != this)
		begin().erase();
	    while (next_by_last != this)
		lbegin().erase();
	}

	// iterator methods
	iterator begin() {
	    return next_by_first;
	}
	iterator end() {
	    return this;
	}
	reverse_iterator rbegin() {
	    return end();
	}
	reverse_iterator rend() {
	    return begin();
	}
	last_iterator lbegin() {
	    return iterator(next_by_last);
	}
	last_iterator lend() {
	    return end();
	}

	void init() {
	    activeMax = rbegin();
	}
	bool inactive() {  
	    return activeMax == rend();
	}
	void setInactive(reverse_iterator it) {
	    activeMax = it;
	}
	reverse_iterator rbeginActive() {
	    return activeMax;
	}

	bool empty() {
	    return begin() == end();
	}

#ifdef DEBUG
	bool inconsistent() {
	    return item != 0;
	}
#endif
    private:
	reverse_iterator activeMax; // end of the active elements of the list
    }; // PP::HitSequenceList

 	    
    /*
     * HitSequenceCollection:
     * this class contains all hit sequences for a given reading frame
     * ordered by block number (one HitSequenceList for each block number
     * to access them ordered by first and last block positions)
     */
    class HitSequenceCollection {
    public:

	HitSequenceCollection(int size=0) {
	    reset(size);
	}
	~HitSequenceCollection() {
	    deleteNodes();
	}
	void addSingle(int b, int offset, const Double& val) {
	    createNode(b, offset, val, end(b), lbegin(b));
	}
	void addSingleFront(int b, int offset, const Double& val) {
	    createNode(b, offset, val, begin(b), lend(b));
	}
	void addHit(int b, int offset, Double val, DistanceType iBD, int blocksize);
	void addHitFront(int b, int offset, Double val, DistanceType iBD);
	    
	void init() {
	    for (int b=0; b<lists.size(); b++)
		lists[b]->init();
	}
	void deleteNodes() {
	    for (int b=0; b<lists.size(); b++)  
		delete lists[b];
	}

    private:
	static HitSequenceList* newlist() {
	    return new HitSequenceList();
	}
    public:
	void reset(int size) {
	    deleteNodes();
	    lists.resize(size);
	    generate_n(lists.begin(), size, &newlist);
	}
	void reset() {
	    deleteNodes();
	    generate(lists.begin(), lists.end(), &newlist);
	}

	// iterator methods
	HitSequenceList::iterator begin(int b) {
	    return list(b).begin();
	}
	HitSequenceList::iterator end(int b) {
	    return list(b).end();
	}
	HitSequenceList::reverse_iterator rbegin(int b) {
	    return list(b).rbegin();
	}
	HitSequenceList::reverse_iterator rend(int b) {
	    return list(b).rend();
	}
	HitSequenceList::last_iterator lbegin(int b) {
	    return list(b).lbegin();
	}
	HitSequenceList::last_iterator lend(int b) {
	    return list(b).lend();
	}

	HitSequenceList& list(int b) {
#ifdef DEBUG
	    if (lists[b]->inconsistent())
		throw ProjectError("root has nonzero item!");
#endif
	    return *lists.at(b);
	}
    private:
	vector<HitSequenceList*> lists;
   }; // PP::HitSequenceCollection

} // namespace PP

#endif
