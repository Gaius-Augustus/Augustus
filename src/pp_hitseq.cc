/**********************************************************************
 * file:    pp_hitseq.cc
 * license: Artistic License, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  protein pattern search extension
 * authors: Oliver Keller
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 20.10.08| Oliver Keller | creation of the file
 * 12.11.08| Oliver Keller | rev. 137, first version in the SVN
 **********************************************************************/

#include "pp_hitseq.hh"


using namespace PP;
	
void HitSequenceCollection::addHit(int b, int offset, Double val, DistanceType iBD, int blocksize) {
    addSingle(b, offset, val);
    if (b==0) return;
    // HitSequenceList& hsList = lists[b-1];
    HitSequenceList::last_iterator 
	it = lbegin(b-1),
	end = lend(b-1);
    if (it == end) return;
    int good_offset = offset-iBD.r.min*3;
    while (it->lastValue() > good_offset) {
	if (++it == end)
	    return;
    }
    good_offset = iBD.has_max ? offset-iBD.r.max*3 : 0;
    while (it->lastValue() >= good_offset) {
	// insert a new node, by_first before the current one, 
	// by_last to the beginning of the list (highest offsets first)
	createNode(*it, it, lbegin(b))->add(offset, val);
	if (++it == end)
	    return;
    }
    good_offset -= 3*blocksize;
    while (it->lastValue() >= good_offset) 
	if (++it == end)
	    return;
    do {
	it.erase();
    } while (it != end);
}

void HitSequenceCollection::addHitFront(int b, int offset, Double val, DistanceType iBD) {
    HitSequenceList* hsSource = lists.at(b+1);
    HitSequenceList::reverse_iterator it = hsSource->rbeginActive(), end = hsSource->rend();
    
    while (it != end && iBD.has_max && 
	   it->firstValue() > offset+ 3*iBD.r.max) ++it;
    hsSource->setInactive(it);
    while (it != end && 
	   it->firstValue() >= offset + 3*iBD.r.min) {
	// insert a new node, by_first to the beginning of the list,
	// by_last at the current position
	createNode(b, offset, val, begin(b), it)->add(*it);
	++it;
    }
}

