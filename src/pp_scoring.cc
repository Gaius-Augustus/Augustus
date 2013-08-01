/**********************************************************************
 * file:    pp_scoring.cc
 * license: Artistic License, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  protein pattern search extension
 * authors: Oliver Keller
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 20.10.08| Oliver Keller | creation of the file
 * 28.01.09| Oliver Keller | rev. 157, final forward-only version
 **********************************************************************/

#include "pp_scoring.hh"

// project includes
#include "properties.hh"
#ifdef DEBUG
#include "statemodel.hh"  // for getStatePair
#endif

// standard C/C++ includes
#include <stdexcept>  // for out_of_range

using namespace PP;


#ifdef DEBUG
Matrix<Double>* ExonScorer::transitions = 0;
int globalThreshUsedCount = 0;
#endif

vector<bool> SubstateModel::stateStrands;
double SubstateModel::fullIntronBonus = 100;
double SubstateModel::minIntronFactor = 0.1;


/*--- class SubstateModel ------------------------------------------*/

SubstateModel::SubstateModel(Profile& p) : 
    prfl(p), 
    allow_truncated(false),
    exhaustive_substates(false)
{
    Properties::assignProperty("/ProteinModel/allow_truncated", allow_truncated);
    Properties::assignProperty("/ProteinModel/exhaustive_substates", exhaustive_substates);

    // initialize collection of hit sequences
    initHitSequences();

    initScores();
    
    maxBlockSize = prfl.blockSize(0);
    for (int blockno = 1; blockno < prfl.blockCount(); blockno++)
	if (prfl.blockSize(blockno) > maxBlockSize) 
	    maxBlockSize = prfl.blockSize(blockno);
}


void SubstateModel::advanceScores(int base) {
    static int base_lookahead[2] = { Constant::ass_end, Constant::dss_start };
    for (int b=0; b < blockCount(); b++) {
	for (int strand=0; strand<=1; strand++) {
	    map<int, Double> newhits;
	    // Constant::base_lookahead;
	    blockScores[strand][b].addBlocksUntil(strand, base + base_lookahead[strand], &newhits); 
	    blockScores[strand][b].removeBlocksUntil(base - base_lookahead[1-strand]);
	    for (map<int,Double>::const_iterator it = newhits.begin();
		 it != newhits.end();
		 ++it) {
		const int& offset = it->first;
		const Double& value = it->second;
		hitSeqColl[strand][offset%3].addHit(b, offset, value,
						    b==0 ? DistanceType(0) :  // ignored 
						    (interBlockDist(b, strand) +
						     Range(blockSize(b-1, strand))),
						    blockSize(b, strand));
	    }
  	}
    }
}

void SubstateModel::scoreAssOrRDss(ViterbiSubmapType& submap, bool complement,
				   int firstBase, const Double& mainProb) {
    ViterbiSubmapType::iterator it = submap.begin(0);
    while (it != submap.end() && it->first.slot < MAX_BLOCKCOUNT) {
 	// int id = it->first;
	Position ppos = it->first;
	Double prob = it->second.p; // * submap.factor;
	if (prob * submap.factor <= mainProb * prfl.getGlobalThresh(complement, ppos)) {
	    ViterbiSubmapType::iterator temp = it;
	    submap.erase(temp, ++it);
	    continue;
	}

#ifdef DEBUG
	if (ppos.i > blockSize(ppos.b, complement) || ppos.i < -interBlockDist(ppos.b, complement).r.max)
	    throw ProjectError("Error in SubstateModel::scoreAss: ppos.i out of range (ppos=(b:" + 
			       itoa(ppos.b) + ", i:" + itoa(ppos.i) +"), strand=" + (complement ? "-" : "+") +")");
#endif

	// if the entry is in a profile slot and the value refers to a block, 
	// create copies containing the suffix factor
	if (ppos.i > 0 && ppos.i < blockSize(ppos.b, complement))  {
	    // blck.addBlocksUntil(firstBase);
	    BlockScoreType& score = blockScores[complement][ppos.b];
#ifdef DEBUG
	    // if (score.addBlocksUntil(complement, firstBase - 3))
	    if (score.end() <=  firstBase - 3 && score.end() <= DNA::len - 3)
		throw ProjectError("SubstateModel::scoreAss: need to call addBlocksUntil!");
#endif
	    try { // BlockScoreType::get throws if firstBase is >= than Block::dnalen
		Double suffixFactor =  score.get(firstBase, ppos.i);
		if (suffixFactor > score.suffixThresh(complement, ppos.i))
		    // store a copy of the viterbi variable containing the suffix factor in the submap
		    // that has a temporary lifetime 
		    // (see [Multi|Single]TargetExonScorer::scoreInternal)
		    submap.push_back(ppos.copyId(), prob * suffixFactor, 0);
	    } catch (out_of_range) {
#ifdef DEBUG
		if (firstBase <= DNA::len-3)
		    throw;
#endif
	    }
	} else
	    submap.push_back(ppos.copyId(), prob, 0);
#ifdef DEBUG
	if (prob * submap.factor > 1) {
	    cout << "Found prob>1." << endl
		 << "factor=" << submap.factor << endl
		 << "beginOfFirstCodon=" << firstBase << endl
		 << "posn=(" << ppos.b << "/" << ppos.i << ")" << endl;
	    throw ProjectError("Error in SubstateModel::scoreAss: prob>1");
	}
#endif
	++it;
    }
#ifdef DEBUG
    if (submap.linkcount() > 1)
	throw ProjectError("ass / reverse dss has linkcount > 1");
#endif
}

void SubstateModel::clearLowScoring(ViterbiSubmapType& submap) {
    if (submap.linkcount() > 1 || submap.empty()) 
	return;

    bool complement = stateStrands[submap.getMainState()];
    bool hasCopies = (submap.begin(MAX_BLOCKCOUNT) != submap.end());
    
    for (int b=0; b<=blockCount(); b++) {
	DistanceType dist = interBlockDist(b, complement);
	ViterbiSubmapType::iterator 
	    begin = submap.begin(b), // bound(SubstateId(b, -dist.r.max-1)),
	    end = submap.bound(b,1), left, right, cand;
#ifdef DEBUG
	if (end != submap.end() && end->first.slot < MAX_BLOCKCOUNT) {
	    Position endpos = end->first;
	    if ((endpos.b == b && endpos.i <= 0) || endpos.b < b)
		throw ProjectError("We end in the wrong block!");
	}
#endif
	if (begin == end) continue;
#ifdef DEBUG
	if (Position(begin->first).b != b) 
	    throw ProjectError("We begin in the wrong block!");
	if (begin != submap.begin()) {
	    cand = begin; cand--;
	    if (Position(cand->first).b >= 0 && Position(cand->first).b == b)
		throw ProjectError("We begin inside the block!");
	}
#endif

	// RELAXATION PHASE (will not remove a lot of substates)
	//          (clear all states closer than RELAX * dist.max
	//           to a better one)
	if (!exhaustive_substates) {
	    left = right = begin;
	    right++;
	    double res = dist.r.max * RELAXATION;
	    while (right != end) {
		if (right->first.value - left->first.value < res) {
		    if (right->second < left->second) {
			cand = right++;
			if (hasCopies)
			    submap.erase(Position(cand->first).copyId());
			submap.eraseAndDecCounts(cand);
			continue;
		    } else {
			if (hasCopies)
			    submap.erase(Position(left->first).copyId());
			if (left==begin) {
			    begin=right;
			    end = submap.bound(b,1); // for the rare case that we
			                             // could just have erased end
			}
			if (left==begin)
			    begin=right;
			
			submap.eraseAndDecCounts(left);
		    }
		}
		left = right++;
	    }
	}
	    
	// PHASE 2: elimination of unneeded states 
	//          (enclosed by two better ones close enough to each other)
	if (dist.has_max) {
	    int width = dist.r.size();
	    SubstateId lowId(b, -width);
	    // width = width >> skipBits[0][b];
	    cand = left = begin;
	    if (++cand == end || width < 3) continue;
	    
	    right = cand; ++right;
	    while (right != end) {
		if (right->first.value - left->first.value <= width &&
		    cand->second <= left->second &&
		    cand->second <= right->second) 
		{
		    if (hasCopies)
			submap.erase(Position(cand->first).copyId());
		    submap.eraseAndDecCounts(cand);
		    if (left->second <= right->second && left != begin) {
			cand = left--;
			continue;
		    } 
		} else
		    left = cand;
		cand = right++;
	    }
	    while (cand != begin &&
		   left->first > lowId && 
		   left->second >= cand->second) {
		if (hasCopies)
		    submap.erase(Position(cand->first).copyId());
		submap.eraseAndDecCounts(cand);
		cand = left--;
	    }
	} else {
	    cand = end;
	    if (--cand == begin) continue;
	    ViterbiSubmapEntry checkVal = cand->second;
	    bool cont;
	    do {
		cont = (--cand != begin);
		if (cand->second <= checkVal) {
		    if (hasCopies)
			submap.erase(Position(cand->first).copyId());
		    submap.eraseAndDecCounts(cand++);
		} else {
		    checkVal = cand->second;
		}
	    } while (cont);
	}
    }
}


/*--- class BonusProbType -------------------------------------------*/

void MultiTargetExonScorer::initBonusProbs(BonusProbType& target, int b, int frame) {
    const BlockScoreType& blscore = blockScores[b];
    for (int l=1-endStateOffset; l<blscore.size(); l++) {
	try {
	    Double newfactor = blscore.get(endOfLastCodon - 3*l +1, 0, l) 
		* mdl.getIntronBonus(b, l, frame, complement);
	    if (l == 0 || newfactor > blscore.prefixThresh(complement, l)) 
		target.push_back(l+endStateOffset, newfactor);
	} catch (out_of_range) {
#ifdef DEBUG
	    if (endOfLastCodon >= blscore.end() && endOfLastCodon <= DNA::len-3 && endOfLastCodon >= 3*l - 1)
		cerr << "BonusProbType::BonusProbType: exception out_of_range caught.\n";
#endif
	}
    }
}

void BonusProbType::update(const Double& p, int begin, int end,
			   ViterbiSubmapType* source, SubstateId leftId) {
#ifdef DEBUG
    if (begin<0)
	throw ProjectError("PP::BonusProbType::update(): called with begin<0!");
#endif
    if (begin >= indexTable.size() ||  // no BonusProbEntry having blockpos>=begin exists
	end <= begin)
	return;    
    int a = indexTable[begin];     // a = first entry with i >= begin
    int b = end < indexTable.size() ? indexTable[end]-1 : indexTable.back();

    for (int j=a; j<=b; j++) {
	if (value(j) < p) v[j].set(p, source, leftId);
    }
}
    
void BonusProbType::mergeIntoSubmap(int blockindex,
				    ViterbiSubmapType& target,
				    SubstateIdMap* leftIds) {
    SubstateId substate(blockindex, 0);
    ViterbiSubmapType::iterator it = target.begin(0);
    for (int j=0; j<v.size(); j++) {
	substate.value = v[j].blockpos;
	Double p = v[j].fullProb() / target.factor;
	SubstateId& leftId = v[j].argmaxId;
	ViterbiSubmapType*& pred = v[j].argmaxMap;
#ifdef DEBUG
	if ( (pred==0) != (leftIds==0 || p==0) || (pred==0) != (leftId.undef()) )
	    throw ProjectError("Predecessor existence mismatch");
#endif
	if (p <= 0) continue;
	while (true) {
	    if (it == target.end() || it->first > substate) {
		target.insert(it, substate, p, pred);
		if (leftIds)
		    (*leftIds)[substate] = leftId;
		break;
	    }
	    if (it->first == substate) {
		if (it->second < p) {
		    it->second.set(p, pred);
		    if (leftIds)
			(*leftIds)[substate] = leftId;
		}
		++it; break;
	    }
	    ++it;
	}
    }
}




/*--- ExonScorer classes --------------------------------------------*/

/*--- class MultiTargetExonScorer -----------------------------------*/

MultiTargetExonScorer::MultiTargetExonScorer(SubstateModel& model,
					     ViterbiSubmapType& substates,
					     StateType type, 
					     int endOfBioExon) :
    ExonScorer(model, endOfBioExon, type), 
    targetVit(substates),
    source(0),
    leftId(-1),
    bestLeftId(0)
{
    for (int b=0; b<blockCount(); b++) {
#ifdef DEBUG
	if (endOfLastCodon >= 2 && blockScores[b].addBlocksUntil(complement, endOfLastCodon-2)) 
	    throw ProjectError("MultiTargetExonScorer: need to call addBlocksUntil!");
#endif
	bonusProbs.push_back(BonusProbType());
	initBonusProbs(bonusProbs.back(), b, endOfBioExon - endOfLastCodon);
    }
    if (!isFirstExon(type)) {
	bestLeftId = targetVit.getPredSubstateMap();
#ifdef DEBUG
	if (bestLeftId == 0)
	    throw ProjectError("bestLeftId == 0");
#endif
    }
    currentHSColl->init();
}

void MultiTargetExonScorer::updateBonusProbs(Position targetPos, const Double& bonusProb) {
    short& prefixScoreBegin = targetPos.i;
    int prefixScoreEnd;
    if (interBlockDist(targetPos.b).has_max) {
	// prefixScoreBegin corresponds to blockStart.max
	// rightmost position of block is the minimal
	// block coordinate at exon end
	prefixScoreEnd = prefixScoreBegin + interBlockRange(targetPos.b).size();
   } else { 
	// prefixScoreBegin is the maximal position
	prefixScoreEnd = prefixScoreBegin + 1;
	if (prefixScoreBegin >= 1)
	    prefixScoreBegin = 0;
    }
    if (prefixScoreEnd > aa_count + endStateOffset) 
	prefixScoreEnd = aa_count + endStateOffset +1;
	
    if (prefixScoreBegin < 1) {
#ifdef DEBUG
	if (prefixScoreBegin < -interBlockRange(targetPos.b).max)
	    throw ProjectError("Error in MTE::updateBonusProbs: updating index out of range");
#endif
	SubstateId substate = targetPos.id();
	if (targetVit.setMax(substate, bonusProb, source) && bestLeftId)
	    (*bestLeftId)[substate] = leftId;
	prefixScoreBegin = 1;
    } 
    if (targetPos.b < blockCount()) {
	bonusProbs[targetPos.b].update(bonusProb, prefixScoreBegin, prefixScoreEnd, source, leftId);
    }
}

bool MultiTargetExonScorer::scoreSourceEntry(int b, int maxBlockStart, Double bonusProb) {
    // returns false if next reachable block exceeds the block start range
    // if there is a block in range, update the bonus probs

    Position targetPos(b, aa_count + endStateOffset - maxBlockStart);
    bool has_max = interBlockDist(b).has_max;
    

    // 1) try to score without block inside
    bool foundTarget = !(has_max && exceedsBlock(targetPos));
    if (foundTarget) 
	updateBonusProbs(targetPos, bonusProb);
    if (b == blockCount()) 
	// we already reached the last block
	return foundTarget;

    // 2) try all the hit sequences that might fit
    HitSequenceList& hsList = currentHSColl->list(b);
    HitSequenceList::reverse_iterator it = hsList.rbeginActive();
    int minBlockStart = interBlockRange(b).alignRight(maxBlockStart).min;
    while (it != hsList.rend()) {
	// iterate over block sequences, in upstream direction
	int dist = (it->firstValue() - beginOfFirstCodon) / 3;
	if (has_max && dist > maxBlockStart) {
	    // firstValue() too far downstream
	    // we discard the block only if it is too far downstream
	    // for any position of the previous block
	    ++it;
	    if (dist > interBlockFullRange(b).max)
		hsList.setInactive(it); // dont check this hit sequence again 
	                                // in this iteration
	} else if (dist < 0 || dist < minBlockStart) 
	    // firstValue() too far upstream: we are done with the hit sequence
	    // but don't delete yet
	    return true;
	else {
	    Position ppos = it->targetPosAt(endOfLastCodon);
	    if (ppos.i >= blockSize(ppos.b)) {
		// ensure that last block ends before exon
		// this is needed as hit sequences are prematurely updated
		ppos = asNextPos(ppos) + endStateOffset;
		if (interBlockDist(ppos.b).has_max && ppos.i > blockSize(ppos.b)) {
		    // lastValue() too far upstream: this hit sequence
		    // will not be needed anymore; we erase it only
		    // if it's not the last block (rinitial/terminal might still
		    // need it, so SingleTargetExonScorer::bestFittingBonusFactor
		    // takes care of its deletion)
		    if (ppos.b < blockCount())
			it.erase();
		    else ++it;
		    continue;
		}
		// block is in range 
		foundTarget = true;
		updateBonusProbs(ppos, bonusProb * it->getFactor());
	    }
	    ++it;
	} 
    }
    // all blocks were downstream of blockStart.min;
    // return false if they were downstream of blockStart.max too
    return foundTarget;
} // MultiTargetExonScorer::scoreSourceEntry

void MultiTargetExonScorer::scoreFirst(ViterbiColumnType& col, int predState, const Double& transEmiProb) {
    if (is_active && !scoreSourceEntry(0, interBlockRange(0).max, col[predState] * transEmiProb))
	is_active = false;
}

void MultiTargetExonScorer::scoreInternal(ViterbiColumnType& col, int predState, const Double& transEmiProb) {
    // here, we evaluate temporary copies of the substate map 
    // the advantage is that we can delete them so there is less to iterate
    // as viterbi/forward moves on
    // the copies have been produced by SubstateModel::scoreAss
    source = &col.getSubstates(predState);
    ViterbiSubmapType::iterator it = source->begin(MAX_BLOCKCOUNT), old_it, end = source->end();
    while (it != end) {
	Position leftPos = it->first;
	leftId = leftPos.id();
	int maxBlockStart;
	int b = leftPos.b;
#ifdef DEBUG
	if (leftPos.i < -interBlockRange(b).max || leftPos.i > blockSize(b))
	    throw ProjectError("Viterbi Submap contains invalid left pos i=" + itoa(leftPos.i) +
			       ".(predState=" + itoa(predState) + ", codon range=[" + itoa(beginOfFirstCodon) +
			       ".." + itoa(endOfLastCodon) + "]). ");
#endif
	if (leftPos.i <= 0) { //
	    maxBlockStart = -leftPos.i;
	} else {
	    int remaining = blockSize(b)-leftPos.i;
	    if (aa_count < remaining) {
		// only relevant if there are multiple introns in block; skip it here
		++it; continue;
	    }
	    maxBlockStart = interBlockRange(++b).max + remaining;
	}
	if (!scoreSourceEntry(b, maxBlockStart, it->second.p * source->factor * transEmiProb)) {
	    old_it = it; 
	    if (b == blockCount()) { // we are done; do not delete 
                it=end; break;       // because terminal/rinitial might 
	    }                        // still need the viterbi variable
	    SubstateId maxId(b + MAX_BLOCKCOUNT,0);
	    do {} while (++it != end && it->first <= maxId);
	    source->erase(old_it, it);
	} else
	    ++it;
    }

    // now, check for short exons with both adjacent introns inside the block
    end = source->begin(MAX_BLOCKCOUNT);
    if (aa_count >= mdl.maxBlockSize || aa_count < MIN_CHECKCOUNT) 
	return;
    for (int b=0; b<blockCount(); b++) {
	int exceeding = blockSize(b)-aa_count;
	if (exceeding <= 1)
	    continue;
	for (it = source->bound(b,1); it != end; ++it) {
	    Position ppos = it->first;
	    if (ppos.b > b || ppos.i >= exceeding)
		break;
	    Double score = blockScores[b].get(beginOfFirstCodon, ppos.i, aa_count);
#ifdef DEBUG
	    if (score <= Block::almostZero) 
		throw ProjectError("Stopcodon in Exon?!");
#endif
	    if (score < getBlock(b).getPartialThresh(complement, ppos.i, ppos.i + aa_count))
		continue;
	    Position rightPos = ppos + aa_count + endStateOffset;
	    if (targetVit.setMax(rightPos.id(), it->second.p * source->factor * transEmiProb * score, source))
		(*bestLeftId)[rightPos.id()] = ppos.id();
	}
    }
} // MultiTargetExonScorer::scoreInternal

void MultiTargetExonScorer::postProcessing(const Double& maxProb) {
#ifdef DEBUG
    if (targetVit.inactive())
	throw ProjectError("Internal Error in ExonScorer::postProcessing");
#endif 

    // merge in bonus probs (i.e., in-block substates)
    for (int b=0; b<blockCount(); b++) 
	bonusProbs[b].mergeIntoSubmap(b, targetVit, bestLeftId);

    // ViterbiSubmapType& addedVit = col.addSubstates(targetVit);

    // apply global thresholds
    ViterbiSubmapType::iterator it = targetVit.begin(0), end = targetVit.begin(MAX_BLOCKCOUNT), temp;
    while (it != end)
	if (it->second.p > maxProb * mdl.prfl.getGlobalThresh(complement, it->first)) 
	    ++it;
	else {
	    temp = it;
	    targetVit.erase(temp, ++it);
	} 
}

/*--- class SingleTargetExonScorer ----------------------------------*/

// constructor for viterbi/forward on terminal/rinitial/(r)single exons
SingleTargetExonScorer::SingleTargetExonScorer(SubstateModel& model,
					       StateType type,
					       int endOfBioExon) :
    ExonScorer(model, endOfBioExon, type), 
    rightPos(model.blockCount(),0),
    predSubstate(-1),
    maxProb(0),
    newmax(false),
    bestSource(0),
    bestHitSeq(0),
    prefixFactor(1)
{
    calcAllowedDist(isLastExon(type));
    currentHSColl->init();
}

// constructor for backward mode
SingleTargetExonScorer::SingleTargetExonScorer(SubstateModel& model, 
					       StateType type, SubstateId substate,
					       int endOfBioExon, int ORFleft) :
    ExonScorer(model, endOfBioExon, type),
    rightPos(substate),
    predSubstate(-1),
    maxProb(0),
    newmax(false),
    bestSource(0),
    bestHitSeq(0),
    prefixFactor(1)

{
#ifdef DEBUG
    currentHSColl = new HitSequenceCollection(blockCount());
#endif
    backward_mode = true;
    if (substate.undef()) // terminal, single, or rinitial exon
	rightPos = Position(model.blockCount(),0);
    else if (rightPos.i > endStateOffset && rightPos.i>0) {
	int i = rightPos.i - endStateOffset;
	const Block& blk = model.getBlock(rightPos.b, complement);
	prefixFactor = blk.scoreFromScratch(complement, endOfLastCodon - 3*i + 1, 0, i);
	if (prefixFactor < blk.getPrefixThresh(complement, i)) {
	    is_active = false;
	    return;
	}
    }
    // moved to SingleTargetExonScorer::addMatches
    // model.addPrefixMatch(rightId() - endStateOffset, endOfLastCodon, endOfBioExon, complement);
    calcAllowedDist(isLastExon(type));
    createHitSequences(ORFleft);
}

void SingleTargetExonScorer::createHitSequences(int ORFleft) {
    // initialize block sequence sets with hits, save ranges for right pos
    HitSequenceCollection& hsColl = *currentHSColl;
    hsColl.reset();
    
    if (rightPos.b == 0)
	return;
    
    // initialize block sequence set for the case that right pos is fixed
    // and all others ignored (backtracking, etype==terminal/rinitial);
    DistanceType distToNext = allowedDist;
    distToNext.r +=  interBlockFullRange(rightPos.b).max - endStateOffset;
    int dnamin = distToNext.has_max ? endOfLastCodon + 1 - 3*distToNext.r.max : ORFleft;
    int b = rightPos.b -1;
    // ensure that block b ends before exon
    int dnamax = endOfLastCodon + 1 - 3*max(distToNext.r.min, blockSize(b));
   
    if (dnamax >= dnamin) 
	while (true) {
	    // we might remove the following later
	    const map<int,Double>& hits = blockScores[b].savedHits(dnamax % 3);
	    map<int,Double>::const_reverse_iterator it, 
		it_end(hits.lower_bound(dnamin)), 
		it_begin(hits.upper_bound(dnamax));
	    if (it_begin == it_end) 
		break; 
	    for (it = it_begin; 
		 it != it_end;
		 ++it) {
		int dnapos = it->first;
		if (b == rightPos.b -1) {
		    hsColl.addSingleFront(b, dnapos, it->second);
		}
		else 
		    hsColl.addHitFront(b, dnapos, it->second, distToNext);
	    }
	    if (b==0)
		break;
	    distToNext = interBlockFullDist(b--);
	    dnamax = it_begin->first - 3*distToNext.r.min;
	    if (ORFleft > dnamax)
		break;
	    if (distToNext.has_max) {
		dnamin = it_end.base()->first - 3*distToNext.r.max;
		if (ORFleft > dnamin) 
		    dnamin = ORFleft;
	    } else 
		dnamin = ORFleft;
#ifdef DEBUG
	    if (dnamin > dnamax)
		throw ProjectError("Fatal error in SingleTargetExonScorer:createHitSequences");
#endif
	}
    // activeMax was used temporarily; is now initialized again
    hsColl.init();
}	    

void SingleTargetExonScorer::scoreFirst(ViterbiColumnType& col, int predState, const Double& transEmiProb) {
    DistanceType firstDist = interBlockDist(0);
    predSubstate.set(0);
    Double predProb = col[predState] * transEmiProb;
    const HitSequence* argmax = 0;
    if (rightPos.b == 0) {
	// range for block start relative to exon start
	newmax = false;
	if (!allowedDist.has(- firstDist.r.max + aa_count +endStateOffset ))
	    return;
    } else 
	predProb *= bestFittingBonusFactor(0, firstDist.r.max, argmax);
    newmax = predProb > maxProb;
    if (newmax) {
	bestHitSeq = argmax;
	maxProb = predProb;
    }
}

void SingleTargetExonScorer::scoreInternal(ViterbiColumnType& col, int predState, const Double& transEmiProb) {
    ViterbiSubmapType& source = col.getSubstates(predState);
    Double bonusProb = 0;
    const HitSequence* argmax = 0;
    newmax = false;
    int len = rightPos.i - endStateOffset;
    if (backward_mode && len > aa_count && rightPos.i > 0) {
	const Block& blk = getBlock(rightPos.b);
	int prefixLen = len - aa_count;
	SubstateId leftId(rightPos.b, prefixLen);
	Double prob = source.get(leftId);
	if (prob == 0)
	    return;
	Double externalScore = blk.scoreFromScratch(complement, beginOfFirstCodon - 3*prefixLen, 0, leftId.value);
#ifdef DEBUG
	if (blk.scoreFromScratch(complement, beginOfFirstCodon, prefixLen, aa_count) <= Block::almostZero) 
	    throw ProjectError("Stop codon in Exon?!");
#endif
	if (prefixFactor / externalScore < blk.getPartialThresh(complement, prefixLen, len))
	    return;
	bonusProb = (prob * transEmiProb) / externalScore;
	if (bonusProb > maxProb) {
	    predSubstate = leftId;
	    newmax = true;
	    bestHitSeq = argmax;
	    bestSource = &source;
	    maxProb = bonusProb;
	}
	return;
    }
    if (!is_active)
	return;
 
    ViterbiSubmapType::iterator it, end;
    if (backward_mode) { // take original viterbi variables
	it = source.begin(0); end = source.bound((rightPos+1).id());
    } else { // in viterbi/sampling, take temporary copies produced by SubstateModel::scoreAss
	it = source.begin(MAX_BLOCKCOUNT); end = source.end();
    }

    while (it != end) {
	Position leftPos = it->first;
	int b = leftPos.b;
#ifdef DEBUG
	if (leftPos.i < -interBlockRange(b).max || leftPos.i > blockSize(b))
	    throw ProjectError("Viterbi Submap contains invalid left pos i=" + itoa(leftPos.i) +
			       ".(predState=" + itoa(predState) + ", codon range=[" + itoa(beginOfFirstCodon) +
			       ".." + itoa(endOfLastCodon) + "]). ");
#endif
	int maxBlockStart;
	if (leftPos.i <= 0)  
	    maxBlockStart = -leftPos.i;
	else {
	    int remaining = blockSize(b)-leftPos.i;
	    if (aa_count < remaining) {
		// can happen only if there are multiple introns in block
		// this is already dealt with above
		++it; continue; 
	    }
	    maxBlockStart = interBlockRange(++b).max + remaining;
	}
	Double leftScore =  backward_mode && leftPos.i>0 ? 
	    it->second.p * getBlock(leftPos.b).checkedSuffixScore(complement, beginOfFirstCodon, leftPos.i) :
	    it->second.p;  // values of temporary copies already contain the score for the block suffix
	leftScore *= source.factor;

	if (b >= rightPos.b) {
	    // no full blocks in exon
	    if (b > rightPos.b) {
#ifdef DEBUG
		// this should be excluded from above
		if (b > rightPos.b + 1 || leftPos.i <= 0 || leftPos.i > rightPos.i)
		    throw ProjectError("Internal Error in SingleTargetExonScorer::scoreInternal: b=" + itoa(b));
#endif
		break; // no more positions to come
	    }
	    // block start relative to exon end
	    int candidateI =  -maxBlockStart + aa_count + endStateOffset;
	    if (candidateI < allowedDist.r.min) {
		// previous block too close
		++it; continue;
	    } else if (allowedDist.has_max && candidateI > allowedDist.r.max) {
		break; // previous block too far
	    }
	    argmax = 0;
	    bonusProb = leftScore * transEmiProb;
	} else {
	    bonusProb = leftScore * transEmiProb * bestFittingBonusFactor(b, maxBlockStart, argmax);
	}
	if (bonusProb > maxProb) {
	    predSubstate = leftPos.id();
	    newmax = true;
	    bestHitSeq = argmax;
	    bestSource = &source;
	    maxProb = bonusProb;
	}
	++it;
    }
}

void SingleTargetExonScorer::postProcessing(const Double& mainProb) {
    if (maxProb <= mainProb) {
	maxProb = 0;
	bestSource = 0;
	predSubstate.set(0);
    }
}

void SingleTargetExonScorer::exportSubstates(ViterbiSubmapType& targetVit) {
    if (backward_mode) {
#ifdef DEBUG
	throw ProjectError("exportSubstates must not be called from backtracking!");
#endif
	return;
    }

    // change the main entry - moved to ExonModel
    // col[state] = this->maxProb;

    // add one entry into the newly submap
    if (maxProb > 0) {
	SubstateId rightId = rightPos.id();
	targetVit.activate();
	targetVit.setMax(rightId, this->maxProb, bestSource);
	// now, make the new predSubstateMap!!!
	if (rightId != predSubstate) 
	    (*targetVit.getPredSubstateMap())[rightId] = predSubstate;
    }

    // and increment the successor counts:
    // REMOVED THIS, because it is dealt with by registerPredecessor
    // if (bestSource)
    // 	bestSource->succCount(predSubstate)++;
}

Double SingleTargetExonScorer::bestFittingBonusFactor(int b, int maxBlockStart, const HitSequence*& argmax) {
    // returns value>0 if new argmax was found
    HitSequenceList& hsList = currentHSColl->list(b);
    if (hsList.empty())
	return 0;
    Double result = 0;
    argmax = 0;
    HitSequenceList::reverse_iterator it = hsList.rbeginActive();
    int minBlockStart = interBlockRange(b).alignRight(maxBlockStart).min;
    if (minBlockStart < 0)
	minBlockStart = 0;
    while (it != hsList.rend()) {
	// iterate over block sequences, in upstream direction
	if (!backward_mode) { 
            // some hit sequences are not compatible with rightPos
#ifdef DEBUG
	    if (rightPos.i != 0 || !is_last) // check if it is a terminal exon!
		throw ProjectError("SingleTargetExonScorer::bestFittingBonusFactor: "
				   "in forward mode, only terminal/rinitial exons are allowed here!");
#endif
	    Position ppos = asNextPos(it->targetPosAt(endOfLastCodon)) + endStateOffset;
	    if (ppos.b != rightPos.b || ppos.i < allowedDist.r.min) { // last block too close
		++it; continue;
	    } 
	    if (allowedDist.has_max && ppos.i > allowedDist.r.max) { 
                // last block too far upstream; delete the hit sequence
		it.erase(); continue;
	    }
	}
#ifdef DEBUG
	// ensure if all hit sequences are compatible with rightPos,
	// otherwise it's a bug
	else {
	    Position ppos = asNextPos(it->targetPosAt(endOfLastCodon)) + endStateOffset;
	    if (rightPos.i > 0 || is_last) {
		if (ppos.b != rightPos.b ||
		    ppos.i <= rightPos.i - interBlockDist(ppos.b).r.size() ||
		    (interBlockDist(ppos.b).has_max && ppos.i > rightPos.i))
		    throw ProjectError("SingleTargetExonScorer: "
				       "we found a hsList not suitable for rightPos\n");
	    } else {
		if (ppos.i > 0 && !interBlockDist(ppos.b).has_max)
		    ppos.i = 0;
		
		// 	ppos.i = mdl.getReferenceI(ppos.b, ppos.i);
		if (ppos != rightPos)
		    throw ProjectError("SingleTargetExonScorer: "
				       "we found a hsList not suitable for rightPos\n");
	    }
	}
#endif
	int dist = (it->firstValue() - beginOfFirstCodon ) / 3;
	Double newfactor = it->getFactor();
	if (interBlockDist(b).has_max && dist > maxBlockStart) {
	    ++it;
	    if (dist > interBlockFullRange(b).max)
		hsList.setInactive(it); // dont check this hit sequence again in this iteration
	} else {
	    if (dist < minBlockStart) 
		return result;
	    else if (newfactor > result) {
		result = newfactor;
		argmax = &(*it);
	    }
	    ++it;
	}
    }
    return result;
}


/*--- only used in DEBUGGING mode -----------------------------------*/

#ifdef DEBUG

Double SingleTargetExonScorer::getFactor() {
    return (bestHitSeq == 0) ? 1 : bestHitSeq->getFactor();
}

void SingleTargetExonScorer::checkResults(const OptionListItem& oli, 
					  const ViterbiMatrixType& viterbi, 
					  int state, int base,
					  Double endPartProb, 
					  Double notEndPartProb) {
    SubstateId substate, rightId = rightPos.id(); 
    if (!is_last) substate = rightId;
    int predState = -1;
    SubstateId predSubstate;
    string error = "";
    Double maxProbNow = 0;
    Double transEmiProb = 0;
    Double leftViterbi = 0;
    Double rightViterbi = viterbi[base].get(state,substate);
    Double suffixFactor = 1;
    Position ppos;
    if (oli.state<0)
	error = "backtracking hasnt found anything valid here!";
    else {
	const ViterbiColumnType& col = viterbi[oli.base<0?0:oli.base];
	getStatePair(oli.state, predState, predSubstate);
	transEmiProb = (*transitions)[predState][state] * endPartProb * notEndPartProb;
	leftViterbi = col.get(predState, predSubstate);
	const ViterbiSubmapType* pred = 0;
	const ViterbiSubmapType* currmap = 0;
	try {pred = &col.getSubstates(predState);} catch(...) {}
	try {currmap = &viterbi[base].getSubstates(state);} catch(...) {}
	ppos = predSubstate;
	newFirstCodon(oli.base - (complement ? Constant::dss_start : Constant::ass_end)+1);
	if (ppos.i > 0) {
	    const Block& blk = getBlock(ppos.b);
	    suffixFactor = ppos.b == rightPos.b ?
		1 / blk.scoreFromScratch(complement, beginOfFirstCodon - 3*ppos.i, 0, ppos.i) :
		blk.checkedSuffixScore(complement, beginOfFirstCodon, ppos.i);
	}
	maxProbNow = leftViterbi * suffixFactor * transEmiProb * getFactor();
	ViterbiSubmapType::const_iterator iter = currmap->bound(rightId);
	if (iter == currmap->end() || iter->first != rightId || iter->second.predMap != pred) {
	    error = "Predecessor submap mismatch!\n";
	} else if (pred != 0 && predSubstate != currmap->getPredecessorSubstate(rightId)) {
	    error = "Predecessor substate mismatch!\n";
	} else if (viterbi[base].get(state) < col[predState] * transEmiProb) {
	    error = "Something is wrong with original viterbi!\n";
	} else if (relative_error(maxProbNow, maxProb)>10e-10) {
	    error = "left Viterbi * suffixFactor * transEmi * hitseqFactor != maxProb";
	} else if (maxProb==0) {
	    error = "maxProb=0";
	} else if (relative_error(maxProb * getPrefixFactor(),rightViterbi)>10e-10) {
	    error = "maxProb * prefixFactor != right viterbi";
	}
#ifndef DEBUG_BACK
	if (!predSubstate.undef() && error=="") 
	    cerr << "substate bonus at base " << base << ": "
		 << viterbi[oli.base].get(predState, predSubstate) / viterbi[oli.base][predState]
		 << endl;
#endif
    }
    if (error != "") {
	cerr << "\n***\nviterbi[" << oli.base << "][" << predState << "](";
	cerr << ppos;
	cerr << ")=" << leftViterbi << "\n"
	     << "viterbi[" << base << "][" << state << "](" << Position(substate) 
	     << ")= " << rightViterbi << "\n";
	cerr << "viterbi[" << base << "][" << getFullStateId(state, substate) << "]\n";
	cerr << "maxProb=" << maxProb << "\n"
	     << "transEmiProb=" << transEmiProb << "\n" 
	     << "factor=" << getFactor() << "\n"
	     << "leftvit * transEmi * factor = " << maxProbNow << "\n"
	     << "prefixFactor=" << getPrefixFactor();
	cerr << "\n***\n";
	throw ProjectError("SingleTargetExonScorer::viterbiForwardAndSampling(backtracking):\n" + error);
    }
}
#endif
    
