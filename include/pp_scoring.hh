/*****************************************************************************\
 * Filename : pp_scoring.hh
 * Author   : Oliver Keller
 * Project  : Gene Prediction with Protein Family Patterns
 *
 * Description: Classes that perform the evaluation of ProteinProfile
 *              emission probabilites
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 20.10.08   | Oliver Keller         | creation of the file
 * 06.11.08   | Oliver Keller         | first version in SVN
 * 28.01.09   | Oliver Keller         | rev. 157, final forward-only version
 \******************************************************************************/

#ifndef __PP_SCORING_HH
#define __PP_SCORING_HH

// project includes
#include "vitmatrix.hh"
#include "pp_hitseq.hh"
#ifdef DEBUG
#include "matrix.hh"
extern int globalThreshUsedCount;
#endif


namespace PP {
    /*
     *  storing maximal bonusProbs for various right sides
     *  a bonusProb is the predProb (predecessor viterbi * transEmiProb)
     *  multiplied with the bonus factor of a hit sequence if any
     */

    struct BonusProbEntry {
	Double value;          // current maximal bonusProb for the blockpos
	Double prefixFactor;   // the factor for the block part [0..i]
	int blockpos;          // the corresponding i for the right side
	ViterbiSubmapType* argmaxMap;
	SubstateId argmaxId;
	BonusProbEntry(int i, Double factor) :
	    value(0), 
	    prefixFactor(factor), blockpos(i),
	    argmaxMap(0), argmaxId() {}
#ifdef DEBUG
	bool operator==(const BonusProbEntry& other) const {
	    return value==other.value 
		&& blockpos==other.blockpos 
		&& prefixFactor==other.prefixFactor;
	}
#endif
	Double fullProb() const {
	    return value * prefixFactor;
	}
	void set(const Double& p, ViterbiSubmapType* source, SubstateId leftPos) {
	    value = p;
	    argmaxMap = source;
	    argmaxId = leftPos;
	}
    }; // struct PP::BonusProbEntry


    struct BonusProbType {
	// BonusProbType(const BlockScoreType&, int base, int offset, bool complement);
#ifdef DEBUG
	string debug_info();
#endif
	void push_back(int i, const Double& factor) {
	    while (indexTable.size() <= i) 
		indexTable.push_back(v.size());
	    BonusProbEntry newval(i, factor);
	    v.push_back(newval);
	}
	void init() {
	}
	Double& value(int j) {
	    return v[j].value;
	}
	const Double& value(int j) const {
	    return v[j].value;
	}
	void update(const Double& p, int begin, int end, 
		    ViterbiSubmapType* source, SubstateId leftId);

	void mergeIntoSubmap(int blockindex,
			     ViterbiSubmapType& target,
			     SubstateIdMap* );

	// the BonusProbEntry, one for each valid right side
	// v[j].blockpos is the corresponding i
	vector<BonusProbEntry> v;

	// for each block pos i, j = indexTable[i] is:
	// - the index of the BonusProbEntry having blockpos=i, or
	// - the index of the first BonusProbEntry having blockpos>i.
	// if neither exists, indexTable has size <= i
	vector<int> indexTable;
#ifdef DEBUG
	bool operator==(const BonusProbType& other) const {
	    return v==other.v && indexTable==other.indexTable;
	}
#endif
    }; // struct PP::BonusProbType


    struct Match {
	Match (string id, int blno, int f, int l, Double d, int off, bool comp) : 
	    blockId(id), blockno(blno), firstBase(f), lastBase(l), score(d), 
	    offset(off), complement(comp) {} 
	// shift is implemented as a function _object_ 
	// since we need it in for_each
	struct shift {
	    shift(int d) : offset(d) {}
	    void operator() (Match& m) const {
		m.firstBase += offset;
		m.lastBase += offset;
	    }
	    int offset;
	};
	int firstCodon() const { // numbering starts with 0
	    return offset/3;
	}
	int lastCodon() const { // equals firstCodon of next match if codon is split
	    return (offset + lastBase - firstBase)/3;
	}
	int phase() const {
	    return mod3(-offset);
	}
	    
	string blockId;
	int blockno;
	int firstBase;
	int lastBase;
	Double score;
	int offset;
	bool complement;
    };


    struct SubstateModel {
	SubstateModel(Profile& p);
	void initScores() {
	    for (int strand=0; strand<2; strand++) {
		blockScores[strand].clear();
		for (int b=0; b<blockCount(); b++) 
		    blockScores[strand].push_back(BlockScoreType( getBlock(b,strand) ));
	    }
	}
	void initHitSequences() {
	    // initialize  hit sequence collections
	    for (int f=0; f<3; f++) {
		hitSeqColl[0][f].reset(blockCount());
		hitSeqColl[1][f].reset(blockCount());
	    }
	}
	void advanceScores(int base);
	void scoreAssOrRDss(ViterbiSubmapType& submap, bool complement, int firstBase, const Double& maxProb);
	void clearLowScoring(ViterbiSubmapType& submap);
	void clearLowScoringSubstates(ViterbiColumnType& vit) {
	    vit.foreachSubstate(this, &SubstateModel::clearLowScoring);
	}
	double getIntronBonus(int b, int i, int frame, bool comp) {
	    const Block& blk = getBlock(b, comp);
	    if (!blk.hasIP())
		return 1.0;
	    double freq;
	    if (comp) {
		if (frame)
		    freq = blk.getIntronFreq(blk.size()-i-1, 3-frame);
		else
		    freq = blk.getIntronFreq(blk.size()-i, 0);
	    } else
		freq = blk.getIntronFreq(i, frame);
	    freq *= fullIntronBonus;
	    return freq > minIntronFactor ? freq : minIntronFactor;
	}
	Profile& prfl;
	int maxBlockSize;
	vector<BlockScoreType> blockScores[2];  //index of blockScores is b, not blockNo!!!
	HitSequenceCollection hitSeqColl[2][3];
	list<Match> currentResult;
	bool allow_truncated;      // allow right-truncated block-scoring exons
	bool exhaustive_substates; // do not use a relaxation heuristic to clear substates
	                           // see SubstateModel::clearLowScoring

	/*---- methods evaluating positions in the profile ----*/ 

	// the global nucleotide coordinate of a position
	// static int getId(int b, int i) {
	//     return Position(b,i).id();
	// }
	// SubstateId toCopyId(SubstateId id) const {
	//     return SubstateId(id.slot + MAX_BLOCKCOUNT);
	// }
	// int fromCopyId(int id) const {
	//     return id + copyOffset;
	// }
	// int getCopyId(Position ppos) const {
	//     return toCopyId(ppos.id());
	// }
	// int getCopyId(int b, int i) const {
	//     return toCopyId(getId(b,i));
	// }

// 	Position getPosOfId(int id) const {
// 	    Position result = Position::fromId(id);
// #ifdef DEBUG
// 	    if (result.b > blockCount()) 
// 		cout << "getB() > #blocks\n";
// #endif
// 	    return result;
// 	}
// 	Position getPosOfCopyId(int id) const {
// 	    return getPosOfId(fromCopyId(id));
// 	}

#ifdef DEBUG
	void checkInRange(Position pos, bool complement) {
	    if (pos.b > blockCount() || pos.b < 0 || 
		(pos.i>0 && pos.i > blockSize(pos.b, complement)))
		throw ProjectError("block position out of range");
	}
#endif 
	int blockCount() const {
	    return prfl.blockCount();
	}
	int blockNoOfB(int b, bool comp) const {
	    if (b == blockCount())
		throw ProjectError("Invalid block no. in SubstateModel::blockNoOfB");
	    return comp ? blockCount() - 1 - b : b;
	}
	int iBlockOfB(int b, bool comp) const {
	    return comp ? blockCount() - b : b;
	}
	DistanceType interBlockDist(int b, bool comp) const {
	    return prfl.interBlockDist(iBlockOfB(b, comp));
	}
	int blockSize(int b, bool comp) const {
	    return b == blockCount() ? 0 : prfl.blockSize(blockNoOfB(b, comp));
	}
	const Block& getBlock(int b, bool comp) const {
	    return prfl[blockNoOfB(b, comp)];
	}

	void addHitSeqMatch(const HitSequence* hs, bool comp) {
	    if (hs)
		for (int b = hs->last(); b >= hs->first(); b--) {
		    int blockno = blockNoOfB(b, comp);
		    int firstBase = (*hs)[b];
		    int lastBase = firstBase + prfl.blockSize(blockno)*3 - 1;
		    Double score = blockScores[comp][b].savedScore(firstBase).getRoot(prfl.blockSize(blockno));
		    currentResult.push_front(Match(prfl[blockno].id, blockno, firstBase, lastBase, score, 0, comp));
		}
	}
	void addPrefixMatch(Position ppos, 
			    int endOfLastCodon, int endOfBioExon, 
			    bool comp) {
 	    int blockno = blockNoOfB(ppos.b, comp);
#ifdef DEBUG
	    checkInRange(ppos, comp);
#endif
	    int startOfMatch = endOfLastCodon - 3*ppos.i + 1;
	    if (startOfMatch <= endOfBioExon) {
		Double score = ppos.i==0 ? 1 : 
		    prfl[blockno].scoreFromScratch(comp, startOfMatch, 0, ppos.i).getRoot(ppos.i);
		int offset = comp ? 
		    3*(prfl.blockSize(blockno) - ppos.i)-endOfBioExon + endOfLastCodon : 0;
		currentResult.push_front(Match(prfl[blockno].id, blockno, 
					       startOfMatch, endOfBioExon, 
					       score, offset, comp));
	    }
#ifdef DEBUG
	    else throw ProjectError("Called addPrefixMatch with invalid value!");
#endif
	}

	void addSuffixMatch(Position ppos, 
			    int beginOfFirstCodon, int beginOfBioExon, 
			    bool comp) {
	    int blockno = blockNoOfB(ppos.b, comp);
#ifdef DEBUG
	    checkInRange(ppos, comp);
#endif
	    int len = prfl.blockSize(blockno) - ppos.i;
	    int endOfMatch = beginOfFirstCodon + 3*len - 1;
	    if (beginOfBioExon <= endOfMatch) {
		Double score = len==0 ? 1 : 
		    prfl[blockno].scoreFromScratch(comp, beginOfFirstCodon, ppos.i).getRoot(len);
		int offset = comp ? 0 : 
		    3*ppos.i-(beginOfFirstCodon-beginOfBioExon);
		currentResult.push_front(Match(prfl[blockno].id, blockno, 
					       beginOfBioExon, endOfMatch, 
					       score, offset, comp));
	    }
#ifdef DEBUG
	    else throw ProjectError("Called addSuffixMatch with invalid value!");
#endif
	}

	void addInternMatch(Position ppos, 
			    int beginOfBioExon, int endOfLastCodon, int endOfBioExon, 
			    bool comp) {
	    int blockno = blockNoOfB(ppos.b, comp);
	    int len = (endOfLastCodon - beginOfBioExon)/3;
	    Double score = len == 0 ? 1 :
		prfl[blockno].scoreFromScratch(comp, endOfLastCodon-3*len+1, ppos.i-len, len).getRoot(len);
	    int offset = comp ? 
		3*(prfl.blockSize(blockno) - ppos.i)-endOfBioExon + endOfLastCodon : 
		3*ppos.i - endOfLastCodon + beginOfBioExon;
	    currentResult.push_front(Match(prfl[blockno].id, blockno, 
					   beginOfBioExon, endOfBioExon, 
					   score, offset, comp));
	}

	void appendMatchesTo(list<Match>& target) {
	    target.splice(target.end(), currentResult);
	}

#ifdef DEBUG
	void outputFullBlocks(int offset=0) {
	    multimap< int, string > output;
	    for (int b=0; b<prfl.blockCount(); b++) 
		for (int f=0; f<6; f++) {
		    bool complement = f%2;
		    const map<int, Double>& hits = blockScores[complement][b].savedHits(f%3);
		    int blockno = blockNoOfB(b, complement);
		    for (map<int, Double>::const_iterator it = hits.begin();
			 it != hits.end(); ++it) {
			ostringstream strm;
			strm << it->first + offset << "\t" 
			     << it->first + offset + 3*prfl.blockSize(blockno) << "\t"
			     << it->second << "\t"
			     << (complement ? '-' : '+') << "\t0\t"
			     << prfl[blockno].id << "\n";
			output.insert(make_pair(it->first, strm.str()));
		    }
		}
	    multimap<int, string>::iterator it;
	    for (it = output.begin(); it != output.end(); ++it) 
		cerr << "\tAUGUSTUS\tblock_hit\t" << it->second ;
	}
	void printGlobalThresh(Double p) {
	    for (int b=0; b<blockCount(); b++) {
		Double p1 = prfl.getGlobalThresh(false, Position(b,0));
		Double p2 = prfl.getGlobalThresh(true, Position(b,0));
		cerr << "(thresh " << b << ") (+):" << p * p1 << " (-):" << p * p2 << "\n";
	    }
	    cerr << "( backward strand states are:";
	    for (int i=0; i<stateStrands.size(); i++) 
		if (stateStrands[i])
		    cerr << " " << i;
	}
#endif
	static void setStateStrands(const vector<StateType>& stateMap) {
	    stateStrands.resize(stateMap.size());
	    for (int i=0; i<stateMap.size(); i++) 
		stateStrands[i] = !isOnFStrand(stateMap[i]);
	}
	    
   private:
	// int copyOffset;
	static vector<bool> stateStrands;
	static double fullIntronBonus;
	static double minIntronFactor;
	    
    }; // struct PP::SubstateModel

    
    class ExonScorer {
    public:
	virtual ~ExonScorer() {
#ifdef DEBUG
	    if (backward_mode)
		delete currentHSColl;
#endif
	}
	void newFirstCodon(int beginOfBioExon) {
	    beginOfFirstCodon = beginOfBioExon + mod3(endOfLastCodon - beginOfBioExon +1);
	    aa_count = (endOfLastCodon - beginOfFirstCodon + 1)/3;
	}

	// score initial / rterminal exons
	virtual void scoreFirst(ViterbiColumnType& col, int predState, const Double& transEmiProb) = 0;
 	// score internal / rinternal exons
	virtual void scoreInternal(ViterbiColumnType& col, int predState, const Double& transEmiProb) = 0;
	void score(ViterbiColumnType& col, int predState, const Double& transEmiProb) {
	    (this->*scoreFunPtr)(col, predState, transEmiProb);
	}
	virtual void postProcessing(const Double& maxProb) = 0;
	virtual void exportSubstates(ViterbiSubmapType&) {}
	virtual SubstateId getPredSubstate() {
	    return SubstateId();
	}
	virtual bool hasNewMax() {
	    return false;
	}
	virtual void addMatches(int, int) {}
	bool validSize() {
	    return aa_count > 0;
	}
#ifdef DEBUG
	virtual Double getMaxProb() {
	    return 0;
	}
	static Matrix<Double>* transitions;
#endif		
    protected:
	ExonScorer(SubstateModel& model, int endOfBioExon, StateType type) : 
	    scoreFunPtr(0),
	    complement(!isOnFStrand(type)),
	    mdl(model),
	    blockScores(model.blockScores[complement]),
	    beginOfFirstCodon(-1),
	    endOfLastCodon(isOnFStrand(type) ? 
			   endOfBioExon - stateReadingFrames[type] :
			   endOfBioExon - (2 - stateReadingFrames[type])),
    	    currentHSColl(&model.hitSeqColl[complement][mod3(endOfLastCodon+1)]),
	    is_active(true),
#ifdef DEBUG
	    is_first(isFirstExon(type)),
	    is_last(isLastExon(type)),
#endif
	    backward_mode(false)
	{
	    if (isFirstExon(type)) 
		scoreFunPtr = &ExonScorer::scoreFirst;
	    else 
		scoreFunPtr = &ExonScorer::scoreInternal;
	    
	    switch (type) {
		case terminal:  case singleG:     // subtract the stop codon
		    endStateOffset=-1; break;
		case initial1: case initial2:     // add the split codon 
		case internal1: case internal2:
		case rinternal0: case rinternal1: 
		case rterminal1: case rterminal0:
		    endStateOffset=1; break;
		case initial0: case internal0:    // splicesite at codon end
		case rinternal2: case rinitial:
		case rterminal2: case rsingleG:   // reverse gene start at exon end
		    endStateOffset=0; break;
		default:
		    throw ProjectError("Internal error: State type " + itoa(type) + " unsuitable for exon scorer.");
	    }
	}


	// accessing prfl
	int blockCount() const {
	    return mdl.blockCount();
	}
	int blockSize(int b) const {
	    return mdl.blockSize(b,complement);
	}
	const Block& getBlock(int b) const {
	    return mdl.getBlock(b, complement);
	}
	bool exceedsBlock(Position ppos) const {
	    return ppos.i > blockSize(ppos.b);
	}
	DistanceType interBlockDist(int b) const {
	    return mdl.interBlockDist(b,complement);
	}
	DistanceType interBlockFullDist(int b) const {
	    if (b==0) 
		return interBlockDist(b);
	    else 
		return interBlockDist(b) + Range(blockSize(b-1));
	}
	Range interBlockRange(int b) const {
	    return interBlockDist(b).r;
	}
	Range interBlockFullRange(int b) const {
	    return interBlockFullDist(b).r;
	}
	Position asNextPos(Position ppos) {
	    ppos.i -= blockSize(ppos.b);
	    ppos.i -= interBlockRange(++ppos.b).max;
	    return ppos;
	}

	// which scoring function to use
	void (ExonScorer::*scoreFunPtr)(ViterbiColumnType&, int predState, const Double&);

	// internal data fields
	bool complement;
	SubstateModel& mdl;
	vector<BlockScoreType>& blockScores;

	/*
	 * codon coordinates of the currently observed exon
	 *
	 * endOfLastCodon is the position of the rightmost nucleotide
	 * with reading frame 2 (forward) resp. 0 (backward) in the
	 * biological exon
	 * endOfLastCodon is the position until which profile scores
	 * are calculated for exons no matter where "base" is in this case
	 */
	int aa_count, beginOfFirstCodon, endOfLastCodon;

	// these are the block hit sequences
	HitSequenceCollection* currentHSColl;

	bool is_active;
#ifdef DEBUG
	bool is_first, is_last;
#endif
	int endStateOffset;   // in reading frames != 0, the state id 
                              // refers to the end of the split codon,
	                      // and not to the last full codon in the exon
	                      // we set endStateOffset:=1 in these cases

	bool backward_mode;  // true if called from sampling or backtracking
    }; // class PP::ExonScorer


    class SingleTargetExonScorer : public ExonScorer {
    public:
 	SingleTargetExonScorer(SubstateModel& model, StateType etype, int endOfBioExon);
 	SingleTargetExonScorer(SubstateModel& model, StateType etype, SubstateId substate,
			       int endofBioExon, int ORFleft);	
	void calcAllowedDist(bool is_last) {
	    allowedDist.has_max = interBlockDist(rightPos.b).has_max || rightPos.i<0;
	    allowedDist.r = 
		(rightPos.i > 0 || is_last) ? 
		interBlockRange(rightPos.b).alignRight(rightPos.i) : rightPos.i;
	}
	void createHitSequences(int ORFleft);

	// score single, backtrack initial/rterminal/single exons
	virtual void scoreFirst(ViterbiColumnType& col, int predState, const Double& transEmiProb);

	// score terminal/rinitial, backtrack terminal/rinitial/internal/rinternal
	virtual void scoreInternal(ViterbiColumnType& col, int predState, const Double& transEmiProb);

	virtual void postProcessing(const Double& mainProb);
	virtual void exportSubstates(ViterbiSubmapType&);
	virtual SubstateId getPredSubstate() {
	    return predSubstate;
	}
	virtual bool hasNewMax() {
	    return newmax;
	}
	virtual void addMatches(int beginOfBioExon, int endOfBioExon) {
	    newFirstCodon(beginOfBioExon);
	    if (rightPos.i - endStateOffset > aa_count) {
		mdl.addInternMatch(rightPos - endStateOffset, beginOfBioExon, endOfLastCodon, endOfBioExon, complement);
		return;
	    } 
	    if (rightPos.i > 0)
		mdl.addPrefixMatch(rightPos - endStateOffset, endOfLastCodon, endOfBioExon, complement);
	    mdl.addHitSeqMatch(bestHitSeq, complement);
	    if (predSubstate.slot>=0) {
		Position ppos = predSubstate;
		if (ppos.i > 0)
		    mdl.addSuffixMatch(ppos, beginOfFirstCodon, beginOfBioExon, complement);
	    }
	}
	// SubstateId rightId() const {
	//     return rightPos.id();
	// }
#ifdef DEBUG
	Double getFactor();
	Double getPrefixFactor() const { 
	    return prefixFactor; 
	}
	virtual Double getMaxProb() {
	    return maxProb;
	}
	void checkResults(const OptionListItem& oli, const ViterbiMatrixType& viterbi, 
			  int, int, Double endPartProb, Double notEndPartProb);
#endif

    private:
	Double bestFittingBonusFactor(int b, int maxBlockStart, const HitSequence*& argmax);
	
	Position rightPos;
	DistanceType allowedDist; // this is the range of values for targetPos.i compatible with rightPos
	SubstateId predSubstate;
	Double maxProb;
	bool newmax;
	ViterbiSubmapType* bestSource;
	const HitSequence* bestHitSeq;
	Double prefixFactor;
    }; // class PP::SingleTargetExonScorer


    class MultiTargetExonScorer : public ExonScorer {
    public:
	MultiTargetExonScorer(SubstateModel& model, ViterbiSubmapType&,
			      StateType etype, int endOfBioExon);
	// no destructor needed: 
	// bestLeftId is owned by the ViterbiSubmap, source just linked to one
	// ~MultiTargetExonScorer() {}

    private:
	void initBonusProbs(BonusProbType&, int b, int frame);
	void updateBonusProbs(Position targetPos, const Double& bonusProb);
	bool scoreSourceEntry(int b, int maxBlockStart, Double bonusProb);
 
	// score initial/rterminal exons
	virtual void scoreFirst(ViterbiColumnType& source, int predState, const Double& transEmiProb);

	// score internal/rinternal exons
	virtual void scoreInternal(ViterbiColumnType& source, int predState, const Double& transEmiProb);

	virtual void postProcessing(const Double& maxProb);
	
	ViterbiSubmapType& targetVit;
	vector<BonusProbType> bonusProbs;
	ViterbiSubmapType* source;
	SubstateId leftId;
	SubstateIdMap* bestLeftId;
    }; // class PP::MultiTargetExonScorer

} // namespace PP
#endif
