/*****************************************************************************\
 * Filename : pp_fastBlockSearcher.hh
 * Author   : Oliver Keller
 * Project  : Gene Prediction with Protein Family Patterns
 *
 * Description: A fast block search class to determine sequence parts
 *              relevant for the ppx extension
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 09.09.09   | Oliver Keller         | creation of the file
 \******************************************************************************/

#ifndef __PP_BLOCKSEARCHER_HH
#define __PP_BLOCKSEARCHER_HH

// project includes
#include "pp_profile.hh"

// standard C/C++ includes
#include <deque>
#include <vector>
#include <map>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

namespace PP {
    
    /*
     * FsSeedCollection: Contains, for each triplet of residues, a 
     *                   set of positions that this tuple is likely
     *                   to appear as part of a block.
     *                   The decision whether a triplet is considered
     *                   a seed for a certain position depends on two
     *                   parameters: 
     *                   - the maximum number of seeds (chosen
     *                     as a subset of 8000 possible triplets) at a
     *                     certain position, reversely proportional to
     *                     the length of the block (so that the total
     *                     number of seeds used for a block is independent
     *                     of its length)
     *                   - if, according to the profile, a high probability
     *                     that a given triplet is a seed can be reached with
     *                     a smaller number of seeds, this reduced number is
     *                     used.
     *                   Given this number, the most likely (i.e. highest
     *                   scoring) triplets are taken into the collection of 
     *                   seeds.
     *
     * The use of a seed collection has the advantage that it can be calculated
     * once and used for filtering during the sequence analysis just by looking
     * at each 9-tuples of nucleotides in the sequence without calculating scores
     * at each point. Furthermore, the decision whether a 9-tuple translates to 
     * a seed can be made simultaneously for all offsets of all blocks at once.
     *
     */
    
    class FsSeedCollection {
//	multimap<string, Position> seeds;
	vector<int> vindex;
	vector<Position> vseeds;
	
	static int expSeedCount;   // maximal number of seeds per block
	static double maxCoverage; // position-specific coverage considered sufficient
//	multimap<string, Position>::iterator current;  // used internally in getNext()
	int current;
	int currend;

#ifdef DEBUG	    
	static double blockP(const PP::Block& blk, string s, int from) {
	    double result = 1;
	    for (int i=0;  i<s.length(); i++) 
		result *= blk[from+i][GeneticCode::get_aa_from_symbol(s[i])];
	    return result;
	}
	static double blockQ(const Block& blk, string s, int from) {
	    double result = 1;
	    for (int i=0;  i<s.length(); i++) 
		result *= blk[from+i].Q(GeneticCode::get_aa_from_symbol(s[i])).doubleValue();
	    return result;
	}
#endif

	static double blockP(const PP::Block& blk, int aas, int from) {

	    double result = 1;
	    for (int i=from+2;  i>=from; i--) {
		result *= blk[i][aas % 20];
		aas /= 20;
	    }
	    return result;
	}
	static double blockQ(const Block& blk, int aas, int from) {
	    double result = 1;
	    for (int i=from+2;  i>=from; i--) {
		result *= blk[i].Q(aas % 20).doubleValue();
		aas /= 20;
	    }
	    return result;
	}

    public:
	// create a seed collection
	FsSeedCollection(const PP::Profile& prfl) :
	    vindex(1,0), current(0), currend(0) 
	{
	    multimap<int, Position> seeds;
	    multimap<int, Position>::iterator it;


	    for (int b = 0; b<prfl.blockCount(); b++) {
		const Block& blk = prfl[b];
#ifdef VERBOSE_SEEDS
		cerr << "Determining seeds for block " << b << " (" << blk.id << ", " << blk.size() << " columns).\n";
#endif
		int maxcount = expSeedCount / blk.size();

		for (int i=0; i<=blk.size()-3; i++) {
		    multimap<double,int> currtriples;
		    string s(3, 0);
		    int val=0;
		    for (int t1 = 0; t1<NUM_AA; t1++) {
			s[0] = GeneticCode::aa_symbols[t1];
			for (int t2 = 0; t2<NUM_AA; t2++) {
			    s[1] = GeneticCode::aa_symbols[t2];
			    for (int t3 = 0; t3<NUM_AA; t3++) {
				s[2] = GeneticCode::aa_symbols[t3];
				double q_val = blockQ(blk, val, i);
				currtriples.insert(make_pair(q_val, val));
				val++;
			    }
			}
		    }
		    int count = 0;
		    double p = 0.0;
		    map<double,int>::reverse_iterator it = currtriples.rbegin();
		    for (count = 0; count < maxcount && p < maxCoverage; count++) {
			p += blockP(blk, it->second, i);
			seeds.insert(make_pair(it->second, Position(b,i)));
			++it;
		    }
#ifdef DEBUG
		    cerr << "[" << i << ".." << (i+2) << "]: coverage=" << p << "   minBlockQ=" << (--it)->first << " count=" << count << endl;
#endif
		}
	    }
	    vseeds.reserve(seeds.size());
	    vindex.reserve(8001);
	    for (it=seeds.begin(); it != seeds.end(); ++it) {
		while (vindex.size() <= it->first) 
		    vindex.push_back(vseeds.size());
		vseeds.push_back(it->second);
	    }
	    while (vindex.size()<=8000) 
		vindex.push_back(vseeds.size());
	}
	
	// return positions for which s is a seed (or null if there are none)
	// the first
	void setFirst(int patt) {
	    current = vindex[patt];
	    currend = vindex[patt+1];
	}
	bool hasNext() {
	    return current < currend;
	}
	// the next
	PP::Position& getNext() {
	    return vseeds[current++];
	}
	
	int size() {
	    return vseeds.size();
	}
    }; // PP::FsSeedCollection

    /*
     * FsHitType
     */
    struct FsHitType;
    typedef deque<FsHitType*> HitQueue;

    struct FsHitType {
	FsHitType(int p, int b, bool r, PartScoreType sc) :
	    pos(p),
	    head(this),
	    blockNo(b),
	    reverse(r),
	    score(sc.score), 
	    blockfrom(sc.from),
	    blockto(sc.to),
	    pathScore(sc.score), 
	    predecessor(0) {}
	int pos;
	int start() {
	    return head->pos;
	}
	FsHitType* head;
	int blockNo;
	bool reverse;
	double score;
	int blockfrom;
	int blockto;
	double pathScore;
	FsHitType* predecessor;
	void linkTo(HitQueue& queue) {
	    while(!queue.empty() && queue.front()->pos < pos - maxIntronLen) 
		queue.pop_front();
	    if (!queue.empty()) {
		predecessor = queue.front();
		head = predecessor->head;
		pathScore = predecessor->pathScore - intronMalus*(pos - predecessor->pos) + score;
	    }
	}
	void pushOn(HitQueue& queue) {
	    while (!queue.empty()) {
		FsHitType* ht = queue.back();
		if (ht->pathScore < pathScore + intronMalus*(pos - ht->pos))
		    queue.pop_back();
		else
		    break;
	    }
	    queue.push_back(this);
	}
	static int maxIntronLen;
	static double intronMalus;
    }; // PP::FsHitType

    /*
     * FsHitCollection
     */
    class FsHitCollection {
	vector<HitQueue> pendingHits[2];
	vector<FsHitType*> finalResult;
	vector<FsHitType*> allHits;
	int size;
    
    public:
	FsHitCollection(int n) : size(n) {
	    pendingHits[0].resize(n);
	    pendingHits[1].resize(n);
	}
	~FsHitCollection() {
	    for (int i=0; i<allHits.size(); i++) 
		delete allHits[i];
	}
	int allHitCount() {
	    return allHits.size();
	}
	int resultCount() {
	    return finalResult.size();
	}
	   
	void newHit(FsHitType*);
	void storeBestResults(multimap<double, FsHitType>&, int mincount, double threshold);
    }; // PP::FsHitCollection



    /*
     * PP::CandidateCollection: collects "candidates" for block hits
     * 
     * Each substring of the input sequence is a potential candidate. 
     * The collection is initialized with the profile (sequence of blocks)
     * and the seed collection (computed once in advance) containing the
     * seeds (see above)
     * 
     */

    // this type represents the results for a particular part of the
    // input sequence when translated and aligned to a block
    // - how many residues of the translation are part of a seed 
    //   (with respect to the corresponding block position)
    // - which is the rightmost residue that is part of seed
    struct ScoringCandidate {
	ScoringCandidate() : inSeedCount(0), lastOffset(-10) {}
	int inSeedCount;
	int lastOffset;
    };

    struct twomer {
	signed char v[2];
	twomer() { v[0]=v[1]=-1; }
	signed char& operator[] (int n) {
	    return v[n];
	}
	int val() { 
	    return (v[0]<0 || v[1]<0) ? -1 : (v[0]*20 + v[1]);
	}
    };
       

    struct CandidateCollection {
	typedef deque<ScoringCandidate> SeedCountsQueue;
	// we have, for both strands and each block a queue
	// of scoring candidates; the index for the queue is
	// given by the distance of the block match start to
	// the last observed nucleotide
	vector<SeedCountsQueue> theCandidates[2];
	FsSeedCollection& seedColl;
	const Profile& prfl;
	
	vector<int> maxseedcount; // maximum of inSeedCounts
	vector<int> seedhitsizes; // sum of inSeedCounts 
	vector<int> candcounts;   // for how many locations bestPartialLogScore is called
	vector<int> realhitcounts; // ...and successful
	
	// pending 2-mers of amino acids, one waiting for each frame/strand
	deque<twomer> patterns;


	CandidateCollection(const Profile& p, FsSeedCollection& coll) :
	    seedColl(coll),
	    prfl(p),
	    maxseedcount(p.blockCount(), 0),
	    seedhitsizes(p.blockCount(), 0),
	    candcounts(p.blockCount(), 0),
	    realhitcounts(p.blockCount(), 0),
	    patterns(6, twomer())
	{
	    // start the computation with a queue of emtpy scoring candidates
	    for (int b=0; b<p.blockCount(); b++) 
		theCandidates[0].push_back(SeedCountsQueue(p.blockSize(b) * 3 -8, 
							ScoringCandidate()));
	    theCandidates[1] = theCandidates[0];
	}
	    
	// put SeedsToEntries: given a 3-aminoacid-motif, add it to those scoring
	// candidates for which it represents a seed at the corresponding position
	void putSeedsToEntries(int currpat, bool reverse) {
	    seedColl.setFirst(currpat);
	    while (seedColl.hasNext()) {
		PP::Position& pos = seedColl.getNext();
		int index = (reverse ? prfl.blockSize(pos.b) - pos.i - 3 : pos.i) * 3;
		    
		ScoringCandidate& pf = theCandidates[reverse][pos.b][index];
#ifdef DEBUG
		if (pf.lastOffset >= 0 && (pf.lastOffset == pos.i || reverse == (pf.lastOffset < pos.i)))
		    cerr << "ERROR: last offset (" << pf.lastOffset << ") " << (reverse ? "<" : ">") 
			 << " current offset (" << pos.i << ") !" << endl;
#endif
		int range = abs(pf.lastOffset-pos.i);
		pf.inSeedCount += range < 3 ? range : 3;
		pf.lastOffset = pos.i;
	    }
	}

	// pushSeeds: move forward the sequence: for each strand, the next pending
	// pair of aminoacids is concatenated to the next residue found and putSeedsToEntries
	// is called with the resulting triple; and new pairs are added to the pending queue
	void pushSeeds(int t) {
	    signed char aa, aab;
	    twomer nextPattern = patterns.back();
	    patterns.pop_back();
	    twomer revPattern = patterns.back();
	    patterns.pop_back();
	    try {
		aa = GeneticCode::map[Seq2Int(3)(PP::DNA::sequence + t-2)];
		aab = GeneticCode::map[Seq2Int(3).rc(PP::DNA::sequence + t-2)];
		int forval = nextPattern.val()*20;
		int revval = revPattern.val();
		if (aa >= 0 && forval >= 0) 
		    putSeedsToEntries(forval + aa, false);
		if (aab >= 0 && revval >= 0)
		    putSeedsToEntries(revval + 400*aab, true);
	    } catch (InvalidNucleotideError e) {
		aa = aab = -2;
	    }
	    nextPattern[0] = nextPattern[1];
	    nextPattern[1] = aa;
	    patterns.push_front(nextPattern);
	    revPattern[1] = revPattern[0];
	    revPattern[0] = aab;
	    patterns.push_front(revPattern);
	}

	// takeCandidates: move forward the candidate queues (one for each block and strand)
	// if a scoring candidate is complete (fully evaluated) and its seed count
	// is above a threshold (here: around 25% of the block width), call it a hit candidate,
	// compute log scores for it, and if there is one non-negative, call it a hit, and
	// add it to the FsHitCollection target
	void takeCandidates(int t, FsHitCollection& target, bool reverse) {
	    for (int b=0; b < prfl.blockCount(); b++) {
		theCandidates[reverse][b].push_front(ScoringCandidate());
		int hc = theCandidates[reverse][b].back().inSeedCount;
		theCandidates[reverse][b].pop_back();
		seedhitsizes[b] += hc;
		if (hc < 0) 
 		    cerr << "ERROR: seedhitcount < 0!\n";
		if (maxseedcount[b] < hc)
		    maxseedcount[b] = hc;
		if (hc > 4 + prfl.blockSize(b)/4) {
		    int blstart = t - prfl.blockSize(b) * 3;
		    if (blstart < 0)
			continue;
		    candcounts[b]++;
		    PP::PartScoreType pt;
		    prfl[b].bestPartialLogScore(reverse, blstart, pt);
		    int width = pt.to - pt.from;
		    if (pt.score >= 0 && 
			width >= MIN_BLOCKSIZE && 
			width >= 0.3*prfl.blockSize(b)) {
			FsHitType* newHit = new FsHitType(blstart, b, reverse, pt);
			target.newHit(newHit);
			realhitcounts[b]++;
		    }
		}
	    }
	}

	// this is the entry point to the candidate collection
	// calls takeCandidates and pushSeeds
	void proceed(int t, FsHitCollection& target) {
	    takeCandidates(t, target, false);
	    takeCandidates(t, target, true);
	    pushSeeds(t);
	}

	void debugging_output() {
	    for (int b=0; b<prfl.blockCount(); b++) {
		cerr << "Block " << b << " (" << prfl[b].id << ")" << endl;
		cerr << "    len: " << prfl.blockSize(b) << endl;
		cerr << "   av c: " << double(seedhitsizes[b])/DNA::len;
		cerr << "  max c: " << maxseedcount[b] << endl;
		cerr << "  % c>t: " << fixed << setw(4) << setprecision(2) << double(candcounts[b]*100)/DNA::len << endl;
	    }
	}
    }; // class PP::CandidateCollection
} // namespace PP
#endif // __PP_BLOCKSEARCHER_HH
