/*****************************************************************************\
 * Filename : pp_profile.hh
 * Author   : Oliver Keller
 * Project  : Gene Prediction with Protein Family Patterns
 *
 * Description: Data structures containing the ProteinProfile /
 *              IntronProfile 
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 22.03.07   | Oliver Keller         | creation of the file
 * 06.11.08   | Oliver Keller         | first version in SVN
 * 28.01.09   | Oliver Keller         | rev. 157, final forward-only version
 \******************************************************************************/

#ifndef __PP_PROFILE_HH
#define __PP_PROFILE_HH

// project includes
#include "geneticcode.hh" // needed to score codons
#include "vitmatrix.hh"   // for SubstateId's

// standard C/C++ includes
#include <deque>
#include <cmath>  // for sqrt
#include <map>
#include <vector>

using namespace std;


// constants used in profile classes
// how many bits are available for i when coding PP::Position as int
#define INTERBLOCKBITS 15
#define MAXINTERBLOCKDIST ((1<<INTERBLOCKBITS)-1)
// this might well be reduced in the future to enable more profile-related
// substates
#define MAX_BLOCKCOUNT 64


// fraction by which the distance range is relaxed
// PP::DistanceType::makeTolerant allows for a range of (1-x)min,(1+x)max
#define RELAXATION 0.05

// minimum size of a partial block that can be discarded
#define MIN_CHECKCOUNT 3
#define MIN_BLOCKSIZE 6

// default values for quantiles (in units of stddev)
// that are transformed into scoring thresholds
#define MIN_SPEC 4.0
#define MIN_SENS 0.4
#define MIN_ANCHOR_SPEC 4.0
#define MIN_ANCHOR_COUNT 1
#define PARTIAL_SPEC 4.5
#define PARTIAL_SENS 2.0
#define GLOBAL_THRESH 2.5


struct Dist {
    double mu;
    double var;
    Dist(double m=0, double v=0) : mu(m), var(v) {}
    Dist& add(double m=0, double v=0) { 
	mu+=m; var+=v; return *this;
    }
    Dist& operator+=(const Dist& other) {
	return add(other.mu, other.var); 
    }
    Dist& operator-=(const Dist& other) {
	return add(-other.mu, -other.var); 
    }
    Dist operator+(const Dist& other) const {
	return Dist(*this) += other;
    }
    Dist operator-(const Dist& other) const {
	return Dist(*this) -= other;
    }

    double stddev() const {
	return sqrt(var>=0?var:-var);
    }

    double normed(double abs) const {
	return (abs - mu)/stddev();
    }
    double abs(double normed) const {
	return normed * stddev() + mu;
    }

};

inline Dist operator-(const Dist& d) {
    return Dist(-d.mu, -d.var);
}

/*
 * A data structure representing the intron profile of a block profile
 *
 */
struct ProfileMap : private map < pair<int, int>, double > {
    void add(int n, int f, double p) {
	insert(make_pair(make_pair(n,f), p));
    }
    double get(int n, int f) const {
	const_iterator it = find(make_pair(n,f));
	return (it == end()) ? 0 : it->second;
    }
}; // class ProfileMap


namespace PP {
    /*
     * Exception class thrown from PP:: classes
     */
    struct ProfileNotFoundError : ProjectError {
 	ProfileNotFoundError(string filename) :
	    ProjectError("PP::Profile: Could not open file \"" + filename + "\".")
	{}
    };
    struct ProfileReadError : ProjectError {
	ProfileReadError(string filename, int lineno) :
	    ProjectError("PP::Profile: Error parsing pattern file\"" + filename + 
			 "\", line " + itoa(lineno) + ".")
	{}
    };

    struct ProfileParseError {
	int lineno; 
	ProfileParseError(int n) : lineno(n) {}
    };

    struct PartParseError {
	int offset;
	PartParseError(int n) : offset(n) {}
    };

    struct ProfileInsigError : ProjectError {
	ProfileInsigError(string message) :
	    ProjectError("PP::Profile: " + message) {}
    };
    
    
    /*
     * A data structure representing an integer interval
     */
    struct Range {
	Range(int x=0) : min(x), max(x) {}
	Range(int l, int h) : min(l), max(h) { 
	    if (min>max) { min = h; max = l; }
	}
	Range operator+= (Range other) {
            min+=other.min;
	    max+=other.max;
            return *this;
        }
        Range operator+ (Range other) const {
	    return Range(min,max)+=other;
        }
	bool has(int elem) const {
	    return min <= elem && elem <= max;
	}
	int size() const { // <= 0 is equivalent to
	    return max-min+1;
	}

	// move interval such that max equals newmax
	Range alignRight(int newmax) const {
	    return Range(newmax-max+min, newmax);
	}
	// move interval such that min equals newmin
	Range alignLeft(int newmin) const {
	    return Range(newmin, newmin+max-min);
	}
	
	// data fields
	int min, max;
    }; // class PP::Range

    inline Range operator- (Range r) {
	return Range(-r.max, -r.min);
    }

    struct DistanceType {
	bool has_max;
	Range r;
	DistanceType() : has_max(true), r() {}
	DistanceType(Range other) : has_max(true), r(other) {}
	friend ostream& operator<< (ostream& strm, DistanceType d) {
	    if (d.has_max)
		return strm <<  d.r.min << "\t" << d.r.max;
  	    else
 		return strm << d.r.min << "\t*";
	}
        friend istream& operator>> (istream& strm, DistanceType& d) {
	    d.has_max = true;
	    strm >> d.r.min >> ws;
	    if (strm.peek() == '*') {
		strm.get();
		d.has_max = false;
		d.r.max = d.r.min;
		return strm;
	    }
	    return strm >> d.r.max;
        }
	DistanceType& operator+= (DistanceType other) {
	    r += other.r;
	    has_max = has_max && other.has_max;
	    return *this;
	}
	DistanceType operator+ (DistanceType other) const {
	    return DistanceType(*this) += other;
	}
	void setInfMax () {
	    r.max = r.min;
	    has_max = false;
	}
	void makeTolerant() {
	    r.min = int(r.min * (1-RELAXATION) + 0.5);
	    if (has_max) {
		r.max = int(r.max * (1+RELAXATION) + 0.5);
		if (r.max >= MAXINTERBLOCKDIST) 
		    setInfMax();
	    } else
		r.max = r.min;
	}
	bool has(int elem) const {
	    return has_max ? r.has(elem) : r.min <= elem;
	}
    };

    /*
     * A data structure representing a distribution of amino acids
     */
    class Column {
    public:
	Column() { values[0]=1.0; }
 	Column(const double* val) {
	    *this = val;
	}
	Column& operator= (const double*);

	// initialize odds-ratios factors 
	void initRatios();

	// return expectation and variance of odds-ratios
	Dist getDist(const Column& model) const;
	Dist getOwnDist() const {
	    return getDist(*this);
	}
	Dist getBackDist() const {
	    return getDist(background);
	}

	// accessing the Column entries
	// by value
	double operator[](int n) const { 
	    return (0 <= n && n<NUM_AA) ? values[n] : -1.0; 
	}
	// by oddRatio
	Double Q(int n) const {
	    return (0 <= n && n<NUM_AA) ? oddRatios[n] : stopCodonScore; 
	}
	Double Q(char c) const {
	    return Q(GeneticCode::get_aa_from_symbol(c));
	}
	// by LogOddRatio
	double L(int n) const { 
	    return Q(n).log(); 
	}

	// stream operations
	friend ostream& operator<<(ostream&, const Column&);
	friend istream& operator>>(istream&, Column& c);
	
	static const Double stopCodonScore;
	static double invalidScore;
	static double weight;
    private:
 	// internal data fields
        double values[NUM_AA];  // these numbers will be positive and
	                        // add up to one; class takes care of the normalisation 
	double oddRatios[NUM_AA]; // the odds against the background distribution

	static const Column background;
	static const double minFreq;   // a frequency of 0.0001 is always ensured by the class;
	                               // but we don't deal here with pseudocounts
    }; // class PP::Column



    struct Position {
	Position() :  b(-1), i(0) {}

	Position(int blockno, int columnNr) :
	    b(blockno), i(columnNr) {}
	Position(const SubstateId& id) : 
	    b(id.slot % MAX_BLOCKCOUNT),
	    i(id.value) {}
	Position operator+ (int l) {
	    return Position(b, i + l);
	}
	Position operator- (int l) {
	    return operator+(-l);
	}
	Position operator++(int) {
	    return Position(b,i++);
	}
	bool operator< (Position other) const {
	    return (b < other.b) || (b == other.b && i < other.i);
	}
	bool operator== (Position other) const {
	    return i==other.i && b==other.b;
	}
	bool operator!= (Position other) const {
	    return !(other == (*this));
	}
	bool operator<= (Position other) const {
	    return !(other < (*this));
	}
	void nextB() { b++; i=0; }

	// static int getB(int id) {
	//     return id >> (INTERBLOCKBITS+1);
	// }
	// static int getI(int id) {
	//     static const int mask = (-1) << (INTERBLOCKBITS+1);
	//     static const int offset = 1 << INTERBLOCKBITS;
	//     return (id | mask) + offset;
	// }
	// static Position fromId(int id) {
	//     if (id < 0) return Position();
	//     return Position(getB(id), getI(id));
	// }

	Position operator= (const SubstateId& id) {
	    b = id.slot % MAX_BLOCKCOUNT;
	    i = id.value;
	    return *this;
	}
	SubstateId id() const {
	    return b<0 ? SubstateId() : SubstateId(b,i);
	}
	SubstateId copyId() const {
	    return b<0 ? SubstateId() : SubstateId(b + MAX_BLOCKCOUNT, i);
	}

     	short b, i;
   }; // class PP::Position
    
    class BlockScoreType;
    
    struct PartScoreType {
	int from, to;
	double score;
    };

    class Block {
//	friend class BlockScoreType;
    public:
	Block(DistanceType, const vector<string>&, string);
	~Block() {
	    delete iP;
	}
	// score(x,b,l) is the emission probability factor for the sequence starting 
	// at dna position x, from the partial block of length l (codons), starting at 
	// offset b
	// BlockScoreType score;
	void initDistributions();
  	bool initThresholds();
	void setIntronProfile(const vector<string>&);
	int size() const { return columns.size(); }
	const Column& operator[] (int columnNr) const { return columns[columnNr]; }
	// bool addBlocksUntil(int newbase, map<int,Double>* result=NULL);
	// void removeBlocksUntil(int newbase) {
	//     newbase -= 3 * size();
	//     while (score.begin() < newbase)
	// 	score.pop_front();
	// }
//	Double saveNewScore();
	// Double savedScore(int offset) {
	//     map<int,Double>::iterator it = hits[offset%3].find(offset);
	//     return it == hits[offset%3].end() ? 0 : it->second;
	// }
	Double scoreFromScratch(bool complement, int dna_offset, int block_offset, int len) const;
	Double scoreFromScratch(bool complement, int dna_offset, int block_offset) const {
	    return scoreFromScratch(complement, dna_offset, block_offset, size()-block_offset);
	}
	Double scoreFromScratch(bool complement, int dna_offset) const {
	    return scoreFromScratch(complement, dna_offset, 0);
	}
	Double checkedSuffixScore(bool complement, int dna_offset, int block_offset) const;
	void bestPartialLogScore(bool complement, int dna_offset, PartScoreType&) const;

	// Double checkedPrefixScore(bool complement, int dna_offset, int len) const {
	//     Double result = scoreFromScratch(complement, dna_offset, 0,len);
	//     return result > getPrefixThresh(complement, len) ? result : 0;
	// }
	
	// const map<int,Double>& savedHits(int f) const {
	//     return hits[f];
	// }
#ifdef DEBUG
	void show_score(int n);
#endif       

	Block reverse(DistanceType dist) {
	    Block result;
	    result.distance = dist;
	    for (int i=size()-1; i>=0; i--) 
		result.columns.push_back(columns[i]);
	    return result;
	}

	friend void initConstants();

	Double getPrefixThresh(bool complement, int i) const;
	Double getSuffixThresh(bool complement, int i) const;
	Double getPartialThresh(bool complement, int from, int to) const {
	    return complement ? 
		thresholdMatrix[size()-from][size()-to] :
		thresholdMatrix[to][from];
	}
	Double getPartialThresh(bool complement, int from) const {
	    return getPartialThresh(complement, from, size());
	}
	Double getThreshold() const {
	    return threshold;
	}
	bool hasIP() const {
	    return iP != 0;
	}
	double getIntronFreq(int i, int frame) const {
	    return iP->get(i,frame);
	}
	bool isAnchor() {
	    return threshold > getSpecThresh(min_anchor_spec); //getBackDist().abs(min_anchor_spec);
	}
// 	Dist getBlockscoreBack() {
// 	    return backDists[0];
// 	}
	// the following methods use logscores
	Dist getOwnDist(int from, int to, bool complement = false) const {
	    return complement ? 
		ownDists[size()-to] - ownDists[size()-from] : 
		ownDists[from] - ownDists[to];
	}
	Dist getOwnDist(int from=0) const {
	    return ownDists[from];
	}
	Dist getBackDist(int from, int to, bool complement = false ) const {
	    return complement ? 
		backDists[size()-to] - ownDists[size()-from]:
		backDists[from] - backDists[to];
	}
	Dist getBackDist(int from=0) const {
	    return backDists[from];
	}
	double ownNormalize(double logscore, int from, int to, bool complement = false) const {
	    return getOwnDist(complement, from, to).normed(logscore);
	}
	double backNormalize(double logscore, int from, int to, bool complement = false) const {
	    return getBackDist(complement, from, to).normed(logscore);
	}
	double getSensThresh(double stddev, int from, int to) const {
	    return getOwnDist(from, to).abs(-stddev);
	}
	double getSensThresh(double stddev, int from=0) const {
	    return getOwnDist(from).abs(-stddev);
	}
	double getSpecThresh(double stddev, int from, int to) const {
	    return getBackDist(from, to).abs(stddev);
	}
	double getSpecThresh(double stddev, int from=0) const {
	    return getBackDist(from).abs(stddev);
	}
	const vector<Dist>& getAllOwnDists() const {
	    return ownDists;
	}

	string id;
	DistanceType distance; // distance range to previous block

	static const Double almostZero;
	static const double invalidScore;

    private:
	Block() : iP(0) {}

	vector<Column> columns;
	ProfileMap* iP;
//	map<int,Double> hits[3]; // will be updated by addBlocksUntil
 
#ifdef DEBUG
	vector<Double> prefixThresh, suffixThresh;
#endif
	vector< vector<Double> > thresholdMatrix;
	vector<Dist> backDists, ownDists;
	// Dist blockscoreBack;
	Double threshold; // threshold for full block scores to be inserted in the hitlist
	static double min_spec;
	static double min_sens;
	static double min_anchor_spec;
	static double partial_spec;
	static double partial_sens;
    }; // class PP::Block

    
    class BlockScoreType {
    public: 
	typedef vector<Double> aligned_type;    // partial (prefix) scores for the block at one dna position
	                                        // back() is the full block score

	typedef deque<aligned_type> array_type; // the aligned_types for all dna positions

	BlockScoreType(const Block& block) : 
	    offset(0), partner(&block) {}
	bool addBlocksUntil(bool complement, int newbase, map<int,Double>* result=NULL);
	void removeBlocksUntil(int newbase) {
	    newbase -= 3 * size();
	    while (begin() < newbase)
		pop_front();
	}
	Double savedScore(int offset) {
	    map<int,Double>::iterator it = hits[offset%3].find(offset);
	    return it == hits[offset%3].end() ? 0 : it->second;
	}
	const map<int,Double>& savedHits(int f) const {
	    return hits[f];
	}
	void init(int newoffset = 0) { 
	    offset = newoffset;
	    values.clear();
	    for (int f=0; f<3; f++) hits[f].clear();
	}
	void pop_front() { // removes the first block from the array and returns its score
	    values.pop_front();
	    ++offset;
 	}
	Double front() {
	    return values.front().back();
	}
// 	void push_back(const aligned_type& new_vals) {
// 	    values.push_back(new_vals);
// 	}
// 	void push_front(const aligned_type& new_vals) {
// 	    values.push_front(new_vals);
// 	    offset--;
// 	}
	const aligned_type& getAlignedScores(int  x) const {
	    return values.at(x-offset);
	}
	Double get(int dna_offset, int block_offset, int len)  const { 
	    if (len==0)
		return 1;
	    if (block_offset == 0)
		return getAlignedScores(dna_offset).at(len-1);
	    const aligned_type& alignedScores = getAlignedScores(dna_offset - 3*block_offset);
	    return alignedScores.at(len+block_offset-1)/alignedScores.at(block_offset-1);
 	}
	Double get(int dna_offset, int block_offset)  const { 
	    if (block_offset == 0)
		return get(dna_offset);
	    const aligned_type& alignedScores = getAlignedScores(dna_offset - 3*block_offset);
	    return alignedScores.back()/alignedScores.at(block_offset-1);
 	}
	Double get(int dna_offset) const {
	    return getAlignedScores(dna_offset).back();
	}
	const int& begin() const {
	    return offset;
	}
	int count() const {
	    return values.size();
	}
	int end() const {
	    return offset + count();
	}
 	int size() const {
	    return partner->size();
	}
	// Double threshold() const {
	//     return partner->threshold;
	// }
	Double prefixThresh(bool complement, int i) const {
	    return partner->getPrefixThresh(complement, i);
	}
	Double suffixThresh(bool complement, int i) const {
	    return partner->getSuffixThresh(complement, i);
	}
   private:
	aligned_type& new_back() {
	    values.push_back(aligned_type());
	    return values.back();
	}
		

	// internal data fields
	array_type values;
	int offset; // dna pos of first base of first evaluated block
	const Block* partner;
	map<int,Double> hits[3]; // will be updated by addBlocksUntil
    }; // class PP::BlockScoreType


	

    class Profile {
	static Column background;
	static double lowest_odds, exp_lowest_odds;

    public:
	// methods of class ProteinProfile
	// initialise pattern from file <filename>
	Profile(string filename);
	
    

// 	static bool idInBlock(int id) {
// 	    return ((id >> INTERBLOCKBITS)) % 2 > 0;
// 	}

// 	Position getPosOfId(int id) const {
// 	    if (id < 0) return Position();
// 	    Position result(getB(id), getI(id));
// #ifdef DEBUG
// 	    if (result.b > blockCount()) {
// 		cout << "getB() > #blocks\n";
// 	    } else if (result.i > blockSize(result.b)) {
// 		cout << "getI() > blocksize\n";
// 	    }
// #endif
// 	    if (result.i < 0) result.i <<= skipBits[result.b];
// 	    return result;
// 	}
// 	Position getPosOfCopyId(int id, bool comp=false) const {
// 	    return getPosOfId(id + copyOffset);
// 	}
	Position firstInRange(int blockno=0) const {
	    return Position(blockno, -interBlockRange(blockno).max);
	}

	int blockCount() const {
	    return blocks.size();
	}
	int blockSize(int blockno) const {
#ifdef DEBUG
	    if (blockno < 0 || blockno >= blocks.size())
		throw ProjectError("prfl.blockSize() called with index out of range");
#endif
	    return blocks[blockno].size();
	}
	DistanceType interBlockDist(int blockno) const {
	    return blockno == blockCount() ? finalDist : blocks[blockno].distance;
	}
	Range interBlockRange(int blockno) const {
	    return interBlockDist(blockno).r;
	}
	// bool has_max(int blockno) const {
	//     return interBlockDist(blockno).has_max;
	// }
	// DistanceType interBlockFullDist(int blockno) const {
	//     if (blockno == 0)
	// 	return blocks[0].distance;
	//     DistanceType result = interBlockDist(blockno);
	//     result.r += blocks[blockno-1].size();
	//     return result;
	// }
// 	int interBlockMaxDist(int blockno) const {
// 	    return interBlockFullRange(blockno).max;
// 	}
// 	bool isInRange(int blockno, int len) const {
// 	    return interBlockRange(blockno).has(len);
// 	}
	Column getColumn(Position pos) const {
	    return blocks[pos.b][pos.i];
	}
	string getName() const {
	    return name;
	}
	string getInfoLine() const {
	    ostringstream result;
	    int blockno=0;
	    while (true) {
		result << "--";
		DistanceType d = interBlockDist(blockno);
		if (d.has_max)
		    result << "[" << d.r.min << ".." << d.r.max << "]";
		else if (d.r.min>0)
		    
		    result << "[" << d.r.min << ",...]";
		result << "--";
		if (blockno == blockCount())
		    return result.str();
		result << "> " << blocks[blockno].id << " (" << blockSize(blockno) << ") <";
		blockno++;
	    }
	}
	ostream& write(ostream&) const;
	
	friend ostream& operator<< (ostream& strm, const Profile& P) {
	    return P.write(strm);
	}
	
	// Position finalPos() const {
	//     return Position(blockCount(),0);
	// }
	// SubstateId getTerminalId() const {
	//     return Position(blockCount(),0).id();
	// }

	Block& operator[] (int blockno) {
	    return blocks[blockno];
	}
	const Block& operator[] (int blockno) const {
	    return blocks[blockno];
	}
	Double getGlobalThresh(bool complement, Position ppos) const {
	    return 
		ppos.b == blockCount() ? 1 :
		globalThresh[complement][ppos.b][ppos.i >= 0 ? ppos.i : 0];
	}
	
	friend void initConstants();

#ifdef DEBUG
	static bool vv_equal(const vector<Double>& v1, const vector<Double>& v2) {
	    return equal(v1.begin(), v1.end(), v2.begin(), almost_equal);
	}
#endif

#ifdef DEBUG
	Profile reverse() {
	    Profile result;
	    vector< vector<Dist> > ownDists;
	    result.name = name + "(rev.)";
	    for (int b=blockCount(); b>0; b--) {
		result.blocks.push_back(blocks[b-1].reverse(interBlockDist(b)));
		result.blocks.back().initThresholds();
		ownDists.push_back(result.blocks.back().getAllOwnDists());
	    }
	    result.calcGlobalThresh(ownDists);
	    if (!equal(result.globalThresh[0].begin(),
		       result.globalThresh[0].end(),
		       globalThresh[1].begin(),
		       vv_equal) ||
		!equal(result.globalThresh[1].begin(),
		       result.globalThresh[1].end(),
		       globalThresh[0].begin(),
		       vv_equal) )
		throw ProjectError("globalThreshs not equal!");
	    for (int b=0; b<blockCount(); b++) 
		for (int i=1; i<result.blockSize(b); i++)
		    if (!almost_equal(result[b].getPrefixThresh(false, i),
				      blocks[blockCount()-b-1].getPrefixThresh(true,i)) ||
			!almost_equal(result[b].getPrefixThresh(true, i),
				      blocks[blockCount()-b-1].getPrefixThresh(false,i)))
			throw ProjectError("prefixThreshs not equal!");
		    else if (!almost_equal(result[b].getSuffixThresh(false, i),
					   blocks[blockCount()-b-1].getSuffixThresh(true,i)) ||
			     !almost_equal(result[b].getSuffixThresh(true, i),
					   blocks[blockCount()-b-1].getSuffixThresh(false,i)))
			throw ProjectError("suffixThreshs not equal!");
	    result.finalDist = interBlockDist(0);
	    return result;
	}
#endif

    private:
	Profile() {}
	void calcGlobalThresh(const vector< vector<Dist> >&); 
	void parse_stream(istream&);

	// data fields 
	string name;
	vector<Block> blocks;  // the blocks
	vector< vector<Double> > globalThresh[2];
	DistanceType finalDist;       // final non-block range
	ProfileMap iP;      // not implemented

	static int    min_anchor_count;
	static double global_thresh;
	static double absolute_malus_threshold;
// 	// moved to SubstateModel
// 	void addSkipBits(Range r) {
// 	    int n=0; 
// 	    while (r.max > (1<<INTERBLOCKBITS)) {
// 		r.max >>= 1;
// 		n++;
// 	    }
// 	    skipBits[0].push_back(n);
// 	}
// 	int copyOffset;

    };  // class PP::Profile		
	    
    struct DNA {
	static const char* sequence;
	static int len;
	static void initSeq(const char* dna, int seqlen) {
	    sequence = dna;
	    len = seqlen;
	}
	static void initSeq(string s) {
	    sequence = s.c_str();
	    len = s.length();
	}
    };

    void initConstants();
} // namespace PP

#ifdef DEBUG
namespace std {
inline ostream& operator<< (ostream& strm, PP::Position ppos) {
    return strm << ppos.b << "/" << ppos.i;
}
}
#endif

inline Double PP::Block::getPrefixThresh(bool complement, int i) const {
#ifdef DEBUG
    if (i<=0 || i>=size()) 
	throw ProjectError("Out of range!");
    Double result = complement ? 
	suffixThresh[size() - i] :
	prefixThresh[i];
    if (!almost_equal(result,getPartialThresh(complement, 0, i)))
	throw ProjectError("partialThresh bug");
#endif
    return getPartialThresh(complement, 0, i);
}

inline Double PP::Block::getSuffixThresh(bool complement, int i) const {
#ifdef DEBUG
    if (i<=0 || i>=size()) 
	throw ProjectError("Out of range!");
    Double result=complement ?
	prefixThresh[size() - i] :
	suffixThresh[i];
    if (!almost_equal(result,getPartialThresh(complement, i)))
	throw ProjectError("partialThresh bug");
#endif
    return getPartialThresh(complement, i);
}
#endif
