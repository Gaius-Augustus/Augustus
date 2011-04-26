/**********************************************************************
 * file:    pp_profile.cc
 * license: Artistic License, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  protein pattern search extension
 * authors: Oliver Keller
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 22.03.07| Oliver Keller | creation of the file
 * 28.01.09| Oliver Keller | rev. 157, final forward-only version
 **********************************************************************/

#include "pp_profile.hh"

// project includes
#include "properties.hh"  // for initConstants

// standard C/C++ includes
#include <iomanip>    // for setw, setprecision
#include <fstream>    // for ifstream
#include <stdexcept>  // for out_of_range
#include <iostream>

using namespace PP;

/* internal constants */
static const double default_amino_frq[NUM_AA] = { 0.07088, 0.05268, 0.06270, 0.05256, 0.05807, 
						  0.04439, 0.04037, 0.07068, 0.05837, 0.07689, 
						  0.06538, 0.09200, 0.05527, 0.03995, 0.03224,
						  0.01312, 0.02253, 0.02353, 0.01793, 0.05046 };



/*--- class Column -------------------------------------------------*/

const Column Column::background = default_amino_frq;
const double Column::minFreq = 0.0001;
const Double Column::stopCodonScore = pow(2.0L,-1000); // 2^(-1000), = almostZero^2

const Double Block::almostZero = pow(2.0L,-500); // upper bound for odds-Score
                                                 // for any sequence containing stop codons

double Block::min_spec = MIN_SPEC;
double Block::min_sens = MIN_SENS;
double Block::min_anchor_spec = MIN_ANCHOR_SPEC;
int    Profile::min_anchor_count = MIN_ANCHOR_COUNT;
double Block::partial_spec = PARTIAL_SPEC;
double Block::partial_sens = PARTIAL_SENS;
double Profile::global_thresh = GLOBAL_THRESH; 
double Profile::absolute_malus_threshold = 0;
double Column::weight = 1.0;
double Column::invalidScore = 1.0;
 
Column& Column::operator= (const double* val) {
    double sum = 0.0;
    for (int a=0; a<NUM_AA; a++) 
	if (val[a]<0)
	    throw out_of_range("Negative value for PP::Column");
	else
	    sum+=val[a]; 
    if (sum <= 0.0) {
	throw out_of_range("No positive value in PP::Column");
    }
    for (int a=0; a<NUM_AA; a++) 
	values[a]=(val[a]/sum)*(1-minFreq*NUM_AA)+minFreq;
    initRatios();
    return *this;
}

void Column::initRatios() {
    if (weight == 1.0)
	for (int a=0; a<NUM_AA; a++)  
	    oddRatios[a] = values[a]/background[a];
    else
	for (int a=0; a<NUM_AA; a++)  
	    oddRatios[a] = std::pow(values[a]/background[a], weight);
}

Dist Column::getDist(const Column& model) const {
    Dist result;
    for (int a=0; a<NUM_AA; a++)  {
	double logfactor = std::log(oddRatios[a]);
	double mu = model[a] * logfactor;
	result.add(mu, mu * logfactor);
    }
    result.var -= result.mu * result.mu;
    return result;
}

ostream& PP::operator<<(ostream& strm, const Column& c) {
    for (int a=0; a<NUM_AA; a++)
	strm << "\t" << c.values[a];
    return strm;
}

istream& PP::operator>>(istream& strm, Column& c) {
    double buf[NUM_AA];
    for (int a=0; a<NUM_AA; a++) 
	if (!(strm >> buf[a]))
	    return strm;
    c = buf;
    return strm;
}


/*--- class IntronProfile ------------------------------------------*/

/*
void IntronProfile::set(const vector<string>& lines) {
    string line;
    for (int lineno=0; lineno<lines.size(); lineno++) {
	istringstream sstrm(lines[lineno]);
	int from, to;
	double prob;
	if ((sstrm >> from >> to >> prob >> ws) && sstrm.eof())
	    (*this)[from][to] = prob;
	else
	    throw PartParseError(lines.size()-lineno);
    }
}

ostream& IntronProfile::pretty_print(ostream& strm) const {
    for (const_iterator it=begin(); it!=end(); it++) {
        strm << "from " << it->first << " to";
        const mapped_type& submap = it->second;
	for (mapped_type::const_iterator it2=submap.begin(); it2!=submap.end(); it2++)
            strm << " " << it2->first << ":"
		 << setw(6) << setprecision(3) << fixed << it2->second;
    }
    return strm << "\n";
}

double IntronProfile::getP(int from) const { 
    const_iterator it = find(from);
    if (it == end())
	return 1.0;
    const mapped_type& submap = it->second;
    return accumulate(submap.begin(), submap.end(), 1.0, sub_second);
}

double IntronProfile::getP(int from, int to) const {
    const_iterator it = find(from);
    if (it == end())
	return 0.0;
    const mapped_type& submap = it->second;
    return submap.count(to)? submap.find(to)->second : 0.0;
}

ostream& PP::operator<< (ostream& strm, const IntronProfile& prfl) {
    for(IntronProfile::const_iterator it = prfl.begin(); it != prfl.end(); ++it) {
	int from = it->first;
	IntronProfile::mapped_type target = it->second;
	for (IntronProfile::mapped_type::iterator it2 = target.begin(); it2 != target.end(); it2++)
	    strm << from << "\t" << it2->first << "\t" << it2->second << "\n";
    }
    return strm;
}

*/ 

/*--- class Block --------------------------------------------------*/

const char* DNA::sequence = 0;
int DNA::len = 0;


Block::Block(DistanceType d, const vector<string>& lines, string default_id) : 
    id(default_id), distance(d), iP(0) {
    // empty lines are ignored; lines saying "name=" set the id; all other lines lead to an exeption 
    // if they cannot be read successfully
    for (int lineno=0; lineno<lines.size(); lineno++) {
	int n;
	Column col = default_amino_frq;
	if (istringstream(lines[lineno]) >> n >> col && n == size())
	    columns.push_back(col);
	else if (lines[lineno].compare(0,5,"name=")==0)
	    id = lines[lineno].substr(5);
	else if (lines[lineno] != "")
	    throw PartParseError(lines.size()-lineno);
    }
}


#ifdef DEBUG
inline Double findThresh(Dist d0, Dist d1, double thr0, double thr1, int variant=0) {
    double
	min_thresh = d0.abs(thr0),
	max_thresh = d1.abs(thr1);
    if (min_thresh <= max_thresh)
	return Double::exp((min_thresh + max_thresh)/2);
    switch (variant) {
	case 1: return Double::exp(max_thresh); // return less
	default: return 0;
    }
}

inline void put_score_line(ostream& strm, double logscore, Dist d0, Dist d1, int size) {
    strm << setprecision(2) << scientific;
    strm << setw(8) << Double::exp(logscore)  << setprecision(4) << fixed << setw(8) << Double::exp(logscore/size)
	 << setw(8) << logscore << setw(8) << logscore/size << setw(8) << d0.normed(logscore) 
	 << setw(8) << d1.normed(logscore) << "\n";
} 
#endif

void Block::initDistributions() {
    ownDists.assign(size()+1, Dist());
    backDists.assign(size()+1, Dist());
    for (int i=size(); i>0; i--) {
	backDists[i-1] = backDists[i] + columns[i-1].getBackDist();
	ownDists[i-1] = ownDists[i] + columns[i-1].getOwnDist();
    }
}    

bool Block::initThresholds() {
    if (size() < MIN_BLOCKSIZE)
	return false;
    thresholdMatrix.resize(size()+1);
    for (int to=0; to <= size(); to++) {
	vector<Double>& current = thresholdMatrix[to];
	current.clear();
	for (int from=0; from<=to-MIN_CHECKCOUNT; from++) {
	    double min_logthresh = getSpecThresh(partial_spec, from, to);
	    double logthresh = getSensThresh(partial_sens, from, to);
	    if (min_logthresh <= logthresh)
		logthresh = (min_logthresh + logthresh)/2;
	    current.push_back(Double::exp(logthresh));
	}
	if (to == size()) 
	    current.resize(size()+1, almostZero);
	else { 
	    if (to < MIN_CHECKCOUNT) 
		current.push_back(almostZero);
	    current.resize(to+1, Double::infinity());
	}
    }
	    
#ifdef DEBUG
    prefixThresh.assign(MIN_CHECKCOUNT, almostZero); // here, allow everything except stop codons
    Dist d0;
    Dist& d1 = ownDists[0];
    for (int i=0; i<size(); i++) {
	if (i>=MIN_CHECKCOUNT)
	    prefixThresh.push_back(findThresh(d0, d1 - ownDists[i], partial_spec, -partial_sens, 1));
	d0 += columns[i].getBackDist();
    }
    if (abs(d0.var/backDists[0].var -1)>1e-5 || abs(d0.mu/backDists[0].mu-1)>1e-5)
	throw ProjectError("partialThresh bug");
#endif
//    blockscoreBack = backDists[0];
    double sens_thresh = getSensThresh(min_sens);
    double spec_thresh = getSpecThresh(min_spec);
#ifdef DEBUG
    Double prefixLast = findThresh(d0, d1, partial_spec, -partial_sens, 1);
    cerr << "--\nBlock: " << id << "\nThreshold calculation\n"
	 << "               score   (av.)   lgscore (av.)  sdev(bk)sdev(blk)\n"
	 << "upper bound: "; put_score_line(cerr, d1.abs(-min_sens), d0, d1, size());
    cerr << "lower bound: "; put_score_line(cerr, d0.abs(min_spec), d0, d1 ,size()); 
    cerr << "threshold:   "; 
    if (spec_thresh <= sens_thresh) {
	put_score_line(cerr, (spec_thresh + sens_thresh)/2, d0, d1 ,size()); 
// 	double dis0 = d0.normed(threshold.log());
// 	double dis1 = d1.normed(threshold.log());
// 	cerr << threshold.getRoot(size()) << " = (back)" << dis0 << "s = (block)" << dis1 << "s\n";
    } else
	cerr << "N/A\n";
    // size()-i is the length of the considered suffix and 
    // must be at least MIN_CHECKCOUNT
    for (int i=0; i<=size()-MIN_CHECKCOUNT; i++) {
	suffixThresh.push_back(findThresh(d0, ownDists[i], partial_spec, -partial_sens, 1));
	d0 -= columns[i].getBackDist();
    }
    suffixThresh.resize(size(), almostZero);
//     cerr << "Prefix thresholds:";
//     for (int i=0; i < size(); i++) 
// 	cerr << " " << prefixThresh[i];
//     cerr << "\nSuffix thresholds:";
//     for (int i=0; i < size(); i++) 
// 	cerr << " " << suffixThresh[i];
//     cerr << "\n";
    if (suffixThresh[0] != prefixLast)
	throw ProjectError("Should be equal!");
#endif    
    if (spec_thresh <= sens_thresh) {
	threshold = Double::exp((spec_thresh + sens_thresh)/2);
	return true;
    }
    return false;
}
    
void Block::setIntronProfile(const vector<string>& lines) {
    iP = new ProfileMap();
    for (int lineno=0; lineno<lines.size(); lineno++) {
	int n,f; double p;
	if (lines[lineno] == "") continue;
	if (istringstream(lines[lineno]) >> n >> f >> p && 
	    0 <= n && n < size() && 0 <= f && f < 3 && 0 < p && p <= 1)
	    iP->add(n,f,p);
	else
	    throw PartParseError(lines.size()-lineno);
    }
}

bool BlockScoreType::addBlocksUntil(bool complement, int newbase, map<int,Double> *result) {
    // add the newly created blocks to (*result)
    // return true if blocks were added
    if (newbase > DNA::len-3)
	newbase=DNA::len-3;
#ifdef DEBUG
    if (newbase < begin())
	cerr << "addBlocksUntil: value for newbase doesn't make much sense!\n";
    if (newbase < end())
	return false;
#endif
    for (int pos = end(); pos <= newbase; pos++) {
	aligned_type& new_vals = new_back();
	new_vals.reserve(size());
	Double value = 1;
	const char* dna = DNA::sequence + pos;
	int count = size();
	if (pos + 3*size() > DNA::len) count = (DNA::len-pos)/3;
	if (complement) {
	    int last = size()-count;
	    for (int i=size()-1; i>=last; i--) {
		try {
		    int aa = GeneticCode::map[Seq2Int(3).rc(dna)];
		    value *= (*partner)[i].Q(aa);
		} catch (InvalidNucleotideError e) {
		    value *= Column::invalidScore;
		}
		new_vals.push_back(value);
		dna += 3;
	    }
	} else {
	    for (int i=0; i<count; i++) {
		try {
		    int aa = GeneticCode::map[Seq2Int(3)(dna)];
		    value *= (*partner)[i].Q(aa);
		} catch (InvalidNucleotideError e) {
		    value *= Column::invalidScore;
		}
		new_vals.push_back(value);
		dna += 3; 
	    }
	}
	new_vals.resize(size(), value * Column::stopCodonScore);
	
#ifdef DEBUG
	if (count==size() && new_vals.back() != value) 
	    throw ProjectError("addBlocksUntil: something is wrong!");
#endif
//	value = new_vals.back();
	if (count==size() && value > partner->getThreshold()) {
	    hits[pos % 3][pos] = value;
	    if (result)
		(*result)[pos] = value;
	}
    }
    return true;
}

// Double Block::saveNewScore() {
//     int offset = score.begin();
//     Double val = score.get(offset); // = score.front().last()
//     score.pop_front();
//     if (val > threshold) {
// 	hits[offset % 3][offset] = val;
// 	return val;
//     } else
// 	return 0;
// }

Double Block::scoreFromScratch(bool complement, int dna_offset, int block_offset, int len) const {
    const char* dna = DNA::sequence + dna_offset;
    if (dna_offset + 3*len > DNA::len)
	return 0;
    Double result = 1;
    if (complement) {
	len = size() - len - block_offset;
	for (int i=size() - block_offset -1; i>=len; i--) {
	    try {
		int aa = GeneticCode::map[Seq2Int(3).rc(dna)];
		result *= columns[i].Q(aa);
	    } catch (InvalidNucleotideError e) {
		result *= Column::invalidScore;
	    }
	    dna +=3;
	}
    } else {
	len += block_offset;
 	for (int i=block_offset; i<len; i++) {
	    try {
		int aa = GeneticCode::map[Seq2Int(3)(dna)];
		result *= columns[i].Q(aa);
	    } catch (InvalidNucleotideError e) {
		result *= Column::invalidScore;
	    }
	    dna += 3; 
	}
    }
    return result;
}

Double Block::checkedSuffixScore(bool complement, int dna_offset, int block_offset) const {
     // if (dna_offset + 3*(size()-block_offset) > DNA::len)
    // 	return 0;
    if (block_offset >= size())
    	return 1;
    Double result = scoreFromScratch(complement, dna_offset, block_offset);
#ifdef DEBUG
    const char* dna = DNA::sequence + dna_offset;
    Double result2 = 1;
    if (dna_offset + 3*(size()-block_offset) > DNA::len)
     	result2 = 0;
    else if (complement)
    	for (int i=size() - block_offset -1; i>=0; i--) {
    	    try {
    		int aa = GeneticCode::map[Seq2Int(3).rc(dna)];
    		result2 *= columns[i].Q(aa);
    	    } catch (InvalidNucleotideError e) {}
    	    dna += 3; 
    	}
    else
    	for (int i=block_offset; i<size(); i++) {
    	    try {
    		int aa = GeneticCode::map[Seq2Int(3)(dna)];
    		result2 *= columns[i].Q(aa);
    	    } catch (InvalidNucleotideError e) {}
    	    dna += 3; 
    	}
    if (result2 != result)
	throw ProjectError("bug in Block::checkedSuffixScore");
#endif
    return result > getSuffixThresh(complement, block_offset) ? 
	result : 0;
}


void getBestPartialProduct(vector<LLDouble>& vec, PartScoreType& result) {
    result.from=0;
    result.to=0;
    int locfrom = 0;
    LLDouble globmax = 1;
    LLDouble locmax = 1;
    for (int i=0; i<vec.size(); i++) {
	locmax *= vec[i];
	if (locmax < 1) {
	    locmax = 1;
	    locfrom = i+1;
	}
	if (globmax < locmax) {
	    globmax = locmax;
	    result.from = locfrom;
	    result.to = i+1;
	}
    }
    result.score = globmax.log();
}


/*
 * bestPartialLogScore: Aligns blockstart with dna_offset and returns
 *                      the log of the best partial score. A sequence part
 *                      of MIN_BLOCKSIZE is considered a
 *                      valid part. Returns from=to, score=0 if block threshold is not
 *                      reached.
 */
void Block::bestPartialLogScore(bool complement, int dna_offset, PartScoreType& result) const {
    string aa_seq;
    const char* seq = DNA::sequence + dna_offset;
    if (complement) 
	seq += 3 * (size()-1);
    while (aa_seq.length() < size()) 
	if (complement) {
	    aa_seq += GeneticCode::revtranslate(seq);
	    seq -= 3;
	} else {
	    aa_seq += GeneticCode::translate(seq);
	    seq += 3;
	}
    int from=0;
    LLDouble locmax=1;
    LLDouble globmax=1;
    result.from=result.to=0;
    for (int t=0; t<aa_seq.length(); t++) {
	locmax *= columns[t].Q(aa_seq[t]);
	if (locmax < 1) {
	    locmax=1;
	    from = t+1;
	}
	if (globmax < locmax) {
	    globmax = locmax;
	    result.from = from;
	    result.to = t+1;
	}
    }
    if (globmax >= getPartialThresh(false, result.from, result.to) ||
	globmax >= getThreshold()) {
	result.score = globmax.log();
    } else {
	result.score = 0;
	result.from = result.to = 0;
    }
}

	
/*--- class Profile ------------------------------------------------*/

Profile::Profile(string filename) {
    ifstream strm(expandHome(filename).c_str());
    if (!strm) 
        throw ProfileNotFoundError(filename);
    try {
#ifdef DEBUG
	cerr << "Reading blocks. Here are the confidence ranges for average column scores:\n";
#endif
	parse_stream(strm);
    } catch (ProfileParseError e) {
	throw ProfileReadError(filename, e.lineno);
    } 
    if (blocks.empty()) 
	throw ProfileInsigError("No usable blocks found in file \"" + filename +"\"");
    if (blockCount() > MAX_BLOCKCOUNT)
	throw ProfileInsigError("More than " + itoa(MAX_BLOCKCOUNT) + " blocks in file \"" 
				+ filename + "\"");
    int anchorCount=0;
    for (int b=0; b<blockCount(); b++) 
	if (blocks[b].isAnchor())
	    anchorCount++;
    if (anchorCount < min_anchor_count)
	throw ProfileInsigError("At least " + itoa(min_anchor_count) + " highly "
				"significant block" + (min_anchor_count==1 ? "" : "s")
				+ "are required in profile, not found in \"" +
				filename + "\"");
		     
    if (name == "") { 
	name = filename.substr(filename.rfind("/")+1);
	name = name.substr(0, name.rfind("."));
    }
#ifdef DEBUG
//    reverse();
#endif
//    invertGlobalThreshs();
//    moved to SubstateModel
//    copyOffset = Position(blockCount()+1,0).id();
//     init();
}

vector<string> readPart(istream& strm, string& next_type, int& lineno) {
    vector<string> result;
    string line;
    next_type = "";
    while (getline(strm,line)) {
        line=line.substr(0, line.find("#"));
	line=line.substr(0, line.find_last_not_of("\t\n\v\f\r ")+1);
	if (line[0]=='[') {
	    lineno++;
	    next_type = line;
	    break;
	}
	result.push_back(line);
    }
    lineno += result.size();
    return result;
}

string readAndConcatPart(istream &strm, string& next_type, int& lineno) {
    vector<string> result = readPart(strm, next_type, lineno);
    if (result.empty()) return "";
    for (int i=1; i<result.size(); i++) 
	if (result[i].find_first_not_of("\t\n\v\f\r ") != string::npos)
	    result[0] += "\n" + result[i];
    return result[0];
}
	
int newlinesFromPos(string s, int posn) {
    int result = 0;
    while ((posn = s.find('\n',posn)+1) > 0)
	result++;
    return result;
}

void Profile::calcGlobalThresh(const vector< vector<Dist> >& ownDists) {
#ifdef DEBUG
    cerr << "Global thresholds for substate bonus:";
#endif
    Dist tail, fullDist;
    globalThresh[0].resize(blockCount());
    globalThresh[1].resize(blockCount());
    for (int b=0; b<blockCount(); b++)
	fullDist += ownDists[b][0];
    for (int b=blockCount()-1; b>=0; b--) {
#ifdef DEBUG
	cerr << "\nBlock " << b << ":";
#endif
	vector<Double>& current = globalThresh[0][b];
	vector<Double>& current_rev = globalThresh[1][blockCount() -1 -b];
	
	current.resize(blockSize(b)+1);
	current_rev.resize(blockSize(b)+1);
	for (int i=0; i<=blockSize(b); i++) {
	    Dist currDist = tail + ownDists[b][i];
	    Dist currDist_rev = fullDist - currDist;
	    current[i] = Double::exp(-currDist.abs(global_thresh));
	    current_rev[blockSize(b)-i] = 
		Double::exp(-currDist_rev.abs(global_thresh));
	    if (current[i] < absolute_malus_threshold)
		current[i] = absolute_malus_threshold;
	    if (current_rev[blockSize(b)-i] < absolute_malus_threshold)
		current_rev[blockSize(b)-i] = absolute_malus_threshold;
	}
	tail += ownDists[b][0];
#ifdef DEBUG
	cerr << scientific << "[0]: " << current[0]  << ", [" 
	     << current.size()/2  << "]: " << current[current.size()/2] 
	     << ", [" << (current.size() -1) << "]: " << current.back();
#endif
    }
#ifdef DEBUG
    cerr << "\n\n";
#endif
}

void Profile::parse_stream(istream & strm) {
    string type;
    vector<string> body; 
    int lineno = 0;  // lineno is the number of lines read 
                     // including the current line containing the type
    vector< vector<Dist> > ownDists;
    char blockName ='A';

    if (readAndConcatPart(strm, type, lineno).find_first_not_of("\t\n\v\f\r ") != string::npos) 
	throw ProfileParseError(lineno); // ignore part before first type id
    if (type == "[name]") 
	name = readAndConcatPart(strm, type, lineno);
    try {
	while (true) {
	    if (type == "[dist]") {
		// read in the allowed distance range
		istringstream lstrm(readAndConcatPart(strm, type, lineno));
		DistanceType addDist;
		if(!(lstrm >> addDist >> ws && lstrm.eof())) 
		    throw ProfileParseError(lineno - newlinesFromPos(lstrm.str(), lstrm.tellg()) -1);
		finalDist += addDist;
	    } else // if dist is not specified, assume arbitrary distance
		finalDist.setInfMax();
	    if (type != "[block]")
		break;
	    blocks.push_back(Block(finalDist, readPart(strm, type, lineno), string("block_")+(blockName++)));
	    if (type == "[intron profile]")
		blocks.back().setIntronProfile(readPart(strm, type, lineno));
	    blocks.back().initDistributions();
	    if (!blocks.back().initThresholds()) {
		// block not statistically significant, so we'll ignore it
		// and just add its size to the allowed distance range
		cerr << "Warning: Block " << blocks.back().id << " is not significant enough, removed from profile.\n";
		finalDist.r += blocks.back().size();
		blocks.pop_back();
	    } else {
		blocks.back().distance.makeTolerant();
		finalDist = DistanceType();
		ownDists.push_back(blocks.back().getAllOwnDists());
	    }
	}
    } catch (PartParseError err) {
	throw ProfileParseError(lineno - err.offset);
    }
    if (type != "")
	throw ProfileParseError(lineno);
    if (!blocks.empty())
	calcGlobalThresh(ownDists);
    finalDist.makeTolerant();
}

ostream& Profile::write(ostream& strm) const {
    strm << "[name]\n" << name << "\n";
    Position pos(0,0);
    while (true) {
	strm << "\n[dist]\n" << interBlockDist(pos.b) << "\n";
	if (pos.b >= blockCount()) break;
	while (pos.i < blockSize(pos.b)) {
	    strm << pos.i << getColumn(pos) << "\n"; pos++;
	}
	pos.nextB();
    }
    // strm << "\n[intron profile]\n" << iP << "\n";
    return strm;
}

void PP::initConstants() {
    Properties::assignProperty("/ProteinModel/block_threshold_spec", Block::min_spec);
    Properties::assignProperty("/ProteinModel/block_threshold_sens", Block::min_sens);
    Properties::assignProperty("/ProteinModel/blockpart_threshold_spec", Block::partial_spec);
    Properties::assignProperty("/ProteinModel/blockpart_threshold_sens", Block::partial_sens);  
    Properties::assignProperty("/ProteinModel/global_factor_threshold", Profile::global_thresh);
    Properties::assignProperty("/ProteinModel/absolute_malus_threshold", 
			       Profile::absolute_malus_threshold);
    Properties::assignProperty("/ProteinModel/invalid_score", Column::invalidScore);
    Properties::assignProperty("/ProteinModel/weight", Column::weight);
}



