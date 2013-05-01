/********************************************************************************
 * 
 * pp_prepare_align.cc:  Remove sequences from a MSA (in FASTA format)
 *                       that cause insertions and deletions
 *
 * In an MSA, to produce PSSMs (position specific scoring matrices), the goal is
 * to find regions, called *blocks*, that are common to all the member sequences, 
 * without insertions and deletions, so the region can be modelled by the matrix. 
 * This program removes sequences from the alignment in order to maximize the 
 * total size of the blocks (columns in the block times number of sequences).

 * If a parameter is specified, it is considered a part of the fasta alignment whose
 * sequences must not be deleted from the alignment
 *

 */


#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <set>
#include <cstdlib> // for getenv
#include <cstdio>  // for printf (debug)
#include <iterator> // for output_iterator (debug)
#include <istream>
#include <fstream>
// #define DEBUG 1
// #define DEF_MIN_COL_COUNT 20
// #define DEF_MAX_COL_COUNT 50
// #define DEF_MAX_DEL .1
#include <iomanip>


double FULL_COL_WEIGHT = 0.8;  // a column is considered "full" if the non-gap ratio exceeds this;
                               // the environment variable PA_FULL_COL_WEIGHT sets this value

double SKIP_COL_WEIGHT = 0.2;  // a column is considered "empty" if the non-gap ratio is below this;
                               // the environment variable PA_SKIP_COL_WEIGHT sets this value

int MINSIZE = 6;               // minimum number of columns that form a block
                               // the environment variable PA_MINSIZE sets this value

int MIN_COL_COUNT = 0;         // set this to ensure a minimum number or size of blocks

using namespace std;

// map<const char*, string> seqnames;
void readDblEnv(double& dbl, const char* envname) {
    const char* envval = getenv(envname);
    cerr << envname << "=";
    if (envval) {
	dbl = atof(envval);
	cerr << dbl << "\n";
    } else {
	cerr << dbl << " (default value)\n";
    }
}

void readIntEnv(int& i, const char* envname) {
    const char* envval = getenv(envname);
    cerr << envname << "=";
    if (envval) {
	i = atoi(envval);
	cerr << i << "\n";
    } else {
	cerr << i << " (default value)\n";
    }
}

void readEnvironment() {
    readDblEnv(FULL_COL_WEIGHT, "PA_FULL_COL_WEIGHT");
    readDblEnv(SKIP_COL_WEIGHT, "PA_SKIP_COL_WEIGHT");
    readIntEnv(MINSIZE, "PA_MINSIZE");
    readIntEnv(MIN_COL_COUNT, "PA_MIN_COL_COUNT");
}

inline bool isAA(char t) {
    return ('A' <= t && t <= 'Z');
}
inline bool nongap(char t) {
    return t != '-' && t != '.';
}

    
bool clash(const map<int, set<int> >& m, int i, int j) {
    // returns true if i is not in the map or contains j
    map<int, set<int> >::const_iterator it = m.find(i);
    return it == m.end() || it->second.count(j) != 0;
}

bool cutAlignment(vector<string>& sequences, int min_size, vector<string>& nodelseqs) {
    ostream_iterator<int>  strmit(cout, "\n");
    int height=0, j=0;
    while (j < sequences.size()) {
	int netlen=0;
	for (int i=0; i<sequences[j].length(); i++)
	    if (isAA(sequences[j][i])) netlen++;
	if (netlen >= MIN_COL_COUNT) {
	    if (height<j)
		sequences[height] = sequences[j];
	    height++;
	}
	j++;
    }
    if (height < j)
	cerr << "Deleted " << (j-height) << " sequences. " << height << " sequences remaining.\n";

    sequences.resize(height);
    if (height == 0)
	return false;
    int width = sequences[0].length();
    if (width < min_size || min_size <= 0) 
	return false;

    /*
     * STEP 1: Calculate column types (depending on number of residues in each column)
     */
    double full_col_bound = FULL_COL_WEIGHT*height;
    double skip_col_bound = SKIP_COL_WEIGHT*height;

//    cerr << "PA_FULL_COL_WEIGHT=" << FULL_COL_WEIGHT << endl;
    // three types of columns:
    // '+' : "block" columns, with a ratio of AA's in column more than FULL_COL_WEIGHT
    //        and containing no lower case or unknown characters
    // '-' : "gap" columns, with a ratio of AA's in column less than SKIP_COL_WEIGHT
    // '.' : "bad" columns, with a ratio of AA's in between, or lower case amino acids

    // A "block" is a succession of columns that contains MIN_SIZE block columns
    // and no bad columns.

    string coltypes = "";
    int blocksize = 0;
    int max_net_width = 0;
    int max_bad_aa_count = 0;
    int badcount=0;

    string blockcols(width, 0);
    string coldepth(width, '.');
#ifdef DEBUG
    vector<int> col_confl_count(width, 0);
#endif
    for (int i=0; i<width; i++) {
	int aa_count = 0;
	bool inactive = false;
	for (int j=0; j<height; j++) {
	    const char& t = sequences[j][i];
	    if (isAA(t)) 
		aa_count++;
	    else if (nongap(t)) {
		inactive = true;
		aa_count++;
	    }
	}
//	cerr << aa_count << endl;
	if (aa_count < skip_col_bound) {
	    bool nodel_ok = true;
	    for (int j=0; j<nodelseqs.size(); j++) {
		const char& t = nodelseqs[j][i];
		if (nongap(t)) {
		    nodel_ok = false;
		    break;
		}
	    }
	    // aa_count below skip_col_bound: make this a gap column 
	    // but no gap columns at the border to bad columns
	    // never put gap columns where one of the nodel sequences
	    // has a residue
	    if (i>0 && coltypes[i-1] != '.' && nodel_ok) {
		coltypes.append(1,'-');
#ifdef DEBUG
		col_confl_count[i] = aa_count;
#endif
	    } else 
		coltypes.append(1,'.');
	} else {
	    for (int j=0; j<nodelseqs.size(); j++) {
		const char& t = nodelseqs[j][i];
		if (!isAA(t)) {
		    inactive = true;
		    break;
		}
	    }
	    
	    if (!inactive && aa_count >= full_col_bound) {
		// aa_count above full_col_bound: make this a block column
		// but only if it is an "aligned" column (ie, upper case)
		// never put block columns where one of the nodel sequences
		// has a gap (then inactive is true)
		coltypes.append(1,'+');
#ifdef DEBUG
		col_confl_count[i] = height - aa_count;
#endif
		blocksize++;
	    } else {
		// if the aa_count did not reach the full_col_bound
		// maximize these aa_counts to give a suggestion for
		// relaxing the bound (lower cases are excluded)
		if (!inactive) {
		    badcount++;
		    if (aa_count > max_bad_aa_count)
			max_bad_aa_count = aa_count;
		}
		coltypes.append(1,'.');
		int j;
		for (j=i-1; j>=0 && coltypes[j] == '-'; j--) {
		    coldepth[j] = coltypes[j] = '.';
#ifdef DEBUG
		    col_confl_count[j] = 0;
#endif
		}
		if (blocksize < min_size) {
		    for (; j>=0 && coltypes[j] != '.'; j--) {
			coldepth[j] = coltypes[j] = '.';
#ifdef DEBUG
			col_confl_count[j] = 0;
#endif
			if (coltypes[j]=='+') badcount++;
		    }
		} else
		    max_net_width += blocksize;
		blocksize = 0;
	    }
	}
	if (coltypes[i] != '.') 
	    coldepth[i] = (aa_count >= height) ? '*' : (aa_count*10)/height + '0';
    } // for j

    // types at right border of alignment
    if (blocksize < min_size) {
	for (int j=width-1; j>=0 && coltypes[j] != '.'; j--) 
	    coldepth[j] = coltypes[j] = '.';
    } else {
	for (int j=width-1; j>=0 && coltypes[j]=='-'; j--)
	    coldepth[j] = coltypes[j] = '.';
	max_net_width += blocksize;
    }
    if (max_net_width < MIN_COL_COUNT || max_net_width == 0) {
	cerr << "ERROR: maximal reachable width is " << max_net_width <<". Reduce PA_FULL_COL_WEIGHT." << endl;
	return false;
    } else if (badcount * 2 > max_net_width) 
	cerr << "WARNING: maximal reachable width is low (" << max_net_width << ")." << endl
	     << "Consider reducing PA_FULL_COL_WEIGHT to " << (float)max_bad_aa_count/height << "." << endl;
    if (max_net_width < MIN_COL_COUNT) 
	MIN_COL_COUNT = max_net_width;
    cerr << "Column depths are: " << endl;
    for (int i=0; i<coldepth.length(); i += 100) {
	cerr << setw(4) << right << i;
	int jmax = coldepth.length();
	if (jmax > i+100) jmax = i+100;
	for (int j=i; j<jmax; j+=10) 
	    cerr << " " << coldepth.substr(j,10);
	cerr << endl;
    }

    // Now, successively delete *conflicting* sequences, which are those 
    // that have either
    //  *   nongaps in gap columns
    //  *   non-aa's in block columns

		
    
    /*
     * STEP 2: Determine sets of sequences conflicting with potential blocks
     */

    // Sequence j belongs to blockstart_conflicts[i], if it is in conflict with
    // one of the columns of the block starting at colum i (i.e., the next MIN_SIZE 
    // block columns and the gap columns between them).
    // If at i there is a gap column or there are less than MIN_SIZE block columns until
    // the next bad column, blockstart_conflicts remains empty

    map<int, set<int> > blockstart_conflicts;
    map<int, set<int> >::iterator conf_it;
    blocksize=0; // the number of block columns from position i to the next block end
                 // (i.e. unused column or end of alignment)
    for (int i=width-1; i>=0; i--) 
	if (coltypes[i] == '+') {
	    blocksize++;
	    if (blocksize >= min_size) 
		blockstart_conflicts[i] = set<int>();
	} else if (coltypes[i] == '.') {
	    blocksize=0;
	}
    for (int j=0; j<height; j++) {
	blocksize=0; // the number of aminoacids in the current sequence from i 
	             // until the next conflict
	map<int, set<int> >::reverse_iterator rev_it = blockstart_conflicts.rbegin();
	for (int i = width-1; i>=0; i--) {
	    if (rev_it == blockstart_conflicts.rend())
		break;
	    switch (coltypes[i]) {
	      case '+' :
		if (isAA(sequences[j][i]))		    
		    // current sequence has an aminoacid at block column
		    blocksize++;
		else  // current sequence is conflicting at position i
		    blocksize=0;
		if (rev_it->first == i) {
		    if (blocksize < min_size) 
			// there is a conflict in the block
			rev_it->second.insert(j);
		    rev_it++;
		}
		break;
	      case '-' :
		// gap column: dont do anything unless current sequence 
		// is in conflict with the gap
		if (nongap(sequences[j][i]))
		    blocksize=0;
	      default:
		break;
	    }
	}
    }
//     int least = 0;
//     int least_pos = -1;
//     for (conf_it = blockstart_conflicts.begin();
// 	 conf_it != blockstart_conflicts.end();
// 	 ++conf_it) {
// 	int cand = conf_it->second.size();
// 	if (cand == 0)
// 	    continue;
// 	if (least == 0 || cand < least) {
// 	    least = cand;
// 	    least_pos = conf_it->first;
// 	}
//     }
//     cerr << "Least conflicts are " << least << " at pos " << least_pos << "." << endl;
	

    blocksize=1;
#ifdef DEBUG
    int alt_net_width=0;
#endif
    for (int i=width-1; i>=0; i--) {
	if (coltypes[i]=='+') {
	    if (blockstart_conflicts.count(i) && blockstart_conflicts[i].empty()) {
		if (blocksize >= min_size) {
#ifdef DEBUG
		    alt_net_width += min_size;
#endif
		    blocksize = min_size;
		    for (int k=i; blocksize>0; k++) 
			if (coltypes[k] == '+') 
			    blockcols[k]=blocksize--;
		} else if (blocksize==0) {
		    blockcols[i]=min_size;
#ifdef DEBUG
		    alt_net_width++;
#endif
		} 
#ifdef DEBUG
		else 
		    cerr << "ERROR: blockstart_conflicts @ impossible position " << i << endl;
#endif
	    } else
		blocksize++;
	} else if (coltypes[i]=='.') 
	    blocksize=1;
    }


    // Sequence j is conflict with blockcontents at position i
    // (i a block or a gap position), if none of the blocks containing i 
    // is compatible with it. This happens not only if the sequence is
    // in conflict with the column i itself but also if the nearest
    // conflicts in both directions are closer than the block length
    vector<string> blockcontent_conflicts(height, string(width, 0));
    
    // conflict_columnsum[i] is the number of sequences that are
    // in conflict with the blockcontent at position i:
    // conflict_columnsum[i] = # { j | blockcontent_conflicts[j][i] == 1 }
    vector<int> conflict_columnsum(width,0);

    // CAVEAT: conflict_columnsum[i] can be 0 even if the column is
    // not part of any valid block (there can be a conflict with a set of
    // columns without a conflict with any of them)

    for (int j = 0; j<height; j++) {
	int dist=min_size; // distance to the next non-conflicting blockstart
	                   // i.e., the last dist blocks had a conflict
	for (int i = 0; i<width; i++) {
	    if (coltypes[i] == '.') 
		dist = min_size;
	    else if (!clash(blockstart_conflicts, i, j))
		dist = 0;
	    else if (dist >= min_size) {
		// this sequence is compatible to none of the min_size
		// preceding blockstarts
		conflict_columnsum[i]++;
		blockcontent_conflicts[j][i] = 1;
	    }
	    if (coltypes[i] == '+')
		dist++;
	}
    }

    // we use the conflict_columnsum for weighing the sequence scores:
    // the score is the total number of conflicting blockcontents (if 
    // all other sequences were fine, this is the number of columns
    // gained by removing the sequence), divided at each column by
    // the total number of conflicted sequences.
    vector<float> seqscores(height,0);
    for (int i=0; i<width; i++) 
	if (coltypes[i]=='+') {
	    for (int j=0; j<height; j++) 
		if (blockcontent_conflicts[j][i])
		    seqscores[j] += 1.0 / conflict_columnsum[i];
	}
    
    int seq_count = height;
    while (true) {
	blocksize=0;
	int net_width = 0;
	// update net_width (the number of block columns without conflicts)
	for (int i = 0; i<width; i++) 
	    if (conflict_columnsum[i] > 0) {
#ifdef DEBUG
		if (blockcols[i]) {
		    cerr << "ERROR: blockcol mismatch at " << i << "." << endl
			 << "have conflict but blockcols[k]>0" <<endl;
		    exit(-1);
		}
#endif
		if (blocksize >= min_size) {
		    // blocksize >= min_size @ conflict: 
		    // all blocksize positions before must have blockcols[.] > 0
		    net_width += blocksize;
#ifdef DEBUG
		    for (int k = i-1; blocksize>0; k--) 
			if (coltypes[k]=='+') {
			    if (!blockcols[k]) {
				cerr << "ERROR: blockcol mismatch (k=" << k << ")." << endl;
				cerr << "blockcol == 0, but we have a block here, starting " << blocksize << " columns before," << endl;
				cerr << "ending before column " << i << endl;
				cerr << string(k, ' ') << "*\n";
				for (k=0; k<width; k++)
				    cerr << k%10;
				cerr << endl;
				for (k=0; k<width; k++)
				    cerr << (col_confl_count[k] ? "!" : " ");
				cerr << endl << coltypes << endl;
				for (k=0; k<width; k++)
				    cerr << (int)blockcols[k];
				cerr << endl;
				exit(-1);
			    }
			    blocksize--;
			}
		} else {
		    // DEBUG
		    for (int k=i-1; blocksize>0; k--)
			if (coltypes[k]=='+'){
			    if (blockcols[k]) {
				cerr << "ERROR: column " << k << " can not be part of a block." << endl
				     << " Not enough columns until the next mismatch @ " << i << "." << endl;
				for (k=0; k<width; k++)
				    cerr << k%10;
				cerr << endl;
				for (k=0; k<width; k++)
				    cerr << (col_confl_count[k] ? "!" : " ");
				cerr << endl << coltypes << endl;
				for (k=0; k<width; k++)
				    cerr << (int)blockcols[k];
				cerr << endl;
				exit(-1);
			    } else 
				blocksize--;
			}
#endif
		}
		//
		blocksize=0;
	    } else if (coltypes[i] == '+') 
		blocksize++;
	if (blocksize >= min_size) 
	    net_width += blocksize;
#ifdef DEBUG
	if (net_width != alt_net_width) {
	    cerr << "ERROR: netwidth mismatch (" << net_width << " vs " << alt_net_width << ")" << endl;
	    cerr << coltypes << endl;
	    for (int i=0; i<blockcols.length(); i++) 
		cerr << int(blockcols[i]);
	    cerr << endl;
	    exit(-1);
	}
#endif
	bool enough_cols = net_width >= MIN_COL_COUNT;
	int currsize = net_width * seq_count;
	cerr << "CURRENT SIZE:    " << net_width << " x " << seq_count << " = " << currsize << ".\n";
	if (net_width == max_net_width) {
	    cerr << "--> Done: Have reached maximal number of columns.\n";
	    break;
	}

	float maxval = enough_cols ?
	    seq_count * net_width: // current size of alignment
	    0; // current size is not accepted (we need more columns)
	int argmax = -1;
	
	// each column defines a set of sequences preventing it from
	// being a block start; we need at least to delete one of these
	// sets if we want to gain columns

	// we maximize over all the sets the estimated gain in the 
	// alignment size 
	set<int>::iterator it;
	for (conf_it = blockstart_conflicts.begin();
	     conf_it != blockstart_conflicts.end(); ++conf_it)  {
	    float added_cols = 0;
	    int i = conf_it->first;
	    set<int>& bad_seqs = conf_it->second;
	    it = bad_seqs.begin();
	    while (it != bad_seqs.end())
		if (sequences[*it] == "")
		    bad_seqs.erase(it++);
		else {
		    added_cols += seqscores[*it];
		    ++it;
		}
	    if (!bad_seqs.empty()) {
		// estimation of the new size (subtract sequences
		// subject to deletion from the current sequence count
		// and add to the net width the estimated number
		// of added columns; as mentioned above, the sum of
		// sequence scores is an estimation for the sequences gained)
		float currval = (seq_count - bad_seqs.size()) * (net_width + added_cols);
		if (currval > maxval) {
		    maxval = currval;
		    argmax = i;
		}
	    }
	}

        if (argmax == -1) {
	    cerr << "--> Done: Have reached optimal size of alignment.\n";
	    break;
	}

	// update matrix and everything
	set<int> deleted_seqs = blockstart_conflicts[argmax];
	cerr << "--> Deleting " << deleted_seqs.size() << " Sequences, making col." << argmax << " a blockstart.\n";
	seq_count -= deleted_seqs.size();
	cerr << "ESTIMATED SIZE: " << maxval/seq_count << "x" << seq_count << " = " << maxval << endl;
	

#ifdef DEBUG
	cerr << " [nodelseqs=" << nodelseqs.size() << "] ";
	cerr << "Column-wise adding scores:\n";
#endif	

	// update seqscores and blockcontent_conflicts
	for (int i = 0; i< width; i++) 
	    if (coltypes[i] == '+' || coltypes[i] == '-') {
		int& oldval = conflict_columnsum[i];
		int newval = oldval;
		for (it = deleted_seqs.begin(); it != deleted_seqs.end(); ++it) 
		    if (blockcontent_conflicts[*it][i]) {
			newval--;
			blockcontent_conflicts[*it][i] = 0;
		    }
		if (newval < oldval) {
		    if (coltypes[i] == '+') {
#ifdef DEBUG
			float x = float(oldval-newval)/oldval;
			cerr << "    col " << i << ": " << x;
#endif
			if (newval > 0) {
			    double seqscoredelta = 1.0/newval - 1.0/oldval;
			    for (int j =0; j< height; j++) {
				if (blockcontent_conflicts[j][i]) 
				    seqscores[j] += seqscoredelta;
			    }
			}
		    }
		    oldval = newval;
		}
	    }

#ifdef DEBUG
	cerr << "Columns actually made block columns: ";
	blocksize=0;
#endif
	
	for (int i=0; i<width; i++) {
	    if (coltypes[i]!='+') 
		continue;
	    conf_it = blockstart_conflicts.find(i);
	    if (conf_it ==  blockstart_conflicts.end())
		continue;
	    set<int>& conflicts=conf_it->second;
	    if (conflicts.empty())
		continue;
	    for (it = deleted_seqs.lower_bound(*conflicts.begin());
		 it != deleted_seqs.end(); ++it) {
		conflicts.erase(*it);
		if (conflicts.empty() || *it > *conflicts.rbegin())
		    break;
	    }
	    if (conflicts.empty()) {
		blocksize=min_size;
		for (int k=i; blocksize > 0; k++) 
		    if (coltypes[k] == '+') {
			if (blockcols[k]==0) {
#ifdef DEBUG
			    cerr << k << ",";
#endif
			    net_width++;
			}
			blockcols[k] = blocksize--;
		    }
	    } 
	}
	cerr << "--\n";
#ifdef DEBUG
	cerr << coltypes << endl;
	for (int i=0; i<width; i++)
	    cerr << (blockcols[i] ? '*' : '.');
	cerr << "\n--\n";
	alt_net_width = net_width;
#endif
	    if (net_width * seq_count <= currsize && enough_cols) {
	    // we dont get bigger - put the sequences back
	    cerr << "--> Done: New size " << net_width * seq_count << " is less than current size.\n";
	    seq_count += deleted_seqs.size();
	    break;
	}

 	// remove sequences
	for (it = deleted_seqs.begin(); it != deleted_seqs.end(); ++it)  {
	    seqscores[*it] = 0;
#ifdef DEBUG
	    for (int i=0; i<width; i++) 
		if (coltypes[i] != '.' && (coltypes[i] == '+') != nongap(sequences[*it][i]))
		    col_confl_count[i]--;
#endif
	    sequences[*it] = "";
	}

    }
    cerr << "Deleted " << (height - seq_count) << " sequences, " << seq_count <<  " sequences remaining.\n";
#ifdef DEBUG
    for (int j=0; j<height; j++)
	if (sequences[j] == "")
	    cerr << j << ",";
#endif
    cerr << endl;
    return true;
}
	
void read_sequences_from_stream(istream& strm, vector<string>& seqs, vector<string>& seqnames) {
    string line;
    while (getline(strm,line)) {
	if (line != "" && line[line.length()-1] == '\r') {
	    line.erase(line.length()-1);
	}
	    
	if (line != "" && line[0]=='>') {
	    seqnames.push_back(line.substr(1));
	    seqs.push_back("");
	} else if (!seqs.empty())
	    seqs.back() += line;
    }
}
    
int main(int argc, char** argv) {
    vector<string> sequences, seqnames;
    vector<string> nodelseqs;
    string buf="";

    readEnvironment();
    read_sequences_from_stream(cin, sequences, seqnames);
    
    if (argc > 1) {
	fstream strm(argv[1]);
	if (strm) 
	    read_sequences_from_stream(strm, nodelseqs, seqnames);
	else {
	    cerr << "Could not read file \"" << argv[1] << "\"." << endl;
	    exit(-1);
	}
    }
#ifdef DEBUG
//    cerr << "[" << nodelseqs.size() << "]" << endl;
#endif

    if (sequences.size() >= 2) {

    for (int i=1; i<sequences.size(); i++)
	if (sequences[i].length() != sequences[0].length()) {
	    cerr << "Sequence " << seqnames[i] << " does not have correct length ("
		 << sequences[i].length() << ", should be " << sequences[0].length() << ")." << endl;
	    exit(-1);
	}
    for (int i=0; i<nodelseqs.size(); i++)
	if (nodelseqs[i].length() != sequences[0].length()) {
	    cerr << "Sequence " << seqnames[i + sequences.size()] << " does not have correct length ("
		 << nodelseqs[i].length() << ", should be " << sequences[0].length() << ")." << endl;
	    exit(-1);
	}
    
    while (!cutAlignment(sequences, MINSIZE, nodelseqs)) {
	FULL_COL_WEIGHT -= 0.05;
	cerr << "Trying PA_FULL_COL_WEIGHT=" << FULL_COL_WEIGHT << endl;
	if (FULL_COL_WEIGHT < SKIP_COL_WEIGHT + 0.05) {
	    cerr << "ERROR: Finally giving up. Look at this file manually!\n";
	    break;
	}
    }

    }
    for (int i=0; i< sequences.size(); i++) 
	if (!sequences[i].empty()) 
	    cout << ">" << seqnames[i] << endl << sequences[i] << endl;
}
	
