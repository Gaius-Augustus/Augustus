/**********************************************************************
 * file:    genomicMSA.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  multiple sequence alignment of genomes for comparative gene prediction
 * authors: Mario Stanke, Alexander Gebauer
 *
 *********************************************************************/

#ifndef _GENOMICMSA
#define _GENOMICMSA

#include "alignment.hh"
#include "exoncand.hh"
#include "geneMSA.hh"
#include "randseqaccess.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

#define NUMCOLNAMES 32
const string colornames[NUMCOLNAMES] = {"aquamarine", "darksalmon", "gainsboro", "gold", "cadetblue", "yellowgreen",
					"crimson", "peru", "cadetblue1", "hotpink1", "lightcyan", "magenta", "mediumseagreen",
					"mintcream", "olivedrab3", "violetred", "grey44", "peachpuff", "chartreuse3",
					"aquamarine2", "darkorange1", "forestgreen", "gray66", "khaki", "olivedrab1", 
					"skyblue4", "maroon2", "grey40", "darkturquoise", "brown4", "seagreen1", "royalblue3"};

struct AliNode {
    Alignment* a;
    int id;
    int weight;
    bool covered;
    int pred;
    int topoIdx;
};

struct AliEdge {
    int weight;
};

struct AliGraph {
    vector<size_t> topo; // topologial order
    int maxWexceptCov;
};

// typedef for an AlignmentGraph: aligments with edges for possible neighborhood in the parental alignment
typedef boost::adjacency_list<boost::setS, // disallow parallel edges
			      boost::vecS, // random access to nodes (deletion inefficient)
			      boost::bidirectionalS,
			      AliNode,
			      boost::property < boost::edge_weight_t, int, AliEdge >,
			      AliGraph> AlignmentGraph;

class AliPath {
public:
    int maxSeqRange(AlignmentGraph &g);
    list<int> path; // the nodes indices on the simple path through the alignment graph
    const MsaSignature *sig; // signature used to construct the graph (some paths differ only through sig)
    set<string> ranges;
    int weights;
    friend ostream& operator<< (ostream& strm, const AliPath &p);
};

typedef AlignmentGraph::vertex_descriptor vertex_descriptor;
typedef AlignmentGraph::edge_descriptor edge_descriptor;

class dfs_time_visitor: public boost::default_dfs_visitor {
    //    typedef typename property_traits < size_t* >::value_type T;
public:
    dfs_time_visitor(size_t *fmap, size_t n) : m_ftimemap(fmap), t(n) { }
    template < typename Vertex, typename Graph >
    void finish_vertex(Vertex u, const Graph & g) { m_ftimemap[--t] = u;}
    size_t *m_ftimemap;
    size_t t;
};



// use funcion overloading on this as the STL list cannot delete from normal and reverse iterators in the same way
void eraseListRange(list<int> L, list<int>::reverse_iterator from, list<int>::reverse_iterator to);
void eraseListRange(list<int> L, list<int>::iterator from, list<int>::iterator to);

class GenomicMSA {
public:
    GenomicMSA(RandSeqAccess *rsa_) : rsa(rsa_){}
   ~GenomicMSA(){}

    void readAlignment(string alignFilename); // reads a multiple species alignment from a *.maf file
    void printAlignment(string outFname); // print alignment in .maf format, to stdout if outFname is empty string
    int numAlignments() { return alignment.size(); }
 
    GeneMSA *getNextGene();
    /** 
     * chromosome lengths as specified in the maf file
     * for each species on hash with sequence names as keys and lengths as values
     */
    vector<map<string,int>> chrLen;
    /**
     * merges pairs of alignments in order to reduce the alignment number in trivial cases
     * without doing any potentially false mergers
     */
    void compactify();
    /**
     * changes alignment list, so that afterwards, each alignment contains a gene range:
     * a single alignment that may contain one or more gene, usually the merger of many .maf alignments
     */
    void findGeneRanges();
    /**
     * write graph of all alignments into a file so graphviz can generate a picture from it
     */
    void writeDot(AlignmentGraph const &g, string fname, MsaSignature const *superSig = NULL);

    void project(AlignmentGraph &g, const MsaSignature *sig);
    AliPath getBestConsensus(AlignmentGraph &g, const MsaSignature *sig, int &numNewCovered);
    int findBestPath(AlignmentGraph &g);
    /**
     * make several small paths from a large one that exceeds maxDNAPieceSize in one species sequence range
     */
    //void chunkyFyPaths(vector<AliPath> &allPaths, AlignmentGraph &g);
    /**
     * cut off redundand ends of path and remove paths with little extra information
     */
    bool prunePaths(vector<AliPath> &allPaths, AlignmentGraph &g);
    template< class Iterator >
    bool prunePathWrt2Other(AliPath &p, Iterator pstart, Iterator pend, 
			    AliPath &other, Iterator ostart, Iterator oend,
			    AlignmentGraph &g, bool forward);
    bool deletePathWrt2Other(AliPath &p, AliPath &other, AlignmentGraph &g);

    static int weight(const Alignment *a, const MsaSignature *sig); // node weight when projecting a to sig
    static int weight(const Alignment *a, const Alignment *b, const MsaSignature *sig); // edge weight after projection to sig

private:
    list<Alignment*> alignment;
    int numSpecies;
    RandSeqAccess *rsa; // the actual data is manages in CompGenePred
    map<string, MsaSignature> signatures;
    static int maxIntronLen;
    static int minGeneLen;
    static int maxGeneLen;
};

#endif  // _GENOMICMSA
