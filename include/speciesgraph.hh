/**********************************************************************
 * file:    speciesgraph.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  builds a directed acyclic graph from a set of sampled genes.
 *          underlying auxiliary structure of the graph consists of seven
 *          neutral lines each representing a type of non-coding segment
 *          In comparative gene prediction for each species an object of
 *          this class is created.
 * authors: Stefanie König
 *
 *********************************************************************/


#ifndef _SPECIESGRAPH
#define _SPECIESGRAPH

//project includes
#include "graph.hh"


class SpeciesGraph : public AugustusGraph {

private:
    AnnoSequence *seqRange;
    list<ExonCandidate*> additionalExons; //exons, which are not sampled
    string speciesname;
    Strand strand;
#ifdef DEBUG
    int count_sampled;               // number of sampled exons
    int count_additional;            // number of additional exons
    int count_overlap;               // overlap between sampled and additional exons
#endif
    double max_weight; // the max weight of a node/edge in the graph, used as an upper/lower bound
    double ec_score; //temp: until there are real scores for exon candidates
    ofstream *sampled_exons;

public:
    SpeciesGraph(list<Status> *states, AnnoSequence *seq, list<ExonCandidate*> &addEx, string name, Strand s, ofstream *se) :
	AugustusGraph(states, seq->length),
	seqRange(seq),
	additionalExons(addEx),
	speciesname(name),
	strand(s),
#ifdef DEBUG
	count_sampled(0),
	count_additional(0),
	count_overlap(0),
#endif
	max_weight(0),
	sampled_exons(se)
    {
	try {
	    ec_score = Properties::getdoubleProperty("/CompPred/ec_score");
	} catch (...) {
	    ec_score = -100;
	}
    }
    ~SpeciesGraph(){
	delete seqRange;
    }
    using AugustusGraph::getKey;

    /*
     * functions to build graph
     */
    
    void buildGraph();
    NodeType fromNeutralLine(Node *node);  // returns type of noncoding segment preceding the exon/node
    NodeType toNeutralLine(Node *node);    // returns type of noncoding segment suceeding the exon/node
    void printGraph(string filename, Node* begin, Node* end, bool only_sampled = false); // prints graph in dot-format
    inline void printGraph(string filename){
	printGraph(filename, head, tail, true);
    }
    /*
     * @getKey(): if 'node' is a neutral node, then key = PosBegin:n_type
     * else key = PosBegin:PosEnd:StateType
     */
    string getKey(Node *n);
    double setScore(Status *st);
    Node* addExon(Status *exon, vector< vector<Node*> > &neutralLines);        // adds a sampled exon to the graph
    void addExon(ExonCandidate *exon, vector< vector<Node*> > &neutralLines);  // adds an exon, which is not sampled
    void addNeutralNodes(Node *node,vector< vector<Node*> > &neutralLines);    // adds neutral nodes and edges to and from an exon
    void addIntron(Node* exon1, Node* exon2, Status *intr);                    // adds a sampled intron
    void printSampledExon(Node *node);

    inline void updateMaxWeight(double weight){
	if(abs(weight) > max_weight){
	    max_weight = abs(weight);
	}
    }
    inline double getMaxWeight() const{                         // upper bound of a maximum weight path
	return 2 * max_weight * nodelist.size();
    }
    inline string getSpeciesname() const{
	return speciesname;
    }
    inline char* getSequence() const{
	return seqRange->sequence;
    }
    inline char* getSeqID() const{
	return seqRange->seqname;
    }
    inline int getSeqOffset() const{
	return seqRange->offset;
    }
    inline int getSeqLength() const{
	return seqRange->length;
    }
    inline Strand getSeqStrand() const{
	return strand;
    }
    inline AnnoSequence* getAnnoSeq() const{
	return seqRange;
    }


    /*
     * maximum weight path problem related functions
     */

    void topSort();                                       // sorts nodelist of graph topologically
    void dfs(Node* node, map<string,Node*> &processed);   // depth first search, subroutine of topSort()                                
    void relax(Node* begin, Node *end);                   // relaxation of all nodes "in between" begin and end ("in between" in terms of the topological ordering)
    inline void relax(){
	relax(head, tail);
    }
    double localChange(MoveObject *move);                 // local path search with calculation of the score difference between new and old local path
    double getScorePath(Node *begin, Node *end);          // calc. the sum of edge weights of the path: begin ~~> end
    Node* getTopSortPred(Node *node);    //get next node in topological order which is on the path
    Node* getTopSortNext(Node *node);    //get previous node in topological order which is on the path
    Node* getNextExonOnPath(Node *node, size_t step);  // returns the i-th next exon on the current path, where i is the step size (size_t step)
    Node* getPredExonOnPath(Node *node, size_t step);  // returns the i-th preceding exon on the current path
    void printNode(Node *node); //print Node with offset

};

struct MoveNode{

    Node *node;
    double weight;

    MoveNode(Node* n, double w) : node(n), weight(w) {}
    ~MoveNode() {}
};

struct MoveEdge{

    Edge *edge;
    double weight;

    MoveEdge(Edge* e, double w) : edge(e), weight(w) {}
    ~MoveEdge() {}
};

class MoveObject{

private:
    SpeciesGraph *graph;
    size_t step_size;
    list<MoveNode> nodes;  // sorted list of nodes
    list<MoveEdge> edges;  // sorted list of edges
    Node* local_head;
    Node* local_tail;

public:
    MoveObject(SpeciesGraph *g, size_t s = 1) :
	graph(g),
	step_size(s),
	local_head(NULL),
	local_tail(NULL)
    {}
    ~MoveObject() {}

    inline Node* getHead() const{
	return local_head;
    }
    inline Node* getTail() const{
	return local_tail;
    }
    inline void addNodeBack(Node* node, double weight){
	nodes.push_back(MoveNode(node, weight));
    }
    inline void addEdgeBack(Edge* edge, double weight){
	edges.push_back(MoveEdge(edge, weight));
    }
    inline void addNodeFront(Node* node, double weight){
	nodes.push_front(MoveNode(node, weight));
    }
    inline void addEdgeFront(Edge* edge, double weight){
	edges.push_front(MoveEdge(edge, weight));
    }
    inline Node* getNodeBack() const{
	return nodes.back().node;
    }
    inline Node* getNodeFront() const{
	return nodes.front().node;
    }
    inline bool nodesIsEmpty() const{
	return nodes.empty();
    }
    inline void shiftHead(size_t step = 1){
	local_head = graph->getPredExonOnPath(local_head,step);
    }
    inline void shiftTail(size_t step = 1){
	local_tail = graph->getNextExonOnPath(local_tail,step);
    }
    void addWeights();
    void undoAddWeights();
    void initLocalHeadandTail();  //determines local_head and local_tail

};

#endif  // _SPECIESGRAPH
 
