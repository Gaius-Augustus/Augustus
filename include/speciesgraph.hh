/**********************************************************************
 * file:    speciesgraph.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  builds a directed acyclic graph from a set of sampled genes.
 *          The underlying auxiliary structure of the graph consists of seven
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

//forward declarations
class Move;

class SpeciesGraph : public AugustusGraph {

private:
    AnnoSequence *seqRange;
    list<ExonCandidate*> *additionalExons; //exons, which are not sampled
    string speciesname;
    Strand strand;
    bool genesWithoutUTRs;
    bool onlyCompleteGenes;
    double max_weight; // the max weight of a node/edge in the graph, used as an upper/lower bound
    double ec_score; //temp: until there are real scores for exon candidates
    ofstream *sampled_exons;        // output file of sampled exons
   
public:
    SpeciesGraph(list<Status> *states, AnnoSequence *seq, list<ExonCandidate*> *addEx, string name, Strand s, bool u, bool o, ofstream *se) :
	AugustusGraph(states, seq->sequence),
	seqRange(seq),
	additionalExons(addEx),
	speciesname(name),
	strand(s),
	genesWithoutUTRs(u),
	onlyCompleteGenes(o),
	max_weight(0),
	sampled_exons(se)
    {
	try {
	    ec_score = Properties::getdoubleProperty("/CompPred/ec_score");
	} catch (...) {
	    ec_score = alpha_e * x0_e - r_be;
	}
    }
    ~SpeciesGraph(){
	// AnnoSequence must be deleted by caller (it is needed at other places)
    }
    using AugustusGraph::getKey;

    /*
     * functions to build graph
     */
    
    void buildGraph();
    /*
     * @getKey(): the key of auxilary nodes is 'PosBegin:n_type'
     * the key of nodes representiong Exons is 'PosBegin:PosEnd:StateType'
     */
    string getKey(Node *n);
    void printGraph(string filename); // prints graph in dot-format
    double setScore(Status *st);
    void printSampledExon(Node *node); //prints sampled exons to display them with gBrowse
    void updateMaxWeight(double weight){
	if(abs(weight) > max_weight)
	    max_weight = abs(weight);
    }
    double getMaxWeight() const {                         // upper bound of a maximum weight path
	return 2 * max_weight * nodelist.size();
    }
    string getSpeciesname() const {return speciesname;}
    char* getSequence() const {return seqRange->sequence;}
    char* getSeqID() const {return seqRange->seqname;}
    int getSeqOffset() const {return seqRange->offset;}
    int getSeqLength() const {return seqRange->length;}
    Strand getSeqStrand() const {return strand;}
    AnnoSequence* getAnnoSeq() const {return seqRange;}

    /*
     * maximum weight path problem related functions
     */

    void topSort();                                       // sorts nodelist of graph topologically
    void dfs(Node* node);   // depth first search, subroutine of topSort()                                
    double relax(Node* begin, Node *end);                   // relaxation of all nodes "in between" begin and end ("in between" in terms of the topological ordering)
    inline double relax(){
	return relax(head, tail);
    }
    double localChange(Move *move);                 // local path search with calculation of the score difference between new and old local path
    double getScorePath(Node *begin, Node *end);          // calc. the sum of edge weights of the path: begin ~~> end
    Node* getTopSortPred(Node *node);    //get next node in topological order which is on the path
    Node* getTopSortNext(Node *node);    //get previous node in topological order which is on the path
    Node* getNextExonOnPath(Node *node, size_t step);  // returns the i-th next exon on the current path, where i is the step size (size_t step)
    Node* getPredExonOnPath(Node *node, size_t step);  // returns the i-th preceding exon on the current path
    void printNode(Node *node); //print Node with offset

private:
    /*
     * subroutines of buildGraph()
     */
    template<class T> Node* addExon(T *exon, vector< vector<Node*> > &neutralLines, map<string,Node*> &auxiliaryNodes); // add an exon to the graph
    void addIntron(Node* pred, Node* succ, Status *intr); // add an intron to the graph
    Node* addLeftSS(Status *exon, vector< vector<Node*> > &neutralLines, map<string,Node*> &auxiliaryNodes);
    Node* addRightSS(Status *exon, vector< vector<Node*> > &neutralLines, map<string,Node*> &auxiliaryNodes);
    Node* addAuxilaryNode(NodeType type, int pos, vector< vector<Node*> > &neutralLines, map<string,Node*> &auxiliaryNodes); // add an auxilary node to the graph if it does not exist already
    Node* getAuxilaryNode(NodeType type, int pos, map<string,Node*> &auxiliaryNodes) const; // get an auxilary node (returns NULL if node does not exists)
    void addAuxilaryEdge(Node *pred, Node *succ); // add an auxilary edge to the graph
    Node* addAuxNodeToLine(NodeType type, int pos, vector< vector<Node*> >&neutralLines); // add an auxilary node to one of the neutral lines if it does not exist already
    NodeType getPredType(StateType type, int begin, int end) ;  // get the type of gene feature that precedes an exon
    NodeType getSuccType(StateType type) ;    // get the type of gene feature that succeeds an exon 
    list<NodeType> getPredTypes(Node *node) ;  // same as above, but includes stop codon types 
    list<NodeType> getSuccTypes(Node *node) ;    // same as above, but includes stop codon types
    bool isGeneStart(Node *exon);
    bool isGeneEnd(Node *exon);
    Node* addNode(Status *exon); // add a sampled CDS or UTR exon
    Node* addNode(ExonCandidate *exon); // add an additional candidate exon
    Node* addNode(NodeType type, int pos); // add an auxilary node
    /*
     * subroutines of printGraph()
     */
    string getDotNodeAttributes(Node *node);
    string getDotEdgeAttributes(Node *pred, Edge *edge);
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

class Move{

private:
    SpeciesGraph *graph;
    size_t step_size;
    list<MoveNode> nodes;  // sorted list of nodes
    list<MoveEdge> edges;  // sorted list of edges
    Node* local_head;
    Node* local_tail;

public:
    Move(SpeciesGraph *g, size_t s = 1) :
	graph(g),
	step_size(s),
	local_head(NULL),
	local_tail(NULL)
    {}
    ~Move() {}

    Node* getHead() const {return local_head;}
    Node* getTail() const {return local_tail;}
    void addNodeBack(Node* node, double weight){nodes.push_back(MoveNode(node, weight));}
    void addEdgeBack(Edge* edge, double weight){edges.push_back(MoveEdge(edge, weight));}
    void addNodeFront(Node* node, double weight){nodes.push_front(MoveNode(node, weight));}
    void addEdgeFront(Edge* edge, double weight){edges.push_front(MoveEdge(edge, weight));}
    Node* getNodeBack() const {return nodes.back().node;}
    Node* getNodeFront() const {return nodes.front().node;}
    bool nodesIsEmpty() const {return nodes.empty();}
    void shiftHead(size_t step = 1){local_head = graph->getPredExonOnPath(local_head,step);}
    void shiftTail(size_t step = 1){local_tail = graph->getNextExonOnPath(local_tail,step);}
    void addWeights();
    void undoAddWeights();
    void initLocalHeadandTail();  //determines local_head and local_tail
};

#endif  // _SPECIESGRAPH
 
