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
    list<ExonCandidate*> additionalExons; //exons, which are not sampled
    string speciesname;
    double max_weight; // the max weight of a node/edge in the graph, used as an upper/lower bound

public:
    SpeciesGraph(list<Status> *states, int dnalength, list<ExonCandidate*> &addEx, string name) :
	AugustusGraph(states, dnalength),
	additionalExons(addEx),
	speciesname(name),
	max_weight(0)
    {}

    using AugustusGraph::getKey;

    //functions to build graph with seven neutral lines
    
    void buildGraph();
    int fromNeutralLine(Node *node);  // determines the type of the neutral line of the preceding neutral node
    int toNeutralLine(Node *node);    // determines the type of the neutral line of the suceeding neutral node

    void printGraph(string filename, Node* begin, Node* end, bool only_sampled = false); // prints graph in dot-format
    inline void printGraph(string filename){
	printGraph(filename, head, tail, true);
    }
    string getKey(Node *n);
    double setScore(Status *st);
    Node* addExon(Status *exon, vector< vector<Node*> > &neutralLines);        // adds a sampled exon to the graph
    Node* addExon(ExonCandidate *exon, vector< vector<Node*> > &neutralLines); // adds an exon, which is not sampled
    void addNeutralNodes(Node *node,vector< vector<Node*> > &neutralLines);    // adds neutral nodes and edges to and from an exon
    void addIntron(Node* exon1, Node* exon2, Status *intr);                    // adds a sampled intron

    inline void updateMaxWeight(double weight){
	if(abs(weight) > max_weight)
	    max_weight = abs(weight);
    }
    inline double getMaxWeight(){                         // upper bound of a maximum weight path
	return 2 * max_weight * nodelist.size();
    }

    // maximum weight path problem related functions

    void topSort();                                       // sorts nodelist of graph topologically
    void dfs(Node* node, map<string,Node*> &processed);   // subroutine of topSort()                                
    void relax(Node* begin, Node *end);                   // relaxation of all nodes "in between" begin and end ("in between" in terms of the topological ordering)
    double localChange(MoveObject *move);                 // local path search with calculation of the score difference between new and old local path
    double getScorePath(Node *begin, Node *end);          // calc. the sum of edge weights of the path: begin ~~> end
    Node* getTopSortPred(Node *node, size_t step = 1); 
    Node* getTopSortNext(Node *node, size_t step = 1);

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
	local_head = graph->getTopSortPred(local_head,step);
    }
    inline void shiftTail(size_t step = 1){
	local_tail = graph->getTopSortNext(local_tail,step);
    }
    void addWeights();
    void undoAddWeights();
    void initLocalHeadandTail();

};

#endif  // _SPECIESGRAPH
