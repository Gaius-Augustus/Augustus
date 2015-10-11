#ifndef _GRAPH_HH
#define _GRAPH_HH

#include <map>
#include <sstream>
#include <queue>
#include <iostream>
#include <stack>
#include "gene.hh"
#include "properties.hh"
#include "exoncand.hh"

using namespace std;

#define NUM_STATENAMES 6

enum Statename{type_unknown=-1, CDS, utr3, utr5, intron, utr3Intron, utr5Intron};

extern string stateNameIdentifiers[NUM_STATENAMES];

/* types of nodes:
 * sampled exons and introns: sampled
 * additional exons, which are not sampled: unsampled_exons
 * neutral nodes: IR, plus0, plus1, plus2, minus0, minus1, minus2 (for each of the 7 neutral lines one type)
 * NOT_KNOWN: default type, for example head and tail
 */
#define NUM_NODETYPES 27

enum NodeType{NOT_KNOWN=-1, IR,
	      plus0, plus1, plus2, minus0, minus1, minus2, T_plus1, TA_plus2, TG_plus2, T_minus1, C_minus1, YY_minus0, // intron types between two CDS exons
	      ncintr, rncintr, // intron between two non-coding exons
	      utr5intr, TLstart, TLstop, utr3intr, rutr5intr, rTLstart, rTLstop, rutr3intr, utrExon, // utr types
	      sampled, unsampled_exon}; // CDS types

extern string nodeTypeIdentifiers[NUM_NODETYPES];

class Status;
class Node;
class Edge;
class Graph;

/*
 * Status stores all the relavant information (for states) from the program specific datastructure
 */

class Status{
public:
    Status(Statename s=type_unknown, int b=0, int e=0, double sc=0.0, const void *it=NULL, Status *n=NULL):
	name(s),
	begin(b),
	end(e),
	score(sc),
	next(n),
	item(it)
    {}
    Statename name;
    int begin, end;
    double score;
    Status *next;
    const void *item;

    bool isIntron() const {return (name >= intron);}
    bool isExon() const {return (name >= CDS && name <= utr5);}
    bool isCDS() const {return (name == CDS);}
    bool isUTR() const {return (name == utr3 || name == utr5);}    
    float getPostProb() const {return (item)? ((State*)item)->apostprob : 0.0;}
    int getLen() const {return end-begin+1;}
    int getFrame() const {return isCDS()? ((State*)item)->frame() : -1;}
    bool hasEvidence() const {return ((State*)item)->evidence;}
    int numEvidence() const {return hasEvidence()? ((State*)item)->evidence->numEvidence : 0;}
};

bool isTlstartOrstop(Status *predExon, Status *succExon);

class Node{
public:
    Node(int s=0, int e=0, float sc=0.0, const void *it=NULL, NodeType t=NOT_KNOWN, Node *p=NULL, bool b=0, Node *nn=NULL, Node *pn=NULL):
	begin(s),
	end(e),
	score(sc),
	item(it),
	n_type(t),
	pred(p),
	label(b),
	nextNontrivialNeutNode(nn),
	prevNontrivialNeutNode(pn)
    {}
    int begin, end;
    float score;
    const void *item;
    NodeType n_type;
    Node *pred;
    bool label;           // label is 1, if node is in path, else label is 0
    list<Edge> edges;
    Node *nextNontrivialNeutNode;   // pointer to Node on neutral line, which is the nearest that can be reached by a nontrivial way (used for new back edges)
    Node *prevNontrivialNeutNode;   // pointer to Node on neutral line, from which the actual node is the nearest that can be reached by a nontrivial way (used for new back edges)
    size_t index;                   // (used for new back edges) (used for tarjan)
    size_t lowlink;                 // (used for tarjan)

    StateType castToStateType(); //casts void* back to State* and returns the StateType
    void addWeight(float weight); //add weight to all outgoing edges
    Edge* getEdge(Node* succ);
    State* getIntron(Node* succ);
    bool isSampled() const {return(n_type == sampled || n_type == utrExon);}

};

class Edge{
public:
    Edge(Node *t=NULL, bool n=true, float sc=0.0, const void *it=NULL):
	to(t),
	score(sc),
	neutral(n),
	item(it)
    {}
    Node *to;
    float score;
    bool neutral;
    const void *item;

    inline bool isSampledIntron(){
	return item;
    }
};

//print functions for Nodes and Edges
ostream& operator<<(ostream& ostrm, Node *node);
ostream& operator<<(ostream& ostrm, const Edge &edge);

/*
 * ordering needed by statelist: 
 * the states of the sampled genes need to be together ordered after their beginning for each gene. 
 * example: 5'UTR-CDS-intron-CDS-intron-CDS-3'UTR-5'UTR-CDS-intron-CDS-3'UTR
 * here we have 2 Genes and states in order of appearance
 */

class Graph{
public:
    Graph(list<Status> *states) : statelist(states) {}
    virtual ~Graph(); 
    void addBackEdges();
    void addBackEdgesComp();
    void tarjan();          // Tarjan's strongly connected components algorithm
    void tarjanAlg(Node* from, stack<Node*> &S, size_t &index);

    list<Node*> nodelist;      //stores all nodes belonging to the graph
    list<Status> *statelist;
    int min, max;
    map<string,Node*> existingNodes;
    Node *head;
    Node *tail;
 
    void buildGraph(); //needs to be called in constructor of derived class

    template<class T> inline bool alreadyProcessed(T *temp){
	return(existingNodes[getKey(temp)]!=NULL);    
    }
    inline bool alreadyProcessed(string key){
	return(existingNodes[key]!=NULL);    
    }
    template<class T> inline Node* getNode(T *temp){
	return existingNodes[getKey(temp)];
    }
    inline Node* getNode(string key){
	return existingNodes[key];
    } 
    // functions needed to build the graph
protected:
    bool edgeExists(Node *e1, Node *e2);
    inline void addToHash(Node *n){
	existingNodes[getKey(n)] = n;
    }
    Node* addExon(Status *exon, vector<Node*> &neutralLine);
    void addPair(Status *exon1, Status *exon2, vector<Node*> &neutralLine);
    void createNeutralLine(vector<Node*> &neutralLine, double weight=0.0, bool onlyComplete=false); 
    void addCompatibleEdges();
    void insertIntron(Node *exon1, Node *exon2);  
    int minInQueue(queue<Node*> *q);
    bool nonneutralIncomingEdge(Node *exon);
    void printGraphToShell();
    void getSizeNeutralLine();
    void addWeightToEdge();

    // program specific functions
    virtual bool exonAtGeneStart(Status *st)=0;
    virtual bool exonAtGeneEnd(Status *st)=0;
    virtual string getKey(Node *n)=0;
    virtual string getKey(Status *st)=0;
    virtual string getKey(State *st)=0;
    virtual string getKey(ExonCandidate *exoncand)=0;
    virtual double getIntronScore(Status *predExon, Status *nextExon)=0;
    virtual void addEdgeFromHead(Status *exon)=0;
    virtual void addEdgeToTail(Status *exon)=0; 
    virtual bool compatible(Node *exon1, Node *exon2)=0;
    virtual double setScore(Status *st)=0;
    virtual void calculateBaseScores()=0;
    virtual void printGraph(string filename)=0;   
    virtual void printGraph2(string filename)=0;  
    virtual bool mergedStopcodon(Node* exon1, Node* exon2)=0;
    virtual bool mergedStopcodon(Status* exon1, Status* exon2)=0;
    virtual bool mergedStopcodon(StateType type1, StateType type2, int end1, int begin2)=0;
    virtual float getAvgBaseProb(Status *st) =0;
};

class AugustusGraph : public Graph{

public:
    AugustusGraph(list<Status> *states, const char* dna) : Graph(states), sequence(dna), seqlength(strlen(dna)){
	try {
	    utr = Properties::getBoolProperty("UTR");
	} catch (...) {
	    utr = false;
	}
	try {
	    alpha_e = Properties::getdoubleProperty("/MeaPrediction/alpha_E");
	} catch (...) {
	    alpha_e = 1;
	}
	try {
	    alpha_i = Properties::getdoubleProperty("/MeaPrediction/alpha_I");
	} catch (...) {
	    alpha_i = 1;
	}
	try {
	    x0_e = Properties::getdoubleProperty("/MeaPrediction/x0_E");
	} catch (...) {
	    x0_e = -10;
	}
	try {
	    x0_i = Properties::getdoubleProperty("/MeaPrediction/x0_I");
	} catch (...) {
	    x0_i = -10;
	}
	try {
	    x1_e = Properties::getdoubleProperty("/MeaPrediction/x1_E");
	} catch (...) {
	    x1_e = 10;
	}
	try {
	    x1_i = Properties::getdoubleProperty("/MeaPrediction/x1_I");
	} catch (...) {
	    x1_i = 10;
	}
	try {
	    y0_e = Properties::getdoubleProperty("/MeaPrediction/y0_E");
	} catch (...) {
	    y0_e = 0.5;
	}
	try {
	    y0_i = Properties::getdoubleProperty("/MeaPrediction/y0_I");
	} catch (...) {
	    y0_i = 0.5;
	}
	try {
	    i1_e = Properties::getdoubleProperty("/MeaPrediction/i1_E");
	} catch (...) {
	    i1_e = 0.25;
	}	
	try {
	    i1_i = Properties::getdoubleProperty("/MeaPrediction/i1_I");
	} catch (...) {
	    i1_i = 0.25;
	}
	try {
	    i2_e = Properties::getdoubleProperty("/MeaPrediction/i2_E");
	} catch (...) {
	    i2_e = 0.75;
	}
	try {
	    i2_i = Properties::getdoubleProperty("/MeaPrediction/i2_I");
	} catch (...) {
	    i2_i = 0.75;
	}
	try {
	    j1_e = Properties::getdoubleProperty("/MeaPrediction/j1_E");
	} catch (...) {
	    j1_e = -5;
	}	
	try {
	    j1_i = Properties::getdoubleProperty("/MeaPrediction/j1_I");
	} catch (...) {
	    j1_i = -5;
	}
	try {
	    j2_e = Properties::getdoubleProperty("/MeaPrediction/j2_E");
	} catch (...) {
	    j2_e = 5;
	}
	try {
	    j2_i = Properties::getdoubleProperty("/MeaPrediction/j2_I");
	} catch (...) {
	    j2_i = 5;
	}
	try {
	    r_be = Properties::getdoubleProperty("/MeaPrediction/r_be"); // threshold for exons on base level
	} catch (...) {
	    r_be = 0.5;
	}
	try {
	    r_bi = Properties::getdoubleProperty("/MeaPrediction/r_bi"); // threshold for introns on base level
	} catch (...) {
	    r_bi = 0.5;
	}
	
	for(int i = 0; i < seqlength*10; i++)
	    baseScore.push_back(0);

    }
    bool exonAtGeneStart(Status *st);
    bool exonAtGeneEnd(Status *st);
    bool exonAtCodingStart(Node *st);
    bool exonAtCodingEnd(Node *st);
    string getKey(Node *n);
    string getKey(Status *st);
    string getKey(State *st);
    string getKey(ExonCandidate *exoncand);
    double getIntronScore(Status *predExon, Status *nextExon);  
    void addEdgeFromHead(Status *exon);
    void addEdgeToTail(Status *exon);
    bool compatible(Node *exon1, Node *exon2);
    bool sameStrand(StateType typeA, StateType typeB);
    bool sameReadingFrame(Node *e1, Node *e2);
    void calculateBaseScores();
    double setScore(Status *st);
    int getBasetype(Status *st, int pos);
    void printGraph(string filename); 
    void printGraph2(string filename);
    bool mergedStopcodon(Node* exon1, Node* exon2);
    bool mergedStopcodon(Status* exon1, Status* exon2);
    bool mergedStopcodon(StateType type1, StateType type2, int end1, int begin2);
    void getPoints(Status *st, double p, double *a1, double *a2, double *b1, double *b2);
    float getAvgBaseProb(Status *st);
    const char* sequence;
    int seqlength;
    vector<double> baseScore;
    bool utr;

  // parameters for scores
  double alpha_e;
  double alpha_i;
  double x0_e;
  double x0_i;
  double x1_e;
  double x1_i;
  double y0_e;
  double y0_i;
  double i1_e;
  double i1_i;
  double i2_e;
  double i2_i;
  double j1_e;
  double j1_i;
  double j2_e;
  double j2_i;
  double r_be;
  double r_bi;  
};

bool compareNodes(Node *first, Node *second);
bool compareEdges(Edge first, Edge second);
bool compareNodeEnds(const Node *first, const Node *second);

template <class T>
inline string to_string (const T& t)
{
    stringstream ss;
    ss << t;
    return ss.str();
}

#endif
