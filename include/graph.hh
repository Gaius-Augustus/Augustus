#ifndef _GRAPH_HH
#define _GRAPH_HH

#include <map>
#include <sstream>
#include <queue>
#include <iostream>
#include "gene.hh"
#include "properties.hh"
#include "exoncand.hh"

using namespace std;

enum Statename{type_unknown=-1, CDS, utr3, utr5, intron, utr3Intron, utr5Intron};


/* types of nodes:
 * sampled exons and introns: sampled
 * additional exons, which are not sampled: unsampled_exons
 * neutral nodes: IR, plus0, plus1, plus2, minus0, minus1, minus2 (for each of the 7 neutral lines one type)
 * NOT_KNOWN: default type, for example head and tail
 */
#define NUM_NODETYPES 10

enum NodeType{NOT_KNOWN=-1, IR, plus0, plus1, plus2, minus0, minus1, minus2, sampled, unsampled_exon};
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
};

class Node{
public:
    Node(int s=0, int e=0, double sc=0.0, const void *it=NULL, NodeType t=NOT_KNOWN, Node *p=NULL, bool b=0, Node *n=NULL, Node *r=NULL):
	begin(s),
	end(e),
	score(sc),
	item(it),
	n_type(t),
	pred(p),
	label(b),
	topSort_next(n),
	topSort_pred(r)
    {}
    int begin, end;
    double score;
    const void *item;
    NodeType n_type;
    Node *pred;
    bool label;           // label is 1, if node is in path, else label is 0
    Node *topSort_next;   // pointer to next node in a topologically sorted list of nodes
    Node *topSort_pred;   // pointer to previous node in a topologically sorted list of nodes
    list<Edge> edges;

    StateType castToStateType(); //casts void* back to State* and returns the StateType

};

class Edge{
public:
    Edge(Node *t=NULL, bool n=true, double sc=0.0, const void *it=NULL):
	to(t),
	score(sc),
	neutral(n),
	item(it)
    {}
    Node *to;
    double score;
    bool neutral;
    const void *item;
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
    template<class T> inline Node* getNode(T *temp){
	return existingNodes[getKey(temp)];
    }

    // functions needed to build the graph
protected:
    bool edgeExists(Node *e1, Node *e2);
    inline void addToHash(Node *n){
	existingNodes[getKey(n)] = n;
    }
    Node* addExon(Status *exon, vector<Node*> &neutralLine);
    void addPair(Status *exon1, Status *exon2, vector<Node*> &neutralLine);
    void createNeutralLine(vector<Node*> &neutralLine); 
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
};

class AugustusGraph : public Graph{

public:
    AugustusGraph(list<Status> *states, int dnalength) : Graph(states), seqlength(dnalength){
	try {
	    utr = Properties::getBoolProperty("UTR");
	} catch (...) {
	    utr = false;
	}
	try {
	    alpha_se = Properties::getdoubleProperty("/MeaPrediction/alpha_SE");
	} catch (...) {
	    alpha_se = 0.745;
	}
	try {
	    alpha_si = Properties::getdoubleProperty("/MeaPrediction/alpha_SI");
	} catch (...) {
	    alpha_si = 0.3775;
	}
	try {
	    alpha_be = Properties::getdoubleProperty("/MeaPrediction/alpha_BE");
	} catch (...) {
	    alpha_be = 0.1;
	}
	try {
	    alpha_bi = Properties::getdoubleProperty("/MeaPrediction/alpha_BI");
	} catch (...) {
	    alpha_bi = 0.01;
	}
	try {
	    r_se = Properties::getdoubleProperty("/MeaPrediction/r_SE");
	} catch (...) {
	    r_se = 0.6;
	}
	try {
	    r_si = Properties::getdoubleProperty("/MeaPrediction/r_SI");
	} catch (...) {
	    r_si = 0.6;
	}
	try {
	    r_be = Properties::getdoubleProperty("/MeaPrediction/r_BE");
	} catch (...) {
	    r_be = 0.6;
	}
	try {
	    r_bi = Properties::getdoubleProperty("/MeaPrediction/r_BI");
	} catch (...) {
	    r_bi = 0.7;
	}
	try {
	    m_se = Properties::getdoubleProperty("m_SE");
	} catch (...) {
	    m_se = 1;
	}
	try {
	    m_si = Properties::getdoubleProperty("m_SI");
	} catch (...) {
	    m_si = 1;
	}
	try {
	    m_be = Properties::getdoubleProperty("m_BE");
	} catch (...) {
	    m_be = 1;
	}
	try {
	    m_bi = Properties::getdoubleProperty("m_BI");
	} catch (...) {
	  m_bi = 1;
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
    int seqlength;
    vector<double> baseScore;
    bool utr;

    // parameters for scores
    double alpha_se;
    double alpha_si;
    double alpha_be;
    double alpha_bi;
    double r_se; 
    double r_si;
    double r_be;
    double r_bi;
  double m_se;
  double m_si;
  double m_be;
  double m_bi;
};

bool compareNodes(Node *first, Node *second);
bool compareEdges(Edge first, Edge second);

template <class T>
inline string to_string (const T& t)
{
    stringstream ss;
    ss << t;
    return ss.str();
}

#endif
