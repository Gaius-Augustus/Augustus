#ifndef _GRAPH_HH
#define _GRAPH_HH

#include <map>
#include <sstream>
#include <queue>

#include "gene.hh"
#include "properties.hh"

using namespace std;

enum Statename{type_unknown=-1, CDS, utr3, utr5, intron, utr3Intron, utr5Intron};
enum Neutral_type{unknown=-1, IR, plus0, plus1, plus2, minus0, minus1, minus2};

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
  Node(int s=0, int e=0, double sc=0.0, const void *it=NULL, Neutral_type t=unknown, Node *p=NULL):
    begin(s),
    end(e),
    score(sc),
    item(it),
    n_type(t),
    pred(p)
  {}
  int begin, end;
  double score;
  const void *item;
  Neutral_type n_type; // helps to identify to which neutral line a neutral Node belongs
  Node *pred;
  list<Edge> edges;
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

/*
 * ordering needed by statelist: 
 * the states of the sampled genes need to be together ordered after their beginning for each gene. 
 * example: 5'UTR-CDS-intron-CDS-intron-CDS-3'UTR-5'UTR-CDS-intron-CDS-3'UTR
 * here we have 2 Genes and states in order of appearance
 */

class Graph{
public:
  Graph(list<Status> *states) : statelist(states){}
  virtual ~Graph(); 
  void addBackEdges();

  list<Node*> nodelist;      //stores all nodes belonging to the graph
  list<Status> *statelist;   
  int min, max;
  map<string,Node*> existingNodes;
  Node *head;
  Node *tail;
 
  void buildGraph(); //needs to be called in constructor of derived class
  void buildGraph7(); // builds graph with seven neutral lines

  protected:	
  // functions needed to build the graph
  inline bool alreadyProcessed(Node *n){
  return(existingNodes[getKey(n)]!=NULL);    
  }
  inline bool alreadyProcessed(Status *st){
  return(existingNodes[getKey(st)]!=NULL); 
  }
  bool edgeExists(Node *e1, Node *e2);
  inline void addToHash(Node *n){
    existingNodes[getKey(n)] = n;
  } 
  Node* getNode(Node *n);
  Node* getNode(Status *st);
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
  virtual double getIntronScore(Status *predExon, Status *nextExon)=0;
  virtual void addEdgeFromHead(Status *exon)=0;
  virtual void addEdgeToTail(Status *exon)=0; 
  virtual bool compatible(Node *exon1, Node *exon2)=0;
  virtual double setScore(Status *st)=0;
  virtual void calculateBaseScores()=0;
  virtual void printGraph(string filename)=0;   

  // additional functions to construct graph with 7 neutral lines
  Node* addExon(Status *exon, vector< vector<Node*> > &neutralLines);
  void addIntron(Node* exon1, Node* exon2, Status *intr);
  virtual int fromNeutralLine(Status *st)=0;
  virtual int toNeutralLine(Status *st)=0;
  virtual void printGraph7(string filename)=0;

};

class AugustusGraph : public Graph{

public:
  AugustusGraph(list<Status> *states, int dnalength) : Graph(states),seqlength(dnalength){
  try {
    utr = Properties::getBoolProperty("UTR");
  } catch (...) {
    utr = false;
  }
  try {
    alpha_se = Properties::getdoubleProperty("/MeaPrediction/alpha_SE");
  } catch (...) {
    alpha_se = 0;
  }
  try {
    alpha_si = Properties::getdoubleProperty("/MeaPrediction/alpha_SI");
  } catch (...) {
    alpha_si = 0;
  }
  try {
    alpha_be = Properties::getdoubleProperty("/MeaPrediction/alpha_BE");
  } catch (...) {
    alpha_be = 0;
  }
  try {
    alpha_bi = Properties::getdoubleProperty("/MeaPrediction/alpha_BI");
  } catch (...) {
    alpha_bi = 0;
  }
  try {
    r_se = Properties::getdoubleProperty("/MeaPrediction/r_SE");
  } catch (...) {
    r_se = 0;
  }
  try {
    r_si = Properties::getdoubleProperty("/MeaPrediction/r_SI");
  } catch (...) {
    r_si = 0;
  }
try {
    r_be = Properties::getdoubleProperty("/MeaPrediction/r_BE");
  } catch (...) {
    r_be = 0;
  }
 try {
    r_bi = Properties::getdoubleProperty("/MeaPrediction/r_BI");
  } catch (...) {
    r_bi = 0;
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
  int seqlength;
  vector<double> baseScore;
  bool utr;

  int fromNeutralLine(Status *st);  
  int toNeutralLine(Status *st);
  void printGraph7(string filename);

  // parameters for scores
  double alpha_se;
  double alpha_si;
  double alpha_be;
  double alpha_bi;
  double r_se; 
  double r_si;
  double r_be;
  double r_bi;
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
