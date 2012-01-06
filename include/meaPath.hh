#ifndef _MEAPATH_HH
#define _MEAPATH_HH

class MEApath{
public:
  MEApath(AugustusGraph *g):graph(g){
    findMEApath();
  }
  ~MEApath(){}

  void findMEApath();
  void getTopologicalOrdering();
  void dfs(Node *n);
  void relax();
  inline list<Node*> getPath(){
    return meaPath;
  }
private:
  AugustusGraph *graph;
  vector<Node*> topSort;
  map<string,Node*> processed;
  list<Node*> meaPath;
};

#endif
