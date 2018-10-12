/*
 * meaPath.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _MEAPATH_HH
#define _MEAPATH_HH

/**
 * 
 * @author Lizzy Gerischer
 */
class MEApath{
public:
  MEApath(AugustusGraph *g):graph(g){}
  ~MEApath(){}

  void findMEApath();
  void getTopologicalOrdering();
  void dfs(Node *n);
  void relax();
  void backtracking();
  inline list<Node*> getPath(){
    return meaPath;
  }
private:
  AugustusGraph *graph;
  vector<Node*> topSort;  // reverse topological ordering of nodelist starting with tail
  map<string,Node*> processed;
  list<Node*> meaPath;
};

#endif
