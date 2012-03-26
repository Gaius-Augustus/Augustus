/**********************************************************************
 * file:    mea7.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  interface for MEA prediction in a graph  with seven neutral lines
 * authors: Stefanie König
 *
 **********************************************************************/

#ifndef _MEA7_HH
#define _MEA7_HH

#include "gene.hh"
#include "graph.hh"

void getMEAtranscripts7(list<Gene> *meaGenes, list<Gene> *alltranscripts, int strlength);
void getMeaGenelist7(list<Node*> meaPath, list<Gene> *meaGenes);


#endif
