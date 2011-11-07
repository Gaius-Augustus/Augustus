/**********************************************************************
 * file:    exoncand.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Alexander Gebauer
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 03.11.11| Mario Stanke  | creation of the file
 **********************************************************************/

#include "exoncand.hh"
#include "intronmodel.hh"
#include "geneticcode.hh"

void getExonCands(const char* dna, float assqthresh, float dssqthresh){
  int n = strlen(dna);
  int max_exon_length = 12000;
  OpenReadingFrame orf(dna, max_exon_length, n);
  
  int frame = 0;
  int base = (n>100)? 100: n;
  int a = orf.leftmostExonBegin(frame, base, false);
  cout << "the leftmost possible begin of a forward-strand exon having frame " << frame << " at position " << base
       << " is " << a << endl;
  if (n > 30){
    bool b = onASS(dna+21);
    cout << "At position 8 on the forward strand there is an AG: " << (b? "yes": "no") << endl;
    b = onRDSS(dna+8);
    cout << "At position 8 on the reverse strand there is a GT: " << (b? "yes": "no") << endl;
  }
  // assqthresh, dssqthresh implemented later
}
