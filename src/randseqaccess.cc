/**********************************************************************
 * file:    randseqaccess.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  random acces to sequence data, e.g. get me chr1:1000-2000 from species 'human'
 * authors: Mario Stanke
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 07.03.12| Mario Stanke  | creation of the file
 **********************************************************************/

#include "randseqaccess.hh"
#include <iostream>

MemSeqAccess::MemSeqAccess(){
  cout << "reading in file names for species from " << Constant::speciesfilenames << endl;
  filenames = getFileNames (Constant::speciesfilenames);
  // read in sequences into memory
  // use Genbank object 
  // GBProcessor gbank(filename);
  // AnnoSequence *inSeq = gbank.getSequenceList();
  // while (inSeq) {add inSeq->sequence to map, testsequence = testsequence->next;}
}

AnnoSequence* MemSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
  cout << "retrieving " << speciesname << " " << chrName << ":" << start << "-" << end << " " << strand << " from memory." << endl;
  // use filenames map and substr(chromosome, start, end-start+1)
  // reverseComplement
  return NULL;
}

DbSeqAccess::DbSeqAccess(){
  cout << "opening database connection using connection data " << Constant::dbaccess << endl;
}

AnnoSequence* DbSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
  cout << "retrieving " << speciesname << " " << chrName << ":" << start << "-" << end << " " << strand << " from database." << endl;
  return NULL;
}

map<string,string> getFileNames (string listfile){
  map<string,string> m;
  return m;
}
