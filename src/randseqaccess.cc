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
 * 21.03.12| Stefanie König| implementation of MemSeqAcess functions
 **********************************************************************/

#include "randseqaccess.hh"
#include "genbank.hh"
#include <iostream>
#include <fstream>

MemSeqAccess::MemSeqAccess(){
  cout << "reading in file names for species from " << Constant::speciesfilenames << endl;
  filenames = getFileNames (Constant::speciesfilenames);
  /*
   * reading in sequences into memory
   */
  for(map<string, string>::iterator it = filenames.begin(); it != filenames.end(); it++){
    GBProcessor gbank(it->second);
    AnnoSequence *inSeq = gbank.getSequenceList();
    while(inSeq){
      string key = it->first + ":" + inSeq->seqname;
      sequences[key] = inSeq->sequence;
      inSeq = inSeq->next;
    }
  }
}

AnnoSequence* MemSeqAccess::getSeq(string speciesname, string chrName, int start, int end){
  cout << "retrieving " << speciesname << " " << chrName << ":" << start << "-" << end << " " << " from memory." << endl;
  AnnoSequence *annoseq = NULL;
  string key = speciesname + ":" + chrName;
  map<string,char*>::iterator it = sequences.find(key);
  if(it != sequences.end()){
    annoseq = new AnnoSequence();
    annoseq->seqname = newstrcpy(chrName);
    annoseq->sequence = newstrcpy(it->second + start, end - start + 1);
    annoseq->length = end-start+1;
    annoseq->offset = start;
  }
  return annoseq;
}

DbSeqAccess::DbSeqAccess(){
  cout << "opening database connection using connection data " << Constant::dbaccess << endl;
}

AnnoSequence* DbSeqAccess::getSeq(string speciesname, string chrName, int start, int end){
  cout << "retrieving " << speciesname << " " << chrName << ":" << start << "-" << end << " " << " from database." << endl;
  return NULL;
}

map<string,string> getFileNames (string listfile){

  map<string,string> filenames;

  ifstream ifstrm(listfile.c_str());
  if (ifstrm.is_open()){
    string line;
    while(getline(ifstrm, line)){
      size_t pos = line.find('\t');
      if (pos != string::npos)
	filenames[line.substr(0,pos)] = line.substr(pos + 1) ;
      else
	throw ProjectError(listfile + " has wrong format in line " + line + ". correct format:\n\n" + 
			   "Homo sapiens <TAB> /dir/to/genome/human.fa\n" + 
			   "Mus musculus <TAB> /dir/to/genome/mouse.fa\n" + 
			   "...\n");
    }
    ifstrm.close();
  }
  else
    throw ProjectError("Could not open input file " + listfile);

  return filenames;
}
