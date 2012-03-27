/**********************************************************************
 * file:    compgenepred.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  comparative gene prediction on multiple species
 * authors: Mario Stanke
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 07.03.12| Mario Stanke  | creation of the file
 **********************************************************************/

#include "compgenepred.hh"
#include "randseqaccess.hh"

CompGenePred::CompGenePred(){
  RandSeqAccess *rsa;
  if (!Constant::speciesfilenames.empty()) {
    rsa = new MemSeqAccess();
  } else {
    rsa = new DbSeqAccess();
  }
  AnnoSequence *annoSeq = rsa->getSeq("dog", "chr31", 39798072, 39798218, minusstrand);
  if (annoSeq){
    cout << "random sequence access test "<< endl;
  }
}

void CompGenePred::start(){
}
