 /**********************************************************************
 * file:    bitmasking.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  usage of bit masks to pack multiple values into a
            single integer
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 **********************************************************************/

#ifndef _BITMASKING_HH
#define _BITMASKING_HH

#include "projectio.hh"

/*
 * class SeqInKey:
 * bit mask for sequence intervals
 * encodes the start position, length, the feature type of a sequence interval in a 64 bit integer
 *   
 * bits            31                      31                 2         
 *     |-------------------------|--------------------|---------------|
 *            start position             length          feature type 
 *
 * the feature type is either "CDS" or "intron"
 */
class SeqIntKey{
     
 public:
    SeqIntKey(uint_fast64_t start, uint_fast64_t len, uint_fast64_t type){
	key = (start << 33)  // 31 bits                                                                                                            
	    + (len << 2 )    // 31 bits                                                                                      
	    + (type);        //  2 bits         

	// verification, only for debugging                      
	if(getStart() != start || getLen() != len || getType() != type)
	    throw ProjectError("internal error in SeqIntKey: packing of start=" + itoa(start)
			       + " len=" + itoa(len) + " type=" + itoa(type) + " failed.\n");
    }
    SeqIntKey(uint_fast64_t k) : key(k) {}
     ~SeqIntKey() {}
    inline uint_fast64_t getStart() const {return (key >> 33);}
    inline uint_fast64_t getLen() const {return ((key>>2) & 0x7FFFFFFF);}  // 0x7FFFFFFF bit mask for retrieving bits 1-31
    inline uint_fast64_t getType() const {return (key & 0x3);}             // 0x3 bit mask for retrieving bits 1-2
    inline uint_fast64_t getEnd() const {return (getStart()+getLen()-1);}
    inline uint_fast64_t getKey() const {return key;}
private:
     uint_fast64_t key;
 };

/*
 * class SeqPosKey:
 * bit mask for sequence positions
 * encodes the sequence ID and position in a 64 bit integer
 *   
 * bits            32                      32           
 *     |-------------------------|--------------------|
 *               seq ID                position        
 *
 */
class SeqPosKey{
     
 public:
    SeqPosKey(uint_fast64_t seqID, uint_fast64_t pos){
	key = (seqID << 32)  // 32 bits                                                                                                          
	    + (pos);         // 32 bits                                                                                                           

	// verification, only for debugging  
	if(getSeqID() != seqID || getPos() != pos )
	    throw ProjectError("internal error in SeqPosKey: packing of seqID=" + itoa(seqID)
			       + " pos=" + itoa(pos) + " failed.\n");
    }
    SeqPosKey(uint_fast64_t k) : key(k) {}
     ~SeqPosKey() {}
    inline uint_fast64_t getSeqID() const {return (key >> 32);}
    inline uint_fast64_t getPos() const {return (key & 0xFFFFFFFF);}   // 0xFFFFFFFF bit mask for retrieving bits 1-32
    inline uint_fast64_t getKey() const {return key;}
private:
     uint_fast64_t key;
 };

#endif   //  _BITMASKING_HH
