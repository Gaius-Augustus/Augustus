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
 * encodes the start position, length, the feature type and strand of a sequence interval in a 64 bit integer
 *   
 * bits            31                      31                 1          1 
 *     |-------------------------|--------------------|---------------|-------|
 *            start position             length          feature type  strand
 *
 * the feature type is either "CDS" or "intron"
 */
class SeqIntKey{
     
 public:
    SeqIntKey(uint_fast64_t start, uint_fast64_t len, uint_fast64_t type, uint_fast64_t strand){
	key = (start << 33)  // 31 bits                                                                                                            
	    + (len << 2 )    // 31 bits                                                                                                                                
	    + (type << 1 )   //  1 bit                                                                                                                                             
	    + (strand);      //  1 bit                                                                                                                             
                  
	// verification, only for debugging                                                                                                                        
	if(getStart() != start || getLen() != len || getType() != type || getStrand() != strand)
	    throw ProjectError("internal error in SeqIntKey: packing of start=" + itoa(start)
			       + " len=" + itoa(len) + " type=" + itoa(type) + " strand=" + itoa(strand) + "failed.\n");
    }
    SeqIntKey(uint_fast64_t k) : key(k) {}
     ~SeqIntKey() {}
    inline uint_fast64_t getStart() const {return (key >> 33);}
    inline uint_fast64_t getLen() const {return ((key>>2) & 0x7FFFFFFF);}  // 0x7FFFFFFF bit mask for retrieving bits 1-31
    inline uint_fast64_t getType() const {return ((key>>1) & 1);}
    inline uint_fast64_t getStrand() const {return (key & 1);}
    inline uint_fast64_t getEnd() const {return (getStart()+getLen()-1);}
    inline uint_fast64_t getKey() const {return key;}
private:
     uint_fast64_t key;
 };

/*
 * class SeqPosKey:
 * bit mask for sequence positions
 * encodes the sequence ID, position and strand in a 64 bit integer
 *   
 * bits            31                      32             1 
 *     |-------------------------|--------------------|-------|
 *               seq ID                position        strand
 *
 */
class SeqPosKey{
     
 public:
    SeqPosKey(uint_fast64_t seqID, uint_fast64_t pos, uint_fast64_t strand){
	key = (seqID << 33)  // 31 bits                                                                                                          
	    + (pos << 1)     // 32 bits                                                                                                           
	    + (strand);      //  1 bit       

	// verification, only for debugging  
	if(getSeqID() != seqID || getPos() != pos || getStrand() != strand)
	    throw ProjectError("internal error in SeqPosKey: packing of seqID=" + itoa(seqID)
			       + " pos=" + itoa(pos) + " strand=" + itoa(strand) + "failed.\n");
    }
    SeqPosKey(uint_fast64_t k) : key(k) {}
     ~SeqPosKey() {}
    inline uint_fast64_t getSeqID() const {return (key >> 33);}
    inline uint_fast64_t getPos() const {return ((key>>1) & 0xFFFFFFFF);}   // 0xFFFFFFFF bit mask for retrieving bits 1-32
    inline uint_fast64_t getStrand() const {return (key & 1);}
    inline uint_fast64_t getKey() const {return key;}
private:
     uint_fast64_t key;
 };

#endif   //  _BITMASKING_HH
