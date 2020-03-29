/*
 * genbank.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _GENBANK_HH
#define _GENBANK_HH

// project includes
#include "types.hh"
#include "gene.hh"

// standard C/C++ includes
#include <list>
#include <fstream>

#ifdef ZIPINPUT
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#endif

#define GBMAXLINELEN 40000

enum FileType {unknown, genbank, fasta};

/**
 * @brief Exception class for the GenBank classes.
 *
 * @author Emmanouil Stafilarakis
 */
class GBError : public ProjectError{
public:
    /**
     * Constructor
     *
     * @param msg The describing error message.
     */
    GBError( string msg ) : ProjectError( msg ) { }
};

/**
 * @brief a genbank feature entry (e.g. CDS, mRNA, TATA_signal)
 *
 * @author Mario Stanke
 */
class GBFeature{
public:
    State *ranges;
    int begin, end;
    string geneid;
    string fkey;
    Strand strand;
    bool complete_l, complete_r; // complete at the left end, right end?
    GBFeature(){
        ranges = NULL;
        begin=end = -1;
        geneid = "";
        strand = plusstrand;
        complete_l = complete_r=true;
    }
    GBFeature(const char *);
    bool checkRange(int len);
    bool operator<(const GBFeature &other) const{
        return (begin<other.begin || (begin==other.begin && end<other.end));
    }
    bool matches(GBFeature &other);
};


//========================================================================

/**
 * @brief A genbank data structure with the information about a gene.
 *
 * @author Emmanouil Stafilarakis
 * @author Mario Stanke
 * @see     GBSplitter
 * @see     GBProcessor
 */
struct GBPositions{
  /// The entire information in GenBank format.
  char*   buffer;
  /// Pointer in 'buffer' on the "ORIGIN" position.
  char*   seqbegin;
  /// The 'buffer' length.
  int     length;
  int     seqlength;
  /// A list of pointers in 'buffer' on the "CDS" positions.
  list<char*> CDSentry;
  /// A list of pointers in 'buffer' on the "mRNA" positions.
  list<char*> mRNAentry;
  list<GBFeature> CDS;
  list<GBFeature> mRNA;
};

//========================================================================




//========================================================================

/**
 * @brief A GenBank database splitter.
 *
 * @author Emmanouil Stafilarakis
 * @see GBProcessor
 */
class GBSplitter{
public:
    GBSplitter( string fname );
    ~GBSplitter( );
    void determineFileType();
    GBPositions* nextData( );
    AnnoSequence *getNextFASTASequence( );
    void clear() {sin.clear(); sin.str(""); ifstrm.close();}
    FileType ftype;
private:
    Boolean     findPositions( GBPositions& pos );
    Boolean     gotoEnd( );
private:
    ifstream    ifstrm;
    std::stringstream sin;
};

//========================================================================

/**
 * @brief A GenBank data processor.
 *
 * @author  Emmanouil Stafilarakis
 * @see     GBSplitter
 */
class GBProcessor{
public:
    /**
     * Constructor
     */
    GBProcessor(string filename);
 
    FileType fileType() {
	return gbs.ftype;
    }

    /**
     *
     */
    GBPositions* nextPosition();
    /**
     * Get the gene information of the current data section.
     * param pos A pointer to a GBPosition object with all needed
     *              information to build a gene object.
     * return A gene object.
     * see Gene, GBPositions
     */
    Gene* getGene( GBPositions* pos );
    AnnoSequence* getAnnoSequence( GBPositions* pos );
    //Gene* getGeneList();
    AnnoSequence* getAnnoSequenceList();
    AnnoSequence* getSequenceList();
private:
    char*   getSequence( GBPositions& pos );
    char*   getJoin( const char* pos, Strand &strand, char *& genename );
private:
    /// The internal GenBank datafile splitter
    GBSplitter  gbs;
    int gbVerbosity;
};

#endif   //  _GENBANK_HH

