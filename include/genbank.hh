/*****************************************************************************\
 * Filename : genbank.hh
 * Authors  : Emmanouil Stafilarakis, Mario Stanke
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 20.09.2001 | Stafilarakis Emm.     | Creation of the file
 * 13.02.2002 | Mario Stanke          | Enabling multiple genes per sequence
 * 05.08.2003 | Mario Stanke          | corrected rInternal frame in GBProcessor::getAnnoSequence
 * 19.09.2005 | Mario Stanke          | GBFeature and reading in UTRs
\******************************************************************************/

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

#define GBMAXLINELEN 10002

enum FileType {unknown, genbank, fasta};

/**
 * Exception class for the GenBank classes.
 *
 * author  Emmanouil Stafilarakis
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
 * GBFeature: a genbank feature entry (e.g. CDS, mRNA, TATA_signal)
 *
 * author  Mario Stanke
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
    ranges=NULL;
    begin=end=-1;
    geneid="";
    strand=plusstrand;
    complete_l=complete_r=true;
  }
  GBFeature(const char *) throw( GBError);
  bool checkRange(int len);
  bool operator<(const GBFeature &other) const{
    return (begin<other.begin || (begin==other.begin && end<other.end));
  }
  bool matches(GBFeature &other);
};


//========================================================================

/**
 * A genbank data structure with the information about a gene.
 *
 * authors Emmanouil Stafilarakis, Mario Stanke
 * see     GBSplitter
 * see     GBProcessor
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
 * A GenBank database splitter.
 *
 * author Emmanouil Stafilarakis
 * see GBProcessor
 */
class GBSplitter{
public:
    GBSplitter( string fname );
    ~GBSplitter( );
    void determineFileType();
    GBPositions* nextData( ) throw( GBError );
    AnnoSequence *getNextFASTASequence( ) throw( GBError );
    void clear() {sin.clear(); sin.str(""); ifstrm.close();}
    FileType ftype;
private:
    Boolean     findPositions( GBPositions& pos ) throw( GBError );
    Boolean     gotoEnd( );
private:
    ifstream    ifstrm;
    std::stringstream sin;
};

//========================================================================

/**
 * A GenBank data processor.
 *
 * author  Emmanouil Stafilarakis
 * see     GBSplitter
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
    GBPositions* nextPosition() throw( GBError );
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
    char*   getSequence( GBPositions& pos ) throw( GBError );
    char*   getJoin( const char* pos, Strand &strand, char *& genename ) throw( GBError );
private:
    /// The internal GenBank datafile splitter
    GBSplitter  gbs;
    int gbVerbosity;
};

#endif   //  _GENBANK_HH

