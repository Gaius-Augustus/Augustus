/**	read, write and access BAM files through bamtools library */

#ifndef _BAMTOOLSACCESS_HH
#define _BAMTOOLSACCESS_HH

#include <iostream>
#include <vector>
#include <memory>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include "bamaccess.hh"
#include "filterBam.h"

class BamToolsAlignmentRecord : public BamAlignmentRecord {
private:
    std::shared_ptr<BamTools::BamAlignment> alignment;
    std::shared_ptr<BamTools::RefVector> refData;
public:
    /**
     * BamToolsAlignmentRecord constructor
     */
    BamToolsAlignmentRecord(std::shared_ptr<BamTools::BamAlignment> alignment, std::shared_ptr<BamTools::RefVector> &refData);

    /** returns the wrapped BamAlignment
     */
    std::shared_ptr<BamTools::BamAlignment> getAlignment() const;

    /** returns the Query template NAME of the alignment, i.e the QNAME Field
     */
    std::string getQueryName() const override final;

    /** returns the length of the query sequence of the alignment, i.e the length of the SEQ Field
     */
    int32_t getQuerySequenceLength() const override final;

    /** returns position (0-based) where alignment starts
     */
    int32_t getStartPosition() const override final;

    /** returns alignment end position
     */
    int32_t getEndPosition() const override final;

    /** returns true if alignment mapped to reverse strand
     */
    bool isReverseStrand() const override final;

    /** returns true if alignment part of paired-end read
     */
    bool isPaired() const override final;

    /** returns true if alignment is mapped
     */
    bool isMapped() const override final;

    /** returns true if alignment's mate is mapped
     */
    bool isMateMapped() const override final;

    /** returns true if alignment is first mate on read
     */
    bool isFirstMate() const override final;

    /** returns true if alignment is second mate on read
     */
    bool isSecondMate() const override final;

    /** returns ID number for reference sequence
     */
    int32_t getRefID() const override final;

    /** returns ID number for reference sequence where alignment's mate was aligned
     */
    int32_t getMateRefID() const override final;

    /** returns the number of equal signs in the query sequence - occur after "samtools calmd -e" was run
     */
    inline uint32_t countEqualSignsInQuerySequence() const override final;

    /** returns the sum of the total length of the M and I cigar operations.
     */
    inline uint32_t sumMandIOperations() const override final;

    /** returns number of insertions wrt Query and Reference through the summation of operations D and I in the CIGAR string
     */
    inline uint32_t sumDandIOperations() const override final;

    /** returns number of soft clippings through the summation of operations S in the CIGAR string
     */
    inline uint32_t sumSOperations() const override final;

    /** returns number of deletions through the summation of operations D in the CIGAR string
     */
    inline uint32_t sumDOperations() const override final;

    /** returns tag data
     *
     * @param tag_name
     * @param value return the tags value
     * @return true if tag exists and contains a valid value of values type
     */
    bool getTagData(const std::string &tag_name, int32_t &value) const override final;

    /** Returns string with name of the reference of an alignment sequence.
     */
    std::string getReferenceName() const override final;
};

class BamToolsWriter : public BamFileWriter {
private:
    BamTools::BamWriter writer;
public:
    /** opens a BAM file for writing and copy headers from the specified reader
     */
    bool openWriter(const std::string &filename, const BamFileReader &reader) override final;

    /** saves the alignment to the alignment archive
     */
    bool saveAlignment(const BamAlignmentRecord_ &alignment) override final;

    /** closes the current BAM file
     */
    void close() override final;
};

class BamToolsReader : public BamFileReader {
private:
    BamTools::BamReader reader;
    BamTools::BamAlignment bamAlignment;
    std::shared_ptr<BamTools::RefVector> refData;

    friend bool BamToolsWriter::openWriter(const std::string &, const BamFileReader &);

public:

    /**
     * BamToolsReader constructor
     */
    BamToolsReader(const globalOptions_t &globalOptions) {
        if (globalOptions.threads > 1) {
        cout << "WARNING: The option \"--threads=" << globalOptions.threads
                << "\" is only valid if SeqLib>=1.2 is used. "
                << "This option is ignored because BamTools is used!" << std::endl;
        }
    }
    
    /** opens a BAM file
     */
    bool openReader(const std::string &filename) override final;

    /** retrieves next available alignment
     */
    bool getNextAlignmentRecord(BamAlignmentRecord_ &alignment) override final;

    /** closes the current BAM file
     */
    bool close() override final;
};

#endif
