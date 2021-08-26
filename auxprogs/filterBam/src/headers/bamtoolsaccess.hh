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
    std::string getQueryName() const override;

    /** returns the length of the query sequence of the alignment, i.e the length of the SEQ Field
     */
    int32_t getQuerySequenceLength() const override;

    /** returns position (0-based) where alignment starts
     */
    int32_t getStartPosition() const override;

    /** returns alignment end position
     */
    int32_t getEndPosition() const override;

    /** returns true if alignment mapped to reverse strand
     */
    bool isReverseStrand() const override;

    /** returns true if alignment part of paired-end read
     */
    bool isPaired() const override;

    /** returns true if alignment is mapped
     */
    bool isMapped() const override;

    /** returns true if alignment's mate is mapped
     */
    bool isMateMapped() const override;

    /** returns true if alignment is first mate on read
     */
    bool isFirstMate() const override;

    /** returns true if alignment is second mate on read
     */
    bool isSecondMate() const override;

    /** returns ID number for reference sequence
     */
    int32_t getRefID() const override;

    /** returns ID number for reference sequence where alignment's mate was aligned
     */
    int32_t getMateRefID() const override;

    /** returns the number of equal signs in the query sequence - occur after "samtools calmd -e" was run
     */
    uint32_t countEqualSignsInQuerySequence() const override;

    /** returns the sum of the total length of the M and I cigar operations.
     */
    uint32_t sumMandIOperations() const override;

    /** returns number of insertions wrt Query and Reference through the summation of operations D and I in the CIGAR string
     */
    uint32_t sumDandIOperations() const override;

    /** returns tag data
     *
     * @param tag_name
     * @param value return the tags value
     * @return true if tag exists and contains a valid value of values type
     */
    bool getTagData(const std::string &tag_name, int32_t &value) const override;

    /** returns tag data
     *
     * @param tag_name
     * @param value return the tags value
     * @return true if tag exists and contains a valid value of values type
     */
    bool getTagData(const std::string &tag_name, std::string &value) const override;

    /** add a new tag of type "Z"
     */
    void addZTag(const std::string &tag_name, const std::string &value) override;

    /** removes a tag
     */
    void removeTag(const std::string &tag_name) override;

    /** Returns string with name of the reference of an alignment sequence.
     */
    std::string getReferenceName() const override;
};

class BamToolsWriter : public BamFileWriter {
private:
    BamTools::BamWriter writer;
public:
    /** opens a BAM file for writing and copy headers from the specified reader
     */
    bool openWriter(const std::string &filename, const BamFileReader &reader) override;

    /** saves the alignment to the alignment archive
     */
    bool saveAlignment(const BamAlignmentRecord_ &alignment) override;

    /** closes the current BAM file
     */
    void close() override;
};

class BamToolsReader : public BamFileReader {
private:
    BamTools::BamReader reader;
    BamTools::BamAlignment bamAlignment;
    std::shared_ptr<BamTools::RefVector> refData;

    friend bool BamToolsWriter::openWriter(const std::string &, const BamFileReader &);

public:
    /** opens a BAM file
     */
    bool openReader(const std::string &filename) override;

    /** retrieves next available alignment
     */
    bool getNextAlignmentRecord(BamAlignmentRecord_ &alignment) override;

    /** closes the current BAM file
     */
    bool close() override;
};

class BamToolsUtils : public BamUtils {
public:

    BamToolsUtils() {
    }

    /**
     * Sort alignments by QueryName in ascending order.
     *
     * @param alignments
     */
    void sortByQueryNameAscending(std::vector<BamAlignmentRecord_> &alignments) const override;

    /**
     * Sort alignments by Position in ascending order.
     *
     * @param alignments
     */
    void sortByPositionAscending(std::vector<BamAlignmentRecord_> &alignments) const override;

    /**
     * Sort alignments by tag name "sc" in descending order.
     *
     * @param alignments
     */
    void sortByscTagDescending(std::vector<BamAlignmentRecord_> &alignments) const override;
};

#endif