
#ifndef _BAMACCESS_HH
#define _BAMACCESS_HH

#include <functional>
#include <iostream>
#include <vector>
#include <memory>

class BamAlignmentRecord;
typedef std::shared_ptr<BamAlignmentRecord> BamAlignmentRecord_;
typedef std::function<bool(const BamAlignmentRecord_, const BamAlignmentRecord_) > BamAlignmentRecordSorter;

class BamAlignmentRecord {
protected:

    BamAlignmentRecord() {} // do not create BamAlignmentRecord object
public:

    virtual ~BamAlignmentRecord() {}; // ensure the destructor of derived classes is invoked

    /** returns the Query template NAME of the alignment, i.e the QNAME Field
     */
    virtual std::string getQueryName() const = 0;

    /** returns the length of the query sequence of the alignment, i.e the length of the SEQ Field
     */
    virtual int32_t getQuerySequenceLength() const = 0;

    /** returns position (0-based) where alignment starts
     */
    virtual int32_t getStartPosition() const = 0;

    /** returns alignment end position
     */
    virtual int32_t getEndPosition() const = 0;

    /** returns true if alignment mapped to reverse strand
     */
    virtual bool isReverseStrand() const = 0;

    /** returns true if alignment part of paired-end read
     */
    virtual bool isPaired() const = 0;

    /** returns true if alignment is mapped
     */
    virtual bool isMapped() const = 0;

    /** returns true if alignment's mate is mapped
     */
    virtual bool isMateMapped() const = 0;

    /** returns true if alignment is first mate on read
     */
    virtual bool isFirstMate() const = 0;

    /** returns true if alignment is second mate on read
     */
    virtual bool isSecondMate() const = 0;

    /** returns ID number for reference sequence
     */
    virtual int32_t getRefID() const = 0;

    /** returns ID number for reference sequence where alignment's mate was aligned
     */
    virtual int32_t getMateRefID() const = 0;

    /** returns the number of equal signs in the query sequence - occur after "samtools calmd -e" was run
     */
    virtual uint32_t countEqualSignsInQuerySequence() const = 0;

    /** returns the sum of the total length of the M and I cigar operations.
     */
    virtual uint32_t sumMandIOperations() const = 0;

    /** returns number of insertions wrt Query and Reference through the summation of operations D and I in the CIGAR string
     */
    virtual uint32_t sumDandIOperations() const = 0;

    /** returns tag data
     *
     * @param tag_name
     * @param value return the tags value
     * @return true if tag exists and contains a valid value of values type
     */
    virtual bool getTagData(const std::string &tag_name, int32_t &value) const = 0;

    /** returns tag data
     *
     * @param tag_name
     * @param value return the tags value
     * @return true if tag exists and contains a valid value of values type
     */
    virtual bool getTagData(const std::string &tag_name, std::string &value) const = 0;

    /** add a new tag of type "Z"
     */
    virtual void addZTag(const std::string &tag_name, const std::string &value) = 0;

    /** removes a tag
     */
    virtual void removeTag(const std::string &tag_name) = 0;

    /** Returns string with name of the reference of an alignment sequence.
     */
    virtual std::string getReferenceName() const = 0;
};

class BamFileReader {
protected:

    BamFileReader() {} // do not create BamReader object
public:

    virtual ~BamFileReader() {}; // ensure the destructor of derived classes is invoked

    /** opens a BAM file
     */
    virtual bool openReader(const std::string &filename) = 0;

    /** retrieves next available alignment
     */
    virtual bool getNextAlignmentRecord(BamAlignmentRecord_ & alignment) = 0;

    /** closes the current BAM file
     */
    virtual bool close() = 0;
};

class BamFileWriter {
protected:

    BamFileWriter() {} // do not create BamReader object
public:

    virtual ~BamFileWriter() {}; // ensure the destructor of derived classes is invoked

    /** opens a BAM file for writing and copy headers from the specified reader
     */
    virtual bool openWriter(const std::string &filename, const BamFileReader &reader) = 0;

    /** saves the alignment to the alignment archive
     */
    virtual bool saveAlignment(const BamAlignmentRecord_ &alignment) = 0;

    /** closes the current BAM file
     */
    virtual void close() = 0;
};

class BamUtils {
protected:

    BamUtils() {} // do not create BamUtils object
public:

    virtual ~BamUtils() {}; // ensure the destructor of derived classes is invoked

    /**
     * Sort alignments by QueryName in ascending order.
     *
     * @param alignments
     */
    virtual void sortByQueryNameAscending(std::vector<BamAlignmentRecord_> &alignments) const = 0;

    /**
     * Sort alignments by Position in ascending order.
     *
     * @param alignments
     */
    virtual void sortByPositionAscending(std::vector<BamAlignmentRecord_> &alignments) const = 0;

    /**
     * Sort alignments by tag name "sc" in descending order.
     *
     * @param alignments
     */
    virtual void sortByscTagDescending(std::vector<BamAlignmentRecord_> &alignments) const = 0;
};

#endif