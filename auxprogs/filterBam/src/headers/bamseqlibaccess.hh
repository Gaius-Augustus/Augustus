/**	read, write and access BAM files through SeqLib library */

#ifndef _BAMSEQLIBSACCESS_HH
#define _BAMSEQLIBACCESS_HH

#include <iostream>
#include <vector>
#include <SeqLib/BamRecord.h>
#include <SeqLib/BamReader.h>
#include <SeqLib/BamWriter.h>
#include <htslib/hts.h>
#include "bamaccess.hh"
#include "filterBam.h"

class BamSeqLibAlignmentRecord final : public BamAlignmentRecord {
private:
    const std::shared_ptr<SeqLib::BamRecord> alignment;
    const std::shared_ptr<const SeqLib::BamHeader> bamHeader;
public:

    /**
     * BamToolsAlignmentRecord constructor
     */
    explicit BamSeqLibAlignmentRecord(std::shared_ptr<SeqLib::BamRecord> alignment, std::shared_ptr<const SeqLib::BamHeader> bamHeader)
    : alignment(std::move(alignment)), bamHeader(std::move(bamHeader)) {
    }

    /** returns the wrapped BamRecord
     */
    std::shared_ptr<SeqLib::BamRecord> getAlignment() const;

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
    inline bool getTagData(const std::string &tag_name, int32_t &value) const override final;

    /** Returns string with name of the reference of an alignment sequence.
     */
    inline std::string getReferenceName() const override final;
};

class BamSeqLibWriter final : public BamFileWriter {
private:
    SeqLib::BamWriter writer = SeqLib::BamWriter(SeqLib::BAM);
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

class BamSeqLibReader final : public BamFileReader {
private:
    SeqLib::BamReader reader;
    SeqLib::BamRecord bamAlignment;
    std::shared_ptr<const SeqLib::BamHeader> bamHeader;
#ifdef SEQLIB_1_2
    std::shared_ptr<const SeqLib::ThreadPool> threadpool;
#endif    
    friend bool BamSeqLibWriter::openWriter(const std::string &, const BamFileReader &);

public:
    /**
     * BamSeqLibReader constructor
     */
    explicit BamSeqLibReader(const globalOptions_t &globalOptions) 
#ifdef SEQLIB_1_2
    : threadpool(std::move(std::make_shared<const SeqLib::ThreadPool>(std::max(1, globalOptions.threads))))
#endif
    {
#ifndef SEQLIB_1_2
        if (globalOptions.threads > 1) {
            cout << "The option \"--threads=" <<  globalOptions.threads 
                 << "\" is only valid if SeqLib>=1.2 is used. "
                 << "This option is ignored because an older SeqLib version is used!" << std::endl;
        }
#endif        
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
