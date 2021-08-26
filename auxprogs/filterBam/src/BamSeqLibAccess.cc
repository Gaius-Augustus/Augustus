/**	read, write and access BAM files through bamtools library */

#include <iostream>
#include <vector>
#include <SeqLib/BamRecord.h>
#include <SeqLib/BamReader.h>
#include <SeqLib/BamWriter.h>
#include <htslib/hts.h>
#include "bamseqlibaccess.hh"
#include "filterBam.h"

// this workaround may be removed in a later version (>1.2) of SeqLib, see: https://github.com/walaj/SeqLib/pull/65
#include <htslib/sam.h>
inline bool LastFlagImpl(const SeqLib::BamRecord &b) {
    return (b.AlignmentFlag() & BAM_FREAD2);
}

/**
 * BamToolsAlignmentRecord constructor
 */
BamSeqLibAlignmentRecord::BamSeqLibAlignmentRecord(std::shared_ptr<SeqLib::BamRecord> alignment, std::shared_ptr<SeqLib::BamHeader> &bamHeader) {
    this->alignment = alignment;
    this->bamHeader = bamHeader;
}

/** returns the wrapped BamRecord
 */
std::shared_ptr<SeqLib::BamRecord> BamSeqLibAlignmentRecord::getAlignment() const {
    return alignment;
}

/** Retrieve the Query template NAME of the alignment, i.e the QNAME Field
 */
std::string BamSeqLibAlignmentRecord::getQueryName() const {
    return alignment->Qname();
}

/** Retrieve the length of the query sequence of the alignment, i.e the length of the SEQ Field
 */
int32_t BamSeqLibAlignmentRecord::getQuerySequenceLength() const {
    return alignment->Length();
}

/** returns position (0-based) where alignment starts
 */
int32_t BamSeqLibAlignmentRecord::getStartPosition() const {
    return alignment->Position();
}

/** returns alignment end position
 */
int32_t BamSeqLibAlignmentRecord::getEndPosition() const {
    return alignment->PositionEnd();
}

/** returns true if alignment mapped to reverse strand
 */
bool BamSeqLibAlignmentRecord::isReverseStrand() const {
    return alignment->ReverseFlag();
}

/** returns true if alignment part of paired-end read
 */
bool BamSeqLibAlignmentRecord::isPaired() const {
    return alignment->PairedFlag();
}

/** returns true if alignment is mapped
 */
bool BamSeqLibAlignmentRecord::isMapped() const {
    return alignment->MappedFlag();
}

/** returns true if alignment's mate is mapped
 */
bool BamSeqLibAlignmentRecord::isMateMapped() const {
    return alignment->MateMappedFlag();
}

/** returns true if alignment is first mate on read
 */
bool BamSeqLibAlignmentRecord::isFirstMate() const {
    return alignment->FirstFlag() && !LastFlagImpl(*alignment);
}

/** returns true if alignment is second mate on read
 */
bool BamSeqLibAlignmentRecord::isSecondMate() const {
    return !alignment->FirstFlag() && LastFlagImpl(*alignment);
}

/** returns ID number for reference sequence
 */
int32_t BamSeqLibAlignmentRecord::getRefID() const {
    return alignment->ChrID();
}

/** returns ID number for reference sequence where alignment's mate was aligned
 */
int32_t BamSeqLibAlignmentRecord::getMateRefID() const {
    return alignment->MateChrID();
}

/** returns the number of equal signs in the query sequence - occur after "samtools calmd -e" was run
 */
uint32_t BamSeqLibAlignmentRecord::countEqualSignsInQuerySequence() const {
    uint32_t numEquals = 0;

    bam1_t* myraw = alignment->raw();
    uint8_t * p = bam_get_seq(myraw);
    int32_t len = myraw->core.l_qseq;
    for (int32_t i = 0; i < len; ++i) {
        if (bam_seqi(p, i) == 0) {
            // zero is defined as '=' - see file hts.c in library htslib: seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
            // compare with file BamRecord.h in library SeqLib:
            //            static const char BASES[16] = {' ', 'A', 'C', ' ',
            //                                           'G', ' ', ' ', ' ',
            //                                           'T', ' ', ' ', ' ',
            //                                           ' ', ' ', ' ', 'N'};
            // -> thats why alignment->Sequence() doesn't contain equal signs
            numEquals++;
        }
    }

    return numEquals;
}

/** returns the sum of the total length of the M and I cigar operations.
 */
uint32_t BamSeqLibAlignmentRecord::sumMandIOperations() const {
    SeqLib::Cigar cigar = alignment->GetCigar();

    uint32_t sumMandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar) { // Length is number of bases
        // Type means operations, i.e. [MIDNSHPX=]
        if (it.Type() == 'M' || it.Type() == 'I') {
            sumMandI += it.Length();
        }
    } // end for

    return sumMandI;
}

/** returns number of insertions wrt Query and Reference through the summation of operations D and I in the CIGAR string
 */
uint32_t BamSeqLibAlignmentRecord::sumDandIOperations() const {
    SeqLib::Cigar cigar = alignment->GetCigar();

    uint32_t sumDandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar) { // Length is number of bases
        // Type means operations, i.e. [MIDNSHPX=]
        if (it.Type() == 'D' || it.Type() == 'I') {
            sumDandI += it.Length();
        }
    } // end for

    return sumDandI;
}

/** returns tag data
 *
 * @param tag_name
 * @param value return the tags value
 * @return true if tag exists and contains a valid value of values type
 */
bool BamSeqLibAlignmentRecord::getTagData(const std::string &tag_name, int32_t &value) const {
#ifdef SEQLIB_1_1_1
    value = alignment->GetIntTag(tag_name);
    return true;
#else
    return alignment->GetIntTag(tag_name, value);
#endif
}

/** returns tag data
 *
 * @param tag_name
 * @param value return the tags value
 * @return true if tag exists and contains a valid value of values type
 */
bool BamSeqLibAlignmentRecord::getTagData(const std::string &tag_name, std::string &value) const {
#ifdef SEQLIB_1_1_1
    value = alignment->GetZTag(tag_name);
    return !value.empty();
#else
    return alignment->GetZTag(tag_name, value);
#endif
}

/** add a new tag of type "Z"
 */
void BamSeqLibAlignmentRecord::addZTag(const std::string &tag_name, const std::string &value) {
    alignment->AddZTag(tag_name, value);
}

/** removes a tag
 */
void BamSeqLibAlignmentRecord::removeTag(const std::string &tag_name) {
    alignment->RemoveTag(tag_name.c_str());
}

/** Returns string with name of the reference of an alignment sequence.
 */
std::string BamSeqLibAlignmentRecord::getReferenceName() const {
#ifdef SEQLIB_1_2
    return alignment->ChrName(*bamHeader);
#else
    // workaround for https://github.com/walaj/SeqLib/issues/35
    // fixed in SeqLib 1.2.0
    int id = alignment->ChrID();
    if (id < 0) {
        return "";
    }
    return bamHeader->IDtoName(id);
#endif
}

/** opens a BAM file
 */
bool BamSeqLibReader::openReader(const std::string &filename) {
    bool opened = reader.Open(filename);
    if (opened) {
        bamHeader = std::make_shared<SeqLib::BamHeader>(reader.Header());
    }
    return opened;
}

/** retrieves next available alignment
 */
bool BamSeqLibReader::getNextAlignmentRecord(BamAlignmentRecord_ &alignment) {
    bool success = reader.GetNextRecord(bamAlignment);
    if (success) {
        alignment = std::make_shared<BamSeqLibAlignmentRecord>(std::make_shared<SeqLib::BamRecord>(bamAlignment), bamHeader);
    }
    return success;
}

/** closes the current BAM file
 */
bool BamSeqLibReader::close() {
    return reader.Close();
}

/** opens a BAM file for writing and copy headers from the specified reader
 */
bool BamSeqLibWriter::openWriter(const std::string &filename, const BamFileReader &reader) {
    const BamSeqLibReader& bslr = dynamic_cast<const BamSeqLibReader&> (reader);
    bool opened = writer.Open(filename);
    if (opened) {
        writer.SetHeader(*(bslr.bamHeader));
        writer.WriteHeader();
    }
    return opened;
}

/** saves the alignment to the alignment archive
 */
bool BamSeqLibWriter::saveAlignment(const BamAlignmentRecord_ &alignment) {
    const BamSeqLibAlignmentRecord& btar = dynamic_cast<const BamSeqLibAlignmentRecord&> (*alignment);
    return writer.WriteRecord(*(btar.getAlignment()));
}

/** closes the current BAM file
 */
void BamSeqLibWriter::close() {
    writer.Close();
}

/**
 * Sort alignments by QueryName in ascending order.
 * 
 * @param alignments
 */
void BamSeqLibUtils::sortByQueryNameAscending(std::vector<BamAlignmentRecord_> &alignments) const {
    static BamAlignmentRecordSorter ascNameSorterFunc = [](const BamAlignmentRecord_ lx, const BamAlignmentRecord_ rx) {
        const BamSeqLibAlignmentRecord& barl = dynamic_cast<const BamSeqLibAlignmentRecord&> (*lx);
        const BamSeqLibAlignmentRecord& barr = dynamic_cast<const BamSeqLibAlignmentRecord&> (*rx);
        std::less<std::string> comp;
        return comp(barl.getQueryName(), barr.getQueryName());
    };

    std::stable_sort(alignments.begin(), alignments.end(), ascNameSorterFunc);
}

/**
 * Sort alignments by Position in ascending order.
 * 
 * @param alignments
 */
void BamSeqLibUtils::sortByPositionAscending(std::vector<BamAlignmentRecord_> &alignments) const {
    static SeqLib::BamRecordSort::ByReadPosition ascPositionSorter;
    static BamAlignmentRecordSorter ascPositionSorterFunc = [](const BamAlignmentRecord_ lx, const BamAlignmentRecord_ rx) {
        const BamSeqLibAlignmentRecord& barl = dynamic_cast<const BamSeqLibAlignmentRecord&> (*lx);
        const BamSeqLibAlignmentRecord& barr = dynamic_cast<const BamSeqLibAlignmentRecord&> (*rx);
        return ascPositionSorter(*(barl.getAlignment()), *(barr.getAlignment()));
    };

    std::stable_sort(alignments.begin(), alignments.end(), ascPositionSorterFunc);
}

/**
 * Sort alignments by tag name "sc" in descending order.
 *
 * @param alignments
 */
void BamSeqLibUtils::sortByscTagDescending(std::vector<BamAlignmentRecord_> &alignments) const {
    static std::string SCORE_TAG = "sc";
    static BamAlignmentRecordSorter descScoreSorterFunc = [](const BamAlignmentRecord_ lx, const BamAlignmentRecord_ rx) {
        const BamSeqLibAlignmentRecord& barl = dynamic_cast<const BamSeqLibAlignmentRecord&> (*lx);
        const BamSeqLibAlignmentRecord& barr = dynamic_cast<const BamSeqLibAlignmentRecord&> (*rx);
        std::string lhsTagValue, rhsTagValue;
        if (!barl.getTagData(SCORE_TAG, lhsTagValue)) return true;
        if (!barr.getTagData(SCORE_TAG, rhsTagValue)) return false;

        return std::stof(lhsTagValue) > std::stof(rhsTagValue); // desc

    };

    std::stable_sort(alignments.begin(), alignments.end(), descScoreSorterFunc);
}
