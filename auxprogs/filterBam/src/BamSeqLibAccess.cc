/**	read, write and access BAM files through SeqLib library */

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

    const bam1_t* myraw = alignment->raw();
    const uint8_t * p = bam_get_seq(myraw);
    const int32_t len = myraw->core.l_qseq;
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

/**
 * Count all CIGAR operation of the specified type.
 * 
 * @param type [MIDNSHPX=]
 * @return summ of all operations
*/
uint32_t BamSeqLibAlignmentRecord::countCigarOperations(const char& type) const {
    uint32_t count = 0;

    // Scanning through all CIGAR operations
    for (const auto& it : alignment->GetCigar()) {
        // Type means operations, i.e. [MIDNSHPX=]
        if (it.Type() == type) {
            count += it.Length();
        }
    }

    return count;
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
#ifdef SEQLIB_1_2
    if (threadpool) {
       reader.SetThreadPool(*threadpool);
    }
#endif
    bool opened = reader.Open(filename);
    if (opened) {
        bamHeader = std::make_shared<const SeqLib::BamHeader>(reader.Header());
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
#ifdef SEQLIB_1_2
    if (bslr.threadpool) {
        writer.SetThreadPool(*(bslr.threadpool));
    }
#endif
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
