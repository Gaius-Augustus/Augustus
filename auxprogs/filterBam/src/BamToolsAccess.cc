/**	read, write and access BAM files through bamtools library */

#include <iostream>
#include <vector>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include "bamtoolsaccess.hh"
#include "filterBam.h"

/**
 * BamToolsAlignmentRecord constructor
 */
BamToolsAlignmentRecord::BamToolsAlignmentRecord(std::shared_ptr<BamTools::BamAlignment> alignment, std::shared_ptr<BamTools::RefVector> &refData) {
    this->alignment = alignment;
    this->refData = refData;
}

/** returns the wrapped BamAlignment
 */
std::shared_ptr<BamTools::BamAlignment> BamToolsAlignmentRecord::getAlignment() const {
    return alignment;
}

/** Retrieve the Query template NAME of the alignment, i.e the QNAME Field
 */
std::string BamToolsAlignmentRecord::getQueryName() const {
    return alignment->Name;
}

/** Retrieve the length of the query sequence of the alignment, i.e the length of the SEQ Field
 */
int32_t BamToolsAlignmentRecord::getQuerySequenceLength() const {
    return alignment->Length;
}

/** returns position (0-based) where alignment starts
 */
int32_t BamToolsAlignmentRecord::getStartPosition() const {
    return alignment->Position;
}

/** returns alignment end position
 */
int32_t BamToolsAlignmentRecord::getEndPosition() const {
    return alignment->GetEndPosition();
}

/** returns true if alignment mapped to reverse strand
 */
bool BamToolsAlignmentRecord::isReverseStrand() const {
    return alignment->IsReverseStrand();
}

/** returns true if alignment part of paired-end read
 */
bool BamToolsAlignmentRecord::isPaired() const {
    return alignment->IsPaired();
}

/** returns true if alignment is mapped
 */
bool BamToolsAlignmentRecord::isMapped() const {
    return alignment->IsMapped();
}

/** returns true if alignment's mate is mapped
 */
bool BamToolsAlignmentRecord::isMateMapped() const {
    return alignment->IsMateMapped();
}

/** returns true if alignment is first mate on read
 */
bool BamToolsAlignmentRecord::isFirstMate() const {
    return alignment->IsFirstMate();
}

/** returns true if alignment is second mate on read
 */
bool BamToolsAlignmentRecord::isSecondMate() const {
    return alignment->IsSecondMate();
}

/** returns ID number for reference sequence
 */
int32_t BamToolsAlignmentRecord::getRefID() const {
    return alignment->RefID;
}

/** returns ID number for reference sequence where alignment's mate was aligned
 */
int32_t BamToolsAlignmentRecord::getMateRefID() const {
    return alignment->MateRefID;
}

/** returns the number of equal signs in the query sequence - occur after "samtools calmd -e" was run
 */
uint32_t BamToolsAlignmentRecord::countEqualSignsInQuerySequence() const {
    std::string alignedBases = alignment->AlignedBases;
    return std::count(alignedBases.begin(), alignedBases.end(), '=');
}

/** returns the sum of the total length of the M and I cigar operations.
 */
uint32_t BamToolsAlignmentRecord::sumMandIOperations() const {
    std::vector<BamTools::CigarOp> cigar = alignment->CigarData;

    uint32_t sumMandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar) { // Length is number of bases
        // Type means operations, i.e. [MIDNSHPX=]
        if (it.Type == 'M' || it.Type == 'I') {
            sumMandI += it.Length;
        }
    } // end for

    return sumMandI;
}

/** returns number of insertions wrt Query and Reference through the summation of operations D and I in the CIGAR string
 */
uint32_t BamToolsAlignmentRecord::sumDandIOperations() const {
    std::vector<BamTools::CigarOp> cigar = alignment->CigarData;

    uint32_t sumDandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar) { // Length is number of bases
        // Type means operations, i.e. [MIDNSHPX=]
        if (it.Type == 'D' || it.Type == 'I') {
            sumDandI += it.Length;
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
bool BamToolsAlignmentRecord::getTagData(const std::string &tag_name, int32_t &value) const {
    char tag_type;
    if (!alignment->GetTagType(tag_name, tag_type)) {
        return false;
    }
    switch (tag_type) {
        case (BamTools::Constants::BAM_TAG_TYPE_INT32):
            return alignment->GetTag(tag_name, value);
        case (BamTools::Constants::BAM_TAG_TYPE_INT8):
            int8_t value8;
            if (!alignment->GetTag(tag_name, value8)) {
                return false;
            }
            value = static_cast<int32_t> (value8);
            return true;
        case (BamTools::Constants::BAM_TAG_TYPE_UINT8):
            uint8_t uvalue8;
            if (!alignment->GetTag(tag_name, uvalue8)) {
                return false;
            }
            value = static_cast<int32_t> (uvalue8);
            return true;
        case (BamTools::Constants::BAM_TAG_TYPE_INT16):
            int16_t value16;
            if (!alignment->GetTag(tag_name, value16)) {
                return false;
            }
            value = static_cast<int32_t> (value16);
            return true;
        case (BamTools::Constants::BAM_TAG_TYPE_UINT16):
            uint16_t uvalue16;
            if (!alignment->GetTag(tag_name, uvalue16)) {
                return false;
            }
            value = static_cast<int32_t> (uvalue16);
            return true;
        case (BamTools::Constants::BAM_TAG_TYPE_UINT32):
            uint32_t uvalue32;
            if (!alignment->GetTag(tag_name, uvalue32)) {
                return false;
            }
            value = static_cast<int32_t> (uvalue32);
            return true;
        default:
            return false;
    }
    return false;
}

/** returns tag data
 *
 * @param tag_name
 * @param value return the tags value
 * @return true if tag exists and contains a valid value of values type
 */
bool BamToolsAlignmentRecord::getTagData(const std::string &tag_name, std::string &value) const {
    return alignment->GetTag(tag_name, value);
}

/** add a new tag of type "Z"
 */
void BamToolsAlignmentRecord::addZTag(const std::string &tag_name, const std::string &value) {
    static std::string ztype = "Z";
    alignment->AddTag(tag_name, ztype, value);
}

/** removes a tag
 */
void BamToolsAlignmentRecord::removeTag(const std::string &tag_name) {
    alignment->RemoveTag(tag_name);
}

/** Returns string with name of the reference of an alignment sequence.
 */
std::string BamToolsAlignmentRecord::getReferenceName() const {
    int32_t refid = getRefID();
    if (refid >= 0 && refid < int32_t((*refData).size())) {
        return (*refData).at(refid).RefName;
    }
    // BAM spec: refID=-1 for a read without a mapping position
    return "[printReferenceName: invalid]";
}

/** opens a BAM file
 */
bool BamToolsReader::openReader(const std::string &filename) {
    bool opened = reader.Open(filename);
    if (opened) {
        refData = std::make_shared<BamTools::RefVector>(reader.GetReferenceData());
    }
    return opened;
}

/** retrieves next available alignment
 */
bool BamToolsReader::getNextAlignmentRecord(BamAlignmentRecord_ &alignment) {
    bool success = reader.GetNextAlignment(bamAlignment);
    if (success) {
        alignment = std::make_shared<BamToolsAlignmentRecord>(std::make_shared<BamTools::BamAlignment>(bamAlignment), refData);
    }
    return success;
}

/** closes the current BAM file
 */
bool BamToolsReader::close() {
    return reader.Close();
}

/** opens a BAM file for writing and copy headers from the specified reader
 */
bool BamToolsWriter::openWriter(const std::string &filename, const BamFileReader &reader) {
    const BamToolsReader& btr = dynamic_cast<const BamToolsReader&> (reader);
    const BamTools::SamHeader header = btr.reader.GetConstSamHeader();

    if (header.HasError()) {
        std::cout << "error in header of '" << filename << "' : " << header.GetErrorString() << std::endl;
        std::cout << "-> data in header section is ignored!" << std::endl;
    }
    return writer.Open(filename, header, *btr.refData);
}

/** saves the alignment to the alignment archive
 */
bool BamToolsWriter::saveAlignment(const BamAlignmentRecord_ &alignment) {
    const BamToolsAlignmentRecord& btar = dynamic_cast<const BamToolsAlignmentRecord&> (*alignment);
    return writer.SaveAlignment(*(btar.getAlignment()));
}

/** closes the current BAM file
 */
void BamToolsWriter::close() {
    writer.Close();
}

