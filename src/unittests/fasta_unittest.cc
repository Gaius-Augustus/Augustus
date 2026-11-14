/*
 * fasta_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for FASTA input
 */

#include "fasta.hh"
#include "genbank.hh"
#include "gtest/gtest.h"

#include <fstream>
#include <sstream>
#include <string>

namespace {
    TEST(FastaTest, ReadFile) {
	char *sequence = NULL, *name = NULL;
	long length;
	std::ifstream ifstrm("test.fa");
	if (!ifstrm){
	    FAIL() << "Could not read unittests/test.fa file.";
	} else {
	    readOneFastaSeq(ifstrm, sequence, name, length);
	    ASSERT_EQ(length, 34);
	    ASSERT_STREQ(name, "aVeryNiceGeneName");
	    ASSERT_STREQ(sequence, "atttgtttgtttgtttgctttgtttgttttgttt");
	}
    }

    // read the next record from a stringstream, return false if there is none
    bool readNext(std::stringstream &ss, std::string &name, std::string &seq) {
	char *sequence = NULL, *cname = NULL;
	long length = 0;
	readOneFastaSeq(ss, sequence, cname, length);
	if (!sequence) {
	    delete [] cname;
	    return false;
	}
	name = cname;
	seq = sequence;
	EXPECT_EQ(length, (long) seq.length());
	delete [] cname;
	delete [] sequence;
	return true;
    }

    TEST(FastaTest, IsFastaCommentLine) {
	EXPECT_TRUE(isFastaCommentLine(";"));
	EXPECT_TRUE(isFastaCommentLine("; a, b"));
	EXPECT_TRUE(isFastaCommentLine(";~schema: https://x.org"));
	EXPECT_TRUE(isFastaCommentLine("  \t;indented"));
	EXPECT_TRUE(isFastaCommentLine(";crlf\r"));
	EXPECT_FALSE(isFastaCommentLine(""));
	EXPECT_FALSE(isFastaCommentLine("   "));
	EXPECT_FALSE(isFastaCommentLine(">seq;comment"));
	EXPECT_FALSE(isFastaCommentLine("acgt;"));
    }

    TEST(FastaTest, LeadingComments) {
	std::stringstream ss(";~schema: https://example.org/fhr.json\n"
			     "; release 1.0, (draft)\n"
			     "\n"
			     "   ;indented comment\n"
			     ">seq1 description\n"
			     "acgt\n"
			     "ACGT\n");
	std::string name, seq;
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "seq1");
	EXPECT_EQ(seq, "acgtACGT");
	EXPECT_FALSE(readNext(ss, name, seq));
    }

    TEST(FastaTest, CommentsInsideAndBetweenRecords) {
	std::stringstream ss(">seq1\n"
			     "; after header: x, y\n"
			     "acgt\n"
			     "  ; inside sequence {z}\n"
			     "tttt\n"
			     "; between records\n"
			     ">seq2\n"
			     "gggg\n"
			     ";at the end\n");
	std::string name, seq;
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "seq1");
	EXPECT_EQ(seq, "acgttttt");
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "seq2");
	EXPECT_EQ(seq, "gggg");
	EXPECT_FALSE(readNext(ss, name, seq));
    }

    TEST(FastaTest, SemicolonInHeaderIsNoComment) {
	std::stringstream ss(">chr1;x=1 desc; more\nacgt\n");
	std::string name, seq;
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "chr1;x=1");
	EXPECT_EQ(seq, "acgt");
    }

    TEST(FastaTest, CRLFComments) {
	std::stringstream ss("; first, comment\r\n>seq1\r\n;inside\r\nacgt\r\nacgt\r\n");
	std::string name, seq;
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "seq1");
	EXPECT_EQ(seq, "acgtacgt");
	EXPECT_FALSE(readNext(ss, name, seq));
    }

    TEST(FastaTest, NoHeaderWithComments) {
	std::stringstream ss(";comment\nacgt\n;comment\nacgt\n");
	std::string name, seq;
	ASSERT_TRUE(readNext(ss, name, seq));
	EXPECT_EQ(name, "unnamed-1");
	EXPECT_EQ(seq, "acgtacgt");
    }

    TEST(FastaTest, ReadFileWithComments) {
	std::ifstream ifstrm("test_comments.fa");
	ASSERT_TRUE(ifstrm) << "Could not read unittests/test_comments.fa file.";
	char *sequence = NULL, *name = NULL;
	long length;
	readOneFastaSeq(ifstrm, sequence, name, length);
	ASSERT_TRUE(sequence != NULL);
	EXPECT_STREQ(name, "aVeryNiceGeneName;version=2");
	EXPECT_STREQ(sequence, "atttgtttgtttgtttgctttgtttgttttgttt");
	EXPECT_EQ(length, 34);
	delete [] name;
	delete [] sequence;
    }

    TEST(FastaTest, GBProcessorWithComments) {
	GBProcessor gbank("test_comments.fa");
	ASSERT_EQ(gbank.fileType(), fasta);
	AnnoSequence *seqlist = gbank.getSequenceList();
	ASSERT_TRUE(seqlist != NULL);
	EXPECT_STREQ(seqlist->seqname, "aVeryNiceGeneName;version=2");
	EXPECT_EQ(seqlist->length, 34);
	EXPECT_STREQ(seqlist->sequence, "atttgtttgtttgtttgctttgtttgttttgttt");
	ASSERT_TRUE(seqlist->next != NULL);
	EXPECT_STREQ(seqlist->next->seqname, "second");
	EXPECT_EQ(seqlist->next->length, 10);
	EXPECT_TRUE(seqlist->next->next == NULL);
	AnnoSequence::deleteSequence(seqlist);
    }
}
