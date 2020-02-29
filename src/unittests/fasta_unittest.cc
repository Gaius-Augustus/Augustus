/*
 * fasta_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for FASTA input
 */

#include "fasta.hh"
#include "gtest/gtest.h"

#include <fstream>
#include <string>

namespace {
    TEST(FastaTest, ReadFile) {
	char *sequence = NULL, *name = NULL;
	int length;
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
}
