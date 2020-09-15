/*
 * alignment_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for StringAlignment (alignment.cc)
 */

#include "alignment.hh"
#include "gtest/gtest.h"


namespace {
    TEST(AlignmentTest, MsaInsertion)
    // inserting unaligned strings into a multiple sequence alignment (MSA)
    {    
        StringAlignment msa(4);
        msa.rows[0] = "a-tt-g";
        msa.rows[1] = "--ctgg";
        msa.rows[2] = "a-ttgc";
        msa.rows[3] = ""; // empty row

        list<MsaInsertion> insList;
        insList.push_back( MsaInsertion(0, 0, string("ggg")) );
        insList.push_back( MsaInsertion(1, 4, string("ttt")) );
        insList.push_back( MsaInsertion(2, 4, string("gaga")) );
        insList.push_back( MsaInsertion(3, 0, string("cc")) );

        msa.insert(insList);
        // std::cout << msa << std::endl;
        
        EXPECT_TRUE(msa.rows[0] == "ggga-tt-----g" &&
                    msa.rows[1] == "-----ctttt-gg" &&
                    msa.rows[2] == "---a-ttgagagc" &&
                    msa.rows[3] == "");
    }
    
    TEST(AlignmentTest, removeGapOnlyCols)
    // removing columns with only gaps
    {    
        StringAlignment msa(3);
        msa.rows[0] = "ggga-tt-----g";
        msa.rows[1] = "-----ctttt-gg";
        msa.rows[2] = "---a-ttgagagc";
        msa.computeLen();
        size_t n = msa.removeGapOnlyCols();
        EXPECT_TRUE(msa.rows[0] == "gggatt-----g" &&
                    msa.rows[1] == "----ctttt-gg" &&
                    msa.rows[2] == "---attgagagc");
        EXPECT_EQ(n, 1);
    }
}
