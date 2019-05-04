/*
 * lldouble_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for LLDouble
 */

#include "lldouble.hh"
#include "gtest/gtest.h"

namespace {
    TEST(LLDoubleTest, Summation) { // 1+3 = 4
	long double a = 1.0, b = 3.0, s = a+b;
	LLDouble lda(a), ldb(b);
	LLDouble lds = lda + ldb;
	EXPECT_DOUBLE_EQ(s, lds.doubleValue());
    }
    TEST(LLDoubleTest, SimpleMultiplication) { // 2*3 = 6
	EXPECT_DOUBLE_EQ(6.0, (LLDouble(2.0) * LLDouble(3.0)).doubleValue());
    }
    TEST(LLDoubleTest, ChainMultiplication) {
	// multiply by a many small factors and then divide by the same factors
	// tests implicitly testPrecision
	unsigned i, n = 10000;
	double* factors = new double[n];
	for (unsigned i = 0; i < n; ++i)
	    factors[i] = (i + 10000.0) / 45678.9;
	
	LLDouble p(1.0);
	for (i = 0; i < n; ++i)
	    p *= factors[i];
	
	for (i = 0; i < n; ++i)
	    p /= factors[i];

	EXPECT_NEAR(1.0, p.doubleValue(), 1e-8);
	delete [] factors;
    }
}
