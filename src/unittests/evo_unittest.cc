/*
 * evo_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for Evo (continuous-time Markov chains, contTimeMC.cc)
 */

#include "contTimeMC.hh"
#include "gtest/gtest.h"

namespace {
    TEST(EvoTest, ExonEvoRateMatrix) 
    {
        int states = 2;
        ExonEvo evo(states);
        evo.setLambda(0.06); // exon gain rate
        evo.setMu(0.02); // exon loss rate
        evo.setPi();
        gsl_matrix *Q = evo.getExonRateMatrix();
       
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 0, 0), -0.06);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 0, 1),  0.06);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 1, 0),  0.02);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 1, 1), -0.02);
        gsl_matrix_free(Q);
    }
    TEST(EvoTest, ExonEvo2StatesTransitionProbs) 
    {
        ExonEvo evo(2);
        EXPECT_EQ(evo.getNumStates(), 2) << "Wrong number of states.";
        vector<double> branchset {0.3, 4.};
        evo.setBranchLengths(branchset);
        evo.setLambda(0.06); // exon gain rate
        evo.setMu(0.02); // exon loss rate
        evo.setPi();
        EXPECT_DOUBLE_EQ(evo.getPi(1), 0.75); // = lambda / (lambda+mu)

        evo.getRateMatrices();
        evo.computeLogPmatrices();
        gsl_matrix *P00 = evo.getSubMatrixP(0, 2.1);
        EXPECT_NEAR(gsl_matrix_get(P00, 0, 1), 0.017785717681568013, 1e-8);
        
        gsl_matrix *P01 = evo.getSubMatrixP(0, 2.3);
        EXPECT_NEAR(gsl_matrix_get(P01, 0, 1), 0.20538822219473185, 1e-8);
        /* verified in Python3 with
           import numpy
           import scipy, scipy.linalg
           Q = numpy.array([[-.06, .06], [.02, -.02]])
           scipy.linalg.expm(Q * .3)[0,1]
           scipy.linalg.expm(Q * 4)[0,1]
        */
    }
}
