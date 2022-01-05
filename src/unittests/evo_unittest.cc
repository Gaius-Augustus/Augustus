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
        double lambda = 0.06;
        double mu = 0.02;
        double ali_error = 0.1;
        ExonEvo evo(states, lambda, mu, ali_error);
        gsl_matrix *Q = evo.getExonRateMatrix();
       
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 0, 0), -0.06);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 0, 1),  0.06);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 1, 0),  0.02);
        EXPECT_DOUBLE_EQ(gsl_matrix_get (Q, 1, 1), -0.02);
        gsl_matrix_free(Q);
    }
    TEST(EvoTest, ExonEvo2StatesTransitionProbs) 
    {
        int states = 2;
        double lambda = 0.06;
        double mu = 0.02;
        double ali_error = 0.1;
        ExonEvo evo(states, lambda, mu, ali_error);
        EXPECT_EQ(evo.getNumStates(), 2) << "Wrong number of states.";
        vector<double> branchset {0.3, 4.};
        evo.setBranchLengths(branchset);
        evo.setPi();
        EXPECT_DOUBLE_EQ(evo.getPi(1), lambda/(lambda+mu));

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
