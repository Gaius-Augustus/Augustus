/*
 * codonevo_unittest.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Unit tests for CodonEvo (codoenevo.cc)
 */

#include "contTimeMC.hh"
#include "codonevo.hh"
#include "geneticcode.hh"
#include "gtest/gtest.h"
#include <stdlib.h>     /* rand */
#include <unistd.h>

int unlink(const char *path) throw ();


namespace {
    TEST(CodonEvoTest, CodonEvoRateReadWrite)
    // writing and reading codon rate matrix to/from file
    {
        // random equilibrium distribution pi
        double *pi = new double[64];
        double normsum = 0.0;
        for (int c=0; c<64; c++){
            pi[c] = 1 + rand() % 10;
            normsum += pi[c];
        }

        for (int c=0; c<64; c++)
            pi[c] /= normsum;
        
        CodonEvo codonevo;
        codonevo.setKappa(4.0);
        codonevo.setPi(pi);
        
        codonevo.setOmegas(3); // k=3
        codonevo.getRateMatrices();
        const char *ratesFname("rates-Q.txt");
        codonevo.writeRateMatrices(ratesFname);
        gsl_matrix *Q0, *Q2;
        // change some values
        Q0 = codonevo.getSubMatrixQ(0);
        Q2 = codonevo.getSubMatrixQ(2);
        gsl_matrix_set (Q0, 0, 0, 0.12345);
        gsl_matrix_set (Q2, 0, 1, 0.012345);
        // check that the values have indeed changed
        EXPECT_NEAR(gsl_matrix_get(Q0, 0, 0), .12345, 1e-6);
        EXPECT_NEAR(gsl_matrix_get(Q2, 0, 1), 0.012345, 1e-6);

        codonevo.readRateMatrices(ratesFname);

        EXPECT_NEAR(gsl_matrix_get(Q0, 0, 0), -1.11186, 1e-6);
        EXPECT_NEAR(gsl_matrix_get(Q2, 0, 1), 0.07978618, 1e-6);

        unlink(ratesFname); // delete temporary file
    }
    TEST(CodonEvoTest, CodonEvoRate) 
    {
        vector<double> branchset {0.3, 4.};
        CodonEvo codonevo; // implicitly uses uniform equilibrium
                           // distribution on sense codons
        codonevo.setKappa(4.0);
        codonevo.setBranchLengths(branchset);
        
        codonevo.setOmegas(20);
        // vector<double> *omegas = codonevo.getOmegas();
        // cout << "omegas: " << *omegas << endl;
        // omegas: [0.166667 , 0.25 , 0.333333 , 0.416667 , 0.5 ,
        // 0.583333 , 0.666667, 0.75 , 0.833333 , 0.916667 , 1 ,
        // 1.09091 ,  1.2 , 1.33333 , 1.5 , 1.71429 , 2 , 2.4 , 3 , 4]

        //codonevo.setPrior(0.5);
        codonevo.getRateMatrices();
        codonevo.computeLogPmatrices();

        gsl_matrix *P00 = codonevo.getSubMatrixP(0, 0.1);
        // printCodonMatrix(P00);
        //gsl_matrix_fprintf(stdout, P00, "%g")
        // below reference values were computed in Python with
        // Mario Stanke: molevol20/uebungen/blatt6/GY94-codonevo.ipynb
        // P('AAA' -> 'AAC'; t=0.3, omega = 1/6)
        EXPECT_NEAR(gsl_matrix_get(P00, 0, 1), 0.0062754453730753945, 1e-9);
        // P('AAA' -> 'AAG'; t=0.3, omega = 1/6)
        EXPECT_NEAR(gsl_matrix_get(P00, 0, 2), 0.12873728943730725, 1e-9);
    }
}
