/*****************************************************************************\
 * Filename : contTimeMC.hh
 * Authors  : Mario Stanke
 *
 * Description: continuous-time Markov chains for sequence evolution
 *              codon models
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 17.02.2013 | Mario Stanke          | creation of the class
\******************************************************************************/
#ifndef _CONTTIMEMC_HH
#define _CONTTIMEMC_HH

// project includes
#include "geneticcode.hh"

// standard C/C++ includes
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

using namespace std;

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa); // transition/transversion ratio, usually >1

/*
 * perform a decompososition of the rate matrix as Q = U * diag(lambda) * U^{-1}
 * returns 0 if sucessful
 */
int eigendecompose(gsl_matrix *Q,       // input
		   double *pi,          // input, must be same as in construction of Q
		   gsl_vector *&lambda, // output, memory for lambda, U and Uinv is allocated within the function
		   gsl_matrix *&U,      // output
		   gsl_matrix *&Uinv);  // output

/*
 * compute the matrix exponential P(t) = exp(Q*t),
 * where Q is given by U, diag(lambda) and U^{-1} as computed by eigendecompose
 */

gsl_matrix *expQt(double t, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv);



void printCodonMatrix(gsl_matrix *M);
#endif    // _CONTTIMEMC_HH
