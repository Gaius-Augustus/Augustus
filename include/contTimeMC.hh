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

using namespace std;

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa); // transition/transversion ratio, usually >1
void printCodonMatrix(gsl_matrix *M);
#endif    // _CONTTIMEMC_HH
