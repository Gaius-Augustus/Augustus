/*****************************************************************************\
 * Filename : contTimeMC.cc
 * licence  : Artistic Licence, see file LICENCE.TXT or 
 *            http://www.opensource.org/licenses/artistic-license.php
 * Authors  : Mario Stanke, mario.stanke@uni-greifswald.de
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 19.02.2013 | Mario Stanke          | creation of the class
\******************************************************************************/

// project includes
#include "contTimeMC.hh"

// standard C/C++ includes
#include <iostream>


/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa){ // transition/transversion ratio, usually >1
    gsl_matrix* Q = gsl_matrix_alloc (64, 64);
    
    return Q;
}
