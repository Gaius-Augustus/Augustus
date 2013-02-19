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
#include "geneticcode.hh"

// standard C/C++ includes
#include <iostream>
#include <iomanip>

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa){ // transition/transversion ratio, usually >1
    gsl_matrix* Q = gsl_matrix_calloc (64, 64); // initialize Q as the zero-matrix
    Seq2Int s2i(3);

    // initialize off-diagonal elements
    // codon i (row)    abc
    // codon j (col)    dbc
    // f                012
    int codoni[3], codonj[3];
    for (int a=0; a<4; a++){
	codoni[0] = a;
	for (int b=0; b<4; b++){
	    codoni[1] = b;
	    for (int c=0; c<4; c++){
		codoni[2] = c;
		int i = 16*codoni[0] + 4*codoni[1] + codoni[2]; // in {0,...63}
		// go through the 3x3 codons that differ from codon i by a single mutation
		for (int f=0;f<3; f++){ // position of mutation
		    for (int d=0; d<4; d++){ // d = substituted base
			if (d != codoni[f]){
			    codonj[0] = codoni[0];
			    codonj[1] = codoni[1];
			    codonj[2] = codoni[2];
			    codonj[f] = d;
			    int j = 16*codonj[0] + 4*codonj[1] + codonj[2]; // in {0,...63}
			    double qij = pi[j];
			    if (GeneticCode::is_purine(codonj[f]) == GeneticCode::is_purine(codoni[f]))
				qij *= kappa; // transition
			    if (GeneticCode::translate(i) != GeneticCode::translate(j))
				qij *= omega; // non-synonymous subst.
			    gsl_matrix_set (Q, i, j, qij);
			}
		    }
		}
	    }
	}
    }

    // initialize diagonal elements (each row must sum to 0)
    

    // normalize Q so that the exptected number of codon mutations per time unit is 1
    
    return Q;
}

void printCodonMatrix(gsl_matrix *M) {    // M is usually Q or P = exp(Qt)
    cout << "      ";
    for (int j=0; j<64; j++)
	cout << "  " << Seq2Int(3).inv(j) << "  "; // head line codon
    cout << endl << "      ";;
    for (int j=0; j<64; j++)
	cout << "   " << GeneticCode::translate(j) << "   "; // head line amino acid
    cout << endl;
    for (int i=0; i<64; i++){
	cout << Seq2Int(3).inv(i) << " " << GeneticCode::translate(i) << " ";
	for (int j=0; j<64; j++){
	    if (gsl_matrix_get(M,i,j) > 0.0)
		cout << setw(6) << fixed << setprecision(4) << gsl_matrix_get(M,i,j) << " ";
	    else 
		cout << setw(7) << " ";
	}
	cout << endl;
    }
}
