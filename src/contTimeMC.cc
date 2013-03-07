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
			    if (!(pi[i] > 0.0))
				qij = 0.0; // rows corresponding to stop codons are also set to 0
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
    double rowsum;
    double scale = 0.0; // factor by which matrix must be scaled to achieve one expected mutation per time unit
    for (int i=0; i<64; i++){
	rowsum = 0.0;
	for (int j=0; j<64; j++)
	    rowsum += gsl_matrix_get(Q, i, j);
	gsl_matrix_set(Q, i, i, -rowsum); // q_{i,i} = - sum_{j!=i} q_{i,j}
	scale += rowsum * pi[i];
    }
    if (scale != 0.0)
	scale = 1.0/scale;
    else
	scale = 1.0; // should not happen
    // normalize Q so that the exptected number of codon mutations per time unit is 1
    gsl_matrix_scale(Q, scale);
    return Q;
}

/*
 * perform a decompososition of the rate matrix as Q = U * diag(lambda) * U^{-1}
 */
int eigendecompose(gsl_matrix *Q,        // input
		    double *pi,          // input, must be same as in construction of Q
		    gsl_vector *&lambda, // output, memory for lambda, U and Uinv is allocated within the function
		    gsl_matrix *&U,      // output
		    gsl_matrix *&Uinv){  // output
    
    // Codon usages pi[j] that vanish are a problem because of division by zero.
    // Therefore, for the purpose of this decomposition, treat any such pi[j] as if 1 instead. This is 
    // not a problem as the j-th column and row of Q is a zero column anyways.
    
    // Compute B = diag(pi^{1/2}) * Q * diag(pi^{-1/2}) which is symmetrical, as Q is time-reversible
    gsl_matrix* B = gsl_matrix_calloc (64, 64);
    for (int i=0; i<64; i++){
	for (int j=0; j<64; j++){
	    if (pi[i] > 0.0 && pi[j] > 0.0){
		double bij = gsl_matrix_get(Q, i, j) * sqrt(pi[i]) / sqrt(pi[j]); // Q vanishes otherwise anyways
		gsl_matrix_set(B, i, j, bij);
	    }
	}
    }

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (64); // work space
    lambda = gsl_vector_calloc (64);
    U = gsl_matrix_calloc (64, 64);
    Uinv = gsl_matrix_calloc (64, 64);
    int status = gsl_eigen_symmv (B, lambda, U, w); // now: B = U * diag(lambda) * U^t, for orthogonal U, i.e. U^t = U^{-1}
    if (status)
	return status;
    gsl_eigen_symmv_free(w);

    // reverse transformation to Q = diag(pi^{-1/2}) * B * diag(pi^{1/2}) = [diag(pi^{-1/2}) * U] * diag(lambda) * [U^t * diag(pi^{1/2})]
    // reuse space of U
    // set Uinv <- U^t * diag(pi^{1/2})
    for (int j=0; j<64; j++){
	double factor = (pi[j]>0.0)? sqrt(pi[j]) : 1;
	for (int i=0; i<64; i++){
	    gsl_matrix_set(Uinv, i, j, factor * gsl_matrix_get(U, j, i)); // (this U is orthogonal)
	}
    }
    // and U <- diag(pi^{-1/2}) * U
    for (int i=0; i<64; i++){
	double factor = (pi[i] > 0.0)? 1/sqrt(pi[i]) : 1;
	for (int j=0; j<64; j++){
	    gsl_matrix_set(U, i, j, factor * gsl_matrix_get(U, i, j));
	}
    }
    // to test whether the decomposition is correct, look at whether Q - U * diag(lambda) * Uinv is close to the 0-matrix
    /*cout << "residue Q - U * diag(lambda) * Uinv" << endl;
    for (int i=0; i<64; i++){
	for (int j=0; j<64; j++){
	    double residue = gsl_matrix_get(Q, i, j);
	    for (int k=0; k<64; k++){
		residue -= gsl_matrix_get(U, i, k) * gsl_vector_get(lambda, k) * gsl_matrix_get(Uinv, k, j);
	    }
	    cout << setw(6) << fixed << setprecision(3) << residue << " ";
	}
	cout << endl;
    }
    */
    return 0;
}

/*
 * compute the matrix exponential P(t) = exp(Q*t),
 * where Q is given by U, diag(lambda) and U^{-1} as computed by eigendecompose
 */
gsl_matrix *expQt(double t, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv){
    gsl_matrix* P = gsl_matrix_alloc (64, 64);    
    double Pij;
    // P(t) = exp(Q*t) = U * diag(exp(lambda_1 * t), ..., exp(lambda_n * t)) * Uinv
    for (int i=0; i<64; i++){
	for (int j=0; j<64; j++){
	    Pij = 0;
	    for (int k=0; k<64; k++)
		Pij += gsl_matrix_get(U, i, k) * exp(gsl_vector_get(lambda, k)*t) * gsl_matrix_get(Uinv, k, j);
	    gsl_matrix_set(P, i, j, Pij);
	}
    }
    return P;
}

/*
 * print a 64x64 codon matrix with neat row and column labels by codon and amino acid
 */
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
	    if (abs(gsl_matrix_get(M,i,j)) > 1e-6)
		cout << setw(6) << fixed << setprecision(3) << gsl_matrix_get(M,i,j) << " ";
	    else 
		cout << setw(7) << " ";
	}
	cout << endl;
    }
}
