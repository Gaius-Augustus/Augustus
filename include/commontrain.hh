/***********************************************************************
 * file:    commontrain.hh
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  used in training of all models
 * authors: Mario Stanke (mario@gobics.de)
 *
 **********************************************************************/

#ifndef _COMMONTRAIN_HH
#define _COMMONTRAIN_HH

// project includes
#include "types.hh"

// standard C/C++ includes
#include <vector>


/**
 * @memo A class for smoothing distributions
 *
 * @doc
 *
 * @author  Mario Stanke
 * @version 0.1
 */

// everything below this is considered 0, for performance
#define SMOOTH_EPSILON Double(1e-20)

class Smooth {
public:
    Smooth (Integer minwindowcount = 0, double slope_of_bandwidth = 0.1) {
	this->minwindowcount = minwindowcount;
	this->slope_of_bandwidth = slope_of_bandwidth;
    }
    ~Smooth() {}

    /*
     * smooth the vector counts of frequencies up to position n
     */
    void smoothCounts(const vector<Integer> &counts, vector<Double>& result, int resultSize=-1);
   
    /*
     * compute an optimal cutoff value d for the start of the geometric
     * distribution of introns
     */
    int geoCutOff(const vector<Integer> &lengths, vector<Double>& result);

private:
    // h: bandwidth, i difference
    // simple triangle kernel
    Double phi_triangle(Integer h, Integer i){
	if (i<0) {
	    i=-i;
	}
	if (i <= h) {
	    return Double((1 - (double) i/h)/h);
	} else 
	    return Double(0.0);
    }
    // h: bandwitdth, standard error=sigma
    // normal distribution kernel
    Double phi_normal(double stderror, Integer i){
	static const Double factor(0.39894228) ;    //  = 1/sqrt(2 pi)
	return factor / stderror * exp(- i / stderror * i / stderror / 2);
    } 

    Integer minwindowcount;      // number of events that need to be at least in each smoothing window
    double slope_of_bandwidth;  // bandwidth of kernel is proportional to index
};

void scaleDblVector(vector<Double>& v, Double sum);

#endif   // _COMMONTRAIN_HH
