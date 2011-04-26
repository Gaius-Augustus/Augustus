/***********************************************************************
 * file:    commontrain.cc
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  used in training of all models
 * authors: Mario Stanke (mario@gobics.de)
 *
 **********************************************************************/

#include "commontrain.hh"

#include <iostream>


void scaleDblVector(vector<Double>& v, Double sum) {
    Double kumSum(0.0), factor;
    for (int i=0; i<v.size(); i++) {
	kumSum += v[i];
    }
    if (kumSum != 0) {
	factor = sum / kumSum;
	for (int i=0; i<v.size(); i++) {
	    v[i] *= factor;
	}
    }
}

/*--- Smooth methods -------------------------------------------------*/

void Smooth::smoothCounts(const vector<Integer> &counts, vector<Double> &result, int resultSize){

    int n = counts.size();
    if (resultSize < 0)
	resultSize = n;
    if (int(result.size()) > resultSize)
	resultSize = result.size();
    result.assign(resultSize, 0.0);
   
    int i, j;
    int bandwidth;
    int cumcountl, cumcountr, numevents=0;
    Boolean negligible;
    
    for (i=0; i<n; i++)
	numevents += counts[i];

    /* for debugging only
    Double sum2 = 0;
    for (i=0; i<n; i++) {
	sum2+=result[i];
    }
    cout << "Summe vorher: " << sum2 << endl; 

    int sum1=0;  
    for (i=0; i<n; i++) {
	sum1 += counts[i];
    }
    cout << "Anzahl vorher: " << sum1 << endl; 
    */

    /*
     *  loop over the positions relevant for the final result vector
     */

    for (i=0; (i<n) && (i < resultSize + 4 * slope_of_bandwidth * resultSize); i++) {
	if (counts[i]>0) {
	    cumcountl= 0;
	    cumcountr= 0;
	    bandwidth = (int) (.01 + slope_of_bandwidth*pow((double)numevents, -.2) * i);
	    if (bandwidth < 1)
		bandwidth = 1; 
	    for (j=i-bandwidth+1; j <= i+bandwidth-1; j++) {
		if (j>= 0 && j<n) {
		    if (j <= i)
			cumcountl += (counts[j]) ? 1 : 0;
		    if (j>= i)
			cumcountr += (counts[j]) ? 1 : 0;		    
		}
	    }
	    // starting from minWindowsSize, enlarge bandwidth until in the 
	    // symmetric window with radius bandwidth are at least minwindowcount 
	    // positions with at least one event
	    while (cumcountl < minwindowcount && cumcountr < minwindowcount && bandwidth < n){
		bandwidth++;
		if (i+bandwidth-1<n) {
		    cumcountl += (counts[i+bandwidth-1]) ? 1 : 0;
		}
		if (i-bandwidth+1>=0) {
		    cumcountr += (counts[i-bandwidth+1]) ? 1 : 0;
		}
	    }
	    if (i<resultSize)
		result[i] += phi_normal(bandwidth, 0) * counts[i];
	    negligible = false;
	    j=1;
	    while (!negligible && (i-j>=0 || i+j<resultSize)){
		Double weight_j = phi_normal(bandwidth, j) * counts[i];
		if (i-j>=0 && i-j<resultSize ) {
		    result[i-j] += weight_j;
		}
		if (i+j<resultSize && i+j>=0) {
		    result[i+j] += weight_j;
		}
		negligible = (weight_j < SMOOTH_EPSILON);
		j++;
	    }
	}
    }
}


/*
 * Smooth::geoCutOff(const vector<Integer> &lencount,  vector<Double>& result)
 * lencount[i] is the weighed number of introns of length i
 * result holds the probabilities of the smoothed length distribution, 
 * memory is allocated in this function
 */
int Smooth::geoCutOff(const vector<Integer> &lencount, vector<Double>& result){
    int d, dtemp;
    int maxpos;
    Double maxvalue, PlenD, PlenDp1, bestRelDiff;
    vector<Double> wholeLenDist;

    smoothCounts(lencount, wholeLenDist);
    scaleDblVector(wholeLenDist, 1.0);

    // search the maximum of the length distribution
    maxvalue = -1;
    maxpos = -1;
  
    for (int i = 0; i < wholeLenDist.size(); i++) {
	if (wholeLenDist[i] > maxvalue) {
	    maxvalue = wholeLenDist[i];
	    maxpos = i;
	}
    }
    cout << "the most probable length is " << maxpos << endl;
    
    /*
     * search for the best place d to cut off
     */

    bestRelDiff = 1000.0;       // = "infinity"
    d = 0;
    Double relDiff;

    /*
     * start with dtemp as the positiion of the maximum of the length distribution
     * increase dtemp and check the size of the jump of the model-distribution
     * if we introduce the cutoff at dtemp
     */
    int numleD = 0, // number of intron with length <= dtemp
	numgD = 0,  // number of intron with length > dtemp
	sumgD = 0;  // sum of intron length > dtemp
    for (int i=0; i < maxpos; i++) {
	numleD += lencount[i];
    } 
    for (int i = maxpos + 1; i < lencount.size(); i++) {
	numgD += lencount[i];
	sumgD += lencount[i] * i;
    }

    for (dtemp = maxpos; dtemp < int(lencount.size()); dtemp++){
	PlenD = wholeLenDist[dtemp]; // we make a very small error here
	PlenDp1 = (double) numgD/(numleD + numgD) /((double) sumgD/numgD-dtemp);
	relDiff = (PlenDp1 > 0) ? PlenD/PlenDp1 : 1000; 

/*	cout << "d= " << dtemp << " ";
	cout << " P(len=D) = " << PlenD << "   P(len=D+1) = " << PlenDp1 <<
	" numleD = " << numleD << endl;*/

	if (relDiff < 1.0 && relDiff > 0) {
	    relDiff = Double(1.0) / relDiff;
	}
	if (relDiff == 0) {
	    relDiff = 1000;
	}
	if (relDiff < bestRelDiff) {
	    bestRelDiff = relDiff;
	    d = dtemp;	
	}
	// TODO: tradeoff between performance and accuracy
	if (relDiff < 1 + (double) dtemp/1000) { 	// relDiff < 1.2  
	    // good cutoff found, don't search any further
	    d = dtemp;
	    break;
	}
	// compute the new counters
	numleD += lencount[dtemp+1]; 
	numgD -= lencount[dtemp+1]; 
	sumgD -= lencount[dtemp+1] * (dtemp+1);
    } 

    result.resize(d+1);
    for (int i=0; i < result.size(); i++) {
	result[i] = wholeLenDist[i];
    }
    //shorte: result.assign(wholeLenDist.begin(), wholeLenDist.begin()+d+1);
    scaleDblVector(result, 1.0);
    return d;
}
