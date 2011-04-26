#include "pp_fastBlockSearcher.hh"

#include <sstream>

int    PP::FsHitType::maxIntronLen = 100000;
double PP::FsHitType::intronMalus  = 1e-5;

int    PP::FsSeedCollection::expSeedCount = 10000;
double PP::FsSeedCollection::maxCoverage  = 0.8;

void PP::FsHitCollection::newHit(FsHitType* ht) {
    allHits.push_back(ht);
    vector<HitQueue>& allQueues = pendingHits[ht->reverse];
    int b = ht->reverse ? ht->blockNo+1 : ht->blockNo-1;
    if (b>=0 && b<size) 
	ht->linkTo(allQueues[b]);
    if (ht->reverse) 
	for (b=ht->blockNo; b>0; b--)
	    ht->pushOn(allQueues[b]);
    else
	for (b=ht->blockNo; b<size-1; b++)
	    ht->pushOn(allQueues[b]);

    // Put this entry in to the final result container
    if (finalResult.empty()) {
	finalResult.push_back(ht);
	return;
    }
    int i = finalResult.size()-1;
    while (i>=0 && finalResult[i]->start() > ht->start()) 
	i--;
    if (i>=0 && finalResult[i]->head == ht->head) {
	if (finalResult[i]->pathScore < ht->pathScore) 
	    finalResult[i] = ht;
    } else if (i == finalResult.size()-1) {
	finalResult.push_back(ht);
    } else  {
	int j = finalResult.size()-1;
	finalResult.push_back(finalResult.back());
	while (--j>i)
	    finalResult[j+1] = finalResult[j];
	finalResult[i+1] = ht;
    }
}

void PP::FsHitCollection::storeBestResults(multimap<double, FsHitType>& result, int mincount, double threshold) {
    int maxcount = finalResult.size();
    result.clear();
	
    if (mincount > maxcount)
	mincount = maxcount;
    int i;
    for (i=0; i<mincount; i++) {
	FsHitType* ht = finalResult[i];
	result.insert(make_pair(ht->pathScore, *ht));
    }
    if (mincount) for(; i<maxcount; i++) {
	FsHitType* ht = finalResult[i];
	double sc = ht->pathScore;
	double listSc = result.begin()->second.pathScore;
	if (listSc > threshold)
	    break;
	if (listSc > sc)
	    continue;
	result.erase(result.begin());
	result.insert(make_pair(sc, *ht));
    }
    for (; i<maxcount; i++) {
	FsHitType* ht = finalResult[i];
	if (ht->pathScore > threshold)
	    result.insert(make_pair(ht->pathScore, *ht));
    }
}
