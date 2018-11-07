/**********************************************************************
 * file:    gene.cc
 * license: Artistic License, see file LICENSE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 * date    |     author      |  changes 
 * --------|-----------------|------------------------------------------
 * 22.01.15| Stefanie Koenig | creation of the file
 **********************************************************************/

// project includes
#include "gene.hh"
#include "projectio.hh"

using namespace std;

string strandIdentifiers[NUM_STRAND_TYPES]=
    {"+", "-"};

string featureTypeIdentifiers[NUM_FEATURE_TYPES]=
    {"CDS", "intron","exon", "UTR", "start_codon", "stop_codon"};

string GeneFeature::getGeneID() const {
    return gene->getGeneID();
}
string GeneFeature::getTxID() const {
    return gene->getTxID();
}
int GeneFeature::getSeqID() const {
    return gene->getSeqID();
}

string GeneFeature::getSource() const {
    return gene->getSource();
}

void GeneFeature::setGene(Gene* gene){
    this->gene = gene;
}
Gene* GeneFeature::getGene(){
    return gene;
}

string GeneFeature::writeFrame() const {
    if(frame >= 0)
	return itoa(frame);
    else
	return ".";
}

bool GeneFeature::sameStrand(Strand other){
    if(this->strand == unknown || other == unknown)
	return true;
    else if(this->strand == other)
	return true;
    return false;
}

bool GeneFeature::sameFrame(int other){
    if(this->frame == -1 || other == -1)
	return true;
    else if(this->frame == other)
	return true;
    return false;
}

void GeneFeature::setEvidence(string e){
    if(!hasEvidence(e))
	extrinsic.push_back(e);
}
string GeneFeature::getEvidence() const {
    string e = "";
    for (int i = 0; i<extrinsic.size(); i++)
	(e.empty())? e=extrinsic[i] : e = e + ":" + extrinsic[i];
    return e;
}
bool GeneFeature::hasEvidence(string e) const {
    for(int i = 0; i < extrinsic.size(); i++){
	if(extrinsic[i] == e)
	    return true;
    }
    return false;
}

Strand getStrand(string token){
    Strand strand = unknown;
    if(token == "+")
	strand = plusstrand;
    if(token == "-")
	strand = minusstrand;
    return strand;
}

int getFrame(string token){
    int frame = -1;
    if(token == "0")
	frame = 0;
    if(token == "1")
	frame = 1;
    if(token == "2")
	frame = 2;  
    return frame;
}

FeatureType getType(std::string token){
    if(token == "CDS")
	return CDS;
    if(token == "intron")
	return intron;
    if(token == "exon")
	return exon;  
    if(token == "start_codon")
	return start;
    if(token == "stop_codon")
	return stop;
    if(token.find("UTR") !=  string::npos)
	return UTR;
    return unkown;
}

int Gene::numGFs(FeatureType t) const {
    int num=0;
    for(std::list<GeneFeature*>::const_iterator it=features.begin(); it!=features.end(); it++){
	if((*it)->isType(t))
	    num++;
    }
    return num;
}

GeneFeature* Gene::getFirstGF(FeatureType t){
    for(std::list<GeneFeature*>::iterator it=features.begin(); it!=features.end(); it++){
	if((*it)->isType(t))
	    return *it;
    }
    return NULL;
}

GeneFeature* Gene::getLastGF(FeatureType t){
    for(std::list<GeneFeature*>::reverse_iterator it=features.rbegin(); it!=features.rend(); it++){
	if((*it)->isType(t))
	    return *it;
    }
    return NULL;
}

void Gene::includeStopInCDS(){
    sortGFs();
    if(tlEnd >= 0){
	if(strand == plusstrand){
	    GeneFeature *last = getLastGF(CDS);
	    if(last){
		long int cdsEnd = last->getEnd();
		if(cdsEnd == tlEnd+2){ // stop already included
		}
		else if( cdsEnd == tlEnd-1){ // include stop
		    last->setLen(last->getLen()+3);
		}
		else if( cdsEnd < tlEnd-1){ // stop codon is a separate CDS exon
		    FeatureType type = CDS;
		    GeneFeature *exon = new GeneFeature(type, tlEnd,tlEnd+2,strand,0);
		    appendFeature(exon);
		    exon->setGene(this);
		}
	    }
	}
	else if(strand == minusstrand){
	    GeneFeature *first = getFirstGF(CDS);
	    if(first){
		long int cdsStart = first->getStart();
		if(cdsStart == tlEnd){ // stop already included

		}
		else if( cdsStart == tlEnd+3){ // include stop
		    first->setStart(cdsStart-3);
		    first->setLen(first->getLen()+3);
		}
		else if( cdsStart > tlEnd+3){ // stop codon is a separate CDS exon
		    FeatureType type = CDS;
		    GeneFeature *exon = new GeneFeature(type, tlEnd,tlEnd+2,strand,0);
		    appendFeature(exon);
		    exon->setGene(this);
		}
	    }
	}
    }
}
void Gene::insertExons(){

    sortGFs();
    GeneFeature* pred_exon = NULL;
    for(list<GeneFeature*>::iterator it=features.begin(); it!=features.end();it++){
        if (!(*it)->isCDS() && !(*it)->isUTR())
            continue;
        if(pred_exon && pred_exon->getEnd() + 1 == (*it)->getStart()){
            pred_exon->setLen((*it)->getEnd()-pred_exon->getStart()+1);
        }
        else{
            FeatureType type = exon;
            GeneFeature *exon = new GeneFeature(type, (*it)->getStart(),(*it)->getEnd(),strand);

            appendFeature(exon);
            exon->setGene(this);

            pred_exon = exon;
        }
    }
}

void Gene::insertMissingGFs(){

    int numIntr = numGFs(intron);
    int numEx = numGFs(exon);

    // insert exons if missing
    if(numEx == 0)
        insertExons();

    // delete UTR GeneFeatures - not needed anymore
    for(list<GeneFeature*>::iterator it=features.begin(); it!=features.end();){
	if ((*it)->isUTR()){
	    delete *it;  
	    it = features.erase(it);
	}
	else {
	    ++it;
	}
    }
    sortGFs();
    numEx = numGFs(exon);
    // insert introns if missing
    if(numIntr == 0 && numEx > 1){
        GeneFeature *pred_exon = NULL;
        for(list<GeneFeature*>::iterator it=features.begin(); it!=features.end();it++){
            if (!(*it)->isExon())
                continue;
            GeneFeature *exon = (*it);
            if(pred_exon){
                FeatureType type = intron;
                GeneFeature *intron = new GeneFeature(type, pred_exon->getEnd(),exon->getStart(),strand);
                appendFeature(intron);
                intron->setGene(this);
            }
            pred_exon = exon;
        }
    }
    // insert start_codon
    if(getTLstart() >= 0){
        FeatureType type = start;
        GeneFeature *start = new GeneFeature(type, getTLstart(), getTLstart()+2, strand);
        appendFeature(start);
        start->setGene(this);
    }
    // insert stop_codon
    if( getTLend() >= 0){
        FeatureType type = stop;
        GeneFeature *stop = new GeneFeature(type, getTLend(), getTLend()+2, strand);
        appendFeature(stop);
        stop->setGene(this);
    }
    sortGFs();
}


bool compareGFs(GeneFeature *a, GeneFeature *b){
    if( a->getStart() != b->getStart() )
	return ( a->getStart() < b->getStart() );
    if ( a->getLen() != b->getLen() )
	return ( a->getLen() >= b->getLen() );
    return a->getFeatureType() <= b->getFeatureType();
}
