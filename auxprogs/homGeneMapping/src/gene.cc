
/**********************************************************************
 * file:    gene.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
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
    {"CDS", "intron"};

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

long int Gene::getStart() const {
    
    long int geneStart = -1;
    if(strand == plusstrand)
	geneStart = (tlStart >= 0)? tlStart : features.front()->getStart();
    else
	geneStart = (tlEnd >=0) ? tlEnd-2 : features.front()->getStart();
    return geneStart;
}

long int Gene::getEnd() const {


    long int geneEnd = -1;
    if(strand == plusstrand)
	geneEnd = (tlEnd >= 0)? tlEnd : features.back()->getEnd();
    else
	geneEnd = (tlStart >=0) ? tlStart+2 : features.back()->getEnd();
    return geneEnd;
}
