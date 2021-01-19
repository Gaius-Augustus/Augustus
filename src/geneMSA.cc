/*
 * geneMSA.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Generation of exon candidates
 */
 
#include "geneticcode.hh"
#include "genomicMSA.hh"
#include "geneMSA.hh"
#include "intronmodel.hh"
#include "namgene.hh"
#include "orthograph.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <unordered_set>

using namespace std;

PhyloTree *GeneMSA::tree = NULL;
CodonEvo *GeneMSA::codonevo = NULL;
int GeneMSA::padding = 1000; // added to seqRange after last alignment block
int GeneMSA::orthoExonID = 1;
int GeneMSA::geneRangeID = 1;
vector<int> GeneMSA::exonCandID;
vector<ofstream*> GeneMSA::exonCands_outfiles;
vector<ofstream*> GeneMSA::orthoExons_outfiles;
vector<ofstream*> GeneMSA::geneRanges_outfiles_gff;
vector<ofstream*> GeneMSA::geneRanges_outfiles_bed;
vector<ofstream*> GeneMSA::omega_outfiles;
unordered_map< bit_vector, PhyloTree*, boost::hash<bit_vector>> GeneMSA::topologies;
map<vector<string>, pair<vector<double>, int> > GeneMSA::computedCumValues;

		

/*
 * constructor of GeneMSA
 */
GeneMSA::GeneMSA(RandSeqAccess *rsa, Alignment *a) {
    int maxDNAPieceSize = Properties::getIntProperty( "maxDNAPieceSize" );
    this->rsa = rsa;
    alignment = a;
    if (!alignment)
	throw ProjectError("Internal error in GeneMSA: Alignment missing.");    
    if (alignment->numRows() != rsa->getNumSpecies())
	throw ProjectError("Error in GeneMSA: Number of species in alignment (" + itoa(alignment->numRows())
			   + ") is not matching the one in the tree (" + itoa(rsa->getNumSpecies()) + ").");
    ltree = NULL; // locus/gene tree, may be different from species tree
    exoncands.resize(alignment->numRows(), NULL);

    /** construct the gene ranges
     * now: simple copy. TODO: extend region when apparently part of the alignment is missing
     * human   ***********---*******************
     * mouse   *******-------******************-
     * chicken ***********----------------------
     *                        ^
     *                        | extend range here
     */
    starts.resize(alignment->numRows(), -1);
    ends.resize(alignment->numRows(), -1);
    offsets.resize(alignment->numRows(), 0);
    for (size_t s=0; s < alignment->numRows(); s++){
	if (alignment->rows[s] && alignment->rows[s]->frags.empty()){
	    delete alignment->rows[s];
	    alignment->rows[s] = NULL;
	}
	if (alignment->rows[s]){
	    int chrLen = rsa->getChrLen(s, alignment->rows[s]->seqID);
	    switch (getStrand(s))
		{
		case plusstrand:            
            starts[s] = alignment->rows[s]->chrStart() - padding;
		    ends[s] = alignment->rows[s]->chrEnd() + padding;
		    // cout << "plus " << chrLen << " " << alignment->rows[s]->seqID << " " << alignment->rows[s]->chrStart() << " " << padding << " " << alignment->rows[s]->chrEnd() << " startEnd " << starts[s] << " " << ends[s] << endl;
		    break;
		case minusstrand:
		    // cout << "minus " << alignment->rows[s]->seqID << " " << chrLen << " " << alignment->rows[s]->chrStart() << " " << padding << " " << alignment->rows[s]->chrEnd() << endl;
		    starts[s] = chrLen - alignment->rows[s]->chrEnd() - padding;
		    ends[s] = chrLen - alignment->rows[s]->chrStart() + padding;
		    break;
		default:
		    throw ProjectError("GeneMSA: Unknown strand in alignment.");
		}
	    // ensure that gene predictions on each region will be done IN ONE PIECE
	    int geneRangeLen = ends[s]-starts[s]+1;
	    if (geneRangeLen > maxDNAPieceSize + 2*padding){
		cerr << "Warning: length of gene range(" << geneRangeLen << ") is species " << rsa->getSname(s) <<  " exceeds maxDNAPieceSize(" 
		     << maxDNAPieceSize << "). Will shorten sequence at both ends to achieve a length of " << maxDNAPieceSize + 2*padding << endl;
		int tooMuch = geneRangeLen - (maxDNAPieceSize + 2*padding);
		starts[s] += (tooMuch + 1)/2;
		ends[s] -= (tooMuch + 1)/2;
		// delete all fragments that do not overlap with the gene range  
                AlignmentRow *row = alignment->rows[s];
                vector<fragment>::iterator first = row->frags.begin();
                while (first != row->frags.end() && first->chrPos + first->len - 1 < starts[s]){
                    ++first;
                }
                row->frags.erase(row->frags.begin(),first);
                vector<fragment>::iterator last = row->frags.begin();
                while (last != row->frags.end() && last->chrPos < ends[s]){
                    ++last;
                }
                row->frags.erase(last,row->frags.end());
                // shorten fragments that are not fully included in the gene range    
                if(!row->frags.empty()){
                    if(row->chrStart() < starts[s]){
                        int diff= starts[s] - row->chrStart();
                        row->frags[0].chrPos += diff;
                        row->frags[0].aliPos += diff;
                        row->frags[0].len -= diff;
                    }
                    if(row->chrEnd() > ends[s]){
                        int diff= row->chrEnd() - ends[s];
                        row->frags[row->frags.size()-1].len -= diff;
		    }
		}
		//TODO: update cumFragLen (probably no longer needed)
	    }
	    if (starts[s] < 0)
		starts[s] = 0;
	    if (ends[s] < 0)
		ends[s] = 0;
	    if (ends[s] >= chrLen)
		ends[s] = chrLen-1;
	    if (starts[s] >= chrLen)
		starts[s] = chrLen-1;
	    offsets[s] = (getStrand(s) == plusstrand)? starts[s] : rsa->getChrLen(s, getSeqID(s)) - 1 - ends[s];
	    /*cout<<  "offset of species "<<s<<" "<<offsets[s]<<"\tstand=";
	    if(getStrand(s) == plusstrand)
	      cout<<"plusstrand";
	    else
	      cout<<"minusstrand";
	    cout<<endl;
	    */
	}
    }
}

/*
 * destructor of GeneMSA
 */
GeneMSA::~GeneMSA(){
    if (alignment)
	delete alignment;
    for (int i=0; i<exoncands.size(); i++) {
	if (exoncands[i]!=NULL) {
	    for(list<ExonCandidate*>::iterator it = exoncands[i]->begin(); it != exoncands[i]->end(); it++){
		delete *it;
	    }
	    exoncands.at(i)->clear();
	    delete exoncands[i];
	}
    }
}

string GeneMSA::getSeqID(int speciesIdx) {
    if (alignment->rows[speciesIdx])
	return alignment->rows[speciesIdx]->seqID;
    else
	return "";
}

Strand GeneMSA::getStrand(int speciesIdx){
    if (alignment->rows[speciesIdx])
	return alignment->rows[speciesIdx]->strand;
    else 
	return STRAND_UNKNOWN;
}

// adds the keys to the map function
map<string,ExonCandidate*>* GeneMSA::getECHash(list<ExonCandidate*> *ec) {
    map<string, ExonCandidate*> *hashCandidates = new  map<string, ExonCandidate*>;
    if (ec->empty())
	return NULL;
	    
    for (list<ExonCandidate*>::iterator lit=ec->begin(); lit!=ec->end(); lit++)
	(*hashCandidates)[(*lit)->key()] = *lit;

    return hashCandidates;
}

// computes and sets the exon candidates for species s
// and inserts them into the hash of ECs if they do not exist already
void GeneMSA::createExonCands(int s, const char *dna, map<int_fast64_t, ExonCandidate*> &ecs, map<int_fast64_t, ExonCandidate*> &addECs){
    double assmotifqthresh = 0.15;
    double assqthresh = 0.3;
    double dssqthresh = 0.7;
    bool onlySampledECs = false;
    int minEClen = 1;
    Properties::assignProperty("/CompPred/onlySampledECs", onlySampledECs);
    if(!onlySampledECs){
	Properties::assignProperty("/CompPred/assmotifqthresh", assmotifqthresh);
	Properties::assignProperty("/CompPred/assqthresh", assqthresh);
	Properties::assignProperty("/CompPred/dssqthresh", dssqthresh);
	// TODO Properties::assignProperty("/CompPred/minExonCandLen", minEClen);
    
	findExonCands(ecs, addECs, dna, minEClen, assmotifqthresh, assqthresh, dssqthresh); 
    }
}

void GeneMSA::setExonCands(vector<map<int_fast64_t, ExonCandidate*> > &ecs){
    for (int s = 0; s < ecs.size(); s++) {
        if (!ecs[s].empty()){
            list<ExonCandidate*> *candidates = new list<ExonCandidate*>;
            for(map<int_fast64_t, ExonCandidate*>::iterator ecit=ecs[s].begin(); ecit!=ecs[s].end(); ecit++){
                candidates->push_back(ecit->second);
            }
            exoncands[s] = candidates;
            ecs[s].clear(); // not needed anymore
            cout << "Found " << exoncands[s]->size() << " ECs on species " << rsa->getSname(s) << endl; 
        }
    }
}

/**
 * createOrthoExons
 */
void GeneMSA::createOrthoExons(list<OrthoExon> &orthoExonsList, map<int_fast64_t, list<pair<int,ExonCandidate*> > > &alignedECs, Evo *evo, float consThres, int minAvLen) {

    //cout << "Creating ortho exons for alignment" << endl << *alignment << endl;
    int k = alignment->rows.size();
    int m = alignment->numFilledRows();; // the number of nonempty rows

    // an ortho exon candidate must have an EC in at least this many species (any subset allowed):
    int minEC = (consThres * m > 2.0)? m * consThres + 0.9999 : 2;
    if(minEC != 2){
	cerr << "Warning: minEC=2, required by current version." << endl;
	minEC = 2;
    }
    cout << "OEs in this gene range must have at least " << minEC << " ECs" << endl;

    map<int_fast64_t, list<pair<int,ExonCandidate*> > >::iterator aec;
    unordered_map<bit_vector,PhyloTree*, boost::hash<bit_vector>>::iterator topit;
    /*
     * Create one ortho exon candidate for each key to which at least minEC exon candidates mapped 
     */
    int numOE = 0;
    for (aec = alignedECs.begin(); aec != alignedECs.end(); ++aec){

	// first remove "absent" ECs, e.g. all tuples (speciesIdx,NULL) from the list        
        bit_vector absent(k,0); // bit 1 for EC absent, but genome aligned
        list<pair<int,ExonCandidate*> >::iterator it = aec->second.begin();
	while (it != aec->second.end()){
	    if(!it->second){
		if(evo->getNumStates() > 2){ // additional state "EC is absent", only if ExonEvo model is initialized with 3 states
		    absent[it->first]=1;
		}
                it = aec->second.erase(it);
            }
            else{
                it++;
            }
        }
	
	// float avLen = 0.0; // not needed anymore
	OrthoExon oe(aec->first, numSpecies());
	bit_vector present(k,0);      // bit 1 for ECs that are present
	/*
	 * leaves: bit 1 if the corresponding leaf node is present in the tree
	 * ExonEvo(2) -> only leaves for "present" ECs
	 * ExonEvo(3) -> leaves for "present" and "absent" ECs
	 * ExonEvo(4) -> all leaves
	 */
	bit_vector leaves = absent; 
	oe.orthoex.resize(k, NULL);
	for (it = aec->second.begin(); it != aec->second.end(); ++it){
	    int s = it->first;
	    ExonCandidate *ec = it->second;
	    // cout << rsa->getSname(s) << "\t" << ec->getStart() + offsets[s] << ".." << ec->getEnd() + offsets[s] << "\t" << *ec << endl;
	    if (oe.orthoex[s])
		throw ProjectError("createOrthoExons: Have two exon candidates from the same species " 
				   + rsa->getSname(s) + " with the same key " + itoa(aec->first));
	    present[s]=1;
	    leaves[s]=1;
	    oe.orthoex[s] = ec;
	    // avLen += ec->len();
	}
	// avLen /= aec->second.size(); // compute average length of exon candidates in oe
	oe.setPresent(present);
	oe.setAbsent(absent);
	// link OE to tree topology
	PhyloTree *t = NULL;
	if(evo->getNumStates() > 3) // additional state "unaligned", only if ExonEvo model is initialized with 4 states
	    leaves = bit_vector(k,1); // original tree, no pruning necessary
	topit = topologies.find(leaves);
	if(topit == topologies.end()){ // insert new topology
	    t = new PhyloTree(*tree);
	    t->prune(leaves,evo); // remove all leaf nodes of "unaligned" species
	    topologies.insert(pair<bit_vector,PhyloTree*>(leaves,t));
	}
	else{
	    t = topit->second;
	}
	if (aec->second.size() >= minEC){ // if the number of exons >=2, create an OrthoExon
	    oe.ID = orthoExonID;
	    oe.setBV(present);
	    oe.setTree(t);	
	    oe.setContainment(0);
	    orthoExonID++;
	    numOE++;
	    orthoExonsList.push_back(oe);
	}
	else { // if the number of exons ==1 , no OrthoExon is needed, the phylogenetic score is constant and can be pre-computed
	    if(t->numSpecies() > 1){ // otherwise the tree has only a single node
		int s = aec->second.begin()->first;
		ExonCandidate *ec = aec->second.begin()->second;
		vector<int> fix_0 = oe.labels; // EC is not predicted
		vector<int> fix_1 = oe.labels; // EC is predicted
		fix_1[s]=1;
		// phylogenetic score, difference bewteen the scores of "EC is predicted" and "EC is not predicted"
		double score = t->MAP(fix_1, oe.weights, evo, true) - t->MAP(fix_0, oe.weights, evo, true);
		ec->setScore(score);
	    }
	}
    }

    alignedECs.clear(); // not needed anymore
    /*
     * Determine 'containment' for each OE: The average number of extra bases (per species) of the largest OE (y) that includes this OE (x) in frame.
     * Idea: A large containment is a sign that OE is actually not true.
     * species
     *  1       xxxxxxxxxxxxxxxxxxxxxxxxx           containment = 5
     *  1     YYyyyyyyyyyyyyyyyyyyyyyyyyyYY
     *  2       xxxxxxxxxxxxxxxxxxxxxxxxx
     *  2     YYyyyyyyyyyyyyyyyyyyyyyyyyyYYYY
     *  3   yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy     
     *  Containing OE y must have ECs for all species where x has ECs, and may have ECs for additional species.
     */
    list<OrthoExon>::iterator oeYit, oeXit, ovlpBeginIt;
    ovlpBeginIt = orthoExonsList.begin();
    for (oeYit = orthoExonsList.begin(); oeYit != orthoExonsList.end();	++oeYit){
	int YaliEnd = oeYit->getAliEnd();
	int YaliStart = oeYit->getAliStart();
	while (ovlpBeginIt != orthoExonsList.end() && ovlpBeginIt->getAliStart() < YaliStart)
	    ++ovlpBeginIt;
	// X ranges from the first OE with the same start as Y to the last OE with the start before Y's end
	oeXit = ovlpBeginIt;
	while (oeXit != orthoExonsList.end() && oeXit->getAliStart() <= YaliEnd){
	    if (oeXit != oeYit && // same OE would yield a containment of 0 anyways
		oeXit->getAliEnd() <= YaliEnd // x is contained in y in the alignment
		&& isOnFStrand(oeXit->getStateType()) == isOnFStrand(oeYit->getStateType())){ // same strand
		size_t containment = 0;
		bool contained = true;
		for (size_t s=0; s < numSpecies() && contained; s++) {
		    if (oeXit->orthoex[s]) {
			if (oeYit->orthoex[s]){	
			    int leftOverhang = oeXit->orthoex[s]->getStart() - oeYit->orthoex[s]->getStart();
			    int rightOverhang = oeYit->orthoex[s]->getEnd() - oeXit->orthoex[s]->getEnd();
			    bool frameCompatible = ((oeYit->orthoex[s]->getFirstCodingBase() - oeXit->orthoex[s]->getFirstCodingBase()) % 3 == 0);
			    if (leftOverhang >= 0 && rightOverhang >= 0 && frameCompatible)
				containment += leftOverhang + rightOverhang;
			    else 
				contained = false;	
			} else {
			    contained = false;
			}
		    }
		}
		if (contained) {
		    // x contains y in above sense
		    containment /= oeXit->numExons(); // average over species present in x
		    if (oeXit->getContainment() < containment)
			oeXit->setContainment(containment);
		}
	    }
	    oeXit++;
	}
    }
    cout << "Found " << numOE << " ortho exons" << endl;
}


void GeneMSA::openOutputFiles(string outdir){
    
    if (Constant::exoncands) // output of exon candidates into a gff file requested
	exonCands_outfiles.resize(tree->numSpecies());
    if(Constant::printOEs)
	orthoExons_outfiles.resize(tree->numSpecies());
    geneRanges_outfiles_gff.resize(tree->numSpecies());
    geneRanges_outfiles_bed.resize(tree->numSpecies());
    omega_outfiles.resize(tree->numSpecies());
    vector<string> species;
    tree->getSpeciesNames(species);
    for (int i=0; i<tree->numSpecies(); i++) {
        string file_exoncand = outdir + "exonCands." + species[i] + ".gff3";
        ofstream *os_ec = NULL;
	if (Constant::exoncands){
	    os_ec = new ofstream(file_exoncand.c_str());
	    if (os_ec) {
		exonCands_outfiles[i] = os_ec;
		(*os_ec) << PREAMBLE << endl;
		(*os_ec) << "#\n#-----  exon candidates  -----" << endl << "#" << endl;
	    }
        }
	string file_geneRanges;
	ofstream *os_gr;
	if (Constant::printGeneRangesBED){
	    file_geneRanges = outdir + "geneRanges." + species[i] + ".bed";
	    os_gr = new ofstream(file_geneRanges.c_str());
	    if (os_gr && (os_gr->rdstate() & std::ofstream::failbit) == 0) {
		geneRanges_outfiles_bed[i] = os_gr;
		(*os_gr) << "#\n#-----  gene ranges  -----" << endl << "#" << endl;
	    } else {
		cerr << "Error writing " << file_geneRanges << endl;
	    }
	}
	// same in gff format
	if (Constant::printGeneRangesGFF){
	    file_geneRanges = outdir + "geneRanges." + species[i] + ".gff";
	    os_gr = new ofstream(file_geneRanges.c_str());
	    if (os_gr && (os_gr->rdstate() & std::ofstream::failbit) == 0) {
		geneRanges_outfiles_gff[i] = os_gr;
		(*os_gr) << "#\n#-----  gene ranges  -----" << endl << "#" << endl;
	    } else {
		cerr << "Error writing " << file_geneRanges << endl;
	    }
	}
	
        string file_orthoexon = outdir + "orthoExons." + species[i] + ".gff3";
	if(Constant::printOEs){
	    ofstream *os_oe = new ofstream(file_orthoexon.c_str());
	    if (os_oe) {
		orthoExons_outfiles[i]=os_oe;
		(*os_oe) << PREAMBLE << endl;
		(*os_oe) << "#\n#----- ortholog exons  -----" << endl << "#" << endl;
	    }
	}
        /*string file_omega = outdir + "omegaExons." + species[i] + ".gff3";
        ofstream *os_omega = new ofstream(file_omega.c_str());
        if (os_omega) {
            omega_outfiles[i]=os_omega;
            (*os_omega) << PREAMBLE << endl;
            (*os_omega) << "#\n#----- exons with a dN/dS ratio smaller than one -----" << endl << "#" << endl;
	}*/
    }
}

void GeneMSA::printStats(){
    for(int s=0; s<numSpecies(); s++){
      cout << rsa->getSname(s) << "." << getSeqID(s) << "\t";
	if (alignment->rows[s])
	    cout << alignment->rows[s]->strand << "\t" << getStart(s) << "\t" << getEnd(s);
	cout << endl;
    }
}

// appends the gene range of each species to the file 'geneRanges.speciesnames.bed'
void GeneMSA::printGeneRanges() {
    for (size_t s=0; s < numSpecies(); s++) {
	if (getStart(s) >= 0) {
	    // output in .bed format
	    if (Constant::printGeneRangesBED) {
		ofstream &fstrm_bed = *geneRanges_outfiles_bed[s]; // write to 'geneRanges.speciesname[s].bed'
		fstrm_bed << getSeqID(s) << "\t" << getStart(s) /* + 1 GM removed this, since BED is 0-based*/ << "\t" << getEnd(s) + 1 << "\t" << alignment->getSignature()
			  << "\t0\t" << getStrand(s) << endl;
	    }
	    // GFF output
	    if (Constant::printGeneRangesGFF) {
		ofstream &fstrm_gff = *geneRanges_outfiles_gff[s]; // write to 'geneRanges.speciesname[s].gff'    
		fstrm_gff << getSeqID(s) << "\tGeneRange\t" << "exon\t" << getStart(s) + 1 << "\t" << getEnd(s) + 1 << "\t0\t"
			  << getStrand(s) << "\t" << ".\t" << "Name=" << geneRangeID;
		fstrm_gff << ";Note=" << alignment->getSignature() << endl;
	    }
	}
    }
    geneRangeID++;
}

// writes the exon candidates of all species into the file 'exonCands.species.gff3'
void GeneMSA::printExonCands() {
    exonCandID.resize(numSpecies(), 1);

    if (!exoncands.empty()) {
        for (int s=0; s < numSpecies(); s++) {
	    ofstream &fstrm = *exonCands_outfiles[s]; // write to 'exonCands.speciesname[i].gff3'
	    list<ExonCandidate*>* sec = exoncands[s];
            if (sec) {
                fstrm << "# sequence:\t" << rsa->getSname(s) << "\t" << getStart(s) + 1 << "-" 
		      << getEnd(s) + 1 << "  " << getEnd(s) - getStart(s) << "bp" << endl;
                for (list<ExonCandidate*>::iterator ecit = sec->begin(); ecit != sec->end(); ++ecit) {
                    fstrm << getSeqID(s)<< "\tEC\t" << "exon\t";
                    if (getStrand(s) == plusstrand) {
                        fstrm << (*ecit)->begin + offsets[s] + 1 << "\t" << (*ecit)->end + offsets[s] + 1 
			      << "\t" << (*ecit)->score << "\t";
                    } else {
			int chrLen = rsa->getChrLen(s, getSeqID(s));
                        fstrm << chrLen - ((*ecit)->end + offsets[s]) << "\t"
			      << chrLen - ((*ecit)->begin+ offsets[s]) << "\t"
			      << (*ecit)->score << "\t";
                    }
		    // the gff strand of the exon is the "strand product" of the alignment strand and exon type strand
		    // e.g. "-" x "-" = "+" << '+' << "\t";
		    fstrm << ((isPlusExon((*ecit)->type) == (getStrand(s) == plusstrand))? '+' : '-');
                    fstrm << "\t" << (*ecit)->gff3Frame() << "\t" << "ID=" << exonCandID[s] << ";"
			  << "Name=" <<stateExonTypeIdentifiers[(*ecit)->type];
		    fstrm << ";upSig=" << (*ecit)->getUpScore() << ";downSig=" << (*ecit)->getDownScore();
		    fstrm << endl;
		    // TODO: adjust type on reverse alignment setrand (compate comment in printSingleOrthoExon)

                    exonCandID[s]++;
		}
            } else {
                fstrm << "#  no exon candidates found " << endl;
            }
        }
    } else {
        cout << "#  no exon candidates found at all" << endl;
    }
}

// writes all ortholog exons of all species in the files 'orthoExons.species.gff3'
// orthoexons are sorted by alignment start coordinate
void GeneMSA::printOrthoExons(list<OrthoExon> &orthoExonsList) {
    if (orthoExonsList.empty())
	return;
    for (list<OrthoExon>::iterator oeit = orthoExonsList.begin(); oeit != orthoExonsList.end(); ++oeit){
        printSingleOrthoExon(*oeit, true);
    }
}


// writes the ortholog exons on one OrthoExon into the files 'orthoExons.species.gff3'
// files: write each one to a file for its species, if false to stdout
void GeneMSA::printSingleOrthoExon(OrthoExon &oe, bool files) {
    bool GBrowseStyle = false; // for viewing in GBrowse use this style
    streambuf *stdout = cout.rdbuf();
    int beginStopOffset, endStopOffset; // to account for excluded stop codons
    for (int s=0; s < numSpecies(); s++) {
	ExonCandidate *ec = oe.orthoex.at(s);
	if (files)
	    cout.rdbuf(orthoExons_outfiles[s]->rdbuf()); // write to 'orthoExons.speciesname[s].gff3'
        if (ec != NULL) {
            // the gff strand of the exon is the "strand product" of the alignment strand and exon type strand
            // e.g. "-" x "-" = "+"
            string gffStrand = ((isPlusExon(ec->type) == (getStrand(s) == plusstrand))? "+" : "-");
            beginStopOffset = endStopOffset = 0;
            if (hasStopCodon(ec->type) && Gene::stopCodonExcludedFromCDS){
                if (gffStrand == "+"){
                    endStopOffset = -3;
                } else {
                    beginStopOffset = 3;
                }
            }
	    cout << getSeqID(s) << "\tOE1\t" << "exon" << "\t";
            if (getStrand(s) == plusstrand){ // strand of alignment
		cout << ec->begin + offsets[s]+1 + beginStopOffset << "\t" << ec->end + offsets[s]+1 + endStopOffset;
            } else {
		int chrLen = rsa->getChrLen(s, getSeqID(s));
                cout << chrLen - (ec->end + offsets[s]) + beginStopOffset << "\t" << chrLen - (ec->begin + offsets[s]) + endStopOffset;
            }
	    cout << "\t" << ec->score << "\t" << gffStrand << "\t" << ec->gff3Frame() << "\t"
                 << "ID=" << oe.ID << ";Name=" << oe.ID << ";Note=" 
		 << stateExonTypeIdentifiers[ec->type]; // TODO: reverse type if alignment strand is "-"
	    if (GBrowseStyle)
		cout << "|" << oe.numExons();
	    else 
		cout << ";n=" << oe.numExons();
	    if (oe.getOmega() >= 0.0){
		if (GBrowseStyle)
		    cout << "|" << oe.getOmega();
		else 
		    cout << ";MLomega=" << oe.getOmega();
	    }
	    if (oe.getEomega() >= 0.0){
		if (GBrowseStyle)
		    cout << "|" << oe.getEomega();
		else
		    cout << ";Eomega=" << oe.getEomega();
	    }
	    if (oe.getVarOmega() >= 0.0){
		if (GBrowseStyle)
		    cout << "|" << oe.getVarOmega();
		else
		    cout << ";VarOmega=" << oe.getVarOmega();
	    }
	    if (oe.getLeftExtOmega() >= 0.0){
	      if (GBrowseStyle)
		cout << "|" << oe.getLeftExtOmega();
	      else
		cout << ";leftBoundaryExtOmega=" << oe.getLeftExtOmega();
            }
	    if (oe.getRightExtOmega() >= 0.0){
              if (GBrowseStyle)
                cout << "|" << oe.getRightExtOmega();
              else
                cout << ";rightBoundaryExtOmega=" << oe.getRightExtOmega();
            }
	    if (oe.getLeftIntOmega() >= 0.0){
              if (GBrowseStyle)
                cout << "|" << oe.getLeftIntOmega();
              else
                cout << ";leftBoundaryIntOmega=" << oe.getLeftIntOmega();
            }
            if (oe.getRightIntOmega() >= 0.0){
              if (GBrowseStyle)
                cout << "|" << oe.getRightIntOmega();
              else
                cout << ";rightBoundaryIntOmega=" << oe.getRightIntOmega();
            }
	    if (oe.getSubst() >= 0){ // number of substitutions
		if (GBrowseStyle)
		    cout << "|" << oe.getSubst();
		else
		    cout << ";subst=" << oe.getSubst();
	    }
	    if (oe.getConsScore() >= 0.0){ // conservation score
		if (GBrowseStyle)
		    cout << "|" << oe.getConsScore();
		else
		    cout << ";cons=" << oe.getConsScore();
	    }
	    if (oe.getLeftConsScore() >= 0.0){ // conservation score of left boundary feature
		if (GBrowseStyle)
		    cout << "|" << oe.getLeftConsScore();
		else
		    cout << ";LeftCons=" << oe.getLeftConsScore();
	    }
	    if (oe.getRightConsScore() >= 0.0){ // conservation score of right boundary feature
		if (GBrowseStyle)
		    cout << "|" << oe.getRightConsScore();
		else
		    cout << ";rightCons=" << oe.getRightConsScore();
	    }
	    if (oe.getDiversity() >= 0.0){ // diversity
		if (GBrowseStyle)
		    cout << "|" << oe.getDiversity();
		else
		    cout << ";div=" << oe.getDiversity();
	    }
	    if (GBrowseStyle)
		cout << "|" << oe.getContainment();
	    else
		cout << ";containment=" << oe.getContainment();
	    if (GBrowseStyle)
                cout << "|" << oe.getPhyleticPattern();
	    else
                cout << ";phyloPat=" << oe.getPhyleticPattern();
	    // output logistic regression score for analysis
	    cout << ";oescore=" << oe.getLogRegScore();
	    // output postprobs of all exons in OE
	    cout << ";ec_postProbs=";
	    for(int i=0; i < numSpecies(); i++){
	      if(oe.orthonode[i] != NULL){
		if(oe.orthonode[i]->n_type == sampled){
		  cout << ((State*)oe.orthonode[i]->item)->apostprob << ",";
		}else{
		  cout << "0,";
		}
	      }
	    }
	    cout << endl;
        }
    }
   cout.rdbuf(stdout); // reset to standard output again 
}

 
void GeneMSA::collect_features(int species, list<OrthoExon> *hects, SpeciesGraph *speciesgraph){

  if (!exoncands.empty() && exoncands[species] != NULL){
    list<ExonCandidate*>* ec = exoncands[species];
    for (list<ExonCandidate*>::iterator ecit = ec->begin(); ecit != ec->end(); ++ecit){
      stringstream key;
      key << "CDS\t" << getSeqID(species)<< "\t";
      if (getStrand(species) == plusstrand) {
	key << (*ecit)->begin + offsets[species] + 1 << "\t" << (*ecit)->end + offsets[species] + 1 << "\t";
      } else {
	int chrLen = rsa->getChrLen(species, getSeqID(species));
	key << chrLen - ((*ecit)->end + offsets[species]) << "\t" << chrLen - ((*ecit)->begin+ offsets[species]) << "\t";
      }
      key << ((isPlusExon((*ecit)->type) == (getStrand(species) == plusstrand))? '+' : '-');
      key << "\t" << (*ecit)->gff3Frame();
      
      string k = key.str();
      //cout << "EC: " << k << endl;
      unordered_map<string, pair<int, vector<double> > >::iterator got = Constant::logReg_feature.find(k);
      if ( got == Constant::logReg_feature.end() ){
	vector<double> feature(14,0);
	feature[0] = (*ecit)->end - (*ecit)->begin + 1;   // exon length
      
	pair<int, vector<double> > p;
	p = make_pair(-1, feature);
	pair<string, pair<int, vector<double> > > entry;
	entry = make_pair(k, p);
	Constant::logReg_feature.insert(entry);
      }
    }
  }
  
  for(list<OrthoExon>::iterator oeit = hects->begin(); oeit != hects->end(); ++oeit){
    ExonCandidate *ec = oeit->orthoex.at(species);
    if(ec == NULL)
      continue;
    stringstream key;
    key << "CDS\t" << getSeqID(species) << "\t" << ec->begin + offsets[species]+1 << "\t" << ec->end + offsets[species]+1 << "\t";
    if(isPlusExon(ec->type))
      key << "+";
    else
      key << "-";
    key << "\t" << ec->gff3Frame();
    //cout << "OE: " << key.str() << endl;
    unordered_map<string, pair<int, vector<double> > >::iterator got = Constant::logReg_feature.find(key.str());
    if ( got == Constant::logReg_feature.end() ){
      throw ProjectError("ortho exon is not an exon candidate!");
    }else{
      vector<double>* feature = &got->second.second;
      if(oeit->getEomega()       > 0){ (*feature)[3] = oeit->getEomega(); }        // omega
      if(oeit->getVarOmega()     > 0){ (*feature)[4] = oeit->getVarOmega(); }      // variance of omega
      if(oeit->getConsScore()    > 0){ (*feature)[5] = oeit->getConsScore(); }     // conservation
      if(oeit->getDiversity()    > 0){ (*feature)[6] = oeit->getDiversity(); }     // diversity
      if(oeit->getContainment()  > 0){ (*feature)[7] = oeit->getContainment(); }   // containment
      if(oeit->numExons()        > 0){ (*feature)[8] = oeit->numExons(); }         // number of species involved in OE
      (*feature)[9] = 1;                                                           // is ortho exon?
      if(oeit->getLeftExtOmega() > 0){ (*feature)[10] = oeit->getLeftExtOmega();}  // LeftExtOmega
      if(oeit->getLeftIntOmega() > 0){ (*feature)[11] = oeit->getLeftIntOmega();}  // LeftIntOmega
      if(oeit->getRightIntOmega()> 0){ (*feature)[12] = oeit->getRightIntOmega();} // RightIntOmega
      if(oeit->getRightExtOmega()> 0){ (*feature)[13] = oeit->getRightExtOmega();} // RightExtOmega
    }
  }
}

/*
 * 
 * This function obtains multiple sequence alignments (MSAs) and their label y=0,1, whether
 * it constitutes a real CDS or not in the reference species.
 */
void GeneMSA::getAllOEMsas(int species, list<OrthoExon> *hects, unordered_map<string,int> *ref_class, vector<AnnoSequence*> const &seqRanges){
    StringAlignment msa(0);
   
    for(list<OrthoExon>::iterator oeit = hects->begin(); oeit != hects->end(); ++oeit){
	ExonCandidate *ec = oeit->orthoex.at(species);
	if (ec == NULL)
	    continue; // can not consider alignments where the reference species has no exon candidate
	stringstream key;
	int beginStopOffset=0, endStopOffset=0; // to account for excluded stop codons
	if (hasStopCodon(ec->type) && Gene::stopCodonExcludedFromCDS){
	    if (isPlusExon(ec->type))
	        endStopOffset = -3;
	    else
                beginStopOffset = 3;
	}
	
	key << "CDS\t" << getSeqID(species) << "\t" << ec->begin + offsets[species] + 1 + beginStopOffset
	    << "\t" << ec->end + offsets[species] + 1 + endStopOffset << "\t";
	if (isPlusExon(ec->type))
	    key << "+";
	else
	    key << "-";
	key << "\t" << ec->gff3Frame();
	unordered_map<string, int>::iterator got = ref_class->find(key.str());
	bool y=0;
	if (got != ref_class->end())
	    y=1;

	try {
	    msa = getMsa(*oeit, seqRanges);
	    cout << "\ny=" << y << "\tOE" << oeit->ID << ": " << key.str() << endl
		 << msa << endl;
	} catch (...) {}
    }
}

StringAlignment GeneMSA::getMsa(OrthoExon const &oe, vector<AnnoSequence*> const &seqRanges, size_t flanking) {
    int k = alignment->numRows();
    int aliStart = oe.getAliStart() - flanking;
    int aliEnd = oe.getAliEnd() + flanking;
    int aliLen = aliEnd - aliStart + 1;
    int gaplen, matchlen, loverhang, prevAliEnd, insertlen;
    StringAlignment msa(k);
    MsaInsertion msains;
    list<MsaInsertion> insList;

    for (size_t s=0; s<k; s++){
	if (alignment->rows[s] == NULL || oe.orthoex[s] == NULL)
	    continue;
	AlignmentRow *row = alignment->rows[s];
	vector<fragment>::const_iterator prev = row->frags.end(),
            from = row->frags.begin(); // this could be more efficient exploiting sortedness

	// search first fragment that is not strictly to the left of the alignment start
	while (from != row->frags.end() && from->aliPos + from->len < aliStart)
	    from++;
        
	prevAliEnd = aliStart - 1;
	while (from != row->frags.end() && from->aliPos <= aliEnd){
            // are there unaligned insertions?
            if (prev != row->frags.end() && from->chrPos > prev->chrPos + prev->len){ // insertion
                msains.s = s;
                msains.insertpos = msa.rows[s].length();
                insertlen = from->chrPos - prev->chrPos - prev->len;
                msains.insert = string(seqRanges[s]->sequence + prev->chrPos + prev->len - offsets[s], insertlen);
                /*
                  cout << "insert of length " << insertlen << " found at "
                     << "s= " << s << " pos " << msains.insertpos << " i.e. "
                     << msains.insert << " inserted after " <<
                     msa.rows[s] << endl;
                */
                insList.push_back(msains);
            }
            // insert gap characters between previous and this fragment
            gaplen = from->aliPos - prevAliEnd - 1;
            if (gaplen > 0)
                msa.rows[s] += string(gaplen, '-');
	    
            loverhang = (from->aliPos < aliStart)? aliStart - from->aliPos : 0;
            matchlen = from->len - loverhang;
            if (matchlen > aliEnd - from->aliPos - loverhang + 1)
                matchlen = aliEnd - from->aliPos - loverhang + 1;
            msa.rows[s] += string(seqRanges[s]->sequence + from->chrPos + loverhang - offsets[s], matchlen);
            prevAliEnd = from->aliPos + loverhang + matchlen - 1;
            prev = from++;
	}
	if (msa.rows[s].size() < aliLen)
            msa.rows[s] += string(aliLen - msa.rows[s].size(), '-');
    }
    if (!insList.empty()){
        // cout << "msa before insertions:\n" << msa << endl;
        // inserts could be long, limit their length to 9
        msa.insert(insList, 9);
        // cout << "msa after  insertions:\n" << msa << endl;
    }
    try {
	msa.computeLen();
    } catch (length_error &e){
	cerr << e.what() << endl << msa << endl;
    }
    msa.removeGapOnlyCols();
    return msa;
}


/** 
 * Two codons are considered aligned, when all 3 of their bases are aligned with each other.
 * Note that not all bases of an ExonCandidate need be aligned.
 * example input (all ECs agree in phase at both boundaries)
 *        
 *                 a|c - - t t|g a t|g t c|g a t|a a 
 *                 a|c - - c t|a a - - - c|a n c|a g
 *                 g|c g - t|t g a|- g t c|g a c|a a
 *                 a|c g t|t t g|a t - t|c g a|c - a
 *                 a|c g - t|t g a|t g t|t g a|- a a
 *                   ^                       ^
 * firstCodonBase    |                       | lastCodonBase (for last species)
 * example output: (stop codons are excluded for single and terminal exons)
 *                 - - -|c t t|- - -|g t c|- - -|g a t
 *                 - - -|c c t|- - -|- - -|- - -|a n c
 *                 c g t|- - -|t g a|g t c|- - -|g a c
 *                 - - -|- - -|- - -|- - -|c g a|- - -
 *                 c g t|- - -|t g a|- - -|t g a|- - -
 *
 */
vector<string> GeneMSA::getCodonAlignment(OrthoExon const &oe, vector<AnnoSequence*> const &seqRanges,
					  const vector<vector<fragment>::const_iterator > &froms, map<unsigned, vector<int> > *alignedCodons, bool generateString, vector<vector<int> > *posStoredCodons, ofstream *codonAli) {
    //printSingleOrthoExon(oe, false);
  //   cout<<"generate codonAlignment for OE("<<oe.ID<<") start: "<<oe.getAliStart()<<" end: "<<oe.getAliEnd()<<" RFC: "<<printRFC(oe.getRFC(offsets))<<endl;
    int k = alignment->rows.size();
    vector<string> rowstrings(k, "");
    // consider only codon columns with a number of codons at least this fraction of the nonempty rows 
   
    /*
     * Store for each aligned codon the triplet of alignment columns encoded in a single long integer
     * key = aliPosOf1stBase * 2^8 + gapsTo2ndBase * 2^4 + gapsTo3rdBase
     * This assumes even on a rare 32 bit machine only that the alignment is shorter than 16,777,216,
     * and that gaps within a codon are at most 15bp. Where this is violated, wrong codon alignments
     * may happen.
     * Values of alignedCodons are pairs of 1) species index s and 2) the chromosomal position of the 
     * first codon base.
     */
    map<unsigned, vector<int> > ac;
    if(alignedCodons == NULL)
	alignedCodons = &ac;
    map<unsigned, vector<int> >::iterator acit;
    map<unsigned, vector<int> > codonAliOE;
    map<unsigned, vector<int> >::iterator oeit;    

    long aliPosOf1stBase, aliPosOf2ndBase, aliPosOf3rdBase;
    int chrCodon1;
    // Map all codons to alignment positions, where possible. This search in LINEAR in the length
    // of all exon candidates and the number of all alignment fragments.
    for (size_t s=0; s<k; s++){
	if (alignment->rows[s] == NULL || oe.orthoex[s] == NULL)
	    continue;
	//	cout<<"species "<<s<<"|"<<getSeqID(s);
	AlignmentRow *row = alignment->rows[s];
	ExonCandidate *ec = oe.orthoex[s];
	int firstCodonBase = oe.getStartInWindow(s) + ( (offsets[s] + ec->getFirstCodingBase() - oe.getStartInWindow(s)) % 3);
	int lastCodonBase = oe.getEndInWindow(s) - ( (oe.getEndInWindow(s) - (offsets[s] + ec->getLastCodingBase())) % 3);
	int firstCodonBaseOE = ec->getFirstCodingBase() + offsets[s];
	int lastCodonBaseOE = ec->getLastCodingBase() + offsets[s];
	//cout << " firstCodingBase: " << ec->getFirstCodingBase() + offsets[s] << " lastCodingBase: " << ec->getLastCodingBase() + offsets[s] << endl;
	//cout<<" firstCodonBase in window: "<<firstCodonBase<<" lastCodonBase in window: "<<lastCodonBase<<" frame: "<<(firstCodonBase % 3)<<endl;
	if ((lastCodonBase - firstCodonBase + 1) % 3 != 0)
	    throw ProjectError("Internal error in getCodonAlignment: frame and length inconsistent.");
	/* TODO: safe computational time by storing already calculated positions
	if(posStoredCodons != NULL){ // calculate only those codons that have not been pocessed yet
	    if((*posStoredCodons)[s][firstCodonBase % 3] >= lastCodonBase )
		continue;
	    else{
		if(firstCodonBase < (*posStoredCodons)[s][firstCodonBase % 3] + 1)
		    firstCodonBase = (*posStoredCodons)[s][firstCodonBase % 3] + 1;
		(*posStoredCodons)[s][firstCodonBase % 3] = lastCodonBase;
	    }
	    }*/
	vector<fragment>::const_iterator from = froms[s];
	for(chrCodon1 = firstCodonBase; chrCodon1 <= lastCodonBase - 2; chrCodon1 += 3){
	    // chrCodon1 is the chromosomal position of the first base in the codon
	    // go the the first fragment that may contain chrCodon1
	    while (from != row->frags.end() && from->chrPos + from->len - 1 < chrCodon1)
		++from;
	    if (from == row->frags.end())
		break; // have searched beyond the codon => finished
	    aliPosOf1stBase = row->getAliPos(chrCodon1, from);
	    if (aliPosOf1stBase >= 0){ // first codon base mappable
		aliPosOf2ndBase = row->getAliPos(chrCodon1 + 1, from);
		if (aliPosOf2ndBase >= 0){
		    aliPosOf3rdBase = row->getAliPos(chrCodon1 + 2, from);
		    if (aliPosOf3rdBase >= 0){
			// all the codon bases were mappable, store codon in the map
			long key = (aliPosOf1stBase << 8) + ((aliPosOf2ndBase - aliPosOf1stBase - 1) << 4)
			    + (aliPosOf3rdBase - aliPosOf2ndBase - 1);
			//cout << "key:" << chrCodon1 << " " << aliPosOf1stBase << " " << (aliPosOf2ndBase - aliPosOf1stBase - 1) << " " << (aliPosOf3rdBase - aliPosOf2ndBase - 1) << " / " << key << endl;
			acit = alignedCodons->find(key);
			if (acit == alignedCodons->end()){ // insert new vector
			  vector<int> cod(k, -1); // -1 missing codon
			    cod[s] = chrCodon1;
			    alignedCodons->insert(pair<unsigned,vector<int> >(key, cod));
			} else {// append new entry to existing list
			  if (acit->second[s] >= 0){ // remove this later
			    //  cout<<"species "<<s<<" : want to insert: "<<chrCodon1<<", already in: "<<acit->second[s]<<endl;
			    //  cout<<"Wow! Another codon has mapped to the same key!"<<endl;
			  }
			    acit->second[s] = chrCodon1;
			}
			if(firstCodonBaseOE >= firstCodonBase && lastCodonBaseOE <= lastCodonBase){
			  oeit = codonAliOE.find(key);
			  if(oeit == codonAliOE.end()){
			    vector<int> cod(k, -1); // -1 missing codon
			    cod[s] = chrCodon1;
			    codonAliOE.insert(pair<unsigned,vector<int> >(key, cod));
			  }else{
			    oeit->second[s] = chrCodon1;
			  }
			}
		    }
		}
	    }
	    // cout<<"codon of OE("<<oe.ID<<"), chrom Pos/RFC : "<< chrCodon1<<" / " <<(chrCodon1 % 3)<<endl;
	}
	//cout<<endl;

    }
    if(codonAli->is_open())
      generateString = true;
 
    if(generateString){
      
      // want to print all ortho exon alignments
      //float minAlignedCodonFrac = 0.3;
      //int m = alignment->numFilledRows();
      //int minAlignedCodons = (m * minAlignedCodonFrac > 2)? m * minAlignedCodonFrac + 0.9999 : 2;
	// Must have at least 'minAlignedCodons' codons in any codon column

	/*
	 * Create one codon alignment column for each key to which at least minAlignedCodons mapped
	 */ 
	for (acit = codonAliOE.begin(); acit != codonAliOE.end(); ++acit){
	  for (size_t s=0; s<k; s++){
	    chrCodon1 = acit->second[s]; // sequence position
	    if (chrCodon1 >= 0)
	      rowstrings[s] += string(seqRanges[s]->sequence + chrCodon1 - offsets[s], 3)+" ";
	    else 
	      rowstrings[s] += "--- ";
	  }
	}

	if (!isOnFStrand(oe.getStateType())){ // reverse complement alignment
	  //	  cout<<"reverse compliment"<<endl;
	  for (size_t s=0; s<k; s++)
		reverseComplementString(rowstrings[s]); 
	}
 
	//cout << "codon alignment:" << endl;

	vector<string> speciesNames;
	tree->getSpeciesNames(speciesNames);

          int maxSnameLen = rsa->getMaxSnameLen();
	  int maxSeqIDLen = alignment->getMaxSeqIdLen();
	  int numSp = 0;
	  for (size_t s=0; s<k; s++)
            if(oe.orthoex[s])
	      numSp++;
	  
	  *codonAli << "\t" << numSp << "\t" << rowstrings[0].size() - rowstrings[0].size()/4 << endl;
	  for (size_t s=0; s<k; s++){
	    if(oe.orthoex[s]){
	      string st = speciesNames[s] + "." + getSeqID(s) + ":" + to_string(oe.orthoex[s]->getStart() + offsets[s] + 1) + "-" + to_string(oe.orthoex[s]->getEnd() + offsets[s] + 1) + "(frame=" + to_string(oe.orthoex[s]->gff3Frame()) + ",type=" + to_string(oe.orthoex[s]->getStateType()) + ",ID=" + to_string(oe.ID) + ")";  
	      *codonAli << setw(maxSeqIDLen + maxSnameLen + 30) << left << st << "\t" << rowstrings[s] << endl;
	    }
	  }
	  *codonAli << endl;
	  
    }
    //cout<<"exiting getCodonAlignment"<<endl;
    return rowstrings;
}


// computes and sets the Omega = dN/dS attribute to all OrthoExons
void GeneMSA::computeOmegas(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, PhyloTree *ctree) {
    // int subst = 0;
    // Initialize for each species the first fragment froms[s] that 
    // is not completely left of the current OrthoExon.
    // This exploits the fact that both fragments and OrthoExons are sorted
    // left-to-right in the alignment.
  

    vector<vector<fragment>::const_iterator > froms(numSpecies());
    for (size_t s=0; s < numSpecies(); s++)
	if (alignment->rows[s])
	    froms[s] = alignment->rows[s]->frags.begin();
    for (list<OrthoExon>::iterator oe = orthoExonsList.begin(); oe != orthoExonsList.end(); ++oe){
	// move fragment iterators to start of exon candidates
	for (size_t s=0; s < numSpecies(); s++)
	    if (alignment->rows[s] && oe->orthoex[s])
		while(froms[s] != alignment->rows[s]->frags.end() 
		      && froms[s]->chrPos + froms[s]->len - 1 < offsets[s] + oe->orthoex[s]->getStart())
		    ++froms[s];
		
	vector<string> rowstrings = getCodonAlignment(*oe, seqRanges, froms);
	double omega;
	// TODO: scale branch lengths to one substitution per codon per time unit
	// cout << "OE" << endl;
	//	printSingleOrthoExon(*oe, false);
	// int subst = 0;
	// omega = codonevo->estOmegaOnSeqTuple(rowstrings, ctree, subst);
	//	cout << "omega=" << omega << endl;
	// oe->setSubst(subst);
	omega = codonevo->graphOmegaOnCodonAli(rowstrings, ctree);
	oe->setOmega(omega);
    }
}


// data to store in the map aliPos. for positions in the alignment store start and end of orthoExons and bitvectors
struct posElements{
    vector<OrthoExon*> oeStart;
    vector<OrthoExon*> oeEnd; 
};

string printBV(bit_vector bv){
    string bvs="";
    for(bit_vector::iterator b = bv.begin(); b != bv.end(); b++){
	if(*b)
	    bvs+= "1|";
	else
	    bvs+= "0|";
    }
    return bvs;
}

string printRFC(vector<int> rfc){
    string rfcs="";
    for(int i=0; i<rfc.size(); i++){
	rfcs+=itoa(rfc[i])+"|";
    }
    return rfcs;
}


void printTest(map<unsigned, vector<int> > alignedCodons, map<int, posElements> aliPos){

  cout<<"--- list of keys for codon algnment ---"<<endl;
    unordered_map<bit_vector, int, boost::hash<bit_vector>> bvCount;
    map<int, posElements>::iterator aliPosIt = aliPos.begin();
    for(map<unsigned, vector<int> >::iterator codonIt = alignedCodons.begin(); codonIt != alignedCodons.end(); codonIt++){
	cout<<"codon position in alignment(1): "<<(codonIt->first >> 8)<<" oe start or end in alignment(2): "<<aliPosIt->first<<endl;
	map<unsigned, vector<int> >::iterator nextCodonIt = codonIt;
	nextCodonIt++;
	while((unsigned)aliPosIt->first >= (codonIt->first >> 8) && (unsigned)aliPosIt->first <= (nextCodonIt->first >> 8)){
	    for(int i=0; i<aliPosIt->second.oeStart.size(); i++)
	      cout<<"oe " << aliPosIt->second.oeStart[i]->ID << " Start: "<<aliPosIt->second.oeStart[i]->getAliStart()<<":"<<aliPosIt->second.oeStart[i]->getAliEnd()<<endl;
	    for(int i=0; i<aliPosIt->second.oeEnd.size(); i++)
	      cout<<"oe " << aliPosIt->second.oeEnd[i]->ID << " End: "<<aliPosIt->second.oeEnd[i]->getAliStart()<<":"<<aliPosIt->second.oeEnd[i]->getAliEnd()<<endl;
	    aliPosIt++;
	}
    }
    cout << "--------------------------------------" << endl;
}

cumValues* GeneMSA::findCumValues(bit_vector bv, vector<int> rfc){
  //cout<<"in findCumValues, Parameters are "<<printBV(bv)<<" : "<<printRFC(rfc)<<endl;
    for(int b = 0; b < bv.size(); b++)
	if(! bv[b])
	    rfc[b] = -1;
    //cout<<"in findCumValues, Parameters after reduction "<<printBV(bv)<<" : "<<printRFC(rfc)<<endl;
    //  vector<cumValues*> allMatchingRFC;

    unordered_map<bit_vector, vector<pair<vector<int>,cumValues>>, boost::hash<vector<bool>>>::iterator it = cumOmega.find(bv);
    if(it == cumOmega.end())
	throw ProjectError("Internal error in findCumValues(): requested bit vector does not exist!");
    //cout<<"bitvector in cumOmega: "<<printBV(it->first)<<endl;
    //cout<<"bitvector has "<<it->second.size()<<" rfcs"<<endl;
    for(int i = 0; i < it->second.size(); i++){
	// cout<<"rfc number "<<i<<endl;
	bool isRFC = true;
	//cout<<"rfc in cumOmega: "<<printRFC(it->second[i].first)<<endl;
	for(int j = 0; j < it->second[i].first.size(); j++){
	    if(it->second[i].first[j] != rfc[j] && rfc[j] != -1){
	      //cout<<"this rfc is not the one we're looking for"<<endl;
		isRFC = false;
		break;
	    }	
	}
	if(isRFC){
	  //cout<<"return cumValues"<<endl;
	    return &it->second[i].second;
	}
    }
    return NULL;
}


void GeneMSA::printCumOmega(){ 
    cout<<"-------------- cumOmega -------------------"<<endl;
    for(unordered_map<bit_vector, vector<pair<vector<int>, cumValues> >, boost::hash<bit_vector> >::iterator coit = cumOmega.begin(); coit != cumOmega.end(); coit++){
	cout<<printBV(coit->first)<<endl;
	for(int i = 0; i < coit->second.size(); i++){
	  cout<<"\t"<<printRFC(coit->second[i].first)<<endl<<"logliks:";
	     for(int j = 0; j < coit->second[i].second.logliks.size(); j++){
	       cout<<"\t"<<coit->second[i].second.logliks[j];
	     }
	     cout<<endl;
	}
    }
    cout<<"--------------------------------------------"<<endl;
}

void GeneMSA::computeOmegasEff(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, PhyloTree *ctree, ofstream *codonAli) {
    cout<<"computing omega for each ortho exon."<<endl;

    // treat forward and reverse strand separately (might be done more efficiently)
    for (int strnd=1; strnd>=0; strnd--){
	bool plusStrand = (bool) strnd;
	if (plusStrand){
	    cout << "--- processing ortho exons on forward strand ---" << endl;
	} else {
	    cout << "--- processing ortho exons on reverse strand ---" << endl;
	}
      
	vector<vector<fragment>::const_iterator > froms(numSpecies());
	for (size_t s=0; s < numSpecies(); s++)
	    if (alignment->rows[s])
		froms[s] = alignment->rows[s]->frags.begin();
        
	cumOmega.clear();
	codonOmega.clear();

	map<int, posElements> aliPos;  // contains all positions where either an orthoExon or a bitvector starts or ends
	map<unsigned, vector<int> > alignedCodons; 
	/*
	 * Store for each aligned codon the triplet of alignment columns encoded in a single long integer
	 * key = aliPosOf1stBase * 2^8 + gapsTo2ndBase * 2^4 + gapsTo3rdBase
	 * This assumes even on a rare 32 bit machine only that the alignment is shorter than 16,777,216,
	 * and that gaps within a codon are at most 15bp. Where this is violated, wrong codon alignments
	 * may happen.
	 * Values of alignedCodons are pairs of 1) species index s and 2) the chromosomal position of the
	 * first codon base.
	 */

	alignedCodons.insert(pair<unsigned, vector<int> >(0,vector<int>(numSpecies(),-1))); // guaratee that codon alignment starts before first OrthoExon
	vector<vector<int> > posStoredCodons(numSpecies(),vector<int>(3,0)); // stores the position of the last codon aligned in getCodonAlignment() for each species and reading frame
    
	int cvID=0;

	cout << "generating codon alignment" << endl;
	//cout << "number of orthoExons: " << orthoExonsList.size() << endl;
	for (list<OrthoExon>::iterator oe = orthoExonsList.begin(); oe != orthoExonsList.end(); ++oe){
	  if(isOnFStrand(oe->getStateType()) != plusStrand)
	    continue;
	  //printSingleOrthoExon(*oe, false);
	  // store start and end information
	  bool aliStart = true;
	  oe->firstAlignedPos.resize(numSpecies());
	  oe->lastAlignedPos.resize(numSpecies());
	  // this iteration over two iterators was more elegant before, but not understood by gcc 4.4.7 in Santa Cruz
	  //	old code equivalent to the next 7 lines:
	  // for (map<int, posElements>::iterator aliPosIt :  { aliPos.find(oe->getAliStart()), aliPos.find(oe->getAliEnd()) }) { 
	  
	  // chop window at borders of alignment
	  int windowStart = min(oe->getAliStart(), Constant::oeExtensionWidth);
	  int windowEnd = min(alignment->aliLen - oe->getAliEnd(), Constant::oeExtensionWidth);
	  map<int, posElements>::iterator start = aliPos.find(oe->getAliStart());
	  map<int, posElements>::iterator end = aliPos.find(oe->getAliEnd());
	  // add  two windows that extend the OE at both bounadries
	  map<int, posElements>::iterator leftBoundaryWindowStart = aliPos.find(oe->getAliStart() - windowStart);
	  map<int, posElements>::iterator leftBoundaryWindowEnd = aliPos.find(oe->getAliStart() - 1);
	  map<int, posElements>::iterator rightBoundaryWindowStart = aliPos.find(oe->getAliEnd() + 1);
	  map<int, posElements>::iterator rightBoundaryWindowEnd = aliPos.find(oe->getAliEnd() + windowEnd);
	  // add two windows inside OE at both boundaries
	  map<int, posElements>::iterator leftBoundaryWindowStartInside = start;
          map<int, posElements>::iterator leftBoundaryWindowEndInside = aliPos.find(oe->getAliStart() + min(oe->getAliLen(), Constant::oeExtensionWidth));
          map<int, posElements>::iterator rightBoundaryWindowStartInside = aliPos.find(oe->getAliEnd() - min(oe->getAliLen(), Constant::oeExtensionWidth));
          map<int, posElements>::iterator rightBoundaryWindowEndInside = end;


	  // order is important! iterative start and end positions are assumed
	  list<pair<map<int, posElements>::iterator, int> > tlist; 
	  tlist.push_back(make_pair(start, oe->getAliStart()));
	  tlist.push_back(make_pair(end, oe->getAliEnd()));
	  tlist.push_back(make_pair(leftBoundaryWindowStart, oe->getAliStart() - windowStart));
	  tlist.push_back(make_pair(leftBoundaryWindowEnd, oe->getAliStart() - 1));
          tlist.push_back(make_pair(rightBoundaryWindowStart, oe->getAliEnd() + 1));
          tlist.push_back(make_pair(rightBoundaryWindowEnd, oe->getAliEnd() + windowEnd));
	  tlist.push_back(make_pair(leftBoundaryWindowStartInside, oe->getAliStart()));
	  tlist.push_back(make_pair(leftBoundaryWindowEndInside, oe->getAliStart() + min(oe->getAliLen(), Constant::oeExtensionWidth)));
	  tlist.push_back(make_pair(rightBoundaryWindowStartInside, oe->getAliEnd() - min(oe->getAliLen(), Constant::oeExtensionWidth)));
	  tlist.push_back(make_pair(rightBoundaryWindowEndInside, oe->getAliEnd()));

	  for (list<pair<map<int, posElements>::iterator, int> >::iterator lit = tlist.begin(); lit != tlist.end(); ++lit){
	    pair<map<int, posElements>::iterator, int> aliPosIt = *lit;

	    if(aliPosIt.first == aliPos.end()){
	      posElements pe;
	      pair<map<int, posElements>::iterator, bool> insertResult;
	      insertResult = aliPos.insert(pair<int, posElements>(aliPosIt.second, pe));
	      aliPosIt.first = insertResult.first;
	    }
	    if(aliStart){
	      aliPosIt.first->second.oeStart.push_back(&(*oe));
	    }else{
	      aliPosIt.first->second.oeEnd.push_back(&(*oe));
	    }
	    aliStart = !aliStart;
	  }

	  
    
	  // generate codon alignments
	  // move fragment iterators to start of exon candidates                                                                        
	  //cout << "+++OE " << oe->ID << " - coordinates in alignment: " << oe->getAliStart() << ":" << oe->getAliEnd() << endl;
	  for (size_t s=0; s < numSpecies(); s++){
	    if (alignment->rows[s] && oe->orthoex[s]){
	      //cout << "EC of species " << s << " = " << oe->orthoex[s]->begin + offsets[s] << ":" << oe->orthoex[s]->end + offsets[s] << endl;
	      //cout << "alignment ends at: " << alignment->rows[s]->aliEnd() << endl;

	      int ww = min(oe->getAliStart(), Constant::oeExtensionWidth);	    
	      //cout << "windowsize left: " << ww << endl;
	      oe->firstAlignedPos[s] = alignment->rows[s]->getChrPos(oe->getAliStart() - ww, froms[s]);

	      while(oe->firstAlignedPos[s] < 0){
		ww--;
		oe->firstAlignedPos[s] = alignment->rows[s]->getChrPos(oe->getAliStart() - ww, froms[s]);
	      }
	      vector<fragment>::const_iterator to = froms[s];
	      ww = min(alignment->rows[s]->aliEnd() - oe->getAliEnd(), Constant::oeExtensionWidth);
	      //cout << "windowsize right: " << ww << endl;
	      oe->lastAlignedPos[s] = alignment->rows[s]->getChrPos(oe->getAliEnd() + ww, to);
              //cout << "chrom. position of first/last aligned base of OE: " << oe->firstAlignedPos[s] << "/" << oe->lastAlignedPos[s] << endl;
	      while(oe->lastAlignedPos[s] < 0){
		ww--;
		oe->lastAlignedPos[s] = alignment->rows[s]->getChrPos(oe->getAliEnd() + ww, to);
	      }
	    }
	    /*if (alignment->rows[s] && oe->orthoex[s])
	      while(froms[s] != alignment->rows[s]->frags.end()
	      && froms[s]->chrPos + froms[s]->len - 1 < lbwPos)
	      ++froms[s];
	    */
	  }
	    
	  getCodonAlignment(*oe, seqRanges, froms, &alignedCodons, false, &posStoredCodons, codonAli);
	}
  
	/* print aliPos
	  cout << "============ aliPos ==============" << endl;
	  for(map<int, posElements>::iterator aliPosIt = aliPos.begin(); aliPosIt != aliPos.end(); aliPosIt++){
	    cout << "## " << aliPosIt->first << endl << "starts" << endl;
	    for(int i=0; i < aliPosIt->second.oeStart.size(); i++){
	      cout << aliPosIt->second.oeStart[i]->ID << "\t";
	    }
	    cout << endl << "ends" << endl;
	    for(int i=0; i < aliPosIt->second.oeEnd.size(); i++){
              cout << aliPosIt->second.oeEnd[i]->ID << "\t";
            }
            cout << endl;
	  }
	*/
  
	cout<<"Merge processing: Traverse alignment left to right"<<endl;

	unordered_map<bit_vector, int, boost::hash<bit_vector>> bvCount;

	map<int, posElements>::iterator aliPosIt = aliPos.begin();
	if(aliPosIt == aliPos.end()){
	  cout<<"No orthoExons on "; 
	  if(plusStrand)
	    cout<<"forward";
	  else
	    cout<<"reverse";
	  cout<<" strand in current gene range!"<<endl;
	  continue;
	}
       
	// walk through codon alignment left to right
	for(map<unsigned, vector<int> >::iterator codonIt = alignedCodons.begin(); codonIt != alignedCodons.end(); codonIt++){
	  /*
	  cout<<"++++codon: "<<(codonIt->first >> 8)<<endl<<"chrom Pos / RFC : "<<endl;
	    for(vector<int>::iterator cit=codonIt->second.begin(); cit!=codonIt->second.end(); cit++){
	      cout<<*cit<<" / "<<(*cit % 3)<<endl;
	    }
	    cout<<"active bit vectors: ";
	  
	    for(unordered_map<bit_vector, int, boost::hash<bit_vector>>::iterator bit=bvCount.begin(); bit!=bvCount.end(); bit++){
	      if(bit->second > 0)
		cout<<printBV(bit->first)<<":"<<bit->second<<"\t";
	    }
	    cout<<endl;
	  */
	  // update bit_vector constellation
	  // ortho Exon starts or ends before current alignment position
	  while((unsigned)aliPosIt->first <= (codonIt->first >> 8) ){
	    //cout<<"next position of aliPos "<<aliPosIt->first<<endl;
	    
	    unordered_map<bit_vector, int, boost::hash<bit_vector>>::iterator bvit;
	    unordered_map<bit_vector, vector<pair<vector<int>, cumValues> >, boost::hash<bit_vector> >::iterator coit;
	    // process all ortho exons that start
	    for(int i=0; i<aliPosIt->second.oeStart.size(); i++){
	      //cout<<"##################ortho exon ("<<aliPosIt->second.oeStart[i]->ID<<") starts: "<<aliPosIt->second.oeStart[i]->getAliStart()<<":"<<aliPosIt->second.oeStart[i]->getAliEnd()<<endl;

	      if (false){
	        cout<<"chromosomal position of each exon:"<<endl;
	        for(int j=0; j<aliPosIt->second.oeStart[i]->orthoex.size(); j++){
	          if(aliPosIt->second.oeStart[i]->orthoex[j]){
		    cout<<"species "<<j<<"\t"<<aliPosIt->second.oeStart[i]->orthoex[j]->begin<<"\toffset: "<<offsets[j]<<"\t"<<aliPosIt->second.oeStart[i]->orthoex[j]->end<<"\toffset: "<<offsets[j]<<"\t"<<aliPosIt->second.oeStart[i]->orthoex[j]->getStateType()<<"\t";

		  ExonCandidate *ec = aliPosIt->second.oeStart[i]->orthoex[j];
		  if (getStrand(j) == plusstrand){ // strand of alignment                                                                   		   cout << "start:" << ec->begin + offsets[j]+1 << "\tend:" << ec->end + offsets[j]+1;
		  } else {
		    int chrLen = rsa->getChrLen(j, getSeqID(j));
		     cout << "start:" << chrLen - (ec->end + offsets[j]) << "\tend:" << chrLen - (ec->begin + offsets[j]);
		  }
		   cout<<endl;

		}else
		    cout<<"species "<<j<<endl;
	        }
	      }
	      //cout<<"---bv of ortho exon:  "<<printBV(aliPosIt->second.oeStart[i]->getBV())<<endl<<"---rfc of ortho exon: "<<printRFC(static_cast<const OrthoExon*>(aliPosIt->second.oeStart[i])->getRFC(offsets))<<endl;
	      if(aliPosIt->second.oeStart[i] == NULL)
		throw ProjectError("Error in posElement.oestart: Pointer to orthoExon is NULL!");
	      bvit = bvCount.find(aliPosIt->second.oeStart[i]->getBV());
	      coit = cumOmega.find(aliPosIt->second.oeStart[i]->getBV());
	      
	      if(bvit == bvCount.end()){
		pair<unordered_map<bit_vector, int, boost::hash<bit_vector>>::iterator,bool> result = bvCount.insert(pair<bit_vector, int>(aliPosIt->second.oeStart[i]->getBV(),0));
		bvit = result.first;
		//if(! result.second)
		//  cerr<<"inserting bit vector "<<printBV(bvit->first)<<" into bvcount failed"<<endl;
		//else
		// cout<<"inserting bit vector "<<printBV(bvit->first)<<" into bvcount done"<<endl;
	        //}else{
		//cout<<"bit vector "<<printBV(bvit->first)<<" already inserted in bvcount"<<endl;
		//if(coit == cumOmega.end())
		//  throw ProjectError("Internal Error in computeOmegaEff(): bit_vector in bvcount but not yet in cumOmega!");
		
	      }
	      bvit->second++;
	      
	      // add reading frame combination if new one occurs
	      int currRFnum; //position in cumOmega vector
	      
	      if(coit == cumOmega.end()){
		vector<pair<vector<int>, cumValues> > vecPair;
		pair<unordered_map<bit_vector, vector<pair<vector<int>, cumValues> >, boost::hash<bit_vector> >::iterator, bool> result = cumOmega.insert(pair<bit_vector, vector<pair<vector<int>, cumValues> > >(aliPosIt->second.oeStart[i]->getBV(), vecPair));
		coit = result.first;
		//if(! result.second)
		//  cerr<<"inserting bit vector "<<printBV(coit->first)<<" into cumOmega failed"<<endl;
		//else
		//  cout<<"inserting bit vector "<<printBV(coit->first)<<" into cumOmega done"<<endl;
		//
	        //}else{
		//cout<<"bit vector "<<printBV(coit->first)<<"already exists in cumOmega"<<endl;

	      }
	      bool rfcIncluded = false;
	      // add new cumValue
	      cumValues cum(cvID);
	      cvID++;
	      //cout << "created cv " << cvID-1 << " with RFC " << printRFC(static_cast<const OrthoExon*>(aliPosIt->second.oeStart[i])->getRFC(offsets)) << endl;
	      pair<vector<int>, cumValues> oeRFC = make_pair(const_cast<const OrthoExon*>(aliPosIt->second.oeStart[i])->getRFC(offsets),cum);
	      for(int rf = 0; rf < coit->second.size(); rf++){
		if(oeRFC.first == coit->second[rf].first){
		  currRFnum = rf;
		  rfcIncluded = true;
		  //cout<<printRFC(oeRFC.first)<<" already in cumOmega ("<<printRFC(coit->second[rf].first)<<")"<<endl;
		  break;
		}
	      }
	      if(! rfcIncluded){
		coit->second.push_back(oeRFC);
		//cout<<"inserting RFC: "<<printRFC(oeRFC.first)<<" to bit_vector "<<printBV(coit->first)<<endl;
		currRFnum = coit->second.size() - 1;
	      }
	      
	      // store cumulative values at the beginning of an OrthoExon
	      cumValues cv = coit->second[currRFnum].second;
	      
	      aliPosIt->second.oeStart[i]->setOmega(&cv.logliks, codonevo, true);
	      if(Constant::computeNumSubs){
		aliPosIt->second.oeStart[i]->setSubst(cv.numSubs, true);
		//cout << "call setSubst() with cv " << cv.id << endl;
	      }
	    }
	    // process all ortho exons that end
	    for(int i=0; i<aliPosIt->second.oeEnd.size(); i++){
	      //cout<<"################ortho exon ("<<aliPosIt->second.oeEnd[i]->ID<<") ends: "<<aliPosIt->second.oeEnd[i]->getAliStart()<<":"<<aliPosIt->second.oeEnd[i]->getAliEnd()<<endl;

	      if (false){
		cout<<"chromosomal position of each exon:"<<endl;
		for(int j=0; j<aliPosIt->second.oeEnd[i]->orthoex.size(); j++){
		  if(aliPosIt->second.oeEnd[i]->orthoex[j])
		    cout<<"species "<<j<<"\t"<<aliPosIt->second.oeEnd[i]->orthoex[j]->begin<<"\toffset: "<<offsets[j]<<"\t"<<aliPosIt->second.oeEnd[i]->orthoex[j]->end<<"\toffset: "<<offsets[j]<<"\t"<<aliPosIt->second.oeEnd[i]->orthoex[j]->getStateType()<<endl;
		  else
		    cout<<"species "<<j<<endl;
		}
	      }
	      bvit = bvCount.find(aliPosIt->second.oeEnd[i]->getBV());
	      if(bvit == bvCount.end())
		throw ProjectError("Error in computeOmegaEff(): bit_vector ends that has never started!");
	      // calculate omega of ortho exon from cumulative sum
	      cumValues *cv = findCumValues(aliPosIt->second.oeEnd[i]->getBV(), const_cast<const OrthoExon*>(aliPosIt->second.oeEnd[i])->getRFC(offsets));
	      if(cv == NULL){
		cerr<<"cum Values has NULL pointer"<<endl;
	      }
	      aliPosIt->second.oeEnd[i]->setOmega(&cv->logliks, codonevo, false);
	      if(Constant::computeNumSubs){
		//cout << "call setSubst() with cv " << cv->id << endl;
		aliPosIt->second.oeEnd[i]->setSubst(cv->numSubs, false);
	      }
	      bvit->second--;
	    }
	    aliPosIt++;
	  }
          
	  // compute omega for current codon alignment
	  // generate array of strings representing one codon alignment
	  vector<string> codonStrings(numSpecies(),"");
	  int numCodons = 0;
	  vector<int> chrCodonPos(numSpecies(),-1);
	  for(size_t s = 0; s < numSpecies(); s++)
	    if (codonIt->second[s] >=0)
	      numCodons++;
	  if(numCodons >= 2){
	    vector<int> rfc(numSpecies(),-1); // reading frame combination of current codon alignment
	    for (size_t s = 0; s < numSpecies(); s++){
	      int chrCodon1 = codonIt->second[s]; // sequence position
	      if (chrCodon1 >= 0){
		codonStrings[s] = string(seqRanges[s]->sequence + chrCodon1 - offsets[s], 3);
		rfc[s] = chrCodon1 % 3;
		chrCodonPos[s] = chrCodon1;
	      }
	      else 
		codonStrings[s] = "---";
	    }
      
	    if (! plusStrand){ // reverse complement alignment
	      for (size_t s = 0; s < numSpecies(); s++)
		reverseComplementString(codonStrings[s]); 
	    }
	    /*
	      cout<<"codon alignment:"<<endl;
	      for (size_t s = 0; s < numSpecies(); s++){
	      cout<<"species "<<s<<"\t"<<codonStrings[s]<<"\t"<<chrCodonPos[s]<<"\t";
	      if (getStrand(s) == plusstrand){ // strand of alignment                                                                   		    
	      cout<<"start:" << chrCodonPos[s] - offsets[s]+1 + starts[s];
	      } else {
	      cout<<"start:" << ends[s] - (chrCodonPos[s] - offsets[s]+1);
	      }
	      cout<<endl;
	      }
	    */
	  
	    // for one alignment position compute omega for all aktive bit_vectors in the correct reading frame combination
	    //cout<<"compute omega for RFC "<<printRFC(rfc)<<endl;
	    //bool foundBV = false;
	    for(unordered_map<bit_vector, int, boost::hash<bit_vector>>::iterator bvit = bvCount.begin(); bvit != bvCount.end(); bvit++){
	      if(bvit->second == 0)
		continue;
	      //cout<<"next Bitvector in bvcount "<<printBV(bvit->first)<<":"<<bvit->second<<endl;
	      cumValues *cv = findCumValues(bvit->first, rfc);    
	      //cout<<"after findCumValues"<<endl;
	      if(cv != NULL){
		// call pruning algo only once for every codonStrings and store omega in map
		int subs = 0; // store number of substitutions
		vector<double> loglik;
		vector<string> cs = pruneToBV(&codonStrings, bvit->first);
		// scipt the next step if cs only consists of "---" entries
		if(cs[0] == "---" && adjacent_find(cs.begin(), cs.end(), not_equal_to<string>()) == cs.end())
		  continue;
		
		map<vector<string>, pair<vector<double>, int> >::iterator oit = computedCumValues.find(cs);
		if(oit==computedCumValues.end()){
		  if(Constant::computeNumSubs)
		    loglik = codonevo->loglikForCodonTuple(cs, ctree, tree, subs);
		  else
		    loglik = codonevo->loglikForCodonTuple(cs, ctree);
		  pair<vector<double>, int> store_cv = make_pair(loglik, subs);
		  computedCumValues.insert(pair<vector<string>,pair<vector<double>, int> >(cs,store_cv));
		}else{
		  loglik = oit->second.first;
		  subs = oit->second.second;
		}
		// calculate columnwise omega and store in appropriate data structure
		/*
		//cout << "calculate omega for codon " << (codonIt->first >> 8) << " ...";
		vector<int> pruned_rfc = pruneToBV(&rfc, bvit->first);
		//cout << "current (reduced) RFC " << printRFC(pruned_rfc) << endl; 
		map<bit_vector, map<vector<int>, vector<double> > >::iterator omegaIt = codonOmega.find(bvit->first);
	        if(omegaIt==codonOmega.end()){
		  map<vector<int>, vector<double> > currRFC;
		  vector<double> o(alignment->aliLen,-1);
		  o[codonIt->first >> 8] = omegaForCodonTuple(&loglik);
		  currRFC.insert(pair<vector<int>, vector<double> >(pruned_rfc, o));
		  codonOmega.insert(pair<bit_vector, map<vector<int>, vector<double> > >(bvit->first, currRFC));
		  int sum_of_rfc = 0;
		  for(map<bit_vector, map<vector<int>, vector<double> > >::iterator oi = codonOmega.begin(); oi != codonOmega.end(); oi++){
		    sum_of_rfc += oi->second.size();
		  }
		  cout << "size of codonOmega: " << codonOmega.size() << "\tsum of RFCs: " << sum_of_rfc << endl;
		}else{
		  map<vector<int>, vector<double> >::iterator rfcIt = omegaIt->second.find(pruned_rfc);
		  if(rfcIt == omegaIt->second.end()){
		    vector<double> o(alignment->aliLen,-1);
		    o[codonIt->first >> 8] = omegaForCodonTuple(&loglik);
		    omegaIt->second.insert(pair<vector<int>, vector<double> >(pruned_rfc, o));
		  }else{
		    if(rfcIt->second[codonIt->first >> 8] != -1)
		      cerr << "Warning: omega was already calculated for same alignment position of first codon, bit_vector and RFC!" << endl;
		    rfcIt->second[codonIt->first >> 8] = omegaForCodonTuple(&loglik);
		  }
		}
		*/
		//cout << "done" << endl;
		//cout<<"loglik of omega: "<<loglik<<endl;
		// store cumulative sum of omega, omega squared and one
		cv->addLogliks(&loglik);
		if(Constant::computeNumSubs)
		  cv->addNumSubs(subs);
		//foundBV=true;
		//printCumOmega();
	      }
	    }
	    //if(! foundBV)
	    // cerr<<"no Bitvector with given RFC found!"<<endl;
	  }
	}
	//cout << "+++ process orthoExon that start or end after the end of the codon alignment!"<<endl;
	int lastPos = aliPosIt->first;
	while(aliPosIt != aliPos.end()){ // process remaining orthoExon ends
	  //if(aliPosIt->second.oeStart.size() > 0)
	  //cerr<<"Warning: there are still orthoexon(s) beginning although codon alignment ended"<<endl; 
	  
	  for(int i=0; i<aliPosIt->second.oeEnd.size(); i++){
	    //cout<<"################ortho exon ("<<aliPosIt->second.oeEnd[i]->ID<<") ends after aliEnd: "<<aliPosIt->second.oeEnd[i]->getAliStart()<<":"<<aliPosIt->second.oeEnd[i]->getAliEnd()<<endl;
	    cumValues *cv;
	    if(aliPosIt->second.oeEnd[i]->getAliStart() >= lastPos){
	      cv = NULL;
	    }else{
	      cv = findCumValues(aliPosIt->second.oeEnd[i]->getBV(), const_cast<const OrthoExon*>(aliPosIt->second.oeEnd[i])->getRFC(offsets));
	    }
	    if(cv != NULL){
	      aliPosIt->second.oeEnd[i]->setOmega(&cv->logliks, codonevo, false);
	      if(Constant::computeNumSubs){
		aliPosIt->second.oeEnd[i]->setSubst(cv->numSubs, false);
                //cout << "call setSubst() with cv " << cv->id << endl;
	      }
	    }
	  }
	  aliPosIt++;
	}
	//printOmegaForCodon(outdir); need to parse outdir
    }
    cout<<"compute omegas done"<<endl;
}

// prune Codon Strings so that only aligned codons of species in bv are used for omega and numSubst calculation 
vector<string> GeneMSA::pruneToBV(vector<string> *cs, bit_vector bv){
  vector<string> cs_pruned(numSpecies(),"---");
  for(int s = 0; s < numSpecies(); s++)
    if(bv[s])
      cs_pruned[s] = (*cs)[s];
  return cs_pruned;
}

// prune RFC so that only active bitvector is used
vector<int> GeneMSA::pruneToBV(vector<int> *rfc, bit_vector bv){
  vector<int> rfc_pruned(numSpecies(),-1);
  for(int s = 0; s < numSpecies(); s++)
    if(bv[s])
      rfc_pruned[s] = (*rfc)[s];
  return rfc_pruned;
}



double GeneMSA::omegaForCodonTuple(vector<double> *loglik){

  if(loglik==NULL)
    cerr << "no likelihood was calculated in this pruning step" << endl;
  double sum = 0;
  int k = codonevo->getK();
  vector<double> postprobs(k, 0.0);
  double maxloglik = *max_element(loglik->begin(), loglik->end());
  for (int u=0; u < k; u++){
    sum += postprobs[u] = exp((*loglik)[u] - maxloglik) * codonevo->getPrior(u);
  }
  
  //cout << "posterior distribution and prior of omega" << endl;
  for (int u=0; u < k; u++){
	  // cout << codonevo->getOmega(u) << "\t" << postprobs[u] <<"\t"<<codonevo->getPrior(u)<< endl;
    postprobs[u] /= sum;
  }
  double omega = 0;
  
  //cout<<"---------------------------------------------------------------------------"<<endl;
  //cout<<"wi\t\tloglikOmegas\tmaxloglik\tpostprobs/sum\tprior\tsum\texp(loglik - maxloglik)"<<endl;
  
  for (int u=0; u < k; u++){
    omega += postprobs[u] * codonevo->getOmega(u);
    //cout<<codonevo->getOmega(u)<<"\t\t"<<loglikOmegas[u]<<"\t\t"<<maxloglik<<"\t\t"<<postprobs[u]<<"\t\t"<<codonevo->getPrior(u)<<"\t\t"<<sum<<"\t\t"<<exp(loglikOmegas[u] - maxloglik)<<endl;                                                     
  }
  
  return omega;
}

void GeneMSA::printOmegaForCodon(string outdir){
  
  for(map<bit_vector, map<vector<int>, vector<double> > >::iterator bvIt = codonOmega.begin(); bvIt != codonOmega.end(); bvIt++){
    int bvToDecimal = 0;
    for(int s = 0; s<numSpecies(); s++)
      bvToDecimal += pow(3,s) * bvIt->first[s]+1;
    string bvStr = "BV" + to_string(bvToDecimal);
    for(map<vector<int>, vector<double> >::iterator rfcIt = bvIt->second.begin(); rfcIt != bvIt->second.end(); rfcIt++){
      int rfcToDecimal = 0;
      for(int s = 0; s<numSpecies(); s++)
	rfcToDecimal += pow(4,s) * rfcIt->first[s]+1;
      string rfcStr = "RFC" + to_string(rfcToDecimal);
      consToWig(rfcIt->second, outdir + "omega." + bvStr + rfcStr);
    }
  }
}



// calculate a columnwise conservation score and output it (for each species) in wiggle format
void GeneMSA::calcConsScore(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, string outdir){

    vector<vector<fragment>::const_iterator > fragsit(numSpecies());
    vector<size_t> seqPos(numSpecies(),0); 
    for (size_t s=0; s < numSpecies(); s++){
	if (alignment->rows[s])
	    fragsit[s] = alignment->rows[s]->frags.begin();
    }
    vector<double> consScore; // vector of conservation scores for each position in the alignment
    for(size_t i = 0; i < alignment->aliLen; i++){
	int a=0,c=0,t=0,g=0; // number of bases of each type in one alignment column
	for(size_t j = 0; j < alignment->rows.size(); j++){
	    AlignmentRow *row = alignment->rows[j];
	    if(row && fragsit[j] != row->frags.end()){
		vector<fragment>::const_iterator it = fragsit[j];
		if( i >= it->aliPos && i <= it->aliPos + it->len - 1){ // character in i-th column, j-th row is not a gap
		    int pos = it->chrPos - offsets[j] + seqPos[j];
                    if(pos < 0 || pos >= seqRanges[j]->length)
                        throw ProjectError("Internal error in GeneMSA::printConsScore: trying to read position" + itoa(pos+1) + "in sequence " + seqRanges[j]->seqname + ".");
		    const char* base = seqRanges[j]->sequence + pos;
		    switch(*base){
		    case 'a': a++; break;
		    case 'c': c++; break;
		    case 't': t++; break;
		    case 'g': g++; break;
		    }
		    seqPos[j]++;
		}
		if(i == it->aliPos + it->len - 1){ // reached the end of the current fragment, move iterator to next fragment
		    ++fragsit[j];
		    seqPos[j]=0;
		}
	    }
	}
	double sc = calcColumnScore(a,c,t,g);
	consScore.push_back(sc);
    }
    // DEBUGGING Mario
    if (consScore.size() != alignment->aliLen)
	cerr << "consScore.size() != alignment->aliLen: " <<  consScore.size() << " != " << alignment->aliLen << endl;
    // calcluate conservation score for each HECT
    for (list<OrthoExon>::iterator oe = orthoExonsList.begin(); oe != orthoExonsList.end(); ++oe){
	double oeConsScore=0.0;
	int oeAliStart = oe->getAliStart();
	int oeAliEnd = oeAliStart + oe->getAliLen(); // the actual end, as oe->getAliLen  = end - start
	for(int pos = oeAliStart; pos <= oeAliEnd; pos++){
	    if (pos > alignment->aliLen || pos < 0)
		throw ProjectError("Internal error in printConsScore: alignment positions of HECTs and geneRanges are inconsistent.");
	    if (pos >= consScore.size()){
		cerr << "invalid access at " << pos << "\t" << consScore.size() <<  endl;
	    }
	    if (pos < consScore.size())
		oeConsScore += consScore[pos];
	}
	oeConsScore/=(oeAliEnd-oeAliStart+1); // average over all alignment columns within a HECT
	oe->setConsScore(oeConsScore);
	// conservation score for left boundary feature
	oeConsScore=0.0;
        int oeLeftBoundAliStart = max(oeAliStart - Constant::oeExtensionWidth, 0);
        int oeLeftBoundAliEnd = max(oeAliStart - 1, 0);
        for(int pos = oeLeftBoundAliStart; pos <= oeLeftBoundAliEnd; pos++){
	  if (pos > alignment->aliLen || pos < 0)
	    throw ProjectError("Internal error in printConsScore: alignment positions of HECTs and geneRanges are inconsistent.");
	  oeConsScore+=consScore[pos];
        }
        oeConsScore/=(oeLeftBoundAliEnd-oeLeftBoundAliStart+1); // average over all alignment columns within a HECT
        oe->setLeftConsScore(oeConsScore);
	// conservation score for right boundary feature
	oeConsScore=0.0;

        int oeRightBoundAliStart = min(oeAliEnd + 1, alignment->aliLen);
        int oeRightBoundAliEnd = min(oeAliEnd + 1 + Constant::oeExtensionWidth, alignment->aliLen - 1);
        for(int pos = oeRightBoundAliStart; pos <= oeRightBoundAliEnd; pos++){
            if (pos >= alignment->aliLen || pos < 0)
                throw ProjectError("Internal error in printConsScore: alignment positions of HECTs and geneRanges are inconsistent.");
            oeConsScore += consScore[pos];
        }
	
        oeConsScore/=(oeRightBoundAliEnd-oeRightBoundAliStart+1); // average over all alignment columns within a HECT

        oe->setRightConsScore(oeConsScore);

    }
    // output for each geneRange and each species a conservation track in wiggle format

    try {
	if (Properties::getBoolProperty( "/CompPred/printConservationWig" ))
	    consToWig(consScore, outdir);
    } catch (...) {}
}

// calculates a conservation score for a single alignment column
double GeneMSA::calcColumnScore(int a, int c, int t, int g){ // input: number of a,c,t and g's in one alignment column
    int N = numSpecies();
    /* frequency of the most frequent nucleotide type divided by the total number of species
     *double max = 0.0;
     *max = a > c ? a : c;
     *max = max > t ? max : t;
     *max = max > g ? max : g;
     *return max/N;
     */
    // entropy (base 4) weighted by the percentage of rows present, scores range from 0 (non-conserved) to 1 (conserved) 
    double sum=a+c+g+t;
    if(sum == 0)
	return 0.0;
    double probA=a/sum;
    double probC=c/sum;
    double probG=g/sum;
    double probT=t/sum;
    double entropy=0.0;
    if(probA > 0.0 && probA < 1.0)
	entropy-=(probA*log2(probA));
    if(probC > 0.0 && probC < 1.0)
	entropy-=(probC*log2(probC));
    if(probT > 0.0 && probT < 1.0)
	entropy-=(probT*log2(probT));
    if(probG > 0.0 && probG < 1.0)
	entropy-=(probG*log2(probG));
    return (1-(0.5*entropy))*(sum/N);
    //return (1-(0.5*entropy)); for getting all ortho exons that have no substitution (consScore = 0) 

}

// print conservation score to wiggle file
void GeneMSA::consToWig(vector<double> &consScore, string outdir){

    for (size_t i = 0; i < alignment->rows.size(); i++){
	AlignmentRow *row = alignment->rows[i];
	if(row){
	    string speciesname=rsa->getSname(i);
	    ofstream outfile(outdir + speciesname + ".wig",fstream::app);
	    if (outfile.is_open()){
		outfile<<"track type=wiggle_0 name=\""<<geneRangeID-1<<"\" description=\""<<alignment->getSignature()<<"\""<<endl;
		outfile<<"variableStep chrom="<<row->seqID<<endl;
		if(consScore.size() == 0)
		  continue;
		for (vector<fragment>::const_iterator it = row->frags.begin(); it != row->frags.end(); ++it){
		  int aliPos=it->aliPos;
		  int chrPos=it->chrPos;
		  while(chrPos < it->chrPos + it->len){
			outfile<<chrPos+1<<" "<<consScore[aliPos]<<endl;
			aliPos++;
			chrPos++;
		    }
		}
	    }
	    outfile.close();
	}
    }
}

void GeneMSA::closeOutputFiles(){
    for (int i=0; i<tree->numSpecies(); i++) {
        if (i < exonCands_outfiles.size() && exonCands_outfiles[i] && exonCands_outfiles[i]->is_open()) {
	    exonCands_outfiles[i]->close();
	    delete exonCands_outfiles[i];
	}
        if (geneRanges_outfiles_gff.size() > i && geneRanges_outfiles_gff[i] && geneRanges_outfiles_gff[i]->is_open()) {
	    geneRanges_outfiles_gff[i]->close();
	    delete geneRanges_outfiles_gff[i];
	}
        if (geneRanges_outfiles_bed.size() > i && geneRanges_outfiles_bed[i] && geneRanges_outfiles_bed[i]->is_open()) {
	    geneRanges_outfiles_bed[i]->close();
	    delete geneRanges_outfiles_bed[i];
	}
        if (i < orthoExons_outfiles.size() && orthoExons_outfiles[i] && orthoExons_outfiles[i]->is_open()) {
	    orthoExons_outfiles[i]->close();
	    delete orthoExons_outfiles[i];
	}
        if (omega_outfiles[i] && omega_outfiles[i]->is_open()) {
	    omega_outfiles[i]->close();
	    delete omega_outfiles[i];
	}
    }
}

void GeneMSA::comparativeSignalScoring(list<OrthoExon> &orthoExonsList){
    cout << "entering comparativeSignalScoring" << endl;
    ExonCandidate *ec;
    for (list<OrthoExon>::iterator oeit = orthoExonsList.begin(); oeit != orthoExonsList.end(); ++oeit){
	cout << "next OrthoExon:" << endl;
	// loop over all exon candidates belonging to this OrthoExon
	for (int s=0; s < numSpecies(); s++) {
	    ec = oeit->orthoex.at(s); 
	    if (ec != NULL) { // otherwise exon candidate is missing in  species s 
		cout << "next ExonCandidate on sequence " << getSeqID(s) << " from species number "  << s
		     <<  " (= " << rsa->getSname(s) << ")\t";
		cout << *ec; // prints begin, end, type, score, assScore, dssScore
		cout << " type as name = " << stateExonTypeIdentifiers[ec->type];
		cout << endl;
		// assScore is only relevant for internal and terminal exons on both strands
		// dssScore is only relevant for internal and initial exons on both strands
		// ec->getAssScore();
		// ec->getDssScore();
		// for now you will not need these methods, but later we may change the scores
		// according to what you find
		// ec->setAssScore(42);
		// ec->setDssScore(42);		
	    }
	}
    }
    cout << "exiting comparativeSignalScoring" << endl;
}

/* constructTree
 * creates, stores are returns the locus tree for the sequence tuple
 * The locus tree is specific to for the alignment and the aligned sequences.
 * It may be different (also in topology) from the species tree of which we
 * consider only one, because of duplications and paralogy.
 */
LocusTree *GeneMSA::constructTree(){
    // Charlotte Janas playground
    ltree = new LocusTree();
    // construct a tree from the alignment
    return ltree;
} 



