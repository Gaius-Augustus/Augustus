
/**********************************************************************
 * file:    genome.cc
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
#include "genome.hh"
#include "bitmasking.hh"

#include <iostream>
#include <iomanip>      // std::setprecision
#include <string.h>     // strstr()
#include <mutex>        // mutual exclusion (mutex) of concurrent execution (threads)

using namespace std;

int Genome::no_genomes=0;
mutex mu;

void Genome::destroyGeneList(){
    for(std::list<Gene*>::iterator it=genes.begin(); it!=genes.end(); it++)
	delete *it;
}
void Genome::destroyHintList(){
    for(std::list<GeneFeature*>::iterator it=hints.begin(); it!=hints.end(); it++)
	delete *it;
} 

void Genome::parseExtrinsicGFF(string gfffilename){
    ifstream ifstrm(gfffilename.c_str());
    if (ifstrm.is_open()){
	int line_no = 0;
	try{
	    while (ifstrm) {
		string line;
		getline(ifstrm, line);
		line_no++;
		if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		    continue;
		vector<string> tokens = splitString(line);
		if(tokens.size() != 9)
		    throw ProjectError("wrong number of columns.\n");
		
		map<string,int>::iterator it = seqnames.find(tokens[0]);
		if(it == seqnames.end()) // no gene on that sequence
		    continue;
		int seqid = it->second;
		long int start = atol(tokens[3].c_str());
		long int end = atol(tokens[4].c_str());
		Strand strand = getStrand(tokens[6]);
		
		/*
		 * find source of extrinsic info, specified in gff as: source=X or src=X
		 */
		const char *spos;
		string esource;
		spos = strstr(tokens[8].c_str(), "source=");
		if (spos)
		    spos += 7;
		if (!spos) {
		    spos = strstr(tokens[8].c_str(), "src=");
		    if (spos)
			spos += 4;
		}
		if (spos && isalpha(*spos)){
		    int skeylen=1;
		    while (isalpha(*(spos+skeylen))){
			skeylen++;
		    }
		    esource = string(spos, skeylen);
		} else {
		    throw ProjectError("No source specified (e.g. by source=M in the last column).\n");
		}
		/*
		 * find all gene features that are supported
		 * by that hint and update their extrinsic source
		 */
		if(tokens[2] == "CDS" || tokens[2] == "intron"){
		    int frame = getFrame(tokens[7]);
		    FeatureType type = (tokens[2] == "CDS")? CDS : intron;
		    if(type == intron){
			start--;
			end++;
		    }
		    SeqIntKey seqInt(start, end-start+1, type);
		    list<GeneFeature*> le = findSeqInt(seqInt.getKey(), seqid, strand);
		    bool isGF = false;
		    for(list<GeneFeature*>::iterator it=le.begin(); it!=le.end(); it++){
			if((*it)->sameFrame(frame)){
			    (*it)->setEvidence(esource); // update extrinsic source
			    isGF = true;
			}
		    }
		    if(!isGF){
			GeneFeature *gf = new GeneFeature(type, start, end, strand, frame);
			gf->setEvidence(esource);
			insertSeqInt(gf, seqid);
			hints.push_back(gf);		 
		    }
		}
	    }	        
	} catch( ProjectError& err ){
	    cerr << "Error in " << gfffilename << " in line " << line_no << endl;
	    throw err;
	}
    }
    else
        throw ProjectError("Could not open input gff file " + gfffilename + ".\n");
}

void Genome::parseGTF(string gtffilename){
    ifstream ifstrm(gtffilename.c_str());
    if (ifstrm.is_open()){
	int line_no = 0;
	try{
	    Gene *gene = NULL;
	    GeneFeature *pred_exon = NULL;

	    while (ifstrm) {
		string line;
		getline(ifstrm, line);
		line_no++;
		if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		    continue;
		vector<string> tokens = splitString(line);
		if(tokens.size() != 9)
		    throw ProjectError("wrong number of columns.\n");

		// determine transcript and geneid
		string attribute = tokens[8];
		stringstream ss(attribute);
		string geneid, transid, token;
		if(tokens[2] != "start_codon" && tokens[2] != "stop_codon" && tokens[2] != "CDS") // skip all features not in {CDS,start_codon,stop_codon}
		    continue;
		while(ss >> token){
		    if( token.compare("gene_id") == 0 ){
			if(ss >> token){
			    for(size_t i = 0; i < token.size(); ++i)
				if(token[i] != '"' && token[i] != ';') geneid += token[i];
			}
			
		    }
		    if( token.compare("transcript_id") == 0 ){
			if(ss >> token){
			    for(size_t i = 0; i < token.size(); ++i)
				if(token[i] != '"' && token[i] != ';') transid += token[i];
			}
		    }
		}
		if(geneid.empty())
		    throw ProjectError("missing gene_id.\n");
		if(transid.empty())
		    throw ProjectError("missing transcript_id.\n");

		if( gene && (gene->getGeneID() != geneid) && (gene->getTxID() == transid) )
		    throw ProjectError("Transcript " + transid + " in multiple genes "
				       + geneid + "," + gene->getGeneID() + ".\n");

		if(!gene || (gene && gene->getTxID() != transid)){ // new transcript starts

		    // insert new sequences both into the seqname AND seqID Hash
		    pair<map<string,int>::iterator,bool> seqid;
		    seqid = seqnames.insert(make_pair(tokens[0], seqnames.size()));
		    if(seqid.second)
			seqIDs.insert(make_pair(seqid.first->second, seqid.first->first));

		    Strand strand = getStrand(tokens[6]);
		    if(gene && gene->hasFeatures())
			genes.push_back(gene);
		    gene= new Gene(geneid,transid, seqid.first->second, strand, tokens[1]);
		    pred_exon = NULL;
		}
		if(tokens[2] == "start_codon"){
		    gene->setTLstart(atol(tokens[3].c_str()));
		}
		if(tokens[2] == "stop_codon"){
		    gene->setTLend(atol(tokens[4].c_str()));
		}
		if(tokens[2] == "CDS"){ // add exon to gene
		    FeatureType type = CDS;
		    long int start = atol(tokens[3].c_str());
		    long int end = atol(tokens[4].c_str());
		    double score = atof(tokens[5].c_str());
		    Strand strand = getStrand(tokens[6]);
		    int frame = getFrame(tokens[7]);
		    GeneFeature *exon = new GeneFeature(type, start,end,strand,frame,score);
		    exon->setGene(gene);

		    // add exon into gfHash
		    insertSeqInt(exon);

		    // add exon start/end into mappedPos
		    insertPos(exon->getSeqID(),start);
		    insertPos(exon->getSeqID(),end);

		    // create an intron if a predecessor exon exists
		    if(pred_exon){
			FeatureType type = intron;
			GeneFeature *intron = new GeneFeature(type, pred_exon->getEnd(),exon->getStart(),strand, -1, 0.0);
			gene->appendFeature(intron);
			intron->setGene(gene);
			// add to gfHash
			insertSeqInt(intron);
		    }
		    gene->appendFeature(exon);
		    pred_exon = exon;
		}
	    }
	    if(gene && gene->hasFeatures())
		genes.push_back(gene);
	} catch( ProjectError& err ){
	    cerr << "Error in " << gtffilename << " in line " << line_no << endl;
	    throw err;
	}
    }
    else
        throw ProjectError("Could not open input gtf file " + gtffilename + ".\n");
}


void Genome::writeGeneFeature(GeneFeature *gf, ofstream &of) const {
    of << getSeqName(gf->getSeqID()) << "\t";
    of << gf->getSource() << "\t";
    of << featureTypeIdentifiers[gf->getFeatureType()] << "\t";
    if(gf->isIntron()){
	of << gf->getStart()+1 << "\t";
	of << gf->getEnd()-1 << "\t";
    }
    else{
	of << gf->getStart() << "\t";
	of << gf->getEnd() << "\t";
    }
    of << gf->getScore() << "\t";
    of << strandIdentifiers[gf->getStrand()] << "\t";
    of << gf->writeFrame() << "\t";
    of << "transcript_id \"" << gf->getTxID() << "\"; ";
    of << "gene_id \"" << gf->getGeneID() << "\";" << endl;
}

void Genome::setTmpDir(string tmpdir){
    string dir = tmpdir + name + "/";
    createDir(dir);
    this->tmpdir=dir;
}

void Genome::printBed(){

    string filename = tmpdir + name + ".bed";
    ofstream of;
    of.open(filename.c_str());

    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
    for(it = mappedPos.begin(); it != mappedPos.end(); it++){
	SeqPosKey key(it->first);
	of << getSeqName(key.getSeqID()) << "\t";
	of << key.getPos()-1 << "\t"; // BED format uses zero-based, half-open coordinates, e.g. [1-25] -> [0,25)
	of << key.getPos() << "\t";  
	of << key.getKey() << endl;
    }
    of.close();
}

void Genome::readBed(Genome &other){

    string bedfilename = tmpdir + other.getName() + ".bed";
    ifstream ifstrm(bedfilename.c_str());
    int line_no = 0;
    if (ifstrm.is_open()){
	while (ifstrm) {
	    string line;
	    getline(ifstrm, line);
	    line_no++;
	    if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		continue;
	    vector<string> tokens = splitString(line);
	    if(tokens.size() != 4)
		throw ProjectError("wrong number of columns in line " + itoa(line_no) + "\n");
	    int seqID = other.getSeqID(tokens[0]);
	    if(seqID < 0) // no gene on that sequence
		continue;
	    long int start = atol(tokens[1].c_str()) + 1;
	    uint_fast64_t key = atoll(tokens[3].c_str()); 

	    // insert mapped position into mappedPos	    
	    SeqPosKey mapped(seqID, start);
	    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
	    it = mappedPos.find(key);                                                                                                 
	    if(it == mappedPos.end())
		throw ProjectError("internal error: unkown SeqIntKey" + itoa(key) + " in genome " + name + "\n");
	    it->second[other.getIdx()].push_back(mapped.getKey());
	}
    }
    else
	throw ProjectError("Could not open temporary bed file " + bedfilename + "\n");

}

// system call to halLiftover (external program) 
void Genome::liftOverTo(Genome &other, string halfile, string halLiftover_exec, string halParam){

    string cmd = halLiftover_exec + " " + halParam + halfile + " " + name + " " + tmpdir + name + ".bed " + other.getName() + " " + tmpdir + other.getName() + ".bed";

    // mu.lock(); // print 'cmd' mutually exclusive
    // cout << "executing " << cmd << endl;
    // mu.unlock();
    string ret = exec(cmd.c_str());
    if(!ret.empty())
	throw ProjectError(ret);
}

string Genome::getSeqName(int seqID) const {

    map<int,string>::const_iterator it = seqIDs.find(seqID);                                                                                                 
    if(it == seqIDs.end())
	throw ProjectError("internal error: unkown sequence ID seqID=" + itoa(seqID) + " in genome " + name + "\n");
    return it->second;

}

int Genome::getSeqID(string seqname) const{

    map<string, int>::const_iterator it = seqnames.find(seqname);                                                                                     
    if(it == seqnames.end())
	return -1;
    return it->second;

}

/*
 * mapping of homologous gene features
 * retrieves for each gene feature gf the mapped start and end positions in each genome j, and
 * assembles all combinations of mapped starts and ends to valid sequence intervals.
 * If an assembled seq interval is part of a gene in genome j,
 * it is appended to the list of homologs of gf
 */
void Genome::mapGeneFeatures(vector<Genome> &genomes){

    for(list<Gene*>::iterator git = genes.begin(); git != genes.end(); git++){	
	list<GeneFeature*> features = (*git)->getFeatureList();

	if(features.empty())
	    continue;
	GeneFeature* first = features.front();
	vector< list< uint_fast64_t > > mappedStarts = findMappedPos(first->getSeqID(), first->getStart());

	for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){
	    
	    // get mapped end positions
	    vector< list< uint_fast64_t > > mappedEnds = findMappedPos((*gfit)->getSeqID(), (*gfit)->getEnd());
	    
	    for(int j = 0; j < genomes.size(); j++){
		string other_name = genomes[j].getName();
		if (other_name == this->name)
		    continue;

		if(mappedStarts[j].empty() || mappedEnds[j].empty()) // at least one boundary is not mappable to genome j
		    continue;
		/*                                                                                       
		 * loop over all combinations of mapped start and end positions and
		 * assemble start position s and an end position e to a seq interval, if they are  
		 * on the same sequence. If s > e
		 * - reverse seq interval                                                                                                      
		 * - set strand to the reverse strand of gene feature gf
		 */
		for(std::list<uint_fast64_t>::iterator msit =mappedStarts[j].begin(); msit != mappedStarts[j].end(); msit++){
		    SeqPosKey s(*msit);
		    for(std::list<uint_fast64_t>::iterator meit =mappedEnds[j].begin(); meit != mappedEnds[j].end(); meit++){
			SeqPosKey e(*meit);
			if( s.getSeqID() == e.getSeqID() ){
			    uint_fast64_t start = s.getPos();
			    uint_fast64_t end = e.getPos();
			    Strand strand = (*gfit)->getStrand();
			    if( start > end ){ // reverse coordinates if they map on different strand
				uint_fast64_t tmp = start;
				start = end;
				end = tmp;
				strand = ((*gfit)->getStrand() == plusstrand)? minusstrand : plusstrand;
			    }
			    /*
			     * if the seq interval is part of a gene in genome j
			     * append it to the list of homologs of the gene feature
			     */
			    SeqIntKey seqInt(start,end-start+1, (*gfit)->getFeatureType());
			    list<GeneFeature*> mapped_features = genomes[j].findSeqInt(seqInt.getKey(), s.getSeqID(), strand);
			    for(list<GeneFeature*>::iterator mapped_gfit = mapped_features.begin(); mapped_gfit !=mapped_features.end();mapped_gfit++){
				if((*gfit)->isIntron()) // append Introns
				    (*gfit)->appendHomolog(*mapped_gfit,j);
				else if((*gfit)->sameFrame((*mapped_gfit)->getFrame()) && ((*gfit)->lenMod3() == (*mapped_gfit)->lenMod3())) // append CDS, if they are in the same frame
				    (*gfit)->appendHomolog(*mapped_gfit,j);				    		    
			    }
			}
		    }
		}
	    }
	    mappedStarts = mappedEnds;
	}
    }
}

void Genome::insertPos(int seqID, long int pos){

    SeqPosKey seqInt(seqID, pos);

    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
    it = mappedPos.find(seqInt.getKey());
    if (it == mappedPos.end()){
	vector< list< uint_fast64_t > > m(no_genomes);
	mappedPos.insert(pair<uint_fast64_t, vector< list< uint_fast64_t > > >(seqInt.getKey(), m));
    }
}


vector< list< uint_fast64_t > > Genome::findMappedPos(int seqID, long int pos){
    
    SeqPosKey seqInt(seqID, pos);
    
    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
    it = mappedPos.find(seqInt.getKey());                                                                                                 
    if(it == mappedPos.end())
	throw ProjectError("internal error: unkown SeqIntKey" + itoa(seqInt.getKey()) + " in genome " + name + "\n");
    return it->second;
}

 
list<GeneFeature*> Genome::findSeqInt(uint_fast64_t key, int seqID, Strand strand) {

    list<GeneFeature*> l;
    map< int, map< uint_fast64_t, list< GeneFeature*> > >::const_iterator ret;
    ret = gfHash.find(seqID);
    if(ret == gfHash.end())
	return l;
    map<uint_fast64_t, list<GeneFeature*> >::const_iterator it;
    it = ret->second.find(key);
    if (it == ret->second.end())
	return l;
    for(list<GeneFeature*>::const_iterator lit=it->second.begin(); lit!=it->second.end(); lit++){
	if((*lit)->sameStrand(strand))
	    l.push_back(*lit);
    }
    return l;
}

void Genome::insertSeqInt(GeneFeature* gf){ // insert exons/introns into gfHash

    SeqIntKey seqInt(gf->getStart(), gf->getLen(), gf->getFeatureType());
        
    pair< map< int, map< uint_fast64_t, list< GeneFeature*> > >::iterator, bool> ret;
    map< uint_fast64_t, list< GeneFeature* > >::iterator it;
    map< uint_fast64_t, list< GeneFeature* > > m;
    
    ret = gfHash.insert(pair< int, map< uint_fast64_t, list< GeneFeature* > > >(gf->getSeqID(),m));
			
    it =ret.first->second.find(seqInt.getKey());
    if (it == ret.first->second.end()){ // insert new exon/intron       
	list<GeneFeature*> l;
	l.push_back(gf);
	ret.first->second.insert(pair<uint_fast64_t,list<GeneFeature*> >(seqInt.getKey(), l));
    } else {                // append exon/intron to existing list                                 
	it->second.push_back(gf);
    }
}

void Genome::insertSeqInt(GeneFeature* gf, int seqID){ // insert exons/introns into gfHash

    SeqIntKey seqInt(gf->getStart(), gf->getLen(), gf->getFeatureType());
        
    pair< map< int, map< uint_fast64_t, list< GeneFeature*> > >::iterator, bool> ret;
    map< uint_fast64_t, list< GeneFeature* > >::iterator it;
    map< uint_fast64_t, list< GeneFeature* > > m;
    
    ret = gfHash.insert(pair< int, map< uint_fast64_t, list< GeneFeature* > > >(seqID,m));
			
    it =ret.first->second.find(seqInt.getKey());
    if (it == ret.first->second.end()){ // insert new exon/intron       
	list<GeneFeature*> l;
	l.push_back(gf);
	ret.first->second.insert(pair<uint_fast64_t,list<GeneFeature*> >(seqInt.getKey(), l));
    } else {                // append exon/intron to existing list                                                                                                 
	it->second.push_back(gf);
    }
}

void Genome::printGFF(string outdir, vector<Genome> &genomes){
    string filename = outdir + name + ".gtf";
    ofstream of;

    of.open(filename.c_str());

    // print header
    for(int i=0; i < genomes.size(); i++){
	of << "# " << genomes[i].getIdx() << "\t" << genomes[i].getName() << endl;
    }
    of << "###" << endl;

    // summary stats (over all genes)
    vector<int> total_mappedStatsG(no_genomes);   // number of transcripts with exact homologs in at least k other genomes
    vector<int> total_mappedStatsE(no_genomes);   // number of CDS with exact homologs in at least k other genomes
    vector<int> total_mappedStatsI(no_genomes);   // number of Intr ...
    vector<int> total_extrinStatsE(no_genomes+1); // number of CDS supported by evidence in at least k other genomes
    vector<int> total_extrinStatsI(no_genomes+1); // number of Intr ...

    
    for(list<Gene*>::iterator git = genes.begin(); git != genes.end(); git++){	
	list<GeneFeature*> features = (*git)->getFeatureList();

	of << "# start transcript " << (*git)->getTxID() << endl;

	writeGene((*git), of);

	vector<map<string,GeneInfo> > ginfo(no_genomes);
	pair<map<string,GeneInfo>::iterator, bool> ret;

	// stats for a single gene
	vector<int> mappedStatsE(no_genomes);   // number of CDS with exact homologs in at least k other genomes
	vector<int> mappedStatsI(no_genomes);   // number of Intr ...
	vector<int> extrinStatsE(no_genomes+1); // number of CDS supported by evidence in at least k other genomes
	vector<int> extrinStatsI(no_genomes+1); // number of Intr ...

	of << "#\n# gene feature level:\n#" << endl;
	for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){

	    printDetailed(*gfit, of);
	    
	    bool hasEvidence = false;
	    bool hasMatch = false;
	    int numMatches = 0;
	    int numEvidence = 0;
	    if((*gfit)->hasEvidence())
		numEvidence++;

	    list<pair<int,GeneFeature*> >homologs = (*gfit)->getHomologs();
	      
	    int pred_idx = 0;

	    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){

		int idx = hit->first;
		GeneFeature* h = hit->second;

		if(pred_idx < idx && hasMatch){ // new genome
		    numMatches++;
		    hasMatch = false; // reset flags
		}
		if(pred_idx < idx && hasEvidence){
		    numEvidence++;
		    hasEvidence = false;
		}

		if(h->hasEvidence()){
		    hasEvidence = true;
		}
		if(h->isPartofGene()){
		    hasMatch = true;
		    GeneInfo gi(h->getGene(),h->isExon(),h->isIntron(),false);
		    ret = ginfo[idx].insert(pair<string,GeneInfo>(h->getTxID(),gi));
		    if(!ret.second){ // transid already in map, only update
			if(h->isExon())
			    ret.first->second.numMatchingEs++;
			else if(h->isIntron())
			    ret.first->second.numMatchingIs++;
		    }
		}
		pred_idx = idx;
	    }
	    if(hasMatch)
		numMatches++;
	    if(hasEvidence)
		numEvidence++;

	    if((*gfit)->isExon()){
		total_mappedStatsE[numMatches]++;
		total_extrinStatsE[numEvidence]++;
		mappedStatsE[numMatches]++;
		extrinStatsE[numEvidence]++;

	    }
	    else{
		total_mappedStatsI[numMatches]++;
		total_extrinStatsI[numEvidence]++;
		mappedStatsI[numMatches]++;
		extrinStatsI[numEvidence]++;
	    }
			   
	}
	of << "#" << endl;
	of << "# number/% of features with exact homologs in at least k other genomes:" << endl;
	of << "#" << endl;
	of << "#---------------------------------------------------------------------------------" << endl;       
	of << "#   k           CDS             Intr             Both" << endl;
	of << "#---------------------------------------------------------------------------------" << endl;
	of << "#" << endl;
	writePictograph(mappedStatsE,mappedStatsI,of);
	of << "#" << endl;
	of << "# number/% of features supported by extrinsic evidence in at least k genomes :" << endl;
	of << "#" << endl;
	of << "#---------------------------------------------------------------------------------" << endl;       
	of << "#   k           CDS             Intr             Both" << endl;
	of << "#---------------------------------------------------------------------------------" << endl;       
	of << "#" << endl;
	writePictograph(extrinStatsE,extrinStatsI,of);
	of << "#\n# transcript level:\n#" << endl;
	
	int numHomologs = 0;

	for(int i=0; i < ginfo.size(); i++){

	    bool hasHomolog = false;
	    for(map<string,GeneInfo>::iterator giit = ginfo[i].begin(); giit != ginfo[i].end(); giit++){
		
		int numE = 1 + (giit->second.gene->numGFs()/2);
		int numI = (giit->second.gene->numGFs()/2);
		
		of << "#";
		of << setw(4) << right << i << " ";
		of << "gene_id="<< setw(12) << left << giit->second.gene->getGeneID();
		of << "tx_id="<< setw(12) << left << giit->first;
		of << setw(3) << right << giit->second.numMatchingEs;
 		of << "/"<< numE << " CDS";
		of <<setw(3) << right << giit->second.numMatchingIs;
		of <<"/"<< numI << " Intr" << endl;
		if(numE == giit->second.numMatchingEs && numI == giit->second.numMatchingIs && (numE+numI) == (*git)->numGFs()){
		    hasHomolog = true;
		    (*git)->appendHomolog(giit->second.gene,i);		    
		}
	    }
	    if(hasHomolog)
		numHomologs++;	
	}
 	total_mappedStatsG[numHomologs]++;

	of << "#" << endl;
	of << "# transcript has an exact homolog in " << numHomologs << " other genomes." <<endl;

	of << "#\n# end transcript " << (*git)->getTxID() << endl;
	of << "###" << endl;

    }
    // print cumulative statistics
    of << "# summary statistic over all genes" << endl;
    of << "#" << endl;
    of << "# gene feature level:" << endl;
    of << "#" << endl;
    of << "# number/% of features with exact homologs in at least k other genomes:" << endl;
    of << "#" << endl;
    of << "#---------------------------------------------------------------------------------" << endl;       
    of << "#   k           CDS             Intr             Both" << endl;
    of << "#---------------------------------------------------------------------------------" << endl;
    of << "#" << endl;
    writePictograph(total_mappedStatsE,total_mappedStatsI,of);
    of << "#" << endl;
    of << "# number/% of features supported by extrinsic evidence in at least k genomes :" << endl;
    of << "#" << endl;
    of << "#---------------------------------------------------------------------------------" << endl;       
    of << "#   k           CDS             Intr             Both" << endl;
    of << "#---------------------------------------------------------------------------------" << endl;       
    of << "#" << endl;
    writePictograph(total_extrinStatsE,total_extrinStatsI,of);
    of << "#" << endl;
    of << "# transcript level:" << endl;
    of << "#" << endl;
    of << "# number/% of transcripts with exact homologs in at least k other genomes:" << endl;
    of << "#" << endl;
    of << "#-----------------------------------------------" << endl;       
    of << "#   k            tx" << endl;
    of << "#-----------------------------------------------" << endl;
    writePictograph(total_mappedStatsG,of);
    of.close();
}

void Genome::printDetailed(GeneFeature *g, std::ofstream &of) const{

    list<pair<int,GeneFeature*> >homologs = g->getHomologs();

    int pred_idx = 0;
    bool hasMatch = false;
    string evidence = "";
    string onlyHint = "*";

    
    of << "# ";
    if(g->isExon())
	of << "CDS  ";
    else
	of << "Intr ";
    if(g->hasEvidence())
	of << g->getEvidence() <<" (";
    else
	of <<"  (";
    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){
	
	int idx = hit->first;
	GeneFeature* h = hit->second;
	
	if(pred_idx < idx && hasMatch){ // next genome
	    of << pred_idx << evidence << onlyHint << " ";
	    hasMatch = false; // reset all flags
	    evidence = "";
	    onlyHint = "*";
	}
	hasMatch = true;
	if(evidence.empty() && h->hasEvidence())
	    evidence = h->getEvidence();
	if(h->isPartofGene()){
	    onlyHint = "";
	}
	pred_idx = idx;
    }   
    if(hasMatch)      
	of << pred_idx << evidence << onlyHint << " ";
    of << ")" << endl;

}

void Genome::writeGene(Gene *g, std::ofstream &of) const{

    // gene line
    /*of << getSeqName(g->getSeqID()) << "\t";
    of << g->getSource() << "\tgene\t";
    of << g->getStart() << "\t";
    of << g->getEnd() << "\t.\t";
    of << strandIdentifiers[g->getStrand()] << "\t.\t";
    of << g->getGeneID() << endl;
    */
    // transcript line
    of << getSeqName(g->getSeqID()) << "\t";
    of << g->getSource() << "\ttranscript\t";
    of << g->getStart() << "\t";
    of << g->getEnd() << "\t.\t";
    of << strandIdentifiers[g->getStrand()] << "\t.\t";
    of << g->getTxID() << endl;
    
    if(g->getStrand() == plusstrand)
	writeTLStart(g, of);
    else
	writeTLEnd(g,of);

    // print exons and introns
    list<GeneFeature*> features = g->getFeatureList();
    for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){
	writeGeneFeature(*gfit, of);
    }
    if(g->getStrand() == plusstrand)
	writeTLEnd(g, of);
    else
	writeTLStart(g,of);
}

void Genome::writeTLStart(Gene *g, std::ofstream &of) const{
  if(g->getTLstart() >= 0){
	of << getSeqName(g->getSeqID()) << "\t";
	of << g->getSource() << "\tstart_codon\t";
	of << g->getTLstart() << "\t";
	of << g->getTLstart()+2 << "\t.\t";
	of << strandIdentifiers[g->getStrand()] << "\t.\t";
	of << "transcript_id \"" << g->getTxID() << "\"; ";
	of << "gene_id \"" << g->getGeneID() << "\";" << endl;
    }
 
}
void Genome::writeTLEnd(Gene *g, std::ofstream &of) const{
    if( g->getTLend() >= 0){
	of << getSeqName(g->getSeqID()) << "\t";
	of << g->getSource() << "\tstop_codon\t";
	of << g->getTLend()-2 << "\t";
	of << g->getTLend() << "\t.\t";
	of << strandIdentifiers[g->getStrand()] << "\t.\t";
	of << "transcript_id \"" << g->getTxID() << "\"; ";
	of << "gene_id \"" << g->getGeneID() << "\";" << endl;
    }
}

#ifdef BOOST
void printHomGeneList(string outfile, vector<Genome> &genomes){
    ofstream of;
    of.open(outfile.c_str());
    if(of.is_open()){

	// print header
	for(int i=0; i < genomes.size(); i++){
	    of << "# " << genomes[i].getIdx() << "\t" << genomes[i].getName() << endl;
	}

	// print list of homologous Tx Ids

	map<string,int> txs; // keys: "(idx,txid)", values: vertex indices in boost Graph
	pair<map<string,int>::iterator,bool> txit;

        Graph G;

	for(int i = 0; i < genomes.size(); i++){
	    for(list<Gene*>::iterator git = genomes[i].genes.begin(); git != genomes[i].genes.end(); git++){
		list<pair<int,Gene*> >homologs = (*git)->getHomologs();
   
		string u_name = itoa(i) + "," + (*git)->getTxID(); 
		txit = txs.insert(make_pair(u_name, txs.size()));
		int u = txit.first->second;
		if(homologs.empty()){
		    boost::add_edge(u,u, G);
		    G[u].name = u_name;
		}
		for(list<pair<int,Gene*> >::iterator hgit = homologs.begin(); hgit != homologs.end(); hgit++){
		    int idx = hgit->first;
		    Gene* g = hgit->second;
		    string v_name = itoa(idx) + "," + g->getTxID(); 
		    txit = txs.insert(make_pair(v_name, txs.size()));
		    int v = txit.first->second;
		    boost::add_edge(u,v,G);
		    G[u].name=u_name;
		    G[v].name=v_name;
		}
	    }
	}
	std::vector<int> component(num_vertices(G));
        connected_components(G, &component[0]);

        int c = 0;
        for (std::vector<int>::size_type i = 0; i != component.size(); ++i){
	    if(component[i] != c)
		of << endl;
	    of << "(" << G[i].name <<") ";
	    c = component[i];
	}
        of << endl;
	of.close();
    }
    else{
	cerr << "Could not open output file " << outfile << endl;
    }
}
#endif
