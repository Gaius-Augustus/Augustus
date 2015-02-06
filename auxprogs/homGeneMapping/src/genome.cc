
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

using namespace std;

int Genome::no_genomes=0;
 

void Genome::parseExtrinsicGFF(string gfffilename){
    ifstream ifstrm(gfffilename.c_str());
    if (ifstrm.is_open()){
	try{
	    int line_no = 0;
	    while (ifstrm) {
		string line;
		getline(ifstrm, line);
		line_no++;
		if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		    continue;
		vector<string> tokens = splitString(line);
		if(tokens.size() != 9)
		    throw ProjectError("wrong number of columns in line " + itoa(line_no) + "\n");

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
		    throw ProjectError("Error in hint line " + itoa(line_no) + "\nNo source specified (e.g. by source=M in the last column)\n");
		}
		/*
		 * find all gene features that are supported
		 * by that hint and update their extrinsic source
		 */
		if(tokens[2] == "CDS"){
		    int frame = getFrame(tokens[7]);
		    FeatureType type = CDS;
		    SeqIntKey seqInt(start, end-start+1,type, strand);
		    list<GeneFeature*> le = findSeqInt(seqInt.getKey(), seqid);
		    for(list<GeneFeature*>::iterator it=le.begin(); it!=le.end(); it++){
			if(frame == (*it)->getFrame())
			    (*it)->setEvidence(esource); // update extrinsic source
		    }
		}
		else if(tokens[2] == "intron"){
		    start--;
		    end++;
		    FeatureType type = intron;
		    if(strand == unknown || strand == plusstrand){
			SeqIntKey seqInt(start, end-start+1, type, plusstrand);
			list<GeneFeature*> li = findSeqInt(seqInt.getKey(), seqid);
			for(list<GeneFeature*>::iterator it=li.begin(); it!=li.end(); it++){
			    (*it)->setEvidence(esource);	
			}
		    }
		    if(strand == unknown || strand == minusstrand){
			SeqIntKey seqInt(start, end-start+1, type, minusstrand);
			list<GeneFeature*> li = findSeqInt(seqInt.getKey(), seqid);
			for(list<GeneFeature*>::iterator it=li.begin(); it!=li.end(); it++){
			    (*it)->setEvidence(esource);	
			}
		    }
		}
	    }	    
	    
	} catch( ProjectError& err ){
	    cerr << "ERROR\n\t" << err.getMessage( ) << endl;
	    exit(1);
	}
    }
    else
        cerr<<"Could not open input gff file "<<gfffilename << endl;
}

void Genome::parseGTF(string gtffilename){
    ifstream ifstrm(gtffilename.c_str());
    if (ifstrm.is_open()){
	//	cout << gtffilename << " is open." << endl;
	try{
	    Gene *gene = NULL;
	    GeneFeature *pred_exon = NULL;

	    int line_no = 0;
	    while (ifstrm) {
		string line;
		getline(ifstrm, line);
		line_no++;
		if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		    continue;
		vector<string> tokens = splitString(line);
		if(tokens.size() != 9)
		    throw ProjectError("wrong number of columns in line " + itoa(line_no) + "\n");

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
		    throw ProjectError("missing gene_id in line  " + itoa(line_no) + "\n");
		if(transid.empty())
		    throw ProjectError("missing transcript_id in line  " + itoa(line_no) + "\n");

		if(!gene || (gene && gene->getGeneID() != geneid)){ // new gene starts

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
		    int frame = getFrame(tokens[7]);
		    GeneFeature *exon = new GeneFeature(type, start,end,frame,score);
		    exon->setGene(gene);

		    // add exon into gfHash
		    insertSeqInt(exon);

		    // add exon start/end into mappedPos
		    insertPos(exon->getSeqID(),start,exon->getStrand());
		    insertPos(exon->getSeqID(),end,exon->getStrand());

		    // create an intron if a predecessor exon exists
		    if(pred_exon){
			FeatureType type = intron;
			GeneFeature *intron = new GeneFeature(type, pred_exon->getEnd(),exon->getStart(), -1, 0.0);
			gene->appendFeature(intron);
			intron->setGene(gene);
			// add to gfHash
			insertSeqInt(intron);
		    }
		    gene->appendFeature(exon);
		    pred_exon = exon;
		}
	    }
	} catch( ProjectError& err ){
	    cerr << "ERROR\n\t" << err.getMessage( ) << endl;
	    exit(1);
	}
    }
    else
        cerr<<"Could not open input gtf file "<<gtffilename << endl;

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
	of << key.getKey() <<"\t0\t";
	of << strandIdentifiers[key.getStrand()] << endl;
    }
    of.close();
}

void Genome::readBed(Genome &other){

    string bedfilename = tmpdir + other.getName() + ".bed";
    ifstream ifstrm(bedfilename.c_str());
    int line_no = 0;
    if (ifstrm.is_open()){
	//	cout << gfffilename << " is open." << endl;
	while (ifstrm) {
	    string line;
	    getline(ifstrm, line);
	    //cout << line << endl;
	    line_no++;
	    if(line.empty() || line[0] == '#') // skip empty lines and comment lines
		continue;
	    vector<string> tokens = splitString(line);
	    if(tokens.size() != 6)
		throw ProjectError("wrong number of columns in line " + itoa(line_no) + "\n");
	    int seqID = other.getSeqID(tokens[0]);
	    if(seqID < 0) // no gene on that sequence
		continue;
	    long int start = atol(tokens[1].c_str()) + 1;
	    uint_fast64_t key = atoll(tokens[3].c_str()); 
	    Strand strand=getStrand(tokens[5]);

	    // insert mapped position into mappedPos	    
	    SeqPosKey mapped(seqID, start, strand);
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
    cout << "executing " << cmd << endl;
    string s = exec(cmd.c_str());
    if(!s.empty()){
	cerr << s << endl;
	exit(1);
    }
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
	vector< list< uint_fast64_t > > mappedStarts = findMappedPos(first->getSeqID(), first->getStart(), first->getStrand());

	for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){
	    
	    // get mapped end positions
	    vector< list< uint_fast64_t > > mappedEnds = findMappedPos((*gfit)->getSeqID(), (*gfit)->getEnd(), (*gfit)->getStrand());
	    
	    for(int j = 0; j < genomes.size(); j++){
		string other_name = genomes[j].getName();
		if (other_name == this->name)
		    continue;

		if(mappedStarts[j].empty() || mappedEnds[j].empty()) // at least one boundary is not mappable to genome j
		    continue;
		/*                                                                                       
		 * loop over all combinations of mapped start and end positions and                                                                          
		 * assemble start position s and an end position e to a seq interval, if they are                                                           
		 * - on the same sequence                                                                                                                                      
		 * - on the same strand                                                                                                                                         
		 * - and s < e                                                                                                                                    
		 */
		for(std::list<uint_fast64_t>::iterator msit =mappedStarts[j].begin(); msit != mappedStarts[j].end(); msit++){
		    SeqPosKey s(*msit);
		    for(std::list<uint_fast64_t>::iterator meit =mappedEnds[j].begin(); meit != mappedEnds[j].end(); meit++){
			SeqPosKey e(*meit);
			if((s.getSeqID() == e.getSeqID()) && (s.getStrand() == e.getStrand()) ){
			    uint_fast64_t start = s.getPos();
			    uint_fast64_t end = e.getPos();
			    if( (*gfit)->getStrand() != s.getStrand() ){ // reverse coordinates if they map on different strand
				uint_fast64_t tmp = start;
				start = end;
				end = tmp;
			    }
			    if(start < end){
				/*
				 * if the seq interval is part of a gene in genome j
				 * append it to the list of homologs of the gene feature
				 */
				SeqIntKey seqInt(start,end-start+1, (*gfit)->getFeatureType(), s.getStrand());
				list<GeneFeature*> mapped_features = genomes[j].findSeqInt(seqInt.getKey(), s.getSeqID());
				for(list<GeneFeature*>::iterator mapped_gfit = mapped_features.begin(); mapped_gfit !=mapped_features.end();mapped_gfit++)
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

void Genome::insertPos(int seqID, long int pos, Strand strand){

    SeqPosKey seqInt(seqID, pos, strand);

    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
    it = mappedPos.find(seqInt.getKey());
    if (it == mappedPos.end()){
	vector< list< uint_fast64_t > > m(no_genomes);
	mappedPos.insert(pair<uint_fast64_t, vector< list< uint_fast64_t > > >(seqInt.getKey(), m));
    }
}


vector< list< uint_fast64_t > > Genome::findMappedPos(int seqID, long int pos, Strand strand){
    
    SeqPosKey seqInt(seqID, pos, strand);
    
    map<uint_fast64_t, vector< list< uint_fast64_t > > >::iterator it;
    it = mappedPos.find(seqInt.getKey());                                                                                                 
    if(it == mappedPos.end())
	throw ProjectError("internal error: unkown SeqIntKey" + itoa(seqInt.getKey()) + " in genome " + name + "\n");
    return it->second;
}

 
list<GeneFeature*> Genome::findSeqInt(uint_fast64_t key, int seqID) {

    list<GeneFeature*> l;
    map< int, map< uint_fast64_t, list< GeneFeature*> > >::const_iterator ret;
    ret = gfHash.find(seqID);
    if(ret == gfHash.end())
	return l;
    map<uint_fast64_t, list<GeneFeature*> >::const_iterator it;
    it = ret->second.find(key);
    if (it == ret->second.end())
	return l;
    return it->second;
}

void Genome::insertSeqInt(GeneFeature* gf){ // insert exons/introns into gfHash
    
    SeqIntKey seqInt(gf->getStart(), gf->getLen(), gf->getFeatureType(), gf->getStrand());
    
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

void Genome::printGFF(string outdir, vector<Genome> &genomes){
    string filename = outdir + name + ".gtf";
    ofstream of;

    of.open(filename.c_str());

    // print header
    for(int i=0; i < genomes.size(); i++){
	of << "# " << genomes[i].getIdx() << "\t" << genomes[i].getName() << endl;
    }
    of << "###" << endl;
    
    for(list<Gene*>::iterator git = genes.begin(); git != genes.end(); git++){	
	list<GeneFeature*> features = (*git)->getFeatureList();

	int numSupportedE = 0;
	int numSupportedI = 0;
	int numMatchingE = 0;
	int numMatchingI = 0;

	of << "# start gene " << (*git)->getGeneID() << endl;

	writeGene((*git), of);

	vector<map<string,GeneInfo> > ginfo(no_genomes);
	pair<map<string,GeneInfo>::iterator, bool> ret;

	
	for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){

	    printDetailed(*gfit, of);
	    
	    bool hasEvidence = (*gfit)->hasEvidence();

	    list<pair<int,GeneFeature*> >homologs = (*gfit)->getHomologs();
	      
	    if(!homologs.empty()){
		if((*gfit)->isExon())
		    numMatchingE++;
		if((*gfit)->isIntron())
		    numMatchingI++;
	    }
	    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){

		bool hasFrameShift = false;

		int idx = hit->first;
		GeneFeature* h = hit->second;
		if(h->hasEvidence()){
		    hasEvidence = true;
		}
		if(h->isExon()){ // check for frame shifts
		    if ( ((*gfit)->getFrame() != h->getFrame()) || ((*gfit)->lenMod3() != h->lenMod3()) )
			hasFrameShift = true;
		}

		GeneInfo gi(h->getGene(),h->isExon(),h->isIntron(),hasFrameShift);
		ret = ginfo[idx].insert(pair<string,GeneInfo>(h->getGeneID(),gi));
		if(!ret.second){ // geneid already in map, only update
		    if(h->isExon())
 			ret.first->second.numMatchingEs++;
		    else if(h->isIntron())
			ret.first->second.numMatchingIs++;
		    if(hasFrameShift)
			ret.first->second.frameshift = true;
		}
	    }
	    if(hasEvidence){
		if((*gfit)->isExon())
		    numSupportedE++;
		if((*gfit)->isIntron())
		    numSupportedI++;
	    }
	}
	int numHomologs = 0;

	for(int i=0; i < ginfo.size(); i++){
	    for(map<string,GeneInfo>::iterator giit = ginfo[i].begin(); giit != ginfo[i].end(); giit++){
		
		int numE = 1 + (giit->second.gene->numGFs()/2);
		int numI = (giit->second.gene->numGFs()/2);
		
		of <<"# "<< i << "\t" << giit->first << "\t" << giit->second.numMatchingEs <<"/"<< numE << " exons " << giit->second.numMatchingIs <<"/"<< numI << " introns ";
		if(giit->second.frameshift)
		    of << "(frameshift)";
		of << endl;
		if(numE == giit->second.numMatchingEs && numI == giit->second.numMatchingIs && (numE+numI) == (*git)->numGFs())
		    numHomologs++;
	    }
		
	}
	if(numHomologs>0)
	    of << "# number of exact homologs: " << numHomologs << endl;

	int numE = 1 + ((*git)->numGFs()/2);
	int numI = ((*git)->numGFs()/2);


	double percSupported = ((double) (numSupportedI + numSupportedE) )/ (*git)->numGFs();
	double percMatching = ((double) (numMatchingI + numMatchingE) )/ (*git)->numGFs();
	
	of << "# % of transcript supported by hints (any source any genome): " << setprecision(3) 
	   << 100.0 * percSupported <<  endl;
	of << "# CDS exons: "<< numSupportedE << "/" << numE << endl;
	of << "# CDS introns: "<< numSupportedI << "/" << numI << endl;

	of << "# % of transcript inconsistent with other genomes: " << setprecision(3) 
	   << (100 * (1 - percMatching)) <<  endl;
	of << "# CDS exons: "<< (numE - numMatchingE) << "/" << numE << endl;
	of << "# CDS introns: "<< (numI - numMatchingI) << "/" << numI << endl;

	of << "# end gene " << (*git)->getGeneID() << endl;
	of << "###" << endl;

    }
    of.close();

}

void Genome::printDetailed(GeneFeature *g, std::ofstream &of) const{

    list<pair<int,GeneFeature*> >homologs = g->getHomologs();
    
    int pred_idx = 0;
    bool hasMatch = false;
    string evidence = "";

    
    of << "# "<<g->getEvidence() <<" (";
    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){
	
	int idx = hit->first;
	GeneFeature* h = hit->second;
	
	if(pred_idx < idx && hasMatch){ // next genome
	    of << pred_idx << evidence << " ";
	    hasMatch = false; // reset all flags
	    evidence = "";
	}
	if(g->isExon()){ // check for frame shifts
	    if ( !(g->getFrame() != h->getFrame()) || (g->lenMod3() != h->lenMod3()) )
		hasMatch = true;
	}
	else
	    hasMatch = true;
	if(evidence.empty() && h->hasEvidence())
	    evidence = h->getEvidence();
	pred_idx = idx;
    }   
    if(hasMatch){ // next genome
	of << pred_idx <<evidence << " ";
    }
    of << ")" << endl;

}

void Genome::writeGene(Gene *g, std::ofstream &of) const{

    // gene line
    of << getSeqName(g->getSeqID()) << "\t";
    of << g->getSource() << "\tgene\t";
    of << g->getStart() << "\t";
    of << g->getEnd() << "\t.\t";
    of << strandIdentifiers[g->getStrand()] << "\t.\t";
    of << g->getGeneID() << endl;

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


