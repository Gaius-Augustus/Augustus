/**********************************************************************
 * file:    genome.cc
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

void Genome::parse(string genefile, string hintsfile, string dbfile){
    parseGTF(genefile);
    if(!hintsfile.empty()) // read hints file if specified
	parseExtrinsicGFF(hintsfile);
#ifdef M_SQLITE
    if(!dbfile.empty()){
	SQLiteDB db = SQLiteDB(dbfile.c_str());
	getDbHints(db);
    }
#endif
    printBed(); // print sequence coordinates, that need to be mapped to the other genomes, to file
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

		string seqname = tokens[0];
		long int start = atol(tokens[3].c_str());
		long int end = atol(tokens[4].c_str());
		Strand strand = getStrand(tokens[6]);
		int frame = getFrame(tokens[7]);

		/*
		 * find source of extrinsic info, specified in gff as: source=X or src=X
		 */
		const char *spos;
		string esource;
		int mult = 1;

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
		spos = strstr(tokens[8].c_str(), "mult=");
		if (spos)
		    spos += 5;
		if (spos)
		    mult = atoi(spos);

		string type = tokens[2];
		if(type == "CDS" || type == "intron" || type == "exon") // currently only CDS , exon and intron hints
		    insertHint(seqname,start,end,strand,esource,mult,frame,type);

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

		FeatureType type = getType(tokens[2]);
		if(type == unkown) // skip all features not in {CDS,start_codon,stop_codon,exon,UTR,intron}
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
		    if(gene && gene->hasFeatures()){
			gene->includeStopInCDS();
			gene->insertMissingGFs();
			genes.push_back(gene);
			insertSeqInts(gene);
		    }
		    gene= new Gene(geneid,transid, seqid.first->second, strand, tokens[1]);
		}
		if(type == start){
		    gene->setTLstart(atol(tokens[3].c_str()));
		}
		if(type == stop){
		    gene->setTLend(atol(tokens[3].c_str()));
		}
		if(type == CDS || type == UTR || type == exon || type == intron){ // add exon to gene
		    long int start = atol(tokens[3].c_str());
		    long int end = atol(tokens[4].c_str());
		    if(type == intron){
			start-=1;
			end+=1;
		    }
		    double score = atof(tokens[5].c_str());
		    Strand strand = getStrand(tokens[6]);
		    int frame = getFrame(tokens[7]);
		    GeneFeature *gf = new GeneFeature(type, start,end,strand,frame,score);
		    gene->appendFeature(gf);
		    gf->setGene(gene);
		}
	    }
	    if(gene && gene->hasFeatures()){
		gene->includeStopInCDS();
		gene->insertMissingGFs();
		genes.push_back(gene);
		insertSeqInts(gene);
	    }
	} catch( ProjectError& err ){
	    cerr << "Error in " << gtffilename << " in line " << line_no << endl;
	    throw err;
	}
    }
    else
        throw ProjectError("Could not open input gtf file " + gtffilename + ".\n");
}


void Genome::writeGeneFeature(GeneFeature *gf, ofstream &of, bool unmapped) const {
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
    of << "gene_id \"" << gf->getGeneID() << "\"; " ;
    if(!gf->isType(stop) && !gf->isType(start)){
	of << "hgm_info \"";
	print_hgm_info(gf,of);
	of  << "\"; ";
	if(unmapped){
	    of << "hgm_mapped \"";
	    print_hgm_unaligned(gf,of);
	    of  << "\"; ";
	}
    }
    of << endl;
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
		throw ProjectError("internal error: unknown SeqIntKey" + itoa(key) + " in genome " + name + "\n");
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
	throw ProjectError("internal error: unknown sequence ID seqID=" + itoa(seqID) + " in genome " + name + "\n");
    return it->second;

}

int Genome::getSeqID(string seqname) const{

    map<string, int>::const_iterator it = seqnames.find(seqname);                                             
    if(it == seqnames.end())
	return -1;
    return it->second;

}

void Genome::write_hgm_gff(vector<Genome> &genomes, string outdir, bool details, bool unmapped){
    for(int j = 0; j < genomes.size(); j++){
	if(idx != j)
	    readBed(genomes[j]); // reading in bed files with mapped coordinates
    }
    mapGeneFeatures(genomes, outdir, details, unmapped);
}

/*
 * mapping of homologous gene features
 * retrieves for each gene feature gf the mapped start and end positions in each genome j, and
 * assembles all combinations of mapped starts and ends to valid sequence intervals.
 * If an assembled seq interval is part of a gene in genome j,
 * it is appended to the list of homologs of gf
 */
void Genome::mapGeneFeatures(vector<Genome> &genomes, string outdir, bool detailed, bool unmapped){
    string filename = outdir + name + ".gtf";
    ofstream of;
    
    of.open(filename.c_str());

    // print header
    for(int i=0; i < genomes.size(); i++){
	of << "# " << genomes[i].getIdx() << "\t" << genomes[i].getName() << endl;
    }
    of << "###" << endl;

    for(list<Gene*>::iterator git = genes.begin(); git != genes.end(); git++){	
	list<GeneFeature*> features = (*git)->features;

	if(features.empty())
	    continue;
	writeTxLine(*git,of);
	GeneInfoCollection gfc(no_genomes);
	for(list<GeneFeature*>::iterator gfit = features.begin(); gfit != features.end(); gfit++){
	    
	    if(!(*gfit)->isMapped()){
		writeGeneFeature(*gfit,of);
		continue;
	    }
	    // get mapped start and end positions
	    vector< list< uint_fast64_t > > mappedStarts = findMappedPos((*gfit)->getSeqID(), (*gfit)->getStart());
	    vector< list< uint_fast64_t > > mappedEnds = findMappedPos((*gfit)->getSeqID(), (*gfit)->getEnd());
	    
	    for(int j = 0; j < genomes.size(); j++){
		string other_name = genomes[j].getName();
		if (other_name == this->name)
		    continue;

		if(mappedStarts[j].empty() || mappedEnds[j].empty()){ // at least one boundary is not mappable to genome j
		    (*gfit)->appendHomolog(NULL,j);
		    continue;
		}
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
				if((*gfit)->isIntron() || (*gfit)->isExon()) // append Introns and Exons
				    (*gfit)->appendHomolog(*mapped_gfit,j);
				else if((*gfit)->sameFrame((*mapped_gfit)->getFrame()) && ((*gfit)->lenMod3() == (*mapped_gfit)->lenMod3())) // append CDS, if they are in the same frame
				    (*gfit)->appendHomolog(*mapped_gfit,j);
			    }
			}
		    }
		}
	    }
	    writeGeneFeature(*gfit, of, unmapped);
	    gfc.createCollection(*gfit);
	    (*gfit)->homologs.clear();
	}
	if(detailed)
	    gfc.printDetailedStats(*git, of);
    }
    mappedPos.clear();
    of.close();
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
	throw ProjectError("internal error: unknown SeqIntKey" + itoa(seqInt.getKey()) + " in genome " + name + "\n");
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

void Genome::insertSeqInts(Gene *g){
    
    for(list<GeneFeature*>::const_iterator it = g->features.begin(); it != g->features.end(); it++){
	if(!(*it)->isMapped())
	    continue;
	// add to gfHash                                                                                                                                             
	insertSeqInt(*it);
	
	// add exon start/end into mappedPos
	if(!(*it)->isIntron()){
	    insertPos((*it)->getSeqID(),(*it)->getStart());
	    insertPos((*it)->getSeqID(),(*it)->getEnd());
	}
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

void GeneInfoCollection::createCollection(GeneFeature *g){
    
    pair<map<string,GeneInfo>::iterator, bool> ret;

    bool hasEvidence = false;
    bool hasMatch = false;
    int numMatches = 0;
    int numEvidence = 0;
    if(g->hasEvidence())
	numEvidence++;
    
    list<pair<int,GeneFeature*> >homologs = g->homologs;
    
    int pred_idx = 0;

    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){

	if(!hit->second)
	    continue;
	
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
	    GeneInfo gi(h->getGene(),h->isCDS(),h->isIntron(),h->isExon(),false);
	    ret = ginfo[idx].insert(pair<string,GeneInfo>(h->getTxID(),gi));
	    if(!ret.second){ // transid already in map, only update
		if(h->isCDS())
		    ret.first->second.numMatchingCs++;
		else if(h->isIntron())
		    ret.first->second.numMatchingIs++;
		else if(h->isExon())
		    ret.first->second.numMatchingEs++;
	    }
	}
	pred_idx = idx;
    }
    if(hasMatch)
	numMatches++;
    if(hasEvidence)
	numEvidence++;
    
    if(g->isCDS()){
	mappedStatsC[numMatches]++;
	extrinStatsC[numEvidence]++;
	
    }
    else if(g->isExon()){
	mappedStatsE[numMatches]++;
	extrinStatsE[numEvidence]++;
    }
    else{
	mappedStatsI[numMatches]++;
	extrinStatsI[numEvidence]++;
    }
}

void GeneInfoCollection::printDetailedStats(Gene *g, std::ofstream &of){

    of << "#\n# gene feature level:\n#" << endl;
    of << "#" << endl;
    of << "# number/% of features with exact homologs in at least k other genomes:" << endl;
    of << "#" << endl;
    of << "#---------------------------------------------------------------------------------------------------" << endl;       
    of << "#   k           CDS             Intr             Exon              Intr+Exon" << endl;
    of << "#---------------------------------------------------------------------------------------------------" << endl;
    of << "#" << endl;
    writePictograph(mappedStatsC,mappedStatsI,mappedStatsE,of);
    of << "#" << endl;
    of << "# number/% of features supported by extrinsic evidence in at least k genomes :" << endl;
    of << "#" << endl;
    of << "#---------------------------------------------------------------------------------------------------" << endl;       
    of << "#   k           CDS             Intr             Exon              Intr+Exon" << endl;
    of << "#---------------------------------------------------------------------------------------------------" << endl;       
    of << "#" << endl;
    writePictograph(extrinStatsC,extrinStatsI,extrinStatsE,of);
    of << "#\n# transcript level:\n#" << endl;
    
    int numHomologs = 0;
    
    for(int i=0; i < ginfo.size(); i++){
	
	bool hasHomolog = false;
	for(map<string,GeneInfo>::iterator giit = ginfo[i].begin(); giit != ginfo[i].end(); giit++){
	    
	    int numC = giit->second.gene->numGFs(CDS);
	    int numI = giit->second.gene->numGFs(intron);
	    int numE = giit->second.gene->numGFs(exon);		
	    
	    of << "#";
	    of << setw(4) << right << i << " ";
	    of << "gene_id="<< setw(30) << left << giit->second.gene->getGeneID();
	    of << "tx_id="<< setw(30) << left << giit->first;
	    of << setw(3) << right << giit->second.numMatchingCs;
	    of << "/"<< numC << " CDS";
	    of <<setw(3) << right << giit->second.numMatchingIs;
	    of <<"/"<< numI << " Intr";
	    of <<setw(3) << right << giit->second.numMatchingEs;
	    of <<"/"<< numE << " exon" << endl;
	    if(numE == giit->second.numMatchingEs && numE == g->numGFs(exon)){
		hasHomolog = true;
		//(*git)->appendHomolog(giit->second.gene,i);		    
	    }
	    if(giit->second.numMatchingEs >=1){ // TODO: make this more reasonable, e.g. half of splice sites coincide
		g->appendHomolog(giit->second.gene,i);		    
	    }
	}
	if(hasHomolog)
	    numHomologs++;	
    }
    of << "#" << endl;
    of << "# transcript has an exact homolog in " << numHomologs << " other genomes.\n#" <<endl;    
}

void Genome::print_hgm_unaligned(GeneFeature *g, std::ofstream &of) const{

    bool isFirst = true;
    for(list<pair<int,GeneFeature*> >::iterator hit = g->homologs.begin(); hit != g->homologs.end(); hit++){
	if(!hit->second){
	    if(isFirst){
		of << hit->first;
		isFirst = false;
	    }
	    else{
		of <<","<< hit->first;
	    }
	}
    }
}	

void Genome::print_hgm_info(GeneFeature *g, std::ofstream &of) const{

    list<pair<int,GeneFeature*> >homologs = g->homologs;

    int pred_idx = 0;
    bool hasMatch = false;
    string evidence = "";
    int mult = 0;
    string onlyHint = "*";
    
    // own evidence
    of << this->idx;
    if(g->hasEvidence())
	of  << g->getEvidence() << "-" << g->getMult();
   	    
    for(list<pair<int,GeneFeature*> >::iterator hit = homologs.begin(); hit != homologs.end(); hit++){
	
	if(!hit->second)
	    continue;

	int idx = hit->first;
	GeneFeature* h = hit->second;
	
	if(pred_idx < idx && hasMatch){ // next genome
	    of <<","<< pred_idx;
	    if(!evidence.empty())
		of << evidence << onlyHint << "-" << mult;
	    hasMatch = false; // reset all flags
	    evidence = "";
	    mult = 0;
	    onlyHint = "*";
	}
	hasMatch = true;
	if(evidence.empty() && h->hasEvidence()){
	    evidence = h->getEvidence();
	    mult = h->getMult();
	}
	if(h->isPartofGene()){
	    onlyHint = "";
	}
	pred_idx = idx;
    }   
    if(hasMatch){      
	of <<","<< pred_idx;
	if(!evidence.empty())
	    of << evidence << onlyHint << "-" << mult;
    }
}

void Genome::writeTxLine(Gene *g, std::ofstream &of) const {

    of << getSeqName(g->getSeqID()) << "\t";
    of << g->getSource() << "\ttranscript\t";
    of << g->getTXstart() << "\t";
    of << g->getTXend() << "\t.\t";
    of << strandIdentifiers[g->getStrand()] << "\t.\t";
    of << g->getTxID() << endl;
}

void Genome::writeGene(Gene *g, std::ofstream &of) const{
    // print gene features
    for(list<GeneFeature*>::iterator gfit = g->features.begin(); gfit != g->features.end(); gfit++)
	writeGeneFeature(*gfit, of);
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
		list<pair<int,Gene*> >homologs = (*git)->homologs;
   
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
	for (std::vector<int>::size_type i = 0; i != component.size(); ++i){
            for (size_t j = 0; j < boost::num_vertices (G); ++j){
                if (component[j] == i)
                    of << "(" << G[j].name <<") ";
            }
            of << endl;
        }
        of.close();
    }
    else{
	cerr << "Could not open output file " << outfile << endl;
    }
}
#endif

void Genome::insertHint(string seqname, long int start, long int end, Strand strand, string esource, int mult, int frame, string f_type){
	    
    // insert new sequences both into the seqname AND seqID Hash
    pair<map<string,int>::iterator,bool> seqid;
    seqid = seqnames.insert(make_pair(seqname, seqnames.size()));
    if(seqid.second)
	seqIDs.insert(make_pair(seqid.first->second, seqid.first->first));
	     
    FeatureType type = getType(f_type);
    if(type == intron){
	start--;
	end++;
    }
    
    SeqIntKey seqInt(start, end-start+1, type);
    list<GeneFeature*> le = findSeqInt(seqInt.getKey(), seqid.first->second, strand);
    bool isGF = false;
    for(list<GeneFeature*>::iterator it=le.begin(); it!=le.end(); it++){
	if((*it)->sameFrame(frame)){
	    (*it)->setEvidence(esource); // update extrinsic source
	    (*it)->setMult(mult);
		    isGF = true;
	}
    }
    if(!isGF){ // add new gene feature
	GeneFeature *gf = new GeneFeature(type, start, end, strand, frame);
	gf->setEvidence(esource);
	gf->setMult(mult);
	insertSeqInt(gf, seqid.first->second);
	hints.push_back(gf);	 
    }    
}

#ifdef M_SQLITE
void Genome::getDbHints(SQLiteDB &db){

    try{
	Statement stmt(&db);
	
	stmt.prepare("SELECT seqname,start,end,typename,strand,frame,mult,esource FROM featuretypes JOIN \
                     (SELECT seqname,start,end,type,strand,frame,mult,esource FROM hints NATURAL JOIN \
                     (SELECT * FROM speciesnames NATURAL JOIN seqnames WHERE speciesname=$1)) ON typeid=type WHERE typename=\"CDS\" OR typename=\"exon\" OR typename=\"intron\";");

	stmt.bindText(1,name.c_str());
	
	while(stmt.nextResult()){  
	    // create new Feature	
	    string seqname = stmt.textColumn(0);
	    long int start=stmt.intColumn(1)+1; // make coordinates 1-based (DB coordinates are 0-based)
	    long int end=stmt.intColumn(2)+1;
	    Strand strand = getStrand(stmt.textColumn(4));
	    string esource = stmt.textColumn(7);
	    int mult = stmt.intColumn(6);
	    string type = stmt.textColumn(3);
	    int frame = stmt.intColumn(5);
	    insertHint(seqname,start,end,strand,esource,mult,frame,type);
	}
    }catch(const char* err) {
	cerr << err << endl;
	throw ProjectError("failed retrieving hints from DB for " + name + "."  + "\n");
    }
}
#endif
		
