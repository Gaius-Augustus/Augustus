/**********************************************************************
 * file:    randseqaccess.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  random acces to sequence data, e.g. get me chr1:1000-2000 from species 'human'
 * authors: Mario Stanke
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 07.03.12| Mario Stanke  | creation of the file
 * 21.03.12| Stefanie König| implementation of MemSeqAcess functions
 * 09.06.12| yuqiulin      | mysql access to sequences
 **********************************************************************/

#include "randseqaccess.hh"
#include "genbank.hh"
#include "intronmodel.hh"
#include "igenicmodel.hh"
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <types.hh>

#ifdef AMYSQL
#include <table_structure.h>
#include <query.h>
#endif

int SpeciesCollection::groupCount = 1;

FeatureCollection* SpeciesCollection::getFeatureCollection(string speciesname){
    int groupID = getGroupID(speciesname);
    if(groupID <= 0)
	return &defaultColl;
    map<int,FeatureCollection>::iterator it = speciesColl.find(groupID);
    if(it != speciesColl.end())
	return &(it->second);
    return &defaultColl;
}

int SpeciesCollection::getGroupID(string speciesname){
    map<string,int>::iterator it = groupIDs.find(speciesname);
    if(it != groupIDs.end())
	return it->second;
    return 0;
}

void SpeciesCollection::addSpeciesToGroup(string skey, int groupID){

    map<string,int>::iterator it=groupIDs.find(skey);
    if(it != groupIDs.end())
	throw ProjectError("SpeciesCollection::readExtrinsicCFGFile: species " +
			   skey +" is assigned to more than one extrinsic config table\n");
    groupIDs.insert(pair<string,int>(skey,groupID));
}

void SpeciesCollection::readExtrinsicCFGFile(vector<string> &speciesNames){
    string filename,skey;
    char buf[1024];
    ifstream datei;

    try {
        filename = Properties::getProperty( "extrinsicCfgFile" );
    } catch(...) {
        cerr << "Could not find parameter 'extrinsicCfgFile'" << endl;
        return;
    }

    try {
        datei.open(filename.c_str());
        if (!datei) {
            throw ProjectError(string("Could not find extrinsic config file ") + filename + ".");
        }
        datei >> comment;
	FeatureCollection fc;
	fc.readSourceRelatedCFG(datei);
	while(datei){
	    FeatureCollection sc(fc);
	    sc.readTypeInfo(datei); // reading in a bonus/malus table
	    if(datei){
		// reading in the species group for which the table is valid
		getline(datei,skey);
		if(skey == "[GROUP]"){
		    cout << "extrinsic group " << groupCount << ":";
		    datei.getline(buf, 1024); // reading in the set of species that belongs to the group
		    stringstream stm(buf);
		    if(stm >> skey){ 
			do{
			    if(skey == "all" || skey == "other"){ // this config table is valid for all other species
				for(int i=0; i<speciesNames.size(); i++){
				    try{
					addSpeciesToGroup(speciesNames[i],groupCount);
					cout <<" "<<speciesNames[i];
				    }catch(ProjectError e){
					if(skey == "all")
					    throw;
				    }
				}
			    }
			    else{
				if(find(speciesNames.begin(), speciesNames.end(), skey) == speciesNames.end()) {
				    throw ProjectError(skey + " is not a valid species identifier.");
				}
				addSpeciesToGroup(skey,groupCount);
				cout <<" "<<skey;
			    }
			}while(stm >> skey);
			cout << endl;
		    }
		    else{
			throw ProjectError("SpeciesCollection::readExtrinsicCFGFile: Please specify a set of species for which config table " + 
					   itoa(groupCount) +" in\n " + filename + "is valid");
		    }
		    speciesColl.insert(pair<int,FeatureCollection>(groupCount,sc));
		    groupCount++;
		    while(datei >> comment >> ws, datei && datei.peek() != '[')
			;
		}
		else{
		    throw ProjectError("SpeciesCollection::readExtrinsicCFGFile: Please specify a set of species for which config table " +
				       itoa(groupCount) +" in\n " + filename + "is valid");
		}
	    }
	    else{
		throw ProjectError("SpeciesCollection::readExtrinsicCFGFile: Please specify a set of species for which config table " +
				   itoa(groupCount) + " in\n " + filename + "is valid");
	    }
	}
        datei.close();
        datei.clear();
    } catch (ProjectError e) {
        cerr << e.getMessage() << endl;
        cerr << "Could not read in file with the configuration of hints: " << filename << endl;
        datei.close();
        datei.clear();
        throw;
    }
}

void SpeciesCollection::readGFFFile(const char *filename){
    /*
     *Read in the configuration file for extrinsic features.
     */
    ifstream datei;
    try {
        datei.open(filename);
        if( !datei ) {
            cerr << "FeatureCollection::readGFFFile( " << filename << " ) : Could not open the file!!!" << endl;
            throw ProjectError();
        }

	/* 
	 * read in line by line
	 */
	Feature f;
	datei >> comment >> ws;
        while (datei) {
            try{
                datei >> f >> comment >> ws;
            } catch (ProjectError e){}
	    // split species name and sequence ID
	    string completeName=f.seqname;
	    string speciesName,seqname;                                                                                                          
	    for (int i=0; i<completeName.length(); i++) {
	       // seperator is the point '.' for example hg19.chr21, has to be changed                         
	       if ((completeName[i] == '-') || (completeName[i] == '.')) {
		  speciesName = completeName.substr(0,i);
		  seqname = completeName.substr(i+1, string::npos);
		  if (i == completeName.length()-1)
		     throw ProjectError("first column in hintfile must be the speciesname and seqname delimited by '.'");
	       }
	    }
	    f.seqname=seqname;
	    // get species specific collection
	    if(withEvidence(speciesName)){
	       FeatureCollection *fc = getFeatureCollection(speciesName);
	       fc->setBonusMalus(f);
	       SequenceFeatureCollection& sfc =  fc->getSequenceFeatureCollection(completeName.c_str());
	       sfc.addFeature(f);
	       fc->hasHintsFile = true;
	    }
	    else{
	       cerr << "Warning: hints are given for species " + speciesName +
		  " but no extrinsic configuration in the extrinsicCfgFile.\n" +
		  " Ignoring all hints for that species." << endl;
			
	    }
        }
        datei.close();
    } catch (ProjectError e) {
        cerr << e.getMessage() << endl;
        throw e;
    } catch(...) {
        cerr << "FeatureCollection::readGFFFile( " << filename << " ) : Could not read the file!!!" << endl;
        datei.close();
    }
}

void RandSeqAccess::setLength(int idx, string chrName, int len){
    map<string,int>::iterator it = chrLen[idx].find(chrName);
    if (it == chrLen[idx].end()){
	chrLen[idx].insert(pair<string,int>(chrName, len));
    } else {
	if (it->second != len)
	    throw ProjectError("Lengths of " + chrName + " inconsistent.");
    }
}

void RandSeqAccess::setSpeciesNames(vector<string> speciesNames){
    numSpecies = speciesNames.size();
    this->speciesNames = speciesNames;
    chrLen.resize(numSpecies, map<string,int>()); // empty initialization
    for (size_t i=0; i < numSpecies; i++){
	if (speciesIndex.find(speciesNames[i]) != speciesIndex.end())
	    throw ProjectError(string("List of species names contains multiple entries: ") + speciesNames[i]);
	speciesIndex.insert(pair<string, size_t>(speciesNames[i], i));
    }
}

int RandSeqAccess::getMaxSnameLen(){
    int maxNameLen = 0;
    for (size_t s=0; s<numSpecies; s++)
	if (speciesNames[s].length() > maxNameLen)
	    maxNameLen = speciesNames[s].length();
    return maxNameLen;
}


int RandSeqAccess::getIdx(string speciesname) {
    map<string,size_t>::iterator it = speciesIndex.find(speciesname);
    if (it == speciesIndex.end())
	return -1;
    else 
	return it->second;
}


int RandSeqAccess::getChrLen(int idx, string chrName){
    map<string,int>::iterator it = chrLen[idx].find(chrName);
    if(it != chrLen[idx].end()) {
	return it->second;
    } else {
	cout << "RandSeqAccess::getChrLen failed on sequence " << chrName << " from species " << speciesNames[idx] << endl;
	return -1; // TODO: error throwing
    }
}

void RandSeqAccess::printStats(){
    cout << "number of species" << numSpecies << endl;
    for (int s=0; s<numSpecies; s++){
	cout << "species " << setw(2) << s << ": " << speciesNames[s] << "\tspeciesIndex= " << getIdx(speciesNames[s]) << endl;
	cout << "sequence lengths for species " << speciesNames[s] << " from alignment file:" << endl;
	for (map<string,int>::iterator it = chrLen[s].begin(); it != chrLen[s].end(); ++it)
	    cout << it->first << " => " << it->second << '\n';
    }
}

MemSeqAccess::MemSeqAccess(vector<string> s){

    setSpeciesNames(s);

    cout << "reading in file names for species from " << Constant::speciesfilenames << endl;
    filenames = getFileNames (Constant::speciesfilenames);
    /*
     * reading in sequences into memory
     */
    for(map<string, string>::iterator it = filenames.begin(); it != filenames.end(); it++){
	GBProcessor gbank(it->second);
	AnnoSequence *inSeq = gbank.getSequenceList();
	while(inSeq){
	    string key = it->first + "." + inSeq->seqname;
	    cout<< "reading in "<<key<<endl;
	    sequences[key] = inSeq->sequence;
	    inSeq = inSeq->next;
	}
    }
    /*
     * reading in exintrinsic evidence into memory
     */
    const char *extrinsicfilename;
    try {
	extrinsicfilename =  Properties::getProperty("hintsfile");
    } catch (...){
	extrinsicfilename = NULL;
	 }
    if (extrinsicfilename || Constant::softmasking) {

	extrinsicFeatures.readExtrinsicCFGFile(speciesNames);
    
	if(extrinsicfilename){
	    cout << "# reading in the file " << extrinsicfilename << " ..." << endl;
	    extrinsicFeatures.readGFFFile(extrinsicfilename);
	    
	    // print Seqs for which hints are given
	    cout << "We have hints for"<<endl;
	    bool seqsWithInfo = false;
	    for(map<string, char*>::iterator it = sequences.begin(); it != sequences.end(); it++){
		string completeName=it->first;
		size_t pos = completeName.find('.');
		string speciesname = completeName.substr(0,pos);
		FeatureCollection *fc = extrinsicFeatures.getFeatureCollection(speciesname);
		if(fc->isInCollections(completeName)){
		    seqsWithInfo = true;
		    cout << completeName << endl;
		}
	    }
	    if(!seqsWithInfo){
		cout << "# WARNING: extrinsic information given but not on any of the sequences in the input set!" << endl;
		cout << "The first column in the hints file must contain the speciesID and seqID separated by '.'"<< endl;
		cout << "(for example 'hg19.chr21')" << endl;
	    }
	}
    }
    else{
	cout << "# No extrinsic information given." << endl;
    }
}

AnnoSequence* MemSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
    AnnoSequence *annoseq = NULL;
    string key = speciesname + "." + chrName;
    map<string,char*>::iterator it = sequences.find(key);
    if(it != sequences.end()){
	annoseq = new AnnoSequence();
	annoseq->seqname = newstrcpy(chrName);
	annoseq->sequence = newstrcpy(it->second + start, end - start + 1);
	annoseq->length = end - start + 1;
	annoseq->offset = start;
	if(strand == minusstrand){
	    char *reverseDNA = reverseComplement(annoseq->sequence);
	    delete [] annoseq->sequence;
	    annoseq->sequence = reverseDNA;
	}
    }
    return annoseq;
}

SequenceFeatureCollection* MemSeqAccess::getFeatures(string speciesname, string chrName, int start, int end, Strand strand){

    string key = speciesname + "." + chrName;
    FeatureCollection *fc = extrinsicFeatures.getFeatureCollection(speciesname);
    SequenceFeatureCollection *c=fc->getSequenceFeatureCollection(key);
    if(!c){
	SequenceFeatureCollection* sfc = new SequenceFeatureCollection(fc); // empty list of hints
	return sfc;
    }
    bool rc=false;
    if(strand == minusstrand)
	rc=true;
    SequenceFeatureCollection* sfc = new SequenceFeatureCollection(*c,start,end,rc); // all hints that end in [start,end]
    return sfc;
}


map<string,string> getFileNames (string listfile){
    map<string,string> filenames;
    ifstream ifstrm(listfile.c_str());
    if (ifstrm.is_open()){
	char buf[256];
	while(ifstrm.getline(buf,255)){
	    stringstream stm(buf);
	    string species, dir;
	    if(stm >> species >> dir){
		dir = expandHome(dir);
		filenames[species] = dir;
	    }
	    else
		throw ProjectError(listfile + " has wrong format in line " + buf + ". correct format:\n\n" + 
				   "Homo sapiens <TAB> /dir/to/genome/human.fa\n" + 
				   "Mus musculus <TAB> /dir/to/genome/mouse.fa\n" + 
				   "...\n");
	}
	ifstrm.close();
    }
    else
        throw ProjectError("Could not open input file " + listfile);
    
    return filenames;
}

DbSeqAccess::DbSeqAccess(vector<string> s){
    setSpeciesNames(s);
  
    if(Constant::dbhints ||  Constant::softmasking ){
	cout << "read in the configuration file for extrinsic features" << endl;
	extrinsicFeatures.readExtrinsicCFGFile(speciesNames);
    }
}

#ifdef AMYSQL
void MysqlAccess::open(){
    dbaccess = Constant::dbaccess;
    split_dbaccess();
    if(Constant::MultSpeciesMode)
	connect_db();
    else
	connect_db(cerr);
}

/*
 * retrieve sequence directly from table genomes(seqid, dnaseq, seqname, start, end, species)
 * arguments and columns start and end are 0-based.
 *
 * database chunks:   |-------------||-------------||-------------||-------------||-------------|
 * requested segment:                      |--------------------------|
 *                                   |   start                       end
 *                             g[0].start                          | 
 */
AnnoSequence* MysqlAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
    mysqlpp::StoreQueryResult store_res;
    string dna, querystr;
    mysqlpp::Query query = con.query();
    query << "SELECT dnaseq,start,end FROM genomes as G,speciesnames as S,seqnames as N WHERE speciesname='" << speciesname << "' AND seqname='"
	  << chrName << "' AND G.speciesid=S.speciesid AND S.speciesid=N.speciesid AND G.seqnr=N.seqnr AND start <= " << end << " AND end >= " << start << " ORDER BY start ASC";
    querystr = query.str();
    //cout << "Executing" << endl << querystr << endl;
    vector<genomes> g;
    query.storein(g);
    
    
    AnnoSequence* annoseq = new AnnoSequence();
    annoseq->seqname = newstrcpy(chrName);
    if (g.empty())
	throw ProjectError("Could not retrieve sequence from database using query:" + querystr);
    else if (g.size() == 1){ // segment overlaps a single dna chunk
        if (!(g[0].start <= start && g[0].end >= end)){
	    if(Constant::MultSpeciesMode){
		return NULL; // temporaryily ignore this problem introduced by Alexander
		throw ProjectError("Tried to retrieve a sequence that is only partially contained in database:"
			       + chrName + ":" + itoa(start) + "-" + itoa(end));
	    }
	    else{
		if(end > g[0].end)
		    end=g[0].end;
	    }
        }
	dna = g[0].dnaseq.substr(start - g[0].start, end-start+1);
    } else {
	vector<genomes>::iterator it = g.begin();
	if (it->end >= end)
	    throw ProjectError("Segment not uniquely represented in database. Have you loaded sequences more than once?");
	dna = it->dnaseq.substr(start - g[0].start);
	while (it != g.end()) {
	    if ((it+1) != g.end() && it->end + 1 != (it+1)->start)
		throw ProjectError("Internal error. Genome sequence not sliced seamlessly into chunks.");
	    it++;
	    if (it != g.end()){
		if ((it+1) != g.end()) {
		    if (it->end >= end)
			throw ProjectError("Segment not uniquely represented in database. Have you loaded sequences more than once?");
		    dna.append(it->dnaseq);
		} else { // last chunk
		    if (it->end < end){
			if(Constant::MultSpeciesMode){
			    throw ProjectError("Tried to retrieve a sequence that is only partially contained in database:" + chrName + ":" + itoa(start) + "-" + itoa(end));
			}
			else{
			    end = it->end;
			}
		    }
		    dna.append(it->dnaseq.substr(0, end - it->start + 1)); // rest of sequence = initial part of last chunk
		}
	    }
	}
    }
    annoseq->sequence = newstrcpy((strand == minusstrand)? reverseComplement(dna.c_str()) : dna.c_str());
    annoseq->length = end - start + 1;
    annoseq->offset = start;
    return annoseq;
}


/*
 * coord_id is an identifier in table 'seq_region'.
 * coord_id==2: this is a contig. There's an entity in 'dna' table, you can retrive sequence directly from it.
 * coord_id==1: this is a chromosome that is consist of more than one entities in 'dna' table. You can't
 * retrieve sequence directly from 'dna' table. The components id and order in which
 * they are assembled can be found in table 'assembly'.
 */
AnnoSequence* MysqlAccess::getSeq2(string speciesname, string chrName, int start, int end, Strand strand){
    mysqlpp::StoreQueryResult store_res;
    AnnoSequence* annoseq = NULL;
    int coord_id, seq_region_id, seq_region_length; 
    vector<assembly> asm_query_region;
    mysqlpp::Query detect_coord_id = con.query();
    detect_coord_id << "select seq_region_id,coord_system_id,length from seq_region where name=\""
		    << chrName << "\"";
    store_res = detect_coord_id.store();
    try {
	if(store_res.num_rows() == 0 ){
	    cerr << "MysqlAccess::getSeq : chrName\"" << chrName
		 << "\" does not exist in database, retrieval of sequence failed." << endl;
	} else {
	    seq_region_id = store_res[0][0];
	    coord_id = store_res[0][1];
	    seq_region_length = store_res[0][2];
/* 
 * In different databases, 'coord_id' are defined differently.
 * get a short substring, if query succeed, the query sequence is a 'contig'
 * otherwise it's a 'chromosome'
 */
	    detect_coord_id << "select substring(sequence from 1 for 10) from dna where seq_region_id="
			    << seq_region_id << endl;
	    store_res = detect_coord_id.store();
	    if(store_res.size()==0){
		coord_id=1;
	    } else {
		coord_id=2;
	    }
	    if(end == -1){// predictionEnd is not defiend.
		end = seq_region_length;
	    }
	    if(start == 0){// predictionStart is not defiend.
		++start;
	    }
	    if(coord_id == 1){
		get_region_coord(seq_region_id,start,end,asm_query_region);
	    }
/*
 * This record can't be found in table 'assembly'.
 * Just for the convenient to put a vector<assembly> to Function:getNextDBSequence.
 */
	    if(coord_id==2){
		assembly row(seq_region_id,seq_region_id,start,end,start,end);
		asm_query_region.push_back(row);
	    }
	}
    }
    catch(const mysqlpp::BadQuery& er){
	cout << "Query error: "<<er.what()<<endl;
    }
    annoseq=getNextDBSequence(chrName,start,end,asm_query_region);
    if(strand == minusstrand){
        char *reverseDNA = reverseComplement(annoseq->sequence);
        delete [] annoseq->sequence;
	annoseq->sequence = reverseDNA;
    }
    return annoseq;
}

SequenceFeatureCollection* MysqlAccess::getFeatures(string speciesname, string chrName, int start, int end, Strand strand){

    FeatureCollection* fc = extrinsicFeatures.getFeatureCollection(speciesname);
    SequenceFeatureCollection* sfc = new SequenceFeatureCollection(fc);
    if(extrinsicFeatures.withEvidence(speciesname) && Constant::dbhints){ // only retrieve hints for the species specified in the extrinsicCfgFile
	mysqlpp::Query query = con.query();
	query << "SELECT source,start,end,score,type,strand,frame,priority,grp,mult,esource FROM hints as H, speciesnames as S,seqnames as N WHERE speciesname='"
	      << speciesname << "' AND seqname='" << chrName << "' AND H.speciesid=S.speciesid AND S.speciesid=N.speciesid AND H.seqnr=N.seqnr AND start <= "
	      << end << " AND end >= " << start;
	//cout << "Executing" << endl << query.str() << endl;
	vector<hints> h;
	query.storein(h);
	
	if(h.empty()){
	    //cout << "no hints retrieved"<<endl;
	}
	else{
	    for(vector<hints>::iterator it = h.begin(); it != h.end(); it++){
		
		// create new Feature
		Feature f;
		f.seqname=chrName;
		f.source=it->source;
		f.type=Feature::getFeatureType(it->type);
		f.feature=featureTypeNames[f.type];
		f.start=it->start;
		f.end=it->end;
		f.score=it->score;
		f.setFrame(it->frame);
		f.setStrand(it->strand);
		f.groupname=it->grp;
		f.priority=it->priority;
		f.mult=it->mult;
		f.esource=it->esource;
		
		// shift positions relative to gene Range
		if(strand == plusstrand)
		    f.shiftCoordinates(start,end);
		else
		    f.shiftCoordinates(start,end,true);
		fc->setBonusMalus(f);
		sfc->addFeature(f);
	    }
	}
	fc->hasHintsFile = true;
    }
    return sfc;
}


int MysqlAccess::split_dbaccess(){
    string::size_type pos1, pos2;
    pos2 = dbaccess.find(','); //string 'dbaccess' is delimited by ',' as default.
    pos1 = 0;        
    while (string::npos != pos2) {
	db_information.push_back(dbaccess.substr(pos1, pos2 - pos1));
	pos1 = pos2 + 1;
	pos2 = dbaccess.find(',', pos1);
    }
    db_information.push_back(dbaccess.substr(pos1));
    return 0;
}

void MysqlAccess::connect_db(ostream& out){
    const  char* db_name = db_information[0].c_str();
    const  char* host = db_information[1].c_str();
    const  char* user = db_information[2].c_str();
    const  char* passwd = db_information[3].c_str();
    try {
	out << "# Opening database connection using connection data \"" << Constant::dbaccess << "\"...\t";
	con.connect(db_name,host,user,passwd);
	con.set_option(new mysqlpp::ReconnectOption(true));
	out << "DB connection OK." << endl;
    }
    catch(const mysqlpp::BadQuery& er){
	out << "Query error: " << er.what() << endl;
    }
}


template<class T>
AnnoSequence* MysqlAccess::getNextDBSequence(string chrName,int start,int end,vector<T>& asm_query_region)
{
    AnnoSequence* annoseq = new AnnoSequence();
    mysqlpp::Query fetchseq_query=con.query();
    mysqlpp::StoreQueryResult mysqlseq; // store the chunk in a StoreQueryResult container.

    //step1: concatenate all the chunks to a single sequence.
    string concat_str;
    int chunkstart,chunkend,fetchlength,gaplength;
    vector<assembly>::iterator iter=asm_query_region.begin();
    int tail=iter->asm_start-1;
    while(iter!=asm_query_region.end()){
	chunkstart=iter->cmp_start;
	chunkend=iter->cmp_end;
	fetchlength=chunkend-chunkstart+1;
	gaplength=iter->asm_start-tail-1;
	concat_str.append(gaplength,'n');//The gaps between chunks are filled with 'n'.
	fetchseq_query<<"select substring(sequence from "
		      <<chunkstart<<" for "<<fetchlength
		      <<") from dna where seq_region_id="
		      <<iter->cmp_seq_region_id;
	try{
	    mysqlseq=fetchseq_query.store();
	    if(mysqlseq.num_rows()==0){
	    cerr<<"get_region_coord Error: No 'dna' corresponds to "<<chrName<<" from "<<chunkstart<<" to "<<chunkend<<endl;
	    }
	    else{
		concat_str.append(mysqlseq[0][0]);
	    }
	}
	catch(const mysqlpp::BadQuery& er){
	    cout << "Query error: "<<er.what()<<endl;
	}
	tail=iter->asm_end;
	iter++;
    }
    //step2: give concatenate sequence to a AnnoSequence object.
    char* tmp_sequence = new char[concat_str.size()+1];
    int pos=-1;
    int i;
    int iters=concat_str.size();
    for (i=0;i<=iters;++i){ 
   	if (isalpha(concat_str[i])){
   	    tmp_sequence[++pos]=tolower(concat_str[i]);
   	}
    }
    tmp_sequence[++pos] = '\0';
    annoseq->sequence=tmp_sequence;
    annoseq->length =end-start+1;
    annoseq->offset=start-1; //predictionStart/End from cmdline start from 1,make it 0-offset here.
    annoseq->seqname=newstrcpy(chrName.c_str());
    return annoseq;
}


/* select first and final segment in 'assembly' table that decide the boundaries of query sequence.
 * The trunks of a query sequence are adjacent non-overlaping dna's segments stored in 'assembly' table.
 *    |atct....|atg........|..|....|......abt|   
 *       START|  query sequence range  |END   
 */
template<class T>
int MysqlAccess::get_region_coord(int seq_region_id,int start,int end,vector<T> &asm_query_region){
    mysqlpp::Query get_region_coord=con.query();
    try{
	get_region_coord<<"select * from assembly where asm_seq_region_id=\""<<seq_region_id<<"\""
			<<" and asm_start <= "<<end
			<<" and asm_end >= "<<start; //assume the trunks store in table 'assembly' are ASC sorted.
	get_region_coord.storein(asm_query_region);

	if(asm_query_region.size()>0){
	    int offset=start-(asm_query_region.begin())->asm_start;
	    if(offset<0){
		cout<<"get_region_coord Warning:chunksize out of range,chunk to "<<(asm_query_region.begin())->asm_start<<"on seq ID:"<<seq_region_id;
	    }
	    (asm_query_region.begin())->asm_start=start;
	    (asm_query_region.begin())->cmp_start += offset;

	    offset=end-(asm_query_region.rbegin())->asm_end;
	    if(offset>0){
		cout<<"get_region_coord Warning:chunksize out of range,chunk to "<<(asm_query_region.rbegin())->asm_end<<"on seq ID:"<<seq_region_id;
	    }
	    (asm_query_region.rbegin())->asm_end=end;
	    (asm_query_region.rbegin())->cmp_end += offset;
	}
	else{
	    cerr<<"get_region_coord Error: No 'dna' corresponds to seq ID:"<<seq_region_id<<" from "<<start<<" to "<<end<<endl;
	}
    }
    catch(const mysqlpp::BadQuery& er){
	cout << "Query error: "<<er.what()<<endl;
    }
    return 0;
}

//template<class T>
//AnnoSequence* MysqlAccess::getDBSequenceList(string chrName,int start,int end,vector<T>& asm_query_region)
//{
//    int chunkcount=1;
//    string tmpstring;
//    char* tempchar;
//    AnnoSequence *seqlist = NULL, *seq, *last = NULL; // seqlist is the header.
//    vector<assembly>::const_iterator asm_iter=asm_query_region.begin(); 
//    int tail=asm_iter->asm_end;
//	if(last){
//	    if(gap>MAX_GAP_LENGTH){
//		string chunkname;
//		chunkname=chrName+".chunk_"+ itoa(start) + "-" + atoi(end);
//		cout<<"************chunk name is "<<chunkname;
//		    //bin    int chunkstart;
//		    //    int chunkend;
////	    seq = getNextDBSequence(chunkname,);
////	seq->next = NULL;
//		    //if (last)
////	  last->next = seq;
//	}
//	else{
//	      seqlist = seq;
//	}
//	last = seq;
//	//++asm_iter;
//    }
//    if (seqlist == NULL)
//	throw ProjectError("No sequences found.");
//    return seqlist;
//}

#endif // AMYSQL

#ifdef SQLITE
AnnoSequence* SQLiteAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){

    int seq_start, seq_end;
    streampos file_start;
    streamsize n=0; // number of characters that are read
    
    AnnoSequence* annoseq = new AnnoSequence();
    annoseq->seqname = newstrcpy(chrName);

    try{
        Statement stmt(&db);
        stmt.prepare("SELECT start,end, offset, streamsize FROM genomes as G,speciesnames as S,seqnames as N WHERE speciesname=$1 AND seqname=$2 \
        AND G.speciesid=S.speciesid AND S.speciesid=N.speciesid AND G.seqnr=N.seqnr AND start <=?3 AND end >=?4 ORDER BY start ASC;");
        stmt.bindText(1,speciesname.c_str());
        stmt.bindText(2,chrName.c_str());
	stmt.bindInt(3,end);
	stmt.bindInt(4,start);
	if(stmt.nextResult()){  
	    seq_start = stmt.intColumn(0);
	    seq_end = stmt.intColumn(1);
	    file_start = (std::streampos)stmt.int64Column(2);
	    n += stmt.intColumn(3);
	    while(stmt.nextResult()){
		seq_end = stmt.intColumn(1);
		n += stmt.intColumn(3);
	    }
	    if(start - seq_start < 0 || seq_end - end < 0){
		if(Constant::MultSpeciesMode){
		    throw ProjectError("failed retrieving sequence " + speciesname + "." + chrName + ":" + itoa(start) + "-" + itoa(end) + "\nout of range error\n");
		}
		else{
		    if(end > seq_end)
			end = seq_end;
		}		    
	    }

	    map<string,string>::iterator it = filenames.find(speciesname);
	    if(it == filenames.end()){
		throw ProjectError("no genome file specified for " + speciesname + " in speciesfilenames");
	    }
	    string file =it->second;
	    ifstream ifstrm(file.c_str());
	    if (ifstrm.is_open()){
		ifstrm.seekg(file_start);
		char* seq = new char[n+1];
		ifstrm.read(seq,n);
		seq[n] = '\0';
		ifstrm.close();

		strip_newlines(seq); // strip of newline characters
		
		annoseq->sequence = newstrcpy(seq + (start-seq_start), end-start+1);
		annoseq->length = end - start + 1;
		annoseq->offset = start;
		delete [] seq;
		if(strand == minusstrand){
		    char *reverseDNA = reverseComplement(annoseq->sequence);
		    delete [] annoseq->sequence;
		    annoseq->sequence = reverseDNA;
		}
		return annoseq;
	    }
	    else
		throw ProjectError("Could not open input file " + file);
	}
	else{
	    throw ProjectError("failed retrieving sequence " + speciesname + "." + chrName + ":" + itoa(start) + "-" + itoa(end) + "\nno such sequence in the db.\n");

	}
    }catch(const char* err){
	cerr<<db.error()<< endl;
	throw ProjectError("failed retrieving sequence " + speciesname + "." + chrName + ":" + itoa(start) + "-" + itoa(end) + "\nsqlite error\n");
    }
    return NULL;
}

SequenceFeatureCollection* SQLiteAccess::getFeatures(string speciesname, string chrName, int start, int end, Strand strand){

    FeatureCollection* fc = extrinsicFeatures.getFeatureCollection(speciesname);
    SequenceFeatureCollection* sfc = new SequenceFeatureCollection(fc);
    if(extrinsicFeatures.withEvidence(speciesname) && Constant::dbhints){ // only retrieve hints for the species specified in the extrinsicCfgFile
	try{
	    Statement stmt(&db);
	    stmt.prepare("SELECT source,start,end,score,type,strand,frame,priority,grp,mult,esource FROM hints as H, speciesnames as S,seqnames as N WHERE \
                      speciesname=?1 AND seqname=?2 AND H.speciesid=S.speciesid AND S.speciesid=N.speciesid AND H.seqnr=N.seqnr AND start <=?3 AND end >=?4;");
	    stmt.bindText(1,speciesname.c_str());
	    stmt.bindText(2,chrName.c_str());
	    stmt.bindInt(3,end);
	    stmt.bindInt(4,start);
	    while(stmt.nextResult()){  
		
		// create new Feature
		Feature f;
		f.seqname=chrName;
		f.source=stmt.textColumn(0);
		f.type=Feature::getFeatureType(stmt.intColumn(4));
		f.feature=featureTypeNames[f.type];
		f.start=stmt.intColumn(1);
		f.end=stmt.intColumn(2);
		f.score=stmt.doubleColumn(3);
		f.setFrame(stmt.textColumn(6));
		f.setStrand(stmt.textColumn(5));
		f.groupname=stmt.textColumn(8);
		f.priority=stmt.intColumn(7);
		f.mult=stmt.intColumn(9);
		f.esource=stmt.textColumn(10);
		
		// shift positions relative to gene Range
		if(strand == plusstrand)
		    f.shiftCoordinates(start,end);
		else
		    f.shiftCoordinates(start,end,true);
		fc->setBonusMalus(f);
		sfc->addFeature(f);
	    }
	}catch(const char* err) {
	    cerr << err << endl;
	    throw ProjectError("failed retrieving hints on " + speciesname + "." + chrName + ":" + itoa(start) + "-" + itoa(end) + "\n");
	}
    } 
    fc->hasHintsFile = true;
    return sfc;
}
#endif
