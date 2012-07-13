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
#include <iostream>
#include <fstream>
#include <types.hh>

#ifdef AMYSQL
#include <table_structure.h>
#include <query.h>
#endif

MemSeqAccess::MemSeqAccess(){
    cout << "reading in file names for species from " << Constant::speciesfilenames << endl;
    filenames = getFileNames (Constant::speciesfilenames);
    /*
     * reading in sequences into memory
     */
    for(map<string, string>::iterator it = filenames.begin(); it != filenames.end(); it++){
	GBProcessor gbank(it->second);
	AnnoSequence *inSeq = gbank.getSequenceList();
	while(inSeq){
	    string key = it->first + ":" + inSeq->seqname;
	    sequences[key] = inSeq->sequence;
	    inSeq = inSeq->next;
	}
    }
}

AnnoSequence* MemSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
    AnnoSequence *annoseq = NULL;
    string key = speciesname + ":" + chrName;
    map<string,char*>::iterator it = sequences.find(key);
    if(it != sequences.end()){
	annoseq = new AnnoSequence();
	annoseq->seqname = newstrcpy(chrName);
	annoseq->sequence = newstrcpy(it->second + start, end - start + 1);
	annoseq->length = end-start+1;
	annoseq->offset = start;
	if(strand == minusstrand){
	    char *reverseDNA = reverseComplement(annoseq->sequence);
	    delete [] annoseq->sequence;
	    annoseq->sequence = reverseDNA;
	}
    }
    return annoseq;
}


map<string,string> getFileNames (string listfile){
    map<string,string> filenames;
    ifstream ifstrm(listfile.c_str());
    if (ifstrm.is_open()){
	string line;
	while(getline(ifstrm, line)){
	    size_t pos = line.find('\t');
	    if (pos != string::npos)
		filenames[line.substr(0,pos)] = line.substr(pos + 1) ;
	    else
		throw ProjectError(listfile + " has wrong format in line " + line + ". correct format:\n\n" + 
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



DbSeqAccess::DbSeqAccess(){
#ifdef AMYSQL
    cout << "opening database connection using connection data\"" << Constant::dbaccess << "\""<<endl;
    dbaccess = Constant::dbaccess;
    split_dbaccess();
    connect_db();
#else
    throw ProjectError("Database access not possible with this compiled version. Please recompile with flag MYSQL.");
#endif
}

#ifndef AMYSQL
AnnoSequence* DbSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
    return NULL;
    // empty dummy for compiler, error message is created in constructor
}
#else // AMYSQL

/*
 * coord_id is a identifier in table 'seq_region'.
 * coord_id==2: this is a contig.there's an entity in 'dna' table,you can retrive sequence directly from it.
 * coord_id==1: this is a chromosome that is consist of more than one entities in 'dna' table.You can't
 * retrieve sequence directly from 'dna' table.The components id and order in which
 * they are assembled can be found in table 'assembly'.
 */
AnnoSequence* DbSeqAccess::getSeq(string speciesname, string chrName, int start, int end, Strand strand){
    mysqlpp::StoreQueryResult store_res;
    AnnoSequence* annoseq=NULL;
    int coord_id,seq_region_id,seq_region_length; 
    vector<assembly> asm_query_region;
    mysqlpp::Query detect_coord_id=con.query();
    detect_coord_id<<"select seq_region_id,coord_system_id,length from seq_region where name=\""
		   <<chrName<<"\"";
    store_res=detect_coord_id.store();
    try{
	if(store_res.num_rows()==0){
	    cerr<<"DbSeqAccess::getSeq : chrName\"" <<chrName
		<<"\" not exist in Database,retrive sequence failed."<<endl;
	}
	else{
	    seq_region_id=store_res[0][0];
	    coord_id=store_res[0][1];
	    seq_region_length=store_res[0][2];
/* 
 * in different datababses,'coord_id' are defined differently.
 * get a short substring,if query succeed，the query sequence is a 'contig'
 * otherwise it's a 'chromosome'
 */
	    detect_coord_id<<"select substring(sequence from 1 for 10) from dna where seq_region_id="<<seq_region_id<<endl;
	    store_res=detect_coord_id.store();
	    if(store_res.size()==0){
		coord_id=1;
	    }
	    else{
		coord_id=2;
	    }
	    if(end==-1){// predictionEnd is not defiend.
		end=seq_region_length;
	    }
	    if(start==0){// predictionEnd is not defiend.
		++start;
	    }
	    if(coord_id==1){
		get_region_coord(seq_region_id,start,end,asm_query_region);
	    }
/*
 *This record can't be found in table 'assembly'.
 *Just for the convenient to put a vector<assembly> to Function:getNextDBSequence.
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

int DbSeqAccess::split_dbaccess(){
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

void DbSeqAccess::connect_db(){
    const  char* db_name = db_information[0].c_str();
    const  char* host = db_information[1].c_str();
    const  char* user = db_information[2].c_str();
    const  char* passwd = db_information[3].c_str();
    try{
	con.connect(db_name,host,user,passwd);
	cout << "DB connection OK:" << user << "\ton " << host << "\tconnects to " << db_name << "\n";
    }
    catch(const mysqlpp::BadQuery& er){
	cout << "Query error: " << er.what() << endl;
    }
}


template<class T>
AnnoSequence* DbSeqAccess::getNextDBSequence(string chrName,int start,int end,vector<T>& asm_query_region)
{
    AnnoSequence* annoseq=new AnnoSequence();
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
    char* tmp_sequence=new char[concat_str.size()+1];
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
int DbSeqAccess::get_region_coord(int seq_region_id,int start,int end,vector<T> &asm_query_region){
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
//AnnoSequence* DbSeqAccess::getDBSequenceList(string chrName,int start,int end,vector<T>& asm_query_region)
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
