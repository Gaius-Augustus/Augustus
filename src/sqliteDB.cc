/**********************************************************************
 * file:    sqliteDB.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  wrapper class around the SQLite interface
 * authors: Stefanie Koenig
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 25.06.14|Stefanie Koenig| creation of the file
 **********************************************************************/


#include "sqliteDB.hh"
#include "hints.hh"

#include <iostream>
#include <stdlib.h>

void SQLiteDB::open(OpenMode mode){
    
    int flag = SQLITE_OPEN_READONLY;
    if(mode == rw)
	flag = SQLITE_OPEN_READWRITE;
    else if(mode == crw)
	flag = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;

    cout <<  "Trying to open database " << dbfile << endl;
    if(sqlite3_open_v2(dbfile, &database, flag, NULL) != SQLITE_OK){ 
	cerr << "Could not open database "<< dbfile << endl;
	cerr << error() << endl;
	exit(1);
    }
}

void SQLiteDB::exec(const char *sql){

    Statement statement(this);
    statement.prepare(sql);
    statement.step();
}

void SQLiteDB::beginTransaction(){
    exec("BEGIN TRANSACTION");
}

void SQLiteDB::endTransaction(){
    exec("END TRANSACTION");
}


void SQLiteDB::close(){
    sqlite3_close(database);   
}

bool SQLiteDB::tableExists(string table_name){

    Statement stmt(this);
    stmt.prepare("SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?1;");
    stmt.bindText(1,table_name.c_str());
    if(stmt.nextResult()){
        if (stmt.boolColumn(0)){
	    return true;
	}
    }
    return false;
}

/*
 * create table `genomes` if it does not already exist
 * the offset is the start position in the file
 * the streamsize is the chunk length in the file (e.g. the
 * length of the sequence segment including newlines)
 */
void SQLiteDB::createTableGenomes(){
    const char *sql="CREATE TABLE IF NOT EXISTS genomes ( \
        seqid INTEGER PRIMARY KEY AUTOINCREMENT, \
        seqnr INTEGER NOT NULL, \
        speciesid INTEGER NOT NULL, \
        start INTEGER NOT NULL, \
        end INTEGER NOT NULL, \
        offset UNSIGNED BIG INT NOT NULL, \
        streamsize INTEGER NOT NULL, \
        FOREIGN KEY (speciesid,seqnr) REFERENCES seqnames(speciesid,seqnr) );";
    exec(sql);
}
/*
 * create table 'speciesnames' if it does not already exist
 */
void SQLiteDB::createTableSpeciesnames(){
    const char* sql = "CREATE TABLE IF NOT EXISTS speciesnames (\
        speciesid INTEGER PRIMARY KEY AUTOINCREMENT,\
        speciesname TEXT UNIQUE NOT NULL);";
    exec(sql);
}
/*
 * create table 'seqnames' if it does not already exist     
 */
void SQLiteDB::createTableSeqnames(){
    const char* sql = "CREATE TABLE IF NOT EXISTS seqnames ( \
         seqnr INTEGER PRIMARY KEY AUTOINCREMENT, \
         speciesid INTEGER NOT NULL REFERENCES speciesnames(speciesid), \
         seqname TEXT NOT NULL, \
         UNIQUE(speciesid,seqname));";
    exec(sql);
}

/*
 * create table 'hints' if it does not already exist
 */
void SQLiteDB::createTableHints(){
    const char* sql = "CREATE TABLE IF NOT EXISTS hints (\
         hintid INTEGER PRIMARY KEY AUTOINCREMENT, \
         speciesid INTEGER NOT NULL, \
         seqnr INTEGER NOT NULL, \
         source TEXT, \
         start INTEGER NOT NULL, \
         end INTEGER NOT NULL, \
         score REAL DEFAULT 0.0, \
         type INTEGER CHECK( type >=0 AND type < 17) NOT NULL, \
         strand TEXT CHECK( strand IN ('+','-','.') ) DEFAULT '.', \
         frame TEXT CHECK( frame IN ('0','1','2','.') ) DEFAULT '.', \
         priority INTEGER DEFAULT -1, \
         grp TEXT DEFAULT '', \
         mult INTEGER CHECK( mult >= 1) DEFAULT 1, \
         esource TEXT NOT NULL, \
         FOREIGN KEY (speciesid,seqnr) REFERENCES seqnames(speciesid,seqnr));";
    exec(sql);
}

void SQLiteDB::createTableFeatureTypes(){

    const char* sql = "CREATE TABLE IF NOT EXISTS featuretypes (\
         typeid INTEGER PRIMARY KEY, \
         typename TEXT CHECK( typename IN ('start','stop','ass','dss','tss','tts','exonpart','exon','intronpart','intron', \
         'irpart','CDS','CDSpart','UTR','UTRpart','nonexonpart','genicpart') ) NOT NULL);";
    exec(sql);

    Statement statement(this);
    
    beginTransaction();
    statement.prepare("INSERT OR IGNORE INTO featuretypes VALUES (?1,?2);");
    for(int i = 0; i < NUM_FEATURE_TYPES; i++){
	const char* name = featureTypeNames[i];
	statement.bindInt(1,i);
	statement.bindText(2,name);
	statement.step();
	statement.reset();
    }
    endTransaction();
}

int SQLiteDB::getSpeciesID(string species){

    Statement stmt(this);
    stmt.prepare("SELECT speciesid FROM speciesnames WHERE speciesname=?1;");
    stmt.bindText(1,species.c_str());
    if(stmt.nextResult()){
        int id = stmt.intColumn(0);
        return id;
    }
    else{
        string sql = "INSERT INTO speciesnames (speciesname) VALUES (\"" + species + "\")";
        exec(sql.c_str());
        return lastInsertID();
    }
}



void Statement::prepare(const char *sql){

    if (sqlite3_prepare_v2(database, sql, -1, &stmt, 0) != SQLITE_OK){
        throw error();
    }
}

void Statement::step(){

    if (sqlite3_step(stmt) != SQLITE_DONE){
        throw error();
    }
}

void Statement::bindInt(int idx, int x){

    if(sqlite3_bind_int(stmt, idx, x) != SQLITE_OK){
	throw error();
    }
}

void Statement::bindInt64(int idx, uint64_t x){

    if(sqlite3_bind_int64(stmt, idx, x) != SQLITE_OK){
	throw error();
    }
}

void Statement::bindDouble(int idx, double d){

    if(sqlite3_bind_double(stmt, idx, d) != SQLITE_OK){
	throw error();
    }
}

void Statement::bindText(int idx, const char* text){

    if(sqlite3_bind_text(stmt, idx, text, strlen(text), NULL) != SQLITE_OK){
	throw error();
    }
}

bool Statement::nextResult(){
    int msg = sqlite3_step(stmt);
    if( msg != SQLITE_DONE && msg != SQLITE_ROW )
	throw error();
    return (msg == SQLITE_ROW);
}
