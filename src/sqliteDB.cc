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

void SQLiteDB::connect(const char* filename){

    cout <<  "Trying to open database " << filename << endl;
    if(sqlite3_open(filename, &database) != SQLITE_OK){ 
	throw error();   
    }
}

void SQLiteDB::exec(const char *sql){

    Statement statement(this);
    statement.prepare(sql);
    statement.step();
    statement.finalize();
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
         mult INTEGER DEFAULT 1, \
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

void Statement::exec(const char *sql){

    prepare(sql);
    step();
    finalize();
}
