/**********************************************************************
 * file:    sqliteDB.cc
 * license: Artistic License, see file LICENSE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  wrapper class around the SQLite interface
 * authors: Stefanie Koenig
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 25.06.14|Stefanie Koenig| creation of the file
 **********************************************************************/


#include "sqliteDB.hh"

#include <iostream>
#include <stdlib.h>
#include <string.h>

void SQLiteDB::open(OpenMode mode){
    
    int flag = SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX;
    if(mode == rw)
	flag = SQLITE_OPEN_READWRITE | SQLITE_OPEN_NOMUTEX;
    else if(mode == crw)
	flag = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_NOMUTEX;

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
