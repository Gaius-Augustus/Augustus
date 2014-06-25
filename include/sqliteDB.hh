/**********************************************************************
 * file:    sqliteDB.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  wrapper class around the SQLite interface
 * authors: Stefanie Koenig
 *
 *********************************************************************/

#ifndef __SQLITEDB_H__
#define __SQLITEDB_H__

#include <string>
#include <vector>
#include <sqlite3.h>

using namespace std;

//forward declarations
class SQLiteDB;
class Statement;

class SQLiteDB
{
 public:
    SQLiteDB() : database(NULL) {}
    ~SQLiteDB(){}
    
    void connect(const char* filename);
    void close();
    void exec(const char *sql);

    // create tables
    void createTableGenomes();
    void createTableSpeciesnames();
    void createTableSeqnames();
    void createTableHints();
    void createTableFeatureTypes();
    void beginTransaction();
    void endTransaction();

    inline int lastInsertID(){return (int)sqlite3_last_insert_rowid(database);}
    inline const char* error(){return sqlite3_errmsg(database);}
    inline int numChanges(){return sqlite3_changes(database);}
    
    friend class Statement;

private:
    sqlite3 *database;
};

class Statement
{
public:
    Statement(SQLiteDB* db) : stmt(NULL), database(db->database) {}
    ~Statement(){}

    void prepare(const char *sql);
    void step();
    void bindInt(int idx, int x);
    void bindDouble(int idx, double d);
    void bindText(int idx, const char* text);
    void exec(const char *sql);
 
    inline void reset(){sqlite3_reset(stmt);}
    inline void finalize(){sqlite3_finalize(stmt);}
    inline bool nextResult(){ return sqlite3_step(stmt) == SQLITE_ROW; }
    inline int numCols(){return sqlite3_column_count(stmt);}
    inline int intColumn(int colNum){return sqlite3_column_int(stmt,colNum);}
    inline double doubleColumn(int colNum){return sqlite3_column_double(stmt,colNum);}
    inline char* textColumn(int colNum){return (char*)sqlite3_column_text(stmt,colNum);}
    inline const char* error(){return sqlite3_errmsg(database);}
    
private:
    sqlite3_stmt *stmt;
    sqlite3 *database;
};

#endif

