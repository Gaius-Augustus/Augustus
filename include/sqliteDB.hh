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

/* database opening flags:
 * ro (read-only mode, throws and error if db does not exist, default flag)  
 * rw (writing only possible if db is not write protected, throws and error if db does not exist) 
 * crw (reading and writing, creates db if it does not exist)  
 */
enum OpenMode{ro,rw,crw};

class SQLiteDB{
public:
    SQLiteDB(const char* f, OpenMode m=ro) : dbfile(f), database(NULL) {
	open(m);
    }
    ~SQLiteDB(){
	close();
    }
    
    void connect(const char* filename);
    void close();
    void exec(const char *sql);
    void open(OpenMode m);
    
    // create tables
    bool tableExists(string table_name);
    void createTableGenomes();
    void createTableSpeciesnames();
    void createTableSeqnames();
    void createTableHints();
    void createTableFeatureTypes();
    int getSpeciesID(string species);

    void beginTransaction();
    void endTransaction();

    inline int lastInsertID(){return (int)sqlite3_last_insert_rowid(database);}
    inline const char* error(){return sqlite3_errmsg(database);}
    inline int numChanges(){return sqlite3_changes(database);}
    
    friend class Statement;

private:
    const char* dbfile;
    sqlite3 *database;
};

class Statement
{
public:
    Statement(SQLiteDB* db) : stmt(NULL), database(db->database) {}
    ~Statement(){finalize();}

    void prepare(const char *sql);
    void step();
    bool nextResult();
    void bindInt(int idx, int x);
    void bindInt64(int idx, uint64_t x);
    void bindDouble(int idx, double d);
    void bindText(int idx, const char* text);
 
    inline void reset(){sqlite3_reset(stmt);}
    inline void finalize(){sqlite3_finalize(stmt);}
    inline int numCols(){return sqlite3_column_count(stmt);}
    inline int intColumn(int colNum){return sqlite3_column_int(stmt,colNum);}
    inline bool boolColumn(int colNum){return sqlite3_column_int(stmt,colNum) != 0;}
    inline uint64_t int64Column(int colNum){return (uint64_t)sqlite3_column_int64(stmt,colNum);}
    inline double doubleColumn(int colNum){return sqlite3_column_double(stmt,colNum);}
    inline char* textColumn(int colNum){return (char*)sqlite3_column_text(stmt,colNum);}
    inline const char* error(){return sqlite3_errmsg(database);}
    
private:
    sqlite3_stmt *stmt;
    sqlite3 *database;
};

#endif

