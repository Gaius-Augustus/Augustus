/*****************************************************************************\
 * Filename : projectio.hh
 * Author   : Emmanouil Stafilarakis
 * Project  : HMM
 * Version  : 0.1
 *
 * Copyright: ©Stafilarakis
 *
 * Description: Some manipulators for streams.
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 11.09.2001 | Stafilarakis Emm.     | Creation of the file
 * 03.09.2002 | Mario Stanke          | change line length to 8192 (before 256)
\******************************************************************************/

#ifndef _PROJECTIO_HH
#define _PROJECTIO_HH

// standard C/C++ includes
#include <istream>
#include <vector>
#include <cstring>


#define MAX_ROW_LEN 8192
using namespace std;

/*
 * istream& comment( istream& strm, char c )
 *--------------------------------------------------------------
 * This function should not be used directly, use the function
 * defined below ('imanip<char> comment( char c )') instead!!
 */
inline istream& comment_c( istream& strm, char c ){
    char buff[MAX_ROW_LEN];
    while( (strm >> ws) && (strm.peek() == c) ){
        strm.getline( buff, MAX_ROW_LEN );
    }
    return strm;
}

/*
 * imanip<char> comment( char c )
 *--------------------------------------------------------------
 * This function removes all comment lines from a stream until
 * a line which not begins with a given comment sign "c".
 * USAGE:
 *        istrm >> comment('CHAR') >> ... ;
 */
/*inline imanip comment( char c ){
    return imanip( comment_c, c );
    }*/


/*
 * istream& comment( istream& strm )
 *--------------------------------------------------------------
 * This function removes all comment lines from a stream until
 * a line which not begins with a comment sign '#'.
 * USAGE:
 *        istrm >> comment >> ... ;
 */
inline istream& comment(istream& strm ){
    return comment_c( strm, '#' );
}



inline istream& find_line_after( istream& strm, const char* str ){
    char buff[MAX_ROW_LEN];

    do{
        strm >> ws;
        strm.getline( buff, MAX_ROW_LEN-1 );
        if ( strncmp( buff, str, strlen(str) ) == 0 )
            return strm;
    } while( strm );
    strm.clear( ios::failbit );
    return strm;
}

struct Goto_line_after { const char *keyword; };
inline Goto_line_after goto_line_after(const char *str){
    Goto_line_after gla;
    gla.keyword = str;
    return gla; 
}

istream& operator>>(istream& is, Goto_line_after gla);

/*
inline imanip<const char*> goto_line_after( const char* str ){
    return imanip<const char*>( find_line_after, str );
}
*/

/*
 *
 */
template <class T>
ostream &operator<<(ostream &output, const vector<T> &v) {
    output << "[";
    if (v.size()>0) {
	output << v[0];
	for(int i=1; i<v.size(); i++) {
	    output << " , " << v[i];
	}
    }
    output << "]" << endl;
    return output;
}


#endif   //  _PROJECTIO_HH
