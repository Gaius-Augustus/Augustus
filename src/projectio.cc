/**********************************************************************
 * file:    projectio.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stafilarakis, Mario Stanke (mario@gobics.de)
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 *         | Stafilarakis  | creation of the class
 * 13.01.03| Mario Stanke  | made goto_line_after conform with g++ 3.2
 **********************************************************************/

#include "projectio.hh"
#include <sstream>

istream& operator>>(istream& is, Goto_line_after gla){ 
    return find_line_after(is, gla.keyword);
}

stringstream& operator>>(stringstream& is, Goto_line_after gla){ 
    return find_line_after(is, gla.keyword);
}
