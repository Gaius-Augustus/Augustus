/*
 * projectio.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#include "projectio.hh"
#include <sstream>

istream& operator>>(istream& is, Goto_line_after gla){ 
    return find_line_after(is, gla.keyword);
}

stringstream& operator>>(stringstream& is, Goto_line_after gla){ 
    return find_line_after(is, gla.keyword);
}
