/**
 * /file global.hpp
 * Use stringstream instead of cerr in Test Mode.
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <sstream>

#include "Test.hpp"

#ifdef TEST_MODE
	extern std::stringstream g_output;
	#define ERR_STRM g_output
#else
	#define ERR_STRM cerr
#endif

#endif /* GLOBAL_H_ */
