/**
 * \file Test.hpp
 * Contains the test environment code and all the test cases.
 */

#ifndef TEST_HPP_
#define TEST_HPP_

//#define TEST_MODE   ///< Flag. Flag set: Test.cpp is executed. Flag not set: main.cpp is executed.

#ifdef TEST_MODE
	#define TEST_CASES ///< Flag. When set, the normal test cases are executed.
	#define TEST_ERRORS ///< Flag. When set, the error message test cases are executed.
	#define TEST_BUGS ///< Flag. When set, the test cases for bug testing are executed.
//	#define TEST_HUMAN19 ///< Flag. When set, the test case with the human chr19 is run.
#endif

#ifdef TEST_HUMAN19
	#define TEST_WITHOUT_Z ///< Flag. When set, the human chromosome 19 is run without the zero coverage UTR detection.
//	#define TEST_WITH_Z ///< Flag. When set, the human chromosome 19 is run with the zero coverage UTR detection.
//	#define TEST_SPLICE_SITE_FILTER ///< Flag. When set, the test case for the splice site filtering is run.
#endif

#endif /* TEST_HPP_ */
