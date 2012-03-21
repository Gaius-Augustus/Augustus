#ifndef _SCANNER_H_
#define _SCANNER_H_

#if ! defined(_SKIP_YYFLEXLEXER_) && ! defined(_SYSINC_FLEXLEXER_H_)
#include <FlexLexer.h>
#define _SYSINC_FLEXLEXER_H_
#endif

#include "../parser/parserbase.h"

class Scanner: public yyFlexLexer
{


 public:
  Parser::STYPE__ *d_val;  //stores semantic values of a token

  // constructor accepts the input and output streams, default stdin and stdout
  Scanner(Parser::STYPE__ *val, std::istream * in = 0, std::ostream * out = 0):
	 yyFlexLexer(in, out),
	 d_val(val)
  { }
  int yylex();

};

#endif



