
/*
*   added by Giovanna Migliorelli 24.07.2019 
*   Tokenizer class can handle more than one separator at the time 
*   aknowledgement : https://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
*   todo : make it a template
*/

#ifndef TOKENIZER_HH
#define TOKENIZER_HH

#include <string>
using namespace std;

class Tokenizer 
{
    public:
        static const std::string DELIMITERS;
        Tokenizer(const std::string& str);
        Tokenizer(const std::string& str, const std::string& delimiters);
        bool NextToken();
        bool NextTokenBlank();
        bool NextTokenBlank(const std::string& delimiters); 
        bool NextToken(const std::string& delimiters);
        string GetToken(){return m_token;};
        void Reset();
    protected:
        const std::string m_string;
        size_t m_offset;
        int m_blanks;
        std::string m_token;
        std::string m_delimiters;
};


#endif