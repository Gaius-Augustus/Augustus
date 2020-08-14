
/*
*   added by Giovanna Migliorelli 24.07.2019 
*   Tokenizer class can handle more than one separator at the time 
*   aknowledgement : https://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
*   todo : make it a template
*/

#include "tokenizer.hh"

const std::string Tokenizer::DELIMITERS(" \t\n\r");

Tokenizer::Tokenizer(const std::string& s) :
    m_string(s), 
    m_offset(0), 
    m_blanks(0),
    m_delimiters(DELIMITERS) {}

Tokenizer::Tokenizer(const std::string& s, const std::string& delimiters) :
    m_string(s), 
    m_offset(0),
    m_blanks(0), 
    m_delimiters(delimiters) {}

bool Tokenizer::NextToken() 
{
    return NextToken(m_delimiters);
}

bool Tokenizer::NextTokenBlank() 
{
    return NextTokenBlank(m_delimiters);
}

bool Tokenizer::NextToken(const std::string& delimiters) 
{
    size_t i = m_string.find_first_not_of(delimiters, m_offset);
    if (std::string::npos == i) 
    {
        m_offset = m_string.length();
        return false;
    }

    size_t j = m_string.find_first_of(delimiters, i);
    if (std::string::npos == j) 
    {
        m_token = m_string.substr(i);
        m_offset = m_string.length();
        return true;
    }

    m_token = m_string.substr(i, j - i);
    m_offset = j;
    return true;
}

// returns an empty string if one delimiters comes close to the next
bool Tokenizer::NextTokenBlank(const std::string& delimiters) 
{
    if(m_blanks>0){
        m_token = "";
        --m_blanks;  
        ++m_offset;
        return true;  
    }
    
    size_t i = m_string.find_first_not_of(delimiters, m_offset);

    if(std::string::npos != i && std::string::npos != m_offset){
        if(i>m_offset){
            m_blanks = i - m_offset;

            if(m_offset>0)
                --m_blanks;
        }
    }
    else if(std::string::npos != m_offset){
        m_blanks = m_string.size() - m_offset;

        if(m_offset==0)
            ++m_blanks;
    }

    if(m_blanks>0){
        m_token = "";
        --m_blanks;  
        ++m_offset;
        return true;  
    }
    
    if (std::string::npos == i) 
    {
        m_offset = m_string.length();
        return false;
    }

    size_t j = m_string.find_first_of(delimiters, i);
    if (std::string::npos == j) 
    {
        m_token = m_string.substr(i);
        m_offset = m_string.length();
        return true;
    }

    m_token = m_string.substr(i, j - i);
    m_offset = j;
    return true;
}
