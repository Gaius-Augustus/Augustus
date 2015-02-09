 /**********************************************************************
 * file:    projectio.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 **********************************************************************/

#ifndef _PROJECTIO_HH
#define _PROJECTIO_HH

#include <vector>
#include <sys/stat.h> // creating directories
#include <sstream>
#include <fstream>

class ProjectError : public std::exception {

public:
    ProjectError(std::string msg) throw() : exception(), message(msg) {}
    ProjectError() throw () : exception() {}
    ~ProjectError() throw() {}    
    std::string getMessage(){return message;}

private:
    std::string message;
};


/*
 * convert int to string
 */
inline std::string itoa(int n) {
    std::ostringstream strm;
    strm << n;
    return strm.str();
}

/*
 * expand the ~ to the $Home directory
 */
inline std::string expandHome(std::string filename){
    if (filename.length()>0 && filename[0]=='~') {
	filename.replace(0,1,getenv("HOME"));
    }
    return filename;
}

/*
 * running external commands
 * returns stdout
 */
inline std::string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
	if(fgets(buffer, 128, pipe) != NULL)
	    result += buffer;
    }
    pclose(pipe);
    return result;
}

/*
 * append slash to the end of a directory name 
 */
inline std::string expandDir(std::string filename){
    if(filename.length()>0 && filename[filename.size()-1] != '/')
	filename += '/';                     // append slash if neccessary 
    filename = expandHome(filename);         // replace "~" by "$HOME
    return filename;
}

/*
 * create directory if not existing
 */
inline void createDir(std::string dir){

    struct stat buffer;
    if( stat(dir.c_str(), &buffer) == -1 || !S_ISDIR(buffer.st_mode)){
        int err = mkdir(dir.c_str(),0755);
        if(err != 0)
	    throw ProjectError("Could not create directory " + dir + "\n");
    }
}

/*
 * split a string into tokens based on a given delimiter
 */
inline std::vector<std::string> splitString(std::string line, char delimiter='\t'){

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while(getline(iss, token, delimiter)){
	if(token.empty())
	    throw ProjectError("Column " + itoa(tokens.size()+1) + " is empty.\n");
	tokens.push_back(token);
    }
    return tokens;
}

#endif   //  _PROJECTIO_HH
