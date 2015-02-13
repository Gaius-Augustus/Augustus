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

/*
 * accumulate all values in a vector in back-to-front order
 */
inline std::vector<int> reverse_accumulate(std::vector<int> v){
    
    int cumValue = 0;
    for(int i=v.size()-1; i >=0; i--){
	cumValue += v[i];
	v[i] = cumValue;
    }
    return v;
}


/*
 * sum elements in two vectors
 */
inline std::vector<int> sum(std::vector<int> u, std::vector<int> v){
    if(u.size() != v.size())
	throw ProjectError("Internal error: sum elements in two vectors of different size.\n");
    for(int i=0; i<u.size(); i++){
	u[i] += v[i];
    }
    return u;
}

/*
 * print a pictograph from values in a vector
 */
inline void writePictograph(std::vector<int> v, std::ofstream &of){

    v = reverse_accumulate(v);
    int sum = v.front();
    
    of.setf(std::ios::showpoint|std::ios::fixed, std::ios::floatfield);
    of.precision(1);

    for(int i=0; i<v.size(); i++){
	double per = (sum>0)? (((double)v[i]*100)/sum) : 0.0;
	double numStars = (sum>0)? ((v[i]*25)/sum) : 0.0;
        of << "#";
	of.width(4); of << i;
	of.width(10); of << v[i];
	of.width(6); of << per << "% ";
        of << std::string(numStars, '*') << std::endl; // pictograph
    }
}

/*
 * print a pictograph from the sum of values of two vectors
 */
inline void writePictograph(std::vector<int> u, std::vector<int> v, std::ofstream &of){

    u = reverse_accumulate(u);
    v = reverse_accumulate(v);
    int sumV = v.front();
    int sumU = u.front();
    std::vector<int> w = sum(u,v);
    
    
    of.setf(std::ios::showpoint|std::ios::fixed, std::ios::floatfield);
    of.precision(1);

    for(int i=0; i<v.size(); i++){
	double percV = (sumV>0)? (((double)v[i]*100)/sumV) : 0.0;
	double percU = (sumU>0)? (((double)u[i]*100)/sumU) : 0.0;
	double percW = ((sumV+sumU)>0)? (((double)w[i]*100)/(sumU+sumV)) : 0.0;
	double numStars = ((sumV+sumU)>0)? ((w[i]*25)/(sumV+sumU)) : 0.0;

        of << "#";
	of.width(4); of << i;
	of.width(10); of << u[i];
	of.width(6); of << percU << "%";
	of.width(10); of << v[i];
	of.width(6); of << percV << "%";
	of.width(10); of << w[i];
	of.width(6); of << percW << "% ";

        of << std::string(numStars, '*') << std::endl; // pictograph
    }
}

#endif   //  _PROJECTIO_HH
