/*****************************************************************************\
 * Filename : lldouble.hh
 * Author   : Emmanouil Stafilarakis
 * Project  : HMM
 * Version  : 0.1
 *
 * authors: Emmanouil Stafilarakis, Mario Stanke, mario@gobics.de
 *
 * Description: This class implements a double object with a very large
 *              range. It is designed to handle very small (or high) floating
 *              point numbers that would otherwise become zero when multiplied
 *              to each other.
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 17.01.2002 | E. Stafilarakis       | Creation of the file.
 * 06.11.2002 | Mario Stanke          | simplify
 * 19.4.2006  | Mario Stanke          | root
 * 20.9.2007  | Oliver Keller         | partial rewrite
 * 26.7.2008  | Oliver Keller         | exponential
 \*****************************************************************************/

#ifndef _LL_DOUBLE_HH
#define _LL_DOUBLE_HH

// standard C/C++ includes
#include <cmath>
#include <sstream>
#ifdef DEBUG
#include <iostream>
#endif

using namespace std;
typedef ios_base::fmtflags fmtflags;

class LLDouble{
    typedef int exponent_type;

    /* class constants follow
     *
     * NOTE:
     * This procedure can lead to problems when LLDoubles are initialized
     * BEFORE the class constants, giving undefined results running testPrecision
     * 
     * Current solution: ensure that all object files initializing LLDoubles are
     *                   mentioned before lldouble.o when calling the linker
     */

    static const double dbl_inf;

    static const double max_val; // = 2^500
    static const double min_val; // = 2^(-500)

    static const double base;    // = 2^1000
    static const double baseinv; // = 2^(-1000)
    static const double logbase; // = log(base) = 693.15

    static const exponent_type max_exponent; 
    static const exponent_type min_exponent;

    static unsigned temperature; // for "heating", heat = (8-temperature)/8, will later often need to compute pow(d,heat) for LLDoubles d
    static double rest[7]; // precomputed values for heating
 public:
    LLDouble(float x=0.0) : value(x), exponent(0) {} // called when no argument is provided
    LLDouble(double d) : value(d), exponent(0) {
	testPrecision();
    }
    LLDouble(long double d);
    LLDouble(int i) : value((double)i), exponent(0) {}
    LLDouble(long i) : value((double)i), exponent(0) {}

    /*
     * conversion to other types
     */
    long double doubleValue() {
	return (long double)value * std::exp((long double)(exponent) * logbase); 
    }
    string toString(int precision=output_precision, 
 		    fmtflags flags=ios::dec) const;
 
    /*
     * arithmetic operators
     */
    LLDouble& operator+=(const LLDouble& other);
    LLDouble& operator-=(const LLDouble& other) {
	return operator+=(-other);
    }
    LLDouble& operator*=( const LLDouble& other ){
	value *= other.value;
	exponent += other.exponent;
	testPrecision();
	return *this;
    }
    LLDouble& operator/=(const LLDouble& other){
	value /= other.value;
	exponent -= other.exponent;
	testPrecision();
	return *this;
    }

    LLDouble operator+( const LLDouble& other ) const {
	return LLDouble(*this) += other;
    }
    LLDouble operator-( const LLDouble& other ) const {
	return LLDouble(*this) -= other;
    }
    LLDouble operator*( const LLDouble& other ) const {
	return LLDouble(*this) *= other;
    }
    LLDouble operator/( const LLDouble& other ) const {
	return LLDouble(*this) /= other;
    }

    friend LLDouble operator-( const LLDouble& dbl ) {
	return LLDouble(-dbl.value, dbl.exponent);
    }
    LLDouble abs() const {
	return LLDouble(std::abs(value), exponent);
    }
    friend LLDouble abs( const LLDouble& dbl ) {
	return dbl.abs();
    }

    /*
     * comparative operators 
     */
    bool operator==(const LLDouble& other) const;
    bool operator>(const LLDouble& other) const;
    bool operator!=(const LLDouble& other) const {
	return !(*this == other);
    }
    bool operator<(const LLDouble& other) const {
	return other > (*this);
    }
    bool operator<=(const LLDouble& other) const {
	return !((*this) > other);
    }
    bool operator>=(const LLDouble& other) const {
	return !(other > (*this));
    }

    /*
     * root and exponential functions
     */
    LLDouble pow(double x) const;
    LLDouble getRoot(int r) const {
	if (value < 0 && r%2)
	    return -pow(-*this,1.0/r);
	return pow(1.0/r);
    }
    double log() const { 
	return std::log(value) + exponent*logbase; 
    }
    double log(int otherbase) const { 
        return log()/std::log((double) otherbase);
    }
    friend double log(const LLDouble& lld) {
	return lld.log();
    }
    friend double log(int otherbase, const LLDouble& lld) {
	return lld.log(otherbase);
    }
    static LLDouble exp(double x);
    static LLDouble pow(const LLDouble& lld, double x) {
	return lld.pow(x);
    }
    LLDouble heated();

    /*
     * I/O stream operators 
     */
    friend istream& operator>>( istream& in, LLDouble& lld ){
	lld.read( in );
	return in;
    }
    friend ostream& operator<<( ostream& out, const LLDouble& lld ){
	int precision = output_precision > 0 ? output_precision : out.precision();
	return out << lld.toString(precision, out.flags());
    }
    
    /*
     * class functions
     */
    static LLDouble getMaxDouble() {
	return LLDouble(max_val, max_exponent);
    }
    static LLDouble getMinDouble() {
	return LLDouble(min_val, min_exponent);
    }
    static void setOutputPrecision(int p){
	output_precision = p;
    };
    static int getOutputPrecision(){
	return output_precision;
    };
    static LLDouble infinity() {
	return LLDouble(dbl_inf, max_exponent);
    }
    static void setTemperature(unsigned t);

private:
    // for internal use: directly set the data fields
    LLDouble(double v, exponent_type e) : 
	value(v), exponent(e) {}
//     void print( ostream& out ) const;
    void read( istream& in );
    void testPrecision( ) {
	// value is 0, or NaN: keep and set exponent=0
	if (value == 0.0 || std::isnan((double) value)) {
	    exponent = 0;
	    return;
	} // value is infinity: set exponent = max_exponent
	else if (std::abs(value) == dbl_inf) {
	    exponent = max_exponent;
	    return;
	}
	// value is too small 
	while( std::abs(value) < min_val) {
	    if (exponent == min_exponent) {
		value = 0; exponent = 0; return;
	    }
	    value *= base;
	    exponent--;
	}
	// value is too large
	while( std::abs(value) > max_val) {
	    if (exponent >= max_exponent) {
		value = value>0? dbl_inf : -dbl_inf; 
		exponent = max_exponent; return;
	    }
	    value *= baseinv;
	    exponent++;
	}
    }
    static int output_precision;
    double       value; // long double : 40% more time, 32% more memory than double, probably no difference
    exponent_type         exponent;
};

/*
 * arithmetic operators for double and LLDouble
 */
inline LLDouble operator/(long double i, const LLDouble& lld ) {
    return LLDouble(i)/lld;
}

inline LLDouble operator*(long double i, const LLDouble& lld) {
    return LLDouble(i)*lld;
}

inline LLDouble operator+(long double i, const LLDouble& lld) {
    return LLDouble(i)+lld;
}

inline LLDouble operator-(long double i, const LLDouble& lld) {
    return LLDouble(i)-lld;
}

#ifdef DEBUG
inline LLDouble relative_error(const LLDouble& d1, const LLDouble& d2) {
    return abs((d1/d2).doubleValue()-1);
}

inline bool relerror_lessthan(const LLDouble& d1, const LLDouble& d2, double rel_error) {
    if (relative_error(d1, d2) >= rel_error) {
	cerr << "relative error: " << relative_error(d1, d2) << "\n";
	return false;
    }
    return true;
}

inline bool almost_equal(const LLDouble& d1, const LLDouble& d2) {
    return relerror_lessthan(d1,d2,0.01);
}
#endif

/*
 * class LogDouble
 *
 * internally stores floating point numbers using their logarithm
 * this saves time when multiplication and division is a frequent operation
 */

class LogDouble{
public:
    LogDouble( double d=0.0 );
    LogDouble( const LogDouble& other ){
	logvalue    = other.logvalue;
    }

    LogDouble operator*( const LogDouble& other ) const;
    LogDouble& operator*=( const LogDouble& other );

    // Assignment operator
    LogDouble& operator=( const LogDouble& other ){
	logvalue    = other.logvalue;
	return *this;
    }
    void print( ostream& out ) const;
private:
    static int outputprecision;
    double  logvalue;
};

ostream& operator<<( ostream& out, const LogDouble& logd );

#endif   //  _LL_DOUBLE_HH
