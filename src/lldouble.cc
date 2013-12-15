/**********************************************************************
 * file:    lldouble.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  floating point number with very high precision
 * authors: Emmanouil Stafilarakis, Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 17.01.02| Stafilarakis  | Creation of the file
 * 2.5.02  | Mario Stanke  | debugging
 * 19.4.06 | Mario Stanke  | root
 * 19.4.06 | Mario Stanke  | root
 * 26.7.08 | Oliver Keller | rewrite; exponential
 **********************************************************************/

#include "lldouble.hh"

#include <iomanip>
#include <limits>
#include <iostream>

/* =====[ LLDouble ]======================================================== */

// constants used only internally
static const double high_exponent = numeric_limits<double>::max_exponent - 24;  // = 1000
static const double log_10 = std::log(10);

// exported constants
const double LLDouble::dbl_inf = numeric_limits<double>::infinity();

// Range for double inside LLDouble
const double LLDouble::max_val = ::pow(2.0, high_exponent/2);  // = 3.2733906078961419e+150
const double LLDouble::min_val = 1.0 / LLDouble::max_val;    // = 3.0549363634996047e-151

// scope of range
const double LLDouble::base = std::pow(2.0, high_exponent);       // = 1.0715086071862673e+301
const double LLDouble::baseinv = 1.0 / LLDouble::base;       // = 9.3326361850321888e-302
const double LLDouble::logbase = std::log(LLDouble::base);   // = 693.14718055994535

const LLDouble::exponent_type LLDouble::max_exponent = numeric_limits<exponent_type>::max();
const LLDouble::exponent_type LLDouble::min_exponent = numeric_limits<exponent_type>::min();
int LLDouble::output_precision = 0;
unsigned LLDouble::temperature = 0;
double LLDouble::rest[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};

/* 
 * =====[ constructors ]=======================================================
 */
LLDouble::LLDouble(long double d) : exponent(0) {
    if (d != 0.0 
	&& std::abs(d) != numeric_limits<long double>::infinity() 
	&& d == d)
    {
	while (std::abs(d) > max_val)  {
	    exponent++; d*=baseinv;
	}
	while (std::abs(d) < min_val) {
	    exponent--; d*=base;
	}
    }
    value = (double) d;
}


/*
 * =====[ string conversion ]===================================================
 */
// print a large float in fixed mode
// exponent must be between 18 and 1000, factor between 1 and 10
inline void print_fixed(ostream& strm, double factor, int exponent) {
    int precision = strm.precision(0);
    strm << factor * 1e17 
	<< setfill('0') << setw(exponent-18) << ""
	<< setprecision(precision) << 0.0;
}
// print in scientific mode
inline void print_scientific(ostream& strm, double factor, double exponent) {
    strm << factor
	 << ((strm.flags() & ios::uppercase) ? 'E' : 'e') 
	 << showpos << fixed << setprecision(0) << exponent;
}

string LLDouble::toString(int precision, fmtflags flags) const {
    ostringstream estrm("");
    estrm.flags(flags);
    estrm.precision(precision);
    if (exponent == 0) 
	estrm << value;
    else {
	long double logvalue = std::log(std::abs(value)) +
	    (long double)(exponent) * logbase;
	double outexp = floor(logvalue/log_10);
	double outval = std::exp(logvalue - outexp * log_10);
	if (outval>=9.5) {
	    outval/=10; outexp++;
	} 
	if (value<0) outval *= -1;
	// = value * pow(10.0, exponent*logbase/log_10 - outexp);
	if (flags & ios::scientific) 
	    estrm << fixed;
	if ((flags & ios::fixed) && outexp < 1000) 
	    if (exponent < 0) 
		estrm << 0.0;
	    else 
		print_fixed(estrm, outval, (int)outexp);
	else
	    print_scientific(estrm, outval, outexp);
    }
    return estrm.str();
}


/*
 * operator +=
 */
LLDouble& LLDouble::operator+=( const LLDouble& other ) {
    if (other.value == 0.0)
	return *this;
    if (value == 0.0 || (other.exponent > exponent + 1)) 
	return *this = other; 
    if  (other.exponent < exponent - 1)
	return *this;
    if ( exponent == other.exponent ) {
	value += other.value;
	testPrecision();  // only case where range should be tested
    } else if (exponent > other.exponent) {
	value += other.value * baseinv; 
    } else {
	value = value * baseinv + other.value;
	exponent++; 
    }
    return (*this);
}


/*
 * =====[ comparative operators ]==============================================
 */
bool LLDouble::operator==(const LLDouble& other) const {
    if (exponent == other.exponent)
	return value==other.value;
    else if (exponent == other.exponent +1)
	return (value == other.value * baseinv);
    else 
	return 
	    exponent == other.exponent -1 &&
	    value == other.value * base;
}

bool LLDouble::operator>(const LLDouble& other) const {
    if ((value >= 0.0 && other.value <= 0.0) ||
	(value <= 0.0 && other.value >= 0.0) ||
	(exponent == other.exponent) )
	return value>other.value;
    int delta = exponent - other.exponent;
    if (delta == -1)
	return value > other.value * base;
    if (delta == 1)
	return value * base > other.value;
    if (value > 0.0)
	return delta > 0;
    else
	return delta < 0;
}


/*
 * =====[ root and exponential functions ]======================================
 */


/*
 * power with LLDouble base
 */
LLDouble LLDouble::pow(double x) const {
    if (x==0) // 0^0 = 1
	return 1;
    if (value == 0)
	return x>0 ? 0 : infinity();
    if (x == 1.0)
	return *this;
    if (value < 0 && // negative base
	(x - 0.25) < x && x < (x + 0.25) && 
	x == floor(x) )  // integer exponent
	return x/2 == floor(x/2) ? // even exponent
	    exp(x * abs().log()) :
	    - exp(x * abs().log());
   return exp(x * log());
}

/*
 * exponential
 */
LLDouble LLDouble::exp(double x) {
    double set_exponent = x / logbase;
    if (set_exponent > max_exponent) 
	return infinity();
    if (set_exponent < min_exponent)
	return 0;
    LLDouble result(std::exp(x - logbase * exponent_type(set_exponent)),
		    exponent_type(set_exponent));
    result.testPrecision();
    return result;
}

void LLDouble::setTemperature(unsigned t){
    temperature = t;
    // precompute base^(1/8), base^(2/8), ..., base^(7/8), so heating goes faster later
    for (size_t i=0; i<7; i++)
	LLDouble::rest[i] = std::pow(base, (double) (i+1)/8);
}

// apply a power function, more efficient through above precomputation and by using exponents that are a fraction of 8
// (3 'std::sqrt' calls take about 52% of the time of one 'std::pow' call)
LLDouble LLDouble::heated(){
    if (temperature == 0 || value == 0.0)
	return *this; // cold: raise to the power of 1
    double heat = (8.0 - temperature) / 8;
    // this is the inefficient but equivalent way:
    //return pow(heat);
    exponent_type newexponent = (exponent_type) (heat * exponent);
    int i = (8 - temperature) * exponent - 8 * newexponent;
    if (exponent < 0 && i != 0){
	newexponent--; // as C++ rounds down for positive numbers and up for negative numbers
	i += 8;
    }
    double r = 1.0;
    if (i)
	r = rest[i-1];
    LLDouble d(r, newexponent);
    d.testPrecision();
    // d needs to be multiplied with pow(value, heat)
    switch (temperature) {
    case 1:
	d *= value; // do this in two steps, as the *= operator implicitly does testPrecision() to prevent underflow
	d /= sqrt(sqrt(sqrt(value)));
	break;
    case 2:
	d *= value;
	d /= sqrt(sqrt(value));
	break;
    case 3:
	{
	    // x^5/8 = sqrt(x) * sqrt(sqrt(sqrt(x)))
	    double x = sqrt(value);
	    d *= x;
	    d *= sqrt(sqrt(x));
	}
	break;
    case 4:
	d *= sqrt(value);
	break;
    case 5:
	{
	    double x = sqrt(sqrt(value));
	    d *= x;
	    d *= sqrt(x);
	}
	break;
    case 6:
	d *= sqrt(sqrt(value));
	break;
    case 7:
	d *= sqrt(sqrt(sqrt(value)));
	break;
    default: // not used
	d *= std::pow(value, heat);
    }
    return d;
}

/*---------------------------------------------------------------------------*\
 |                                                                           |
 |  private methods                                                          |
 |                                                                           |
\*---------------------------------------------------------------------------*/


/*
 * =====[ read ]===============================================================
 */
void LLDouble::read( istream& in ){
    // to do: enable higher exponents
    if (in >> value) {
	exponent = 0;
	testPrecision();
    }
}


/* =====[ LogDouble ]======================================================= */

int LogDouble::outputprecision = 3;

LogDouble::LogDouble( double d ){
    if (d==0){
	logvalue = - numeric_limits<double>::infinity();
    } else {
	logvalue = log10(d);
    }
}

ostream& operator<<( ostream& out, const LogDouble& logd ){
    logd.print( out );
    return out;
}

void LogDouble::print( ostream& out ) const{
    out.precision(outputprecision);
    out << logvalue << " ";
    if (logvalue > -3 && logvalue < 3)
	out << pow(10.0, logvalue);
    else {
	int main = (int) logvalue;
	double rest = logvalue-main;
	out << pow(10.0, rest) << "e" << main;
    }
}

LogDouble LogDouble::operator*( const LogDouble& other ) const{
    LogDouble newdouble;
    newdouble.logvalue = logvalue + other.logvalue;
    return newdouble;
}

LogDouble& LogDouble::operator*=( const LogDouble& other ){
    logvalue += other.logvalue;
    return *this;
}

/* ========================================================================= */

