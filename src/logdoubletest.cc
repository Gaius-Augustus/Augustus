#include "lldouble.hh"

#include <iostream>

int main() {
    LogDouble ld, lf=0.9999;
    LLDouble d, f=0.9999;
    ld = 0.25;
    d = 0.25;
    clock_t anfang, ende;
    cout << "sizeof(LLDouble) " << sizeof(LLDouble) << endl;
    cout << "sizeof(LogDouble) " << sizeof(LogDouble) << endl;
    cout << "sizeof(double) " << sizeof(double) << endl;
    cout << "sizeof(long double) " << sizeof(long double) << endl;
  
 
    anfang = clock();
    for (int i=0; i<100000000; i++) {
	d = d * f;
	ld = ld * lf;
    }
    ende = clock();
    cout << "both time " << ende - anfang << endl;

    anfang = clock();
    for (int i=0; i<100000000; i++) {
	d = d * f;
    }
    ende = clock();
    cout << "LLDouble " << d << " time " << ende - anfang << endl;

    anfang = clock();
    for (int i=0; i<100000000; i++) {
	ld = ld * lf;
    }
    ende = clock();
    cout << "LogDouble " << ld << " time " << ende - anfang << endl;


    return 0;
}
// saves about 1/3 of the time
