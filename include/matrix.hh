/*
 * matrix.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _MATRIX_HH
#define _MATRIX_HH

// standard C/C++ includes
#include <vector>
#include <ostream>  // for operator<<


using namespace std;

/**
 * @brief A simple matrix class. Base class for all mathematical matrix objects.
 * @details This class provides the functionality of a matrix. It is designed
 *          to easily allocate space without the need to free this explicitily.
 *
 * @author  Emmanouil Stafilarakis
 * @version 0.1
 */
template <class T>
class Matrix {
public:
    /// The type of the stored values
    typedef T         value_t;
public:
    /**
     * Constructor
     */
    Matrix(int n = 0, int m = 0) {
	init(n,m);
    }
    /**
     * Destructor
     */
    ~Matrix( ){ 
	delete[] data;
    }
   /**
     * this is a row: just a pointer
     */
    value_t* operator[] ( int i ) {
	return data + i*width;
    }
    /**
     *
     */
    const value_t* operator[] ( int i ) const {
	return data + i*width;
    }
    /**
     *
     */
    value_t& operator() ( int i, int j ) {
	return (*this)[i][j];
    }
    /**
     *
     */
    const value_t& operator() ( int i, int j ) const {
	return (*this)[i][j];
    }
    /**
     * this clones the row
     */
    vector<T> getRow(int i) const {
	return vector<T>(data + i*width, data + (i+1)*width);
    }
    /**
     *
     */
    vector<T> getColumn(int j) const {
	vector<T> result;
	result.reserve(height);
	for (value_t* p = data + j; p < data + size; p+=width)
	    result.push_back(*p);
	return result;
    }
    /**
     *
     */
    int getColSize( ) const {
	return height;
    }
    /**
     *
     */
    int getRowSize( ) const {
	return width;
    }
    /**
     *
     */
    void operator=(const Matrix<value_t>& mat) {
	resize(mat.height, mat.width);
	for (int i=0; i<size; i++) 
	    data[i] = mat.data[i];
    }
    void assign(int n, int m, value_t t = value_t()) {
	resize(n,m);
	for (int i=0; i<size; i++) 
	    data[i] = t;
    }
    friend ostream& operator<<(ostream& out, const Matrix<value_t>& mat) {
	for (int i = 0; i < mat.size; i++) {
	    out << mat.data[i];
	    if (i+1 % mat.width)
		out << "\t";
	    else
		out << endl;
	}
	return out;
    }

private:
    void resize(int n = 0, int m = 0) {
	if (n == height && m == width) 
	    return;
	delete[] data;
	init(n,m);
    }
    void init(int n = 0, int m = 0) {
	if (n>0 && m>0) {
	    height=n; width=m; size=n*m;
	    data  = new value_t[size];
	} else {
	    size=height=width=0; data=0;
	}
    }
    int     height;  // n
    int     width;  // m
    int     size;    // n*m;
    value_t   *data;  // the nxm matrix
};

#endif   //  _MATRIX_HH

