/**
 * \file flex_vec.hpp
 */

#ifndef FLEX_VEC_HPP_
#define FLEX_VEC_HPP_

#include <vector>

/**
 * @brief Special vector in which the value of its first index can be chosen individually
 */
template<class T> class Flex_Vec : public std::vector<T>
{
public:
	Flex_Vec(int first_pos, int last_pos) : std::vector<T>(last_pos - first_pos + 1), m_first_pos(first_pos) { }

	using std::vector<T>::at;
	using std::vector<T>::operator[];
	using std::vector<T>::size;
	using std::vector<T>::clear;

	// TODO: at -> []
    T& operator[](int i) { return at(i - m_first_pos); }
    const T& operator[](int i) const { return at(i - m_first_pos); }
    T& operator[](unsigned i) { return at(i - m_first_pos); }
    const T& operator[](unsigned i) const { return at(i - m_first_pos); }
private:
    int m_first_pos;
};


#endif /* FLEX_VEC_HPP_ */
