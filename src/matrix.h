#pragma once

#include <cassert>
#include <vector>

namespace tmd {
// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain
#define TMD_MATRIX_DEFINE_OPERATORS \
	const T& operator()(int i) const { return m_data[i]; } \
	      T& operator()(int i)       { return m_data[i]; } \
	const T& operator()(int i, int j) const { return m_data[index(i, j)]; } \
	      T& operator()(int i, int j)       { return m_data[index(i, j)]; }

////////////////////////////////////////////////////////////////////////
template<typename T>
class Matrix {
	std::vector<T> m_data;
	int m_i, m_j;
public:
	// column-major
	Matrix() : m_i(0), m_j(0) {}
	Matrix(int i, int j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {
		#ifdef DEBUG
		assert(m_j >=0);
		assert(m_i >=0);
		#endif
	}
	int index(int i, int j) const {//return the index in m_data, column-major
		#ifdef DEBUG
		assert(j < m_j && j >=0);
		assert(i < m_i && i >=0);
		#endif
		return i + m_i*j; //column-major
	}
	void resize(int m, int n, const T& filler_val) {//resize to mxn, new sizes should be the same or greater than the old, preserves original data
		#ifdef DEBUG
		assert(m >= dim_1());//new sizes should be the same or greater than the old
		assert(n >= dim_2());
		#endif
		if(m == dim_1() && n == dim_2()) return; // no-op
		std::vector<T> tmp(m*n, filler_val);
		for(int i = 0; i < m_i; i++) {
			for(int j = 0; j < m_j; j++) {
				tmp[i+m*j] = (*this)(i, j);
			}
		}
		m_data = tmp;
		m_i = m;
		m_j = n;
	}
	void append(const Matrix<T>& x, const T& filler_val) {//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
		int m = dim_1();
		int n = dim_2();
		resize(m + x.dim_1(), n + x.dim_2(), filler_val);
		for(int i = 0; i < x.dim_1(); i++) {
			for(int j = 0; j < x.dim_2(); j++) {
				(*this)(i+m, j+n) = x(i, j);
			}
		}
	}
	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	int dim_1() const { return m_i; }
	int dim_2() const { return m_j; }
};
//////////////////////////////////////////////////////////////////////////




// inline int triangular_matrix_index(int n, int i, int j) {
// 	assert(j < n);
// 	assert(i <= j);
// 	return i + j*(j+1)/2;
// }
// //permissive
// inline int triangular_matrix_index_permissive(int n, int i, int j) {
// 	return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
// }
//only upper-right triangular{i,j}
// 0  1  3  6
//    2  4  7
//       5  8
//          9
template<typename T>
class Triangular_Matrix {
	std::vector<T> m_data;
	int m_dim;
public:
	Triangular_Matrix() : m_dim(0) {}
	Triangular_Matrix(int n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {
		#ifdef DEBUG
		assert(n > 0);
		#endif
	}
	int index(int i, int j) const {//return corresponding index in vector m_data
		#ifdef DEBUG
		assert(i >= 0 && j >= 0);
		assert(j < m_dim);
		assert(i <= j);
		#endif
		return i + j*(j+1)/2;
	}
	int index_permissive(int i, int j) const {
		return (i < j) ? index(i, j) : index(j, i);
	}
	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	int dim() const {//return m_dim, dimension of the Matrix
		return m_dim;
	}
};
///////////////////////////////////////////////////////////////////////////





//no diagonal elements
template<typename T>
class Strictly_Triangular_Matrix {
	std::vector<T> m_data;
	int m_dim;
public:
	Strictly_Triangular_Matrix() : m_dim(0) {}
	Strictly_Triangular_Matrix(int n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {
		#ifdef DEBUG
		assert(n > 0);
		#endif
	}
	int index(int i, int j) const {//return corresponding index in vector m_data
		#ifdef DEBUG
		assert(i >= 0 && j > 0);
		assert(j < m_dim);
		assert(i < j);
		assert(j >= 1); // by implication, really
		#endif
		return i + j*(j-1)/2;
	}
	int index_permissive(int i, int j) const {
		return (i < j) ? index(i, j) : index(j, i);
	}
	void resize(int n, const T& filler_val) {//resize m_data to has only the upper right Matrix elements, new sizes should be the same or greater than the old, preserves original data
		#ifdef DEBUG
		assert(n >= m_dim);
		#endif
		if(n == m_dim) return; // no-op
		m_dim = n;
		m_data.resize(n*(n-1)/2, filler_val); // preserves original data
	}
	void append(const Strictly_Triangular_Matrix<T>& m, const T& filler_val) {//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
		int n = dim();
		resize(n + m.dim(), filler_val);
		for(int i = 0; i < m.dim(); ++i) {
			for(int j = i+1; j < m.dim(); ++j) {
				(*this)(i+n, j+n) = m(i, j);
			}
		}
	}
	void append(const Matrix<T>& rectangular, const Strictly_Triangular_Matrix<T>& triangular) {//arguments rectangular and triangular must have the same dim
		//if(rectangular.dim_2() == 0) do nothing
		//if(rectangular.dim_1() == 0) assign argument triangular to this object
		//else append these three matrices like the following
		// i  i  i  r  r  r  r
		//    i  i  r  r  r  r
		//       i  r  r  r  r
		//          t  t  t  t
		//             t  t  t
		//                t  t
		//                   t
		//where i,r,t are this, rectangular and triangular, respectively.
		#ifdef DEBUG
		assert(dim() == rectangular.dim_1());
		assert(rectangular.dim_2() == triangular.dim());
		#endif
		// a filler value is needed by append or resize
		// we will use a value from rectangular as the filler value
		// but it can not be obtained if dim_1 or dim_2 is 0
		// these cases have to be considered separately
		if(rectangular.dim_2() == 0) return;
		if(rectangular.dim_1() == 0) {
			(*this) = triangular;
			return;
		}
		const T& filler_val = rectangular(0, 0); // needed by 'append below'
		int n = dim();
		append(triangular, filler_val);
		for(int i = 0; i < rectangular.dim_1(); ++i) {
			for(int j = 0; j < rectangular.dim_2(); ++j) {
				(*this)(i, n + j) = rectangular(i, j);
			}
		}
	}
	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	int dim() const {
		return m_dim;
	}
};
#undef TMD_MATRIX_DEFINE_OPERATORS
////////////////////////////////////////////////////////////////////////////////////////////

}