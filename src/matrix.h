#pragma once

#include <cassert>
#include <vector>

namespace tmd {
// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain
#define TMD_MATRIX_DEFINE_OPERATORS \
	const T& operator()(Size_Type i) const { return m_data[i]; } \
	      T& operator()(Size_Type i)       { return m_data[i]; } \
	const T& operator()(Size_Type i, Size_Type j) const { return m_data[index(i, j)]; } \
	      T& operator()(Size_Type i, Size_Type j)       { return m_data[index(i, j)]; }

////////////////////////////////////////////////////////////////////////
template<typename T>
class Matrix {
	std::vector<T> m_data;
	Size_Type m_i, m_j;
public:
	// column-major
	Matrix() : m_i(0), m_j(0) {}
	Matrix(Size_Type i, Size_Type j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {}
	Size_Type index(Size_Type i, Size_Type j) const {//return the index in m_data, column-major
		assert(j < m_j);
		assert(i < m_i);
		return i + m_i*j; //column-major
	}
	void resize(Size_Type m, Size_Type n, const T& filler_val);//resize to mxn, new sizes should be the same or greater than the old, preserves original data
	void append(const Matrix<T>& x, const T& filler_val);//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	Size_Type dim_1() const { return m_i; }
	Size_Type dim_2() const { return m_j; }
};
// //////////////////////////////////////////////////////////////////////////
template<typename T>
class Triangular_Matrix {
	std::vector<T> m_data;
	Size_Type m_dim;
public:
	Triangular_Matrix();//m_dim(0)
	Triangular_Matrix(Size_Type n, const T& filler_val);//m_data(n*(n+1)/2, filler_val), m_dim(n)
	Size_Type index(Size_Type i, Size_Type j) const;//return corresponding index in vector m_data
	//only upper-right triangular{i,j}
	// 0  1  3  6
	//    2  4  7
	//       5  8
	//          9
	Size_Type index_permissive(Size_Type i, Size_Type j) const;//return (i < j) ? index(i, j) : index(j, i);
	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	Size_Type dim() const;//return m_dim, dimension of the Matrix
};
// ///////////////////////////////////////////////////////////////////////////
// //no diagonal elements
// template<typename T>
// class strictly_triangular_matrix {
// 	std::vector<T> m_data;
// 	Size_Type m_dim;
// public:
// 	strictly_triangular_matrix();//m_dim(0)
// 	strictly_triangular_matrix(Size_Type n, const T& filler_val);//m_data(n*(n-1)/2, filler_val), m_dim(n)
// 	Size_Type index(Size_Type i, Size_Type j) const;//return corresponding index in vector m_data
// 	Size_Type index_permissive(Size_Type i, Size_Type j) const;//return (i < j) ? index(i, j) : index(j, i);
// 	void resize(Size_Type n, const T& filler_val);//resize m_data to has only the upper right Matrix elements, new sizes should be the same or greater than the old, preserves original data
// 	void append(const strictly_triangular_matrix<T>& m, const T& filler_val);//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
// 	void append(const Matrix<T>& rectangular, const strictly_triangular_matrix<T>& triangular);//arguments rectangular and triangular must have the same dim
// 	//if(rectangular.dim_2() == 0) do nothing
// 	//if(rectangular.dim_1() == 0) assign argument triangular to this object
// 	//else append these three matrices like the following
// 	// i  i  i  r  r  r  r
// 	//    i  i  r  r  r  r
// 	//       i  r  r  r  r
// 	//          t  t  t  t
// 	//             t  t  t
// 	//                t  t
// 	//                   t
// 	//where i,r,t are this, rectangular and triangular, respectively.
// 	TMD_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
// 	Size_Type dim() const;//return m_dim, dimension of the Matrix
// };
#undef TMD_MATRIX_DEFINE_OPERATORS


///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//Matrix
template<typename T>
void Matrix<T>::resize(Size_Type m, Size_Type n, const T& filler_val) {
	if(m == dim_1() && n == dim_2()) return; // no-op
	assert(m >= dim_1());//new sizes should be the same or greater than the old
	assert(n >= dim_2());
	std::vector<T> tmp(m*n, filler_val);
	for(Size_Type i = 0; i < m_i; i++) {
		for(Size_Type j = 0; j < m_j; j++) {
			tmp[i+m*j] = (*this)(i, j);
		}
	}
	m_data = tmp;
	m_i = m;
	m_j = n;
}
template<typename T>
void Matrix<T>::append(const Matrix<T>& x, const T& filler_val) {
	Size_Type m = dim_1();
	Size_Type n = dim_2();
	resize(m + x.dim_1(), n + x.dim_2(), filler_val);
	for(Size_Type i = 0; i < x.dim_1(); i++) {
		for(Size_Type j = 0; j < x.dim_2(); j++) {
			(*this)(i+m, j+n) = x(i, j);
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////
//only upper-right triangular{i,j}
// 0  1  3  6
//    2  4  7
//       5  8
//          9
inline Size_Type triangular_matrix_index(Size_Type n, Size_Type i, Size_Type j) {
	assert(j < n);
	assert(i <= j);
	return i + j*(j+1)/2;
}
//permissive
inline Size_Type triangular_matrix_index_permissive(Size_Type n, Size_Type i, Size_Type j) {
	return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
}
//Triangular_Matrix
template<typename T>
Triangular_Matrix<T>::Triangular_Matrix() : m_dim(0) {}
template<typename T>
Triangular_Matrix<T>::Triangular_Matrix(Size_Type n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {}
template<typename T>
Size_Type Triangular_Matrix<T>::index(Size_Type i, Size_Type j) const { return triangular_matrix_index(m_dim, i, j); }
template<typename T>
Size_Type Triangular_Matrix<T>::index_permissive(Size_Type i, Size_Type j) const { return (i < j) ? index(i, j) : index(j, i); }
template<typename T>
Size_Type Triangular_Matrix<T>::dim() const { return m_dim; }



// ////////////////////////////////////////////////////////////////////////////////////////////
// //strictly_triangular_matrix
// template<typename T>
// strictly_triangular_matrix<T>::strictly_triangular_matrix() : m_dim(0) {}
// template<typename T>
// strictly_triangular_matrix<T>::strictly_triangular_matrix(Size_Type n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {}
// template<typename T>
// Size_Type strictly_triangular_matrix<T>::index(Size_Type i, Size_Type j) const {
// 	assert(j < m_dim);
// 	assert(i < j);
// 	assert(j >= 1); // by implication, really
// 	return i + j*(j-1)/2;
// }
// template<typename T>
// Size_Type strictly_triangular_matrix<T>::index_permissive(Size_Type i, Size_Type j) const { return (i < j) ? index(i, j) : index(j, i); }

// template<typename T>
// void strictly_triangular_matrix<T>::resize(Size_Type n, const T& filler_val) {
// 	if(n == m_dim) return; // no-op
// 	assert(n > m_dim);
// 	m_dim = n;
// 	m_data.resize(n*(n-1)/2, filler_val); // preserves original data
// }
// template<typename T>
// void strictly_triangular_matrix<T>::append(const strictly_triangular_matrix<T>& m, const T& filler_val) {
// 	Size_Type n = dim();
// 	resize(n + m.dim(), filler_val);
// 	for(Size_Type i = 0; i < m.dim(); ++i) {
// 		for(Size_Type j = i+1; j < m.dim(); ++j) {
// 			(*this)(i+n, j+n) = m(i, j);
// 		}
// 	}
// }
// template<typename T>
// void strictly_triangular_matrix<T>::append(const Matrix<T>& rectangular, const strictly_triangular_matrix<T>& triangular) {
// 	assert(dim() == rectangular.dim_1());
// 	assert(rectangular.dim_2() == triangular.dim());
// 	// a filler value is needed by append or resize
// 	// we will use a value from rectangular as the filler value
// 	// but it can not be obtained if dim_1 or dim_2 is 0
// 	// these cases have to be considered separately
// 	if(rectangular.dim_2() == 0) return;
// 	if(rectangular.dim_1() == 0) {
// 		(*this) = triangular;
// 		return;
// 	}
// 	const T& filler_val = rectangular(0, 0); // needed by 'append below'
// 	Size_Type n = dim();
// 	append(triangular, filler_val);
// 	for(Size_Type i = 0; i < rectangular.dim_1(); ++i) {
// 		for(Size_Type j = 0; j < rectangular.dim_2(); ++j) {
// 			(*this)(i, n + j) = rectangular(i, j);
// 		}
// 	}
// }
// template<typename T>
// Size_Type strictly_triangular_matrix<T>::dim() const { return m_dim; }


}