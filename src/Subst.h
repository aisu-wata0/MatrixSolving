#ifndef SUBST_H
#define SUBST_H
/**
@file Subst.h
*/

#include "Matrix.h"

namespace std {

#define max(x,y) ((x) > (y) ? x : y)

enum SubstDirection {
	SubstForwards,
	SubstBackwards,
};

enum SubstDiagonal {
	DiagonalUnit,
	DiagonalValue,
};

enum SubstPermutation {
	SubstPermute,
	SubstNoPermute,
};

template<typename T, typename A>
inline double& at(std::vector<T,A>& vec, long i, long j) {
	return vec.at(i);
}
template<typename T, typename A>
inline const double& at(std::vector<T,A> const& vec, long i, long j) {
	return vec.at(i);
}

/**
 * @brief Subtitution method for linear sistems, forward or backward, uses P from pivoting on the LU decomposition
 * @param T Triangular coeficient Matrix
 * @param X Variables to be found
 * @param I Independant terms
 * @param forward true is forward subst, false is backward.
 * @param P from LU Decomposition
 * @param unit_diagonal true if diagonal is equal to 1
 * @param col Column of the matrix I to be used as B
 */
template<SubstDirection direction, SubstDiagonal diagonal, SubstPermutation Permutation, class TMatrix, class XMatrix, class IMatrix>
void subst(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	long i, j;
	int step;
	long size = T.size;
	
	if(direction == SubstForwards){
		i = 0;
		step = +1;
	} else {
		i = size-1;
		step = -1;
	}

	for (; i >= 0 && i < size; i += step) {
		if(Permutation == SubstPermute)
			{ at(X, i,col) = at(I, P.at(i),col); }
		else{ at(X, i,col) = at(I, i,col); }
		if(direction == SubstForwards)
			{ j = 0; }
		else{ j = size-1; }
		
		for (; j != i; j += step) {
			at(X, i,col) -= at(X, j,col) * T.at(i,j);
		}
		
		if(diagonal != DiagonalUnit)
			{ at(X, i,col) /= T.at(i,i); }
	}
}


}
#endif