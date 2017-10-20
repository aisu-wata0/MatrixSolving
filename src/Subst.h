#ifndef SUBST_H
#define SUBST_H
/**
@file Subst.h
*/

#include "Matrix.h"

namespace std {


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
template<SubstDirection direction, SubstDiagonal diagonal, class TMatrix>
void subst_P(TMatrix& T, vector<double>& X, MatrixColMajor& I, vector<long>& P, long col){
	double sum;
	long i, j;
	int step;
	long size = T.size;
	
	if(direction == SubstForwards){
		i = 1;
		step = +1;
		X.at(0) = I.at(P.at(0),col);
		if(diagonal != DiagonalUnit){
			X.at(i) /= T.at(0, 0);
		}
	} else {
		i = size-2;
		step = -1;
		X.at(size-1) = I.at(P.at(size-1),col);
		if(diagonal != DiagonalUnit){
			X.at(i) /= T.at(size-1, size-1);
		}
	}

	for (; i >= 0 && i < size; i += step) {
		sum = I.at(P.at(i), col);
		if(direction == SubstForwards) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j) * T.at(i, j);
		}
		X.at(i) = sum;
		if(diagonal != DiagonalUnit){
			X.at(i) /= T.at(i, i);
		}
	}
}
/**
 * @brief Subtitution method for linear sistems, forward or backward
 * @param T Triangular coeficient Matrix
 * @param X Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 * @param col Column of the matrix to be used as B
 */
template<SubstDirection direction, class TMatrix>
void subst(TMatrix& T, MatrixColMajor& X, vector<double>& B, long col) {
	long i, j;
	int step;
	long size = T.size;

	if(direction == SubstForwards){
		i = 1;
		step = +1;
		X.at(0,col) = B.at(0) / T.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		X.at(size-1,col) = B.at(size-1) / T.at(size-1, size-1);
	}
	
	for (; i >= 0 && i < size; i += step) {
		X.at(i,col) = B.at(i);
		if(direction == SubstForwards) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			X.at(i,col) -= X.at(j,col) * T.at(i, j);
		}
		X.at(i,col) = X.at(i,col) / T.at(i, i);
	}
	
//	int bstep = CACHE_LINE_SIZE * step;
//	// forwards
//	for(long bi = 0; bi < size; bi += bstep){
//		X.at(i,col) = B.at(i);
//		for(long bj = 0; bj != bi; bj += bstep){
//			// go though current block
//			for(long i = bi; i < bi+bstep; i += step){
//				for(long j = bj; j < bj+bstep; j += step){
//					X.at(i,col) -= X.at(j,col) * T.at(i, j);
//				}
//			}
//		}
//		X.at(i,col) = X.at(i,col) / T.at(i, i);
//	}
}


}
#endif