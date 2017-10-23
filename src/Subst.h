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
void substR(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	long i, j;
	long bi, bj;
	int step;
	int bstep;
	long size = T.size;
	long endi, endj;
	
	if(direction == SubstForwards){
		bi = 0;
		step = +1;
		bstep = step * CACHE_LSZ;
	} else {
		bi = size-1;
		step = -1;
		bstep = step * CACHE_LSZ;
	}
	
	for (; bi >= 0 && bi < size ; bi += bstep) {
		if(direction == SubstForwards)
			{ bj = 0;		endi = min(bi + bstep, size); }
		else{ bj = size-1;	endi = max(bi + bstep, -1); }
		
		for(int i = bi; i != endi; i += step)
			if(Permutation == SubstPermute)
				{ at(X, i,col) = at(I, P.at(i),col); }
			else{ at(X, i,col) = at(I, i,col); }
		
		for (; bj != bi; bj += bstep) {
			for (i = bi; i != endi; i += step) {
//				if(direction == SubstForwards)
//					{ endj = min(bj + bstep, i); }
//				else{ endj = max(bj + bstep, i); }
				endj = bj + bstep;
				for (j = bj; j != endj; j += step) {
					at(X, i,col) = at(X, i,col) - T.at(i,j) * at(X, j,col);
				}
			}
		}
		// Remainder diagonal
		for (i = bi; i != endi; i += step) {
			if(direction == SubstForwards)
				{ endj = min(bj + bstep, i); }
			else{ endj = max(bj + bstep, i); }
			
			for (j = bj; j != endj; j += step) {
				at(X, i,col) = at(X, i,col) - T.at(i,j) * at(X, j,col);
			}
			if(diagonal != DiagonalUnit)
				{ at(X, i,col) /= T.at(i,i); }
		}
	}
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
void substUnroll(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	long i, j;
	int step;
	int bstep;
	long size = T.size;
	
	if(direction == SubstForwards){
		i = 0;
		step = +1;
		bstep = step * CACHE_LSZ;
	} else {
		i = size-1;
		step = -1;
	}
	
	#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
	for (; i >= 0 && i < size; i += step) {
		at(X, i,col) = at(I, i,col);
		for (j = 0; j < ROUND_DOWN(i,CACHE_LSZ); j += bstep) {
			for(int c = 0; c < bstep; ++c)
				at(X, i,col) -= at(X, j+c,col) * T.at(i,j+c);
		}
		// Remainder
		for(; j < i; ++j)
			{ at(X, i,col) -= at(X, j,col) * T.at(i,j); }
		
		if(diagonal != DiagonalUnit)
			{ at(X, i,col) /= T.at(i,i); }
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
	long size = T.size;
	long bi; long bj;
	long i; long j;
	long iend; long jend;
	long step;

	// TODO test iterate blocks by col (bj) instead of as currently by row (bi)
	if(direction == SubstForwards){
		bj = 0;
		step = +1;
	} else {
		bj = size-1;
		step = -1;
	}
	long bstep = step * CACHE_LSZ;
	
	for(; bj >= 0 && bj < size; bj += bstep){
		if(direction == SubstForwards) {
			jend = bj + bstep > size ? size : bj + bstep;
		} else {
			jend = bj + bstep < 0 ? -1 : bj + bstep;
		}
		// go though diagonal block
		for(i = bj; i != jend; i += step){
			if((bj == 0 && direction == SubstForwards) || (bj == size-1 && direction == SubstBackwards)){
				X.at(i,col) = B.at(i);
			}
			for(j = bj; j != jend && j != i; j += step){
				X.at(i,col) -= X.at(j,col) * T.at(i, j);
			}
			X.at(i,col) = X.at(i,col) / T.at(i, i);
		}
		for(bi = bj+bstep; bi >= 0 && bi < size; bi += bstep){
			if(direction == SubstForwards) {
				iend = bi + bstep > size ? size : bi + bstep;
			} else {
				iend = bi + bstep < 0 ? -1 : bi + bstep;
			}
			if((bj == 0 && direction == SubstForwards) || (bj == size-1 && direction == SubstBackwards))
				for(i = bi; i != iend; i += step){
					X.at(i,col) = B.at(i);
				}
			// go though current block
			for(i = bi; i != iend; i += step){
				for(j = bj; j != jend && j != i; j += step){
					X.at(i,col) -= X.at(j,col) * T.at(i, j);
				}
			}
		}
		// remainder
//		for(i = size - mod(size, bstep); i < size; i += step){
//			for(j = bj; j < bj+bstep; j += step){
//				cout << "(" << i << "," << j << ")" << endl; 
//			}
//		}
	}
	// remainder
//	for(i = size - mod(size, bstep); i < size; i += step){
//		for(j = size - mod(size, bstep); j != i; j += step){
//			cout << "(" << i << "," << j << ")" << endl; 
//		}
//		cout << "divided by (" << i << "," << i << ") dia" << endl;
//	}
}


}
#endif