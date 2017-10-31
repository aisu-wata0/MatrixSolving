#ifndef SUBST_H
#define SUBST_H

#include "Matrix.hpp"

namespace gm {

#define max(x,y) ((x) > (y) ? x : y)

enum class Direction {
	Forwards,
	Backwards,
};

enum class Diagonal {
	Unit,
	Value,
};

enum class Permute {
	True,
	False,
};

template<class Elem>
Elem& at(varray<Elem>& arr, size_t i, size_t j){
	return arr.at(i);
}
template<class Elem>
const Elem& at(varray<Elem> const& arr, size_t i, size_t j){
	return arr.at(i);
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
template<Direction direction, Diagonal diagonal, Permute permute, class TMatrix, class XMatrix, class IMatrix>
void subst(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	size_t i, j;
	int step;
	size_t size = T.size();
	
	if(direction == Direction::Forwards){
		i = 0;
		step = +1;
	} else {
		i = size-1;
		step = -1;
	}

	for (; i >= 0 && i < size; i += step) {
		if(permute == Permute::True)
			{ at(X, i,col) = at(I, P.at(i),col); }
		else{ at(X, i,col) = at(I, i,col); }
		if(direction == Direction::Forwards)
			{ j = 0; }
		else{ j = size-1; }
		
		for (; j != i; j += step) {
			at(X, i,col) -= at(X, j,col) * T.at(i,j);
		}
		
		if(diagonal != Diagonal::Unit)
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
template<Direction direction, Diagonal diagonal, Permute permute, class TMatrix, class XMatrix, class IMatrix>
void substR(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	long i, j;
	long bi, bj;
	int step;
	int bstep;
	long size = T.size();
	long endi, endj;
	
	if(direction == Direction::Forwards){
		bi = 0;
		step = +1;
		bstep = step * BL1;
	} else {
		bi = size-1;
		step = -1;
		bstep = step * BL1;
	}
	
	for (; bi >= 0 && bi < size ; bi += bstep) {
		if(direction == Direction::Forwards)
			{ bj = 0;		endi = min(bi + bstep, size); }
		else{ bj = size-1;	endi = max(bi + bstep, -1); }
		
		for(int i = bi; i != endi; i += step)
			if(permute == Permute::True)
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
			if(direction == Direction::Forwards)
				{ endj = min(bj + bstep, i); }
			else{ endj = max(bj + bstep, i); }
			
			for (j = bj; j != endj; j += step) {
				at(X, i,col) = at(X, i,col) - T.at(i,j) * at(X, j,col);
			}
			if(diagonal != Diagonal::Unit)
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
template<Direction direction, Diagonal diagonal, Permute permute, class TMatrix, class XMatrix, class IMatrix>
void substUnroll(TMatrix& T, XMatrix& X, IMatrix& I, vector<long>& P, long col){
	long i, j;
	int step;
	int bstep;
	long size = T.size;
	
	if(direction == Direction::Forwards){
		i = 0;
		step = +1;
		bstep = step * BL1;
	} else {
		i = size-1;
		step = -1;
	}
	
	#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
	for (; i >= 0 && i < size; i += step) {
		at(X, i,col) = at(I, i,col);
		for (j = 0; j < ROUND_DOWN(i,BL1); j += bstep) {
			for(int c = 0; c < bstep; ++c)
				at(X, i,col) -= at(X, j+c,col) * T.at(i,j+c);
		}
		// Remainder
		for(; j < i; ++j)
			{ at(X, i,col) -= at(X, j,col) * T.at(i,j); }
		
		if(diagonal != Diagonal::Unit)
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
template<Direction direction, class TMatrix>
void subst(TMatrix& T, MatrixColMajor<double>& X, vector<double>& B, long col) {
	long size = T.size;
	long bi; long bj;
	long i; long j;
	long iend; long jend;
	long step;

	// TODO test iterate blocks by col (bj) instead of as currently by row (bi)
	if(direction == Direction::Forwards){
		bj = 0;
		step = +1;
	} else {
		bj = size-1;
		step = -1;
	}
	long bstep = step * BL1;
	
	for(; bj >= 0 && bj < size; bj += bstep){
		if(direction == Direction::Forwards) {
			jend = bj + bstep > size ? size : bj + bstep;
		} else {
			jend = bj + bstep < 0 ? -1 : bj + bstep;
		}
		// go though diagonal block
		for(i = bj; i != jend; i += step){
			if((bj == 0 && direction == Direction::Forwards) || (bj == size-1 && direction == Direction::Backwards)){
				X.at(i,col) = B.at(i);
			}
			for(j = bj; j != jend && j != i; j += step){
				X.at(i,col) -= X.at(j,col) * T.at(i, j);
			}
			X.at(i,col) = X.at(i,col) / T.at(i, i);
		}
		for(bi = bj+bstep; bi >= 0 && bi < size; bi += bstep){
			if(direction == Direction::Forwards) {
				iend = bi + bstep > size ? size : bi + bstep;
			} else {
				iend = bi + bstep < 0 ? -1 : bi + bstep;
			}
			if((bj == 0 && direction == Direction::Forwards) || (bj == size-1 && direction == Direction::Backwards))
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

#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax, kmax;
	size_t bstep[5];
	//const size_t unr = 2;
	//double acc[unr*unr];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;/**/
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size);
		jmax = min(bj[0]+bstep[0] , size);
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < (bk[0]+bstep[0]); k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < i; k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
}
#undef ind

#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;	
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);
	
	for(i = 0; i < size; i += 1){
		for(j = 0; j < size; j += 1){
			for(k = 0; k != i; k += 1){
				ind(X,i,j) = ind(X,i,j) - ind(LU,i,k) * ind(X,k,j);
			}
			if(diagonal == Diagonal::Value)
				ind(X,i,j) /= ind(LU,i,i);
		}
	}
}
#undef ind

#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU10(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax, kmax;
	size_t bstep[5];
	//const size_t unr = 2;
	//double acc[unr*unr];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;/**/
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);
	
	asm("LUTiled1.0");
	for (bi[1] = 0; bi[1] < size; bi[1] += bstep[1])
	for (bj[1] = 0; bj[1] < size; bj[1] += bstep[1])
	for (bk[1] = 0; bk[1] < (bi[1]+bstep[1]); bk[1] += bstep[1]) {
		bimax[0] = min(bi[1]+bstep[1], size);
		bjmax[0] = min(bj[1]+bstep[1], size);
		for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
		for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
			bkmax[0] = min(bk[1]+bstep[1], bi[0]);
			imax = min(bi[0]+bstep[0] , size);
			jmax = min(bj[0]+bstep[0] , size);
			for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
				for (i = bi[0]; i < imax; i += 1)
				for (j = bj[0]; j < jmax; j += 1) {
					for (k = bk[0]; k < (bk[0]+bstep[0]); k += 1)
						ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				}
			} // Last block in K, diagonal
			if(bk[0] == bi[0] && bk[1] == bi[1])
			for (bk[0] = bi[0]; bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
				for (i = bi[0]; i < imax; i += 1)
				for (j = bj[0]; j < jmax; j += 1) {
					for (k = bk[0]; k < i; k += 1)
						ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
					if(diagonal == Diagonal::Value)
						ind(X, i, j) /= ind(LU, i, i);
				}
			}
		}
	}
	asm("END LUTiled1.0");
}
#undef ind


}
#endif