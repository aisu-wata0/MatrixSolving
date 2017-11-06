#ifndef SUBST_H
#define SUBST_H

#include <assert.h>

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
void subst(TMatrix& T, XMatrix& X, IMatrix& I, vector<size_t>& P, size_t col){
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
 * @brief Solves LU*X = B, B having LU.size() columns of independant terms.
 * Resulting in LU.size() columns of X, each being the solution to LU*x = B.col(j).
 * Tiling on L0, SSE, and Unrolling on i,j
 * @param LU triangular matrix, Lower or Upper
 * @param X Matrix of solutions
 * @param B Matrix of independent terms
 * @param P Permutation obtained in GaussEl()
 */
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0AU(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
// Defines to index Matrices, if direction is backwards, access is reversed
#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-(i), (size-1)-(j)))
#define indvi(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)/vn-(i), (size-1)-(j)))
#define indvj(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)-(i), (size-1)/vn-(j)))

	size_t size = X.sizeMem();
	size_t i, j, k, kv;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax;
	size_t bstep[5];
	size_t isrt;
	const size_t iunr = 2;
	const size_t junr = 4;
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${32} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/

	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);

	size_t vn = X.vecN(); // number of doubles in vec
	vec<double> acc[iunr*junr];
	
#define vect(v) for(size_t v=0; v < vn; ++v) // ease vectorization
#define unrll(u,step) for(size_t u = 0; u < step; ++u) // ease unrolling
#define unr(iu,iunr,ju,junr) unrll(iu,iunr) unrll(ju,junr) // unroll 2 dimensions
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0]) // L1 tiling
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size); // setting tile limits
		jmax = min(bj[0]+bstep[0] , size);
		if(direction == Direction::Forwards)
			isrt = bi[0];
		else isrt = max(bi[0], X.pad());
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax -(iunr-1); i += iunr) { // i unroll
				for (j = bj[0]; j < jmax -(junr-1); j += junr) { // j unroll
assert(((direction == Direction::Backwards) && (((size-1-bk[0])-(vn-1)) % 4 == 0))
|| ((direction == Direction::Forwards) && (bk[0] % 4 == 0))); // So that vectorization doesn't give segfault -O2
// Multiply current tile: i,j = A krow * IA kcol
// For (i,j): from i to i+iunr; from j to j+junr
#define kloop(iunr, junr)	\
					unr(iu,iunr,ju,junr) vect(v) acc[iu*junr + ju][v] = 0;	\
					/*vectorized loop*/	\
					for (kv = bk[0]/vn; kv < (bk[0]+bstep[0])/vn; ++kv)	\
						unr(iu,iunr,ju,junr)	\
						acc[iu*junr+ju].v += indvj(LU, i+iu, kv).v * indvi(X, kv, j+ju).v;	\
					unr(iu,iunr,ju,junr) /*vect result sum*/	\
					vect(v) ind(X, i+iu, j+ju) -= acc[iu*junr+ju][v];
// end define
					kloop(iunr,junr)
				}
				for (j = j; j < jmax; ++j) { // j unroll reminder
					kloop(iunr,1)
				}
			}
			for (i = i; i < imax; ++i) { // i unroll remainder
				for (j = bj[0]; j < jmax -(junr-1); j += junr) { // j unroll
					kloop(1,junr)
				}
				for (j = j; j < jmax; ++j) { // j unroll reminder
					kloop(1,1)
				}
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = isrt; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				for (k = bk[0]; k < i; ++k)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
#undef vect
#undef unrll
#undef unr
#undef kloop
#undef ind
#undef indvi
#undef indvj
}

/**
 * @brief Solves LU*X = B, B having LU.size() columns of independant terms.
 * Resulting in LU.size() columns of X, each being the solution to LU*x = B.col(j).
 * Tiling on L0 and SSE
 * @param LU triangular matrix, Lower or Upper
 * @param X Matrix of solutions
 * @param B Matrix of independent terms
 * @param P Permutation obtained in GaussEl()
 */
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0A(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
// Defines to index Matrices, if direction is backwards, access is reversed
#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-(i), (size-1)-(j)))
#define indvi(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)/vn-(i), (size-1)-(j)))
#define indvj(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)-(i), (size-1)/vn-(j)))
	
	size_t size = X.sizeMem();
	size_t i, j, k, kv;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax;
	size_t bstep[5];
	size_t isrt;
	//const size_t unr = 2;
	//double acc[unr*unr];
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${32} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);
	size_t vn = X.vecN();
	vec<double> acc;
	
#define vect(v) for(size_t v = 0; v < vn; ++v)
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size);
		jmax = min(bj[0]+bstep[0] , size);
		if(direction == Direction::Forwards)
			isrt = bi[0];
		else isrt = max(bi[0], X.pad());
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				assert( ((direction == Direction::Backwards) && (((size-1-bk[0])-(vn-1)) % 4 == 0))
				|| ((direction == Direction::Forwards) && (bk[0] % 4 == 0)));
				vect(v)
					acc[v] = 0;
				for (kv = bk[0]/vn; kv < (bk[0]+bstep[0])/vn; ++kv)
					acc.v = acc.v - indvj(LU, i, kv).v * indvi(X, kv, j).v;
				vect(v)
					ind(X,i,j) += acc[v];
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = isrt; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				for (k = bk[0]; k < i; ++k)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
#undef vect
#undef ind
#undef indvi
#undef indvj
}

/**
 * @brief Solves LU*X = B, B having LU.size() columns of independant terms.
 * Resulting in LU.size() columns of X, each being the solution to LU*x = B.col(j).
 * Tiling on L0
 * @param LU triangular matrix, Lower or Upper
 * @param X Matrix of solutions
 * @param B Matrix of independent terms
 * @param P Permutation obtained in GaussEl()
 */
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
// Defines to index Matrices, if direction is backwards, access is reversed
#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-(i), (size-1)-(j)))
	
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax, kmax;
	size_t bstep[5];
	bstep[0] = 24;
	bstep[1] = bstep[0]*3;
	
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
#undef ind
}


/**
 * @brief Solves LU*X = B, B having LU.size() columns of independant terms.
 * Resulting in LU.size() columns of X, each being the solution to LU*x = B.col(j).
 * Tiling on L0 and L1
 * @param LU triangular matrix, Lower or Upper
 * @param X Matrix of solutions
 * @param B Matrix of independent terms
 * @param P Permutation obtained in GaussEl()
 */
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU10(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
// Defines to index Matrices, if direction is backwards, access is reversed
#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-(i), (size-1)-(j)))
	
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
#undef ind
}

/**
 * @brief Solves LU*X = B, B having LU.size() columns of independant terms.
 * Resulting in LU.size() columns of X, each being the solution to LU*x = B.col(j).
 * @param LU triangular matrix, Lower or Upper
 * @param X Matrix of solutions
 * @param B Matrix of independent terms
 * @param P Permutation obtained in GaussEl()
 */
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
// Defines to index Matrices, if direction is backwards, access is reversed
#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-(i), (size-1)-(j)))
	
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
#undef ind
}


}
#endif