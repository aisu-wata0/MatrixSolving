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


}
#endif