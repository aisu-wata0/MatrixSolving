
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctgmath>
//#include <likwid.h>
#include <unistd.h>

#include "Matrix.hpp"
#include "GaussEl.hpp"
#include "Subst.hpp"
#include "Chronometer.hpp"

namespace gm {
using namespace std;

double total_time_iter = 0.0;
double total_time_residue = 0.0;
double lu_time = 0.0;

template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLU0(LUMatrix& LU, IAMatrix& X, IMatrix& B, varray<size_t>& P){
	static
	MatrixColMajor<double> Z(X.size());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0AU<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0AU<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
}

/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \n
 * Tiling on L0, SSE, unrolling on i,j
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AUIJ(AMatrix& A, IAMatrix& IA, IMatrix& R){
	ssize_t size = A.size();
	ssize_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	ssize_t bstep[5];
	/**/
	const ssize_t iunr = 2;
	const ssize_t junr = 4;
	/* export GCC_ARGS=" -D IUNRLL=${2} -D JUNRLL=${4}"* const size_t iunr = IUNRLL; const size_t junr = JUNRLL;/**/
	vec<double> acc[iunr*junr];
	/**/
	bstep[0] = B2L1;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"* bstep[0] = L0; bstep[1] = bstep[0]*L1M;/**/
	ssize_t i, j, k, kv;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// Multiply A*IA and subtract from R
	ssize_t vn = R.vecN(); // number of elements on the register (vectorization)
	
#define vect(v) for(ssize_t v=0; v < vn; ++v) // ease vectorization
#define unrll(u,step) for(size_t u = 0; u < step; ++u) // ease unrolling
#define unr(iu,iunr,ju,junr) unrll(iu,iunr) unrll(ju,junr) // unroll 2 dimensions
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0]) // L1 tiling
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		ssize_t imax = min(bi[0]+bstep[0], size); // setting tile limits
		ssize_t jmax = min(bj[0]+bstep[0], size);
		ssize_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax -(iunr-1); i += iunr) { // i unroll
			for (j = bj[0]; j < jmax -(junr-1); j += junr) { // j unroll
// Multiply current tile: i,j = A krow * IA kcol
// For (i,j): from i to i+iunr; from j to j+junr
#define kloop(iunr, junr)	\
				unr(iu,iunr,ju,junr) vect(v) acc[iu*junr + ju][v] = 0;	\
				for (kv = bk[0]/vn; kv < kmax/vn; ++kv) /*vectorized loop*/	\
					unr(iu,iunr,ju,junr)	\
					acc[iu*junr+ju].v += A.atv(i+iu, kv).v * IA.atv(kv, j+ju).v;	\
				for(k = kv*vn; k < kmax; ++k) /*vect remainder*/	\
					unr(iu,iunr,ju,junr)	\
					R.at(i+iu, j+ju) -= A.at(i+iu, k) * IA.at(k, j+ju);	\
				unr(iu,iunr,ju,junr) /*vect result sum*/	\
				vect(v) R.at(i+iu, j+ju) -= acc[iu*junr+ju][v];
// end define
				kloop(iunr, junr)
			}
			for(j = j; j < jmax; ++j){ // j unroll reminder
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
	}
#undef unrll
#undef kloop
	// Calculate norm error from R
	ssize_t iv;
	double errNorm = 0;
	vec<double> errNormV{0};
	for(j = 0; j < size; ++j){
		for(iv = 0; iv < R.sizeVec(); ++iv) // vect loop
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(i = R.remStart(); i < R.size(); ++i) // vect remainder
			errNormV[R.vecN()-1] += R.at(i,j)*R.at(i,j);
	}
	vect(v) errNorm += errNormV[v]; // vect result sum
	
	return sqrt(errNorm);
#undef vect
}
/**
 * @brief Calculates inverse of A into IA
 * @param LU decomposition of A
 * @param IA return value, no init needed
 * @param P LU pivot permutation
 * @param iter_n
 */
template<class AMatrix, class LUMatrix, class IAMatrix>
void inverse_refining(AMatrix& A, LUMatrix& LU, IAMatrix& IA, varray<size_t>& P, long iter_n){
	long i=0;
	// number of digits of iter_n, for pretty printing
	long digits = (long)log10((double) iter_n) + 1;
	double c_residue;
	//double l_residue;
	size_t size = A.size();
	// Optm: iterating line by line
	MatrixColMajor<double> W(A.size()), R(A.size());
	
	for(size_t j = 0; j < size; ++j){
		for(size_t i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(size_t i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}

	//LIKWID_MARKER_START("INV");
	
	//solveMLU(LU, IA, R, P);
	solveMLU0(LU, IA, R, P);
	
	//LIKWID_MARKER_STOP("INV");
	//LIKWID_MARKER_START("RES");
	
	c_residue = residue0AUIJ(A, IA, R);
	
	//LIKWID_MARKER_STOP("RES");
	
	cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	while(i < iter_n){
		// (abs(l_residue - c_residue)/c_residue > EPSILON) && (l_residue > c_residue)
		// relative approximate error
		i += 1;
		// R: residue of IA

		timer.start();
		//LIKWID_MARKER_START("INV");
		
		//solveMLU(LU, W, R, P);
		solveMLU0(LU, W, R, P);
		
		//LIKWID_MARKER_STOP("INV");
		// W: residues of each variable of IA
		// adjust IA with found errors
		//LIKWID_MARKER_START("SUM");
		
		//add(IA, W);
#define unrll(u,step) for(ssize_t u = 0; u < step; ++u) // ease unrolling
		ssize_t i, j, junr = 8, size = IA.size();
		for(j=0; j < size -(junr-1); j += junr)
			for(i=0; i < size; ++i)
				unrll(ju,junr)
				IA.at(i,j+ju) += W.at(i,j+ju);
		for(j = j; j < size; ++j)
			for(i=0; i < size; ++i)
				IA.at(i,j) += W.at(i,j);
#undef unrll
		
		//LIKWID_MARKER_STOP("SUM");
		total_time_iter += timer.tickAverage();
		
		//l_residue = c_residue;
		timer.start();
		//LIKWID_MARKER_START("RES");
		
		c_residue = residue0AUIJ(A, IA, R);
		
		//LIKWID_MARKER_STOP("RES");
		total_time_residue += timer.tick();
		
		cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	}
}


}