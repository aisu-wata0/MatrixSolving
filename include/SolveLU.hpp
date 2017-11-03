
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

/**
 * @brief Solves LU system using subst functions.
 * @param LU Matrix to find solution
 * @param X Variables to be found
 * Output: Value of found variables is stored in X
 * @param B Independent term matrix
 * @param P Permutation vector resulting of the pivoting
 * @param col Column of the matrix to be used as B
 */
template<class LUMatrix, class XMatrix, class BMatrix>
inline void solveLU(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P, long col){
	static
	varray<double> Z(LU.sizeMem());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substR<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P, col);
	// find X; Ux=Z
	substR<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P, col);
}
/**
 * @brief Finds inverse matrix,
 * also solves linear sistems A*IA = I where each of I's columns are different Bs
 * @param LU Matrix to find inverse
 * @param IA Inverse matrix to be found
 * Output: Inverse matrix is stored in IA
 * @param I Identity matrix
 * @param P Permutation vector resulting of the pivoting
 */
template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLU(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
	for(long j = 0; j < X.size(); j++){
		// for each X col solve SL to find the X col values
		static
		varray<double> Z(LU.sizeMem());
		if(Z.size() != X.size()){ Z.alloc(X.size()); }
		// find Z; LZ=B
		subst<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P, j);
		// find X; Ux=Z
		subst<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P, j);
	}
}

template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLU0(LUMatrix& LU, IAMatrix& X, IMatrix& B, varray<size_t>& P){
	static
	MatrixColMajor<double> Z(X.size());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0A<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0A<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
}

/**
 * @brief  Calculates residue into I, A*IA shuold be close to Identity
 * @param A original coef matrix
 * @param IA solution to A*IA = I
 * @param I Output residue, no init needed
 * @return Norm of the residue
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue(AMatrix& A, IAMatrix& IA, IMatrix& I){
	double err_norm = 0.0;
	size_t size = A.size();

	for(size_t j = 0; j < size; ++j){
		// for each column of the inverse
		for(size_t i = 0; i < size; i++){
			// for each line of A
			I.at(i,j) = 0;
			// multiply A line to the current inverse col
			for(size_t k = 0; k < size; ++k){
				I.at(i,j) = I.at(i,j) - A.at(i,k) * IA.at(k,j);
			}
			if(i == j){
				I.at(i,j) += 1;
			}
			
			err_norm += I.at(i,j)*I.at(i,j);
		}
	}
	return sqrt(err_norm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \
 * Uses tiling on L0.
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = 24;
	bstep[1] = bstep[0]*4;
	bstep[2] = bstep[1]*4;
	bstep[3] = bstep[2]*5;
	
	size_t i, j, k, kv;
	
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// multiply
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j)
		for (k = bk[0]; k < kmax; ++k){
			R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
		}
	}
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv)
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.vecEnd(); i < R.size(); ++i)
			errNormV[R.regEN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.regEN(); ++v) errNorm += errNormV[v];
	
	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \
 * Uses tiling on L0, SSE.
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0A(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${32} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	size_t i, j, k, kv;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	size_t nv = R.regEN();
#define vect(v) for(size_t v=0; v < nv; ++v)
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j) {
			vec<double> acc;
			vect(v) acc[v] = 0;
			for (kv = bk[0]/nv; kv < kmax/nv; ++kv)
				acc.v = acc.v - A.atv(i, kv).v * IA.atv(kv, j).v;
			for(k = kv*nv; k < kmax; ++k)
				R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
			vect(v) R.at(i, j) += acc[v];
		}
	}
#undef vect
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv)
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.vecEnd(); i < R.size(); ++i)
			errNormV[R.regEN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.regEN(); ++v) errNorm += errNormV[v];
	
	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \
 * Uses tiling on L0, SSE, Unrolling on k.
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AU(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	const size_t unr = 8;
	vec<double> acc;
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/

	size_t i, j, k, kv, rem;

	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	size_t nv = R.regEN();
#define vect(v) for(size_t v=0; v < nv; ++v)
#define unr(v) for(size_t v=0; v < unr; ++v)
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j) {
			vect(u) vect(v) acc[v] = 0;
			for (kv = bk[0]/nv; kv < kmax/nv -(unr-1); kv += unr)
				unr(u) acc.v += A.atv(i, kv+u).v * IA.atv(kv+u, j).v;
//			for (rem = 0; rem < kmax/nv % unr; ++rem) // unroll remainder
//				acc.v += A.atv(i, kv+rem).v * IA.atv(kv+rem, j).v;
			for(k = kv*nv; k < kmax; ++k) // vect remainder
				R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
			vect(v) R.at(i, j) -= acc[v]; // vect result sum
		}
	}
#undef vect
#undef unr
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv)
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.vecEnd(); i < R.size(); ++i)
			errNormV[R.regEN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.regEN(); ++v) errNorm += errNormV[v];

	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \
 * Uses tiling on L0, SSE, Unrolling on i,j. (Doen't care for the remainder of the unrolling, unnacurate)
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AUU(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	const size_t unr = 2;
	vec<double> acc[unr*unr];
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	size_t i, j, k, kv, rem;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	size_t nv = R.regEN();
#define vect(v) for(size_t v=0; v < nv; ++v)
#define unr(u,n) for(size_t u = 0; u < n; ++u)
#define unr2(iu,ju,n) unr(iu,n) unr(ju,n)
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; i += unr)
		for (j = bj[0]; j < jmax; j += unr) {
			unr2(iu,ju,unr) vect(v) acc[iu*unr + ju][v] = 0;
			for (kv = bk[0]/nv; kv < kmax/nv; ++kv)
				unr2(iu,ju,unr)
					acc[iu*unr+ju].v += A.atv(i+iu, kv).v * IA.atv(kv, j+ju).v;
//			for (rem = 0; rem < kmax/nv % unr; ++rem) // unroll remainder
//				acc[0].v -= A.atv(i, kv+rem).v * IA.atv(kv+rem, j).v;
			for(k = kv*nv; k < kmax; ++k) // vect remainder
				unr2(iu,ju,unr)
					R.at(i+iu, j+ju) -= A.at(i+iu, k) * IA.at(k, j+ju);
			unr2(iu,ju,unr)
				vect(v) R.at(i+iu, j+ju) -= acc[iu*unr+ju][v]; // vect result sum
		}
	}
#undef vect
#undef unr
#undef unr2
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv)
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.vecEnd(); i < R.size(); ++i)
			errNormV[R.regEN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.regEN(); ++v) errNorm += errNormV[v];

	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \
 * Uses tiling on L0, SSE, unrolling on i,j.
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AUIJ(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	const size_t iunr = 2;
	const size_t junr = 4;
	vec<double> acc[iunr*junr];
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	size_t i, j, k, kv, rem;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	size_t nv = R.regEN();
#define vect(v) for(size_t v=0; v < nv; ++v)
#define unri(iu) for(size_t iu = 0; iu < iunr; ++iu)
#define unrj(ju) for(size_t ju = 0; ju < junr; ++ju)
#define unr(iu,ju) unri(iu) unrj(ju)
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax -(iunr-1); i += iunr) { // i unroll
			for (j = bj[0]; j < jmax -(junr-1); j += junr) { // j unroll
				unr(iu,ju) vect(v) acc[iu*junr + ju][v] = 0;
				for (kv = bk[0]/nv; kv < kmax/nv; ++kv) // vectorized loop
					unr(iu,ju)
					acc[iu*junr+ju].v += A.atv(i+iu, kv).v * IA.atv(kv, j+ju).v;
				for(k = kv*nv; k < kmax; ++k) // vect remainder
					unr(iu,ju)
					R.at(i+iu, j+ju) -= A.at(i+iu, k) * IA.at(k, j+ju);
				unr(iu,ju) // vect result sum
				vect(v) R.at(i+iu, j+ju) -= acc[iu*junr+ju][v];
			}
			for(j = j; j < jmax; ++j){ // j unroll reminder
				unri(iu) vect(v) acc[iu*junr][v] = 0;
				for (kv = bk[0]/nv; kv < kmax/nv; ++kv) // vectorized loop
					unri(iu)
					acc[iu*junr].v += A.atv(i+iu, kv).v * IA.atv(kv, j).v;
				for(k = kv*nv; k < kmax; ++k) // vect remainder
					unri(iu)
					R.at(i+iu, j) -= A.at(i+iu, k) * IA.at(k, j);
				unri(iu) // vect result sum
				vect(v) R.at(i+iu, j) -= acc[iu*junr][v]; 
			}
		}
		for(i = i; i < imax; ++i){ // i unroll remainder
			for (j = bj[0]; j < jmax -(junr-1); j += junr) { // j unroll
				unrj(ju) vect(v) acc[ju][v] = 0;
				for (kv = bk[0]/nv; kv < kmax/nv; ++kv) // vectorized loop
					unrj(ju)
					acc[ju].v += A.atv(i, kv).v * IA.atv(kv, j+ju).v;
				for(k = kv*nv; k < kmax; ++k) // vect remainder
					unrj(ju)
					R.at(i, j+ju) -= A.at(i, k) * IA.at(k, j+ju);
				unrj(ju) // vect result sum
				vect(v) R.at(i, j+ju) -= acc[ju][v];
			}
			for(j = j; j < jmax; ++j){  // j unroll reminder
				unri(iu) vect(v) acc[iu*junr][v] = 0;
				for (kv = bk[0]/nv; kv < kmax/nv; ++kv) // vectorized loop
					unri(iu)
						acc[iu*junr].v += A.atv(i+iu, kv).v * IA.atv(kv, j).v;
				for(k = kv*nv; k < kmax; ++k) // vect remainder
					unri(iu)
						R.at(i+iu, j) -= A.at(i+iu, k) * IA.at(k, j);
				unri(iu)
					vect(v) R.at(i+iu, j) -= acc[iu*junr][v]; // vect result sum
			}
		}
	}
#undef vect
#undef unr
#undef unr2
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv)
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.vecEnd(); i < R.size(); ++i)
			errNormV[R.regEN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.regEN(); ++v) errNorm += errNormV[v];

	return sqrt(errNorm);
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
		
		for(size_t j=0; j < IA.size(); ++j)
			for(size_t i=0; i < IA.size(); ++i)
				IA.at(i,j) += W.at(i,j);
		//add(IA, W);
		
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