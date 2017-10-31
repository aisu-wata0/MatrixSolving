
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
inline void solveMLU0(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
	static
	MatrixColMajor<double> Z(X.size());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
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

template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0(AMatrix& A, IAMatrix& IA, IMatrix& I){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = 8;
	bstep[1] = bstep[0]*4;
	bstep[2] = bstep[1]*4;
	bstep[3] = bstep[2]*5;
	
	for(size_t j = 0; j < size; ++j){
		for(size_t i = 0; i < j; ++i)
			I.at(i,j) = 0;
		I.at(j,j) = 1;
		for(size_t i = j+1; i < size; ++i)
			I.at(i,j) = 0;
	}
	
	for (bi[3] = 0; bi[3] < size; bi[3] += bstep[3])
	for (bj[3] = 0; bj[3] < size; bj[3] += bstep[3])
	for (bk[3] = 0; bk[3] < size; bk[3] += bstep[3]){
		bimax[3] = min(bi[3]+bstep[3], size);
		bjmax[3] = min(bj[3]+bstep[3], size);
		bkmax[3] = min(bk[3]+bstep[3], size);
		for (bi[2] = bi[3]; bi[2] < bimax[3]; bi[2] += bstep[2])
		for (bj[2] = bj[3]; bj[2] < bjmax[3]; bj[2] += bstep[2])
		for (bk[2] = bk[3]; bk[2] < bkmax[3]; bk[2] += bstep[2]){
			bimax[2] = min(bi[2]+bstep[2], size);
			bjmax[2] = min(bj[2]+bstep[2], size);
			bkmax[2] = min(bk[2]+bstep[2], size);
			for (bi[1] = bi[2]; bi[1] < bimax[2]; bi[1] += bstep[1])
			for (bj[1] = bj[2]; bj[1] < bjmax[2]; bj[1] += bstep[1])
			for (bk[1] = bk[2]; bk[1] < bkmax[2]; bk[1] += bstep[1]){
				bimax[1] = min(bi[1]+bstep[1], size);
				bjmax[1] = min(bj[1]+bstep[1], size);
				bkmax[1] = min(bk[1]+bstep[1], size);
				for (bi[0] = bi[1]; bi[0] < bimax[1]; bi[0] += bstep[0])
				for (bj[0] = bj[1]; bj[0] < bjmax[1]; bj[0] += bstep[0])
				for (bk[0] = bk[1]; bk[0] < bkmax[1]; bk[0] += bstep[0]){
					size_t imax = min(bi[0]+bstep[0], size);
					size_t jmax = min(bj[0]+bstep[0], size);
					size_t kmax = min(bk[0]+bstep[0], size);
					for (size_t i = bi[0]; i < imax; ++i)
					for (size_t j = bj[0]; j < jmax; ++j) {
						for (size_t k = bk[0]; k < kmax; ++k) {
							I.at(i, j) = I.at(i, j) - A.at(i, k) * IA.at(k, j);
						}
					}
				}
			}
		}
	}
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < I.sizeVec(); ++iv)
			errNormV.v += I.atv(iv,j).v*I.atv(iv,j).v;
		for(size_t i = I.vecEnd(); i < I.size(); ++i)
			errNormV[I.regEN()-1] += I.at(i,j)*I.at(i,j);
	}
	for(size_t v=0; v < I.regEN(); ++v) errNorm += errNormV[v];
	
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
void inverse_refining(AMatrix& A, LUMatrix& LU, IAMatrix& IA, vector<long>& P, long iter_n){
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
	
	c_residue = residue0(A, IA, R);
	
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
		
		c_residue = residue0(A, IA, R);
		
		//LIKWID_MARKER_STOP("RES");
		total_time_residue += timer.tick();
		
		cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	}
}


}