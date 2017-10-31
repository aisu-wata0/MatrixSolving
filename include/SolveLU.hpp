
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
inline void solve_lu(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P, long col){
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
inline void substMLUZ(LUMatrix& LU, IAMatrix& Z, IMatrix& B, vector<long>& P){
	for(size_t j = 0; j < B.size(); j++){
		subst<Direction::Backwards, Diagonal::Unit, Permute::True>(LU, Z, B, P, j);
	}
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

	for(long icol=0; icol < size; icol++){
		// for each column of the inverse
		for(long i=0; i < size; i++){
			// for each line of A
			I.at(i,icol) = 0;
			// multiply A line to the current inverse col
			for(long j=0; j < size; j++){
				I.at(i,icol) -= A.at(i,j)*IA.at(j,icol);
			}
			if(i == icol){
				I.at(i,icol) += 1;
			}

			err_norm += I.at(i,icol)*I.at(i,icol);
		}
	}
	return sqrt(err_norm);
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
	
	// Optm: iterating line by line
	MatrixColMajor<double> W(A.size()), R(A.size());
	identity(R);

	//LIKWID_MARKER_START("INV");
	
	inverse(LU, IA, R, P);
	
	//LIKWID_MARKER_STOP("INV");
	//LIKWID_MARKER_START("RES");
	
	c_residue = residue(A, IA, R);
	
	//LIKWID_MARKER_STOP("RES");
	
	cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	while(i < iter_n){
		// (abs(l_residue - c_residue)/c_residue > EPSILON) && (l_residue > c_residue)
		// relative approximate error
		i += 1;
		// R: residue of IA

		timer.start();
		//LIKWID_MARKER_START("INV");
		
		inverse(LU, W, R, P);
		
		//LIKWID_MARKER_STOP("INV");
		// W: residues of each variable of IA
		// adjust IA with found errors
		//LIKWID_MARKER_START("SUM");
		
		add(IA, W);	//TOptm: add at the same time its calculating W
		
		//LIKWID_MARKER_STOP("SUM");
		total_time_iter += timer.tickAverage();
		
		//l_residue = c_residue;
		timer.start();
		//LIKWID_MARKER_START("RES");
		
		c_residue = residue(A, IA, R);
		
		//LIKWID_MARKER_STOP("RES");
		total_time_residue += timer.tick();
		
		cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	}
}


}