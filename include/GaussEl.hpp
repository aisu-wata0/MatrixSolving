#ifndef GAUSSEL_H
#define GAUSSEL_H

#include "Double.h"
#include "Matrix.hpp"

namespace gm {
using namespace std;

/**
 @brief For the matrix LU finds its LU decomposition overwriting it
 Has partial pivoting, stores final indexes in P
 @param LU Matrix to be decomposed Output: lower triangle of this matrix will store L 1 diagonal implicit, upper triangle stores U
 @param P Permutation vector resulting of the pivoting
 */
inline void GaussEl(const Matrix<double>& A, Matrix<double>& LU, varray<size_t>& P) {
	// copy A to LU
	set(LU, A);
	// initializing permutation vector
	for(size_t i = 0; i < A.sizeMem(); i++){
		P.at(i) = i;
	}

	// for each pivot
	for(size_t p = 0; p < A.size(); p++){
		/* partial pivoting */
		size_t maxRow = p;
		for(size_t i = p+1; i < A.size(); i++){
			// for each value below the p pivot
			if(abs(LU.at(i,p)) > abs(LU.at(maxRow,p))) maxRow = i;
		} // finds max value
		// pivots rows of U
		swap_rows(LU, maxRow, p);
		swap(P.at(p), P.at(maxRow));

		if(close_zero(LU.at(p,p))){
			fprintf(stderr, "Found a pivot == 0, system is not solvable with partial pivoting");
			exit(EXIT_FAILURE);
		}
		// LU.at(p,p) = 1; implicit
		//for(size_t i = p+1; i < A.size; i++){	// going from pivot+1 to end
		for(size_t i = A.size()-1; i >= p+1; i--){	// going from end to pivot+1 Optm: If no pivoting ocurred: 1 less cache miss
			if(not close_zero(LU.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
				// subtract pivot row U.at(p, _) from current row LU.at(i, _)
				for(size_t k = p+1; k < A.size(); k++){
					LU.at(i, k) -= LU.at(p, k) * LU.at(i, p);
					// mulitply pivot line value to multiplier
				}
			} else {
				// pivot not subtracted from line
				LU.at(i, p) = 0.0;
			}
		}
	}
}


}
#endif