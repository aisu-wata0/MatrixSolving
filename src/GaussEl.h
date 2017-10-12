#ifndef GAUSSEL_H
#define GAUSSEL_H
/**
@file GaussEl.h
*/
namespace std {

#include "Double.h"
#include "Matrix.h"

/**
 @brief For the matrix LU finds its LU decomposition overwriting it
 Has partial pivoting, stores final indexes in P
 @param LU Matrix to be decomposed Output: lower triangle of this matrix will store L 1 diagonal implicit, upper triangle stores U
 @param P Permutation vector resulting of the pivoting
 */
void GaussEl(Matrix A, Matrix& LU, vector<long>& P) {
	// copy A to LU
	LU.resize(A.size);
	LU.set(A);
	// initializing permutation vector
	for(long i = 0; i < LU.size; i++){
		P.at(i) = i;
	}

	for(long p = 0; p < LU.size; p++){
		// for each pivot
		/* partial pivoting */
		long maxRow = p;
		for(long i = p+1; i < LU.size; i++){
			// for each value below the p pivot
			if(abs(LU.at(i,p)) > abs(LU.at(maxRow,p))) maxRow = i;
		} // finds max value
		// pivots rows of U
		LU.swap_rows(maxRow, p);
		swap(P.at(p), P.at(maxRow));

		// LU.at(p,p) = 1; implicit
		if(close_zero(LU.at(p,p))){
			fprintf(stderr, "Found a pivot == 0, system is not solvable with partial pivoting");
			exit(EXIT_FAILURE);
		}
		//for (long i = p+1; i < LU.size; i++) {	//going from below pivot to end
		for (long i = LU.size-1; i >= p+1; i--) {	//going from end to pivot
			// for each line below pivot
			if (!close_zero(LU.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
				// subtract pivot from current line (in U)
				for (long k = p+1; k < LU.size; k++) {
					// for each collumn starting from pivot's
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