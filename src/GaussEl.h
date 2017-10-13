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
template<class MLower, class MUpper>
void GaussEl(Matrix& A, MLower& L, MUpper& U, vector<long>& P) {
	// copy A to LU
	/* Optm: test *
	set(L, U, A);
	/**/
	L.set(A);
	U.set(A);
	/**/
	
	// initializing permutation vector
	for(long i = 0; i < A.size; i++){
		P.at(i) = i;
	}

	// for each pivot
	for(long p = 0; p < A.size; p++){
		/* partial pivoting */
		long maxRow = p;
		for(long i = p+1; i < A.size; i++){
			// for each value below the p pivot
			if(abs(L.at(i,p)) > abs(L.at(maxRow,p))) maxRow = i;
		} // finds max value
		// pivots rows of U
		swap_rows(L, U, maxRow, p);
		swap(P.at(p), P.at(maxRow));

		if(close_zero(U.at(p,p))){
			fprintf(stderr, "Found a pivot == 0, system is not solvable with partial pivoting");
			exit(EXIT_FAILURE);
		}
		L.at(p,p) = 1;
		
		// for each line below pivot
		//for (long i = A.size-1; i >= p+1; i--) {	//going from end to pivot
		for (long i = p+1; i < A.size; i++) {	//going from below pivot to end
			if (!close_zero(L.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				L.at(i, p) = L.at(i, p)/U.at(p, p);
				
				// subtract pivot row U.at(p, _) from current row LU.at(i, _)
				// L side
				for (long k = p+1; k < i; k++) {
					// for each collumn starting from pivot's
					L.at(i, k) -= U.at(p, k) * L.at(i, p);
					// mulitply pivot line value to multiplier
				}
				// U side
				for (long k = i; k < A.size; k++) {
					// for each collumn starting from pivot's
					U.at(i, k) -= U.at(p, k) * L.at(i, p);
					// mulitply pivot line value to multiplier
				}
			} else {
				// pivot not subtracted from line
				L.at(i, p) = 0.0;
			}
		}
	}
}


}
#endif