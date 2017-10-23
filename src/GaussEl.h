#ifndef GAUSSEL_H
#define GAUSSEL_H
/**
@file GaussEl.h
*/

#include "Double.h"
#include "Matrix.h"

namespace std {
double total = 0.0;
double average = 0.0;
Timer gauss_timer;
/**
 @brief For the matrix LU finds its LU decomposition overwriting it
 Has partial pivoting, stores final indexes in P
 @param LU Matrix to be decomposed Output: lower triangle of this matrix will store L 1 diagonal implicit, upper triangle stores U
 @param P Permutation vector resulting of the pivoting
 */
void GaussEl(const Matrix& A, Matrix& LU, vector<long>& P) {
	long p, i, j;
	long bi, bj;
	long endi, endj;
	int bstep = CACHE_LSZ;
	long size = A.size;
	// copy A to LU
	set(LU, A);
	// initializing permutation vector
	for(long i = 0; i < A.size; i++){
		P.at(i) = i;
	}
	total = 0.0;
	// for each pivot
	for(p = 0; p < A.size; p++){
		/* partial pivoting */
		long maxRow = p;
		for(i = p+1; i < A.size; i++){
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
		
		// Calc pivot multipliers
		//for(long i = p+1; i < LU.size; i++){	// going from pivot+1 to end
		for(i = LU.size-1; i != p; i--){	// going from end to pivot+1 Optm: If no pivoting ocurred: 1 less cache miss
			LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
		}
		
		/**
		for(i = p+1; i < A.size; i++){   // going from pivot+1 to end
		//for(long i = A.size-1; i >= p+1; i--){	// going from end to pivot+1 Optm: If no pivoting ocurred: 1 less cache miss
			// subtract pivot row U.at(p, _) from current row LU.at(i, _)
			for(j = p+1; j < A.size; j++){
				LU.at(i, j) = LU.at(i, j) - LU.at(p, j) * LU.at(i, p);
				// mulitply pivot line value to multiplier
			}
		}
		/**/
		// TEST
		if((p+1) % bstep == 0) {
			gauss_timer.start();
			/* Streaming *
			for(i = p+1; i < A.size; i++)
				for(j = p+1; j < A.size; j++)
					LU.at(i, j) = LU.at(i, j) - LU.at(p, j) * LU.at(i, p);
			/* Blocking */
			for (bi = p+1; bi < size ; bi += bstep)
				for(j = p+1; j != size; j += 1)
					for(int c = 0; c != bstep; c += 1)
						LU.at(bi+c, j) = LU.at(bi+c, j) - LU.at(p, j) * LU.at(bi+c, p);
			/**/
			double el = gauss_timer.elapsed();
			if((size-(p+1))>0)
				total += el;
		}
	}
	average += total;
}


}
#endif