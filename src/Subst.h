#ifndef SUBST_H
#define SUBST_H
/**
@file Subst.h
*/
namespace std {

/**
 * @brief Subtitution method for linear sistems, forward or backward, uses P from pivoting on the LU decomposition
 * @param L Coeficient Matrix
 * @param X Variables to be found
 * @param I Independant terms
 * @param forward true is forward subst, false is backward.
 * @param P from LU Decomposition
 * @param unit_diagonal true if diagonal is equal to 1
 * @param col Column of the matrix I to be used as B
 */
void subst_P(Matrix& L, vector<double>& X, MatrixColMajor& I, bool forward, vector<long>& P, bool unit_diagonal, long col){
	double sum;
	long i, j;
	int step;
	long size = L.size;

	if (forward){
		i = 1;
		step = +1;
		X.at(0) = I.at(P.at(0),col);
		if(!unit_diagonal){
			X.at(i) /= L.at(0, 0);
		}
	} else {
		i = size-2;
		step = -1;
		X.at(size-1) = I.at(P.at(size-1),col);
		if(!unit_diagonal){
			X.at(i) /= L.at(size-1, size-1);
		}
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = I.at(P.at(i), col);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j) * L.at(i, j);
		}
		X.at(i) = sum;
		if(!unit_diagonal){
			X.at(i) /= L.at(i, i);
		}
	}
}

/**
 * @brief Subtitution method for linear sistems, forward or backward
 * @param A Coeficient Matrix
 * @param X Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 * @param col Column of the matrix to be used as B
 */
void subst(Matrix& A, MatrixColMajor& X, vector<double>& B, bool forward, long col) {
	double sum;
	long i, j;
	int step;
	long size = A.size;

	if (forward){
		i = 1;
		step = +1;
		X.at(0,col) = B.at(0) / A.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		X.at(size-1,col) = B.at(size-1) / A.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = B.at(i);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j,col) * A.at(i, j);
		}
		X.at(i,col) = sum / A.at(i, i);
	}
}

}
#endif