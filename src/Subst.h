#ifndef SUBST_H
#define SUBST_H
/**
@file Subst.h
*/
namespace std {

/**
 * @brief Subtitution method for linear sistems, forward or backward, uses P from pivoting on the LU decomposition
 * @param T Triangular coeficient Matrix
 * @param X Variables to be found
 * @param I Independant terms
 * @param forward true is forward subst, false is backward.
 * @param P from LU Decomposition
 * @param unit_diagonal true if diagonal is equal to 1
 * @param col Column of the matrix I to be used as B
 */
template<bool forward, bool unit_diagonal, class TMatrix>
void subst_P(TMatrix& T, vector<double>& X, MatrixColMajor& I, vector<long>& P, long col){
	double sum;
	long i, j;
	int step;
	long size = T.size;
	
	if(forward){
		i = 1;
		step = +1;
		X.at(0) = I.at(P.at(0),col);
		if(not unit_diagonal){
			X.at(i) /= T.at(0, 0);
		}
	} else {
		i = size-2;
		step = -1;
		X.at(size-1) = I.at(P.at(size-1),col);
		if(not unit_diagonal){
			X.at(i) /= T.at(size-1, size-1);
		}
	}

	for (; i >= 0 && i < size; i += step) {
		sum = I.at(P.at(i), col);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j) * T.at(i, j);
		}
		X.at(i) = sum;
		if(not unit_diagonal){
			X.at(i) /= T.at(i, i);
		}
	}
}
/**
 * @brief Subtitution method for linear sistems, forward or backward
 * @param T Triangular coeficient Matrix
 * @param X Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 * @param col Column of the matrix to be used as B
 */
template<bool forward, class TMatrix>
void subst(TMatrix& T, MatrixColMajor& X, vector<double>& B, long col) {
	double sum;
	long i, j;
	int step;
	long size = T.size;

	if(forward){
		i = 1;
		step = +1;
		X.at(0,col) = B.at(0) / T.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		X.at(size-1,col) = B.at(size-1) / T.at(size-1, size-1);
	}

	for (; i >= 0 && i < size; i += step) {
		sum = B.at(i);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j,col) * T.at(i, j);
		}
		X.at(i,col) = sum / T.at(i, i);
	}
}


}
#endif