#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctgmath>

#include "Timer.h"
#include "Matrix.h"

using namespace std;

Timer timer;
double total_time_iter = 0.0;
double total_time_residue = 0.0;

/**
 * @brief For the matrix LU finds its LU decomposition overwriting it
 * Has partial pivoting, stores final indexes in P
 * @param LU Matrix to be decomposed
 * Output: lower triangle of this matrix will store L 1 diagonal implicit, upper triangle stores U
 * @param P Permutation vector resulting of the pivoting
 */
void GaussEl(Matrix& LU, vector<long>& P) {
	// initializing permutation vector
	for(long i = 0; i <= LU.size-1; i++){
		P.at(i) = i;
	}

	for(long p = 0; p <= LU.size-1; p++){
		// for each pivot
		/* partial pivoting */
		long maxRow = p;
		for(long i = p+1; i <= LU.size-1; i++){
			// for each value below the p pivot
			if(abs(LU.at(i,p)) > abs(LU.at(maxRow,p))) maxRow = i;
		} // finds max value
		// pivots rows of U
		LU.swap_rows_from(p, maxRow, p);
		swap(P.at(p), P.at(maxRow));

		// LU.at(p,p) = 1; implicit
		if(close_zero(LU.at(p,p))){
			fprintf(stderr, "Found a pivot == 0, system is not solvable with partial pivoting");
			exit(1);
		}
		//for (long i = p+1; i <= LU.size-1; i++) {	//going from below pivot to end
		for (long i = LU.size-1; i >= p+1; i--) {	//going from end to pivot
			// for each line below pivot
			if (!close_zero(LU.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
				// subtract pivot from current line (in U)
				for (long k = p+1; k <= LU.size-1; k++) {
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

void subst_P(Matrix& L, vector<double>& X, Matrix& I, bool forward, vector<long>& P, bool unit_diagonal, long col){
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
		sum = I.at(i, col);
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
 * @brief Subtitution method for linear sistems, forward or backward, uses P from pivoting on the LU decomposition
 * @param A Coeficient Matrix
 * @param X Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 * @param P from LU Decomposition
 * @param unit_diagonal true if diagonal is equal to 1
 */
void subst_P(Matrix& A, vector<double>& X, vector<double>& B, bool forward, vector<long>& P, bool unit_diagonal) {
	double sum;
	long i, j;
	int step;
	long size = A.size;

	if (forward){
		i = 1;
		step = +1;
		X.at(0) = B.at(P.at(0));
		if(!unit_diagonal){
			X.at(i) /= A.at(0, 0);
		}
	} else {
		i = size-2;
		step = -1;
		X.at(size-1) = B.at(P.at(size-1));
		if(!unit_diagonal){
			X.at(i) /= A.at(size-1, size-1);
		}
	}

	for(; i >= 0 && i <= size-1 ; i += step){
		sum = B.at(P.at(i));
		if(forward) {j = 0;} else {j = size-1;}
		for(; j != i; j += step){
			// from the start to diagonal, subst values and sub them from the sum
			sum -= X.at(j) * A.at(i, j);
		}
		X.at(i) = sum;
		if(!unit_diagonal){
			X.at(i) /= A.at(i, i);
		}
	}
}

void subst(Matrix& A, Matrix& X, vector<double>& B, bool forward, long col) {
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

/**
 * @brief Subtitution method for linear sistems, forward or backward
 * @param A Coeficient Matrix
 * @param X Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 */
void subst(Matrix& A, vector<double>& X, vector<double>& B, bool forward) {
	double sum;
	long i, j;
	int step;
	long size = A.size;

	if (forward){
		i = 1;
		step = +1;
		X.at(0) = B.at(0) / A.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		X.at(size-1) = B.at(size-1) / A.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = B.at(i);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= X.at(j) * A.at(i, j);
		}
		X.at(i) = sum / A.at(i, i);
	}
}

void solve_lu(Matrix LU, Matrix& X, Matrix& B, vector<long>& P, long col){
	vector<double> Z(LU.size);
	// find Z; LZ=B
	subst_P(LU, Z, B, true, P, true, col);

	// find X; Ux=Z
	subst(LU, X, Z, false, col);
}

void solve_lu(Matrix LU, vector<double>& X, vector<double>& B, vector<long>& P){
	vector<double> Z(X.size());
	// find Z; LZ=B
	subst_P(LU, Z, B, true, P, true);

	// find X; Ux=Z
	subst(LU, X, Z, false);
}

void inverse(Matrix& LU, Matrix& IA, Matrix& I, vector<long>& P){
	int j;
	long size = LU.size;
	vector<double> X(size), B(size), lhs(size);

	for(j = 0; j < size; j++){
		// for each IA col solve SL to find the IA col values
		solve_lu(LU, IA, I, P, j);
	}
}

void JacobiIt(Matrix& A, vector<double>& indTerms, vector<double>& X){
	int total;
	vector<double> new_x;
	// for each line
	for(long i = 0; i <= A.size; ++i){
		total = 0;
		// substitute each X
		for(long j = 0; j <= X.size(); j++) {
			if(j != i){
				total = total + A.at(i, j) * X.at(j);
			}
		}
		new_x.at(i) = (indTerms.at(i) - total) / A.at(i, i);
	}
}

template <class T>
vector<T> add_vec(vector<T> a, vector<T> b, int sign = 1){
	vector<T> x(a.size());

	for(long i=0; i <= ((long)x.size() - 1); i++){
		x.at(i) = a.at(i) + sign*b.at(i);
	}

	return x;
}

vector<double> lhs_value(Matrix A, vector<double> X){
	vector<double> lhs(X.size());

	for(long i=0; i <= X.size()-1; i++){
		lhs.at(i) = 0;
		for(long j=0; j <= X.size()-1; j++){
			lhs.at(i) += A.at(i, j)*X.at(j);
		}
	}
	return lhs;
}

void lu_refining(Matrix& A, Matrix& LU, vector<double>& X, vector<double>& B, vector<long>& P){
	vector<double> lhs(B.size()), r(B.size()), w(B.size());

	cout<<"\n\nRefining "<< endl;
	for(long i=0; i <= 3; i++) {
		lhs = lhs_value(A, X);
		// r: residue of X
		r = add_vec(B, lhs, -1);
		cout<<"\n\n System Residue "<< i << endl;
		printv(r);

		// w: residues of each variable of X
		solve_lu(LU, w, r, P);

		cout<<"\n Variables Residue "<< i << endl;
		printv(w);
		
		// adjust X with found errors
		X = add_vec(X, w);
	}
}

/**
 * @brief Calculates residue into I, A*IA supposed to be Identity
 * @param A original matrix
 * @param IA solution to A*IA = I
 * @param I return value, no init needed
 * @return Max residue value found in abs()
 */
double residue(Matrix& A, Matrix& IA, Matrix& I){
	double norm = 0.0;
	
	for(long icol=0; icol <= A.size-1; icol++){
		// for each column of the inverse
		for(long i=0; i <= IA.size-1; i++){
			// for each line of A
			I.at(i,icol) = 0;
			// multiply A line to the current inverse col
			for(long j=0; j <= IA.size-1; j++){
				I.at(i,icol) -= A.at(i,j)*IA.at(j,icol);
			}
			if(i == icol)
				I.at(i,icol) += 1;
			norm += I.at(i,icol)*I.at(i,icol);
		}
	}
	return sqrt(norm);
}

/**
 * @brief Calculates inverse of A into IA
 * @param LU decomposition of A
 * @param IA return value, no init needed
 * @param P LU pivot permutation
 */
void inverse_refining(Matrix& A, Matrix& LU, Matrix& IA, vector<long>& P, long k){
	long i=0;
	long digits = k > 0 ? (long) log10((double) k) + 1 : 1;
	double c_residue, l_residue;
//	Matrix W(A.size, 0), R(A.size, 0), I(A.size, 0);
	// Optm: iterating line by line
	Matrix W(A.size, BY_LINE), R(A.size, BY_LINE), I(A.size, BY_LINE);

	identity(I);
	// TODO: inverse_id(LU, IA, P); that doesn't need Identity matrix in the memory
	inverse(LU, IA, I, P);
	c_residue = residue(A, IA, R);
	l_residue = c_residue*2; // enter condition at least once
	
	cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	while(i < k){
		// (abs(l_residue - c_residue)/c_residue > EPSILON) && (l_residue > c_residue)
		// relative approximate error
		i += 1;
		// R: residue of IA
		
		timer.start();
		inverse(LU, W, R, P);
		// W: residues of each variable of IA
		// adjust IA with found errors
		IA.add(W);	//TODO: Optm: add at the same time its calculating W
		total_time_iter += timer.elapsed();

		l_residue = c_residue;
		
		timer.start();
		c_residue = residue(A, IA, R);
		total_time_residue += timer.elapsed();
		cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	}
}

int main(int argc, char **argv) {
	cout.precision(17);
	cout << scientific;
	srand(20172);

	long size, k = 64;
	
	// if(has_in_file_argument)
		ifstream in_f(IODIR "in.txt");
		cin.rdbuf(in_f.rdbuf()); //redirect

	// if(has_out_file_argument)
		ofstream o_f(IODIR "out.txt");
		cout.rdbuf(o_f.rdbuf()); //redirect

//	read matrix size
	cin>> size;

	Matrix A(size, BY_COL), LU(size, BY_COL);
	// Matrix IA(size, BY_COL);
	// Optm: iterating line by line
	Matrix IA(size, BY_LINE);
	vector<double> B(size), z(size);
	vector<long> P(size);
	vector<double> X;
	
//	read matrix coef
	/**/
	for(long i=0; i<=size-1; i++){
		for(long j=0; j<=size-1; j++){
			cin>> A.at(i,j);
			LU.at(i,j) = A.at(i,j);
		}
	}
	/**
	read B values
	for(i=0; i<=size-1; i++)
		cin>> B.at(i);
	/**/
	timer.start();
	double lu_time = 0.0;
	
	GaussEl(LU, P);
	lu_time = timer.elapsed();

	cout<<"#\n";
	inverse_refining(A, LU, IA, P, k);

	
	cout<< defaultfloat;
	cout<<"# Tempo LU: "<< lu_time <<"\n"; // TODO
	cout<<"# Tempo iter: "<< total_time_iter/(double)k <<"\n";
	cout<<"# Tempo residuo: "<< total_time_residue/(double)k <<"\n#\n";
	cout << scientific;
	printm(IA);
	/**
	find first iteration of X
	solve_lu(LU, X, B, P);

	lu_refining(A, LU, X, B, P);

	cout<<"\nSet of solution is"<<endl;
	printv(X);
	cout << endl;
	/**/

	/**
	ofstream o_f;
	cout.precision(17);
	cout << scientific;
	for(int n = 2; n <= 2048; n <<= 1){
		o_f.open(IODIR + to_string(n) + ".txt");
		cout.rdbuf(o_f.rdbuf()); //redirect
		
		generateSquareRandomMatrix(n);
		
		o_f.close();
	}
	/**/
	
	return 0;
}
