#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

#include "TestTimer.cpp"
#include "Matrix.h"

using namespace std;

/**
 * @brief For the matrix LU finds its LU decomposition overwriting it
 * Has partial pivoting, stores final indexes in P
 * @param LU Matrix to be decomposed
 * Output: lower triangle of this matrix will store L 1 diagonal implicit, upper triangle stores U
 * @param P Permutation vector resulting of the pivoting
 */
void GaussEl(Matrix& LU, vector<long>& P) {
	//TestTimer t("GaussEl");

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
		
		// LU.at(p,p) = 1; // change subst method
		close_zero(LU.at(p,p));
		if(close_zero(LU.at(p,p))){
			cout<<"Found a pivot == 0, system is not solvable with partial pivoting"<< endl;
			exit(1);
		}
		for (long i = p+1; i <= LU.size-1; i++) {
			// for each line below pivot
			if (!close_zero(LU.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
				// subtract pivot from current line (in U)
				for (long k = p+1; k <= LU.size-1; k++) {
					// for each collumn stating from pivot's
					LU.at(i, k) -= LU.at(p, k) * LU.at(i, p);
					// mulitply pivot line value to multiplier
				}
			}
		}
	}
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
		X.at(size-1) = B.at(P.at(size-1)) / A.at(size-1, size-1);
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

void solve_lu(Matrix LU, vector<double>& X, vector<double>& B, vector<long>& P){
	vector<double> Z(X.size());
	// find Z; LZ=B
    subst_P(LU, Z, B, true, P, true);
	cout<<"\nZ"<< endl;
	printv(Z);
	// find X; Ux=Z
    subst(LU, X, Z, false);
}

void inverse(Matrix& I, Matrix& A, Matrix& LU, vector<long>& P){  
	int i, j, n;  
	long size = LU.size;
	vector<double> X(size), B(size), lhs(size);

	// Creating identity matrix
	for(i = 0; i < size; i++){  
		for(j = 0; j < size; j++) I.at(i, j) = 0;  
		I.at(i, i) = 1;  
	}

	for(i = 0; i < size; i++){			// for each row
		for(j = 0; j < size; j++){		// copying row elements
			X.at(j) = A.at(j, i);		// from the original matrix
			B.at(j) = I.at(j, i);		// from identity
		}
		solve_lu(LU, X, B, P);		// solve this row system
		
		for(n = 0; n < size; n++)
			I.at(n, i) = X.at(n);
		// Copying solution to inverse
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

    for(long i=0; i <= x.size()-1; i++){
        x.at(i) = a.at(i) + sign*b.at(i);
    }

    return x;
}

void lu_refining(Matrix A, Matrix LU, vector<vector<double>>& X, vector<double>& B, vector<long>& P){
	vector<double> lhs(B.size()), r(B.size()), w(B.size());

	cout<<"\n\nRefining "<< endl;
	for(long i=0; i <= 3; i++) {
    	lhs = lhs_value(A, X.at(i));
    	// r: residue of X
    	r = add_vec(B, lhs, -1);
		cout<<"\n\n System Residue "<< i << endl;
		printv(r);

    	// w: residues of each variable of X
    	solve_lu(LU, w, r, P);

		cout<<"\n Variables Residue "<< i << endl;
		printv(w);

		X.push_back(vector<double>(B.size()));
    	// adjust X with found errors
    	X.at(i+1) = add_vec(X.at(i), w);
	}
}

void inverse_refining(Matrix A, Matrix LU, vector<vector<double>>& X, vector<double>& B, vector<long>& P){
	vector<double> lhs(B.size()), r(B.size()), w(B.size());

	cout<<"\n\nRefining "<< endl;
	for(long i=0; i <= 3; i++) {
    	lhs = lhs_value(A, X.at(i));
    	// r: residue of X
    	r = add_vec(B, lhs, -1);
		cout<<"\n\n System Residue "<< i << endl;
		printv(r);

    	// w: residues of each variable of X
    	solve_lu(LU, w, r, P);

		cout<<"\n Variables Residue "<< i << endl;
		printv(w);

		X.push_back(vector<double>(B.size()));
    	// adjust X with found errors
    	X.at(i+1) = add_vec(X.at(i), w);
	}
}

union my_double {
	double f;
	long long exp : 11;
	long long mant : 53;
} my_double_t; 


//	my_double g; // TODO
//	g.f = 1.0;
//	g.mant = g.mant + 1;
//	cout<<"\nX =  "<<g.f - 1.0<<"  ";

int main(int argc, char **argv) {
	fstream file;
	long size,i,j;
	file.open("in.txt");
    
	cout<<"Enter the order of matrix = ";
	file>> size;

	Matrix A(size, 0), LU(size, 0), I(size, 0);
	vector<double> B(size), z(size);
	vector<long> P(size);
	vector<vector<double>> X;
	X.push_back(vector<double>(size));

	cout<<"Enter all coefficients of matrix : ";
	for(i=0; i<=size-1; i++){
		cout<<"\nRow "<<i<<"  ";
		for(j=0; j<=size-1; j++){
			file>> A.at(i,j);
			LU.at(i,j) = A.at(i,j);
		}
	}
//	cout<<"Enter elements of B matrix"<<endl;
//	for(i=0; i<=size-1; i++)
//		file>> B.at(i);

	GaussEl(LU, P);

	cout<<"\n\nLU matrix is "<< endl;
	printm(LU);
	cout<<"\nPivoting Permutaton is "<< endl;
	printv(P);
	
    inverse(I, A, LU, P);
    cout<<"\nInv matrix is "<< endl;
	printm(I);

	// find first iteration of X
//	solve_lu(LU, X.at(0), B, P);
//
//	lu_refining(A, LU, X, B, P);
//
//	cout<<"\nSet of solution is"<<endl;
//	printv(X.at(X.size()-1));
//	cout << endl;

	return 0;
}
