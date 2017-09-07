#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

#include "TestTimer.cpp"

using namespace std;

#define EPSILON 1e-7

bool close_zero(double x){
	return fabs(x) < EPSILON;
}

class Matrix
{
public:
	long size;
	vector<double> matrix;

	Matrix(long size, int type) : size(size){
		if (type == 0){
			matrix.resize(size*size);
		}
	}

	double& at(long i, long j) {
		return matrix.at(i*size + j);
	}
	
	void swap_rows_from(long row0, long row1, long start){
		if(row0 == row1)
			return;
		for(long j = start; j <= size-1; j++){
			// for each collumn
			swap(this->at(row0, j), this->at(row1, j));
		}
	}
	
	void swap_rows(long row0, long row1){
		swap_rows_from(row0, row1, 0);
	}
	
	void print(){
		for(long i = 0; i < size; i++){
			for(long j = 0; j < size; j++){
				cout << this->at(i, j) <<'\t';
			}
			cout << endl;
		}
	}

	void print(vector<double>& indTerms){
		for (long i = 0; i < size; i++) {
			for (long j = 0; j < size; j++) {
				cout << this->at(i, j) <<"x"<< j <<'\t';
			}
			cout <<"=\t"<< indTerms[i] << endl;
		}
	}
};

void subst(Matrix& coefs, vector<double>& var_x, vector<double>& indTerms, bool forward, vector<long>& P) {
	double sum;
	long i, j;
	int step;
	long size = coefs.size;

	if (forward){
		i = 1;
		step = +1;
		var_x.at(0) = indTerms.at(P.at(0)) / coefs.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x.at(size-1) = indTerms.at(P.at(size-1)) / coefs.at(size-1, size-1);
	}
	
	for(k = 0; k <= size-1; k++){
		if()
		swap(indTerms.at(), indTerms.at());
	}

	for(; i >= 0 && i <= size-1 ; i += step){
		sum = indTerms.at(P.at(i));
		if(forward) {j = 0;} else {j = size-1;}
		for(; j != i; j += step){
			sum -= var_x.at(j) * coefs.at(i, j);
		}
		var_x.at(i) = sum / coefs.at(i, i);
	}
}

/**
 * @brief Subtitution method for linear sistems, forward or backward
 * @param A Coeficient Matrix
 * @param x Variables to be found
 * @param B Independant terms
 * @param forward true is forward subst, false is backward.
 */
void subst(Matrix& coefs, vector<double>& var_x, vector<double>& indTerms, bool forward) {
	double sum;
	long i, j;
	int step;
	long size = coefs.size;

	if (forward){
		i = 1;
		step = +1;
		var_x.at(0) = indTerms.at(0) / coefs.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x.at(size-1) = indTerms.at(size-1) / coefs.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = indTerms.at(i);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= var_x.at(j) * coefs.at(i, j);
		}
		var_x.at(i) = sum / coefs.at(i, i);
	}
}

void printm(Matrix& matrix){
	for(long i = 0; i <= matrix.size-1; i++){
		for(long j = 0; j <= matrix.size-1; j++)
			cout<< matrix.at(i,j) <<'\t';
		cout<< '\n';
	}
}

template <class T>
void printv(vector<T>& x){
	for(T& vx : x)
		cout << vx <<'\t';
}

/**
 * @brief For the matrix A finds its LU decomposition without overwriting.
 * Has partial pivoting, stores final indexes in P
 * @param A Matrix to be decomposed
 * @param LU Decomposition, lower triangle of this matrix will store L without the 1 diagonal, upper triangle stores U
 * @param P Permutation vector resulting of the pivoting
 */
void GaussEl(Matrix& coefs, Matrix& LU, vector<long>& P) {
	TestTimer t("GaussEl");
	long size = coefs.size;

	// initializing permutation vector
	for(long i = 0; i <= size-1; i++){
		P.at(i) = i;
	}

	for(long p = 0; p <= size-1; p++){
		// for each pivot
		/* partial pivoting */
		long maxRow = p;
		for(long i = p+1; i <= size-1; i++){
			// for each value below the p pivot
			if(abs(LU.at(i,p)) > abs(LU.at(maxRow,p))) maxRow = i;
		} // finds max value
		// pivots rows of U
		LU.swap_rows_from(p, maxRow, p);
		swap(P.at(p), P.at(maxRow));
		
		// LU.at(p,p) = 1; // change subst method
		close_zero(LU.at(p,p))
		if(close_zero(LU.at(p,p))){
			cout<<"Found a pivot == 0, system is not solvable with partial pivoting"<< endl;
			exit(1);
		}
		for (long i = p+1; i <= size-1; i++) {
			// for each line below pivot
			if (!close_zero(LU.at(i,p))){
				// only subtract pivot line if coeficient is not null
				// find pivot multiplier, store in L
				LU.at(i, p) = LU.at(i, p)/LU.at(p, p);
				// subtract pivot from current line (in U)
				for (long k = p+1; k <= size-1; k++) {
					// for each collumn stating from pivot's
					LU.at(i, k) -= LU.at(p, k) * LU.at(i, p);
					// mulitply pivot line value to multiplier
				}
			}
		}
	}
}

void JacobiIt(Matrix& coefs, vector<double>& indTerms, vector<double>& x){
    int total;
    vector<double> new_x;
    // for each line
    for(long i = 0; i <= coefs.size; ++i){
        total = 0;
        // substitute each x
        for(long j = 0; j <= x.size(); j++) {
            if(j != i){
                total = total + coefs.at(i, j) * x.at(j);
            }
        }
        new_x.at(i) = (indTerms.at(i) - total) / coefs.at(i, i);
    }
}

vector<double> found_values(Matrix coefs, vector<double> x){
    vector<double> values(x.size());

    for(long i=0; i <= x.size()-1; i++){
        values.at(i) = 0;
        for(long j=0; j <= x.size()-1; j++){
            values.at(i) += coefs.at(i, j)*x.at(j);
        }
    }
    return values;
}

void solve_lu(Matrix l, Matrix u, vector<double>& x, vector<double>& b, vector<long>& P){
	vector<double> z(x.size());
	// find z; LZ=b
    subst(l, z, b, true, P);
	// find x; Ux=z
    subst(u, x, z, false);
}

template <class T>
vector<T> add_vec(vector<T> a, vector<T> b, int sign = 1){
    vector<T> x(a.size());

    for(long i=0; i <= x.size()-1; i++){
        x.at(i) = a.at(i) + sign*b.at(i);
    }

    return x;
}

void lu_refining(Matrix coef, Matrix l, Matrix u, vector<vector<double>>& x, vector<double>& b, vector<long>& P){
	vector<double> values(b.size()), r(b.size()), w(b.size());

	cout<<"\n\nRefining "<< endl;
	for(long i=0; i <= 3; i++) {
    	values = found_values(coef, x.at(i));
    	// r: residue of x
    	r = add_vec(b, values, -1);
		cout<<"\n\n System Residue "<< i << endl;
		printv(r);

    	// w: residues of each variable of x
    	solve_lu(l, u, w, r, P);

		cout<<"\n Variables Residue "<< i << endl;
		printv(w);

		x.push_back(vector<double>(b.size()));
    	// adjust x with found errors
    	x.at(i+1) = add_vec(x.at(i), w);
	}
}

union my_double {
	double f;
	long long exp : 11;
	long long mant : 53;
} my_double_t; 

int main(int argc, char **argv) {
	my_double g; // TODO
	g.f = 1.0;
	g.mant = g.mant + 1;
	cout<<"\nX =  "<<g.f - 1.0<<"  ";
	
	
	fstream file;
	file.open("in.txt");

	long size,i,j;

	cout<<"Enter the order of matrix ! ";
	file>> size;

	Matrix coef(size), lu(size);
	vector<double> b(size), z(size);
	vector<long> P(size);
	vector<vector<double>> x;
	x.push_back(vector<double>(size));

	cout<<"Enter all coefficients of matrix : ";
	for(i=0; i<=size-1; i++){
		cout<<"\nRow "<<i<<"  ";
		for(j=0; j<=size-1; j++)
			file>> coef.at(i,j);
			lu.at(i,j) = coef.at(i,j);
	}
	cout<<"Enter elements of b matrix"<<endl;
	for(i=0; i<=size-1; i++)
		file>> b.at(i);

	//crout(coef, l, u);
	//Cholesky(coef, l, u);
	GaussEl(coef, lu, P);

	//* Displaying LU matrix *
	cout<<"\n\nL matrix is "<< endl;
	printm(l);
	cout<<"\nU matrix is "<< endl;
	printm(u);
	cout<<"\nPivoting Permutaton is "<< endl;
	printv(P);

	//* find first iteration of x *
	solve_lu(l, u, x.at(0), b, P);

	lu_refining(coef, l, u, x, b, P);

	cout<<"\nSet of solution is"<<endl;
	printv(x.at(x.size()-1));
	cout << endl;

	return 0;
}
