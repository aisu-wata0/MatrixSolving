#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

#include "TestTimer.cpp"

using namespace std;

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

	void swap_rows(long row0, long row1){
		if(row0 == row1)
			return;
		for(long j = 0; j <= size-1; j++){
			swap(this->at(row0, j), this->at(row1, j));
		}
	}

	void print(){
		for (long i = 0; i < size; i++) {
			for (long j = 0; j < size; j++) {
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

void subst(Matrix& A, vector<double>& var_x, vector<double>& indTerms, bool forward, vector<long>& perm) {
	double sum;
	long i, j;
	int step;
	long size = A.size;

	if (forward){
		i = 1;
		step = +1;
		var_x.at(0) = indTerms.at(perm.at(0)) / A.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x.at(size-1) = indTerms.at(perm.at(size-1)) / A.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = indTerms.at(perm.at(i));
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= var_x.at(j) * A.at(i, j);
		}
		var_x.at(i) = sum / A.at(i, i);
	}
}
/* */
void subst(Matrix& A, vector<double>& var_x, vector<double>& indTerms, bool forward) {
	double sum;
	long i, j;
	int step;
	long size = A.size;

	if (forward){
		i = 1;
		step = +1;
		var_x.at(0) = indTerms.at(0) / A.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x.at(size-1) = indTerms.at(size-1) / A.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = indTerms.at(i);
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= var_x.at(j) * A.at(i, j);
		}
		var_x.at(i) = sum / A.at(i, i);
	}
}
/**/
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

/**/
void GaussEl(Matrix& A, Matrix& L, Matrix& U, vector<long>& P) {
	TestTimer t("GaussEl");
	long size = A.size;

	for(long i = 0; i <= size-1; i++){
		P.at(i) = i;
		for(long j = 0; j <= size-1; j++){
			U.at(i, j) = A.at(i,j);
		}
	}

	for(long j = 0; j <= size-1; j++){
		// for each pivot
		/* partial pivoting */
		long maxRow = j;
		for(long p = j+1; p <= size-1; p++){
			if(abs(U.at(p,j)) > abs(U.at(maxRow,j))) maxRow = p;
		}
		U.swap_rows(j, maxRow);
		swap(P.at(j), P.at(maxRow));
		/**/
		L.at(j,j) = 1;
		for (long i = j+1; i <= size-1; i++) {
			// for each line below pivot
			if (U.at(i, j) != 0){ // TODO: comparison with 0
				// find pivot multiplier
				L.at(i, j) = U.at(i, j)/U.at(j, j);
				// subtract pivot from current line
				for (long k = 0; k <= size-1; k++) {
					U.at(i, k) -= U.at(j, k) * L.at(i, j);
				}
			} else {
				cout<<"Found a pivot == 0, system is not solvable with partial pivoting"<< endl;
				L.at(i,j) = 0;
			}
		}
	}
}

void solve_lu(Matrix l, Matrix u, vector<double>& x, vector<double>& b, vector<long>& P){
	vector<double> z(x.size());
	//* find z; LZ=b *
    subst(l, z, b, true, P);
	//* find x; Ux=z *
    subst(u, x, z, false);
}

void inverse(Matrix& I, Matrix A, Matrix L, Matrix U, vector<long>& P){  
    int i, j, n;  
    long size = L.size;
    vector<double> x(size), b(size);
    
    //Creating identity matrix
    for(i = 0; i < size; i++){  
        for(j = 0; j < size; j++) I.at(i, j) = 0;  
        I.at(i, i) = 1;  
    }
    
    for(i = 0; i < size; i++){          //for each row
        for(j = 0; j < size; j++){      //copying row elements
            x.at(j) = A.at(j, i);      	//from the original matrix
            b.at(j) = I.at(j, i);       //from identity
        }
        solve_lu(L, U, x, b, P);        //solve this row system
        for(n = 0; n < size; n++) I.at(n, i) = x.at(n);        //Copying solution to inverse
    }
   
}  
   
/* */
void Cholesky(Matrix& a, Matrix& l, Matrix& u) {
	long i, j, p, size = a.size;
	double sum, sum2;
	/**/
	for(i=0; i <= size-1; i++){
		sum = 0;
		for(j=0; j <= i-1; j++){
			sum2 = 0;
			for(p = 0; p <= j-1; p++)
				sum2 += l.at(i,p) * l.at(j,p);
			l.at(i,j) = (a.at(i,j) - sum2) / l.at(j,j);
			u.at(j,i) = l.at(i,j);

			sum += l.at(i,j)*l.at(i,j);
		}
		l.at(i,i) = sqrt(a.at(i,i) - sum);
		u.at(i,i) = l.at(i,i);
	}
	/**/
}

////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_Decomposition(double *A, int n)                              //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to decompose the n x n matrix A       //
//     into a lower triangular matrix L and a unit upper triangular matrix U  //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     This routine is suitable for those classes of matrices which when      //
//     performing Gaussian elimination do not need to undergo partial         //
//     pivoting, e.g. positive definite symmetric matrices, diagonally        //
//     dominant band matrices, etc.                                           //
//     For the more general case in which partial pivoting is needed use      //
//                    Crout_LU_Decomposition_with_Pivoting.                   //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is Crout_LU_Solve          //
//     (see below).                                                           //
//                                                                            //
//     The Crout method is given by evaluating, in order, the following       //
//     pair of expressions for k = 0, ... , n - 1:                            //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
/* */
int crout(Matrix& a, Matrix& l, Matrix& u) {
	long k, i, j, p, size = a.size;
	double sum;

	/**
	for(k = 0; k <= size-1; k++){
		u.at(k,k) = 1;
		for(i = k; i <= size-1; i++){
			sum = 0;
			for(p = 0; p <= k-1; p++){
				sum += l.at(i,p)*u.at(p,k);
			}
			l.at(i,k) = a.at(i,k)-sum;
		}

		for(j = k+1; j <= size-1 ;j++){
			sum = 0;
			for(p = 0; p <= k-1; p++){
				sum += l.at(k,p)*u.at(p,j);
			}
			u.at(k,j) = (a.at(k,j)-sum)/l.at(k,k);
		}
	}
	/**
	for(k=0; k <= size-1; k++){
		for(i=0; i <= size-1; i++){
			if(i<=k){
				u.at(i,k)=a.at(i,k);
				for(p=0; p<i-1; p++)
				  u.at(i,k)-=l.at(i,p)*u.at(p,k);
				if(i==k)
				  l.at(i,k)=1;
				else
				  l.at(i,k)=0;
			}else{
				l.at(i,k)=a.at(i,k);
				for(p=0; p<=k-1; p++)
				  l.at(i,k)-=l.at(i,p)*u.at(p,k);
				l.at(i,k)/=u.at(k,k);
				u.at(i,k)=0;
			}
		}
	}
	/*For each row and column, k = 0, ..., n-1,
	find the lower triangular matrix elements for column k
	and if the matrix is non-singular (nonzero diagonal element).
	find the upper triangular matrix elements for row k. */
	/**
	int row;
	double *p_k, *p_row, *p_col, *A = &a.matrix[0];
	long n = a.size;
	
	for (k = 0, p_k = A; k < n; p_k += n, k++) {
		for (i = k, p_row = p_k; i < n; p_row += n, i++) {
			for (p = 0, p_col = A; p < k; p_col += n, p++)
				*(p_row + k) -= *(p_row + p) * *(p_col + k);
		}  
		if ( *(p_k + k) == 0.0 ) return -1;
		for (j = k+1; j < n; j++) {
			for (p = 0, p_col = A; p < k; p_col += n,  p++)
				*(p_k + j) -= *(p_k + p) * *(p_col + j);
			*(p_k + j) /= *(p_k + k);
		}
	}
	/* *
	for (k = 0; k <= size-1; k++){
		for (i = 0; i <= size-1; i++){
			if (i < k){
				l.at(i,k) = 0;
			}else{
				l.at(i,k) = a.at(i,k);
				for (p = 0; p < k; p++){
					l.at(i,k) = l.at(i,k) - l.at(i,p) * u.at(p,k);
				}
			}
		}
		for (i = 0; i <= size-1; i++){
			if (i < k){
				u.at(k,i) = 0;
			}else if (i == k){
				u.at(k,i) = 1;
			}else{
				u.at(k,i) = a.at(k,i) / l.at(k,k);
				for (p = 0; p < k; p++){
					u.at(k,i) = u.at(k,i) - ((l.at(k,p) * u.at(p,i)) / l.at(k,k));
				}
			}
	  }
	}
	/**/
	return 0;
}

void JacobiIt(Matrix& A, vector<double>& indTerms, vector<double>& x){
    int total;
    vector<double> new_x;
    // for each line
    for(long i = 0; i <= A.size; ++i){
        total = 0;
        // substitute each x
        for(long j = 0; j <= x.size(); j++) {
            if(j != i){
                total = total + A.at(i, j) * x.at(j);
            }
        }
        new_x.at(i) = (indTerms.at(i) - total) / A.at(i, i);
    }
}

vector<double> found_values(Matrix A, vector<double> x){
    vector<double> values(x.size());

    for(long i=0; i <= x.size()-1; i++){
        values.at(i) = 0;
        for(long j=0; j <= x.size()-1; j++){
            values.at(i) += A.at(i, j)*x.at(j);
        }
    }
    return values;
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

int main(int argc, char **argv) {
	fstream file;
	file.open("in.txt");

	long size,i,j;
    
	cout<<"Enter the order of matrix ! ";
	file>> size;

	Matrix coef(size, 0), l(size, 0), u(size, 0), I(size, 0);
	vector<double> b(size), z(size);
	vector<long> P(size);
	vector<vector<double>> x;
	x.push_back(vector<double>(size));

	cout<<"Enter all coefficients of matrix : ";
	for(i=0; i<=size-1; i++){
		cout<<"\nRow "<<i<<"  ";
		for(j=0; j<=size-1; j++)
			file>> coef.at(i,j);
	}
	cout<<"Enter elements of b matrix"<<endl;
	for(i=0; i<=size-1; i++)
		file>> b.at(i);

	//crout(coef, l, u);
	//Cholesky(coef, l, u);
	GaussEl(coef, l, u, P);
    inverse(I, coef, l, u, P);
	//* Displaying LU matrix *
	cout<<"\n\nL matrix is "<< endl;
	printm(l);
	cout<<"\nU matrix is "<< endl;
	printm(u);
    cout<<"\nInv matrix is "<< endl;
	printm(I);
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
