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

void subst(Matrix& coefs, vector<double>& var_x, vector<double>& indTerms, bool forward, vector<long>& perm) {
	double sum;
	long i, j;
	int step;
	long size = coefs.size;
    
	if (forward){
		i = 1;
		step = +1;
		var_x.at(0) = indTerms.at(perm.at(0)) / coefs.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x.at(size-1) = indTerms.at(perm.at(size-1)) / coefs.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = indTerms.at(perm.at(i));
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= var_x.at(j) * coefs.at(i, j);
		}
		var_x.at(i) = sum / coefs.at(i, i);
	}
}
/* */
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
void GaussEl(Matrix& coefs, Matrix& L, Matrix& U, vector<long>& P) {
	TestTimer t("GaussEl");
	long size = coefs.size;

	for(long i = 0; i <= size-1; i++){
		P.at(i) = i;
		for(long j = 0; j <= size-1; j++){
			U.at(i, j) = coefs.at(i,j);
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
/* */
void crout(Matrix& a, Matrix& l, Matrix& u) {
	long k, i, j, p, size = a.size;
	double sum;
	/**/
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
	/**/
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
	//* find z; LZ=b *
    subst(l, z, b, true, P);

	cout<< endl <<"Z is"<< endl;
	printv(z);
	cout<< endl;
	
    subst(u, x, z, false);
}

template <class T>
vector<T> add_vec(vector<T> a, vector<T> b, int sign = 1){
    vector<T> x(a > b ? a.size() : b.size());
    
    for(long i=0; i <= x.size(); i++){
        x.at(i) = a.at(i) + sign*b.at(i);
    }
    
    return x;
}

int main(int argc, char **argv) {
	fstream file;
	file.open("in.txt");

	long size,i,j;

	cout<<"Enter the order of matrix ! ";
	file>> size;

	Matrix coef(size, 0), l(size, 0), u(size, 0);
	vector<double> b(size), z(size);
	vector<long> P(size);
	vector<vector<double>> x(5);
	x.at(0).resize(size);

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
	/* *
	for (i = 0; i <= size-1; i++){
		for (j = 0; j <= size-1; j++){
			if (j < i){
				l.at(j,i) = 0;
			}else{
				l.at(j,i) = coef.at(j,i);
				for (k = 0; k < i; k++){
					l.at(j,i) = l.at(j,i) - l.at(j,k) * u.at(k,i);
				}
			}
		}
		for (j = 0; j <= size-1; j++){
			if (j < i){
				u.at(i,j) = 0;
			}else if (j == i){
				u.at(i,j) = 1;
			}else{
				u.at(i,j) = coef.at(i,j) / l.at(i,i);
				for (k = 0; k < i; k++){
					u.at(i,j) = u.at(i,j) - ((l.at(i,k) * u.at(k,j)) / l.at(i,i));
				}
			}
	  }
	}
	/**/

	//* Displaying LU matrix *
	cout<< endl << endl <<"L matrix is "<< endl;
	printm(l);
	cout<< endl <<"U matrix is "<< endl;
	printm(u);
	cout<< endl <<"Pivoting Permutaton is "<< endl;
	printv(P);

	//* find first iteration of x *
	solve_lu(l, u, x.at(0), b, P);
	
	vector<double> values(x.size()), r(x.size()), w(x.size());
	for(long i=1; i <= 3 ; i++) {
    	values = found_values(coef, x.at(i));
    	// r: residue of x
    	r = add_vec(b, values, -1);
    	
    	// w: residues of each variable of x
    	solve_lu(l, u, w, r, P);
    	
    	x.at(i+1).resize(size);
    	// adjust x with found errors
    	x.at(i+1) = add_vec(x.at(i), w);
	}
	

	cout<<endl<<"Set of solution is"<<endl;
	//printv(x.at(x.size()-1));
	cout << endl;

	return 0;
}