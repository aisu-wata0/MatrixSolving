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
	vector<double> aux;

	Matrix(long size, int type) : size(size){
		if (type == 0){
			matrix.resize(size*size);
			aux.resize(size);
		}
	}

	double& at (long i, long j) {
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

void subst(Matrix& coefs, vector<double>& indTerms, vector<double>& var_x, bool forward) {
	double sum;
	long i, j;
	int step;
	long size = coefs.size;

	if (forward){
		i = 1;
		step = +1;
		var_x[0] = indTerms[0] / coefs.at(0, 0);
	} else {
		i = size-2;
		step = -1;
		var_x[size-1] = indTerms[size-1] / coefs.at(size-1, size-1);
	}

	for (; i >= 0 && i <= size-1 ; i += step) {
		sum = indTerms[i];
		if(forward) {j = 0;} else {j = size-1;}
		for (; j != i; j += step) {
			sum -= var_x[j] * coefs.at(i, j);
		}
		var_x[i] = sum / coefs.at(i, i);
	}
}
/* */

/**/
void printm(Matrix& matrix){
	for(long i = 0; i <= matrix.size-1; i++){
		for(long j = 0; j <= matrix.size-1; j++)
			cout<< matrix.at(i,j) <<'\t';
		cout<< '\n';
	}
}
void printv(vector<double>& x){
	for(double& vx : x)
		cout << vx <<'\t';
}
#include <chrono>
#include <thread>
/**/
void GaussEl(Matrix& coefs, Matrix& L, Matrix& U, vector<double>& P) {
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
				// subtract pivot
				for (long k = 0; k <= size-1; k++) {
					U.at(i, k) -= U.at(j, k) * L.at(i, j);
				}
			} else {
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
/* */

/* */
int main() {
	fstream file;
	file.open("in.txt");

	long size,i,j,p;
	double sum;

	cout<<"Enter the order of matrix ! ";
	file>> size;

	// double a[10][10],l[10][10]={0},u[10][10]={0},sum,b[10],z[10]={0},x[10]={0};
	Matrix coef(size, 0), l(size, 0), u(size, 0);
	vector<double> b(size), z(size), x(size), P(size);

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
	//double auxx = b.at()
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

	/* Displaying LU matrix */
	cout<< endl << endl <<"LU matrix is "<< endl;
	printm(l);
	cout<< endl;
	printm(u);
	cout<< endl;
	printv(P);

	/* FINDING Z; LZ=b */
	for(i=0; i<=size-1; i++){
		//forward subtitution method
		sum = 0;
		for(p=0; p<i; p++)
		sum += l.at(i,p)*z[p];
		z[i] = (b.at(P.at(i))-sum)/l.at(i,i);
	}

	cout<<endl<<"Z is"<<endl;
	printv(z);

	cout << endl;
	/* FINDING X; UX=Z */
	for(i=size-1; i>=0; i--){
		sum=0;
		for(p=size-1; p>=i; p--)
			sum += u.at(i,p)*x[p];
		x[i] = (z[i]-sum)/u.at(i,i);
	}

	/* DISPLAYING SOLUTION */
	cout<<endl<<"Set of solution is"<<endl;
	printv(x);
	cout << endl;

	return 0;
}
/* */
/* *
int main(int argc, char **argv)
{
	fstream file;
	file.open("in.txt");

	long size;

	file >> size;

	vector<double> indTerms(size);
	vector<double> var_x(size);

	Matrix coefs(size, 0);

	for(long i = 0; i <= size-1; i++)
		for(long j = 0;j <= size-1; j++)
			file >> coefs.at(i, j);

	for(long i = 0; i <= size-1; i++)
		file >> indTerms[i];

	coefs.print(indTerms);

	subst(coefs, indTerms, var_x, false);

	for(long i = 0; i <= size-1; i++)
		cout<< var_x[i] <<'\t';
	cout<< endl;

	return 0;
}
/* */
