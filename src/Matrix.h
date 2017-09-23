#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

namespace std {

#define EPSILON 1e-15
#define E_FACTOR 16

union Double {
	double f;
	long long exp : 11;
	long long mant : 53;
}; 

bool close_zero(double x){
	return fabs(x) <= EPSILON;
}

double inc(double x, int i){
	Double y;
	y.f = x;
	y.mant += i;
	return y.f;
}

/**
 * @brief  uses a factor of the epsilon of 'a' to compare if b is in it's range
 * @return true if they are close enough,
 */
bool near(double a, double b){
	double min_a = a - (a - inc(a, -1)) * E_FACTOR;
	double max_a = a + (inc(a, +1) - a) * E_FACTOR;

	return (min_a <= b && max_a >= b);
}

enum Matrix_type { BY_COL, BY_LINE };

/**
 * @brief Stores values from specified type of matrix, can be acessed using conventional indexing .at(i,j)
 */
class Matrix
{
public:
	long size;
	Matrix_type type;
	vector<double> matrix;

	Matrix(long size, Matrix_type type) : size(size), type(type), matrix(size*size){
	}

	double& at(long i, long j) {
		if(type == BY_COL){
			return matrix.at(i*size + j);
		} else {
			return matrix.at(j*size + i);
		}
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
	
	void add(Matrix& b, double sign = 1){
		for(long i=0; i <= size-1; i++){
			for(long j=0; j <= size-1; j++){
				this->at(i,j) += sign*b.at(i,j);
			}
		}
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

void identity(Matrix& I){
	for(long i = 0; i < I.size; i++){
		for(long j = 0; j < I.size; j++){
			I.at(i, j) = 0;
		}
		I.at(i, i) = 1;
	}
}

void generateSquareRandomMatrix(long n){
	long i, j;
	double invRandMax = 1.0/(double)RAND_MAX;
	
	cout<< n <<"\n";
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			cout<< (double)rand() * invRandMax <<"\t";
		}
		cout<<"\n";
	}
}

void randomMatrix(Matrix& M, Matrix& LU){
	long i, j;
	double invRandMax = 1.0/(double)RAND_MAX;
	
	for(i = 0; i < M.size; i++){
		for(j = 0; j < M.size; j++){
			M.at(i,j) = (double)rand() * invRandMax;
			LU.at(i,j) = M.at(i,j);
		}
	}
}

void printm(Matrix& matrix){
	cout<<  matrix.size <<"\n";
	for(long i = 0; i <= matrix.size-1; i++){
		for(long j = 0; j <= matrix.size-1; j++)
			cout<< matrix.at(i,j) <<"\t";
		cout<< endl;
	}
}

template <class T>
void printv(vector<T>& x){
	for(T& vx : x)
		cout << vx <<'\t';
}

}

#endif // MATRIX_H