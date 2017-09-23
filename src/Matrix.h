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

/**
 * @brief stores a double
 * f: double
 * exp: exponent
 * mant: mantissa
 */
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
 * Optimization: The next in memory value can be decided on the constructor
 */
class Matrix
{
public:
	long size;
	Matrix_type type;
	vector<double> matrix;

	/**
	 * @param size of matrix, total number of lines
	 * @param type optional: if the next in memory value is on the following column or line
	 * default: column
	 */
	Matrix(long size, Matrix_type type = BY_COL) : size(size), type(type), matrix(size*size){
	}

	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 */
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
	
	/**
	 * @brief Increments b into the matrix
	 * @param b
	 * @param sign -1 with you want to add -b
	 */
	void add(Matrix& b, double sign = 1){
		for(long i=0; i <= size-1; i++){
			for(long j=0; j <= size-1; j++){
				this->at(i,j) += sign*b.at(i,j);
			}
		}
	}
	
	/**
	 * @brief prints matrix coeficients
	 */
	void print(){
		for(long i = 0; i < size; i++){
			for(long j = 0; j < size; j++){
				cout << this->at(i, j) <<'\t';
			}
			cout << endl;
		}
	}

	/**
	 * @brief prints matrix coeficients together with B
	 * @param B
	 */
	void print(vector<double>& B){
		for (long i = 0; i < size; i++) {
			for (long j = 0; j < size; j++) {
				cout << this->at(i, j) <<"x"<< j <<'\t';
			}
			cout <<"=\t"<< B[i] << endl;
		}
	}
};

/**
 * @brief sets I to identity
 * @param I
 */
void identity(Matrix& I){
	for(long i = 0; i < I.size; i++){
		for(long j = 0; j < I.size; j++){
			I.at(i, j) = 0;
		}
		I.at(i, i) = 1;
	}
}

/**
 * @brief prints random matrix
 * @param n size
 */
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

/**
 * @brief Assigns the same random matrix to M and LU,
 * they need to be the same size and have been initialized with some size
 * @param M
 * @param LU
 */
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

/**
 * @brief prints matrix
 * @param matrix
 */
void printm(Matrix& matrix){
	cout<<  matrix.size <<"\n";
	for(long i = 0; i <= matrix.size-1; i++){
		for(long j = 0; j <= matrix.size-1; j++)
			cout<< matrix.at(i,j) <<"\t";
		cout<< endl;
	}
}

/**
 * @brief prints vector x
 * @param x
 */
template <class T>
void printv(vector<T>& x){
	for(T& vx : x)
		cout << vx <<'\t';
}

}

#endif // MATRIX_H