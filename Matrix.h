#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

namespace std {

#define EPSILON 1e-7
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

/**
 * @brief Stores values from specified type of matrix, can be acessed using conventional indexing .at(i,j)
 */
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
	
	void add(Matrix& b){
		for(long i=0; i <= size-1; i++){
			for(long j=0; j <= size-1; j++){
				this->at(i,j) += b.at(i,j);
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

}

#endif // MATRIX_H