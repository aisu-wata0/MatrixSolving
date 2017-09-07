#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctgmath>

namespace std {

#define EPSILON 1e-7

bool close_zero(double x){
	return fabs(x) < EPSILON;
}

/**
 * @class Matrix
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