#ifndef MATRIX_H
#define MATRIX_H
/**
@file Matrix.h
*/
namespace std {

long div_down(long n, long d) {
    return n / d - (((n > 0) ^ (d > 0)) && (n % d));
}

#define mod(X,Y) ((((X) % (Y)) + (Y)) % Y)
#define CACHE_LINE_SIZE 16
#define PAD(X) (div_down((X),CACHE_LINE_SIZE)*(CACHE_LINE_SIZE*(CACHE_LINE_SIZE-1))/2)
// Optm: test switching, the below doesnt work probably
//#define PAD(X) ((long)floor((X)/(double)CACHE_LINE_SIZE)*(CACHE_LINE_SIZE*(CACHE_LINE_SIZE-1))/2)

#define PADDING true

/**
 * @brief Stores values of matrix in a vector, Row Major Order
 */
class Matrix
{
public:
	long size;
	vector<double> matrix;
	/**
	 * @param size of matrix, total number of lines
	 */
	Matrix(long size) : size(size), matrix(size*size){
	}
	
	Matrix(){
	}
	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 */
	double& at(long i, long j) {
		return matrix.at(i*size + j);
	}
	
	void resize(long new_size){
		 size = new_size;
		 matrix.resize(size*size);
	}
};

/**
 * @brief Stores values of matrix in a vector, Column Major Order
 */
class MatrixColMajor : public Matrix
{
public:
	using Matrix::Matrix;
	
	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 */
	double& at(long i, long j) {
		return matrix.at(j*size + i);
	}
};

template<class Mat>
void swap_rows(Mat& M, long row0, long row1){
	if(row0 == row1)
		return;
	for(long j = 0; j < M.size; j++){
		swap(M.at(row0, j), M.at(row1, j));
	}
}
/**
 * @brief Increments b into the matrix
 * @param b
 * @param sign -1 with you want to add -b
 */
template<class Mat>
void add(Mat& M, Matrix& b, double sign = 1){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < M.size; j++){
			M.at(i,j) += sign*b.at(i,j);
		}
	}
}
/**
 * @brief copy matrix A to yourself
 */
template<class Mat>
void set(Mat& M, Matrix& A){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < M.size; j++){
			M.at(i,j) = A.at(i,j);
		}
	}
}

template<class Mat>
void print(Mat& M){
	for(long i = 0; i < M.size; i++){
		for(long j = 0; j < M.size; j++){
			cout << M.at(i, j) <<'\t';
		}
		cout << endl;
	}
}
/**
 * @brief sets I to identity
 */
template<class Mat>
void identity(Mat& I){
	for(long i = 0; i < I.size; i++){
		for(long j = 0; j < I.size; j++){
			I.at(i, j) = 0;
		}
		I.at(i, i) = 1;
	}
}
/**
 * @brief Assigns random matrix to M
 * @param M needs to have been allocated
 */
template<class Mat>
void randomMatrix(Mat& M){
	long i, j;
	double invRandMax = 1.0/(double)RAND_MAX;
	
	for(i = 0; i < M.size; i++){
		for(j = 0; j < M.size; j++){
			M.at(i,j) = (double)rand() * invRandMax;
		}
	}
}
/**
 * @brief prints matrix with size in the first line
 */
template<class Mat>
void printm(Mat& matrix){
	cout<<  matrix.size <<"\n";
	print(matrix);
}
/**
 * @brief prints vector x
 */
template <class T>
void printv(vector<T>& x){
	for(T& vx : x)
		cout << vx <<'\t';
}


}
#endif
