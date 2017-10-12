#ifndef MATRIX_H
#define MATRIX_H
/**
@file Matrix.h
*/
namespace std {

#define mod(X,Y) (((X) % (Y)) < 0 ? ((X) % (Y)) + (Y) : ((X) % (Y)))
#define CACHE_LINE_SIZE 8
#define PAD(X) ((X)+1)/CACHE_LINE_SIZE*(CACHE_LINE_SIZE*(CACHE_LINE_SIZE-1))/2
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
	/**
	 * @brief swaps rows starting from col 'start'
	 */
	void swap_rows_from(long row0, long row1, long start){
		if(row0 == row1)
			return;
		for(long j = start; j < size; j++){
			// for each collumn
			swap(this->at(row0, j), this->at(row1, j));
		}
	}

	void swap_rows(long row0, long row1){
		swap_rows_from(row0, row1, 0);
	}

	void resize(long new_size){
		 size = new_size;
		 matrix.resize(size*size);
	}
	/**
	 * @brief Increments b into the matrix
	 * @param b
	 * @param sign -1 with you want to add -b
	 */
	void add(Matrix& b, double sign = 1){
		for(long i=0; i < size; i++){
			for(long j=0; j < size; j++){
				this->at(i,j) += sign*b.at(i,j);
			}
		}
	}
	/**
	 * @brief copy matrix M to yourself
	 */
	void set(Matrix& M){
		for(long i=0; i < size; i++){
			for(long j=0; j < size; j++){
				this->at(i,j) = M.at(i,j);
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

/**
 * @brief sets I to identity
 * @param I
 */
void identity(MatrixColMajor& I){
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
void randomMatrix(Matrix& M){
	long i, j;
	double invRandMax = 1.0/(double)RAND_MAX;

	for(i = 0; i < M.size; i++){
		for(j = 0; j < M.size; j++){
			M.at(i,j) = (double)rand() * invRandMax;
		}
	}
}
/**
 * @brief prints matrix
 * @param matrix
 */
void printm(Matrix& matrix){
	cout<<  matrix.size <<"\n";
	for(long i = 0; i < matrix.size; i++){
		for(long j = 0; j < matrix.size; j++)
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


/**
 * @brief Stores values of a triangular matrix in a array, each line goes only until elements of the diagonal
 * ex: having ​​CACHE_LINE_SIZE = ​3​
 *   0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2
 * [00, -​, -​|10,11​, -​|20,21,22|30,31,32,33​, -, -​|40,41,42,43,44, -|​50,...​]
 * pad=​2​;​   pad=​1​;   pad=​0;   ​pad=​2​;​            ​pad=​1​;​            pad=0;​
 * +2      ​ +1​       +0       +2​ ​               +1​​                +0
 * // having i​​ = 4;
 * // 3 = ​(5/3 -2/3)*  3*2/2​; // ok​
 * padding_before = (i+1)/​​​CACHE_LINE_SIZE *(​​​​​CACHE_LINE_SIZE *(​CACHE_LINE_SIZE-1)/2);
 * a1 = ​​​CACHE_LINE_SIZE​ - mod(i,​​CACHE_LINE_SIZE​)​;
 * sum = (​CACHE_LINE_SIZE​ ​- a1)*(a1+CACHE_LINE_SIZE​-1​)/2;
 * sum = mod(i,​​CACHE_LINE_SIZ​​E​)*(​​2*CACHE_LINE_SIZE​ ​​-1 - mod(i,​​CACHE_LINE_SIZE​)​​)/2​;
 * // 2 = mod(4,3)*(2*3 -1 - mod(4,3))/2; // ok
 * ​pad_btotal = ​padding_before​+​​sum 
 * [i*(i+1)/2​ + ​​pad_btota​l​​ + j​​]​
 */
class MatrixTriSeq
{
public:
	double* matrix;
	long size;
	/**
	 * @param size of matrix, total number of lines
	 */
	MatrixTriSeq(long size){
		malloc(size);
	}

	MatrixTriSeq(){
	}
	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 */
	double& at(long i, long j) {
		long sum = mod(i, CACHE_LINE_SIZE)*(2*CACHE_LINE_SIZE -1 - mod(i, CACHE_LINE_SIZE))/2;
		long pad_btotal = PAD(i) + sum;
		return matrix[i*(i+1)/2 + pad_btotal + j];
		//i*(i+1)/2​ + (i+1)*​​​(x-1)/2 ​+ mod(i,​x​)*(​​2*x ​​-1 - mod(i,​x​)​​)/2​ + j​;
	}

	void swap_rows(long row0, long row1){
		swap_rows_from(row0, row1, 0);
	}

	void malloc(long new_size){
		size = new_size;
		matrix = malloc((size*size + PAD(size))*sizeof(double));
	}
	/**
	 * @brief Increments b into the matrix
	 * @param b
	 * @param sign -1 with you want to add -b
	 */
	void add(Matrix& b, double sign = 1){
		for(long i=0; i < size; i++){
			for(long j=0; j < i+1; j++){
				this->at(i,j) += sign*b.at(i,j);
			}
		}
	}
	/**
	 * @brief copy matrix M to yourself
	 */
	void set(Matrix& M){
		for(long i=0; i < size; i++){
			for(long j=0; j < i+1; j++){
				this->at(i,j) = M.at(i,j);
			}
		}
	}

	/**
	 * @brief prints matrix
	 */
	void print(){
		for(long i = 0; i < size; i++){
			for(long j = 0; j < i+1; j++){
				cout << this->at(i, j) <<'\t';
			}
			for(long j = i+1; j < size; j++){
				cout << (double)0.0 <<'\t';
			}
			cout << endl;
		}
	}
};

/**
 * @brief swaps rows starting from col 'start'
 */
void swap_rows(MatrixTriSeq& Low, MatrixTriSeq& Upp, long row0, long row1){
	if(row0 == row1)
		return;
	
	if(row0 > row1){swap(row0, row1);}
	/* 6 1 2
	 * 3 4 5
	 * 0 7 8 */
	for(long j = 0; j < row0+1; j++){
		// for each collumn
		swap(Low.(row1, j), Low.at(row0, j));
	}
	/* 6 7 2
	 * 3 4 5
	 * 0 1 8 */
	for(long j = row0+1; j < row1+1; j++){
		// for each collumn
		swap(Upp.at(row0, j), Low.at(row1, j));
	}
	/* 6 7 8
	 * 3 4 5
	 * 0 1 2 */
	for(long j = row1; j < size; j++){
		// for each collumn
		swap(Upp.at(row0, j), Upp.at(row1, j));
	}
}


}
#endif
