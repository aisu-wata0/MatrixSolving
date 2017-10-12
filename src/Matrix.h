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
#define CACHE_LINE_SIZE 4
#define PAD(X) ((long)floor(((X)+1)/(double)CACHE_LINE_SIZE)*(CACHE_LINE_SIZE*(CACHE_LINE_SIZE-1))/2)
// Optm: test
//#define PAD(X) (div_down(((X)+1),CACHE_LINE_SIZE)*(CACHE_LINE_SIZE*(CACHE_LINE_SIZE-1))/2)

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
 * @brief prints matrix with size in the first line
 */
template <class T>
void printm(T& matrix){
	cout<<  matrix.size <<"\n";
	for(long i = 0; i < matrix.size; i++){
		for(long j = 0; j < matrix.size; j++)
			cout<< matrix.at(i,j) <<"\t";
		cout<< endl;
	}
}
/**
 * @brief prints vector x
 */
template <class T>
void printv(vector<T>& x){
	for(T& vx : x)
		cout << vx <<'\t';
}

/**
 * @brief Stores values of a triangular matrix in a array, each line goes only until elements of the diagonal
 */
class MatrixTri
{
public:
	double* matrix;
	long size;
	long m_size;
	
	void reserve(long new_size){
		size = new_size;
		m_size = size*size + PAD(size);
		matrix = (double*)malloc(m_size*sizeof(double));
	}
};

/**
 * @class MatrixTriUpp
 * @brief Stores values of a upper triangular matrix in a array, each line goes only until elements of the diagonal
 */
class MatrixTriLow : public MatrixTri
{
public:
	long desl;
	long first_pad_seq;
	
	MatrixTriLow(){
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	MatrixTriLow(long size){
		reserve(size);
	}
	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 * ex: having ​​CACHE_LINE_SIZE = ​3​
	 *   0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2
	 * [00, -​, -​|10,11​, -​|20,21,22|30,31,32,33​, -, -​|40,41,42,43,44, -|​50,...​]
	 * pad=​2​;​   pad=​1​;   pad=​0;   ​pad=​2​;​            ​pad=​1​;​            pad=0;​
	 * +2      ​ +1​       +0       +2​ ​               +1​​                +0
	 * // having i​​ = 4;
	 * // 3 = ​floor(5/3) * 3*2/2​; // ok​
	 * padding_before = (i+1)/​​​CACHE_LINE_SIZE *(​​​​​CACHE_LINE_SIZE *(​CACHE_LINE_SIZE-1)/2);
	 * a1 = ​​​CACHE_LINE_SIZE​ - mod(i,​​CACHE_LINE_SIZE​)​;
	 * sum = (​CACHE_LINE_SIZE​ ​- a1)*(a1+CACHE_LINE_SIZE​-1​)/2;
	 * sum = mod(i,​​CACHE_LINE_SIZ​​E​)*(​​2*CACHE_LINE_SIZE​ ​​-1 - mod(i,​​CACHE_LINE_SIZE​)​​)/2​;
	 * // 2 = mod(4,3)*(2*3 -1 - mod(4,3))/2; // ok
	 * ​pad_btotal = ​padding_before​+​​sum 
	 * [i*(i+1)/2​ + ​​pad_btota​l​​ + j​​]​		*/
	double& at(long i, long j) {
		long pos = mod(i, CACHE_LINE_SIZE);
		long sum = pos*(2*CACHE_LINE_SIZE -1 - pos)/2;
		long pad_btotal = PAD(i) + sum;
		return matrix[i*(i+1)/2 + pad_btotal + j];
		//i*(i+1)/2​ + (i+1)*​​​(x-1)/2 ​+ mod(i,​x​)*(​​2*x ​​-1 - mod(i,​x​)​​)/2​ + j​;
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
 * @class MatrixTriUpp
 * @brief Stores values of a upper triangular matrix in a array, each line goes only until elements of the diagonal
 */
class MatrixTriUpp : public MatrixTri
{
public:
	long desl;
	long first_pad_seq;
	
	MatrixTriUpp(){
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	MatrixTriUpp(long size){
		reserve(size);
	}
	/**
	 * @param i
	 * @param j
	 * @return element of position i,j
	 * ex: having ​​CACHE_LINE_SIZE = ​4; size % CACHE_LINE_SIZE == 2;
	 * row = 0  1  2  3  4  5  6  7  8  9  0  1  2
	 * pad = ​2  3  0  1  2  3  0  1  2  3  0  1  2
	 *       ----     -------------------     ----
	 * first_pad_seq       PAD(i-desl+1)       sum
	 * pad_btotal = first_pad_seq + PAD(i-desl+1) + sum
	 * at(i,j) = [i*(i+1)/2​ + ​​pad_btota​l​​ + j​​ -i]​ */
	double& at(long i, long j) {
		long pos = mod(i-desl,CACHE_LINE_SIZE);
		long sum = pos*(pos+1)/2;
		long pad_btotal = first_pad_seq + PAD(i-desl-1) + sum;
		return matrix[ i*(2*size-i+1)/2 + pad_btotal + j -i];
		//i*(i+1)/2​ + (i+1)*​​​(x-1)/2 ​+ mod(i,​x​)*(​​2*x ​​-1 - mod(i,​x​)​​)/2​ + j​;
	}
	void reserve(long new_size){
		desl = mod((new_size+1),CACHE_LINE_SIZE);
		long a1 = mod(new_size,CACHE_LINE_SIZE);
		long an = CACHE_LINE_SIZE-1;
		long n = an -a1 +1;
		first_pad_seq = n*(a1 + an)/2;
		MatrixTri::reserve(new_size);
	}
	/**
	 * @brief copy matrix M to yourself
	 */
	void set(Matrix& M){
		for(long i=0; i < size; i++){
			for(long j=i; j < size; j++){
				this->at(i,j) = M.at(i,j);
			}
		}
	}
	
	void print(){
		for(long i = 0; i < size; i++){
			for(long j = 0; j < i; j++){
				cout << (double)0.0 <<'\t';
			}
			for(long j = i; j < size; j++){
				cout << this->at(i, j) <<'\t';
			}
			cout << endl;
		}
	}
};

/**
 * @brief copy matrix M to LU
 */
template<class MLower, class MUpper>
void set(MLower& Low, MUpper& Upp, Matrix& M){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < i+1; j++){
			Low.at(i,j) = M.at(i,j);
		}
		for(long j=i; j < M.size; j++){
			Upp.at(i,j) = M.at(i,j);
		}
	}
}

/**
 * @brief swaps rows starting from col 'start'
 */
 template<class MLower, class MUpper>
void swap_rows(MLower& Low, MUpper& Upp, long row0, long row1){
	if(row0 == row1)
		return;
	/* 1 0 0  0 1 2  row0 = 1, row1 = 2
	 * 3 1 0  0 4 5
	 * 6 7 1  0 0 8 */
	if(row0 > row1){swap(row0, row1);}
	/* 1 0 0  0 1 2
	 * 6 1 0  0 4 5
	 * 3 7 1  0 0 8 */
	for(long j = 0; j < row0; j++){
		// for each collumn
		swap(Low.at(row0, j), Low.at(row1, j));
	}
	/* 1 0 0  0 1 2
	 * 6 1 0  0 7 5
	 * 3 4 1  0 0 8 */
	for(long j = row0; j < row1; j++){
		// for each collumn
		swap(Upp.at(row0, j), Low.at(row1, j));
	}
	/* 1 0 0  0 1 2
	 * 6 1 0  0 7 8
	 * 3 4 1  0 0 5 */
	for(long j = row1; j < Upp.size; j++){
		// for each collumn
		swap(Upp.at(row0, j), Upp.at(row1, j));
	}
}


}
#endif
