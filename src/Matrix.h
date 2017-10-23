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

#define CACHE_LINE_SIZE 64 // likwid-topology: Cache line size:	64
//#define CACHE_LSZ CACHE_LINE_SIZE/sizeof(double)
#define CACHE_LSZ 16 // how many doubles in a line

#define PAD(X) (div_down((X),CACHE_LSZ)*(CACHE_LSZ*(CACHE_LSZ-1))/2)
// Optm: test switching, the below doesnt work probably
//#define PAD(X) ((long)floor((X)/(double)CACHE_LSZ)*(CACHE_LSZ*(CACHE_LSZ-1))/2)

#define PADDING true

// Non member access functions
template<class Mat>
double& at(Mat& M, long i, long j){
	return M.at(i,j);
}
template<class Mat>
const double& at(Mat const& M, long i, long j){
	return M.at(i,j);
}

/**
 * @brief Stores values of matrix in a vector, Row Major Order
 */
class Matrix
{
public:
	double* arr;
	long size;
	long m_size;
	
	void mem_alloc(long size){
		arr = (double*)malloc((size*size)*sizeof(double));
	}
	
	void set_size(long new_size){
		size = new_size;
		if(PADDING){
			m_size = size + mod(CACHE_LSZ - size, CACHE_LSZ);
			if(mod(m_size/CACHE_LSZ, 2) == 0){
				m_size = m_size + CACHE_LSZ; // make sure m_size is odd multiple of cache line
			}
		} else {
			m_size = size;
		}
		mem_alloc(m_size);
	}
	
	Matrix(){
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	Matrix(long size){
		set_size(size);
	}
	
	~Matrix(){
		free(arr);
	}
	
	inline long m_pos(long i, long j) const {
		return i*m_size + j;
	}
	double& at(long i, long j) {
		return arr[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr[m_pos(i,j)];
	}
};

/**
 * @brief Stores values of matrix in a vector, Column Major Order
 */
class MatrixColMajor : public Matrix
{
public:
	using Matrix::Matrix;
	
	inline long m_pos(long i, long j) const {
		return j*m_size + i;
	}
	double& at(long i, long j) {
		return arr[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr[m_pos(i,j)];
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
 * @brief M += B
 * @param sign -1 with you want to add -b
 */
template<class Mat>
void add(Mat& M, Mat& B, double sign = 1){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < M.size; j++){
			M.at(i,j) += sign*B.at(i,j);
		}
	}
}
/**
 * @brief copy matrix A to yourself
 */
template<class Mat>
void set(Mat& M, const Matrix& A){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < M.size; j++){
			M.at(i,j) = A.at(i,j);
		}
	}
}
/**
 * @brief sets all matrix to parameter
 */
template<class Mat>
void set(Mat& M, double x){
	for(long i=0; i < M.size; i++){
		for(long j=0; j < M.size; j++){
			M.at(i,j) = x;
		}
	}
}

template<class Mat>
void print(const Mat& M){
	for(long i = 0; i < M.size; i++){
		for(long j = 0; j < M.size; j++){
			cout << M.at(i, j) <<' ';
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
void printm(const Mat& M){
	cout<<  M.size <<"\n";
	print(M);
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
 * @class MatrixTriUpp
 * @brief Lower triangular implementation of MatrixTri
 */
class MatrixTriLow : public Matrix
{
public:
	/**
	 * ex: having CACHE_LSZ = 3;
	 * [00, -, -|10,11, -|20,21,22|30,31,32,33, -, -|40,41,42,43,44, -|50,...]
	 * pad=2;    pad=1;   pad=0;   pad=2;            pad=1;           pad=0;
	 * pad_total(i = 4) == +2 +1 +0 +2 = 5
	 * https://docs.google.com/spreadsheets/d/1y2ffNz4jD6LwmL2JnWPCmYlUoSuY7Q04Ey6REYM7DXg/edit?usp=sharing */
	inline long pad_total(long i) const {
		if(not PADDING) return 0;
		long pos = mod(i, CACHE_LSZ);
		long sum = pos*(2*CACHE_LSZ -1 - pos)/2;
		return (PAD(i) + sum);
	}
	
	inline long m_pos(long i, long j) const {
		return ( i*(i+1)/2 + j + pad_total(i) );
	}
	double& at(long i, long j) {
		return arr[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr[m_pos(i,j)];
	}
	
	void set_size(long new_size){
		size = new_size;
		m_size = m_pos(size-1, size-1) +1;
		mem_alloc(m_size);
	}
	
	MatrixTriLow(){
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	MatrixTriLow(long size){
		set_size(size);
	}
	
	/**
	 * @brief copy matrix M to yourself
	 */
	void set(const Matrix& M){
		for(long i=0; i < size; i++){
			for(long j=0; j < i+1; j++){
				this->at(i,j) = M.at(i,j);
			}
		}
	}
	
	void test(){
		for(long i = 0; i < size; i++){
			cout<<"row "<< i << " padding = "<< pad_total(i) <<endl;
			for(long j = 0; j < i+1; j++){
				this->at(i, j) = 1;
			}
		}
		for(long i = 0; i < size; i++){
			for(long j = 0; j < i+1; j++){
				if(this->at(i, j) == 0){
					cerr<<"TWO POSITIONS ACESSING SAME MEMORY: "<<"at("<< i <<","<< j <<")"<<endl;
				}
				this->at(i, j) = 0;
			}
		}
	}
};

void print(const MatrixTriLow& L){
	long size = L.size;
	for(long i = 0; i < size; i++){
		for(long j = 0; j < i+1; j++){
			cout << L.at(i, j) <<'\t';
		}
		for(long j = i+1; j < size; j++){
			cout << (double)0.0 <<'\t';
		}
		cout << endl;
	}
}

/**
 * @class MatrixTriUpp
 * @brief Upper triangular implementation of Matrix
 */
class MatrixTriUpp : public Matrix
{
public:
	long desl;
	long first_pad_seq;
	/**
	 * ex: having CACHE_LSZ = 4; (size % CACHE_LSZ) == 2; (meaning first row will have 2 pads
	 * row = 0  1  2  3  4  5  6  7  8  9  0  1  2
	 * pad = 2  3  0  1  2  3  0  1  2  3  0  1  2
	 *       ----     -------------------     ----
	 * first_pad_seq       PAD(i-desl+1)       sum
	 * pad_btotal = first_pad_seq + PAD(i-desl+1) + sum 
	 * https://docs.google.com/spreadsheets/d/1y2ffNz4jD6LwmL2JnWPCmYlUoSuY7Q04Ey6REYM7DXg/edit?usp=sharing */
	inline long pad_total(long i) const {
		if(not PADDING) return 0;
		long pos = mod(i-desl,CACHE_LSZ);
		long sum = pos*(pos+1)/2;
		return (first_pad_seq + PAD(i-desl) + sum);
	}
	
	inline long m_pos(long i, long j) const {
		return ( i*(2*size-i+1)/2 + j -i + pad_total(i) );
	}
	double& at(long i, long j) {
		return arr[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr[m_pos(i,j)];
	}
	
	void set_size(long new_size){
		size = new_size;
		desl = mod((size+1),CACHE_LSZ);
		if(mod(size,CACHE_LSZ) != CACHE_LSZ-1){
			long a1 = CACHE_LSZ - mod(size,CACHE_LSZ);
			long an = CACHE_LSZ-1;
			long n = an -a1 +1;
			first_pad_seq = n*(a1 + an)/2;
		} else {
			first_pad_seq = 0;
		}
		m_size = m_pos(size-1, size-1) +1;
		mem_alloc(m_size);
	}
	
	MatrixTriUpp(){
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	MatrixTriUpp(long size){
		set_size(size);
	}
	/**
	 * @brief copy matrix M to yourself
	 */
	void set(const Matrix& M){
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
	
	void test(){
		for(long i = 0; i < size; i++){
			cout<<"row "<< i << " padding = "<< pad_total(i) <<endl;
			for(long j = i; j < size; j++){
				this->at(i, j) = 1;
			}
		}
		for(long i = 0; i < size; i++){
			for(long j = i; j < size; j++){
				if(this->at(i, j) == 0){
					cerr<<"TWO POSITIONS ACESSING SAME MEMORY"<<endl;
				}
				this->at(i, j) = 0;
			}
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

template<class MLower, class MUpper>
void printLU(MLower& Low, MUpper& Upp){
	long size = Low.size;
	for(long i = 0; i < size; i++){
			for(long j = 0; j < i; j++){
					cout << Low.at(i,j) <<'\t';
			}
			for(long j = i; j < size; j++){
					cout << Upp.at(i,j) <<'\t';
			}
			cout << endl;
   }
}

}
#endif
