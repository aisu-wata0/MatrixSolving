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
#define L1_LINE_DN (CACHE_LINE_SIZE/sizeof(double)) // how many doubles in a line

#define CACHE_L1_SIZE 32*1024/2 // likwid-topology: Size:	 32 kB
// divided by 2 because we wont be able to fill L1 completely without throwing
// useful values out
#define L1_DN (CACHE_L1_SIZE/sizeof(double)) // how many doubles in L1 cache
// 2048
// size of the block that should fit into L1
#define BL1 (((size_t)sqrt(L1_DN)) - mod(((size_t)sqrt(L1_DN)), L1_LINE_DN))
// 40 % 8 == 0, (40*40 < 2048)

#define REG_SZ 32 // how many bytes in a register
#define dn (REG_SZ/sizeof(double)) // how many doubles is a register
typedef double v4df __attribute__ ((vector_size (REG_SZ))); // vector of four doubles

typedef union vdouble
{
  v4df v;
  double d[dn];
} vdouble;

typedef union vdoublep
{
  vdouble* v;
  double* d;
} vdoublep;


#define PAD(X) (div_down((X),L1_LINE_DN)*(L1_LINE_DN*(L1_LINE_DN-1))/2)
// Optm: test switching, the below doesnt work probably
//#define PAD(X) ((long)floor((X)/(double)L1_LINE_DN)*(L1_LINE_DN*(L1_LINE_DN-1))/2)

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
	vdoublep arr;
	long size;
	long m_size;
	long m_sizev;
	
	void mem_alloc(long size){
		arr.d = (double*)malloc((size*size)*sizeof(double));
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	Matrix(long size)
	:	size(size){
		if(PADDING){
			m_size = size + mod(L1_LINE_DN - size, L1_LINE_DN);
			if(mod(m_size/L1_LINE_DN, 2) == 0){
				m_size = m_size + L1_LINE_DN; // make sure m_size is odd multiple of cache line
			}
		} else {
			m_size = size;
		}
		m_sizev = m_size/dn;
		mem_alloc(m_size);
	}
	
	long m_posv(long i, long j) const {
		return i*m_sizev + j;
	}
	vdouble& atv(long i, long j) {
		return arr.v[m_posv(i,j)];
	}
	const vdouble& atv(long i, long j) const {
		return arr.v[m_posv(i,j)];
	}
	long m_pos(long i, long j) const {
		return i*m_size + j;
	}
	double& at(long i, long j){
		return arr.d[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr.d[m_pos(i,j)];
	}
};

/**
 * @brief Stores values of matrix in a vector, Column Major Order
 */
class MatrixColMajor : public Matrix
{
public:
	using Matrix::Matrix;
	
	long m_posv(long i, long j) const {
		return j*m_sizev + i;
	}
	vdouble& atv(long i, long j) {
		return arr.v[m_posv(i,j)];
	}
	const vdouble& atv(long i, long j) const {
		return arr.v[m_posv(i,j)];
	}
	long m_pos(long i, long j) const {
		return j*m_size + i;
	}
	double& at(long i, long j){
		return arr.d[m_pos(i,j)];
	}
	const double& at(long i, long j) const {
		return arr.d[m_pos(i,j)];
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
void print(Mat& M){
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
void printm(Mat& M){
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


}
#endif
