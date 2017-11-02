#ifndef BYTES_H
#define BYTES_H

#define mod(X,Y) ((((X) % (Y)) + (Y)) % Y)
#define Lower_Multiple(SZ,D) ((SZ) - ((SZ) % D))

#define CACHE_LINE_SIZE (64) // likwid-topology: Cache line size:	64
#define L1_LINE_DN (CACHE_LINE_SIZE/sizeof(double)) // how many doubles in a line
// likwid-topology: Size
#define L1KiB 32
#define CACHE_L1_SIZE (L1KiB*1024/2)
#define CACHE_L2_SIZE (256*1024/2)
#define CACHE_L3_SIZE (3*1024*1024/2)
// divided by 2 because we wont be able to fill L1 completely without throwing
// useful values out

// aproximate minimum number of lines L1 cache has
// (for this capacity)(min lines mean max associativ)
#define L1LINE_N ((L1KiB*1024/8)/64)

// 2048
#define L1_DN (CACHE_L1_SIZE/sizeof(double)) // how many doubles in L1 cache
#define L2_DN (CACHE_L2_SIZE/sizeof(double)) // how many doubles in L1 cache
#define L3_DN (CACHE_L3_SIZE/sizeof(double)) // how many doubles in L1 cache

// size of the block to fit one matrix in L1
#define MAX_BL1 ((long)sqrt(L1_DN))
// align to cache line
#define BL1 (Lower_Multiple(MAX_BL1, L1_LINE_DN))
// 40 % 8 == 0, (40*40 < 2048)

// to fit 3 matrixes
#define MAX_B3L1 ((long)sqrt(L1_DN/3))
#define B3L1 (Lower_Multiple(MAX_B3L1, L1_LINE_DN))

#define MAX_B3L2 ((long)sqrt(L2_DN/3))
#define B3L2 (Lower_Multiple(MAX_B3L2, B3L1))

#define MAX_B3L3 ((long)sqrt(L3_DN/3))
#define B3L3 (Lower_Multiple(MAX_B3L3, B3L2))

#define REG_SZ (32) // how many bytes in a register

#define regDN (REG_SZ/sizeof(double)) // how many doubles is a register
#define regFN (REG_SZ/sizeof(float)) // how many floats is a register

#endif


#ifndef VARRAY_HPP
#define VARRAY_HPP

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>


#define unroll(v,n) for(size_t v = 0; v < n; v++)
#define unroll2(vi,vj,n) unroll(vi,n)unroll(vj,n)
#define vec(v) for(size_t v = 0; v < regDN; v++)

namespace gm
{
using namespace std;

template <typename num>
bool isPowerOfTwo (num n) {
	return (n > 0) && ((n == 0) || ((n & (n - 1)) == 0));
}

template<typename elem>
inline size_t regN() { return (REG_SZ/sizeof(elem)); }

/**
 * @class vec
 * @brief Vectorized type
 * the number of elem in the vec is regN<elem>()*/
template<typename elem>
struct vec
{
	elem __attribute__ ((vector_size (REG_SZ)))  v; // vectorization of elems

	inline const elem& operator[] (size_t i) const {
		assert(i < regN<elem>() && "Vectorized elem out of register access");
		return v[i];
	}
	inline elem& operator[] (size_t i) {
		assert(i < regN<elem>() && "Vectorized elem out of register access");
		return v[i];
	}
};

template<typename elem>
union vecp
{
	vec<elem>*  v;
	elem* p;
};

/**
 * @brief Allocates size*sizeof(elem) and aligns it into a boundary-bit boundary
 * @param ppElement pointer to array (pointer)
 * @param size number of elems in your resulting pointer
 * @param boundary : Power of two, else the behavior is undefined
 * @return The pointer you should free, don't lose it
 */
template<typename elem>
void* al_allloc(elem** ppElement, size_t size, size_t boundary){
	*ppElement = (elem*)malloc((size)*sizeof(elem) + boundary-1);
	if(*ppElement == NULL){
		cerr <<"failed to malloc "<< (size)*sizeof(elem)/1024 <<" KiB"<< endl;
		exit(0);
	}
	void* pMem = *ppElement;
	*ppElement = (elem*)(((uintptr_t)pMem + (boundary-1)) & (uintptr_t)( ~ (boundary-1)));
	return pMem;
}

template<class elem>
class varray
{
public:
	vecp<elem> arr;
	size_t mSize; //  n of types
	size_t mSizeVec; // n of vec<elem>s

	size_t mSizeMem; // n of elem in memory
	size_t mSizeVecMem; // n of vec<elem>s in memory

	size_t mEndVec; // elem index where vectorization ends
	void* mpMem; // pointer to be freed

	/** @brief elem number in a register */
	inline size_t regEN(){ return regN<elem>(); }

	void memAlloc(size_t mSizeMem){
		if(mpMem != NULL) { free(mpMem); }
		mpMem = al_allloc(&arr.p, mSizeMem, CACHE_LINE_SIZE);
	}

	void alloc(size_t size){
		mSize = size;
		mSizeVec = mSize/regEN();

		mSizeMem = mSize;
		mSizeMem = mSize;
		// add to make it multiple of L1_LINE_DN
		mSizeMem += mod(-mSize,(ptrdiff_t)L1_LINE_DN);
		//mSizeMem += mod(-mSize, (ptrdiff_t)L1_LINE_DN); TODO test

		//if(((mSizeMem/L1_LINE_DN) % 2) == 0)
		if(mSizeMem > L1LINE_N-1 && isPowerOfTwo(mSizeMem))
			mSizeMem = mSizeMem + L1_LINE_DN; // make sure mSizeMem is odd multiple

		mSizeVecMem = mSizeMem/regEN();

		mEndVec = Lower_Multiple(mSize, regEN());

		memAlloc(mSizeMem);
	}

	varray() {
		mpMem = NULL;
	}
	varray(size_t size) {
		mpMem = NULL;
		alloc(size);
	}
	~varray(){
		if(mpMem != NULL) { free(mpMem); }
	}

	/** @brief size of the varray */
	size_t size(){ return mSize; }
	/** @brief size of the vectorized varray */
	size_t sizeVec(){ return mSizeVec; }

	/** @brief vectorization end index */
	size_t vecEnd(){ return mEndVec; }
	/**
	 * @brief Vectorized access
	 * @return return vec in the position i,
	 * it will contain elems (i*regEN()) to (i*regEN() + regEN() -1)
	 */
	inline vec<elem>& atv(size_t i) {
		assert(i < mSizeVec && "varray vec access out of bounds");
		return arr.v[i];
	}
	inline const vec<elem> & atv(size_t i) const {
		assert(i < mSizeVec && "varray vec access out of bounds");
		return arr.v[i];
	}

	inline elem& at(size_t i){
		assert(i < mSize && "varray access out of bounds");
		return arr.p[i];
	}
	inline const elem& at(size_t i) const {
		assert(i < mSize && "varray access out of bounds");
		return arr.p[i];
	}

	elem* begin() { return &arr.p[0]; }
	elem* end() { return &arr.p[size()]; }

	vec<elem>* beginVec() { return &arr.v[0]; }
	vec<elem>* endVec() { return &arr.v[sizeVec()]; }

	const elem& operator[] (size_t i) const { return at(i); }
	elem& operator[] (size_t i) {return at(i); }
};

/**
 * @brief prints vector
 */
template<class Cont>
void printv(Cont& C){
	for(size_t i = 0; i < C.size(); i++)
		cout << C.at(i) <<" ";
}
/**
 * @brief prints vector
 */
template<template<class> class Cont, class Elem>
void printv(Cont<Elem>& C){
	for(Elem& x : C)
		cout << x <<" ";
}


}
#endif


#ifndef MATRIX_HPP
#define MATRIX_HPP


namespace gm
{
using namespace std;

#define PAD(X) (div_down((X),L1_LINE_DN)*(L1_LINE_DN*(L1_LINE_DN-1))/2)
// Optm: test switching, the below doesnt work probably
//#define PAD(X) ((size_t)floor((X)/(double)L1_LINE_DN)*(L1_LINE_DN*(L1_LINE_DN-1))/2)

#define PADDING true

// Non member access functions

template<template<class> class Cont, class Elem>
Elem& at(Cont<Elem>& C, size_t i, size_t j){
	return C.at(i,j);
}
template<template<class> class Cont, class Elem>
const Elem& at(Cont<Elem> const& M, size_t i, size_t j){
	return M.at(i,j);
}

/**
 * @brief Stores values of matrix in a vector, Row Major Order
 */
template<class Elem>
class Matrix
{
public:
	varray<Elem> varr;
	size_t mSize;
	size_t mSizeVec; // n of vec<element>s

	size_t mSizeMem;
	size_t mSizeVecMem;
	size_t mPad;

	size_t mEndVec;

	/** @return Number of elements in register*/
	size_t regEN() { return varr.regEN(); }

	void memAlloc(size_t mSize){
		varr.alloc(mSize*mSize);
	}
	void alloc(size_t size){
		mSize = size;
		mSizeVec = mSize/regEN();
		mSizeMem = mSize;
		if(PADDING){
			mSizeMem += mod(-mSize,(ptrdiff_t)L1_LINE_DN);
			// make sure m_size is odd multiple of cache line
			//if(((mSizeMem/L1_LINE_DN) % 2) == 0)
			if(mSizeMem > L1LINE_N-1 && isPowerOfTwo(mSizeMem))
				mSizeMem = mSizeMem + L1_LINE_DN;
		}
		mPad = mSizeMem - mSize;
		mSizeVecMem = mSizeMem/regEN();
		mEndVec = Lower_Multiple(mSize, regEN());
		memAlloc(mSizeMem);
	}
	/**
	 * @param size of matrix, total number of lines
	 */
	Matrix(size_t size){
		alloc(size);
	}

	Matrix(){}

	/** @brief size of the varray */
	size_t size() const { return mSize; }
	/** @brief size of the vectorized varray */
	size_t sizeVec() const { return mSizeVec; }
	/** @brief size of the vectorized varray */
	size_t sizeMem() const { return mSizeMem; }
	/** @brief vectorization end index */
	size_t vecEnd(){ return mEndVec; }
	/** @brief vectorization end index */
	size_t pad(){ return mPad; }

	size_t indVecMem(size_t i, size_t j) const {
		assert(i < mSizeMem && j < mSizeVecMem);
		return i*mSizeVecMem + j;
	}
	vec<Elem>& atv(size_t i, size_t j) {
		return varr.atv(indVecMem(i,j));
	}
	const vec<Elem>& atv(size_t i, size_t j) const {
		return varr.atv(indVecMem(i,j));
	}
	size_t indMem(size_t i, size_t j) const {
		assert(i < mSizeMem && j < mSizeMem);
		return i*mSizeMem + j;
	}
	Elem& at(size_t i, size_t j){
		return varr.at(indMem(i,j));
	}
	const Elem& at(size_t i, size_t j) const {
		return varr.at(indMem(i,j));
	}
};

/**
 * @brief Stores values of matrix in a vector, Column Major Order
 */
template<class Elem>
class MatrixColMajor : public Matrix<Elem>
{
public:
	using Matrix<Elem>::Matrix;
	using Matrix<Elem>::varr;
	using Matrix<Elem>::mSizeMem;
	using Matrix<Elem>::mSizeVecMem;

	size_t indVecMem(size_t i, size_t j) const {
		assert(i < mSizeVecMem && j < mSizeMem);
		return j*mSizeVecMem + i;
	}
	vec<Elem>& atv(size_t i, size_t j) {
		return varr.atv(indVecMem(i,j));
	}
	const vec<Elem>& atv(size_t i, size_t j) const {
		return varr.atv(indVecMem(i,j));
	}
	size_t indMem(size_t i, size_t j) const {
		assert(i < mSizeMem && j < mSizeMem);
		return j*mSizeMem + i;
	}
	Elem& at(size_t i, size_t j){
		return varr.at(indMem(i,j));
	}
	const Elem& at(size_t i, size_t j) const {
		return varr.at(indMem(i,j));
	}
	/**/
};

template<class Mat>
void swap_rows(Mat& M, size_t row0, size_t row1){
	if(row0 == row1)
		return;
	for(size_t j = 0; j < M.size(); j++){
		swap(M.at(row0, j), M.at(row1, j));
	}
}
/**
 * @brief M += B
 * @param sign -1 with you want to add -b
 */
template<class Mat>
void add(Mat& M, Mat& B, double sign = 1){
	for(size_t i=0; i < M.size(); i++){
		for(size_t j=0; j < M.size(); j++){
			M.at(i,j) += sign*B.at(i,j);
		}
	}
}
/**
 * @brief copy matrix A to yourself
 */
template<class Mat, class elem>
void set(Mat& M, const Matrix<elem>& A){
	for(size_t i=0; i < M.size(); i++){
		for(size_t j=0; j < M.size(); j++){
			M.at(i,j) = A.at(i,j);
		}
	}
}
/**
 * @brief sets all matrix to parameter
 */
template<class Mat>
void set(Mat& M, double x){
	for(size_t i=0; i < M.size(); i++){
		for(size_t j=0; j < M.size(); j++){
			M.at(i,j) = x;
		}
	}
}

template<class Mat>
void print(Mat& M){
	for(size_t i = 0; i < M.size(); i++){
		for(size_t j = 0; j < M.size(); j++){
			cout << M.at(i, j) <<" ";
		}
		cout << endl;
	}
}
/**
 * @brief sets I to identity
 */
template<class Mat>
void identity(Mat& I){
	for(size_t i = 0; i < I.size(); i++){
		for(size_t j = 0; j < I.size(); j++){
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
	size_t i, j;
	double invRandMax = 1.0/(double)RAND_MAX;

	for(i = 0; i < M.size(); i++){
		for(j = 0; j < M.size(); j++){
			M.at(i,j) = (double)rand() * invRandMax;
		}
	}
}
/**
 * @brief prints matrix with size in the first line
 */
template<class Mat>
void printm(Mat& M){
	cout<<  M.size() <<"\n";
	print(M);
}


}
#endif


using namespace std;
using namespace gm;



//#define max(x,y) ((x) > (y) ? x : y)

enum class Direction {
	Forwards,
	Backwards,
};

enum class Diagonal {
	Unit,
	Value,
};

enum class Permute {
	True,
	False,
};


#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax, kmax;
	size_t bstep[5];
	//const size_t unr = 2;
	//double acc[unr*unr];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;/**/
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);

	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size);
		jmax = min(bj[0]+bstep[0] , size);
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				for (k = bk[0]; k < (bk[0]+bstep[0]); ++k)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				for (k = bk[0]; k < i; ++k)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
}
#undef ind



#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
#define indvi(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)/nv-i, (size-1)-j))
#define indvj(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)-i, (size-1)/nv-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0A(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
	size_t size = X.sizeMem();
	size_t i, j, k, kv;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t imax, jmax, kmax;
	size_t bstep[5];
	size_t isrt;
	//const size_t unr = 2;
	//double acc[unr*unr];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;/**/
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);
	size_t nv = X.regEN();
	vec<double> acc;

#define vect(v) for(size_t v = 0; v < nv; ++v)
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size);
		jmax = min(bj[0]+bstep[0] , size);
		if(direction == Direction::Forwards)
			isrt = bi[0];
		else isrt = max(bi[0], X.pad());
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				assert( ((direction == Direction::Backwards) && (((size-1-bk[0])-(nv-1)) % 4 == 0))
				|| ((direction == Direction::Forwards) && (bk[0] % 4 == 0)));
				vect(v)
					acc[v] = 0;
				for (kv = bk[0]/nv; kv < (bk[0]+bstep[0])/nv; ++kv)
					acc.v = acc.v - indvj(LU, i, kv).v * indvi(X, kv, j).v;
				vect(v)
					ind(X,i,j) += acc[v];
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = isrt; i < imax; ++i)
			for (j = bj[0]; j < jmax; ++j) {
				for (k = bk[0]; k < i; ++k)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
#undef vect
#undef ind
#undef indv
}



#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
	size_t size = X.size();
	size_t i, j, k;
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);

	for(i = 0; i < size; i += 1){
		for(j = 0; j < size; j += 1){
			for(k = 0; k != i; k += 1){
				ind(X,i,j) = ind(X,i,j) - ind(LU,i,k) * ind(X,k,j);
			}
			if(diagonal == Diagonal::Value)
				ind(X,i,j) /= ind(LU,i,i);
		}
	}
}
#undef ind


#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))

#define indv(M,i,j) (direction == Direction::Forwards ? \
	M.atv(i, j) : \
	M.atv((size-1)-i, (size-1)-j))
template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLUA(LUMatrix& LU, XMatrix& X, BMatrix& B, varray<size_t>& P){
	size_t size = X.size();
	size_t i, j, k, kv;
	vec<double> acc;
	for(j = 0; j < size; ++j)
		for(i = 0; i < size; ++i)
			if(permute == Permute::True)
				ind(X, i, j) = ind(B, P.at(i), j);
			else
				ind(X, i, j) = ind(B, i, j);

#define vect(v) for(size_t v = 0; v < 4; v++)
	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j){
			vect(v)
				acc[v] = 0;
			for(kv = 0; kv < i/LU.regEN(); ++kv){
				acc.v = acc.v - indv(LU,i,kv).v * indv(X,kv,j).v;
			}
			vect(v)
				ind(X,i,j) += acc[v];
			for(k = kv*4; k < i; ++k){
				ind(X,i,j) = ind(X,i,j) - ind(LU,i,k) * ind(X,k,j);
			}
			if(diagonal == Diagonal::Value)
				ind(X,i,j) /= ind(LU,i,i);
		}
	}
}
#undef vect
#undef ind



template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLU0(LUMatrix& LU, IAMatrix& X, IMatrix& B, varray<size_t>& P){
	static
	MatrixColMajor<double> Z(X.size());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0A<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0A<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
}

#include <cmath>
#include <ctgmath>

template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0(AMatrix& A, IAMatrix& IA, IMatrix& I){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = 8;
	bstep[1] = bstep[0]*4;
	bstep[2] = bstep[1]*4;
	bstep[3] = bstep[2]*5;

	for(size_t j = 0; j < size; ++j){
		for(size_t i = 0; i < j; ++i)
			I.at(i,j) = 0;
		I.at(j,j) = 1;
		for(size_t i = j+1; i < size; ++i)
			I.at(i,j) = 0;
	}

	for (bi[3] = 0; bi[3] < size; bi[3] += bstep[3])
	for (bj[3] = 0; bj[3] < size; bj[3] += bstep[3])
	for (bk[3] = 0; bk[3] < size; bk[3] += bstep[3]){
		bimax[3] = min(bi[3]+bstep[3], size);
		bjmax[3] = min(bj[3]+bstep[3], size);
		bkmax[3] = min(bk[3]+bstep[3], size);
		for (bi[2] = bi[3]; bi[2] < bimax[3]; bi[2] += bstep[2])
		for (bj[2] = bj[3]; bj[2] < bjmax[3]; bj[2] += bstep[2])
		for (bk[2] = bk[3]; bk[2] < bkmax[3]; bk[2] += bstep[2]){
			bimax[2] = min(bi[2]+bstep[2], size);
			bjmax[2] = min(bj[2]+bstep[2], size);
			bkmax[2] = min(bk[2]+bstep[2], size);
			for (bi[1] = bi[2]; bi[1] < bimax[2]; bi[1] += bstep[1])
			for (bj[1] = bj[2]; bj[1] < bjmax[2]; bj[1] += bstep[1])
			for (bk[1] = bk[2]; bk[1] < bkmax[2]; bk[1] += bstep[1]){
				bimax[1] = min(bi[1]+bstep[1], size);
				bjmax[1] = min(bj[1]+bstep[1], size);
				bkmax[1] = min(bk[1]+bstep[1], size);
				for (bi[0] = bi[1]; bi[0] < bimax[1]; bi[0] += bstep[0])
				for (bj[0] = bj[1]; bj[0] < bjmax[1]; bj[0] += bstep[0])
				for (bk[0] = bk[1]; bk[0] < bkmax[1]; bk[0] += bstep[0]){
					size_t imax = min(bi[0]+bstep[0], size);
					size_t jmax = min(bj[0]+bstep[0], size);
					size_t kmax = min(bk[0]+bstep[0], size);
					for (size_t i = bi[0]; i < imax; ++i)
					for (size_t j = bj[0]; j < jmax; ++j) {
						for (size_t k = bk[0]; k < kmax; ++k) {
							I.at(i, j) = I.at(i, j) - A.at(i, k) * IA.at(k, j);
						}
					}
				}
			}
		}
	}
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < I.sizeVec(); ++iv)
			errNormV.v += I.atv(iv,j).v*I.atv(iv,j).v;
		for(size_t i = I.vecEnd(); i < I.size(); ++i)
			errNormV[I.regEN()-1] += I.at(i,j)*I.at(i,j);
	}
	for(size_t v=0; v < I.regEN(); ++v) errNorm += errNormV[v];

	return sqrt(errNorm);
}

/**
 * @brief Calculates inverse of A into IA
 * @param LU decomposition of A
 * @param IA return value, no init needed
 * @param P LU pivot permutation
 * @param iter_n
 */
template<class AMatrix, class LUMatrix, class IAMatrix>
void inverse_refining(AMatrix& A, LUMatrix& LU, IAMatrix& IA, varray<size_t>& P, size_t iter_n){
	long i=0;
	// number of digits of iter_n, for pretty printing
	long digits = (long)log10((double) iter_n) + 1;
	double c_residue;
	//double l_residue;
	size_t size = A.size();
	// Optm: iterating line by line
	MatrixColMajor<double> W(A.size()), R(A.size());

	for(size_t j = 0; j < size; ++j){
		for(size_t i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(size_t i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}

	//LIKWID_MARKER_START("INV");

	//solveMLU(LU, IA, R, P);
	solveMLU0(LU, IA, R, P);

	//LIKWID_MARKER_STOP("INV");
	//LIKWID_MARKER_START("RES");

	c_residue = residue0(A, IA, R);

	//LIKWID_MARKER_STOP("RES");

	while(i < iter_n){
		// (abs(l_residue - c_residue)/c_residue > EPSILON) && (l_residue > c_residue)
		// relative approximate error
		i += 1;
		// R: residue of IA

		//solveMLU(LU, W, R, P);
		solveMLU0(LU, W, R, P);

		//LIKWID_MARKER_STOP("INV");
		// W: residues of each variable of IA
		// adjust IA with found errors
		//LIKWID_MARKER_START("SUM");

		for(size_t j=0; j < IA.size(); ++j)
			for(size_t i=0; i < IA.size(); ++i)
				IA.at(i,j) += W.at(i,j);
		//add(IA, W);



		c_residue = residue0(A, IA, R);

		//LIKWID_MARKER_STOP("RES");
	}
}

void testVec(size_t size){
	Matrix<double> LU(size);
	MatrixColMajor<double> B(LU.size()), X(LU.size());
	size_t i,j,kv;

    size_t step = 8;
	vec<double> acc[step];

#define vect(v) for(size_t v = 0; v < step; v++)
    vect(s)
        acc[s][0] += 4;
	for(i=0; i < size; ++i)
	for(j=0; j < size; ++j)
	for(kv=0; kv < size/4; kv += 1){
		for(size_t s = 0; s < step; ++s)
            acc[s].v += B.atv(i,kv+s).v * X.atv(kv+s,j).v;
        for(size_t o = 0; o < 4; o++)
            acc[o][0] += o;
    }

	cout << acc[0][0] << endl;
#undef vect
}

int main(int argc, char **argv) {
	size_t size = 1024;
	size_t iter_n = 20;
	Matrix<double> LU(size);
	Matrix<double> A(size);
	MatrixColMajor<double> IA(LU.size());
	varray<size_t> P(LU.size());

	inverse_refining(A, LU, IA, P, iter_n);

	//testVec(1024);
    return 0;
}
