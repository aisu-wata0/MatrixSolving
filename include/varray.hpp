#ifndef VARRAY_HPP
#define VARRAY_HPP

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "bytes.h"

#define unroll(v,n) for(size_t v = 0; v < n; v++)
#define unroll2(vi,vj,n) unroll(vi,n)unroll(vj,n)
#define vec(v) for(size_t v = 0; v < regDN; v++)

namespace gm
{
using namespace std;


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
	
	const elem& operator[] (size_t i) const {
		assert(i < regN<elem>() && "Vectorized elem out of register access");
		return v[i];
	}
	elem& operator[] (size_t i) {
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
		mSizeMem += mod(-mSize,(long)L1_LINE_DN);
		//mSizeMem += mod(-mSize, (ptrdiff_t)L1_LINE_DN); TODO test
		
		if(((mSizeMem/L1_LINE_DN) % 2) == 0){
			mSizeMem = mSizeMem + L1_LINE_DN; // make sure mSizeMem is odd multiple
		}
		
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
	vec<elem>& atv(size_t i) {
		assert(i < mSizeVec && "varray vec access out of bounds");
		return arr.v[i];
	}
	const vec<elem> & atv(size_t i) const {
		assert(i < mSizeVec && "varray vec access out of bounds");
		return arr.v[i];
	}
	
	elem& at(size_t i){
		assert(i < mSize && "varray access out of bounds");
		return arr.p[i];
	}
	const elem& at(size_t i) const {
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