#ifndef VARRAY_H
#define VARRAY_H
/**
@file varray.h
*/

#include <stdint.h>

namespace std {

long div_down(long n, long d) {
    return n / d - (((n > 0) ^ (d > 0)) && (n % d));
}

#define mod(X,Y) ((((X) % (Y)) + (Y)) % Y)

#define CACHE_LINE_SIZE (64) // likwid-topology: Cache line size:	64
#define L1_LINE_DN (CACHE_LINE_SIZE/sizeof(double)) // how many doubles in a line

#define CACHE_L1_SIZE (32*1024/2) // likwid-topology: Size:	 32 kB
// divided by 2 because we wont be able to fill L1 completely without throwing
// useful values out
#define L1_DN (CACHE_L1_SIZE/sizeof(double)) // how many doubles in L1 cache
// 2048
// size of the block that should fit into L1
#define BL1 ((long)((sqrt(L1_DN)) - mod(((long)sqrt(L1_DN)), L1_LINE_DN)))
// 40 % 8 == 0, (40*40 < 2048)

#define REG_SZ (32) // how many bytes in a register
#define dn (REG_SZ/sizeof(double)) // how many doubles is a register
// 4

template<typename type>
struct vec
{
  type __attribute__ ((vector_size (REG_SZ)))  v; // vectorization of doubles
};

template<typename type>
union vp
{
  vec<type>*  v;
  type* d;
};

#define unroll(v,n) for(size_t v = 0; v < n; v++)
#define unroll2(vi,vj,n) unroll(vi,n)unroll(vj,n)
#define vec(v) for(size_t v = 0; v < dn; v++)

/**
 * @brief Allocates size*sizeof(type) and aligns it into a boundary-bit boundary
 * @param tpp pointer to set the array
 * @param size number of elements in your resulting pointer
 * @param boundary : Power of two, else the behavior is undefined
 * @return The pointer you should free, don't lose it
 */
template<typename type>
void* al_allloc(type** tpp, size_t size, size_t boundary){
	*tpp = (type*)malloc((size)*sizeof(type) + boundary-1);
	if(*tpp == NULL){
		cerr <<"failed to malloc "<< (size)*sizeof(type)/1024 <<" KiB"<< endl;
		exit(0);
	}
	void* mem_p = *tpp;
	*tpp = (type*)(((uintptr_t)mem_p + (boundary-1)) & (uintptr_t)( ~ (boundary-1)));
	return mem_p;
}

template<typename type>
class varray
{
public:
	vp<type> arr;
	size_t size; //  size of doubles
	size_t m_size; // size of doubles in memory
	size_t m_size_v; // size of vectorized doubles in memory
	void* mem_p;
	
	void mem_alloc(long m_size){
		mem_p = al_allloc(&arr.d, m_size, CACHE_LINE_SIZE);
	}
	
	void alloc_size(size_t size){
		this->size = size;
		m_size = size;
		
		m_size += mod(-size,(long)BL1); // add to make it multiple of BL1
		
		if(((m_size/L1_LINE_DN) % 2) == 0){
			m_size = m_size + BL1; // make sure m_size is odd multiple
		}
		m_size_v = m_size/dn;
		mem_alloc(m_size);
	}
	
	varray() {
		mem_p = NULL;
	}
	
	varray(size_t size) {
		mem_p = NULL;
		alloc_size(size);
	}
	
	~varray(){
		if(mem_p != NULL)
			{ free(mem_p); }
	}
	
	vec<type>& atv(long i) {
		return arr.v[i];
	}
	const vec<type> & atv(long i) const {
		return arr.v[i];
	}
	type& at(long i){
		return arr.d[i];
	}
	const type& at(long i) const {
		return arr.d[i];
	}
};


}
#endif