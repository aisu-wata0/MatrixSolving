/**
 * @brief Solves LU system using subst functions.
 * @param LU Matrix to find solution
 * @param X Variables to be found
 * Output: Value of found variables is stored in X
 * @param B Independent term matrix
 * @param P Permutation vector resulting of the pivoting
 * @param col Column of the matrix to be used as B
 */
template<class LUMatrix, class XMatrix, class BMatrix>
inline void solveLU(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P, long col){
	static
	varray<double> Z(LU.sizeMem());
	if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	subst<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P, col);
	// find X; Ux=Z
	subst<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P, col);
}
/**
 * @brief Finds inverse matrix,
 * also solves linear sistems A*IA = I where each of I's columns are different Bs
 * @param LU Matrix to find inverse
 * @param IA Inverse matrix to be found
 * Output: Inverse matrix is stored in IA
 * @param I Identity matrix
 * @param P Permutation vector resulting of the pivoting
 */
template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLU(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
	for(long j = 0; j < X.size(); j++){
		// for each X col solve SL to find the X col values
		static
		varray<double> Z(LU.sizeMem());
		if(Z.size() != X.size()){ Z.alloc(X.size()); }
		// find Z; LZ=B
		subst<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P, j);
		// find X; Ux=Z
		subst<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P, j);
	}
}

/**
 * @brief  Calculates residue into I, A*IA shuold be close to Identity
 * @param A original coef matrix
 * @param IA solution to A*IA = I
 * @param I Output residue, no init needed
 * @return Norm of the residue
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue(AMatrix& A, IAMatrix& IA, IMatrix& I){
	double err_norm = 0.0;
	size_t size = A.size();

	for(size_t j = 0; j < size; ++j){
		// for each column of the inverse
		for(size_t i = 0; i < size; i++){
			// for each line of A
			I.at(i,j) = 0;
			// multiply A line to the current inverse col
			for(size_t k = 0; k < size; ++k){
				I.at(i,j) = I.at(i,j) - A.at(i,k) * IA.at(k,j);
			}
			if(i == j){
				I.at(i,j) += 1;
			}
			
			err_norm += I.at(i,j)*I.at(i,j);
		}
	}
	return sqrt(err_norm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \n
 * Tiling on L0
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = 24;
	bstep[1] = bstep[0]*4;
	bstep[2] = bstep[1]*4;
	bstep[3] = bstep[2]*5;
	
	size_t i, j, k, kv;
	// set R to identity
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// Multiply A*IA and subtract from R
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j)
		for (k = bk[0]; k < kmax; ++k){
			R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
		}
	}
	// Calculate norm error from R
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv) // vect loop
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.remStart(); i < R.size(); ++i) // vect remainder
			errNormV[R.vecN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.vecN(); ++v) // vect result sum
		errNorm += errNormV[v];
	
	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \n
 * Tiling on L0, SSE
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0A(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${32} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	size_t i, j, k, kv;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// Multiply A*IA and subtract from R
	size_t vn = R.vecN();
	
#define vect(v) for(size_t v=0; v < vn; ++v)
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j) {
			vec<double> acc;
			vect(v) acc[v] = 0;
			for (kv = bk[0]/vn; kv < kmax/vn; ++kv)
				acc.v = acc.v - A.atv(i, kv).v * IA.atv(kv, j).v;
			for(k = kv*vn; k < kmax; ++k)
				R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
			vect(v) R.at(i, j) += acc[v];
		}
	}
#undef vect
	// Calculate norm error from R
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv) // vect loop
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.remStart(); i < R.size(); ++i) // vect remainder
			errNormV[R.vecN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.vecN(); ++v) // vect result sum
		errNorm += errNormV[v];
	
	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \n
 * Tiling on L0, SSE, Unrolling on k
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AU(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	const size_t kunr = 8;
	vec<double> acc;
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/

	size_t i, j, k, kv, rem;

	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// Multiply A*IA and subtract from R
	size_t vn = R.vecN();
	
#define vect(v) for(size_t v=0; v < vn; ++v)
#define unr(u,n) for(size_t u = 0; u < n; ++u)
	
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; ++i)
		for (j = bj[0]; j < jmax; ++j) {
			vect(u) vect(v) acc[v] = 0;
			for (kv = bk[0]/vn; kv < kmax/vn -(kunr-1); kv += kunr)
				unr(u,kunr) acc.v += A.atv(i, kv+u).v * IA.atv(kv+u, j).v;
			for(k = kv*vn; k < kmax; ++k) // vect remainder
				R.at(i, j) = R.at(i, j) - A.at(i, k) * IA.at(k, j);
			vect(v) R.at(i, j) -= acc[v]; // vect result sum
		}
	}
#undef vect
#undef unr
	// Calculate norm error from R
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv) // vect loop
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.remStart(); i < R.size(); ++i) // vect remainder
			errNormV[R.vecN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.vecN(); ++v) // vect result sum
		errNorm += errNormV[v];

	return sqrt(errNorm);
}
/**
 * @brief Given a matrix A and it's inverse, calculates residue into R. \n
 * Tiling on L0, SSE, Unrolling on i,j (Doen't care for the remainder of the unrolling, unnacurate)
 */
template<class AMatrix, class IAMatrix, class IMatrix>
inline double residue0AUU(AMatrix& A, IAMatrix& IA, IMatrix& R){
	size_t size = A.size();
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	const size_t unr = 2;
	vec<double> acc[unr*unr];
	/**/
	bstep[0] = B2L1;
	bstep[1] = bstep[0]*3;
	/* export GCC_ARGS=" -D L0=${24} -D L1M=${3}"*
	bstep[0] = L0;
	bstep[1] = bstep[0]*L1M;/**/
	size_t i, j, k, kv, rem;
	for(j = 0; j < size; ++j){
		for(i = 0; i < j; ++i)
			R.at(i,j) = 0;
		R.at(j,j) = 1;
		for(i = j+1; i < size; ++i)
			R.at(i,j) = 0;
	}
	// Multiply A*IA and subtract from R
	size_t vn = R.vecN();
	
#define vect(v) for(size_t v=0; v < vn; ++v)
#define unr(u,n) for(size_t u = 0; u < n; ++u)
#define unr2(iu,ju,n) unr(iu,n) unr(ju,n)

	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0])
	for (bk[0] = 0; bk[0] < size; bk[0] += bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		for (i = bi[0]; i < imax; i += unr)
		for (j = bj[0]; j < jmax; j += unr) {
			unr2(iu,ju,unr) vect(v) acc[iu*unr + ju][v] = 0;
			for (kv = bk[0]/vn; kv < kmax/vn; ++kv)
				unr2(iu,ju,unr)
					acc[iu*unr+ju].v += A.atv(i+iu, kv).v * IA.atv(kv, j+ju).v;
			for(k = kv*vn; k < kmax; ++k) // vect remainder
				unr2(iu,ju,unr)
					R.at(i+iu, j+ju) -= A.at(i+iu, k) * IA.at(k, j+ju);
			unr2(iu,ju,unr)
				vect(v) R.at(i+iu, j+ju) -= acc[iu*unr+ju][v]; // vect result sum
		}
	}
#undef vect
#undef unr
#undef unr2
	// Calculate norm error from R
	double errNorm = 0;
	vec<double> errNormV{0};
	for(size_t j = 0; j < size; ++j){
		for(size_t iv = 0; iv < R.sizeVec(); ++iv) // vect loop
			errNormV.v += R.atv(iv,j).v*R.atv(iv,j).v;
		for(size_t i = R.remStart(); i < R.size(); ++i) // vect remainder
			errNormV[R.vecN()-1] += R.at(i,j)*R.at(i,j);
	}
	for(size_t v=0; v < R.vecN(); ++v) // vect result sum
		errNorm += errNormV[v];

	return sqrt(errNorm);
}