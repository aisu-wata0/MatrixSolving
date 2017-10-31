
#include <cstring>
#include "Chronometer.hpp"
#include "Subst.hpp"

#define block(i,i0,imax, j,j0,jmax, k,k0,kmax, bstep) \
for (i = i0; i < imax; i += bstep) \
	for (j = j0; j < jmax; j += bstep) \
		for (k = k0; k < kmax; k += bstep)

//inline void test(Matrix& LU, MatrixColMajor& B, Matrix& X){
void t_matrix_mult(size_t size){
	Matrix<double> LU(size);
	MatrixColMajor<double> X(size);
	MatrixColMajor<double> B(size);
	size_t i, j, k;
	
	printf("%llu\n", size);
	
	asm("SETUP");
	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			if(true){
				LU.at(i,j) = i*1000 + j;
				B.at(j,i) = j*1000 + i;
			} else {
				LU.at(i,j) = j;
				B.at(j,i) = (double)i/size;
			}
		}
	}
	size_t repetitionN = 3;
	size_t reps;
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	/**
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/**/
	bstep[0] = B3L1;
	bstep[1] = bstep[0]*4;
	bstep[2] = bstep[1]*4;
	bstep[3] = bstep[2]*5;
	/**
	bstep[1] = bstep[3]*2;
	bstep[2] = bstep[3]*2;
	bstep[3] = bstep[3]*2;
	/**
	bstep[0] = 2;
	bstep[1] = bstep[0]*2;
	bstep[2] = bstep[1]*2;
	bstep[3] = bstep[2]*2;
	/**/
	// export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} "
	//bstep[1] = bstep[0]*L1M;
	//bstep[2] = bstep[1]*L2M;
	//bstep[2] = bstep[1]*L3M;
	
	printf("%llu, %llu, %llu, %llu\n", bstep[0], bstep[1], bstep[2], bstep[3]);
	const bool PRINT_MATRIX = false;
	//const size_t unr = 2;
	//double acc[unr*unr];
	
	Chronometer<128> timer;
	/**/
	// warmup
	set(X,0);
	block(bi[0],0,size, bj[0],0,size, bk[0],0,size, bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
			X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
		}
	}
	// end warmup
	
	
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("Tiled0");
	block(bi[0],0,size, bj[0],0,size, bk[0],0,size, bstep[0]){
		size_t imax = min(bi[0]+bstep[0], size);
		size_t jmax = min(bj[0]+bstep[0], size);
		size_t kmax = min(bk[0]+bstep[0], size);
		block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
			X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
		}
	}
	asm("END Tiled0");
	}
	printf("Tiled0  \t%f sec\n", timer.tick()/repetitionN);
	//if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("Tiled1.0");
	block(bi[1],0,size, bj[1],0,size, bk[1],0,size, bstep[1]){
		bimax[1] = min(bi[1]+bstep[1], size);
		bjmax[1] = min(bj[1]+bstep[1], size);
		bkmax[1] = min(bk[1]+bstep[1], size);
		block(bi[0],bi[1],bimax[1], bj[0],bj[1],bjmax[1], bk[0],bk[1],bkmax[1], bstep[0]){
			size_t imax = min(bi[0]+bstep[0], size);
			size_t jmax = min(bj[0]+bstep[0], size);
			size_t kmax = min(bk[0]+bstep[0], size);
			block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
				X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
			}
		}
	}
	asm("END Tiled1.0");
	}
	printf("Tiled1.0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("Tiled2.1.0");
	block(bi[2],0,size, bj[2],0,size, bk[2],0,size, bstep[2]){
		bimax[2] = min(bi[2]+bstep[2], size);
		bjmax[2] = min(bj[2]+bstep[2], size);
		bkmax[2] = min(bk[2]+bstep[2], size);
		block(bi[1],bi[2],bimax[2], bj[1],bj[2],bjmax[2], bk[1],bk[2],bkmax[2], bstep[1]){
			bimax[1] = min(bi[1]+bstep[1], size);
			bjmax[1] = min(bj[1]+bstep[1], size);
			bkmax[1] = min(bk[1]+bstep[1], size);
			block(bi[0],bi[1],bimax[1], bj[0],bj[1],bjmax[1], bk[0],bk[1],bkmax[1], bstep[0]){
				size_t imax = min(bi[0]+bstep[0], size);
				size_t jmax = min(bj[0]+bstep[0], size);
				size_t kmax = min(bk[0]+bstep[0], size);
				block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
					X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
				}
			}
		}
	}
	asm("END Tiled2.1.0");
	}
	printf("Tiled2.1.0 \t%f sec\n", timer.tick()/repetitionN);
	//if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("Tiled3.2.1.0");
	block(bi[3],0,size, bj[3],0,size, bk[3],0,size, bstep[3]){
		bimax[3] = min(bi[3]+bstep[3], size);
		bjmax[3] = min(bj[3]+bstep[3], size);
		bkmax[3] = min(bk[3]+bstep[3], size);
		block(bi[2],bi[3],bimax[3], bj[2],bj[3],bjmax[3], bk[2],bk[3],bkmax[3], bstep[2]){
			bimax[2] = min(bi[2]+bstep[2], size);
			bjmax[2] = min(bj[2]+bstep[2], size);
			bkmax[2] = min(bk[2]+bstep[2], size);
			block(bi[1],bi[2],bimax[2], bj[1],bj[2],bjmax[2], bk[1],bk[2],bkmax[2], bstep[1]){
				bimax[1] = min(bi[1]+bstep[1], size);
				bjmax[1] = min(bj[1]+bstep[1], size);
				bkmax[1] = min(bk[1]+bstep[1], size);
				block(bi[0],bi[1],bimax[1], bj[0],bj[1],bjmax[1], bk[0],bk[1],bkmax[1], bstep[0]){
					size_t imax = min(bi[0]+bstep[0], size);
					size_t jmax = min(bj[0]+bstep[0], size);
					size_t kmax = min(bk[0]+bstep[0], size);
					block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
						X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
					}
				}
			}
		}
	}
	asm("END Tiled3.2.1.0");
	}
	printf("Tiled3.2.1.0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	////
	// LU solving
	////
	/**/
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	set(X,B);
	asm("LUTiled0");
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		//bkmax[0] = bi[0]+bstep[0];
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			size_t imax = min(bi[0]+bstep[0], size);
			size_t jmax = min(bj[0]+bstep[0], size);
			size_t kmax = bk[0]+bstep[0];
			block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
				X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
			}
		} // Last block in K, diagonal
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			size_t imax = min(bi[0]+bstep[0], size);
			size_t jmax = min(bj[0]+bstep[0], size);
			for (size_t i = bi[0]; i < imax; ++i)
			for (size_t j = bj[0]; j < jmax; ++j) {
				for (size_t k = bk[0]; k < i; ++k) {
					X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
				}
				X.at(i, j) /= LU.at(i, i);
			}
		}
	}
	asm("END LUTiled0");
	}
	printf("LUTiled0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	set(X,B);
	asm("LUTiled1.0");
	block (bi[1],0,size, bj[1],0,size, bk[1],0,(bi[1]+bstep[1]), bstep[1]) {
		bimax[0] = min(bi[1]+bstep[1], size);
		bjmax[0] = min(bj[1]+bstep[1], size);
		for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
		for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
			bkmax[0] = min(bk[1]+bstep[1], bi[0]);
			for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
				size_t imax = min(bi[0]+bstep[0], size);
				size_t jmax = min(bj[0]+bstep[0], size);
				size_t kmax = bk[0]+bstep[0];
				block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
					X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
				}
			} // Last block in K, diagonal
			if(bk[0] == bi[0] && bk[1] == bi[1]){
				for (; bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
					size_t imax = min(bi[0]+bstep[0], size);
					size_t jmax = min(bj[0]+bstep[0], size);
					for (size_t i = bi[0]; i < imax; ++i)
					for (size_t j = bj[0]; j < jmax; ++j) {
						for (size_t k = bk[0]; k < i; ++k) {
							X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
						}
						X.at(i, j) /= LU.at(i, i);
					}
				}
			}
		}
	}
	asm("END LUTiled1.0");
	}
	printf("LUTiled1.0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,B);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("LUTiled2.1.0");
	block (bi[2],0,size, bj[2],0,size, bk[2],0,(bi[2]+bstep[2]), bstep[2]) {
		bimax[1] = min(bi[2]+bstep[2], size);
		bjmax[1] = min(bj[2]+bstep[2], size);
		for (bi[1] = bi[2]; bi[1] < bimax[1]; bi[1] += bstep[1])
		for (bj[1] = bj[2]; bj[1] < bjmax[1]; bj[1] += bstep[1]) {
			bkmax[1] = min(bk[2]+bstep[2], bi[1]+bstep[1]);
			for (bk[1] = bk[2]; bk[1] < bkmax[1]; bk[1] += bstep[1]) {
				bimax[0] = min(bi[1]+bstep[1], size);
				bjmax[0] = min(bj[1]+bstep[1], size);
				for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
				for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
					bkmax[0] = min(bk[1]+bstep[1], bi[0]);
					for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
						size_t imax = min(bi[0]+bstep[0], size);
						size_t jmax = min(bj[0]+bstep[0], size);
						size_t kmax = bk[0]+bstep[0];
						block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
							X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
						}
					} // Last block in K, diagonal
					if(bk[0] == bi[0] && bk[1] == bi[1] && bk[2] == bi[2]){
						for (; bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
							size_t imax = min(bi[0]+bstep[0], size);
							size_t jmax = min(bj[0]+bstep[0], size);
							for (size_t i = bi[0]; i < imax; ++i)
							for (size_t j = bj[0]; j < jmax; ++j) {
								for (size_t k = bk[0]; k < i; ++k) {
									X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
								}
								X.at(i, j) /= LU.at(i, i);
							}
						}
					}
				}
			}
		}
	}
	asm("END LUTiled2.1.0");
	}
	printf("LUTiled2.1.0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,B);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
	asm("LUTiled3.2.1.0");
	block (bi[3],0,size, bj[3],0,size, bk[3],0,(bi[3]+bstep[3]), bstep[3]) {
		for (bi[2] = bi[3]; bi[2] < bimax[2]; bi[2] += bstep[2])
		for (bj[2] = bj[3]; bj[2] < bjmax[2]; bj[2] += bstep[2]) {
			bkmax[2] = min(bk[3]+bstep[3], bi[2]+bstep[2]);
			for (bk[2] = bk[3]; bk[2] < bkmax[2]; bk[2] += bstep[2]) {
				bimax[1] = min(bi[2]+bstep[2], size);
				bjmax[1] = min(bj[2]+bstep[2], size);
				for (bi[1] = bi[2]; bi[1] < bimax[1]; bi[1] += bstep[1])
				for (bj[1] = bj[2]; bj[1] < bjmax[1]; bj[1] += bstep[1]) {
					bkmax[1] = min(bk[2]+bstep[2], bi[1]+bstep[1]);
					for (bk[1] = bk[2]; bk[1] < bkmax[1]; bk[1] += bstep[1]) {
						bimax[0] = min(bi[1]+bstep[1], size);
						bjmax[0] = min(bj[1]+bstep[1], size);
						for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
						for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
							bkmax[0] = min(bk[1]+bstep[1], bi[0]);
							for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
								size_t imax = min(bi[0]+bstep[0], size);
								size_t jmax = min(bj[0]+bstep[0], size);
								size_t kmax = bk[0]+bstep[0];
								block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
									X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
								}
							} // Last block in K, diagonal
							if(bk[0] == bi[0] && bk[1] == bi[1] && bk[2] == bi[2]){
								for (; bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
									size_t imax = min(bi[0]+bstep[0], size);
									size_t jmax = min(bj[0]+bstep[0], size);
									for (size_t i = bi[0]; i < imax; ++i)
									for (size_t j = bj[0]; j < jmax; ++j) {
										for (size_t k = bk[0]; k < i; ++k) {
											X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
										}
										X.at(i, j) /= LU.at(i, i);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	asm("END LUTiled3.2.1.0");
	}
	printf("LUTiled3.2.1.0 \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	/**
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		set(X,B);
		asm("LU_Naive");
		for(i=0; i < size; i += step){
			for(j=0; j < size; j += step){
				for(k=0; k < i; k += 1){
					X.at(i,j) = X.at(i,j) - LU.at(i,k) * X.at(k,j);
				}
				X.at(i,j) /= LU.at(i,i);
			}
		}
		asm("END LU_Naive");
	}
	printf("LU_Naive: \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	
	/**
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		asm("Naive");
		for(i=0; i < size; i += step){
			for(j=0; j < size; j += step){
				for(k=0; k < size; k += 1){	
					X.at(i,j) = X.at(i,j) + LU.at(i,k) * B.at(k,j);
				}
			}
		}
		asm("END Naive");
	}
	printf("Naive: \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	
	/**
	set(X,0);
	timer.tick();
	asm("Tile_BI");
	for (bi = 0; bi < size; bi += bstep[0]){
		for (j = 0; j < size; j += unr){
			for (i = bi; i < (bi+bstep[0]); i += unr){
				unroll2(ci,cj,unr) acc[ci*unr + cj] = 0;
				for (k = 0; k < size; k++){
					unroll2(ci,cj,unr) acc[ci*unr + cj] += LU.at(i+ci, k) * at(B, k, j+cj);
				}
				unroll2(ci,cj,unr) at(X, i+ci, j+ci) = acc[ci*unr + cj];
			}
		}
	}
	asm("END Tile_BI");
	printf("Tile_BI: \t%f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
}

#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))

template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;	
	if(permute == Permute::True)
		set<direction>(X,B,P);
	else set<direction>(X,B);
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

template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU0(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
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
	
	asm("LUTiled0");
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		imax = min(bi[0]+bstep[0] , size);
		jmax = min(bj[0]+bstep[0] , size);
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < (bk[0]+bstep[0]); k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < i; k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
	asm("END LUTiled0");
}


template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
inline void substMLU10(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
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
	
	asm("LUTiled1.0");
	block (bi[1],0,size, bj[1],0,size, bk[1],0,(bi[1]+bstep[1]), bstep[1]) {
		bimax[0] = min(bi[1]+bstep[1], size);
		bjmax[0] = min(bj[1]+bstep[1], size);
		for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
		for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
			bkmax[0] = min(bk[1]+bstep[1], bi[0]);
			imax = min(bi[0]+bstep[0] , size);
			jmax = min(bj[0]+bstep[0] , size);
			for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
				for (i = bi[0]; i < imax; i += 1)
				for (j = bj[0]; j < jmax; j += 1) {
					for (k = bk[0]; k < (bk[0]+bstep[0]); k += 1)
						ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				}
			} // Last block in K, diagonal
			if(bk[0] == bi[0] && bk[1] == bi[1])
			for (bk[0] = bi[0]; bk[0] < (bi[0]+bstep[0]); bk[0] += bstep[0]) {
				for (i = bi[0]; i < imax; i += 1)
				for (j = bj[0]; j < jmax; j += 1) {
					for (k = bk[0]; k < i; k += 1)
						ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
					if(diagonal == Diagonal::Value)
						ind(X, i, j) /= ind(LU, i, i);
				}
			}
		}
	}
	asm("END LUTiled1.0");
}
#undef ind

template<class LUMatrix, class IAMatrix, class IMatrix>
inline void solveMLUNew(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
	//static
	MatrixColMajor<double> Z(X.size());
	//if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
}