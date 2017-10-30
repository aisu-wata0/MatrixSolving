
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
	size_t step;
	
	size = LU.size();
	cout << size << ":" << LU.sizeMem() << endl;
	
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
	
	step = +1;
	//bstep = step * BL1/3; // cache block size
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
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
	size_t block_n[16];
	
	const size_t unr = 2;
	const bool PRINT_MATRIX = false;
	double acc[unr*unr];
	
	Chronometer<128> timer;
	/**
	// warmup
	set(X,0);
	timer.tick();
	for (size_t i0 = 0; i0 < size; i0 += bstep[0]) {
		for (size_t j0 = 0; j0 < size; j0 += bstep[0]) {
			for (size_t k0 = 0; k0 < size; k0 += bstep[0]) {
				size_t imax = i0 + bstep[0] > size ? size : i0 + bstep[0];
				size_t jmax = j0 + bstep[0] > size ? size : j0 + bstep[0];
				size_t kmax = k0 + bstep[0] > size ? size : k0 + bstep[0];
				for (size_t j = j0; j < jmax; ++j) {
					for (size_t i = i0; i < imax; ++i) {
						for (size_t k = k0; k < kmax; ++k) {
							X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
						}
					}
				}
			}
		}
	}
	// end warmup
	
	
	memset(block_n, 0, sizeof(block_n));
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		asm("Tiled_0");
		block(bi[0],0,size, bj[0],0,size, bk[0],0,size, bstep[0]){
			block_n[0]++;
			size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
			size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
			size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];
			block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
				block_n[4]++;
				X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
			}
		}
		asm("END Tiled_0");
	}
	printf("Tiled_0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	memset(block_n, 0, sizeof(block_n));
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		asm("Tiled_1.0");
		block(bi[1],0,size, bj[1],0,size, bk[1],0,size, bstep[1]){
		block_n[1]++;
			size_t i1max = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
			size_t j1max = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
			size_t k1max = bk[1] + bstep[1] > size ? size : bk[1] + bstep[1];
			block(bi[0],bi[1],i1max, bj[0],bj[1],j1max, bk[0],bk[1],k1max, bstep[0]){
				block_n[0]++;
				size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
				size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
				size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];
				block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
					block_n[4]++;
					X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
				}
			}
		}
		asm("END Tiled_1.0");
	}
	printf("Tiled_1.0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	memset(block_n, 0, sizeof(block_n));
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		asm("Tiled_2.1.0");
		block(bi[2],0,size, bj[2],0,size, bk[2],0,size, bstep[2]){
			block_n[2]++;
			size_t i2max = bi[2] + bstep[2] > size ? size : bi[2] + bstep[2];
			size_t j2max = bj[2] + bstep[2] > size ? size : bj[2] + bstep[2];
			size_t k2max = bk[2] + bstep[2] > size ? size : bk[2] + bstep[2];
			block(bi[1],bi[2],i2max, bj[1],bj[2],j2max, bk[1],bk[2],k2max, bstep[1]){
				block_n[1]++;
				size_t i1max = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
				size_t j1max = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
				size_t k1max = bk[1] + bstep[1] > size ? size : bk[1] + bstep[1];
				block(bi[0],bi[1],i1max, bj[0],bj[1],j1max, bk[0],bk[1],k1max, bstep[0]){
					block_n[0]++;
					size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
					size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
					size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];
					block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
						block_n[4]++;
						X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
					}
				}
			}
		}
		asm("END Tiled_2.1.0");
	}
	printf("Tiled_2.1.0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	
	memset(block_n, 0, sizeof(block_n));
	set(X,0);
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		asm("Tiled_3.2.1.0");
		block(bi[3],0,size, bj[3],0,size, bk[3],0,size, bstep[3]){
			block_n[3]++;
			size_t i3max = bi[3] + bstep[3] > size ? size : bi[3] + bstep[3];
			size_t j3max = bj[3] + bstep[3] > size ? size : bj[3] + bstep[3];
			size_t k3max = bk[3] + bstep[3] > size ? size : bk[3] + bstep[3];
			block(bi[2],bi[3],i3max, bj[2],bj[3],j3max, bk[2],bk[3],k3max, bstep[2]){
				block_n[2]++;
				size_t i2max = bi[2] + bstep[2] > size ? size : bi[2] + bstep[2];
				size_t j2max = bj[2] + bstep[2] > size ? size : bj[2] + bstep[2];
				size_t k2max = bk[2] + bstep[2] > size ? size : bk[2] + bstep[2];
				block(bi[1],bi[2],i2max, bj[1],bj[2],j2max, bk[1],bk[2],k2max, bstep[1]){
					block_n[1]++;
					size_t i1max = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
					size_t j1max = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
					size_t k1max = bk[1] + bstep[1] > size ? size : bk[1] + bstep[1];
					block(bi[0],bi[1],i1max, bj[0],bj[1],j1max, bk[0],bk[1],k1max, bstep[0]){
						block_n[0]++;
						size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
						size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
						size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];
						block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
							block_n[4]++;
							X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
						}
					}
				}
			}
		}
		asm("END Tiled_3.2.1.0");
	}
	printf("Tiled_3.2.1.0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	
	size_t test[5];
	////
	// LU solving
	////
	memset(test, 0, sizeof(test));
	memset(block_n, 0, sizeof(block_n));
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		set(X,B);
		asm("LU_Tiled_0");
		for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
		for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
			test[0]++;
			//bkmax[0] = bi[0] + bstep[0];
			for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
				block_n[0]++;
				size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
				size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
				size_t kmax = bk[0] + bstep[0];
				block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
					block_n[6]++;
					X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
				}
			}
			// Last block in K, diagonal
			for (bk[0] = (bi[0]); bk[0] < (bi[0] + bstep[0]); bk[0] += bstep[0]) {
				block_n[2]++;
				size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
				size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
				for (size_t i = bi[0]; i < imax; ++i)
				for (size_t j = bj[0]; j < jmax; ++j) {
					for (size_t k = bk[0]; k < i; ++k) {
						block_n[8]++;
						X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
					}
					X.at(i, j) /= LU.at(i, i);
				}
			}
		}
		asm("END LU_Tiled_0");
	}
	printf("LU_Tiled_0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	memset(test, 0, sizeof(test));
	memset(block_n, 0, sizeof(block_n));
	timer.tick();
	for(reps = 0; reps < repetitionN; reps++){
		set(X,B);
		asm("LU_Tiled_1.0");
		block (bi[1],0,size, bj[1],0,size, bk[1],0,(bi[1]+bstep[1]), bstep[1]) {
			test[0] = 0;
			block_n[15]++;
			bimax[0] = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
			bjmax[0] = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
			for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
			for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
				test[0]++;
				bkmax[0] = min(bk[1] + bstep[1], bi[0]);
				for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
					block_n[0]++;
					size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
					size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
					size_t kmax = bk[0] + bstep[0];
					block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
						block_n[6]++;
						X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
					}
				}
				// Last block in K, diagonal
				if(bk[0] == bi[0] && bk[1] == bi[1]){
					for (; bk[0] < (bi[0] + bstep[0]); bk[0] += bstep[0]) {
						block_n[2]++;
						size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
						size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
						for (size_t i = bi[0]; i < imax; ++i)
						for (size_t j = bj[0]; j < jmax; ++j) {
							for (size_t k = bk[0]; k < i; ++k) {
								block_n[8]++;
								X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
							}
							X.at(i, j) /= LU.at(i, i);
						}
					}
				}
			}
		}
		asm("END LU_Tiled_1.0");
	}
	printf("LU_Tiled_1.0: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	/**/
	
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
	printf("LU_Naive: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	
//	memset(test, 0, sizeof(test));
//	memset(block_n, 0, sizeof(block_n));
//	set(X,B);
//	timer.tick();
//	for(reps = 0; reps < repetitionN; reps++){
//		asm("LU_Tiled_2.1.0");
//		//block (bi[2],0,size, bj[2],0,size, bk[2],0,(bi[2]+bstep[2]), bstep[2]) {
//			block (bi[1],0,size, bj[1],0,size, bk[1],0,(bi[1]+bstep[1]), bstep[1]) {
//				test[0] = 0;
//				block_n[15]++;
//				bimax[0] = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
//				bjmax[0] = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
//				for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
//				for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
//					test[0]++;
//					bkmax[0] = min(bk[1] + bstep[1], bi[0]);
//					for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
//						block_n[0]++;
//						size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
//						size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
//						size_t kmax = bk[0] + bstep[0];
//						block(i,bi[0],imax, j,bj[0],jmax, k,bk[0],kmax, 1){
//							block_n[6]++;
//							X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
//						}
//					}
//					// Last block in K, diagonal
//					if(bk[0] == bi[0] && bk[1] == bi[1]){
//						for (; bk[0] < (bi[0] + bstep[0]); bk[0] += bstep[0]) {
//							block_n[2]++;
//							size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
//							size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
//							for (size_t i = bi[0]; i < imax; ++i)
//							for (size_t j = bj[0]; j < jmax; ++j) {
//								for (size_t k = bk[0]; k < i; ++k) {
//									block_n[8]++;
//									X.at(i, j) = X.at(i, j) - LU.at(i, k) * X.at(k, j);
//								}
//								X.at(i, j) /= LU.at(i, i);
//							}
//						}
//					}
//				}
//			}
//		//}
//		asm("END LU_Tiled_2.1.0");
//	}
//	printf("LU_Tiled_2.1.0: %f sec\n", timer.tick()/repetitionN);
//	if(PRINT_MATRIX) { printm(X); cout << endl; }
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
	printf("Naive: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	/**
	set(X,0);
	timer.tick();
	asm("Tile_BI");
	for (bi = 0; bi < size; bi += bstep[0]){
		for (j = 0; j < size; j += unr){
			for (i = bi; i < (bi + bstep[0]); i += unr){
				unroll2(ci,cj,unr) acc[ci*unr + cj] = 0;
				for (k = 0; k < size; k++){
					unroll2(ci,cj,unr) acc[ci*unr + cj] += LU.at(i+ci, k) * at(B, k, j+cj);
				}
				unroll2(ci,cj,unr) at(X, i+ci, j+ci) = acc[ci*unr + cj];
			}
		}
	}
	asm("END Tile_BI");
	printf("Tile_BI: %f sec\n", timer.tick()/repetitionN);
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
}

#define ind(M,i,j) (direction == Direction::Forwards ? \
	M.at(i, j) : \
	M.at((size-1)-i, (size-1)-j))

template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
void substMLU(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;
	/**/
	const size_t unr = 2;
	double acc[unr*unr];
	
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
void substMLU0(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	//size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;
	/**/
	//const size_t unr = 2;
	//double acc[unr*unr];
	
	if(permute == Permute::True)
		set<direction>(X,B,P);
	else set<direction>(X,B);
	asm("LU_Tiled_0");
	for (bi[0] = 0; bi[0] < size; bi[0] += bstep[0])
	for (bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
		size_t imax = min(bi[0]+bstep[0] , size);
		size_t jmax = min(bj[0]+bstep[0] , size);
		for (bk[0] = 0; bk[0] < (bi[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < (bk[0] + bstep[0]); k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
			}
		} // Last block in K, diagonal, divide by pivot
		for (bk[0] = (bi[0]); bk[0] < (bi[0] + bstep[0]); bk[0] += bstep[0]) {
			for (i = bi[0]; i < imax; i += 1)
			for (j = bj[0]; j < jmax; j += 1) {
				for (k = bk[0]; k < i; k += 1)
					ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				if(diagonal == Diagonal::Value)
					ind(X, i, j) /= ind(LU, i, i);
			}
		}
	}
	asm("END LU_Tiled_0");
}


template<Direction direction, Diagonal diagonal, Permute permute,
	class LUMatrix, class XMatrix, class BMatrix>
void substMLU10(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P){
	size_t size = X.size();
	size_t i, j, k;
	size_t bi[5], bj[5], bk[5];
	size_t bimax[5], bjmax[5], bkmax[5];
	size_t bstep[5];
	/**/
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	/* export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} *
	bstep[1] = bstep[0]*L1M;
	bstep[2] = bstep[1]*L2M;
	bstep[2] = bstep[1]*L3M;
	/**/
	//const size_t unr = 2;
	//double acc[unr*unr];
	
	if(permute == Permute::True)
		set<direction>(X,B,P);
	else set<direction>(X,B);
	asm("LU_Tiled_1.0");
	block (bi[1],0,size, bj[1],0,size, bk[1],0,(bi[1]+bstep[1]), bstep[1]) {
		bimax[0] = bi[1] + bstep[1] > size ? size : bi[1] + bstep[1];
		bjmax[0] = bj[1] + bstep[1] > size ? size : bj[1] + bstep[1];
		for (bi[0] = bi[1]; bi[0] < bimax[0]; bi[0] += bstep[0])
		for (bj[0] = bj[1]; bj[0] < bjmax[0]; bj[0] += bstep[0]) {
			bkmax[0] = min(bk[1] + bstep[1], bi[0]);
			size_t imax = min(bi[0]+bstep[0] , size);
			size_t jmax = min(bj[0]+bstep[0] , size);
			for (bk[0] = bk[1]; bk[0] < bkmax[0]; bk[0] += bstep[0]) {
				for (i = bi[0]; i < imax; i += 1)
				for (j = bj[0]; j < jmax; j += 1) {
					for (k = bk[0]; k < (bk[0] + bstep[0]); k += 1)
						ind(X, i, j) = ind(X, i, j) - ind(LU, i, k) * ind(X, k, j);
				}
			} // Last block in K, diagonal
			if(bk[0] == bi[0] && bk[1] == bi[1])
			for (bk[0] = bi[0]; bk[0] < (bi[0] + bstep[0]); bk[0] += bstep[0]) {
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
	asm("END LU_Tiled_1.0");
}
#undef ind

template<class LUMatrix, class IAMatrix, class IMatrix>
void solveMLUNew(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
	//static
	MatrixColMajor<double> Z(X.size());
	//if(Z.size() != X.size()){ Z.alloc(X.size()); }
	// find Z; LZ=B
	substMLU0<Direction::Forwards, Diagonal::Unit, Permute::True>(LU, Z, B, P);
	// find X; Ux=Z
	substMLU0<Direction::Backwards, Diagonal::Value, Permute::False>(LU, X, Z, P);
}