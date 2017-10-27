

#define block(i,i0,imax, j,j0,jmax, k,k0,kmax, bstep) \
for (size_t i = i0; i < imax; i += bstep) \
	for (size_t j = j0; j < jmax; j += bstep) \
		for (size_t k = k0; k < kmax; k += bstep)

//inline void test(Matrix& LU, MatrixColMajor& B, Matrix& X){
void t_matrix_mult(size_t size){
	Matrix LU(size);
	MatrixColMajor X(size);
	MatrixColMajor B(size);
	
	size_t i, j, k;
	size_t step;
	
	size = LU.size;
	const size_t b_size = LU.b_size;
	cout << size << ":" << b_size << ":" << LU.m_size << endl;
	cout <<"L1_LINE_DN = "<< L1_LINE_DN <<";  L1_DN = "<< L1_DN <<";  BL1 = "<< BL1 <<";  B3L1 = "<< B3L1 <<";  "<< endl;
	
	asm("SETUP");
	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			if(true){
				LU.at(i,j) = i*size + j;
				B.at(j,i) = j*size + i;
			} else {
				LU.at(i,j) = j;
				B.at(j,i) = (double)i/size;
			}
		}
	}
	
	step = +1;
	//bstep = step * BL1/3; // cache block size
	size_t bi[5], bj[5], bk[5];
	size_t bstep[5];
	bstep[0] = 8;
	bstep[1] = bstep[0]*3;
	bstep[2] = bstep[1]*3;
	bstep[3] = bstep[2]*4;
	// export GCC_ARGS=" -D L1M=${3} -D L2M=${3} L3M=${4} "
	//bstep[1] = bstep[0]*L1M;
	//bstep[2] = bstep[1]*L2M;
	//bstep[2] = bstep[1]*L3M;
	
	const size_t unr = 2;
	const bool PRINT_MATRIX = false;
	double acc00, acc01, acc10, acc11;
	double acc[unr*unr];
	
	Chronometer<128> timer;
	double cl;
	
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
	
	
	set(X,0);
	timer.tick();
	asm("Tiled0");
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
	asm("END Tiled0");
	printf("Tiled0: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**
	
	
	set(X,0);
	timer.tick();
	asm("Tiled0");
	for (size_t bi[0] = 0; bi[0] < size; bi[0] += bstep[0]) {
		for (size_t bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
			for (size_t bk[0] = 0; bk[0] < size; bk[0] += bstep[0]) {
				size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
				size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
				size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];

				for (size_t j = bj[0]; j < jmax; ++j) {
					for (size_t i = bi[0]; i < imax; ++i) {
						for (size_t k = bk[0]; k < kmax; ++k) {
							X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
						}
					}
				}
			}
		}
	}
	asm("END Tiled0");
	printf("Tiled0: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
	
	
	set(X,0);
	timer.tick();
	asm("Tiled1");
	block(i2,0,size, j2,0,size, k2,0,size, bstep[1]){
		size_t i1max = i2 + bstep[1] > size ? size : i2 + bstep[1];
		size_t j1max = j2 + bstep[1] > size ? size : j2 + bstep[1];
		size_t k1max = k2 + bstep[1] > size ? size : k2 + bstep[1];
		block(i1,i2,i1max, j1,j2,j1max, k1,k2,k1max, bstep[0]){
			size_t imax = i1 + bstep[0] > size ? size : i1 + bstep[0];
			size_t jmax = j1 + bstep[0] > size ? size : j1 + bstep[0];
			size_t kmax = k1 + bstep[0] > size ? size : k1 + bstep[0];
			block(i,i1,imax, j,j1,jmax, k,k1,kmax, 1){
				X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
			}
		}
	}
	asm("END Tiled1");
	printf("Tiled1: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	set(X,0);
	timer.tick();
	asm("Tiled2");
	block(i3,0,size, j3,0,size, k3,0,size, bstep[2]){
		size_t i2max = i3 + bstep[2] > size ? size : i3 + bstep[2];
		size_t j2max = j3 + bstep[2] > size ? size : j3 + bstep[2];
		size_t k2max = k3 + bstep[2] > size ? size : k3 + bstep[2];
		block(i2,i3,i2max, j2,j3,j2max, k2,k3,k2max, bstep[1]){
			size_t i1max = i2 + bstep[1] > size ? size : i2 + bstep[1];
			size_t j1max = j2 + bstep[1] > size ? size : j2 + bstep[1];
			size_t k1max = k2 + bstep[1] > size ? size : k2 + bstep[1];
			block(i1,i2,i1max, j1,j2,j1max, k1,k2,k1max, bstep[0]){
				size_t imax = i1 + bstep[0] > size ? size : i1 + bstep[0];
				size_t jmax = j1 + bstep[0] > size ? size : j1 + bstep[0];
				size_t kmax = k1 + bstep[0] > size ? size : k1 + bstep[0];
				block(i,i1,imax, j,j1,jmax, k,k1,kmax, 1){
					X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
				}
			}
		}
	}
	asm("END Tiled2");
	printf("Tiled2: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	
	set(X,0);
	timer.tick();
	asm("Tiled3");
	block(i4,0,size, j4,0,size, k4,0,size, bstep[3]){
		size_t i3max = i4 + bstep[3] > size ? size : i4 + bstep[3];
		size_t j3max = j4 + bstep[3] > size ? size : j4 + bstep[3];
		size_t k3max = k4 + bstep[3] > size ? size : k4 + bstep[3];
		block(i3,i4,i3max, j3,j4,j3max, k3,k4,k3max, bstep[2]){
			size_t i2max = i3 + bstep[2] > size ? size : i3 + bstep[2];
			size_t j2max = j3 + bstep[2] > size ? size : j3 + bstep[2];
			size_t k2max = k3 + bstep[2] > size ? size : k3 + bstep[2];
			block(i2,i3,i2max, j2,j3,j2max, k2,k3,k2max, bstep[1]){
				size_t i1max = i2 + bstep[1] > size ? size : i2 + bstep[1];
				size_t j1max = j2 + bstep[1] > size ? size : j2 + bstep[1];
				size_t k1max = k2 + bstep[1] > size ? size : k2 + bstep[1];
				block(i1,i2,i1max, j1,j2,j1max, k1,k2,k1max, bstep[0]){
					size_t imax = i1 + bstep[0] > size ? size : i1 + bstep[0];
					size_t jmax = j1 + bstep[0] > size ? size : j1 + bstep[0];
					size_t kmax = k1 + bstep[0] > size ? size : k1 + bstep[0];
					block(i,i1,imax, j,j1,jmax, k,k1,kmax, 1){
						X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
					}
				}
			}
		}
	}
	asm("END Tiled3");
	printf("Tiled3: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	size_t bkmax[5];
	set(X,0);
	timer.tick();
	asm("Tiled0");
	for (size_t bi[0] = 0; bi[0] < size; bi[0] += bstep[0]) {
		for (size_t bj[0] = 0; bj[0] < size; bj[0] += bstep[0]) {
			bkmax[0] = bi[0];
			bkmax[0] = bi[0] + bstep[0];
			for (size_t bk[0] = 0; bk[0] < kmax[0]; bk[0] += bstep[0]) {
				size_t imax = bi[0] + bstep[0] > size ? size : bi[0] + bstep[0];
				size_t jmax = bj[0] + bstep[0] > size ? size : bj[0] + bstep[0];
				size_t kmax = bk[0] + bstep[0] > size ? size : bk[0] + bstep[0];

				for (size_t j = bj[0]; j < jmax; ++j) {
					for (size_t i = bi[0]; i < imax; ++i) {
						for (size_t k = bk[0]; k < kmax; ++k) {
							X.at(i, j) = X.at(i, j) + LU.at(i, k) * B.at(k, j);
						}
					}
				}
			}
		}
	}
	asm("END Tiled0");
	printf("Tiled0: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
	
	set(X,0);
	timer.tick();
	asm("Naive");
	for(i=0; i < size; i += step){
		for(j=0; j < size; j += step){
			for(k=0; k < size; k += 1){	
				X.at(i,j) = X.at(i,j) + LU.at(i,k) * B.at(k,j);
			}
		}
	}
	asm("END Naive");
	printf("Naive: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	
	
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
	printf("Tile_BI: %f sec\n", timer.tick(&cl));
	if(PRINT_MATRIX) { printm(X); cout << endl; }
	/**/
}