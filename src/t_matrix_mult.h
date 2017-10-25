
//inline void test(Matrix& LU, MatrixColMajor& B, Matrix& X){
void test(long size){
	Matrix LU(size);
	MatrixColMajor B(size);
	Matrix X(size);
	
	long bi, bj, bk;
	long i, j, k;
	long iend; long jend;
	long step;
	
	size = LU.m_size;
	cout << size <<endl;
	
	asm("SETUP");
	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			LU.at(i,j) = i*size + j;
			B.at(i,j) = j*size + i;
		}
	}
	
	SubstDirection direction = SubstForwards;
	if(true){
		bj = 0;
		step = +1;
	} else {
		bj = size-1;
		step = -1;
	}
	long bstep = step * BL1; // cache block size
	
	set(X,0);
	timer.start();
	asm("NO BLOCKING");
	for(i=0; i < size; i += step){
		for(j=0; j < size; j += step){
			X.at(i,j) = 0;
			for(k=0; k < size; k += 1){	
				at(X, i,j) = at(X, i,j) + LU.at(i,k) * at(B, k,j);
			}
		}
	}
	asm("END NO BLOCKING");
	cout<<"# "<< timer.elapsed() <<"\n";
	// printm(X);
	
	double acc00, acc01, acc10, acc11;
	/**
	set(X,0);
	timer.start();
	asm("BLOCKING BI");
	for (bi = 0; bi < size; bi += bstep){
		for (j = 0; j < size; j += 2){
			for (i = bi; i < bi + bstep; i += 2){
				acc00 = acc01 = acc10 = acc11 = 0;
				for (k = 0; k < size; k++){
					acc00 += LU.at(i+0, k) * at(B, k, j+0);
					acc01 += LU.at(i+0, k) * at(B, k, j+1);
					acc10 += LU.at(i+1, k) * at(B, k, j+0);
					acc11 += LU.at(i+1, k) * at(B, k, j+1);
				}
				at(X, i+0,j +0) = acc00;
				at(X, i+0,j +1) = acc01;
				at(X, i+1,j +0) = acc10;
				at(X, i+1,j +1) = acc11;
			}
		}
	}
	asm("END BLOCKING BI");
	cout<<"# "<< timer.elapsed() <<"\n";
	// printm(X);
	/**/
	
	#define unr 2
	#define unroll(i,j,unr) for(long i = 0; i < unr; i++) for(long j = 0; j < unr; j++)
	
	double acc[unr*unr];
	long ci, cj;
	set(X,0);
	timer.start();
	asm("BLOCKING BI");
	for (bi = 0; bi < size; bi += bstep){
		for (j = 0; j < size; j += unr){
			for (i = bi; i < (bi + bstep); i += unr){
				unroll(ci,cj,unr) acc[ci*unr + cj] = 0;
				for (k = 0; k < size; k++){
					unroll(ci,cj,unr) acc[ci*unr + cj] += LU.at(i+ci, k) * at(B, k, j+cj);
				}
				unroll(ci,cj,unr) at(X, i+ci, j+ci) = acc[ci*unr + cj];
			}
		}
	}
	asm("END BLOCKING BI");
	cout<<"# "<< timer.elapsed() <<"\n";
	// printm(X);
	
	
	set(X,0);
	timer.start();
	asm("BLOCKING BK");
	for(bi = 0; bi < size; bi += bstep){
		for(bk = 0; bk < size; bk += bstep){
			for(j=0; j < size; j += 2){
				for(i = bi; i < bi + bstep; i += 2 ){
					if(bk == 0){
						acc00 = acc01 = acc10 = acc11 = 0;
					} else {
						acc00 = at(X, i+0,j +0);
						acc01 = at(X, i+0,j +1);
						acc10 = at(X, i+1,j +0);
						acc11 = at(X, i+1,j +1);
					}
					for(k = bk; k < bk + bstep; k++){
						acc00 += LU.at(i+0, k) * at(B, k, j+0);
						acc01 += LU.at(i+0, k) * at(B, k, j+1);
						acc10 += LU.at(i+1, k) * at(B, k, j+0);
						acc11 += LU.at(i+1, k) * at(B, k, j+1);
					}
					at(X, i+0,j +0) = acc00;
					at(X, i+0,j +1) = acc01;
					at(X, i+1,j +0) = acc10;
					at(X, i+1,j +1) = acc11;
				}
			}
		}
	}
	asm("END BLOCKING BK");
	cout<<"# "<< timer.elapsed() <<"\n";
	// printm(X);
}