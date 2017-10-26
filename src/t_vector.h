
void t_vector(long size){
	Matrix LU(size);
	vector<double> X(LU.m_size);
	MatrixColMajor B(size);
	
	long bi, bk;
	long i, j, k;
	size = LU.m_size;
	
	cout << size << endl;
	
	asm("Setup");
	for(i = 0; i < size; i++){
		X.at(i) = 0;
		for(j = 0; j != size/dn; j += 1){
			vec(dj) LU.at(i,j*dn+dj) = (i*size + j*dn+dj);
			B.atv(j,i) = LU.atv(i,j);
		}
	}
	asm("END Setup");
	bool PRINT_MATRIX = false;
	
	size_t unr = 4;
	int c = 0;
	clock_t t[128];
	t[c++] = clock();
	asm("VECmult_Tilroll");
	for(bi = 0; bi < size; bi += BL1){
		unroll(ui,BL1) X.at(bi+ui) = 0;
		for(bk = 0; bk < size/dn; bk += BL1/dn){
			asm("VECmult_Tilroll inner");
			for(i = bi; i < (bi + BL1); i += unr){
				vdouble acc[unr];
				//vdouble acc[unr];
				//memset(acc, 0, sizeof(acc));
				unroll(ui,unr) vec(d) acc[ui][d] = 0;
				for(k = bk; k < (bk + BL1/dn); k++){
					unroll(ui,unr) acc[ui] = acc[ui] + LU.atv(i+ui,k) * B.atv(k,0);
				}
				unroll(ui,unr)
					vec(d) X.at(i+ui) += acc[ui][d];
			}
			asm("END VECmult_Tilroll inner");
		}
	}
	asm("END VECmult_Tilroll");
	t[c++] = clock();
	printf("VECmult_Tilroll: %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	asm("VECmult_Tiling");
	for(bi = 0; bi < size; bi += BL1){
		unroll(ui,BL1) X.at(bi+ui) = 0;
		for(bk = 0; bk < size/dn; bk += BL1/dn){
			for(i = bi; i < (bi + BL1); i++){
				vdouble acc = {0};
				for(k = bk; k < (bk + BL1/dn); k++){
					acc = acc + LU.atv(i,k) * B.atv(k,0);
				}
				vec(d) X.at(i) += acc[d];
			}
		}
	}
	asm("END VECmult_Tiling");
	t[c++] = clock();
	printf("VECmult_Tiling:  %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	asm("VECmult_Unroll");
	for(i = 0; i < size; i += unr){
		vdouble acc[unr];
		unroll(ui,unr)
			vec(d) acc[ui][d] = 0;
		for(k = 0; k < size/dn; k++){
			unroll(ui,unr) acc[ui] = acc[ui] + LU.atv(i+ui,k) * B.atv(k,0);
		}
		unroll(ui,unr) X.at(i+ui) = 0;
		unroll(ui,unr)
			vec(d) X.at(i+ui) += acc[ui][d];
	}
	asm("END VECmult_Unroll");
	t[c++] = clock();
	printf("VECmult_Unroll:  %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	asm("VECmult");
	for(i = 0; i < size; i++){
		vdouble acc = {0};
		for(k = 0; k < size/dn; k++){
			acc = acc + LU.atv(i,k) * B.atv(k,0);
		}
		X.at(i) = 0;
		vec(d) X.at(i) += acc[d];
	}
	asm("END VECmult");
	t[c++] = clock();
	printf("VECmult:         %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	asm("Acc");
	for(i = 0; i < size; i++){
		double acc[dn] = {0};
		for(k = 0; k < size/dn; k++){
			vec(d) acc[d] = acc[d] + LU.at(i, k*4+d) * B.at(k*4+d, 0);
		}
		X.at(i) = 0;
		vec(d) X.at(i) += acc[d];
	}
	asm("END Acc");
	t[c++] = clock();
	printf("Acc:             %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	asm("Naive");
	for(i = 0; i < size; i++){
		X.at(i) = 0;
		for(k = 0; k < size; k++){
			X.at(i) = X.at(i) + LU.at(i, k) * B.at(k, 0);
		}
	}
	asm("END Naive");
	t[c++] = clock();
	printf("Naive:           %f sec\n", (double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	cout << X.at(0) <<":"<< X.at(X.size()-1) <<"\tend size "<< size << endl;
}
