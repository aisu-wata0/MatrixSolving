
void t_vector(long size){
	Matrix LU(size);
	varray<double> X(LU.m_size);
	varray<double> B(LU.m_size);
	
	long bi, bk;
	long i, j, k;
	long vj, vk;
	size = LU.m_size;
	
	cout << size << endl;
	
	asm("Setup");
	for(i = 0; i < size; i++){
		X.at(i) = 0;
		B.at(i) = i;
		for(vj = 0; vj < size/dn; vj += 1){
			vec(dj) LU.at(i, vj*dn+dj) = (i*size + vj*dn+dj);
		}
	}
	asm("END Setup");
	bool PRINT_MATRIX = false;
	size_t repetionN = 100;
	size_t reps;
	
	size_t unr = 4;
	int c = 0;
	clock_t t[128];
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("VECmult_Tilroll");
		for(bi = 0; bi < size; bi += BL1){
			unroll(ui,BL1) X.at(bi+ui) = 0;
			for(bk = 0; bk < size; bk += BL1){
				//asm("VECmult_Tilroll inner");
				for(i = bi; i < (bi + BL1); i += unr){
					vec<double> acc[unr];
					//memset(acc, 0, sizeof(acc));
					unroll(ui,unr) vec(d) acc[ui].v[d] = 0;
					for(vk = bk/dn; vk < (bk + BL1)/dn; vk++){
						unroll(ui,unr) acc[ui].v = acc[ui].v + LU.atv(i+ui,vk).v * B.atv(vk).v;
					}
					unroll(ui,unr)
						vec(d) X.at(i+ui) += acc[ui].v[d];
				}
				//asm("END VECmult_Tilroll inner");
			}
		}
		asm("END VECmult_Tilroll");
	}
	t[c++] = clock();
	printf("VECmult_Tilroll: %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("VECmult_Tiling");
		for(bi = 0; bi < size; bi += BL1){
			unroll(ui,BL1) X.at(bi+ui) = 0;
			for(bk = 0; bk < size; bk += BL1){
				for(i = bi; i < (bi + BL1); i++){
					vec<double> acc = {0};
					for(vk = bk/dn; vk < (bk + BL1)/dn; vk++){
						acc.v = acc.v + LU.atv(i,vk).v * B.atv(vk).v;
					}
					vec(d) X.at(i) += acc.v[d];
				}
			}
		}
		asm("END VECmult_Tiling");
	}
	t[c++] = clock();
	printf("VECmult_Tiling:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("VECmult_Unroll");
		for(i = 0; i < size; i += unr){
			vec<double> acc[unr];
			unroll(ui,unr)
				vec(d) acc[ui].v[d] = 0;
			for(k = 0; k < size/dn; k++){
				unroll(ui,unr) acc[ui].v = acc[ui].v + LU.atv(i+ui,k).v * B.atv(k).v;
			}
			unroll(ui,unr) X.at(i+ui) = 0;
			unroll(ui,unr)
				vec(d) X.at(i+ui) += acc[ui].v[d];
		}
		asm("END VECmult_Unroll");
	}
	t[c++] = clock();
	printf("VECmult_Unroll:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("VECmult");
		for(i = 0; i < size; i++){
			vec<double> acc = {0};
			for(k = 0; k < size/dn; k++){
				acc.v = acc.v + LU.atv(i,k).v * B.atv(k).v;
			}
			X.at(i) = 0;
			vec(d) X.at(i) += acc.v[d];
		}
		asm("END VECmult");
	}
	t[c++] = clock();
	printf("VECmult:         %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("Acc");
		for(i = 0; i < size; i++){
			double acc[dn] = {0};
			for(k = 0; k < size/dn; k++){
				vec(d) acc[d] = acc[d] + LU.at(i, k*4+d) * B.at(k*4+d);
			}
			X.at(i) = 0;
			vec(d) X.at(i) += acc[d];
		}
		asm("END Acc");
	}
	t[c++] = clock();
	printf("Acc:             %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetionN; reps++){
		asm("Naive");
		for(i = 0; i < size; i++){
			X.at(i) = 0;
			for(k = 0; k < size; k++){
				X.at(i) = X.at(i) + LU.at(i, k) * B.at(k);
			}
		}
		asm("END Naive");
	}
	t[c++] = clock();
	printf("Naive:           %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	cout << X.at(0) <<":"<< X.at(X.size-1) <<"\tend size "<< size << endl;
}
