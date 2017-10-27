
void t_vector(size_t size) {
	Matrix LU(size);
	varray<double> X(LU.m_size);
	varray<double> B(LU.m_size);
	
	size_t bi, bk;
	size_t i, j, k;
	size_t vj, vk;
	
	cout << size << ":" << LU.m_size << endl;
	
	asm("Setup");
	for(i = 0; i < size; i++){
		X.at(i) = 0;
		B.at(i) = 1;
		for(vj = 0; vj < size/dn; vj += 1){
			vec(dj) LU.at(i, vj*dn+dj) = (i);
		}
	}
	asm("END Setup");
	bool PRINT_MATRIX = false;
	if(PRINT_MATRIX) { print(LU); cout << endl; }
	if(PRINT_MATRIX) { printv(B); cout << endl; }
	
	size_t repetitionN = 100;
	size_t reps;
	
	const size_t bstep = BL1;
	size_t unr = 4;
	int c = 0;
	clock_t t[128];
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		asm("VECmult_Tilroll");
		for(bi = 0; bi < size; bi += bstep){
			unroll(ui,bstep) X.at(bi+ui) = 0;
			for(bk = 0; bk < size; bk += bstep){
				asm("VECmult_Tilroll_inner");
				for(i = bi; i < (bi + bstep); i += unr){
					vec<double> acc[unr];
					//memset(acc, 0, sizeof(acc));
					unroll(ui,unr) vec(d) acc[ui].v[d] = 0;
					for(vk = bk/dn; vk < (bk + bstep)/dn; vk++){
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
	printf("VECmult_Tilroll: %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		asm("VECmult_Tiling");
		for(bi = 0; bi < size; bi += bstep){
			unroll(ui,bstep) X.at(bi+ui) = 0;
			for(bk = 0; bk < size; bk += bstep){
				for(i = bi; i < (bi + bstep); i++){
					vec<double> acc = {0};
					for(vk = bk/dn; vk < (bk + bstep)/dn; vk++){
						acc.v = acc.v + LU.atv(i,vk).v * B.atv(vk).v;
					}
					vec(d) X.at(i) += acc.v[d];
				}
			}
		}
		asm("END VECmult_Tiling");
	}
	t[c++] = clock();
	printf("VECmult_Tiling:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		asm("VECmult_Unroll");
		for(i = 0; i < size; i += unr){
			vec<double> acc[unr];
			unroll(ui,unr) vec(d) acc[ui].v[d] = 0;
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
	printf("VECmult_Unroll:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
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
	printf("VECmult:         %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
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
	printf("Acc:             %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
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
	printf("Naive:           %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	cout << X.at(0) <<":"<< X.at(X.size-1) <<"\tend size "<< size << endl;
}
