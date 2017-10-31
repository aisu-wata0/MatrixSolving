
namespace gm
{
using namespace std;


void t_vector(size_t size) {
	Matrix<double> LU(size);
	varray<double> X(LU.sizeMem());
	varray<double> B(LU.sizeMem());
	
	size_t bi, bk;
	size_t i, k;
	size_t vj, vk;
	size_t dn = LU.regEN();
	
	cout << size << ":" << LU.sizeMem() << endl;
	
	for(i = 0; i < size; i++){
		X.at(i) = 0;
		B.at(i) = 1;
		for(vj = 0; vj < LU.sizeVec(); vj += 1){
			vec(dj) LU.at(i, vj*dn+dj) = (i);
		}
	}
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
		for(bi = 0; bi < size; bi += bstep){
			unroll(ui,bstep) X.at(bi+ui) = 0;
			for(bk = 0; bk < size; bk += bstep){
				for(i = bi; i < (bi + bstep); i += unr){
					vec<double> acc[unr];
					memset(acc, 0, sizeof(acc));
					for(vk = bk/dn; vk < (bk + bstep)/dn; vk++){
						unroll(ui,unr) acc[ui].v = acc[ui].v + LU.atv(i+ui,vk).v * B.atv(vk).v;
					}
					unroll(ui,unr)
						vec(d) X.at(i+ui) += acc[ui].v[d];
				}
			}
		}
	}
	t[c++] = clock();
	printf("VECmult_Tilroll: %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
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
	}
	t[c++] = clock();
	printf("VECmult_Tiling:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
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
	}
	t[c++] = clock();
	printf("VECmult_Unroll:  %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		for(i = 0; i < size; i++){
			vec<double> acc = {0};
			for(k = 0; k < size/dn; k++){
				acc.v = acc.v + LU.atv(i,k).v * B.atv(k).v;
			}
			X.at(i) = 0;
			vec(d) X.at(i) += acc.v[d];
		}
	}
	t[c++] = clock();
	printf("VECmult:         %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		for(i = 0; i < size; i++){
			double acc[dn] = {0};
			for(k = 0; k < size/dn; k++){
				vec(d) acc[d] = acc[d] + LU.at(i, k*4+d) * B.at(k*4+d);
			}
			X.at(i) = 0;
			vec(d) X.at(i) += acc[d];
		}
	}
	t[c++] = clock();
	printf("Acc:             %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	
	for(reps = 0; reps < repetitionN; reps++){
		for(i = 0; i < size; i++){
			X.at(i) = 0;
			for(k = 0; k < size; k++){
				X.at(i) = X.at(i) + LU.at(i, k) * B.at(k);
			}
		}
	}
	t[c++] = clock();
	printf("Naive:           %f sec\n", ((double)(t[c-1] - t[c-2]) / CLOCKS_PER_SEC) / repetitionN);
	if(PRINT_MATRIX) { printv(X); cout << endl; }
	t[c++] = clock();
	
	cout << X.at(0) <<":"<< X.at(X.size()-1) <<"\tend size "<< size << endl;
}

}