
void test2(long size){
	Matrix LU(size);
	MatrixColMajor B(size);
	MatrixColMajor X(size);
	
	long bi, bj, bk;
	long i, j, k;
	long iend; long jend;
	long step;
    #define vec(v) for(size_t v = 0; v < dn; v++)
	size = LU.m_size;
	
	cout << size << endl;
	
	asm("VECTOR");
	for(i = 0; i < size; i++){
		for(j = 0; j != size/dn; j += 1){
			vec(vj) LU.at(i,j*dn+vj) = (i*size + j*dn+vj);
			X.atv(j,i).v = LU.atv(i,j).v;
		}
	}
	asm("END VECTOR");
	
	printm(LU);
	printm(X);
	/**/
	asm("VMULT");
	for(i = 0; i < size; i++){
		vdouble acc = {{0}};
		for(k = 0; k < size/dn; k++){
			acc.v = acc.v + LU.atv(i,k).v * X.atv(k,0).v;
		}
		B.at(i,0) = 0;
		vec(v) B.at(i,0) += acc.d[v];
	}
	asm("END VMULT");
	/**
	asm("MULT");
	for(i = 0; i < size; i++){
		double acc[dn] = {0};
		for(k = 0; k < size/dn; k++){
			vec(v) acc[v] = acc[v] + LU.at(i, k*4+v) * X.at(k*4+v, 0);
		}
		B.at(i,0) = 0;
		vec(v) B.at(i,0) += acc[v];
	}
	asm("END MULT");
	/**
	asm("M");
	for(i = 0; i < size; i++){
		B.at(i,0) = 0;
		for(k = 0; k < size; k++){
			B.at(i,0) = B.at(i,0) + LU.at(i, k) * X.at(k, 0);
		}
	}
	asm("END M");
	/**/
	printm(B);
}
