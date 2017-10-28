/**
@mainpage

Inverts input matrix using LU decomposition by Gauss Elimination and refining
Usage: %s [-e inputFile] [-o outputFile] [-r randSize] -i Iterations

@authors Bruno Freitas Serbena
@authors Luiz Gustavo Jhon Rodrigues
*/
/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctgmath>
//#include <likwid.h>
#include <unistd.h>

#include "Timer.h"
#include "Matrix.h"
#include "GaussEl.h"
#include "Subst.h"

using namespace std;

Timer timer;
double total_time_iter = 0.0;
double total_time_residue = 0.0;
double lu_time = 0.0;

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
void solve_lu(LUMatrix& LU, XMatrix& X, BMatrix& B, vector<long>& P, long col){
	static vector<double> Z(LU.m_size);
	// find Z; LZ=B
	substR<SubstForwards, DiagonalUnit, SubstPermute>(LU, Z, B, P, col);
	// find X; Ux=Z
	substR<SubstBackwards, DiagonalValue, SubstNoPermute>(LU, X, Z, P, col);
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
void inverse(LUMatrix& LU, IAMatrix& IA, IMatrix& I, vector<long>& P){
	for(long j = 0; j < IA.size; j++){
		// for each IA col solve SL to find the IA col values
		solve_lu(LU, IA, I, P, j);
	}
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
//template<class LUMatrix, class IAMatrix, class IMatrix>
//void inverseB(LUMatrix& LU, IAMatrix& X, IMatrix& B, vector<long>& P){
//	long bi, bj, bk;
//	long i, j, k;
//	long iend; long jend;
//	long step;
//	
//	size = LU.m_size;
//	
//	SubstDirection direction = SubstForwards;
//	if(true){
//		bj = 0;
//		step = +1;
//	} else {
//		bj = size-1;
//		step = -1;
//	}
//	long bstep = step * L1_DN;
//	
//	#define unroll(i,j,unr) for(i = 0; i < unr; i++) for(j = 0; j < unr; j++)
//	
//	long unr = 2;
//	vector<double> acc(unr*unr);
//	long ci, cj;
//	for (bi = 0; bi < size; bi += bstep){
//		for (j = 0; j < size; j += unr){
//			for(i = bi; i != bi + bstep; i += unr){
////				if(direction == SubstForwards)
////					{ endk = min(bk + bstep, i); }
////				else{ endk = max(bk + bstep, i); }
//				// Reset acc
//				unroll(ci,cj,unr) acc.at(ci*unr + cj) = 0;
//				// Iterate to calculate at(X, i, j) .. at(X, i+unr-1, j+unr-1)
//				for(k = 0; k != i; k += step){
//					unroll(ci,cj,unr) acc.at(ci*unr + cj) += LU.at(i+ci, k) * at(X, k, j+cj);
//				}
//				// Remaining in unr
//				unroll(ci,cj,unr) {
//					for(ck = 0; i+ck != i+ci; ck++)
//						acc.at(ci*unr + cj) += LU.at(i+ci, k+ck) * acc.at(ck*unr + cj);
//						// acc.at(ci*unr + cj) += LU.at(i+ci, k+ck) * at(X, i+ck, j+cj);
//					acc.at(ci*unr + cj) /= LU.at(i+ci, i+ci);
//				}
//				// Atribute acc to at(X, i, j) .. at(X, i+unr-1, j+unr-1)
//				unroll(ci,cj,unr) at(X, i+ci, j+cj) = acc.at(ci*unr + cj);
//			}
//		}
//	}
//	///////////////////////////
//			for(i = bi; i != bi + bstep; i += unr){
////				if(direction == SubstForwards)
////					{ endk = min(bk + bstep, i); }
////				else{ endk = max(bk + bstep, i); }
//				// Reset acc
//				unroll(ci,cj,unr) acc.at(ci*unr + cj) = 0;
//				// Iterate to calculate at(X, i, j) .. at(X, i+unr-1, j+unr-1)
//				for(k = 0; k != i; k += step){
//					unroll(ci,cj,unr) acc.at(ci*unr + cj) += LU.at(i+ci, k) * at(X, k, j+cj);
//				}
//				// Remaining in unr
//				unroll(ci,cj,unr) {
//					for(ck = 0; i+ck != i+ci; ck++)
//						acc.at(ci*unr + cj) += LU.at(i+ci, k+ck) * acc.at(ck*unr + cj);
//						// acc.at(ci*unr + cj) += LU.at(i+ci, k+ck) * at(X, i+ck, j+cj);
//					acc.at(ci*unr + cj) /= LU.at(i+ci, i+ci);
//				}
//				// Atribute acc to at(X, i, j) .. at(X, i+unr-1, j+unr-1)
//				unroll(ci,cj,unr) at(X, i+ci, j+cj) = acc.at(ci*unr + cj);
//			}
//			// Remainder diagonal
//			for (i = bi; i != endi; i += step) {
//				if(direction == SubstForwards)
//					{ endk = min(bk + bstep, i); }
//				else{ endk = max(bk + bstep, i); }
//				
//				for (k = bk; k != endk; k += step) {
//					at(X, i,col) = at(X, i,col) - T.at(i,k) * at(X, k,col);
//				}
//				if(diagonal != DiagonalUnit)
//					{ at(X, i,col) /= T.at(i,i); }
//			}
//	///////////////////////////
//	for (; bi >= 0 && bi < size ; bi += bstep) {
//		if(direction == SubstForwards)
//			{ bk = 0;		endi = min(bi + bstep, size); }
//		else{ bk = size-1;	endi = max(bi + bstep, -1); }
//		
//		for(int i = bi; i != endi; i += step)
//			if(Permutation == SubstPermute)
//				{ at(X, i,col) = at(I, P.at(i),col); }
//			else{ at(X, i,col) = at(I, i,col); }
//		
//		for (; bk != bi; bk += bstep) {
//			for (i = bi; i != endi; i += step) {
////				if(direction == SubstForwards)
////					{ endk = min(bk + bstep, i); }
////				else{ endk = max(bk + bstep, i); }
//				endk = bk + bstep;
//				for (k = bk; k != endk; k += step) {
//					at(X, i,col) = at(X, i,col) - T.at(i,k) * at(X, k,col);
//				}
//			}
//		}
//		// Remainder diagonal
//		for (i = bi; i != endi; i += step) {
//			if(direction == SubstForwards)
//				{ endk = min(bk + bstep, i); }
//			else{ endk = max(bk + bstep, i); }
//			
//			for (k = bk; k != endk; k += step) {
//				at(X, i,col) = at(X, i,col) - T.at(i,k) * at(X, k,col);
//			}
//			if(diagonal != DiagonalUnit)
//				{ at(X, i,col) /= T.at(i,i); }
//		}
//	}
//}
/**
 * @brief  Calculates residue into I, A*IA shuold be close to Identity
 * @param A original coef matrix
 * @param IA solution to A*IA = I
 * @param I Output residue, no init needed
 * @return Norm of the residue
 */
template<class AMatrix, class IAMatrix, class IMatrix>
double residue(AMatrix& A, IAMatrix& IA, IMatrix& I){
	double err_norm = 0.0;

	for(long icol=0; icol < A.size; icol++){
		// for each column of the inverse
		for(long i=0; i < A.size; i++){
			// for each line of A
			I.at(i,icol) = 0;
			// multiply A line to the current inverse col
			for(long j=0; j < A.size; j++){
				I.at(i,icol) -= A.at(i,j)*IA.at(j,icol);
			}
			if(i == icol){
				I.at(i,icol) += 1;
			}

			err_norm += I.at(i,icol)*I.at(i,icol);
		}
	}
	return sqrt(err_norm);
}
/**
 * @brief Calculates inverse of A into IA
 * @param LU decomposition of A
 * @param IA return value, no init needed
 * @param P LU pivot permutation
 * @param iter_n
 */
template<class AMatrix, class LUMatrix, class IAMatrix>
void inverse_refining(AMatrix& A, LUMatrix& LU, IAMatrix& IA, vector<long>& P, long iter_n){
	long i=0;
	// number of digits of iter_n, for pretty printing
	long digits = (long)log10((double) iter_n) + 1;
	double c_residue;
	//double l_residue;
	
	// Optm: iterating line by line
	MatrixColMajor W(A.size), R(A.size);
	identity(R);

	//LIKWID_MARKER_START("INV");
	
	inverse(LU, IA, R, P);
	
	//LIKWID_MARKER_STOP("INV");
	//LIKWID_MARKER_START("RES");
	
	c_residue = residue(A, IA, R);
	
	//LIKWID_MARKER_STOP("RES");
	
	cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	while(i < iter_n){
		// (abs(l_residue - c_residue)/c_residue > EPSILON) && (l_residue > c_residue)
		// relative approximate error
		i += 1;
		// R: residue of IA

		timer.start();
		//LIKWID_MARKER_START("INV");
		
		inverse(LU, W, R, P);
		
		//LIKWID_MARKER_STOP("INV");
		// W: residues of each variable of IA
		// adjust IA with found errors
		//LIKWID_MARKER_START("SUM");
		
		add(IA, W);	//TOptm: add at the same time its calculating W
		
		//LIKWID_MARKER_STOP("SUM");
		total_time_iter += timer.elapsed();

		//l_residue = c_residue;
		timer.start();
		//LIKWID_MARKER_START("RES");
		
		c_residue = residue(A, IA, R);
		
		//LIKWID_MARKER_STOP("RES");
		total_time_residue += timer.elapsed();
		
		cout<<"# iter "<< setfill('0') << setw(digits) << i <<": "<< c_residue <<"\n";
	}
}
/**
 * @brief Assigns matrix from cin to M
 * @param A needs to have been allocated
 */
template <class Mat>
void readMatrix(Mat& A){
	for(long i=0; i < A.size; i++){
		for(long j=0; j < A.size; j++){
			cin>> A.at(i,j);
		}
	}
}

int mainBAK(int argc, char **argv) {
	//LIKWID_MARKER_INIT;
	
	cout.precision(17);
	cout << scientific;
	srand(20172);
	
	bool input = true;
	long size = 0, iter_n = -1;

	int c;
	char* inputFile = NULL;
	char* outputFile = NULL;
	ifstream in_f;
	ofstream o_f;
	streambuf* coutbuf = cout.rdbuf(); //save old buf;

	while (( c = getopt(argc, argv, "e:o:r:i:")) != -1){
		switch (c){
			case 'e':
				inputFile = optarg;
				in_f.open(inputFile);
				cin.rdbuf(in_f.rdbuf()); //redirect
				break;
			case 'o':
				outputFile = optarg;
				o_f.open(outputFile);
				cout.rdbuf(o_f.rdbuf()); //redirect
				break;
			case 'r':	//Generate random matrix
				input = false;
				size = stol(optarg);
				break;
			case 'i':
				iter_n = stol(optarg);
				break;
			case ':':
			// missing option argument
				fprintf(stderr, "%s: option '-%c' requires an argument\n", argv[0], optopt);
				break;
			default:
				fprintf(stderr, "Usage: %s [-e inputFile] [-o outputFile] [-r randSize] -i Iterations\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if(iter_n == -1){
		fprintf(stderr, "Usage: %s [-e inputFile] [-o outputFile] [-r randSize] -i Iterations\n", argv[0]);
		fprintf(stderr, "-i Iterations is not optional\n");
		exit(EXIT_FAILURE);
	}
	
	if(input){
		cin>> size;
	}
	
	Matrix A(size);
	
	if(input){
		readMatrix(A);
		in_f.close();
	}else {
		randomMatrix(A);
	}
	
	Matrix LU(size);
	vector<long> P(A.m_size);
	
	timer.start();
	//LIKWID_MARKER_START("LU");
	
	GaussEl(A, LU, P);
	
	//LIKWID_MARKER_STOP("LU");
	lu_time = timer.elapsed();

	// Optm: iterating line by line
	MatrixColMajor IA(size);
	
	cout<<"#\n";
	inverse_refining(A, LU, IA, P, iter_n);

	cout<< defaultfloat;
	cout<<"# Tempo LU: "<< lu_time <<"\n";
	cout<<"# Tempo iter: "<< total_time_iter/(double)iter_n <<"\n";
	cout<<"# Tempo residuo: "<< total_time_residue/(double)iter_n <<"\n#\n";
	printm(IA);
	
	//LIKWID_MARKER_CLOSE;
	in_f.close();
	cout.rdbuf(coutbuf); //redirect
	o_f.close();
	return 0;
}







#define asm(x) ;

#include "t_matrix_mult.h"

#include "t_vector.h"

int main()
{
	//LIKWID_MARKER_INIT;
	//cout.precision(4);
	//cout << scientific;
	srand(20172);
	size_t size = 8192*2;
	/**
	vector<size_t> V_sz = {8192/4,8192/2};
	//vector<size_t> V_sz = {BL1*2};
	for (auto& sz : V_sz){
		sz = Lower_Multiple(sz,BL1);
	}
	for (auto sz : V_sz){
		for(size = sz; size < sz+1; size++){
			t_vector(size);
			cout << endl;
		}
	}
	/**/
	//vector<size_t> V_sz = {144,288,288*2,288*3,288*4};
	vector<size_t> V_sz = {8};
	for (auto sz : V_sz){
		for(size = sz; size < sz+1; size++){
			t_matrix_mult(size);
			cout << endl;
		}
	}
	/**/
	//LIKWID_MARKER_CLOSE;
	return 0;
}