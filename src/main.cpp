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
void solve_lu(Matrix& LU, MatrixColMajor& X, MatrixColMajor& B, vector<long>& P, long col){
	vector<double> Z(LU.size);
	// find Z; LZ=B
	subst_P<SubstForwards, DiagonalUnit>(LU, Z, B, P, col);
	// find X; Ux=Z
	subst<SubstBackwards>(LU, X, Z, col);
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
void inverse(Matrix& LU, MatrixColMajor& IA, MatrixColMajor& I, vector<long>& P){
	for(long j = 0; j < IA.size; j++){
		// for each IA col solve SL to find the IA col values
		solve_lu(LU, IA, I, P, j);
	}
}
/**
 * @brief  Calculates residue into I, A*IA shuold be close to Identity
 * @param A original coef matrix
 * @param IA solution to A*IA = I
 * @param I Output residue, no init needed
 * @return Norm of the residue
 */
double residue(Matrix& A, MatrixColMajor& IA, MatrixColMajor& I){
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
void inverse_refining(Matrix& A, Matrix& LU, MatrixColMajor& IA, vector<long>& P, long iter_n){
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

int main(int argc, char **argv) {
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
	vector<long> P(size);
	
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


int mainTEST()
{
	SubstDirection direction = SubstForwards;
	long size = 3;
	Matrix A(size);
	set(A, 0);
	// TODO test iterate blocks by col (bj) instead of as currently by row (bi)
	long bi; long bj;
	long i; long j;
	long iend; long jend;
	long step;
	// TODO test iterate blocks by col (bj) instead of as currently by row (bi)
	// forwards
	if(direction == SubstForwards){
		bj = 0;
		step = +1;
	} else {
		bj = size-1;
		step = -1;
	}
	long bstep = step * CACHE_LSZ;
	
	for(; bj >= 0 && bj < size; bj += bstep){
		if(direction == SubstForwards) {
			jend = bj + bstep > size ? size : bj + bstep;
		} else {
			jend = bj + bstep < 0 ? -1 : bj + bstep;
		}
		// go though diagonal block
		for(i = bj; i != jend; i += step){
			if((bj == 0 && direction == SubstForwards) || (bj == size-1 && direction == SubstBackwards)){
				cout << "started line (" << i << "," << bj << ") dia" << endl;
				if(A.at(i,i) != 0) { cout << "ERROR: already started" << endl; }
				A.at(i,i) = 1;
			}
			for(j = bj; j != jend && j != i; j += step){
				cout << "(" << i << "," << j << ") dia" << endl; 
				if(A.at(j,j) != -1 || A.at(i,i) < 1) cout << "ERROR: X not found yet or not inited" << endl; 
			}
			cout << "divided by (" << i << "," << i << ") dia" << endl;
			if(A.at(i,i) < 1) cout << "ERROR: dividend not inited" ;
			A.at(i,i) = -1;
		}
		for(bi = bj+bstep; bi >= 0 && bi < size; bi += bstep){
			if(direction == SubstForwards) {
				iend = bi + bstep > size ? size : bi + bstep;
			} else {
				iend = bi + bstep < 0 ? -1 : bi + bstep;
			}
			if((bj == 0 && direction == SubstForwards) || (bj == size-1 && direction == SubstBackwards))
				for(i = bi; i != iend; i += step){
					cout << "started line (" << i << "," << bj << ")" << endl;
					if(A.at(i,i) != 0) { cout << "ERROR: already started" << endl; }
					A.at(i,i) = 1;
				}
			// go though current block
			for(i = bi; i != iend; i += step){
				for(j = bj; j != jend && j != i; j += step){
					cout << "(" << i << "," << j << ")" << endl; 
				}
			}
		}
		// remainder
//		for(i = size - mod(size, bstep); i < size; i += step){
//			for(j = bj; j < bj+bstep; j += step){
//				cout << "(" << i << "," << j << ")" << endl; 
//			}
//		}
	}
	// remainder
//	for(i = size - mod(size, bstep); i < size; i += step){
//		for(j = size - mod(size, bstep); j != i; j += step){
//			cout << "(" << i << "," << j << ")" << endl; 
//		}
//		cout << "divided by (" << i << "," << i << ") dia" << endl;
//	}
	return 0;
}