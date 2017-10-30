/**
@mainpage

Inverts input matrix using LU decomposition by Gauss Elimination and refining
Usage: %s [-e inputFile] [-o outputFile] [-r randSize] -i Iterations

@authors Bruno Freitas Serbena
@authors Luiz Gustavo Jhon Rodrigues
*/

#include "main.hpp"

/**
 * @brief Assigns matrix from cin to M
 * @param A needs to have been allocated
 */
template <class Mat>
void readMatrix(Mat& A){
	for(long i=0; i < A.size(); i++){
		for(long j=0; j < A.size(); j++){
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
	
	Matrix<double> A(size);
	
	if(input){
		readMatrix(A);
		in_f.close();
	}else {
		randomMatrix(A);
	}
	
	Matrix<double> LU(size);
	vector<long> P(A.sizeMem());
	
	timer.start();
	//LIKWID_MARKER_START("LU");
	
	GaussEl(A, LU, P);
	
	//LIKWID_MARKER_STOP("LU");
	lu_time = timer.tick();
	
	
	MatrixColMajor<double> IA(size);
	/**
	cout<<"#\n";
	inverse_refining(A, LU, IA, P, iter_n);

	cout<< defaultfloat;
	cout<<"# Tempo LU: "<< lu_time <<"\n";
	cout<<"# Tempo iter: "<< total_time_iter/(double)iter_n <<"\n";
	cout<<"# Tempo residuo: "<< total_time_residue/(double)iter_n <<"\n#\n";
	printm(IA);
	/**/
	
	MatrixColMajor<double> I(size);
	identity(I);
	
	//inverseLUT(LU, IA, I, P);
	inverseLU(LU, IA, I, P);
	printm(IA);
	
	
	//LIKWID_MARKER_CLOSE;
	in_f.close();
	cout.rdbuf(coutbuf); //redirect
	o_f.close();
	return 0;
}









int mainTests()
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
	vector<size_t> V_sz = {144,288,288*2,288*3,288*4};
	//vector<size_t> V_sz = {8};
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