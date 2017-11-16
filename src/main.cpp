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
	LIKWID_MARKER_INIT;
	cout.precision(8);
	cout << scientific;
	srand(20172);
	
	ifstream in_f;
	ofstream o_f;
	streambuf* coutbuf = cout.rdbuf(); //save old buf;
	
	bool input;
	size_t size, iter_n;
	// redirects cout & cin
	parseArgs(argc, argv, input, size, iter_n, in_f, o_f);
	
	Matrix<double> A;
	
	if(input){
		cin>> size;
		A.alloc(size);
		readMatrix(A);
		in_f.close();
	}else {
		A.alloc(size);
		randomMatrix(A);
	}
	
	Matrix<double> LU(size);
	varray<size_t> P(A.sizeMem());
	
	double lu_time = 0.0;
	timer.start();
	//LIKWID_MARKER_START("LU");
	
	GaussEl(A, LU, P);
	
	//LIKWID_MARKER_STOP("LU");
	lu_time = timer.tick();
	
	MatrixColMajor<double> IA(size);
	
	//cout<<"#\n";
	inverse_refining(A, LU, IA, P, iter_n);
	
	cout<< defaultfloat;
	cout<< lu_time <<"\n";
	cout<< timer_iter.averageTotal() <<"\n";
	cout<< timer_residue.averageTotal() << endl;
	//printm(IA);
	
	LIKWID_MARKER_CLOSE;
	in_f.close();
	cout.rdbuf(coutbuf); //redirect
	o_f.close();
	return 0;
}

void parseArgs(int& argc, char**& argv,
bool& input, size_t& size, size_t& iter_n, ifstream& in_f, ofstream& o_f){
	int c;
	input = true;
	size = 0; iter_n = -1;
#define errMsg "Usage: %s [-e inputFile] [-o outputFile] [-r randSize] -i Iterations\n"
	while ((c = getopt(argc, argv, "e:o:r:i:")) != -1){
		switch (c){
			case 'e':
				// inputFile
				in_f.open(optarg);
				cin.rdbuf(in_f.rdbuf()); //redirect
				break;
			case 'o':
				// outputFile
				o_f.open(optarg);
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
				fprintf(stderr, errMsg, argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	
	if(iter_n == -1){
		fprintf(stderr, errMsg, argv[0]);
		fprintf(stderr, "-i Iterations is not optional\n");
		exit(EXIT_FAILURE);
	}
#undef errMsg
}





int mainT(int argc, char **argv)
{
	LIKWID_MARKER_INIT;
	//cout.precision(4);
	//cout << scientific;
	srand(20172);
	
	ifstream in_f;
	ofstream o_f;
	streambuf* coutbuf = cout.rdbuf(); //save old buf;
	
	bool input;
	size_t size, iter_n;
	
	parseArgs(argc, argv, input, size, iter_n, in_f, o_f);
	
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
	vector<size_t> V_sz{144,288,288*2,288*3,288*4};
	//vector<size_t> V_sz = {144};
	for (auto sz : V_sz){
		for(size = sz; size < sz+1; size++){
			//t_matrix_mult(size);
			cout << endl;
		}
	}
	/**/
	LIKWID_MARKER_CLOSE;
	cout.rdbuf(coutbuf); //redirect
	return 0;
}