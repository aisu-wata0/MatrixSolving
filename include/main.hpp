
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctgmath>
//#include <likwid.h>
#include <unistd.h>

#include "Matrix.hpp"
#include "GaussEl.hpp"
#include "Subst.hpp"
#include "Chronometer.hpp"
#include "SolveLU.hpp"

using namespace std;
using namespace gm;

void parseArgs(int& argc, char**& argv, bool& input, size_t& size, size_t& iter_n, ifstream& in_f, ofstream& o_f);
