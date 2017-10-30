
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

#define asm(x) ;
#include "t_matrix_mult.hpp"
#include "t_vector.hpp"
