
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

using namespace std;
using namespace gm;

#define ChronometerHistoryMax 16
Chronometer<ChronometerHistoryMax> timer;
double total_time_iter = 0.0;
double total_time_residue = 0.0;
double lu_time = 0.0;
