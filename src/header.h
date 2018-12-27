#ifndef HEADER_H
#define HEADER_H

#include <Rcpp.h>
#include <iostream>
#include <random>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
// Shifted hill function
extern double Hs_Racipe(double A, double AB0, size_t n_ab, double lambda_ab);
// change hyperparameter values according to parameter range
extern void applyParameterRange(const double& parameter_range, double& g_min,
                                double& g_max, double& k_min,
                         double& k_max, double lambda_min, double& lambda_max);
//uniformly distributed random number generator in (0,1) range
extern std::mt19937_64 u_generator;
extern std::uniform_real_distribution<double> u_distribution;
extern std::mt19937_64 g_generator;
extern std::normal_distribution<double> g_distribution;
// All definitions in header.cpp
#endif
