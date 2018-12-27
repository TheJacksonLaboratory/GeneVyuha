#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

int test( Rcpp::List simulationData)
{
  NumericVector simulationParameters = as<NumericVector>(simulationData["simParams"]);
  NumericVector stochasticParameters = as<NumericVector>(simulationData["stochParams"]);
  NumericVector hyperParameters = as<NumericVector>(simulationData["hyperParams"]);
  double g_min = hyperParameters[0];
  double g_max = hyperParameters[1];
  double k_min = hyperParameters[2];
  double k_max = hyperParameters[3];
  double lambda_min = hyperParameters[4];
  double lambda_max = hyperParameters[5];
  double n_min = hyperParameters[6];
  double n_max = hyperParameters[7];
  int possible_interactions = static_cast<int>(hyperParameters[8]);
  int threshold_max = static_cast<int>(hyperParameters[9]);
  double standard_deviation_factor = hyperParameters[10];
  Rcout<<g_min<<"\t"<<g_max<<"\t"<<possible_interactions<<"\t"<<n_max<<"\n";
  possible_interactions = 5;
  hyperParameters[0] = 15.5;
  Rcout<<hyperParameters[8]<<"\t"<<possible_interactions<<"\n";
  Rcout<<hyperParameters[0]<<"\t"<<g_min<<"\n";
  return 0;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
test(simulationData)
*/
