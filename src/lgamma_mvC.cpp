#include <Rcpp.h>
using namespace Rcpp;

// C++ implementation of the logarithm of the gamma function
double lgamma_mvC(double x, double p){

  double res = (p*(p-1.0)/4.0)*log(M_PI);

  for(int i=0; i < p; i++){
    res = res + lgamma(x - (double)i/2.0);
  }

  return (res);
}