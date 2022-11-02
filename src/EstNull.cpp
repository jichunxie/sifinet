#include <Rcpp.h>
using namespace Rcpp;

//' EstNull
//' This function is a Rcpp version of Wenguang Sun and Tony T. Cai's 
//' EstNull.func R function, estimating null distribution from data.
//' Wenguang Sun and Tony T. Cai. "Oracle and adaptive compound decision rules for false discovery rate control". 
//' In: Journal of the American Statistical Association 102 (479 2007), pp. 901–912. doi: 10.1198/016214507000000545.
//' @author Qi Gao
//' @param x Input vector of all coexpression values 
//' @param gamma Parameter setting the stopping threshold 
//' @return List of mean and std
//' @export
// [[Rcpp::export]]

List EstNull(NumericVector x, double gamma = 0.1) {
  int n = x.length();
  double gan = std::pow(n, -gamma);
  
  double s;
  double phiplus;
  double phiminus;
  double phi;
  double dphiplus;
  double dphiminus;
  double shat = 0;
  double uhat = 0;
  
  for (int t = 1; t <= 1000; t++){
    s = t * 1.0 / 200;
    phiplus = mean(cos(s * x));
    phiminus = mean(sin(s * x));
    phi = std::sqrt(std::pow(phiplus, 2) + std::pow(phiminus, 2));
    if (phi <= gan){
      dphiplus = -mean(x*sin(s*x));
      dphiminus = mean(x*cos(s*x));
      shat = std::sqrt(-(phiplus * dphiplus + phiminus * dphiminus) / (s * phi * phi));
      uhat = -(dphiplus * phiminus - dphiminus * phiplus) / (phi * phi);
      break;
    }
  }
  
  return List::create(Rcpp::Named("mean") = uhat,
                      Rcpp::Named("std") = shat);
}

