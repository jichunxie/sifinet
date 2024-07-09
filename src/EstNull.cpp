#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' EstNull
//' This function is a Rcpp version of Wenguang Sun and Tony T. Cai's
//' EstNull.func R function, estimating null distribution from data.
//' Sun, W., & Cai, T. T. (2007). Oracle and Adaptive Compound Decision
//' Rules for False Discovery Rate Control.
//' Journal of the American Statistical Association,
//' 102(479), 901–912.
//' @author Qi Gao
//' @param x Input vector of all coexpression values
//' @param gamma Parameter setting the stopping threshold
//' @return List of mean and std
//' @export
// [[Rcpp::export]]

Rcpp::List EstNull(arma::vec x, double gamma = 0.1) {
  int n = x.n_elem;
  double gan = std::pow(n, -gamma);

  double s;
  double phiplus;
  double phiminus;
  double phi;
  double dphiplus;
  double dphiminus;
  double shat = 0;
  double uhat = 0;

  for (int t = 1; t <= 1000; t++) {
    s = t * 1.0 / 200;
    phiplus = mean(cos(s * x));
    phiminus = mean(sin(s * x));
    phi = std::sqrt(std::pow(phiplus, 2) + std::pow(phiminus, 2));
    if (phi <= gan) {
      dphiplus = -arma::as_scalar(x.t() * sin(s * x)) / n;
      dphiminus = arma::as_scalar(x.t() * cos(s * x)) / n;
      shat = std::sqrt(-(phiplus * dphiplus + phiminus * dphiminus) /
                       (s * phi * phi));
      uhat = -(dphiplus * phiminus - dphiminus * phiplus) / (phi * phi);
      break;
    }
  }

  return Rcpp::List::create(Rcpp::Named("mean") = uhat,
                            Rcpp::Named("std") = shat);
}
