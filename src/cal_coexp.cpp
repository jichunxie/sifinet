#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' cal_coexp
//' This function calculates the coexpression patterns between genes
//' and returns the coexpression matrix.
//' @author Qi Gao
//' @param X Input binarized cell (row) by gene (column) matrix
//' @return Coexpression matrix
//' @export
// [[Rcpp::export]]
arma::mat cal_coexp(arma::mat X) {
  int p = X.n_cols;
  int n = X.n_rows;
  arma::vec q(p);
  for (int i = 0; i < p; i++) {
    q(i) = mean(X.col(i));
  }
  arma::vec mq = 1 - q;
  arma::mat c = X.t() * X - q * q.t() * n;
  arma::mat d = sqrt(n * q * q.t() % (mq * mq.t()));

  return (c / d);
}
