#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//' cal_conn
//' This function calculates the first 3 order connectivities for each gene
//' and returns the list of vectors of connectivities.
//' @author Qi Gao
//' @param data Input gene by gene coexpression matrix 
//' @param thres Gene pairs with coexpression exceed thres would be assigned an edge between them in the coexpression network
//' @param m Sample size used for the calculation of 3rd order connectivities
//' @param abso Whether to calculate connectivities in absolute network (TRUE) or positive network (FALSE)
//' @param niter Number of sample used for the calculation of 3rd order connectivities
//' @return List of connectivities C1, C2, and C3
//' @export
// [[Rcpp::export]]

List cal_conn(arma::mat data, double thres = 3, int m = 10, bool abso = 1, int niter = 100){
  data.diag().zeros();
  int n = data.n_rows;
  
  arma::mat data_thres(n, n);
  
  if (abso){
    data_thres = arma::conv_to<arma::mat>::from(abs(data) >= thres);
  } else {
    data_thres = arma::conv_to<arma::mat>::from(data >= thres);
  }
  
  arma::vec D1(n);
  arma::vec weight(n);
  for(int i = 0; i < n; i++){
    D1(i) = sum(data_thres.row(i));
  }
  
  arma::vec D2(n);
  arma::vec D3(n);
  arma::vec D3_temp(niter);
  arma::uvec cur_set(m);
  int nz = 0;
  
  arma::uvec allid = arma::linspace<arma::uvec>(0, n-1, n);
  for (int i = 0; i < n; i++){
    if (D1[i] > 2){
      arma::uvec idx = find(data_thres.row(i) == 1);
      arma::mat data_temp = data_thres.rows(idx);
      
      for(int k = 0; k < n; k++){
        weight(k) = sum(data_temp.col(k));
      }
      
      D2[i] = accu(weight.elem(idx)) / D1[i] / (D1[i]-1);
      
      
      nz = sum(weight > 0);
      
      if (nz > m){
        for (int j = 0; j < niter; j++){
          cur_set = Rcpp::RcppArmadillo::sample(allid, m, FALSE, weight);
          
          arma::mat data_temp1 = data_thres.submat(cur_set, cur_set);
          
          D3_temp[j] = accu(data_temp1) / m / (m-1);
        }
        D3[i] = mean(D3_temp);
      } else{
        cur_set = find(weight > 0);
        arma::mat data_temp1 = data_thres.submat(cur_set, cur_set);
        
        D3[i] = accu(data_temp1) / nz / (nz-1);
      }
      
      
    } 
  }
  
  
  
  return List::create(Rcpp::Named("C1") = D1,
                      Rcpp::Named("C2") = D2,
                      Rcpp::Named("C3") = D3);
}
