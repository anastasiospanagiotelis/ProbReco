// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]

#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp14)]]

using Eigen::Map;  
using Eigen::Dynamic;
using Eigen::Matrix;  

//The following template needs to be defined to do AD

struct scoretemp {
  //Following inputs are all constant w.r.t. differentiation
  const Matrix<double, Dynamic, Dynamic> S_; // Linear constraints matrix
  const Matrix<double, Dynamic, 1> y_;  // Realisation
  const Matrix<double, Dynamic, Dynamic> x_; // Draws from probabilistic forecast
  const Matrix<double, Dynamic, Dynamic> xs_; // Draws from probabilistic forecast
  int score_;//Score to be used
  double alpha_; //Additional parameters
  scoretemp(const Matrix<double, Dynamic, Dynamic>& S,
     const Matrix<double, Dynamic, 1>& y,
     const Matrix<double, Dynamic, Dynamic>& x, 
     const Matrix<double, Dynamic, Dynamic>& xs, 
     int score,
     double alpha) : S_(S),y_(y),x_(x),xs_(xs),score_(score),alpha_(alpha) { }
  template <typename T>
  T operator()(const Matrix<T, Dynamic, 1>& Gvec) const{
    int n = S_.rows();  // Number of series
    int m = S_.cols();  // Number of base series 
    int Q = x_.cols();  // Number of draws from probabilistic forecast
    Matrix<T, Dynamic, 1> g1 = Gvec.topRows(m); //First m parameters for translation
    Matrix<T, Dynamic, 1> g2 = Gvec.bottomRows(n*m); //Remaining parameters for matrix
    Matrix<T, Dynamic, Dynamic> G = to_matrix(g2 , m , n); //Put into matrix
    Matrix<T, Dynamic, Dynamic> SG = stan::math::multiply(S_,G); // S times G (needs stan::math::multiply)
    Matrix<T, Dynamic, Dynamic> xr = stan::math::multiply(SG,x_); // Draw from reconciled forecast (needs stan::math::multiply)
    Matrix<T, Dynamic, Dynamic> xrs = stan::math::multiply(SG,xs_); // Draw from reconciled forecast (needs stan::math::multiply)
    
    Matrix<T, Dynamic, 1> d = stan::math::multiply(S_,g1);  // Translation
    
    T out; // output
    if (score_==1){
      
      T term1 = 0; // Intialise terms for cumulative sum
      T term2 = 0; // Intialise terms for sumulative sum    
      
      
      Matrix<T, Dynamic, 1> yd = y_-d;  // y-d
      
      Matrix<T, Dynamic, 1> dif1; // Initialise vectors
      Matrix<T, Dynamic, 1> dif2; // Initialise vectors
      

      
      for (int q=0; q<Q; q++){
        dif1 = xr.col(q) - xrs.col(q); // Difference for first term (translation cancels)
        dif2 = (yd-xr.col(q)); // Difference for second term (translation cancels)
        term1 += stan::math::pow(dif1.norm(),alpha_); // Update sum
        term2 += stan::math::pow(dif2.norm(),alpha_); // Update sum
        
      }
      out = ( term2 - (0.5 * term1 ) ) / Q ;
    }
    if (score_==2){
      //Variogram to go here
      out = 0;
      T dif;
      T term1 = 0 ;
      T term2 = 0 ;
      for (int i=0; i<= (n -1); i++){
        for(int j=(i+1);j<n; j++){
          dif = stan::math::abs( y_(i) - y_(j) ); //Difference in y
          term1 = stan::math::pow( dif, alpha_ );
          term2 = 0;
          for (int q=0; q<Q; q++){
            dif = stan::math::abs( d(i) - d(j) + xr(i,q) - xrs(j,q) );
            term2 += stan::math::pow( dif, alpha_ );
          }
          term2 = term2 /Q; //Divide by sample size
          out += stan::math::pow(( term1 - term2 ),2);
        }
      }
      out=2*out;
  
    }
    
    return (out); // Energy score
  }
};

// Computes energy score contribution with unconstrained G
// 
// Sin Matrix encoding linear constraints
// yin Vector valued realisation
// xin Vector valued realisation from base predictive
// xsin Copy of vector valued realisation from base predictive
// Gin Reconciliation matrix (vectorised)
// Contribution to energy score and gradient w.r.t G
// [[Rcpp::export(name=".score")]]
Rcpp::List score(Rcpp::NumericMatrix Sin, // Constraints matrix
                Rcpp::NumericVector yin, // Realisations
                Rcpp::NumericMatrix xin,// Draws from probabilistic forecasts
                Rcpp::NumericMatrix xsin, // Draws from probabilistic forecasts
                Rcpp::NumericVector Gin,
                int scorein,
                double alphain) {
  
   //Convert from R objects to Eigen objects
   Map<Matrix<double,Dynamic,Dynamic>> S(Rcpp::as<Map<Matrix<double,Dynamic,Dynamic>> >(Sin));
   Map<Matrix<double,Dynamic, 1>> y(Rcpp::as<Map<Matrix<double,Dynamic, 1>> >(yin));
   Map<Matrix<double,Dynamic, Dynamic>> x(Rcpp::as<Map<Matrix<double,Dynamic, Dynamic>> >(xin));
   Map<Matrix<double,Dynamic, Dynamic>> xs(Rcpp::as<Map<Matrix<double,Dynamic, Dynamic>> >(xsin));
   Map<Matrix<double,Dynamic, 1>> G(Rcpp::as<Map<Matrix<double,Dynamic, 1>> >(Gin));
   
   scoretemp f(S,y,x,xs,scorein,alphain); // Instantiate template
   double fx; // Value of function
   Matrix<double, Dynamic, 1> grad_fx; // Initialise gradient
   stan::math::gradient(f, G, fx, grad_fx); //Automatic differentiation
   // Convert back to list
   return Rcpp::List::create(Rcpp::Named("grad") = Rcpp::wrap(grad_fx), 
                             Rcpp::Named("val") = fx);
 }


