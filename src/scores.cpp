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


struct ip {
  const Matrix<double, Dynamic, Dynamic> S_;
  const Matrix<double, Dynamic, 1> y_;
  const Matrix<double, Dynamic, Dynamic> x_;
  const Matrix<double, Dynamic, Dynamic> xs_;
  ip(const Matrix<double, Dynamic, Dynamic>& S,
     const Matrix<double, Dynamic, 1>& y,
     const Matrix<double, Dynamic, Dynamic>& x, 
     const Matrix<double, Dynamic, Dynamic>& xs) : S_(S),y_(y),x_(x),xs_(xs) { }
  template <typename T>
  T operator()(const Matrix<T, Dynamic, 1>& Gvec) const{
    int n = S_.rows();
    int m = S_.cols();
    int Q = x_.cols();
    Matrix<T, Dynamic, 1> g1 = Gvec.topRows(m);
    Matrix<T, Dynamic, 1> g2 = Gvec.bottomRows(n*m);
    Matrix<T, Dynamic, Dynamic> G = to_matrix(g2 , m , n);
    Matrix<T, Dynamic, Dynamic> SG = stan::math::multiply(S_,G);
    Matrix<T, Dynamic, Dynamic> xr = stan::math::multiply(SG,x_);
    Matrix<T, Dynamic, Dynamic> xrs = stan::math::multiply(SG,xs_);
    
    Matrix<T, Dynamic, 1> d = stan::math::multiply(S_,g1);
    Matrix<T, Dynamic, 1> yd = y_-d;
    
    Matrix<T, Dynamic, 1> dif1;
    Matrix<T, Dynamic, 1> dif2;
    
    T term1 = 0;
    T term2 = 0;
    
    T ip1;
    T ip2;
    
    for (int j=0; j<Q; j++){
      dif1 = xr.col(j) - xrs.col(j);
      dif2 = (yd-xr.col(j));
      ip1 = dif1.transpose() * dif1;
      ip2 = dif2.transpose() * dif2;
      term1 += sqrt(ip1);
      term2 += sqrt(ip2);
    }
    
    
    return (((0.5 * term1 ) - term2) / Q );
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
// [[Rcpp::export(name=".energy_i")]]
Rcpp::List energy_i(Rcpp::NumericMatrix Sin, 
                Rcpp::NumericVector yin,
                Rcpp::NumericMatrix xin, 
                Rcpp::NumericMatrix xsin,
                Rcpp::NumericVector Gin){
   Map<Matrix<double,Dynamic,Dynamic>> S(Rcpp::as<Map<Matrix<double,Dynamic,Dynamic>> >(Sin));
   Map<Matrix<double,Dynamic, 1>> y(Rcpp::as<Map<Matrix<double,Dynamic, 1>> >(yin));
   Map<Matrix<double,Dynamic, Dynamic>> x(Rcpp::as<Map<Matrix<double,Dynamic, Dynamic>> >(xin));
   Map<Matrix<double,Dynamic, Dynamic>> xs(Rcpp::as<Map<Matrix<double,Dynamic, Dynamic>> >(xsin));
   Map<Matrix<double,Dynamic, 1>> G(Rcpp::as<Map<Matrix<double,Dynamic, 1>> >(Gin));
   ip f(S,y,x,xs);
   double fx;
   Matrix<double, Dynamic, 1> grad_fx;
   stan::math::gradient(f, G, fx, grad_fx);
   return Rcpp::List::create(Rcpp::Named("grad") = Rcpp::wrap(grad_fx), 
                             Rcpp::Named("val") = fx);
 }


