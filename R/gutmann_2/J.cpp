#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

double h(arma::rowvec x, arma::rowvec theta);
double pm(arma::rowvec x, arma::rowvec theta);
double pn(arma::rowvec x);

// [[Rcpp::export]]
double J(arma::rowvec theta, arma::mat X, arma::mat Y) {
    int N = X.n_rows;
    int P = X.n_cols;
    
    double t=0;
    arma::rowvec x, y;
    for (int i=0; i < N; ++i) {
	x = X.row(i);
	y = Y.row(i);
	t += std::log(h(x, theta)) + std::log(1 - h(y, theta));
    }
    return -t / (2.0 * N);
}

double h(arma::rowvec x, arma::rowvec theta) {
    double Gu = std::log(pm(x, theta)) - std::log(pn(x));
    return 1.0 / ( 1 + std::exp(-Gu));
}

// [[Rcpp::export]]
double pm(arma::rowvec x, arma::rowvec theta) {
    double x_ = x(0);
    double m = 0;
    double b = 1 / std::sqrt(2);
    double c = theta(0);
    double t = c * 1 / (2 * b) * std::exp(- std::abs(x_ - m) / b);
    return t;
}

// [[Rcpp::export]]
double pn(arma::rowvec x) {
    int P = x.n_rows;
    double res = std::exp(dot(x, x))/std::sqrt( std::pow(2 * PI, P));
    return  res;
}

