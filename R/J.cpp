#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

double h(arma::rowvec x, arma::rowvec theta);
double pm(arma::rowvec x, arma::rowvec theta);
double pn(arma::rowvec x);
double pn2(arma::rowvec x, arma::mat S);
double pm2(arma::rowvec x, arma::mat B, double c);
double h2(arma::rowvec x, arma::mat S, arma::mat B, double c);

// [[Rcpp::export]]
double J(arma::rowvec theta, arma::mat X, arma::mat Y, arma::mat S) {
    int N = X.n_rows;
    int P = X.n_cols;

    arma::mat B(P, P);
    for (int i = 0; i < P; ++i) {
	for (int j = 0; j < P; ++j) {
	    // Rcpp::Rcout << theta(i*P + j) << "\t";
	    B(i, j) = theta(i*P + j);
	}
    }
    // Rcpp::Rcout << "\n" << A;
    double c = theta(16); 

    double t = 0;
    arma::rowvec x, y;
    for (int i=0; i < N; ++i) {
	x = X.row(i);
	y = Y.row(i);
	// t += std::log(h(x, theta)) + std::log(1 - h(y, theta));
	t += std::log(h2(x, S, B, c)) + std::log(1 - h2(y, S, B, c));
    }
    return -t / (2*N);
}

// [[Rcpp::export]]
double h2(arma::rowvec x, arma::mat S, arma::mat B, double c) {
    double Gu = std::log(pm2(x, B, c)) - std::log(pn2(x, S));
    return 1.0 / ( 1 + std::exp(-Gu));
}

// [[Rcpp::export]]
double pm2(arma::rowvec x, arma::mat B, double c) {
    // Rcpp::Rcout << arma::abs(B * x.t());
    double t = std::exp(-std::sqrt(2) * arma::sum(arma::abs(B * x.t())) + c);
    return t;
}

// [[Rcpp::export]]
double pn2(arma::rowvec x, arma::mat S) {
    double P = x.n_cols;
    arma::mat Sinv = S.i();
    double Sdet = arma::det(S);
    arma::mat xSx = (x * Sinv * x.t());
    double res = std::exp(- xSx(0, 0) /2.0) / 
	( std::pow(2*arma::datum::pi, P/2.0) * std::sqrt(Sdet) );
    return res;
}

double h(arma::rowvec x, arma::rowvec theta) {
    double Gu = std::log(pm(x, theta)) - std::log(pn(x));
    return 1.0 / ( 1 + std::exp(-Gu));
}

// [[Rcpp::export]]
double pm(arma::rowvec x, arma::rowvec theta) {
    double x_ = x(0);
    double b = 1 / std::sqrt(2);
    double c = theta(0);
    double t = std::exp(-std::sqrt(2) * arma::sum(arma::abs(x))) * c;
    return t;
}

// [[Rcpp::export]]
double pn(arma::rowvec x) {
    double P = x.n_cols;
    double res = std::exp(-dot(x, x)/2.0) / 
	std::pow(2*arma::datum::pi, P/2.0);
    return res;
}
