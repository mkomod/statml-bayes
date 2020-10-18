#include <cmath>
#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double lnpm0(arma::rowvec x, arma::rowvec b);
double lnpn(arma::rowvec y, arma::mat S_ing, double log_NormalisingConst);


// [[Rcpp::export]]
double J(NumericVector theta, NumericMatrix X, NumericMatrix Y, 
	arma::mat S_inv, double log_NormalisingConst) {
    int P = X.ncol();
    int N = X.nrow();
    double c = theta(std::pow(P, 2));
    arma::mat A(P, P);
    for (int i = 0; i < P; ++i) {
	for (int j = 0; j < P; ++j) {
	    A(i, j) = theta(i*P + j);		// fill
	}
    }    
    arma::mat B = A.i();			// Invert

    arma::rowvec x(P);
    arma::rowvec y(P);  
    arma::rowvec b(P);
    double t = 0, qx, qy, px, py, hx, hy;
    for (int i = 0; i < N; ++i) {
	x = X.row(i);
	y = Y.row(i);
	for (int j = 0; j < P; ++j) {
	    b = B.row(j);
	    px += lnpm0(x, b);
	    py += lnpm0(y, b);
	}
	px += c; py += c;			// add on c
	px = exp(px); py = exp(py);
	qx = exp(lnpn(x, S_inv, log_NormalisingConst));
	qy = exp(lnpn(y, S_inv, log_NormalisingConst));
	hx = px / (px - qx);
	hy = py / (py - qy);
	Rcout << "hx: " << hx << " hy: " << hy;
	t += std::log(hx) + std::log(1.0 - hy);
    }
    return -t / (2 * N);
}


// [[Rcpp::export]]
double
lnpm0(arma::rowvec x, arma::rowvec b) {
    double t = -sqrt(2) * arma::dot(x, b);	// dot product of a, b
    return t;
}

// [[Rcpp::export]]
double
lnpn(arma::rowvec y, arma::mat S_inv, double log_NormalisingConst) {
    arma::mat yt = y.t();
    arma::mat a = (y * S_inv * yt);
    double t = -1.0/2.0 * a(0, 0) - log_NormalisingConst;
    return t;
}
