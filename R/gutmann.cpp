#include <cmath>
#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double lnpm0(arma::vec x, arma::rowvec b);
double lnpn(arma::vec y, arma::mat S_ing, double log_NormalisingConst);


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

    arma::vec x(P);
    arma::vec y(P);  
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
	qx = lnpn(x, S_inv, log_NormalisingConst);
	qy = lnpn(y, S_inv, log_NormalisingConst);
	hx = 1.0 / (1.0 + exp(qx - px));
	hy = 1.0 / (1.0 + exp(qy - py));
	t += std::log(hx) + std::log(1.0 - hy);
    }
    return -t;
}


double
lnpm0(arma::vec x, arma::rowvec b) {
    double t = -sqrt(2) * arma::dot(x, b);	// dot product of a, b
    return t;
}

double
lnpn(arma::vec y, arma::mat S_inv, double log_NormalisingConst) {
    arma::mat yt = y.t();
    arma::mat a = (yt * S_inv * y);
    double t = a(0, 0) - log_NormalisingConst;
    return t;
}
