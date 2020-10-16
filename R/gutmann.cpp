#include <cmath>
#include <RcppArmadillo.h>

using namespace Rcpp;

double lnpm(NumericVector x, NumericVector b, double c);
double lnpn(NumericVector y);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double J(NumericVector theta, NumericMatrix X, NumericMatrix Y) {
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

    double t = 0;
    for (int i = 0; i < N; ++i) {
	NumericVector x = X.row(i);
	NumericVector y = Y.row(i);
	NumericVector b(4);
	double px = 0;
	double py = 0;
	for (int j = 0; j < P; ++j) {
	    for (int k = 0; k < P; ++k) {
		b[k] = B(j, k);
	    }
	    px = lnpm(X(i, _), b, c);
	    py = lnpm(Y(i, _), b, c);
	}
	double qx = lnpn(x);
	double qy = lnpn(y);
	double hx = 1.0 / (1.0 + exp(qx - px));
	double hy = 1.0 / (1.0 + exp(qy - py));
	t = t + std::log(hx) + std::log(1.0 - hy);
    }
    return -t;
}


double
lnpm(NumericVector x, NumericVector b, double c) {
    double t = 0;
    for (int i = 0; i < x.size(); ++i) {
	t += std::abs(x[i] * b[i]);
    }
    return -sqrt(2) * t + c;
}

double
lnpn(NumericVector y) {
    double t = sum(dnorm(y, true));
    return t;
}
