#include <Rcpp.h>
#include <cmath>


double laplaceInvCdf(double p, double m, double b);

// [[Rcpp::export]]
Rcpp::NumericVector
rlaplace(const int n, const double m=0, const double sd = 1) {
    Rcpp::NumericVector unif_vals = Rcpp::runif(n);
    Rcpp::NumericVector vals(n);
    double b = sd / sqrt(2);
    for (int i = 0; i < unif_vals.size(); ++i) {
	vals[i] = laplaceInvCdf(unif_vals[i], m, b);
    }
    return vals;
}

double
laplaceInvCdf(double p, double m, double b) {
    return m - b * R::sign(p - 0.5) * std::log(1 - 2 * std::abs(p - 0.5));
}

