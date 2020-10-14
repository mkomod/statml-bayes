#include <Rcpp.h>
#include <cmath>


double laplaceInvCdf(double p, double m, double s);

// [[Rcpp::export]]
Rcpp::NumericVector
rlaplace(const int n, const double m=0, const double sd = 1) {
    Rcpp::NumericVector unif_vals = Rcpp::runif(n);
    Rcpp::NumericVector vals(n);
    for (int i = 0; i < unif_vals.size(); ++i) {
	vals[i] = laplaceInvCdf(unif_vals[i], m, sd);
    }
    return vals;
}

double
laplaceInvCdf(double p, double m, double s) {
    return m - s * R::sign(p - 0.5) * std::log(1 - 2 * std::abs(p - 0.5));
}

