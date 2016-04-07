#include "negative_binomial_distribution.h"


NegativeBinomialDistribution::NegativeBinomialDistribution() {
}


NegativeBinomialDistribution::~NegativeBinomialDistribution() {

}


double NegativeBinomialDistribution::probability(double p, double k, double x) {
	// pre	: x >= 0
	// ret	: exception \/ P[X = x], given X~NegativeBinomial(p,k)

	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	// clear floating point exceptions
	std::feclearexcept(FE_ALL_EXCEPT);

	// y over x, use logarithmic gamma function and transform using exponential
	double binCoef = exp(lgammal(x+k-1 + 1) - lgammal(x + 1) - lgammal(x+k-1 - x + 1));

	result = binCoef * pow(p, x) * pow(1 - p, k);

	// detect possible floating point operation exception
	if (std::fetestexcept(FE_ALL_EXCEPT)) {
		throw std::runtime_error("floating point operation failure in NegativeBinomialDistribution.probability");
	} // if

	return result;

} // probability


double NegativeBinomialDistribution::probabilityByNormalApproximation(double p, double k, double x) {
	// pre	: True
	// ret	: P[X-normal = x], given X~NegativeBinomial(p,k) and X-normal~Normal(mu=k/p, sigma=sqrt(k*(1-p)/(p*p)))

	// pre-conditions satisfied

	double result = 0.0;

	double mu = k / p;
	double sigma = sqrt(k * ((1 - p) / (p*p)));

	NormalDistribution kees;

	result = kees.probability(mu, sigma, x);

	return result;

} // probabilityByNormalApproximation