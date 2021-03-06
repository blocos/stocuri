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

	double px = pow(p, x);
	double p1x = pow(1 - p, k);


	result = binCoef * px * p1x;

	
	// detect possible floating point operation exception
	if (std::fetestexcept(FE_OVERFLOW)) {
		throw std::overflow_error("overflow in NegativeBinomialDistribution.probability");
	}

	if (std::fetestexcept(FE_UNDERFLOW)) {
		throw std::underflow_error("underflow in NegativeBinomialDistribution.probability");
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

	//std::cout << mu << std::endl;
	//std::cout << sigma << std::endl;

	NormalDistribution kees;

	// check z
	double z = (x - mu) / sigma;

	if (z <= -30) {
		result = 0.0;
	}
	else {

		try{
			result = kees.probability(mu, sigma, x);
		}
		catch (std::exception &me) {
			std::cout << me.what() << std::endl;

			std::cout << "z: " << (x - mu) / sigma << std::endl;

			result = 0.0;
		}
	}

	return result;

} // probabilityByNormalApproximation