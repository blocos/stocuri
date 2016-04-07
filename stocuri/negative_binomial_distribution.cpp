#include "negative_binomial_distribution.h"


NegativeBinomialDistribution::NegativeBinomialDistribution() {
}


NegativeBinomialDistribution::~NegativeBinomialDistribution() {

}


double NegativeBinomialDistribution::probability(double p, double k, double x) {
	double result = 0.0;

	// clear floating point exceptions
	std::feclearexcept(FE_ALL_EXCEPT);


	// y over x

	// exp(lgamma(r + k) - lgamma(r) - lgamma(k + 1)) * pow(p, r) * pow((1 - p), k)
	//double over = exp(lgammal(x + k) - lgammal(x) - lgammal(k + 1));

	// y over x
	double over = exp(lgammal(x+k-1 + 1) - lgammal(x + 1) - lgammal(x+k-1 - x + 1));

	result = over * pow(p, x) * pow(1 - p, k);

	std::cout << "o:" << over << std::endl;

	// detect possible underflow on exponential with negative lambda
	if (std::fetestexcept(FE_OVERFLOW)) {
		//std::cout << over << std::endl;
		throw std::overflow_error("floating point operation failure in NegativeBinomialDistribution.probability");
	} // if

	return result;

} // probability


double NegativeBinomialDistribution::probabilityByNormalApproximation(double p, double k, double x) {
	double result = 0.0;

	double mu = k / p;
	double sigma = sqrt(k * ((1 - p) / (p*p)));

	NormalDistribution kees;

	result = kees.probability(mu, sigma, x);

	//	std::cout << "result N.pdf: " << result << std::endl;
	//std::cout << "mu: " << mu << std::endl;
	//std::cout << "sigma: " << sigma << std::endl;

	if (result < 0) {
		result = 0;
	}

	return result;
} // probabilityByNormalApproximation