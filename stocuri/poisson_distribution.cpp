#include "poisson_distribution.h"

PoissonDistribution::PoissonDistribution() {
}


PoissonDistribution::~PoissonDistribution() {
}


double PoissonDistribution::recursiveProbability(double lambda, double x){
	// pre	: True
	// ret	: [ (lambda/x)*recursiveProbability(lambda, x-1) | x > 0 ] && [1 | x == 0]

	// pre-conditions satisfied

	double result;

	if (x == 0) {
		result = 1.0;
	} else {
		result = (lambda / x) * recursiveProbability(lambda, x - 1);
	} // eif

	return result;

} // recursiveProbability


double PoissonDistribution::probability(double lambda, double x) {
	// pre	: lambda >= 0 /\ x >= 0
	// ret	: exception \/ P[X = x], given X~Poisson(lambda)

	assert(lambda >= 0);
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	// clear floating point exceptions
	std::feclearexcept(FE_ALL_EXCEPT);

	result = exp(-lambda);

	// detect possible underflow on exponential with negative lambda
	if (std::fetestexcept(FE_UNDERFLOW)) {
		throw std::underflow_error("underflow in PoissonDistribution.pdf");
	} // if

	result = recursiveProbability(lambda, x) * result;

	// detect possible underflow after recursion
	if (std::fetestexcept(FE_UNDERFLOW)) {
		throw std::underflow_error("underflow in PoissonDistribution.pdf");
	} // if

	return result;

} // pdf


double PoissonDistribution::probabilityBySterlingApproximation(double lambda, double x) {
	// pre	: True
	// ret	: PDF based on Sterling's Approximation

	// pre-conditions satisfied

	double result = exp(x * log(lambda) - lgamma(x + 1.0) - lambda);
	
	return result;

} // pdfSterlingApproximation


double PoissonDistribution::probabilityByNormalApproximation(double lambda, double x) {
	// pre	: True
	// ret	: P[X-normal = x], given X~Poisson(lambda) and X-normal~Normal(mu=lambda, sigma=sqrt(lambda))
	
	double result = 0.0;

	NormalDistribution kees;
	result = kees.probability(lambda, sqrt(lambda), x);

	return result;

} // probabilityByNormalApproximation
