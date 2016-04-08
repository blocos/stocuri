#include "normal_distribution.h"


NormalDistribution::NormalDistribution()
{
}


NormalDistribution::~NormalDistribution()
{
}

double NormalDistribution::probability(double mu, double sigma, double x) {
	// pre	: True
	// ret	: P[X = x], given X~Normal(mu, sigma)

	// pre-conditions satisfied

	double result = 0.0;

	double pi = 3.1415926535897;

	std::feclearexcept(FE_ALL_EXCEPT);

	double z = (x - mu) / sigma;

	result = 1.0/sigma * (1.0 / sqrt(2.0 * pi)) * exp(-0.5*z*z);
	
	// detect possible floating point operation error
	if (std::fetestexcept(FE_UNDERFLOW) || std::fetestexcept(FE_OVERFLOW)) {
		result = 0.0;
		//throw(std::runtime_error("floating point operation error in Normal.pdf"));
	} // if

	return result;

} // probability
