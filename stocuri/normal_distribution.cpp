#include "normal_distribution.h"


NormalDistribution::NormalDistribution()
{
}


NormalDistribution::~NormalDistribution()
{
}

double NormalDistribution::probability(double mu, double sigma, double x) {

	//std::cout << "mu: " << mu << std::endl;
	//std::cout << "sigma: " << sigma << std::endl;
	//std::cout << "x: " << x << std::endl;

	double result = 0.0;

	double pi = 3.1415926535897;

	double z = (x - mu) / sigma;

	//std::cout << "z: " << z << std::endl;
	//std::cout << "1.0/sigma: " << 1.0 / sigma << std::endl;
	///std::cout << "(1.0 / sqrt(2.0 * pi)): " << (1.0 / sqrt(2.0 * pi)) << std::endl;
	//std::cout << "-0.5*x*x: " << -0.5*z*z << std::endl;
	//std::cout << "exp(-0.5*x*x): " << exp(-0.5*z*z) << std::endl;

	result = 1.0/sigma * (1.0 / sqrt(2.0 * pi)) * exp(-0.5*z*z);
	
	std::feclearexcept(FE_ALL_EXCEPT);


	// detect possible underflow
	if (std::fetestexcept(FE_UNDERFLOW)) {
		result = 0.0;
	}

	return result;
}
