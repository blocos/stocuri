#ifndef NEGATIVE_BINOMIAL_DISTRIBUTION_H
#define NEGATIVE_BINOMIAL_DISTRIBUTION_H

#include <iostream>
#include <cassert>
#include <cfenv>
#include <exception>
#include <cmath>

class NegativeBinomialDistribution {

	public:
		NegativeBinomialDistribution();
		~NegativeBinomialDistribution();

		double probability(double p, double k, double x);

};

#endif