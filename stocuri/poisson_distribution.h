#ifndef POISSON_DISTRIBUTION_H
#define POISSON_DISTRIBUTION_H


#include <iostream>
#include <cassert>
#include <cfenv>
#include <exception>
#include <cmath>

#include "normal_distribution.h"


class PoissonDistribution {

	private:
		double recursiveProbability(double lambda, double x);

	public:
		PoissonDistribution();
		~PoissonDistribution();

		double probability(double lambda, double x);

		double probabilityBySterlingApproximation(double lambda, double x);

		double probabilityByNormalApproximation(double lambda, double x);
};

#endif

