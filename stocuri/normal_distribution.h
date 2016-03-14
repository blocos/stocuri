#ifndef NORMAL_DISTRIBUTION_H
#define NORMAL_DISTRIBUTION_H

#include <iostream>
#include <cassert>
#include <cfenv>
#include <exception>
#include <cmath>


class NormalDistribution {

	public:
		NormalDistribution();
		~NormalDistribution();

		double probability(double mu, double sigma, double x);

};

#endif

