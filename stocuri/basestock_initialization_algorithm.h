#ifndef BASESTOCK_INITIALIZATION_ALGORITHM
#define BASESTOCK_INITIALIZATION_ALGORITHM

#include <iostream>
#include <cassert>

#include "two_echelon_distribution_network.h"
#include "poisson_distribution.h"

#include <QtCore>


class BasestockInitializationAlgorithm {
	private:

		TwoEchelonDistributionNetwork *network;

		double calculateBetaJ(double BigM, int retailer);
	public:
		BasestockInitializationAlgorithm();
		~BasestockInitializationAlgorithm();

		int run(TwoEchelonDistributionNetwork *network, QList<double> *targetAggregateFillRates);
};

#endif

