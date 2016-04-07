#ifndef GREEDY_ALGORITHN_H
#define GREEDY_ALGORITHN_H

#include <iostream>
#include <cmath>
#include <cassert>

#include <QtCore>

#include "poisson_distribution.h"
#include "normal_distribution.h"
#include "negative_binomial_distribution.h"

#include "two_echelon_distribution_network.h"

class GreedyAlgorithm {
	private:	

		bool GRAVES = true;
		double BOCUTOFF = 0.00000001;
		bool debug = true;

		int n = 0;
		int w = 0;

		TwoEchelonDistributionNetwork *network;

		// warehouse ------------------------------------------------------------------------------------------------------------
		
		// probabilities
		double pPartsOnOrderAtWarehouse(int product, int x);
		double pPartsOnHandAtWarehouse(int product, int x);
		double pPartsOnBackorderAtWarehouse(int product, int x);
		
		// expected values
		double ePartsOnHandAtWarehouse(int product);
		double ePartsOnBackorderAtWarehouse(int product);

		// special
		double pPartsOnBackorderAtWarehouseFromRetailer(int product, int retailer, int x);
		double ePartsOnBackorderAtWarehouseFromRetailer(int product, int retailer);


		// retailer -------------------------------------------------------------------------------------------------------------
		
		// here
		double pPartsOnHandAtRetailer(int product, int retailer, int x);
		double pPartsOnBackorderAtRetailer(int product, int retailer, int x);

		double ePartsOnHandAtRetailer(int product, int retailer);
		double ePartsOnBackorderAtRetailer(int product, int retailer);
		
		double vPartsOnBackorderAtWarehouse(int product);
		double vPartsOnBackorderAtWarehouseFromRetailer(int product, int retailer);

		double ePartsOnOrderAtRetailer(int product, int retailer);
		double vPartsOnOrderAtRetailer(int product, int retailer);
		
		// support --------------------------------------------------------------------------------------------------------------

		bool isTargetFillRatesSatisfied(TwoEchelonDistributionNetwork *network, QList<double> *EBOj, QList<double> *targetAggregateFillRates);
		unsigned long long binomialCoefficient(unsigned long long n, unsigned long long k);
		double calculateDeltaEBO(int product, int j, QList<double> *targetAggregateFillRates);


	public:
		GreedyAlgorithm();
		~GreedyAlgorithm();

		double pPartsOnOrderAtRetailer(int product, int retailer, int x);
		double pPartsOnOrderAtRetailer2Moment(int product, int retailer, int x);

		void setNetwork(TwoEchelonDistributionNetwork *network){ this->network = network; debug = true; };

		QList<double> evaluateNetwork(TwoEchelonDistributionNetwork *network);
		int optimizeNetwork(TwoEchelonDistributionNetwork *network, QList<double> *targetAggregateFillRates);

};

#endif

