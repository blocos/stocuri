#ifndef GREEDY_ALGORITHN_H
#define GREEDY_ALGORITHN_H

#include <iostream>
#include <cmath>
#include <cassert>

#include <QtCore>

#include "poisson_distribution.h"
#include "normal_distribution.h"

#include "two_echelon_distribution_network.h"

class GreedyAlgorithm {
	private:	

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
		
		double pPartsOnOrderAtRetailer(int product, int retailer, int x);
		double pPartsOnHandAtRetailer(int product, int retailer, int x);
		double pPartsOnBackorderAtRetailer(int product, int retailer, int x);
		
		/*
		

		double ePartsOnHandAtWJ(int i, int j);
		double ePartsOnBackorderAtWJ(int i, int j);

		void initializeGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ);
		void initializeDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ);

		void clearS(QList<QList<double>*> *aS);
		void clearGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ);
		void clearDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ);

		bool stopingCriterionMet();
		
		
		//unsigned long long int factorial(unsigned long long int x);
		//unsigned binomialCoef(unsigned n, unsigned k);
		*/

	public:
		GreedyAlgorithm();
		~GreedyAlgorithm();

		void evaluateNetwork(TwoEchelonDistributionNetwork *network);
		void optimizeNetwork(TwoEchelonDistributionNetwork *network);

};

#endif

