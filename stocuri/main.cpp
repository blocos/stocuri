#include <QtCore/QCoreApplication>

#include <iostream>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <exception>
#include <cassert>
#include <stdexcept>

#include "basestock_initialization_algorithm.h"
#include "greedy_algorithm.h"
#include "poisson_distribution.h"
#include "normal_distribution.h"

// ensure access to floating point operations registers, needed to detect underflows, overflows and so on
#pragma fenv_access (on)

int main ( int argc, char *argv[] ) {

	// settings -----------------------------------------------------------------------------------------------------------------

	bool	USE_MULTI_INCREMENT = true;
	int		MULTI_INCREMENT_MAX_PRODUCTS = 10;
	double	MULTI_INCREMENT_MIN_EBO = 0.0;
	bool	USE_GRAVES = true;
	double	BO_CUT_OFF = 0.000000001;
	double	APPROX = 30;

	QString el = "94"; // el in {94, 100, 7038435, 7038847, 7045042, 7132473, 1060024472, 3010101808 }
	QString custom = "mimp-10";

	int nProducts = 10;

	QString demandFile = "demand-" + el + ".csv";
	QString settingsFile = "preprocessesed-settings-" + el + ".csv";


	// core ---------------------------------------------------------------------------------------------------------------------

	QCoreApplication a ( argc, argv );

	QList<double> *targetAggregateFillRates = new QList<double>();
	targetAggregateFillRates->append(0.95);
	targetAggregateFillRates->append(0.95);
	targetAggregateFillRates->append(0.95);

	// retailers, products
	TwoEchelonDistributionNetwork *network = new TwoEchelonDistributionNetwork(3, nProducts);
	network->loadFromFile(settingsFile, demandFile);
	std::cout << "data loaded" << std::endl;
	
	/*
	// set arrival rates
	network->setArrivalRateAtWarehouse(1, 300.3); 
	network->setArrivalRateAtRetailer(1, 1, 100.1);
	network->setArrivalRateAtRetailer(1, 2, 200.2);

	// set lead times
	network->setLeadTimeToWarehouse(1, 4);
	network->setLeadTimeToRetailer(1, 1, 1);
	network->setLeadTimeToRetailer(1, 2, 1);

	// set inventory holding cost
	network->setInventoryHoldingCostAtWarehouse(1, 1);
	network->setInventoryHoldingCostAtRetailer(1, 1, 1);
	network->setInventoryHoldingCostAtRetailer(1, 2, 1);
	*/	

	// ----------------------------------------------------------------------------------------------------------- optimizatia --

	BasestockInitializationAlgorithm *gaia = new BasestockInitializationAlgorithm();
	std::cout << "initializing base-stock levels..." << std::endl;
	gaia->run(network, targetAggregateFillRates);
	std::cout << "base-stock levels initialized" << std::endl;


	// define relevant times
	time_t start;
	time_t stop;

	time(&start);

	GreedyAlgorithm *gerrit = new GreedyAlgorithm(USE_MULTI_INCREMENT, MULTI_INCREMENT_MAX_PRODUCTS, MULTI_INCREMENT_MIN_EBO, USE_GRAVES, BO_CUT_OFF, APPROX);
	std::cout << "optimizing base-stock levels..." << std::endl;
	int result = gerrit->optimizeNetwork(network, targetAggregateFillRates);

	time(&stop);

	int duration = difftime(stop, start);

	std::cout << "base-stock levels optimized in " << duration << " seconds" << std::endl;

	// ------------------------------------------------------------------------------------------------------------- evaluation --

	std::cout << "********************************" << std::endl;

	/*
	// set base-stock levels
	network->setBaseStockLevelAtWarehouse(1, 2);
	network->setBaseStockLevelAtRetailer(1, 1, 1);
	network->setBaseStockLevelAtRetailer(1, 2, 1);
	*/

	QList<double> EBOj = gerrit->evaluateNetwork(network);

	for (int j = 1; j <= network->sizeRetailers(); j++){
		double Mj = 0.0;
		for (int i = 1; i <= network->sizeProducts(); i++){
			Mj = Mj + network->getArrivalRateAtRetailer(i, j);
		} // for
		std::cout << "demand at j=" << j << " equals: " << Mj << std::endl;
		std::cout << "EBOj for j=" << j << " equals: " << EBOj[j - 1] << std::endl;
		std::cout << "beta star for j=" << j << " equals: " << (1 - (EBOj[j - 1] / Mj)) << std::endl;
	} // for

	// -------------------------------------------------------------------------------------------------------- write to files --

	network->writeBaseStockLevelsToFile("base-stock-levels-" + el + "-" + custom + ".txt");

	// -------------------------------------------------------------------------------------------------------------- clean up --

	delete network;
	network = 0;

	delete targetAggregateFillRates;
	targetAggregateFillRates = 0;

	delete gerrit;
	gerrit = 0;

	// --------------------------------------------------------------------------------------------------- execute application --

	return a.exec();
	
} // main
