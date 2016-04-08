#include <QtCore/QCoreApplication>

#include <iostream>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <exception>
#include <cassert>
#include <stdexcept>
#include <ctime>

#include "basestock_initialization_algorithm.h"
#include "greedy_algorithm.h"
#include "poisson_distribution.h"
#include "normal_distribution.h"

// ensure access to floating point operations registers, needed to detect underflows
#pragma fenv_access (on)

int main ( int argc, char *argv[] ) {
	
	QCoreApplication a ( argc, argv );

	QList<double> *targetAggregateFillRates = new QList<double>();
	targetAggregateFillRates->append(0.95);
	targetAggregateFillRates->append(0.95);
	targetAggregateFillRates->append(0.95);
	
	// retailers = 3, products = 5
	TwoEchelonDistributionNetwork *network = new TwoEchelonDistributionNetwork(3, 22); // 2, 1
	network->loadFromFile("preprocessesed-settings.csv", "demand.csv");
	
	//qDebug() << network->getArrivalRateAtWarehouse(1);
	//qDebug() << network->getArrivalRateAtRetailer(1, 1);
	//qDebug() << network->getArrivalRateAtRetailer(1, 2);
	//qDebug() << network->getArrivalRateAtRetailer(1, 3);


	// set arrival rates
	/*network->setArrivalRateAtWarehouse(1, 300.3); 
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

	//network->setBaseStockLevelAtWarehouse(1, 300);*/
	

	// ----------------------------------------------------------------------------------------------------------- optimizatia --

	GreedyAlgorithm *gerrit = new GreedyAlgorithm();
	
	BasestockInitializationAlgorithm *bia = new BasestockInitializationAlgorithm();
	bia->run(network, targetAggregateFillRates);

	std::cout << "initialized" << std::endl;

	time_t start;
	time(&start);

	int result = gerrit->optimizeNetwork(network, targetAggregateFillRates);

	time_t stop;
	time(&stop);

	int duration = difftime(stop, start);

	int seconds = duration % 60;
	int minutes = floor(duration / 60);
	int hours = floor(duration / 3600 );


	std::cout << "optimization duration: " << duration << "s" << std::endl;

	// ------------------------------------------------------------------------------------------------------------- evaluation --

	// set base-stock levels

	//setBaseStockLevelAtRetailer(1, 3, 1);

	//network->setBaseStockLevelAtWarehouse(1, 2);
	//network->setBaseStockLevelAtRetailer(1, 1, 1);
	//network->setBaseStockLevelAtRetailer(1, 2, 1);
	/*
	network->setBaseStockLevelAtWarehouse(2, 2);
	network->setBaseStockLevelAtRetailer(2, 1, 3);
	network->setBaseStockLevelAtRetailer(2, 2, 5);*/

	QList<double> EBOj = gerrit->evaluateNetwork(network);
	///network->setBaseStockLevelAtRetailer(5, 3, 2);

	///EBOj = gerrit->evaluateNetwork(network);

	std::cout << "********************" << std::endl;

	for (int j = 1; j <= network->sizeRetailers(); j++){
		double Mj = 0.0;
		for (int i = 1; i <= network->sizeProducts(); i++){
			Mj = Mj + network->getArrivalRateAtRetailer(i, j);
		}
		std::cout << "demand at j=" << j << " equals: " << Mj << std::endl;
		std::cout << "EBOj for j=" << j << " equals: " << EBOj[j - 1] << std::endl;
		std::cout << "beta star for j=" << j << " equals: " << (1 - (EBOj[j - 1] / Mj)) << std::endl;
	} // for

	// -------------------------------------------------------------------------------------------------------- write to files --

	network->writeBaseStockLevelsToFile("base-stock-levels.txt");

	network->writeBaseStockLevelsToExcel("base-stock-levels.xls");


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
