#include <QtCore/QCoreApplication>

#include <iostream>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <exception>
#include <cassert>
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
	TwoEchelonDistributionNetwork *network = new TwoEchelonDistributionNetwork(3, 100); // 2, 1
	network->loadFromFile("preprocessesed-settings.csv", "demand.csv");

	// ----------------------------------------------------------------------------------------------------------- optimizatia --

	GreedyAlgorithm *gerrit = new GreedyAlgorithm();
	
	BasestockInitializationAlgorithm *bia = new BasestockInitializationAlgorithm();
	bia->run(network, targetAggregateFillRates);

	std::cout << "initialized" << std::endl;


	time_t start;
	time(&start);

	//int result = gerrit->optimizeNetwork(network, targetAggregateFillRates);

	time_t stop;
	time(&stop);

	int duration = difftime(stop, start);

	int seconds = duration % 60;
	int minutes = floor(duration / 60);
	int hours = floor(duration / 3600 );


	std::cout << "optimization duration: " << duration << "s" << std::endl;

	// ------------------------------------------------------------------------------------------------------------- evaluation --

	QList<double> EBOj = gerrit->evaluateNetwork(network);

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
