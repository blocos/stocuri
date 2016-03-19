#include <QtCore/QCoreApplication>

#include <iostream>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <exception>
#include <cassert>

#include "greedy_algorithm.h"

#include "poisson_distribution.h"
#include "normal_distribution.h"

#include "BasicExcel.hpp"

using namespace YExcel;

// ensure access to floating point operations registers, needed to detect underflows
#pragma fenv_access (on)

int main ( int argc, char *argv[] ) {
	
	QCoreApplication a ( argc, argv );

	QList<double> *targetAggregateFillRates = new QList<double>();
	targetAggregateFillRates->append(0.95);
	targetAggregateFillRates->append(0.95);
	
	// retailers = 2, products = 2
	TwoEchelonDistributionNetwork *network = new TwoEchelonDistributionNetwork(2, 2);

	// set arrival rates
	network->setArrivalRateAtWarehouse(1, 3); 
	network->setArrivalRateAtRetailer(1, 1, 1);
	network->setArrivalRateAtRetailer(1, 2, 2);

	network->setArrivalRateAtWarehouse(2, 4);
	network->setArrivalRateAtRetailer(2, 1, 2);
	network->setArrivalRateAtRetailer(2, 2, 2);
	
	// set lead times
	network->setLeadTimeToWarehouse(1, 6);
	network->setLeadTimeToRetailer(1, 1, 1);
	network->setLeadTimeToRetailer(1, 2, 1);

	network->setLeadTimeToWarehouse(2, 4);
	network->setLeadTimeToRetailer(2, 1, 1);
	network->setLeadTimeToRetailer(2, 2, 1);

	// set inventory holding cost
	network->setInventoryHoldingCostAtWarehouse(1, 1);
	network->setInventoryHoldingCostAtRetailer(1, 1, 1);
	network->setInventoryHoldingCostAtRetailer(1, 2, 1);

	network->setInventoryHoldingCostAtWarehouse(2, 1);
	network->setInventoryHoldingCostAtRetailer(2, 1, 1);
	network->setInventoryHoldingCostAtRetailer(2, 2, 1);

	// ----------------------------------------------------------------------------------------------------------- optimizatia --

	GreedyAlgorithm *gerrit = new GreedyAlgorithm();
	int result = gerrit->optimizeNetwork(network, targetAggregateFillRates);


	// ------------------------------------------------------------------------------------------------------------- evaluation --

	// set base-stock levels
	/*network->setBaseStockLevelAtWarehouse(1, 2);
	network->setBaseStockLevelAtRetailer(1, 1, 1);
	network->setBaseStockLevelAtRetailer(1, 2, 1);
	
	network->setBaseStockLevelAtWarehouse(2, 2);
	network->setBaseStockLevelAtRetailer(2, 1, 3);
	network->setBaseStockLevelAtRetailer(2, 2, 5);*/

	QList<double> EBOj = gerrit->evaluateNetwork(network);

	for (int j = 1; j <= network->sizeRetailers(); j++){
		double Mj = 0.0;
		for (int i = 1; i <= network->sizeProducts(); i++){
			Mj = Mj + network->getArrivalRateAtRetailer(i, j);
		}
		std::cout << "demand at j=" << j << " equals: " << Mj << std::endl;
		std::cout << "EBOj for j=" << j << " equals: " << EBOj[j - 1] << std::endl;
		std::cout << "beta star for j=" << j << " equals: " << (1 - (EBOj[j - 1] / Mj)) << std::endl;
	} // for

	 // ------------------------------------------------------------------------------------------------------- write to excel --
	


	/*
	BasicExcel e;
	e.New(1);

	BasicExcelWorksheet* sheet = e.GetWorksheet("Sheet1");
	BasicExcelCell* cell;

	if (sheet) {

		int row = 0;

		for (int j = 0; j <= network->sizeRetailers(); j++){

			for (int i = 1; i <= network->sizeProducts(); i++) {
				
				cell = sheet->Cell(row, 0); // location
				cell->Set((int)j);

				cell = sheet->Cell(row, 1); // product
				cell->Set((int)i);

				cell = sheet->Cell(row, 2); // base-stock level
				
				if (j == 0) {
					cell->Set((double)network->getBaseStockLevelAtWarehouse(i));
				} else {
					cell->Set((double)network->getBaseStockLevelAtRetailer(i, j));
				} // eif

				row = row + 1;

			} // for

		} // for

	} // if
	
	e.SaveAs("base-stock-levels.xls");

	std::cout << "wrote results to Excel" << std::endl;
*/

	network->writeBaseStockLevelsToFile("base-stock-levels.txt");

	return a.exec();
	
	
} // main
