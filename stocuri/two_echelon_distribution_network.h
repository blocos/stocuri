#ifndef TWO_ECHELON_DISTRIBUTION_NETWORK
#define TWO_ECHELON_DISTRIBUTION_NETWORK

#include "supply_chain_network.h"

#include <QtCore>

#include <cassert>
#include <fstream>
#include <iostream>

#include "BasicExcel.hpp"

using namespace YExcel;

class TwoEchelonDistributionNetwork : SupplyChainNetwork {

	private:
		QList<QList<double>*> *arrivalRates;
		QList<QList<double>*> *leadTimes;
		QList<QList<double>*> *baseStockLevels;
		QList<QList<double>*> *inventoryHoldingCosts;

		QVector<QString> *productsIdToNumber;

		int nRetailers;
		int nProducts;

	public:
		TwoEchelonDistributionNetwork();
		TwoEchelonDistributionNetwork(int nRetailers, int nProducts);
		~TwoEchelonDistributionNetwork();

		// warehouse related 
		void setArrivalRateAtWarehouse(int product, double arrivalRate);
		double getArrivalRateAtWarehouse(int product);

		void setLeadTimeToWarehouse(int product, double leadTime);
		double getLeadTimeToWarehouse(int product);

		void setBaseStockLevelAtWarehouse(int product, double baseStockLevel);
		double getBaseStockLevelAtWarehouse(int product);

		void setInventoryHoldingCostAtWarehouse(int product, double inventoryHoldingCost);
		double getInventoryHoldingCostAtWarehouse(int product);

		// retailer related
		void setArrivalRateAtRetailer(int product, int retailer, double arrivalRate);
		double getArrivalRateAtRetailer(int product, int retailer);

		void setLeadTimeToRetailer(int product, int retailer, double leadTime);
		double getLeadTimeToRetailer(int product, int retailer);

		void setBaseStockLevelAtRetailer(int product, int retailer, double baseStockLevel);
		double getBaseStockLevelAtRetailer(int product, int retailer);

		void setInventoryHoldingCostAtRetailer(int product, int retailer, double inventoryHoldingCost);
		double getInventoryHoldingCostAtRetailer(int product, int retailer);

		// general
		int sizeRetailers();
		int sizeProducts();

		// overriding virtual methods

		int loadFromFile(QString fileNameSettings, QString filenNameDemand);
		int writeBaseStockLevelsToFile(QString fileName);
		int writeBaseStockLevelsToExcel(QString fileName);

};

#endif

