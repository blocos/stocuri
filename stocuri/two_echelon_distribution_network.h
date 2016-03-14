#ifndef TWO_ECHELON_DISTRIBUTION_NETWORK
#define TWO_ECHELON_DISTRIBUTION_NETWORK

#include "supply_chain_network.h"

#include <QtCore>

#include <cassert>

class TwoEchelonDistributionNetwork : SupplyChainNetwork {

	/*private:
		QList<QList<double>*> *m;
		QList<QList<double>*> *t;
		QList<QList<double>*> *S;*/

	private:
		QList< QList< double > *> *arrivalRates;

	public:
		TwoEchelonDistributionNetwork();
		TwoEchelonDistributionNetwork(int nRetailers, int nProducts);
		~TwoEchelonDistributionNetwork();

		// warehouse related 
		void setArrivalRateAtWarehouse(int product);
		void getArrivalRateAtWarehouse(int product);

		void setLeadTimeToWarehouse(int product);
		void getLeadTimeToWarehouse(int product);

		void setBaseStockLevelAtWarehouse(int product);
		void setBaseStockLevelAtWarehouse(int product);

		// retailer related
		void setArrivalRateAtWarehouse(int retailer, int product);
		void getArrivalRateAtWarehouse(int retailer, int product);

		void setLeadTimeToRetailer(int retailer, int product);
		void getLeadTimeToWarehouse(int retailer, int product);

		void setBaseStockLevelAtRetailer(int retailer, int product);
		void setBaseStockLevelAtWarehouse(int retailer, int product);


		void sizeRetailers();
		void sizeProducts();


};

#endif

