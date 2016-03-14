#include "two_echelon_distribution_network.h"

// ---------------------------------------------------------------------------------------------- CONSTRUCTORS AND DESTRUCTORS --

TwoEchelonDistributionNetwork::TwoEchelonDistributionNetwork() {
}


TwoEchelonDistributionNetwork::TwoEchelonDistributionNetwork(int nWarehouses, int nProducts) {
	// pre	: nWarehouse > 0 && nProducts > 0
	// post	: 

	assert(nWarehouses > 0);
	assert(nProducts > 0);

	// pre-conditions satisfied

	QList<QList<double>*> *arrivalRates = new QList<QList<double>*>();

	for (int j = 0; j < nWarehouses; j++){

		QList<double>* products = new QList<double>();

		for (int i = 0; i < nProducts; i++){
			products->append(0);
		} // for

		arrivalRates->append(products);

	} // for

} // TwoEchelonDistributionNetwork


TwoEchelonDistributionNetwork::~TwoEchelonDistributionNetwork() {
}

// ------------------------------------------------------------------------------------------------------- GETTERS AND SETTERS --

void TwoEchelonDistributionNetwork::setArrivalRateAtWarehouse(int product){
}

void TwoEchelonDistributionNetwork::getArrivalRateAtWarehouse(int product){
}

void TwoEchelonDistributionNetwork::setLeadTimeToWarehouse(int product){
}

void TwoEchelonDistributionNetwork::getLeadTimeToWarehouse(int product){
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtWarehouse(int product){
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtWarehouse(int product){
}

void TwoEchelonDistributionNetwork::setArrivalRateAtWarehouse(int retailer, int product){
}

void TwoEchelonDistributionNetwork::getArrivalRateAtWarehouse(int retailer, int product){
}

void TwoEchelonDistributionNetwork::setLeadTimeToRetailer(int retailer, int product){
}

void TwoEchelonDistributionNetwork::getLeadTimeToWarehouse(int retailer, int product){
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtRetailer(int retailer, int product){
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtWarehouse(int retailer, int product){
}