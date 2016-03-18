#include "two_echelon_distribution_network.h"

// ---------------------------------------------------------------------------------------------- CONSTRUCTORS AND DESTRUCTORS --

TwoEchelonDistributionNetwork::TwoEchelonDistributionNetwork() {
}


TwoEchelonDistributionNetwork::TwoEchelonDistributionNetwork(int nRetailers, int nProducts) {
	// pre	: nWarehouse > 0 /\ nProducts > 0
	// post	: arrivalRates = leadTimes = baseStockLevels = inventoryHoldingCosts = [ [ 0, ..., 0 ], ..., [ 0, ..., 0 ] ]

	assert(nRetailers > 0);
	assert(nProducts > 0);

	// pre-conditions satisfied

	this->nRetailers = nRetailers;
	this->nProducts = nProducts;

	arrivalRates = new QList<QList<double>*>();
	leadTimes = new QList<QList<double>*>();
	baseStockLevels = new QList<QList<double>*>();
	inventoryHoldingCosts = new QList<QList<double>*>();

	// retailers + 1 warehouse
	for (int j = 0; j < nRetailers + 1; j++){

		QList<double>* arrivalRatesAtJ = new QList<double>();
		QList<double>* leadTimesAtJ = new QList<double>();
		QList<double>* baseStockLevelsAtJ = new QList<double>();
		QList<double>* inventoryHoldingCostAtJ = new QList<double>();

		for (int i = 0; i < nProducts; i++){
			arrivalRatesAtJ->append(0);
			leadTimesAtJ->append(0);
			baseStockLevelsAtJ->append(0);
			inventoryHoldingCostAtJ->append(0);
		} // for

		arrivalRates->append(arrivalRatesAtJ);
		leadTimes->append(leadTimesAtJ);
		baseStockLevels->append(baseStockLevelsAtJ);
		inventoryHoldingCosts->append(inventoryHoldingCostAtJ);

	} // for

} // TwoEchelonDistributionNetwork


TwoEchelonDistributionNetwork::~TwoEchelonDistributionNetwork() {

}

// ------------------------------------------------------------------------------------------------------- GETTERS AND SETTERS --

void TwoEchelonDistributionNetwork::setArrivalRateAtWarehouse(int product, double arrivalRate) {
	(*(*arrivalRates)[0])[product - 1] = arrivalRate;
}

double TwoEchelonDistributionNetwork::getArrivalRateAtWarehouse(int product) {
	return arrivalRates->at(0)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setLeadTimeToWarehouse(int product, double leadTime) {
	(*(*leadTimes)[0])[product - 1] = leadTime;
}

double TwoEchelonDistributionNetwork::getLeadTimeToWarehouse(int product) {
	return leadTimes->at(0)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtWarehouse(int product, double baseStockLevel) {
	(*(*baseStockLevels)[0])[product - 1] = baseStockLevel;
}

double TwoEchelonDistributionNetwork::getBaseStockLevelAtWarehouse(int product) {
	return baseStockLevels->at(0)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setInventoryHoldingCostAtWarehouse(int product, double inventoryHoldingCost){
	(*(*inventoryHoldingCosts)[0])[product - 1] = inventoryHoldingCost;
}

double TwoEchelonDistributionNetwork::getInventoryHoldingCostAtWarehouse(int product) {
	return inventoryHoldingCosts->at(0)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setArrivalRateAtRetailer(int product, int retailer, double arrivalRate) {
	(*(*arrivalRates)[retailer])[product - 1] = arrivalRate;
}

double TwoEchelonDistributionNetwork::getArrivalRateAtRetailer(int product, int retailer) {
	return arrivalRates->at(retailer)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setLeadTimeToRetailer(int product, int retailer, double leadTime) {
	(*(*leadTimes)[retailer])[product - 1] = leadTime;
}

double TwoEchelonDistributionNetwork::getLeadTimeToRetailer(int product, int retailer) {
	return leadTimes->at(retailer)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setBaseStockLevelAtRetailer(int product, int retailer, double baseStockLevel) {
	(*(*baseStockLevels)[retailer])[product - 1] = baseStockLevel;
}

double TwoEchelonDistributionNetwork::getBaseStockLevelAtRetailer(int product, int retailer) {
	return baseStockLevels->at(retailer)->at(product - 1);
}

void TwoEchelonDistributionNetwork::setInventoryHoldingCostAtRetailer(int product, int retailer, double inventoryHoldingCost) {
	(*(*inventoryHoldingCosts)[retailer])[product - 1] = inventoryHoldingCost;
}

double TwoEchelonDistributionNetwork::getInventoryHoldingCostAtRetailer(int product, int retailer) {
	return inventoryHoldingCosts->at(retailer)->at(product - 1);
}


// ------------------------------------------------------------------------------------------------------------ PUBLIC METHODS --

int TwoEchelonDistributionNetwork::sizeProducts() {
	return nProducts;
}

int TwoEchelonDistributionNetwork::sizeRetailers() {
	return nRetailers;
}