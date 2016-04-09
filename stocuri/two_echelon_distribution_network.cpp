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

	productsIdToNumber = new QVector<QString>(nProducts, "");

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
	for (int j = 0; j < nRetailers + 1; j++){
		delete (*arrivalRates)[j];
		delete (*leadTimes)[j];
		delete (*baseStockLevels)[j];
		delete (*inventoryHoldingCosts)[j];
	}

	delete arrivalRates;
	arrivalRates = 0;

	delete leadTimes;
	leadTimes = 0;

	delete baseStockLevels;
	baseStockLevels = 0;

	delete inventoryHoldingCosts;
	inventoryHoldingCosts = 0;

	delete productsIdToNumber;

} // ~TwoEchelonDistributionNetwork


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


// -------------------------------------------------------------------------------------------------------- OVERRIDDEN METHODS --

int TwoEchelonDistributionNetwork::loadFromFile(QString fileNameSettings, QString filenNameDemand){
	int result = 0;

	bool sfl = true;
	int lineCounter = 0;

	// parse settings -----------------------------------------------------------------------------------------------------------

	std::ifstream fileInputSettings(fileNameSettings.toStdString());
	
	std::string line;
	while (std::getline(fileInputSettings, line, '\n')) {
		if (sfl) {
			sfl = false;
			continue;
		} 

		// convert
		QString qLine = QString::fromStdString(line);

		// split in parts
		QStringList parts = qLine.split(",");

		// id (0), number (1), price (2), H-RSHQ (3), H-CSC (4), L-CSC (5), L-USA (6), L-Singapore (7), L-China (8), L-USA+1 (9), L-Singapore+1 (10), L-China+1 (11)

		// retrieve id
		int id = parts[0].toInt();

		// got one additional line
		lineCounter = lineCounter + 1;

		// store relation id->number
		(*productsIdToNumber)[id-1] = parts[1];

		// store lead times
		double leadTimeCSC = parts[5].toDouble();
		this->setLeadTimeToWarehouse(id, leadTimeCSC);

		double leadTimeUSA = parts[9].toDouble();
		this->setLeadTimeToRetailer(id, 1, leadTimeUSA);

		double leadTimeSingapore = parts[9].toDouble();
		this->setLeadTimeToRetailer(id, 2, leadTimeSingapore);

		double leadTimeChina = parts[9].toDouble();
		this->setLeadTimeToRetailer(id, 3, leadTimeChina);

		// store inventory holding costs
		double inventoryHoldingCostCSC = parts[4].toDouble();
		this->setInventoryHoldingCostAtWarehouse(id, inventoryHoldingCostCSC);

		double inventoryHoldingCostRSHQ = parts[3].toDouble();
		this->setInventoryHoldingCostAtRetailer(id, 1, inventoryHoldingCostRSHQ);
		this->setInventoryHoldingCostAtRetailer(id, 2, inventoryHoldingCostRSHQ);
		this->setInventoryHoldingCostAtRetailer(id, 3, inventoryHoldingCostRSHQ);

		if (lineCounter == nProducts) {
			break;
		}

	} // while

	fileInputSettings.close();

	qDebug() << *productsIdToNumber;

	// parse demand -------------------------------------------------------------------------------------------------------------
	
	sfl = true;
	lineCounter = 0;

	std::ifstream fileInputDemand(filenNameDemand.toStdString());

	while (std::getline(fileInputDemand, line, '\n')) {
		if (sfl) {
			sfl = false;
			continue;
		} // if

		// convert
		QString qLine = QString::fromStdString(line);

		// read one additional line
		lineCounter = lineCounter + 1;

		// split in parts
		QStringList parts = qLine.split(",");

		// region (0), number (1), rate (2)

		// retrieve name
		QString name = parts[0];
		QString number = parts[1];
		
		double rate = parts[2].toDouble();

		int id = productsIdToNumber->indexOf(number.trimmed()) + 1;

		if (name == "USA") {
			this->setArrivalRateAtRetailer(id, 1, rate);
		} else if (name == "Singapore") {
			this->setArrivalRateAtRetailer(id, 2, rate);
		} else if (name == "China"){
			this->setArrivalRateAtRetailer(id, 3, rate);
		} else if (name == "World"){
			this->setArrivalRateAtWarehouse(id, rate);
		} else {
			// undefined
		} // eif
		
		if ((lineCounter / 4) == nProducts) {
			break;
		} // if

	} // while

	fileInputDemand.close();

	return result;

} // loadFromFile


int TwoEchelonDistributionNetwork::writeBaseStockLevelsToFile(QString fileName){
	bool result = 0;

	std::ofstream fileOutput(fileName.toStdString(), std::ios::out);

	QVector<QString> regions = QVector<QString>() << "CSC" << "USA" << "Singapore" << "China";

	if (!fileOutput.is_open()) {
		result = 1;
	} else {
		for (int j = 0; j <= sizeRetailers(); j++) {
			for (int i = 1; i <= sizeProducts(); i++) {
				double SiX = 0.0;

				if (j == 0){
					SiX = getBaseStockLevelAtWarehouse(i);
				} else{
					SiX = getBaseStockLevelAtRetailer(i, j);
				} // eif

				QString number = productsIdToNumber->at(i - 1);
				fileOutput << regions[j].toStdString() << ", " << number.toStdString() << ", " << SiX << std::endl;
			} // for
		} // for

		fileOutput.close();

	} // eif

	return result;

} // writeBaseStockLevelsToFile