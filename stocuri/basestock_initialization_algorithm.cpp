#include "basestock_initialization_algorithm.h"


BasestockInitializationAlgorithm::BasestockInitializationAlgorithm()
{
}


BasestockInitializationAlgorithm::~BasestockInitializationAlgorithm()
{
}

double BasestockInitializationAlgorithm::pPartsOnOrderAtRetailer(int product, int retailer, double x) {
	// pre	: ( 1 <= product /\ product <= |I| ) /\ ( 1 <= retailer /\ retailer <= |J| ) /\ x >= 0
	// ret	: P[X = x], X~Poisson(lambda) with lambda = arrival rate * lead time for item

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeProducts());
	assert(x >= 0);

	PoissonDistribution tetrodotoxin;

	double result = 0.0;

	double mij = network->getArrivalRateAtRetailer(product, retailer);
	double Lij = network->getLeadTimeToRetailer(product, retailer);

	double lambda = mij * Lij;

	try{
		result = tetrodotoxin.probabilityBySterlingApproximation(lambda, x);
	}
	catch (std::exception& me) {
		std::cout << me.what() << std::endl;

		result = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
	} // catch me if you can

	return result;

} // pPartsOnOrderAtRetailer


double BasestockInitializationAlgorithm::calculateBetaJ(double BigM, int retailer) {

	double result = 0.0;

	PoissonDistribution tetrodotoxin;

	for (int i = 1; i <= network->sizeProducts(); i++){
		double sum = 0.0;

		double mij = network->getArrivalRateAtRetailer(i, retailer);

		for (int x = 0; x <= network->getBaseStockLevelAtRetailer(i, retailer); x++){
			double pxx = pPartsOnOrderAtRetailer(i, retailer, x);

			sum = sum + pxx;
		} // for

		result = result + ((mij / BigM)*sum);

	} // for

	return result;
}

int BasestockInitializationAlgorithm::run(TwoEchelonDistributionNetwork *network, QList<double> *targetAggregateFillRates) {
	int result = 0;

	this->network = network;

	

	for (int j = 1; j <= network->sizeRetailers(); j++) {
		double BigM = 0;
		for (int i = 1; i <= network->sizeProducts(); i++) {
			BigM = BigM + network->getArrivalRateAtRetailer(i, j);

			double mij = network->getArrivalRateAtRetailer(i, j);
			double Lij = network->getLeadTimeToRetailer(i, j);

			double f = ceil(mij*Lij);

			double max = 0.0;

			if (f > max) {
				max = f;
			}

			network->setBaseStockLevelAtRetailer(i, j, max);
		} // for


		// 
		double betaJ = calculateBetaJ(BigM, j);

		if (betaJ < (*targetAggregateFillRates)[j - 1]) {
			
			QVector<double> delta = QVector<double>(network->sizeProducts(), 0);

			while (true) {
				
				double argmax = -1;
				double max = 0;

				for (int i = 1; i <= network->sizeProducts(); i++){
					double mij = network->getArrivalRateAtRetailer(i, j);
					double hij = network->getInventoryHoldingCostAtRetailer(i, j);
					double Sij = network->getBaseStockLevelAtRetailer(i, j);

					double pxx = pPartsOnOrderAtRetailer(i, j, Sij);

					delta[i - 1] = (mij*pxx) / (BigM * hij);

					if (max < delta[i - 1]) {
						max = delta[i - 1];
						argmax = i;
					}
				} // for

				double Skj = network->getBaseStockLevelAtRetailer(argmax, j);

				Skj = Skj + 1;

				network->setBaseStockLevelAtRetailer(argmax, j, Skj);

				betaJ = calculateBetaJ(BigM, j);

				if (betaJ >= (*targetAggregateFillRates)[j - 1]) {
					break;
				} // if
			} // while
		} // if
	} // for


	return result;
}
