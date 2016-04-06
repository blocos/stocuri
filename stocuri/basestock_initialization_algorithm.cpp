#include "basestock_initialization_algorithm.h"


BasestockInitializationAlgorithm::BasestockInitializationAlgorithm()
{
}


BasestockInitializationAlgorithm::~BasestockInitializationAlgorithm()
{
}


double BasestockInitializationAlgorithm::calculateBetaJ(double BigM, int retailer) {

	double result = 0.0;

	PoissonDistribution tetrodotoxin;

	for (int i = 1; i <= network->sizeProducts(); i++){
		double sum = 0.0;

		double mij = network->getArrivalRateAtRetailer(i, retailer);
		double Lij = network->getLeadTimeToRetailer(i, retailer);

		double lambda = mij * Lij;

		for (int x = 0; x <= network->getBaseStockLevelAtRetailer(i, retailer); x++){
			double pxx = 0.0;
			try{
				pxx = tetrodotoxin.probability(lambda, x);
			}
			catch (std::exception& me) {
				std::cout << me.what() << std::endl;

				pxx = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
			} // catch me if you can

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
		} // for


		// 
		double betaJ = calculateBetaJ(BigM, j);

		

		std::cout << "beta for j=" << j << ": " << betaJ << std::endl;



	} // for


	return result;
}
