#include "greedy_algorithm.h"

// ---------------------------------------------------------------------------------------------- CONSTRUCTORS AND DESTRUCTORS --

GreedyAlgorithm::GreedyAlgorithm() {
}

GreedyAlgorithm::~GreedyAlgorithm() {
}

// ------------------------------------------------------------------------------------------------------- GETTERS AND SETTERS --

// none


// ---------------------------------------------------------------------------------------------- PRIVATE METHODS -- warehouse --

double GreedyAlgorithm::pPartsOnOrderAtWarehouse(int product, int x) {
	// pre	: 1 <= product /\ product <= |I| /\ x >= 0
	// ret	: P[X = x], X~Poisson(lambda) with lambda = arrival rate * lead time for item

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	// determine mean for Poisson distribution
	double lambda = network->getArrivalRateAtWarehouse(product) * network->getLeadTimeToWarehouse(product);

	// initialize Poisson distribution
	PoissonDistribution tetrodotoxin;

	try {
		// try regular calculation
		//result = tetrodotoxin.probabilityBySterlingApproximation(lambda, x);
		result = tetrodotoxin.probability(lambda, x);
	} catch (std::exception& me) {
		std::cout << me.what() << std::endl;

		// catched me, try calculation with Normal approximation
		result = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
	} // catch-me-if-you-can

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnOrderAtWarehouse


double GreedyAlgorithm::pPartsOnHandAtWarehouse(int product, int x) {
	// pre	: 1 <= product <= |I| /\ x >= 0
	// ret	: 0

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnHandAtWarehouse


double GreedyAlgorithm::pPartsOnBackorderAtWarehouse(int product, int x) {
	// pre	: 1 <= product <= |I| /\ x >= 0
	// ret	: P[Bi0(Si0) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0;

	double Si0 = network->getBaseStockLevelAtWarehouse(product);

	if (x == 0) {
		double sigma = 0;
		for (int y = 0; y <= Si0; y++){
			sigma = sigma + pPartsOnOrderAtWarehouse(product, y);
		} // for
		result = sigma;
	} else {
		return pPartsOnOrderAtWarehouse(product, Si0 + x);
	} // eif

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnBackorderAtWarehouse


double GreedyAlgorithm::ePartsOnHandAtWarehouse(int product) {
	// pre	: 1 <= product <= |I|
	// ret	: E[OHi0(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());

	// pre-conditions satisfied

	double result = 0.0;

	double Si0 = network->getBaseStockLevelAtWarehouse(product);

	for (int x = 0; x <= Si0; x++){
		result = result + ((Si0 - (double)x)*pPartsOnOrderAtWarehouse(product, x));
	} // for

	assert(result >= 0);

	return result;

} // ePartsOnHandAtWarehouse


double GreedyAlgorithm::ePartsOnBackorderAtWarehouse(int product) {
	// pre	: 1 <= product <= |I|
	// ret	: E[BOi0(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());

	// pre-conditions satisfied

	double result = 0.0;

	double ai0 = network->getArrivalRateAtWarehouse(product);
	double Li0 = network->getLeadTimeToWarehouse(product);
	double Si0 = network->getBaseStockLevelAtWarehouse(product);

	result = (ai0 * Li0) - Si0 + ePartsOnHandAtWarehouse(product);

	assert(result >= 0);

	return result;

} // ePartsOnBackorderAtWarehouse


double GreedyAlgorithm::pPartsOnBackorderAtWarehouseFromRetailer(int product, int retailer, int x) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: P[BOi0j(Si0) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());
	assert(x >= 0);
	
	// pre-conditions satisfied

	double result = 0;

	double mi0 = network->getArrivalRateAtWarehouse(product);
	double mij = network->getArrivalRateAtRetailer(product, retailer);
	
	double previous = 0;

	int y = x;

	while (true) {

		// y over x, expressed in the logarithmic gamma function, translate using exponential
		double binCoef = exp(lgammal(y + 1) - lgammal(x + 1) - lgammal(y - x + 1));
		double powerX = pow(mij / mi0, x);
		double powerYminX = pow(1.0 - (mij / mi0), y - x);
		double pBO = pPartsOnBackorderAtWarehouse(product, y);

		double value = binCoef * powerX * powerYminX * pBO;

		result = result + value;

		y = y + 1;

		if (result - previous < BOCUTOFF) {
			break;
		} // if

		previous = result;

	} // while

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnBackorderAtWarehouseFromRetailer


double GreedyAlgorithm::ePartsOnBackorderAtWarehouseFromRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: E[BOi0j(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mi0 = network->getArrivalRateAtWarehouse(product);
	double mij = network->getArrivalRateAtRetailer(product, retailer);

	result = (mij / mi0)*ePartsOnBackorderAtWarehouse(product);

	assert(result >= 0);
	
	return result;

} // ePartsOnBackorderAtWarehouseFromRetailer


// ----------------------------------------------------------------------------------------------- PRIVATE METHODS -- retailer --

double GreedyAlgorithm::pPartsOnOrderAtRetailer(int product, int retailer, int x) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: P[Xij(Si0) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	for (int k = 0; k <= x; k++) {
		// P[Y=k]*P[BOj=x-k]

		// 1. P[Y=k]

		double mij = network->getArrivalRateAtRetailer(product, retailer);
		double Lij = network->getLeadTimeToRetailer(product, retailer);
		double lambda = mij * Lij;
		double p1 = 0.0;

		PoissonDistribution tetrodotoxin;

		try {
			// try regular calculation
			p1 = tetrodotoxin.probability(lambda, k);
		} catch (std::exception& me) {
			std::cout << me.what() << std::endl;

			// catched me, try calculation with Normal approximation
			p1 = tetrodotoxin.probabilityByNormalApproximation(lambda, k);
		} // try to catch me

		// 2. P[BOj = x - k]
		double p2 = pPartsOnBackorderAtWarehouseFromRetailer(product, retailer, x-k);

		result = result + (p1 * p2);

	} // while

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnOrderAtRetailer


double GreedyAlgorithm::vPartsOnBackorderAtWarehouse(int product) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: Var[BOi0j(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());

	// pre-conditions satisfied

	double result = 0.0;
	
	double EBO = ePartsOnBackorderAtWarehouse(product);
	
	double Si0 = network->getBaseStockLevelAtWarehouse(product);
	double mi0 = network->getArrivalRateAtWarehouse(product);
	double Li0 = network->getLeadTimeToWarehouse(product);

	double EX = mi0*Li0;
	double EX2 = mi0*Li0*((mi0*Li0) + 1);

	double sum = 0.0;
	for (int x = 0; x <= Si0; x++){
		sum = sum + (((Si0 - (double)x)*(Si0 - (double)x))*pPartsOnOrderAtWarehouse(product, x));
	} // for

	double EBO2 = EX2 - (2 * Si0 * EX) + (Si0*Si0) - sum;

	result = EBO2 - (EBO * EBO);

	return result;

} // vPartsOnBackorderAtWarehouseFromRetailer


double GreedyAlgorithm::vPartsOnBackorderAtWarehouseFromRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: Var[BOi0j(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mi0 = network->getArrivalRateAtWarehouse(product);
	double mij = network->getArrivalRateAtRetailer(product, retailer);
	double fij = mij / mi0;

	double EBO = ePartsOnBackorderAtWarehouse(product);
	double VarBO = vPartsOnBackorderAtWarehouse(product);

	result = (fij * fij * VarBO) + (fij * (1 - fij) * EBO);

	return result;

} // vPartsOnBackorderAtWarehouseFromRetailer


double GreedyAlgorithm::ePartsOnOrderAtRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: E[Xij(Si0,Sij)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mij = network->getArrivalRateAtRetailer(product, retailer);
	double Lij = network->getLeadTimeToRetailer(product, retailer);
	double EBOj = ePartsOnBackorderAtWarehouseFromRetailer(product, retailer);

	result = mij*Lij + EBOj;

	return result;

} // ePartsOnOrderAtRetailer


double GreedyAlgorithm::vPartsOnOrderAtRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: Var[Xij(Si0,Sij)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mij = network->getArrivalRateAtRetailer(product, retailer);
	double Lij = network->getLeadTimeToRetailer(product, retailer);
	double VarBOj = vPartsOnBackorderAtWarehouseFromRetailer(product, retailer);

	result = mij*Lij + VarBOj;

	return result;
} // vPartsOnOrderAtRetailer


double GreedyAlgorithm::pPartsOnHandAtRetailer(int product, int retailer, int x){
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: P[OHij(Si0, Sij) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnHandAtRetailer


double GreedyAlgorithm::pPartsOnBackorderAtRetailer(int product, int retailer, int x) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: P[BOij(Si0, Sij) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

	double Sij = network->getBaseStockLevelAtRetailer(product, retailer);

	if (x == 0){
		for (int y = 0; y <= Sij; y++) {
			result = result + pPartsOnOrderAtRetailer(product, retailer, y);
		} // for
	} else{
		result = pPartsOnOrderAtRetailer(product, retailer, Sij + x);
	} // eif

	assert(result >= 0);
	assert(result <= 1);

	return result;

} // pPartsOnBackorderAtRetailer


double GreedyAlgorithm::ePartsOnHandAtRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc|
	// ret	: E[BOij(Si0, Sij)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double Sij = network->getBaseStockLevelAtRetailer(product, retailer);

	for (int x = 0; x <= Sij; x++){
		if (GRAVES) {
			result = result + ((Sij - x)*pPartsOnOrderAtRetailer2Moment(product, retailer, x));

			if (debug){
				std::cout << "t n2m: " << pPartsOnOrderAtRetailer(product, retailer, x) << std::endl;
				std::cout << "t 2m: " << pPartsOnOrderAtRetailer2Moment(product, retailer, x) << std::endl;
			}
		}
		else {
			result = result + ((Sij - x)*pPartsOnOrderAtRetailer(product, retailer, x));
		}
		
	} // for

	assert(result >= 0);

	return result;

} // ePartsOnHandAtRetailer


double GreedyAlgorithm::ePartsOnBackorderAtRetailer(int product, int retailer){
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc|
	// ret	: E[BOij(Si0, Sij)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mij = network->getArrivalRateAtRetailer(product, retailer);
	double Lij = network->getLeadTimeToRetailer(product, retailer);
	double Sij = network->getBaseStockLevelAtRetailer(product, retailer);

	result = (mij*Lij) - Sij + ePartsOnBackorderAtWarehouseFromRetailer(product, retailer) + ePartsOnHandAtRetailer(product, retailer);	

	assert(result >= 0);	

	return result;

} // ePartsOnBackorderAtWarehouse


// ------------------------------------------------------------------------------------------------------------ PUBLIC METHODS --

QList<double> GreedyAlgorithm::evaluateNetwork(TwoEchelonDistributionNetwork *network) {
	// pre	: True
	// post	: [EBO1, EBO2, ..., EBOj], 1 <= j <= |Jloc|

	// pre-conditions satisfied

	this->network = network;

	std::cout << "VarBOi0: " << vPartsOnBackorderAtWarehouse(1) << std::endl;

	std::cout << "EXij: " << ePartsOnOrderAtRetailer(1, 1) << std::endl;
	std::cout << "VarXij: " << vPartsOnOrderAtRetailer(1, 1) << std::endl;

	double mi0 = network->getArrivalRateAtWarehouse(1);
	double Li0 = network->getLeadTimeToWarehouse(1);

	double EX2 = vPartsOnOrderAtRetailer(1, 1) + ePartsOnOrderAtRetailer(1,1)*ePartsOnOrderAtRetailer(1,1);

	std::cout << "EXij squared: " << EX2 << std::endl;

	double p = (vPartsOnOrderAtRetailer(1, 1) - ePartsOnOrderAtRetailer(1, 1)) / vPartsOnOrderAtRetailer(1, 1);
	std::cout << "p: " << p << std::endl;

	double k = ((1 - p) / p) *ePartsOnOrderAtRetailer(1, 1);

	std::cout << "k: " << k << std::endl;

	std::cout << "------------------------" << std::endl;



	std::cout << "EBOj: " << ePartsOnBackorderAtRetailer(1, 1) << std::endl;


	GRAVES = false;

	std::cout << "EBOj 2: " << ePartsOnBackorderAtRetailer(1, 1) << std::endl;
	

	for (int x = 0; x <= 6; x++) {
		//std::cout << pPartsOnBackorderAtWarehouse(1, x) << std::endl;
	}

	for (int x = 0; x <= 3; x++) {
		//std::cout << pPartsOnBackorderAtRetailer(1, 1, x) << std::endl;
	}

	QList<double> EBOj = QList<double>();

	for (int j = 1; j <= network->sizeRetailers(); j++) {
		EBOj.append(0);


		for (int i = 1; i <= network->sizeProducts(); i++) {

			double eboj = ePartsOnBackorderAtRetailer(i, j);

			EBOj[j - 1] = EBOj[j - 1] + eboj;

		
			std::cout << "product: " << i << ", location: " << j << ": E[BOj]=" << eboj << std::endl;
		
		} // for
	} // for

	return EBOj;

} // evaluateNetwork

int GreedyAlgorithm::optimizeNetwork(TwoEchelonDistributionNetwork *network, QList<double> *targetAggregateFillRates) {
	// pre	: network->sizeRetailers = |targetAggregateFillRates|
	// post	: []

	assert(network->sizeRetailers() == targetAggregateFillRates->size());

	// pre-conditions satisfied

	this->network = network;

	// Lines 2 - 3
	// implicitly expected from argument 'network'

	// Lines 4 - 5

	QList<double> *EBOj = new QList<double>();
	
	for (int j = 1; j <= network->sizeRetailers(); j++){

		EBOj->append(0);

		for (int i = 1; i <= network->sizeProducts(); i++) {

			double mij = network->getArrivalRateAtRetailer(i, j);
			double Li0 = network->getLeadTimeToWarehouse(i);
			double Lij = network->getLeadTimeToRetailer(i, j);

			(*EBOj)[j - 1] = (*EBOj)[j - 1] + (mij*(Li0 + Lij));

		} // for

	} // for
	

	// Line 6

	double costS = 0;

	// Line 7 

	if (!isTargetFillRatesSatisfied(network, EBOj, targetAggregateFillRates)) {
		
		// Line 8
		
		while ( true ) {

			// Lines 9 - 10

			QVector<QVector<double>> gamma = QVector<QVector<double>>(network->sizeProducts(), QVector<double>(network->sizeRetailers() + 1, 0));

			int k = -1;
			int l = -1;

			double deltaMax = -1000.0;

			

			for (int i = 1; i <= network->sizeProducts(); i++) {

				// for each location, including the warehouse
				for (int j = 0; j <= network->sizeRetailers(); j++) {

					double deltaEBO = calculateDeltaEBO(i, j, targetAggregateFillRates);
				
					double hiX = 0.0;

					if (j == 0) {
						hiX = network->getInventoryHoldingCostAtWarehouse(i);
					} else{
						hiX = network->getInventoryHoldingCostAtRetailer(i, j);
					} // eif

					gamma[i - 1][j] = deltaEBO / hiX;

					// Line 11
					if (deltaMax < gamma[i - 1][j]) {
						k = i;
						l = j; 
						deltaMax = gamma[i - 1][j];
					} // if

				} // for

			} //

			// Line 12

			double SiX = 0.0;

			if (l == 0) {
				SiX = network->getBaseStockLevelAtWarehouse(k);
				network->setBaseStockLevelAtWarehouse(k, SiX + 1);
			} else {
				SiX = network->getBaseStockLevelAtRetailer(k, l);
				network->setBaseStockLevelAtRetailer(k, l, SiX + 1);
			} // eif

			// Line 13

			for (int j = 1; j <= network->sizeRetailers(); j++){
				(*EBOj)[j - 1] = 0;
				for (int i = 1; i <= network->sizeProducts(); i++) {
					(*EBOj)[j - 1] = (*EBOj)[j - 1] + ePartsOnBackorderAtRetailer(i, j);
				} // for
			} // for

			if (isTargetFillRatesSatisfied(network, EBOj, targetAggregateFillRates)) {
				break;
			} // if
		} // while
	} // if


	std::cout << "n: " << n << std::endl;
	std::cout << "w: " << w << std::endl;

	// Line 14;

	return 0;

} // optimizeNetwork


// ------------------------------------------------------------------------------------------------ PRIVATE METHODS -- support --

bool GreedyAlgorithm::isTargetFillRatesSatisfied(TwoEchelonDistributionNetwork *network, QList<double> *EBOj, QList<double> *targetAggregateFillRates){
	// pre	: True
	// ret	: not Line 7 or Line 13

	bool result = true;

	for (int j = 1; j <= network->sizeRetailers(); j++) {

		double arrivalRateOverI = 0.0;
		for (int i = 1; i <= network->sizeProducts(); i++) {
			arrivalRateOverI = arrivalRateOverI + network->getArrivalRateAtRetailer(i, j);
		} // for

		double actualFillRate = (1 - ((*EBOj)[j - 1] / arrivalRateOverI));

		if (actualFillRate < targetAggregateFillRates->at(j - 1)) {
			result = false;
		} // if

	} // for 

	return result;
}


unsigned long long GreedyAlgorithm::binomialCoefficient(unsigned long long n, unsigned long long k) {
	unsigned long long c = 1, i;

	if (k > n - k) {
		k = n - k;  // take advantage of symmetry
	} // if

	for (i = 1; i <= k; i++, n--) {
		if (c / i > UINT_MAX / n) {
			return 0;  // return 0 on overflow
			assert(false);
		} // if
		c = c / i * n + c%i * n / i;  // split c*n/i into (c/i*i + c%i)*n/i
	} // for

	return c;

} // binomialCoefficient


double GreedyAlgorithm::calculateDeltaEBO(int product, int j, QList<double> *targetAggregateFillRates) {
	// pre	: 1 <= product <= |I| /\ 0 <= location <= |J|
	// ret	: deltaEBOij

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(0 <= j);
	assert(j <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	for (int l = 1; l <= network->sizeRetailers(); l++){

		double Mj = 0.0;
		double EBOlS = 0.0;

		for (int i = 1; i <= network->sizeProducts(); i++){
			Mj = Mj + network->getArrivalRateAtRetailer(i, l);
			EBOlS = EBOlS + ePartsOnBackorderAtRetailer(i, l);
		} // for

		double m1 = 0;

		m1 = (EBOlS / Mj) - (1.0 - (*targetAggregateFillRates)[l-1]);

		if (m1 <= 0) {
			m1 = 0;
		} // if

		double EBOlSplus1 = 0.0;

		double SiX = 0.0;

		if (j == 0) {
			SiX = network->getBaseStockLevelAtWarehouse(product);
			network->setBaseStockLevelAtWarehouse(product, SiX + 1);
		} else {
			SiX = network->getBaseStockLevelAtRetailer(product, j);
			network->setBaseStockLevelAtRetailer(product, j, SiX + 1);
		}

		for (int i = 1; i <= network->sizeProducts(); i++){
			EBOlSplus1 = EBOlSplus1 + ePartsOnBackorderAtRetailer(i, l);
		} // for

		if (j == 0) {
			network->setBaseStockLevelAtWarehouse(product, SiX);
		}
		else {
			network->setBaseStockLevelAtRetailer(product, j, SiX);
		}

		double m2 = 0;

		m2 = (EBOlSplus1 / Mj) - (1.0 - (*targetAggregateFillRates)[l - 1]);

		if (m2 <= 0) {
			m2 = 0;
		} // if

		result = result + (m1 - m2);

	} // for

	return result;

} // calculateDeltaEBO