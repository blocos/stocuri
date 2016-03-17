#include "greedy_algorithm.h"

// ---------------------------------------------------------------------------------------------- CONSTRUCTORS AND DESTRUCTORS --

GreedyAlgorithm::GreedyAlgorithm() {
}

GreedyAlgorithm::~GreedyAlgorithm() {
}

// ------------------------------------------------------------------------------------------------------- GETTERS AND SETTERS --
/*
unsigned long long int TwoEchelonSparePart::factorial(unsigned long long int x) {
  return (x == 1 || x == 0) ? 1 : factorial(x - 1) * x;
}

unsigned TwoEchelonSparePart::binomialCoef(unsigned n, unsigned k) {
	unsigned c = 1, i;
	if (k > n - k) {
		k = n - k;  // take advantage of symmetry
	}
	for (i = 1; i <= k; i++, n--) {
		if (c / i > UINT_MAX / n) {
			return 0;  // return 0 on overflow 
		}
		c = c / i * n + c%i * n / i;  // split c*n/i into (c/i*i + c%i)*n/i
	}
	return c;
}

*/





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
		result = tetrodotoxin.probability(lambda, x);
	} catch (std::exception& me) {
		std::cout << me.what() << std::endl;

		// catched me, try calculation with Normal approximation
		result = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
	} // tryCatchMe

	return result;

} // pPartsOnOrderAtWarehouse


double GreedyAlgorithm::pPartsOnHandAtWarehouse(int product, int x) {
	// pre	: 1 <= product <= |I| /\ x >= 0
	// ret	: -1

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(x >= 0);

	// pre-conditions satisfied

	double result = 0.0;

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
		result = result + ((Si0 - x)*pPartsOnOrderAtWarehouse(product, x));
	} // for

	return result;

} // ePartsOnHandAtWarehouse


double GreedyAlgorithm::ePartsOnBackorderAtWarehouse(int product) {
	// pre	: 1 <= product <= |I|
	// ret	: E[BOi0(Si0)]

	assert(1 <= product);
	assert(product <= network->sizeProducts());

	// pre-conditions satisfied

	double result = 0.0;

	double arrivalRatei0 = network->getArrivalRateAtWarehouse(product);
	double leadTimei0 = network->getLeadTimeToWarehouse(product);
	double Si0 = network->getBaseStockLevelAtWarehouse(product);

	result = (arrivalRatei0 * leadTimei0) - Si0 + ePartsOnHandAtWarehouse(product);

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


	/*
	double prev = 0;

	int y = x;

	// not safe!

	while (true) {
		double term1 = binomialCoef(y, x);
		double term2 = pow(m->at(j)->at(i - 1) / m->at(0)->at(i - 1), x);
		double term3 = pow(1 - (m->at(j)->at(i - 1) / m->at(0)->at(i - 1)), y - x);
		double term4 = pPartsOnBackorderAtW0(i, y);

		double value = term1 * term2 * term3 * term4;

		result = result + value;

		y = y + 1;

		if (result - prev < 0.0000001) {
			break;
		}

		prev = result;

	}*/

	return result;

} // pPartsOnBackorderAtWarehouseFromRetailer


double GreedyAlgorithm::ePartsOnBackorderAtWarehouseFromRetailer(int product, int retailer) {
	// pre	: 1 <= product <= |I| /\ 1 <= retailer <= |Jloc| /\ x >= 0
	// ret	: E[BOi0j(Si0) = x]

	assert(1 <= product);
	assert(product <= network->sizeProducts());
	assert(1 <= retailer);
	assert(retailer <= network->sizeRetailers());

	// pre-conditions satisfied

	double result = 0.0;

	double mi0 = network->getArrivalRateAtWarehouse(product);
	double mij = network->getArrivalRateAtRetailer(product, retailer);

	result = (mij / mi0)*ePartsOnBackorderAtWarehouse(product);

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
			p1 = tetrodotoxin.probability(lambda, x);
		}
		catch (std::exception& me) {
			std::cout << me.what() << std::endl;

			// catched me, try calculation with Normal approximation
			p1 = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
		} // tryCatchMe

		// 2. P[BOj = x - k]
		double p2 = pPartsOnBackorderAtWarehouseFromRetailer(product, retailer, x - k);

		result = result + (p1 * p2);
	}

	return result;

} // pPartsOnOrderAtRetailer


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

	return result;

} // pPartsOnBackorderAtRetailer


/*




double GreedyAlgorithm::ePartsOnBackorderAtWJ(int i, int j){
	double result = 0.0;

	result = (m->at(j)->at(i - 1)*t->at(j)->at(i - 1)) + ePartsOnBackorderAtW0FromWJ(i, j) - S->at(j)->at(i-1) + ePartsOnHandAtWJ(i, j);

	return result;
}

void GreedyAlgorithm::evaluate(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS){
	m = am;
	t = at;
	S = aS; 
	
	std::cout << "evaluate function invoked" << std::endl;

	for (int x = 0; x <= 3; x++) {
		std::cout << pPartsOnBackorderAtWJ(1, 1, x) << std::endl;
	}

	//std::cout << pPartsOnBackorderAtW0FromWJ(1, 1, 0) << std::endl;

	std::cout << (m->at(1)->at(0) * t->at(1)->at(0)) << std::endl;

	std::cout << "expected backorders at location 1: " << ePartsOnBackorderAtWJ(1, 1) << std::endl;
	
}

void GreedyAlgorithm::initializeGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		QList<double> *tmp = new QList<double>();
		for (int i = 0; i < cardyI; ++i){
			tmp->append(0);
		} // for
		gamma->append(tmp);
	} // for
} // initializeGamma

void GreedyAlgorithm::initializeDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		QList<double> *tmp = new QList<double>();
		for (int i = 0; i < cardyI; ++i){
			tmp->append(0);
		} // for
		delta->append(tmp);
	} // for
} // initializeDelta

void GreedyAlgorithm::clearS(QList<QList<double>*> *aS) {
	for (int j = 0; j < aS->size(); ++j) {
		for (int i = 0; i < aS->at(j)->size(); ++i){
			(*(*aS)[j])[i] = 0;
		} // for
	} // for
} // clearS

void GreedyAlgorithm::clearGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		for (int i = 0; i < cardyI; ++i){
			(*(*gamma)[j])[i] = 0;
		} // for
	} // for
} // initializeGamma

void GreedyAlgorithm::clearDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		for (int i = 0; i < cardyI; ++i){
			(*(*delta)[j])[i] = 0;
		} // for
	} // for
} // initializeDelta

bool GreedyAlgorithm::stopingCriterionMet(){
	bool result = false;
	return result;
}

void GreedyAlgorithm::greedyProcedure(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS){
	
	m = am;
	t = at;
	S = aS;

	int cardyJ = aS->size();
	int cardyJlocal = aS->size() - 1;
	int cardyI = aS->at(0)->size();

	std::cout << cardyJ << std::endl;
	std::cout << cardyJlocal << std::endl;

	// Greedy Algorithm according to Basten and Van Houtum (2014)

	// step 1

	clearS(S);

	int costS = 0;

	QList<double> EBOj = QList<double>();

	for (int j = 0; j < cardyJlocal; ++j){
		EBOj.append(0);
		std::cout << "location j=" << j + 1 << " " << EBOj[j] << std::endl;

		for (int i = 0; i < cardyI; ++i){
			std::cout << "item " << i + 1 << std::endl;
			EBOj[j] = EBOj[j] + m->at(j+1)->at(i) * (t->at(0)->at(i) + t->at(j+1)->at(i));
		} // for

		std::cout << "expected backorders at location j=" << j + 1 << " equals " << EBOj[j] << std::endl;

	} // for
	

	// step 2

	QList<QList<double>*> *gamma = new QList<QList<double>*>();
	initializeGamma(gamma, cardyI, cardyJ);

	QList<QList<double>*> *delta = new QList<QList<double>*>();
	initializeDelta(delta, cardyI, cardyJ);
	
	while (!stopingCriterionMet()){
		
		// calculate EBO for each location J-local

		for (int j = 0; j < cardyJlocal; ++j){
			EBOj[j] = 0;
			
			for (int i = 0; i < cardyI; ++i){
				EBOj[j] = EBOj[j] + ePartsOnBackorderAtWJ(i+1, j+1);
			}
		}


		// for all locations
		for (int j = 0; j < cardyJ; ++j) {
			// for all items
			for (int i = 0; i < cardyI; ++i){
				
				for (int l = 0; l < cardyJlocal; ++l) {
					
					double term1 = ((EBOj[l] / m->at(l + 1)->at(i)) - 0.05 <= 0) ? 0.05 : (EBOj[l] / m->at(l + 1)->at(i)) - 0.05;
					
					// increase S(i,j)
					(*(*S)[j])[i] = (*(*S)[j])[i] + 1;
					
					double EBOSpe = 0.0;
					for (int v = 0; v < cardyI; ++v){
						EBOSpe = EBOSpe + ePartsOnBackorderAtWJ(v+1, l + 1);
					}

					double term2 = ((EBOSpe / m->at(l + 1)->at(i)) - 0.05 <= 0) ? 0.05 : (EBOSpe / m->at(l + 1)->at(i)) - 0.05;

					// decrease S(i,j)
					(*(*S)[j])[i] = (*(*S)[j])[i] - 1;


					(*(*delta)[j])[i] = (*(*delta)[j])[i] + ( term1 - term2 );

				} // for

				std::cout << "delta for (" << i << "," << j << ")" << (*(*delta)[j])[i] << std::endl;

			} // for
		} // for


		break;

		
	}

} // greedyProcedure

*/


// ------------------------------------------------------------------------------------------------------------ PUBLIC METHODS --

void GreedyAlgorithm::evaluateNetwork(TwoEchelonDistributionNetwork *network) {

}

void GreedyAlgorithm::optimizeNetwork(TwoEchelonDistributionNetwork *network) {
	this->network = network;
}