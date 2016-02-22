#include "two_echelon_spare_part.h"

TwoEchelonSparePart::TwoEchelonSparePart() {
}

TwoEchelonSparePart::~TwoEchelonSparePart() {
}

unsigned long long int TwoEchelonSparePart::factorial(unsigned long long int x) {
  return (x == 1 || x == 0) ? 1 : factorial(x - 1) * x;
}

unsigned TwoEchelonSparePart::binomialCoef(unsigned n, unsigned k) {
	unsigned c = 1, i;
	if (k > n - k) {
		k = n - k;  /* take advantage of symmetry */
	}
	for (i = 1; i <= k; i++, n--) {
		if (c / i > UINT_MAX / n) {
			return 0;  /* return 0 on overflow */
		}
		c = c / i * n + c%i * n / i;  /* split c*n/i into (c/i*i + c%i)*n/i */
	}
	return c;
}

double TwoEchelonSparePart::pPartsInRepairAtW0(int i, int x) {
	//return pow((m->at(0)->at(i - 1) * t->at(0)->at(i - 1)), x) / factorial(x))*exp(-m->at(0)->at(i - 1)*t->at(0)->at(i - 1));

	return exp(x * log((m->at(0)->at(i - 1) * t->at(0)->at(i - 1))) - lgamma(x + 1.0) - (m->at(0)->at(i - 1) * t->at(0)->at(i - 1)));
}

double TwoEchelonSparePart::pPartsOnHandAtW0(int i, int x) {
	if (x == 0) {
		//
	}
	else {
		//return pPartsInRepairAtW0(i, (*S)[0][i - 1] - x);
	}
	return 0.0;
}

double TwoEchelonSparePart::pPartsOnBackorderAtW0(int i, int x) {
	if (x == 0) {
		double sigma = 0;
		for (int y = 0; y <= S->at(0)->at(i-1); y++){
			sigma = sigma + pPartsInRepairAtW0(i, y);
		}
		return sigma;
	}
	else {
		return pPartsInRepairAtW0(i, S->at(0)->at(i - 1) + x);
	}
}

double TwoEchelonSparePart::pPartsOnBackorderAtW0FromWJ(int i, int j, int x) {

	double result = 0;

	double prev = 0;

	int y = x;

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
		
	}

	return result;
}


double TwoEchelonSparePart::ePartsOnHandAtW0(int i) {
	double result = 0.0;

	for (int x = 0; x <= S->at(0)->at(i - 1); x++){
		result = result + ((S->at(0)->at(i - 1) - x)*pPartsInRepairAtW0(i, x));
	}

	return result;
}

double TwoEchelonSparePart::ePartsOnBackorderAtW0(int i) {
	double result = 0.0;

	result = (m->at(0)->at(i-1)*t->at(0)->at(i-1)) - S->at(0)->at(i - 1) + ePartsOnHandAtW0(i);

	return result;
}

double TwoEchelonSparePart::ePartsOnBackorderAtW0FromWJ(int i, int j) {
	double result = 0.0;

	result = (m->at(j)->at(i - 1) / m->at(0)->at(i - 1) )*ePartsOnBackorderAtW0(i);

	return result;
}

double TwoEchelonSparePart::pPartsInRepairAtWJ(int i, int j, int x){
	double result = 0.0;

	for (int k = 0; k <= x; k++) {
		// P[Y=k]*P[BOj=x-k]

		// 1. P[Y=k]
		double pyk = exp((double)k * log((m->at(j)->at(i - 1) * t->at(j)->at(i - 1))) - lgamma((double)k + 1.0) - (m->at(j)->at(i - 1) * t->at(j)->at(i - 1)));

		// 2. P[BOj = x - k]
		double pbok = pPartsOnBackorderAtW0FromWJ(i, j, x - k);

		result = result + (pyk * pbok);
	}

	return result;
}

double TwoEchelonSparePart::pPartsOnBackorderAtWJ(int i, int j, int x) {
	double result = 0.0;

	if (x == 0){
		for (int y = 0; y <= S->at(j)->at(i - 1); y++) {
			result = result + pPartsInRepairAtWJ(i, j, y);
		}
	}
	else{
		result = pPartsInRepairAtWJ(i, j, S->at(j)->at(i-1) + x);
	}

	return result;

}

double TwoEchelonSparePart::ePartsOnHandAtWJ(int i, int j) {
	double result = 0.0;

	for (int x = 0; x <= S->at(j)->at(i - 1); x++){
		result = result + ((S->at(j)->at(i-1)-x)*pPartsInRepairAtWJ(i, j, x));
	} // for

	return result;
}

double TwoEchelonSparePart::ePartsOnBackorderAtWJ(int i, int j){
	double result = 0.0;

	result = (m->at(j)->at(i - 1)*t->at(j)->at(i - 1)) + ePartsOnBackorderAtW0FromWJ(i, j) - S->at(j)->at(i-1) + ePartsOnHandAtWJ(i, j);

	return result;
}

void TwoEchelonSparePart::evaluate(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS){
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

void TwoEchelonSparePart::initializeGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		QList<double> *tmp = new QList<double>();
		for (int i = 0; i < cardyI; ++i){
			tmp->append(0);
		} // for
		gamma->append(tmp);
	} // for
} // initializeGamma

void TwoEchelonSparePart::initializeDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		QList<double> *tmp = new QList<double>();
		for (int i = 0; i < cardyI; ++i){
			tmp->append(0);
		} // for
		delta->append(tmp);
	} // for
} // initializeDelta

void TwoEchelonSparePart::clearS(QList<QList<double>*> *aS) {
	for (int j = 0; j < aS->size(); ++j) {
		for (int i = 0; i < aS->at(j)->size(); ++i){
			(*(*aS)[j])[i] = 0;
		} // for
	} // for
} // clearS

void TwoEchelonSparePart::clearGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		for (int i = 0; i < cardyI; ++i){
			(*(*gamma)[j])[i] = 0;
		} // for
	} // for
} // initializeGamma

void TwoEchelonSparePart::clearDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ) {
	for (int j = 0; j < cardyJ; ++j) {
		for (int i = 0; i < cardyI; ++i){
			(*(*delta)[j])[i] = 0;
		} // for
	} // for
} // initializeDelta

bool TwoEchelonSparePart::stopingCriterionMet(){
	bool result = false;
	return result;
}

void TwoEchelonSparePart::greedyProcedure(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS){
	
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