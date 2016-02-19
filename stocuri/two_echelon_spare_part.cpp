#include "two_echelon_spare_part.h"

TwoEchelonSparePart::TwoEchelonSparePart() {
}

TwoEchelonSparePart::TwoEchelonSparePart(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS){
	m = am;
	t = at;
	S = aS;

}

TwoEchelonSparePart::~TwoEchelonSparePart() {
}

int TwoEchelonSparePart::factorial(int x) {
  return (x == 1 || x == 0) ? 1 : factorial(x - 1) * x;
}

double TwoEchelonSparePart::pPartsInRepairAtW0(int i, int x) {
	return (pow((m->at(0)->at(i - 1) * t->at(0)->at(i - 1)), x) / factorial(x))*exp(-m->at(0)->at(i - 1)*t->at(0)->at(i - 1));
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

void TwoEchelonSparePart::evaluate() {
	std::cout << "evaluate function invoked" << std::endl;

	std::cout << pPartsOnBackorderAtW0(1, 1) << std::endl;
}