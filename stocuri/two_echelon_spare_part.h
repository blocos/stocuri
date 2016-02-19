#ifndef TWOECHELONSPAREPART_H
#define TWOECHELONSPAREPART_H

#include <iostream>
#include <cmath>

#include <QtCore>

class TwoEchelonSparePart {
	private:	

		QList<QList<double>*> *m;
		QList<QList<double>*> *t;
		QList<QList<double>*> *S;

		double pPartsInRepairAtW0(int i, int x);
		double pPartsOnHandAtW0(int i, int x);
		double pPartsOnBackorderAtW0(int i, int x);
		int factorial(int x);

	public:

		TwoEchelonSparePart();
		TwoEchelonSparePart(QList<QList<double>*> *m, QList<QList<double>*> *t, QList<QList<double>*> *S);
		~TwoEchelonSparePart();

		void evaluate();
};

#endif

