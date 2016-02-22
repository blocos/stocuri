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

		unsigned long long int factorial(unsigned long long int x);
		unsigned binomialCoef(unsigned n, unsigned k);

		double pPartsInRepairAtW0(int i, int x);
		double pPartsOnHandAtW0(int i, int x);
		double pPartsOnBackorderAtW0(int i, int x);
		double pPartsOnBackorderAtW0FromWJ(int i, int j, int x);

		double ePartsOnHandAtW0(int i);
		double ePartsOnBackorderAtW0(int i);
		double ePartsOnBackorderAtW0FromWJ(int i, int j);
		
		double pPartsInRepairAtWJ(int i, int j, int x);
		double pPartsOnBackorderAtWJ(int i, int j, int x);

		double ePartsOnHandAtWJ(int i, int j);
		double ePartsOnBackorderAtWJ(int i, int j);

		void initializeGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ);
		void initializeDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ);

		void clearS(QList<QList<double>*> *aS);
		void clearGamma(QList<QList<double>*> *gamma, int cardyI, int cardyJ);
		void clearDelta(QList<QList<double>*> *delta, int cardyI, int cardyJ);

		bool stopingCriterionMet();

	public:

		TwoEchelonSparePart();
		TwoEchelonSparePart(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS);
		~TwoEchelonSparePart();

		void evaluate(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS);
		void greedyProcedure(QList<QList<double>*> *am, QList<QList<double>*> *at, QList<QList<double>*> *aS);

};

#endif

