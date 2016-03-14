#include <QtCore/QCoreApplication>

#include <iostream>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <exception>
#include <cassert>

#include "two_echelon_spare_part.h"

#include "poisson_distribution.h"
#include "normal_distribution.h"

#pragma fenv_access (on)

int main ( int argc, char *argv[] ) {
	

	double lambda = 4000;
	double x = 4000;
	double result = 0.0;

	PoissonDistribution tetrodotoxin;

	try {
		result = tetrodotoxin.probability(lambda, x);
		std::cout << "probability by Poisson: " << result << std::endl;
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		result = tetrodotoxin.probabilityByNormalApproximation(lambda, x);
		std::cout << "probability by Normal approximation of Poisson distribution: " << result << std::endl;
	}



	QCoreApplication a ( argc, argv );
	/*
	QList<QList<double>*> *m = new QList<QList<double>*>();
	QList<QList<double>*> *t = new QList<QList<double>*>();
	QList<QList<double>*> *S = new QList<QList<double>*>();

	QList<double> *m0 = new QList<double>();
	m0->append(0.3);

	QList<double> *m1 = new QList<double>();
	m1->append(0.1);

	QList<double> *m2 = new QList<double>();
	m2->append(0.2);

	m->append(m0);
	m->append(m1);
	m->append(m2);


	QList<double> *t0 = new QList<double>();
	t0->append(4);

	QList<double> *t1 = new QList<double>();
	t1->append(1);

	QList<double> *t2 = new QList<double>();
	t2->append(1);

	t->append(t0);
	t->append(t1);
	t->append(t2);

	QList<double> *S0 = new QList<double>();
	S0->append(2);

	QList<double> *S1 = new QList<double>();
	S1->append(19);

	QList<double> *S2 = new QList<double>();
	S2->append(1);

	S->append(S0);
	S->append(S1);
	S->append(S2);



	TwoEchelonSparePart *sChain = new TwoEchelonSparePart();
	
	
	sChain->evaluate(m, t, S);
	
	
	sChain->greedyProcedure(m, t, S);

	*/





	
	
	return a.exec();
}
