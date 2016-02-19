
#include <QtCore/QCoreApplication>

#include "two_echelon_spare_part.h"

#include <iostream>


int main ( int argc, char *argv[] ) {
	
	QCoreApplication a ( argc, argv );

	QList<QList<double>*> *m = new QList<QList<double>*>();
	QList<QList<double>*> *t = new QList<QList<double>*>();
	QList<QList<double>*> *S = new QList<QList<double>*>();

	QList<double> *m0 = new QList<double>();
	m0->append(0.3);

	QList<double> *m1 = new QList<double>();
	m1->append(0.1);

	QList<double> *m2 = new QList<double>();
	m2->append(0.1);

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
	S1->append(1);

	QList<double> *S2 = new QList<double>();
	S2->append(1);

	S->append(S0);
	S->append(S1);
	S->append(S2);



	TwoEchelonSparePart *sChain = new TwoEchelonSparePart(m, t, S);
	sChain->evaluate();


	std::cout << "hello world" << std::endl;

	

	return a.exec();
}
