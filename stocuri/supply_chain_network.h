#ifndef SUPPLY_CHAIN_NETWORK_H
#define SUPPLY_CHAIN_NETWORK_H

#include <QtCore\qstring.h>

class SupplyChainNetwork {
	public:
		SupplyChainNetwork();
		~SupplyChainNetwork();

		virtual int loadFromFile(QString fileName);
		virtual int writeBaseStockLevelsToFile(QString fileName);

};

#endif

