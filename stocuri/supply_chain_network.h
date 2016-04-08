#ifndef SUPPLY_CHAIN_NETWORK_H
#define SUPPLY_CHAIN_NETWORK_H

#include <QtCore>

class SupplyChainNetwork {
	public:
		SupplyChainNetwork();
		~SupplyChainNetwork();

		virtual int loadFromFile(QString fileNameSettings, QString filenNameDemand);
		virtual int writeBaseStockLevelsToFile(QString fileName);

};

#endif

