#ifndef PYDOF_H
#define PYDOF_H

#include "pydofData.h"

class Pydof:PydofData{

	public:
		int getNoFilters(); 

		int getNoMaxCoeffs();

		int getNoCoeffs(int);

		const double* getCoeffs(int);

};

Pydof pydof; 

#endif
