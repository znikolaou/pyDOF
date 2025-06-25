#include "pydof.h"

int Pydof::getNoFilters(){ return noFilters; }

int Pydof::getNoMaxCoeffs(){ return noMaxCoeffs; }

int Pydof::getNoCoeffs(int filterId){ return noCoeffs[filterId]; }

const double* Pydof::getCoeffs(int filterId){ return coeffs[filterId]; }
