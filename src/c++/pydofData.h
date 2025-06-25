#ifndef PYDOF_DATA_H
#define PYDOF_DATA_H

//todo: implement python library to write this class!
class PydofData{
	protected: 
	       static const int noFilters=3, noMaxCoeffs=4;
	       static constexpr int noCoeffs[noFilters]={
		       2,
		       2,
		       2
	       };
	       static constexpr double coeffs[noFilters][noMaxCoeffs]={
		       {10.0,
			20.0,
			0.0,
			0.0},
		       {30.0,
			40.0,
			0.0,
			0.0},
		       {50.0,
			60.0,
			0.0,
			0.0}
	       };
	       
}; 

#endif 
