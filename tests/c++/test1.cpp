#include <iostream>
#include "pydof.h"
#include "pydofFilter.h"
using namespace std; 

int main(){
	Pydof pydf;
	pydf.printFilterInfo(); 
	pydf.printFilterCoeffs(0);

	PydofFilter pdFilter(1,0,0,"ab","cd","ef"); 

}
