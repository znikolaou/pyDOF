#include <iostream>
#include "pydof.h"
using namespace std; 

int main(){
	for(int i=0;i<pydof.getNoFilters();i++){
		cout<<"Filter id="<<i<<endl;
		for(int j=0;j<pydof.getNoCoeffs(i);j++){
			cout<<"j="<<j<<" "<<pydof.getCoeffs(i)[j]<<endl;
		}
	}

	return 0;
}
