#include "mcSampling.h"

void mcSampling::initialise(int TPAR, double *PARMIN, double *PARMAX){

tPar = TPAR;
parMin = PARMIN;
parMax = PARMAX; 
par = new double [tPar];

//srand(time(NULL)); //initialises the random sampling from the current clock time.

} //end of initialise();

void mcSampling::sample(){

	for (int i = 0; i < tPar; i++){
		ran = rand();
		ran /= RAND_MAX;
		par[i] = (ran*(parMax[i]-parMin[i]))+parMin[i];
	}

} //end of function sample();