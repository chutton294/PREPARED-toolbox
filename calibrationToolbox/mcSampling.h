#include <time.h>
#include <stdlib.h>


#ifndef MCSAMPLING_HPP_INCLUDED
#define MCSAMPLING_HPP_INCLUDED

class mcSampling{

public:

int tPar; //total number of model parameters to calibration. note: not error model parameters
double *parMax; //array to store the max range of the uniform prior distribution
double *parMin; //array to store the minimum range of the uniform prior distribution
double *par; //vectore storing the latest mc samples from the respective uniform prior distributions.

void initialise(int TPAR, double *PARMIN, double *PARMAX);
//note: only pointers to arrays are copied in function initialise to fascilitate code re-use. hence, any array passed may be modified

void sample(); //samples randomly from the uniform prior distributions 

private: 

	double ran;

}; //end of class mcsampling

#endif //MCSAMPLING_HPP_INCLUDED

