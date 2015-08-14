#include "genericFunctions.h"

#ifndef MODELCALIBRATION_HPP_INCLUDED
#define MODELCALIBRATION_HPP_INCLUDED

class modelCalibration: public genericFunctions{

public:

int tPar; //total number of model parameters to calibration. note: not error model parameters
int tObs; //total number of observations used in calibration.
int tSamples; //the number of samples to retain during the sampling for posterior analysis
double **par; //parameter store for top tSamples
double **pred; //predictions store for the top samples, at the observation locations used in calibration.
double *like; //store fthe likelihoods associated with each sample

void initialise(int TMODELPAR, int TOBS, int TSAMPLES);
//note: only pointers to arrays are copied in function initialise to fascilitate code re-use. hence, any array passed may be modified

void addSample(double likelihood, double * parameters, double *predictions);

private:

	int filled;
	double minimum;
	int rank;


}; //end of class model calibration

#endif //MODELCALIBRATION_HPP_INCLUDED

