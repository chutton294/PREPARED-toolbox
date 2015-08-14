#include "stdlib.h"
#include <cstdlib>
#include "math.h"
#include <string.h>
#include "time.h"
#include "genericFunctions.h"

#ifndef PFSAMPLING_HPP_INCLUDED
#define PFSAMPLING_HPP_INCLUDED

class pfSampling{ //class of generic, commonly used functions in data analysis.

public: 

int tObs; //total number of observations used in calibration.
int tSamples; //the number of samples to retain during the sampling for posterior analysis
int tPar; //total number of parameters in the sample.
double **par; //parameter store for top tSamples
double **pred; //predictions store for the top samples, at the observation locations used in calibration.
double *like; //store fthe likelihoods associated with each sample
double dirac; //dirac function used in paramter perturbation.
genericFunctions *gen; //initialise a pointer to a generic functions object.

void initialise(int TMODELPAR, int TOBS, int TSAMPLES);

void resample();

void perturbParameters(double *parMin, double *parMax);

private:

void copySample(int a, int b); //sample b (parameters and states) are copied to a.

void mhSampling(int a, int b); //sample b is copied to sample a.

}; //end of class pfSampling.h

#endif //PFSAMPLING_HPP_INCLUDED