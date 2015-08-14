#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>
#include "genericFunctions.h"

#ifndef PARAMETERANALYSIS_HPP_INCLUDED
#define PARAMETERANALYSIS_HPP_INCLUDED

class parameterAnalysis: public genericFunctions{

public:

	double **CDF;
	double **PDF;
	double **binCent;
	double *maxProb; //parameter set giving the maximum probability.
	double *meanProb; //parameter set giving the mean probability. expected value.
	double *stdProb; //std of the mean probability. 
	double *cdfDiff; 
	double *pMin;
	double *pMax;
	double *corrCoeff;
	double *ciU;
	double *ciL;
	int bins;
	int tSamples;
	int tPar;
	int coeffs;
	int *pc1; 
	int *pc2; //stores the index of the parameter in the correlation.

	void runAnalysis(double **par, double *prob, int & TSAMPLES, int & TPAR, double * parmin, double * parmax, int BINS, double CI);
	
};

#endif	//PARAMETERANALYSIS_HPP_INCLUDED