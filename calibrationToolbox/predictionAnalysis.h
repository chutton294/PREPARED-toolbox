#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>
#include "genericFunctions.h"

#ifndef PREDICTIONANALYSIS_HPP_INCLUDED
#define PREDICTIONANALYSIS_HPP_INCLUDED

class predictionAnalysis: public genericFunctions{

public:

	int tObs;
	int tSamples;
	double *ciU;
	double *ciL;
	int bins;
	double **PDF;
	double **CDF;
	double **binCent;
	double *predMin;
	double *predMax;
	double **pars;
	double **pred;
	double *prob;

	void initialise(double **PRED, double * PROB, int & TSAMPLES, int & TOBS);

	void runAnalysis(double & CI, int BINS);

	void predictionIntervals(double & CI, int BINS, void (*f)(double *errPar,double & result), int tErrPar, int tErrSamp, double **par, int tPar);
		
};

#endif	//PREDICTIONANALYSIS_HPP_INCLUDED