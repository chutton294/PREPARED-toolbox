#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>
#include <iomanip>

#include "genericFunctions.h"
#include "parameterAnalysis.h"
#include "predictionAnalysis.h"

#ifndef INFORMALBAYESANALYSIS_HPP_INCLUDED
#define INFORMALBAYESANALYSIS_HPP_INCLUDED

class informalBayesAnalysis{ //class of functions for informal bayesian analysis.

public: 

int tSamples;
int tObs; 
int tPar;
int tThresh;

double **pars;
double *obs;
double *like;
double *parmax;
double *parmin;
double *thresholds; //an array of fractions of the total number of runs, or samples.
double * mlePar; 
double maxLike;
double **pred; //predictions from which to derive confidence intervals from GLUE
int *parRank;
parameterAnalysis * pA;
predictionAnalysis *prA;

double ci;

void initialise(double **PARS, double *LIKE, int TSAMPLES, int  TPAR, double *PARMIN, double * PARMAX, double **PRED, int & TOBS, double * OBS, double CI);

void runAnalysis(int TTHRESH, double *thrs);

void outputTables(std::string filename);

void outputPDFCDF(std::string filename);

void outputPredInt(std::string filename);

private:

void sortPar();

void sortPred();


}; //end of class genericFunctions

#endif	//INFORMALBAYESANALYSIS_HPP_INCLUDED

