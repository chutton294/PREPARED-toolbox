#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>
#include "time.h"

#ifndef GENERICFUNCTIONS_HPP_INCLUDED
#define GENERICFUNCTIONS_HPP_INCLUDED

class genericFunctions{ //class of generic, commonly used functions in data analysis.

public: 

//calculate the mean absolute error between two vectors.
void MAE(double *vector1, double * vector2, int & vectorlength, double & result);

//sorts values from largest (MLE) to smallest.
void bubbleSort(double *values,int length);

//produces a list "rank" of "values" from highest to lowest. so if rank[0] = 63, then value[63] is the largest in the list. e.g. the MLE. useful for then sorting associated parameters.
void bubbleSortRank(double *values, int *rank, int length); 

void normalise(double *values, double *norms, int length);

//computes a probDist for a point probabilities over a distribution, defined by par.
void probDist(double *par, double * prob, int & length, double & min, double & max, int & bins, double *PDF, double * binCent);

void cumDist(double *PDF, int & bins, double * CDF);

void CDFDiff(double *CDF, double *binCent, int & bins, double & maxDiff, double &parmin, double &parmax);

//square of pearsons product moment correlation coefficient
void coeffDeterm(double *p1, double * p2, int &length, double & R2);

void max(double *p1, int & length, double & maxVal, int & rank);

void min(double *p1, int & length, double & minVal, int & rank);

void confInt(double *bin, double *CDF, int & length, double & ciL, double & ciU, double & ci);

//initialises the random function using the computer clock
void initialiseRand();

void sampleNormDist(double std, double & result);

void calcMean(double *vect1, int &length, double & mean);

void calcStd(double *vect1, int &length, double & std);

void randInt(double lower, double upper, double & result);

void perturbKernSmooth(double *vect1, int length, double dirac);

void variogram(double *x, double *y, double *z, double *variance, double * bincent, double *tObs, int length, int bins);

}; //end of class genericFunctions

#endif	//GENERICFUNCTIONS_HPP_INCLUDED

