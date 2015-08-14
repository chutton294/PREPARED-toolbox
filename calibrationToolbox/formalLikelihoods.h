#include "math.h"
#include "stdlib.h"

#ifndef FORMALLIKELIHOODS_HPP_INCLUDED
#define FORMALLIKELIHOODS_HPP_INCLUDED

class formalLikelihoods{

public:

	double *observations;
	double *std;
	int tObs;

	void initialise(double *OBSERVATIONS, int TOBS);

	void (*sg)(double *, double &);   

	void negLogGaussLF(double *predictions, double & result);

	static void sampGaussLF(double *par,double & result);

	void normLogProb(double *normprob, double *loglike, int tSamples);

	void setStd(double *predictions, double STD, double code);

}; //end of class formal likelihoods

#endif //FORMALLIKELIHOODS


