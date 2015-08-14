#include "math.h"


#ifndef INFORMALLIKELIHOODS_HPP_INCLUDED
#define INFORMALLIKELIHOODS_HPP_INCLUDED

class informalLikelihoods{

public:

	double obsMean;
	double denom;
	double *observations;
	int tObs;
	double exp;
	double numerator;

	void initialiseNSE(double *OBSERVATIONS, int & TOBS, double EXP);

	void initialiseNSSE(double *OBSERVATIONS, int & TOBS);

	void runNSE(double *predictions, double & result);

	void runNSSE(double *predictions, double & result);

private:
	void calculateDenom();

}; //end of class informal likelihoods

#endif //INFORMALLIKELIHOODS

