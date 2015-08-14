#include "informalLikelihoods.h"


void informalLikelihoods::initialiseNSE(double *OBSERVATIONS, int & TOBS, double EXP){

	tObs = TOBS;
	observations = OBSERVATIONS;
	exp = EXP;

	calculateDenom();
}

void informalLikelihoods::initialiseNSSE(double *OBSERVATIONS, int & TOBS){

	tObs = TOBS;
	observations = OBSERVATIONS;

	for (int i = 0; i < tObs; i++){
		obsMean += observations[i];
	}

	denom = obsMean*obsMean;
}

void informalLikelihoods::calculateDenom(){

	denom = 0;
	obsMean = 0;

	for (int i = 0; i < tObs; i++){
		obsMean += observations[i];
	}
	obsMean /= tObs;

	for (int i = 0; i < tObs; i++){
		denom += pow(abs(observations[i] - obsMean), exp);
	}
	
} //end of calculateDenom

void informalLikelihoods::runNSE(double *predictions, double & result){

	numerator = 0; 
	
	for (int i = 0; i < tObs; i++){
		numerator += pow(abs(observations[i]-predictions[i]),exp);
	}

	result = 1-(numerator/denom);

	if (result < 0){result = 0;}

} //end of runNSE

void informalLikelihoods::runNSSE(double *predictions, double & result){

	numerator = 0; 
	
	for (int i = 0; i < tObs; i++){
		numerator += pow(abs(observations[i]-predictions[i]),2);
	}

	result = 1-(numerator/denom);

	if (result < 0){result = 0;}

} //end of runNSE
