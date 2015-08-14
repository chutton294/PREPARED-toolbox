#include "modelCalibration.h"

void modelCalibration::initialise(int TPAR, int TOBS, int TSAMPLES){

	tPar = TPAR;
	tObs = TOBS;
	tSamples = TSAMPLES;
	filled = 0;

	par = new double *[tPar];
	pred = new double *[tObs];
	like = new double[tSamples];

	for (int i = 0; i < tPar; i++){
		par[i] = new double[tSamples];
	}
	for (int i = 0; i < tObs; i++){
		pred[i] = new double[tSamples];
	}

}

void modelCalibration::addSample(double likelihood, double * parameters, double *predictions){

	if(filled < tSamples){

		like[filled] = likelihood;

		for (int i = 0; i < tPar; i++){
			par[i][filled] = parameters[i];
		}
		for (int i = 0; i < tObs; i++){
			pred[i][filled] = predictions[i];
		}

		if (filled == (tSamples -1)){
			min(like,tSamples,minimum,rank);
		}
	}

	if (filled >= tSamples){
		
		if (likelihood > minimum){
			like[rank] = likelihood;

			for (int i = 0; i < tPar; i++){
				par[i][rank] = parameters[i];
			}
			for (int i = 0; i < tObs; i++){
				pred[i][rank] = predictions[i];
			}

			min(like,tSamples,minimum,rank);
		}
	}

	filled += 1;

} //end of addSample
