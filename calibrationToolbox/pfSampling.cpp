#include "pfSampling.h"

void pfSampling::initialise(int TPAR, int TOBS, int TSAMPLES){

	tPar = TPAR;
	tObs = TOBS;
	tSamples = TSAMPLES;
	gen = new genericFunctions[1];
	gen[0].initialiseRand();

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

void pfSampling::mhSampling(int a, int b){

	double ratio = like[a]-like[b];
		double rand;
		gen[0].randInt(0,1,rand);

		double test = log(rand);

		if (test > ratio){
			copySample(a,b);
		}

}

void pfSampling::resample(){

	double temp;
	gen[0].randInt(0,tSamples-1,temp);
	int initSample = temp;

	for (int i = initSample; i < tSamples-1; i++){

		mhSampling(i+1,i);
	}

	mhSampling(0,tSamples-1);

	for (int i = 0; i < initSample-1; i++){

		mhSampling(i+1,i);

	}

}

void pfSampling::copySample(int a, int b){

	//copy parameters of the sample.
	for (int i = 0; i < tPar; i++){
		par[i][a]  = par[i][b];
	}

	//copy states of the sample.
	for (int i = 0; i < tObs; i++){
		pred[i][a]  = pred[i][b];
	}
	
}

void pfSampling::perturbParameters(double *parMin, double *parMax){

	for (int i = 0; i < tPar; i++){
		gen[0].perturbKernSmooth(par[i],tSamples,0.8);	

		//check to ensure all samples are within the range of the uniform prior distributions.
		for (int j  = 0; j < tSamples; j++){
			if (par[i][j] < parMin[i]){par[i][j] ==parMin[i];}
			if (par[i][j] > parMax[i]){par[i][j] ==parMax[i];}
		}
	}	

} //end of function perturb parameters.