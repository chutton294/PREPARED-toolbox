#include "parameterAnalysis.h"


void parameterAnalysis::runAnalysis(double **par, double *prob, int & TSAMPLES, int & TPAR, double * parmin, double * parmax, int BINS, double CI){

	//dimension arrays for storage of results.
	bins = BINS;
	tSamples = TSAMPLES;
	tPar = TPAR;
	pMin = new double[tPar];
	pMax = new double[tPar];
	cdfDiff = new double [tPar];
	maxProb = new double[tPar];
	meanProb = new double[tPar];
	stdProb = new double[tPar];
	PDF = new double *[tPar]; 
	CDF = new double *[tPar];
	binCent = new double *[tPar];
	ciU = new double [tPar];
	ciL = new double [tPar];
	
	coeffs = 0; 

	for (int i = 0; i < tPar-1; i++){
		for (int j = i+1; j < tPar; j++){
			coeffs += 1;
		}
	}
	corrCoeff = new double[coeffs];
	pc1 = new int[coeffs];
	pc2 = new int[coeffs];
	

	for (int i = 0; i < tPar; i++){
		PDF[i] = new double[bins];
		CDF[i] = new double[bins];
		binCent[i] = new double[bins];

		int junk = 0;
		//define parameter min and max values for each parameter.
		min(par[i],tSamples,pMin[i],junk);
		max(par[i],tSamples,pMax[i],junk);

		//calculate parameter PDF.
		probDist(par[i],prob,tSamples,pMin[i],pMax[i],bins,PDF[i],binCent[i]);

		//calculate parameter CDF.
		cumDist(PDF[i],bins,CDF[i]);

		//calculate difference in CDF from the uniform distributions.
		CDFDiff(CDF[i],binCent[i],bins,cdfDiff[i],parmin[i],parmax[i]);

		//calculate the specified confidence intervals on the parameter values.
		confInt(binCent[i],CDF[i],bins,ciL[i],ciU[i],CI);
	}

	//calculate correlation matrix between the parameters.
		int count = 0;
		for (int i = 0; i < tPar-1; i++){
			for (int j = i+1; j < tPar; j++){
				coeffDeterm(par[i], par[j], tSamples, corrCoeff[count]);
				pc1[count] = i;
				pc2[count] = j;
				count += 1;
			}
		}
	
	//calculate the maximum probability parameter set.
		double junk;
		int set;
		max(prob,tSamples,junk,set);

		for (int i = 0; i < tPar; i++){
			maxProb[i] = par[i][set];
		}
		for (int i = 0; i < tPar; i++){
			meanProb[i] = 0;
			stdProb[i] = 0;
			for (int j = 0; j < tSamples; j++){
				meanProb[i] += par[i][j]*prob[j];
				stdProb[i] += par[i][j]*par[i][j]*prob[j];
			}
			stdProb[i] -= (meanProb[i]*meanProb[i]);
			stdProb[i] = sqrt(stdProb[i]);
		}


	} //endof runAnalysis