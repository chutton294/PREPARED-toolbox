// model calibration.cpp example application

#include "stdafx.h"

//header files of calibration classes.
#include "modelCalibration.h"
#include "mcSampling.h"
#include "informalLikelihoods.h"
#include "informalBayesAnalysis.h"
#include "genericFunctions.h"


//linear model function, with gradient (a) and intercept (b).
void linearModel(double a, double b, double x, double & y);

//function to generate "data" using the linear model. 
//Observations perturbed with random normal noise, with standard deviation (std).
void generateData(int length, double a, double b, double std);

double *observations; //vector to store the observations generated for use in model calibration
double *predictions; //vector to store the latest model predictions

int tObs; //number of observations used for calibration
int parameters; //total number of parameters.
int retainedParameterSets; //number of parameter sets retained for posterior analysis.


int _tmain(int argc, _TCHAR* argv[])
{	
	tObs = 1000;
	parameters = 3;
	retainedParameterSets = 15000;
	generateData(tObs, 0.6, 0.5, 2);

	//set prior parameter distributions for the model.
	double * minpar; minpar = new double[parameters];
	double *maxpar; maxpar = new double[parameters];
	minpar[0] = 0.5; maxpar[0] = 0.7;
	minpar[1] = 0.35; maxpar[1] = 0.65;
	minpar[2] = 1; maxpar[2] = 3;
	
	//initialise the calibration class.
	modelCalibration cal;
	cal.initialise(parameters, tObs, retainedParameterSets);
	
	//initialise the formal likelihood class.
	informalLikelihoods like;
	like.initialiseNSE(observations,tObs,2);

	//initialise the sampling class.
	mcSampling samp;
	samp.initialise(parameters, minpar, maxpar);

	informalBayesAnalysis iba;

	//loop through making model runs.
	for (int i = 0; i < 200000; i++){

		std::cout<<i<<std::endl;

		//sample from the specified prior distributions.
		samp.sample();

		//run the model.
		for (int i = 0; i < tObs; i++){
			linearModel(samp.par[0], samp.par[1], i+1, predictions[i]);
		}

		//calculate the result of the predictions.
		double result;
		like.runNSE(predictions,result);		
		
		//if sufficiently good, store in calibration master storage.
		cal.addSample(result, samp.par, predictions);
	}
	
	//intitialise the formal bayesian analysis class.
	iba.initialise(cal.par, cal.like, cal.tSamples, cal.tPar, samp.parMin, samp.parMax, cal.pred, cal.tObs,observations,0.95);
				
	//define the thresholds for running the informal bayesian analysis.
	int tThresh = 5;
	double *thresholds; thresholds = new double[tThresh];
	thresholds[0] = 0.1;
	thresholds[1] = 0.2;
	thresholds[2] = 0.4;
	thresholds[3] = 0.7;
	thresholds[4] = 1;

	//run model performance evaluation.
	iba.runAnalysis(tThresh,thresholds);

	//output summary tables of model parameter sensitivity to file.
	iba.outputTables("C:/Users/desktop/toolbox_demo/output tables.txt");

	//output parameter probability (and cumulative) distribution functions to file.
	iba.outputPDFCDF("C:/Users/desktop/toolbox_demo/output pdfcdf.txt");

	//output prediction intervals and confidence intervals to file.
	iba.outputPredInt("C:/Users/desktop/toolbox_demo/output confidence intervals.txt");
	
	return 0;
}

void linearModel(double  a, double  b, double  x, double & y){

	 y = (a*x) + b;

} //end of linear

void generateData(int length, double a, double b, double std){

	genericFunctions gen;
	gen.initialiseRand();

	//generate data for the linear model calibration
	observations = new double[length];
	predictions = new double[length];
	double *residuals; residuals = new double[length];

	for (int i = 0; i < length; i++){ 
		linearModel(a, b,i+1,observations[i]);

		gen.sampleNormDist(std,residuals[i]);
		observations[i] += residuals[i];
	}

	double actualStd = 0;

	gen.calcStd(residuals,length,actualStd);

	delete [] residuals;
}