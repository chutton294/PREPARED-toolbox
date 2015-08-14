#include "formalBayesAnalysis.h"

void formalBayesAnalysis::initialise(double **PARS, double *LIKE, int TSAMPLES, int  TPAR, double *PARMIN, double * PARMAX, double **PRED, double *OBS, int & TOBS, double CI){

tPar = TPAR;
like = LIKE;
tSamples = TSAMPLES;
obs = OBS;
pars = PARS;
parRank = new int[tSamples];
parmax = PARMAX;
parmin = PARMIN;
pred = PRED;
tObs = TOBS;
ci =  CI;

}

void formalBayesAnalysis::runAnalysis(void (*f)(double *errPar,double & result),int tErrPar, int tErrSamp){

	genericFunctions gen;
		
	//calculate max likelihood parameter set
	int rank;
	gen.max(like, tSamples, maxLike, rank);
	mlePar = new double[tPar];
	for (int i = 0; i < tPar; i++){
		mlePar[i] = pars[i][rank];
	}

	//sort parameters smallest to largest
	gen.bubbleSortRank(like,parRank,tSamples); //buuble sort rank no longer sorts the likelihoods too...

	//now sort the parameters, predictions and likelihoods accordingly.
	sortPar();
	sortPred();

	gen.bubbleSort(like, tSamples);

	//loop through each specified threshold, calculating pdfs, cdfs and correlations. 
	//class for each thresholds. in each store: CDF, PDF, parsens by deviation, correlation matrix. 
	pA = new parameterAnalysis[1];
	prC = new predictionAnalysis[1];
	prP = new predictionAnalysis[1];
	prC[0].initialise(pred,like,tSamples,tObs);
	prP[0].initialise(pred,like,tSamples,tObs);

	pA[0].runAnalysis(pars,like,tSamples,tPar,parmin,parmax,1000,ci);

	prC[0].runAnalysis(ci,1000);

	//now calculate prediction intervals, in addition to the confidence intervals ran in prA[0].run analysis.
	prP[0].predictionIntervals(ci,1000,f,tErrPar,tErrSamp,pars,tPar);


} //end of run analysis

void formalBayesAnalysis::sortPar(){

	//create temp store;
	double **parTemp; parTemp = new double * [tPar]; 
	for (int i = 0; i < tPar;i++){ parTemp[i] = new double[tSamples];}

	for (int i = 0; i < tSamples; i++){
		for (int j = 0; j < tPar;j++){
			parTemp[j][i] = pars[j][parRank[i]];
		}
	}

	//copy back to original store;
	for (int i = 0; i < tSamples; i++){
		for (int j = 0; j < tPar;j++){
			pars[j][i] = parTemp[j][i];
		}
	}
	
	//delete temp store;
	for (int i = 0; i < tPar;i++){ delete [] parTemp[i];}
	delete [] parTemp;

}

void formalBayesAnalysis::outputTables(std::string filename){

	//output all summary information here for all thresholds in two summary tables.
	std::ofstream outa(filename,std::ios::out);

	outa<<"Summary_information_on_informal_bayesian_calibration"<<std::endl<<std::endl;

	outa<<"	";
	outa<<" "<<"Parameter"<<std::endl;

	outa<<" "<<"Score"<<' ';

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
	}
	outa<<std::endl;
	
	outa<<' '<<"MLE"<<' ';


	for (int i = 0; i < tPar; i++){
		outa<<pA[0].maxProb[i]<<' ';
	}
	outa<<std::endl;

		outa<<" "<<"Mean"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[0].meanProb[j]<<' ';
		}
		outa<<std::endl;
	
		outa<<" "<<"SD"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[0].stdProb[j]<<' ';
		}
		outa<<std::endl;
	
		outa<<" "<<"CDFdiff"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[0].cdfDiff[j]<<' ';
		}
		outa<<std::endl;


	outa<<std::endl;
	outa<<std::endl;

	outa<<"alternative tables mean and standard deviation in brackets"<<std::endl<<std::endl;


	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
			outa<<(pA[0].meanProb[i])<<"("<<pA[0].stdProb[i]<<")"<<' ';
		outa<<std::endl;
	}
	

	outa<<std::endl;
	outa<<std::endl;

	outa<<"tables of model sensitivity of SD first, and second, cdf diff "<<std::endl<<std::endl;


	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
		outa<<pA[0].stdProb[i]<<' ';
		
		outa<<std::endl;
	}

	outa<<std::endl;
	outa<<std::endl;


	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
		outa<<pA[0].cdfDiff[i]<<' ';
		
		outa<<std::endl;
	}

	outa<<std::endl;
	outa<<std::endl;

	outa<<"parameter interaction sensitivity"<<std::endl;

	outa<<std::endl;


	int count = 0;
	for (int ii = 0; ii < tPar-1; ii++){
		for (int j = ii+1; j < tPar; j++){
			outa<<"P"<<(ii+1)<<'-'<<"P"<<(j+1)<<' ';
			outa<<pA[0].corrCoeff[count]<<' ';
			count +=1;
			outa<<std::endl;
		}
	}

	outa<<std::endl;
	outa<<std::endl;

	outa<<"sensitivity analysis tables"<<std::endl<<std::endl;

		outa<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<"P"<<(j+1)<<' ';
		}
		outa<<std::endl;
		
		count = 0;
		for (int ii = 0; ii < tPar-1; ii++){
			outa<<"P"<<(ii+1)<<' ';
			for (int j = 0; j < ii+1; j++){
				outa<<"--"<<' ';
			}
			for (int j = ii+1; j < tPar; j++){
				outa<<pA[0].corrCoeff[count]<<' ';
				count +=1;
			}
			outa<<std::endl;
		}
		outa<<std::endl<<std::endl;



	outa.close();
	

}

void formalBayesAnalysis::outputPDFCDF(std::string filename){

	std::ofstream outa(filename,std::ios::out);
	
	


	for (int i = 0; i < tPar; i++){
		outa<<"parameter: "<<(i+1)<<std::endl;
		outa<<"parameter_value"<<' '<<"PDF"<<' '<<"CDF"<<std::endl;

			for (int j = 0; j < pA[0].bins; j++){
					outa<<pA[0].binCent[i][j]<<' '<<pA[0].PDF[i][j]<<' '<<pA[0].CDF[i][j]<<' ';	
					outa<<std::endl;
			}
		outa<<std::endl;
	}

	outa.close();

}

void formalBayesAnalysis::outputPredInt(std::string filename){

	std::ofstream outa(filename,std::ios::out);
	
	outa<<"lower_pi"<<' '<<"lower_ci"<<' '<<"observation"<<' '<<"upper_ci"<<' '<<"upper_pi"<<std::endl;
	
	for (int i = 0; i < tObs; i++){
			outa<<prP[0].ciL[i]<<' '<<prC[0].ciL[i]<<' '<<obs[i]<<' '<<prC[0].ciU[i]<<' '<<prP[0].ciU[i];
		outa<<std::endl;
	}
	
	outa.close();

} //end of output prediction interval.

void formalBayesAnalysis::sortPred(){

	//create temp store;
	double **predTemp; predTemp = new double * [tObs];
	for (int i = 0; i < tObs; i++){predTemp[i] = new double [tSamples];}
	
	for (int k = 0; k < tSamples; k++){
		for (int i = 0; i < tObs; i++){
				predTemp[i][k] = pred[i][parRank[k]];
		}
	}
		
	for (int k = 0; k < tSamples; k++){
		for (int i = 0; i < tObs; i++){
			 pred[i][k] = predTemp[i][k];
		}
	}
			
	
	//delete temp store;
	for (int i = 0; i < tObs;i++){ delete [] predTemp[i];}
	delete [] predTemp;

} //end of sort pred.