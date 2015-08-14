#include "informalBayesAnalysis.h"

void informalBayesAnalysis::initialise(double **PARS, double *LIKE, int TSAMPLES, int  TPAR, double *PARMIN, double * PARMAX, double **PRED, int & TOBS, double *OBS, double CI){

pars = PARS;
like = LIKE;
tSamples = TSAMPLES;
obs = OBS;
tPar = TPAR;
parRank = new int[tSamples];
parmax = PARMAX;
parmin = PARMIN;
pred = PRED;
tObs = TOBS;
ci =  CI;

}

void informalBayesAnalysis::runAnalysis(int TTHRESH, double *thrs){

	genericFunctions gen;
	tThresh = TTHRESH;
	
	thresholds = new double[tThresh];
	pA = new parameterAnalysis[tThresh];
	prA = new predictionAnalysis[tThresh];


	for (int i = 0; i < tThresh; i++){
		thresholds[i] = thrs[i];
	}

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
	
	//now loop through each thresold, selecting the appropriate information threshold, then running pA analysis.
	for (int i = 0; i < tThresh; i++){

		//create input information.
		int sets = thresholds[i]*tSamples;
		double *prob; prob = new double[sets];
		double **p; p = new double * [tPar]; 
		double **pr; pr = new double *[tObs];

		for (int i = 0; i < tPar; i++){ 
			p[i] = new double[sets];
		}
		for (int i = 0; i < tObs; i++){
			pr[i] = new double[sets];
		}

		//do we not need to create a new store here for the prediction intervals too????????
		for (int i = 0; i < sets;i++){
			prob[i] = like[i];
			for (int j = 0; j < tPar; j++){
				p[j][i] = pars[j][i];
			}
			for (int j = 0; j < tObs; j++){
				pr[j][i] = pred[j][i];
			}
		}

		//normalise the likelihoods to create pseudo probabilities.
		gen.normalise(prob,prob,sets);


		pA[i].runAnalysis(p,prob,sets,tPar,parmin,parmax,1000,ci);



		//now create the prediction intervals for each observation.
		prA[i].initialise(pr,prob,sets,tObs);
		prA[i].runAnalysis(ci,1000);


		//no delete arrays created for that threshold.
		for (int i = 0; i < tPar; i++){ 
			delete [] p[i];
		}
		for (int i = 0; i < tObs; i++){
			delete [] pr[i] ;
		}
		delete [] p;
		delete [] pr;

	}

}

void informalBayesAnalysis::sortPar(){

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

void informalBayesAnalysis::outputTables(std::string filename){


	int prec = 3;


	//output all summary information here for all thresholds in two summary tables.
	std::ofstream outa(filename,std::ios::out);

		outa<<"Summary_information_on_informal_bayesian_calibration"<<std::endl<<std::endl;

	outa<<"	";
	outa<<" "<<"Parameter"<<std::endl;

	outa<<"Threshold"<<' '<<"Score"<<' ';

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
	}
	outa<<std::endl;
	
	outa<<' '<<"MLE"<<' ';


	for (int i = 0; i < tPar; i++){
		outa<<pA[0].maxProb[i]<<' ';
	}
	outa<<std::endl;

	for (int i = 0; i < tThresh; i++){
		outa<<thresholds[i]*tSamples<<" "<<"Mean"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[i].meanProb[j]<<' ';
		}
		outa<<std::endl;
		outa<<thresholds[i]*tSamples<<' ';
		outa<<"SD"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[i].stdProb[j]<<' ';
		}
		outa<<std::endl;
		outa<<thresholds[i]*tSamples<<' ';
		outa<<"CDFdiff"<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<pA[i].cdfDiff[j]<<' ';
		}
		outa<<std::endl;
	}

	outa<<std::endl;
	outa<<std::endl;

	outa<<"alternative tables mean and standard deviation in brackets"<<std::endl<<std::endl;

	for (int i = 0; i < tThresh; i++){
		outa<<' '<<thresholds[i]*tSamples<<' ';
	}

	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
		for (int j = 0; j < tThresh; j++){
			outa<<std::setprecision(prec)<<(pA[j].meanProb[i])<<"("<<std::setprecision(prec)<<pA[j].stdProb[i]<<")"<<' ';
		}
		outa<<std::endl;
	}
	

	outa<<std::endl;
	outa<<std::endl;

	outa<<"tables of model sensitivity of SD first, and second, cdf diff "<<std::endl<<std::endl;

	for (int i = 0; i < tThresh; i++){
		outa<<' '<<thresholds[i]*tSamples<<' ';
	}

	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
		for (int j = 0; j < tThresh; j++){
			outa<<pA[j].stdProb[i]<<' ';
		}
		outa<<std::endl;
	}

	outa<<std::endl;
	outa<<std::endl;

	
	for (int i = 0; i < tThresh; i++){
		outa<<' '<<thresholds[i]*tSamples<<' ';
	}

	outa<<std::endl;

	for (int i = 0; i < tPar; i++){
		outa<<"P"<<(i+1)<<' ';
		for (int j = 0; j < tThresh; j++){
			outa<<pA[j].cdfDiff[i]<<' ';
		}
		outa<<std::endl;
	}

	outa<<std::endl;
	outa<<std::endl;

	outa<<"tables of model mean parameter values, and specified confidence intervals "<<std::endl<<std::endl;


	outa<<"confidence interval: "<<ci*100<<"%"<<std::endl;
	outa<<std::endl;

	outa<<"	"<<"_";
	
	for (int i = 0; i < tThresh; i++){
		outa<<' '<<thresholds[i]*tSamples<<' ';
	}

	outa<<std::endl;
	
		for (int i = 0; i < tPar; i++){
			outa<<"P"<<(i+1)<<' '<<"UpperCI"<<' ';
			for (int j = 0; j < tThresh; j++){
				outa<<pA[j].ciU[i]<<' ';
			}
			outa<<std::endl;
			outa<<"P"<<(i+1)<<' '<<"Mean"<<' ';
			for (int j = 0; j < tThresh; j++){
				outa<<pA[j].meanProb[i]<<' ';
			}
			outa<<std::endl;
			outa<<"P"<<(i+1)<<' '<<"LowerCI"<<' ';
			for (int j = 0; j < tThresh; j++){
				outa<<pA[j].ciL[i]<<' ';
			}
			outa<<std::endl;
		}
	

	outa<<std::endl;
	outa<<std::endl;



	outa<<"parameter interaction sensitivity as a function of thresholds"<<std::endl;

	outa<<"threshold: ";
	for (int i = 0; i < tThresh; i++){
		outa<<thresholds[i]*tSamples<<' ';
	}
	outa<<std::endl;


	int count = 0;
	for (int ii = 0; ii < tPar-1; ii++){
		for (int j = ii+1; j < tPar; j++){
			outa<<"P"<<(ii+1)<<'-'<<"P"<<(j+1)<<' ';
			for (int i = 0; i < tThresh; i++){
			outa<<pA[i].corrCoeff[count]<<' ';
			}
			count +=1;
			outa<<std::endl;
		}
	}

	outa<<std::endl;
	outa<<std::endl;

	outa<<"sensitivity analysis tables for each thresold"<<std::endl<<std::endl;

	for (int i = 0; i < tThresh; i++){
		outa<<"threshold: "<<thresholds[i]*tSamples<<std::endl;
		outa<<std::endl;

		outa<<' ';
		for (int j = 0; j < tPar; j++){
			outa<<"P"<<(j+1)<<' ';
		}
		outa<<std::endl;
		
		int count = 0;
		for (int ii = 0; ii < tPar-1; ii++){
			outa<<"P"<<(ii+1)<<' ';
			for (int j = 0; j < ii+1; j++){
				outa<<"--"<<' ';
			}
			for (int j = ii+1; j < tPar; j++){
				outa<<pA[i].corrCoeff[count]<<' ';
				count +=1;
			}
			outa<<std::endl;
		}
		outa<<std::endl<<std::endl;

	}

	outa.close();
	

}

void informalBayesAnalysis::outputPDFCDF(std::string filename){

	std::ofstream outa(filename,std::ios::out);
	
	for (int i = 0; i < tPar; i++){
		outa<<"parameter: "<<(i+1)<<std::endl;
		
		for (int a = 0; a < tThresh; a++){
			for (int aa = 0; aa< 3; aa++){
				outa<<"thresh_"<<thresholds[a]*tSamples<<' ';
			}
		}
		outa<<std::endl;


		for (int a = 0; a < tThresh; a++){
			outa<<"parameter_value"<<' '<<"PDF"<<' '<<"CDF"<<' ';
		}
		outa<<std::endl;

		for (int a = 0; a < tThresh; a++){
			for (int j = 0; j < pA[0].bins; j++){
				for (int a = 0; a < tThresh; a++){
					outa<<pA[a].binCent[i][j]<<' '<<pA[a].PDF[i][j]<<' '<<pA[a].CDF[i][j]<<' ';				
				}
				outa<<std::endl;
			}
		outa<<std::endl;
		}

	}
	outa.close();

}

void informalBayesAnalysis::outputPredInt(std::string filename){

	std::ofstream outa(filename,std::ios::out);
	
	for (int j = 0; j < tThresh; j++){
			outa<<"lower_"<<thresholds[j]*tSamples<<' ';
			outa<<"upper_"<<thresholds[j]*tSamples<<' ';
		}

	outa<<"observations";
	
	outa<<std::endl;

	for (int i = 0; i < tObs; i++){
		for (int j = 0; j < tThresh; j++){
			outa<<prA[j].ciL[i]<<' '<<prA[j].ciU[i]<<' ';
		}
		outa<<obs[i]<<' ';
		outa<<std::endl;
	}
	
	outa.close();

} //end of output prediction interval.

void informalBayesAnalysis::sortPred(){

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