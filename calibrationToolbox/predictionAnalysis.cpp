#include "predictionAnalysis.h"

void predictionAnalysis::initialise(double **PRED, double * PROB, int & TSAMPLES, int & TOBS){

	tSamples = TSAMPLES;
	tObs = TOBS;
	pred = PRED;
	prob = PROB;

} //end of function initialise.

void predictionAnalysis::runAnalysis(double & CI, int BINS){

	
		bins = BINS;
		PDF = new double *[tObs];
		CDF = new double *[tObs];
		binCent = new double *[tObs];
		ciU = new double [tObs];
		ciL = new double [tObs];
		predMin = new double [tObs];
		predMax = new double [tObs];


		for (int i = 0; i < tObs; i++){
			PDF[i] = new double[bins];
			CDF[i] = new double[bins];
			binCent[i] = new double[bins];

			int junk = 0;
			//define prediction min and max values for each parameter.
			min(pred[i],tSamples,predMin[i],junk);
			max(pred[i],tSamples,predMax[i],junk);

			//calculate parameter PDF.
			probDist(pred[i],prob,tSamples,predMin[i],predMax[i],bins,PDF[i],binCent[i]);

			//calculate parameter CDF.
			cumDist(PDF[i],bins,CDF[i]);

			//calculate confidence intervals in CDF from the uniform distributions.
			confInt(binCent[i],CDF[i],bins,ciL[i],ciU[i],CI);
		
		}


}//endof runAnalysis

void predictionAnalysis::predictionIntervals(double & CI, int BINS, void (*f)(double *errPar,double & result), int tErrPar, int tErrSamp, double **par, int tPar){

		bins = BINS;
		PDF = new double *[tObs];
		CDF = new double *[tObs];
		binCent = new double *[tObs];
		ciU = new double [tObs];
		ciL = new double [tObs];
		predMin = new double [tObs];
		predMax = new double [tObs];

		int tLen = tSamples*tErrSamp;
		int count = 0;
		int count1 = 0;
		double sum = 0;

		//loop through all of the observations calculating the desired prediction intervals
		for (int i = 0; i < tObs; i++){
			std::cout<<i<<' ';
			PDF[i] = new double[bins];
			CDF[i] = new double[bins];
			binCent[i] = new double[bins];
			int junk = 0;
			sum = 0;
			double sum1 = 0;
			count1 = 0;

			double *predStore; predStore = new double[tLen];
			double *probability; probability = new double[tLen];

			//sample errors from the error model for each prediction.
			for (int j = 0; j < tSamples; j++){

				count = 0;
				//obtain parameters from total parameter set to derive error model sample.
				double *p1; p1 = new double[tErrPar];
				for (int k = tPar - tErrPar; k < tPar; k++){
					p1[count] = par[k][j];
					count += 1;
				}

				//sum += prob[j];


				double *error; error = new double[tErrSamp];
			
				for (int k = 0; k < tErrSamp; k++){; //as this is the variance, yet we want to input the std
					f(p1,error[k]);
					sum += error[k];
				}
							
				//std::cout<<pred[i][j]<<std::endl;
				for (int k = 0; k < tErrSamp; k++){
					predStore[count1] = error[k] + pred[i][j];
					probability[count1] = prob[j]/static_cast<double>(tErrSamp);
					sum1 += probability[count1];
					count1 += 1;
				}

				delete [] error;
				delete [] p1;
			}
					
			//define parameter min and max values for each parameter.			
			min(predStore,tLen,predMin[i],junk);
			max(predStore,tLen,predMax[i],junk);

			/*std::ofstream outa("C:/Users/ch294/Documents/output confidence intervals.txt",std::ios::out);
	
				for (int k = 0; k < tLen; k++){
					outa<<predStore[k]<<' '<<probability[k]<<std::endl;
				}

			outa.close();*/



			//calculate parameter PDF.
			probDist(predStore,probability,tLen,predMin[i],predMax[i],bins,PDF[i],binCent[i]);

			cumDist(PDF[i],bins,CDF[i]);

			confInt(binCent[i],CDF[i],bins,ciL[i],ciU[i],CI);

			//std::cout<<sum<<' '<<sum1<<' '<<CDF[0][bins-1]<<' '<<ciL[i]<<' '<<ciU[i]<<std::endl;

			delete []predStore;
			delete []probability;
		
		}
		

};