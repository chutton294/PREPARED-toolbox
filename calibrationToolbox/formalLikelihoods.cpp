#include "formalLikelihoods.h"


void formalLikelihoods::negLogGaussLF(double *predictions, double & result){

	//returns the log gauss sum of the final terms in the log gaussian likelihood function.
	result = ((-tObs*0.5)*log(2*3.1415926535897));
	
	double one = 1;

	for (int i = 0; i< tObs; i++){
			double var_used = std[i]*std[i];
			if (var_used == 0){ var_used = 0.0000000001;} //to prevent errors in the data through zero division.
			result += -0.5*log(var_used);
			result += (predictions[i]-observations[i])*(predictions[i]-observations[i])* -0.5 * (one/var_used);
	}
	
	//to prevent errors in calculation of one parameter set affecting the whole sampling analysis later on
	if ((result/result)!=1){
		result = -100000000;
	}
	
} //end of negLogGaussLF

void formalLikelihoods::normLogProb(double *normprob, double *loglike, int samples){

	//converts a vector of log likelihoods to a vector of normalised probabilities.

	//first find MLE.
	double max = 0;
	int max_prob = 0;
	for (int j = 0; j<samples; j++){
		if (j == 0){ 
			max = loglike[j];
		}
		if (loglike[j] > max){ 
			max = loglike[j];
			max_prob = j; 
		}
	}
	
	double sum = 0;

	for (int j = 0; j<samples; j++){
		sum += exp(loglike[j]-max);
	}

	double one = 1;
	normprob[max_prob] = one/sum;  //calculat probability of MLE.
	
	//use MLE to calculate all other probabilities in the vector.
	for (int j = 0; j<samples; j++){
		if(j != max_prob){
			normprob[j] = exp(loglike[j]-max)*normprob[max_prob];
		}
	}

} //end of function normLogProb

void formalLikelihoods::setStd(double *predictions, double STD, double code){

	//allows for consideration of heteroscadestic error in the error variance of the model.
	//if code == 0, then normal variance for all predictions
	//if code == 1, then heteroscedastic variance for all predictions

	if (code == 0){

		for (int i = 0; i < tObs; i++){
			std[i] = STD;
		}

	}

	if (code == 1){

		for (int i = 0; i < tObs; i++){
			std[i] = predictions[i]*STD;
		}

	}

}

void formalLikelihoods::initialise(double *OBSERVATIONS, int TOBS){

	observations = OBSERVATIONS;
	tObs = TOBS;
	std = new double[tObs];
	sg = &sampGaussLF;

} //end of function initialise.

void formalLikelihoods::sampGaussLF(double *par, double & result){

	//*par should be a vector of size 1, storing a standard deviation for gaussian distribution 
	//double pi =  3.14159265;
	//double denom = par[0]*sqrt(2*pi);
	//double gaussMax = static_cast<double>(1)/(denom);

	//double d;
	//
	//do{
	//	d= (rand()/static_cast<double>(RAND_MAX))*gaussMax;
	//} while (d == 0);

	//result = -2*par[0]*par[0]*log(d*denom);

	//if (result < 0){result *=-1;}

	//sqrt(result);

	//double sign = rand()/static_cast<double>(RAND_MAX);

	//if (sign > 0.5){result*=-1;}

	double rand1 = 0;
	double rand2 = 0;

	do{
		rand1 =  rand()/static_cast<double>(RAND_MAX);
	} while (rand1 == 0);

	do{
		rand2 =  rand()/static_cast<double>(RAND_MAX);
	} while (rand2 == 0);
		
	result = par[0]*sqrt(-2*log(rand1))*cos(2*3.14159265*(rand2));

}