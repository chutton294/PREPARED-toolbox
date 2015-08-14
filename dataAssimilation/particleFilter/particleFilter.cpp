#include "stdlib.h"
#include <iostream>
#include "time.h"
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>

#include "particleFilter.h"

void particleFilter::initialise(int TPARTICLES, int TSTATES, int TPARAMETERS, int TPREDICTIONS){

	tParticles = TPARTICLES;
	tStates = TSTATES;
	tParameters = TPARAMETERS;
	tPredictions = TPREDICTIONS;
	sampCount = 0;
	delta = 0.8;
	parPert = 0;
	statePert = 0;

	weight = new double [tParticles];
	cumWeight = new double[tParticles];
	copies = new double[tParticles];

	for (int i = 0; i < tParticles; i++){
		weight[i] = 0;
		cumWeight[i] = 0;
	}


	//dimension the main storages in the model, depending on whether there needs to be a state storage or not.
			
		states = new double *[tParticles];
		parameters = new double *[tParticles];
		predictions = new double *[tParticles];

		//summary information stores
		parMean = new double [tParameters];
		stateMean = new double [tStates];
		parRanMin = new double[tParameters]; 
		parRanMax = new double [tParameters];
		parRan95L = new double[tParameters];
		parRan95U = new double[tParameters];

		stateRanMin = new double[tStates]; 
		stateRanMax = new double[tStates];
		stateRan95L = new double[tStates];
		stateRan95U = new double[tStates];
		
		for (int i = 0; i < tParticles; i++){
			states[i] = new double[tStates];
			parameters[i] = new double[tParameters];
			predictions[i] = new double[tPredictions];
		}




}

void particleFilter::calculateStatistics(){

	for (int i = 0; i < tParameters; i++){
		parMean[i] = 0;
		for (int j = 0; j < tParticles; j++){
			parMean[i] += cumWeight[j]*parameters[j][i];
		}

		double *par; par = new double[tParticles];

		for (int j = 0; j < tParticles; j++){
			par[j] = parameters[j][i];
		}
		int junk;
		gen.min(par,tParticles,parRanMin[i],junk);
		gen.max(par,tParticles,parRanMax[i],junk);
		
		//calculate the 95% confidence intervals for the parameter.
		

		int bins = 1000;
		double upper = 0; double lower = 0;
		double *PDF; PDF = new double[bins];
		double *CDF; CDF = new double[bins];
		double *binCent; binCent = new double[bins];
		gen.probDist(par,cumWeight,tParticles,parRanMin[i],parRanMax[i],bins,PDF,binCent);
		gen.cumDist(PDF,bins,CDF);
		double CI = 0.95;
		gen.confInt(binCent,CDF,bins,parRan95L[i],parRan95U[i],CI);
		
		delete [] PDF;
		delete [] CDF;
		delete [] binCent;
		delete [] par;
	}

	for (int i = 0; i < tStates; i++){
		stateMean[i] = 0;
		for (int j = 0; j < tParticles; j++){
			stateMean[i] += weight[j]*states[j][i];
		}

		double *st; st = new double[tParticles];
		for (int j = 0; j < tParticles; j++){
			st[j] = states[j][i];
		}
		int junk;
		gen.min(st,tParticles,stateRanMin[i],junk);
		gen.max(st,tParticles,stateRanMax[i],junk);

			
		int bins = 1000;
		double upper = 0; double lower = 0;
		double *PDF; PDF = new double[bins];
		double *CDF; CDF = new double[bins];
		double *binCent; binCent = new double[bins];
		gen.probDist(st,weight,tParticles,stateRanMin[i],stateRanMax[i],bins,PDF,binCent);
		gen.cumDist(PDF,bins,CDF);
		double CI = 0.95;
		gen.confInt(binCent,CDF,bins,stateRan95L[i],stateRan95U[i],CI);
		
		delete [] PDF;
		delete [] CDF;
		delete [] binCent;
		delete [] st;
	}

}

void particleFilter::initialiseStates(double *parmin, double *parmax){

	for (int i = 0; i < tParticles; i++){
		for (int j = 0; j < tStates; j++){
			gen.randInt(parmin[j],parmax[j],states[i][j]);
		}
	}

}

void particleFilter::initialiseParameters(double *parmin, double *parmax){

	for (int i = 0; i < tParticles; i++){
		for (int j = 0; j < tParameters; j++){
			gen.randInt(parmin[j],parmax[j],parameters[i][j]);
		}
	}
	
}

void particleFilter::copyStates(){

	//copy particle states.
	for (int i=0; i< tParticles; i++){

		//find a particle with more than 1 copy to duplicate
		if (copies[i] > 1){ 

			do{
				//fina a particle with no copies to over-write the states of that copy.
				int itt = -1;
				int check = 0;

				do{
					itt +=1;
					if (copies[itt] == 0){
						check = 1;
					}
				} while (check != 1);

				//copy the state vector of the model
				for (int j = 0; j < tStates; j++){
					states[itt][j] = states[i][j];
				}

				copies[i] -= 1;
				copies[itt] += 1;

			} while (copies[i] > 1);

		}
	}

}

void particleFilter::copyParameters(){

	//copy particle states.
	for (int i=0; i< tParticles; i++){

		//find a particle with more than 1 copy to duplicate
		if (copies[i] > 1){ 

			//loop through, copying all excess copies to other particles.
			do {
				//fina a particle with no copies to over-write the states of that copy.
				int itt = -1;
				int check = 0;

				do{
					itt +=1;
					if (copies[itt] == 0){
						check = 1;
					}
				} while (check != 1);

				//copy the state vector of the model
				for (int j = 0; j < tParameters; j++){
					parameters[itt][j] = parameters[i][j];
				}

				copies[i] -= 1;
				copies[itt] += 1;

			} while (copies[i] > 1);

		}
	}

}

void particleFilter::perturbStates(double *stateMin, double *stateMax){

	double *st; st = new double[tParticles]; 

	for (int i = 0; i < tStates; i++){
		for(int j = 0; j < tParticles; j++){
			st[j] = states[j][i];
		}
		
		double pert = (stateMax[i] - stateMin[i])*statePert;
		
		gen.perturbKernSmooth(st,tParticles,delta,pert);
		
		for(int j = 0; j < tParticles; j++){
			states[j][i] = st[j];
			if (states[j][i] < stateMin[i]){states[j][i] = stateMin[i];}
			if (states[j][i] > stateMax[i]){states[j][i] = stateMax[i];}
		}
	}

	delete [] st;
}

void particleFilter::perturbParameters(double *parMin, double *parMax){

	double *par; par = new double[tParticles]; 

	for (int i = 0; i < tParameters; i++){
		for(int j = 0; j < tParticles; j++){
			par[j] = parameters[j][i];
		}
		double pert = (parMax[i] - parMin[i])*parPert;

		gen.perturbKernSmooth(par,tParticles,0.8,pert);

		//re-assign the values back to the parameters
		for(int j = 0; j < tParticles; j++){
			parameters[j][i] = par[j];
			if (parameters[j][i] < parMin[i]){parameters[j][i] = parMin[i];}
			if (parameters[j][i] > parMax[i]){parameters[j][i] = parMax[i];}
		}
	}

	delete [] par;
}

void particleFilter::resetCumWeight(){

	for(int i = 0; i < tParticles; i++){
		cumWeight[i] = 0;
	}

	sampCount = 0;

}

void particleFilter::parameterUpdate(double *parMin, double *parMax){

		if (sampCount == sampFreq){

				stochUniResampling(cumWeight,copies);

				copyParameters();

				perturbParameters(parMin, parMax);

				resetCumWeight();
			}
			
}

void particleFilter::calculateWeights(double *observations,double std){

	//calculate the particle weights using the negative log likelihood function...	
	for (int i = 0; i < tParticles; i++){

		//calculate the negative log likelihood for the ensemble member.
		weight[i] = ((-tPredictions*0.5)*log(2*3.1415926535897));
	
		double one = 1;

		for (int j = 0; j< tPredictions; j++){
			double var_used = std*std;
			if (var_used == 0){ var_used = 0.0000000001;} //to prevent errors in the data through zero division.
			weight[i] += -0.5*log(var_used);
			weight[i] += (predictions[i][j]-observations[j])*(predictions[i][j]-observations[j])* -0.5 * (one/var_used);
		}
	
		//to prevent errors in calculation of one parameter set affecting the whole sampling analysis later on
		if ((weight[i]/weight[i])!=1){
			weight[i] = -100000000;
		}

	} //end of particle loop

	//normalise the weights......
	//first find MLE.
	double max = 0;
	int max_prob = 0;
	for (int j = 0; j<tParticles; j++){
		if (j == 0){ 
			max = weight[j];
		}
		if (weight[j] > max){ 
			max = weight[j];
			max_prob = j; 
		}
	}
	
	double sum = 0;

	for (int j = 0; j<tParticles; j++){
		sum += exp(weight[j]-max);
	}

	double one = 1;
	weight[max_prob] = one/sum;  //calculat probability of MLE.
	
	//use MLE to calculate all other probabilities in the vector.
	for (int j = 0; j<tParticles; j++){
		if(j != max_prob){
			weight[j] = exp(weight[j]-max)*weight[max_prob];
		}
	}

	/*for (int i = 0; i <tParticles;i++){
		std::cout<<weight[i]<<' ';
		for (int j = 0; j < 3; j++){
			std::cout<<states[i][j]<<' ';
		}
		std::cout<<std::endl;
	}
*/
}

void particleFilter::stochUniResampling(double *weight, double *copies){


	//put all particles on a line from 0 to 1.
	double *CDF; CDF = new double[tParticles];

	CDF[0] = weight[0];
	for (int i = 1; i < tParticles; i++){
		CDF[i]= weight[i] + CDF[i-1];
		//std::cout<<CDF[i]<<std::endl;
	}


	//calculate the random length of the line.
	double rando;
	do {
		rando = rand();
	} while (rando <= 0);
	double ranmax = RAND_MAX;
	double inter = static_cast<double>(1)/tParticles;
	double line = rando*inter/ranmax;
	double sum = 0;

	//determine number of copies.
	for (int i = 0; i < tParticles; i++){
		copies[i] = 0;
		do{
			if (line < CDF[i]){
				copies[i] += 1;
				line += inter;
			}
		} while (line < CDF[i]);
		sum += copies[i];
	}

	
	////now sort the particles and copies, and use bubble sort to assign copies monotonically from smallest to largest.
	//////////////////calculate CDF and then use bubble sort to move the particle weights from lowest to highest...

	//double *temp; temp = new double[tParticles];
	//int *index; index = new int[tParticles];
	//double *cops; cops = new double[tParticles];

	//for (int i = 0; i < tParticles; i++){
	//	temp[i] = weight[i];
	//	cops[i] = copies[i];
	//	index[i] = i;
	//}


	////sort particle weights, and particle copies from smallest to largest.
	//int count = 0;
	//do {
	//	count = 0;
	//	for (int i = 0; i < tParticles -1; i++){
	//		if (temp[i] > temp[i+1]){
	//			float one = temp[i];
	//			int on = index[i];
	//			float two = temp[i+1];
	//			int tw = index[i+1];
	//			float cop1 = cops[i];
	//			float cop2 = cops[i+1];

	//			temp[i] = two;
	//			index[i] = tw;
	//			temp[i+1] = one;
	//			index[i+1] = on;
	//			cops[i] = cop2;
	//			cops[i+1] = cop1;
	//			count += 1;
	//		}

	//	}

	//} while ( count > 0);

		
	//if(susSort == 1){

	////now bubble sort the copies, but maintain the weights in their sorted position.
	//count = 0;
	//do {
	//	count = 0;
	//	for (int i = 0; i < tParticles -1; i++){
	//		if (cops[i] > cops[i+1]){
	//			float cop1 = cops[i];
	//			float cop2 = cops[i+1];
	//			cops[i] = cop2;
	//			cops[i+1] = cop1;
	//			count += 1;
	//		}

	//	}

	//} while ( count > 0);


	//}



	//int store = 0;
	////now reassign the new copies based on the sorted weights.
	//for (int i = 0; i < tParticles; i++){
	//	int ind = index[i];
	//	copies[ind] = cops[i];
	//	store += copies[ind];
	//}


	

		/*for (int i = 0; i < tParticles; i++){
			std::cout<<copies[i]<<' '<<weight[i]<<' '<<states[i][0]<<' '<<parameters[i][0]<<std::endl;
		}*/

}

void particleFilter::updateCumWeight(){

	for (int i = 0; i < tParticles; i++){
		cumWeight[i] = (cumWeight[i]*sampCount)+ weight[i];
		cumWeight[i] /= (sampCount+1);
	}

	sampCount += 1;

}
