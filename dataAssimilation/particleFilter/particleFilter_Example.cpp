#include "stdafx.h"
#include "iostream"
#include "string.h"
#include "genericFunctions.h"
#include "particleFilter.h"

//prior ranges for the parameter value estimates.
double *parmin; double *parmax; 

//prior ranges from which the initial state conditions are initialised.
double *statemin; double *statemax; 

//function to generate the data used in assimilation from the model.
void generateData(double noise); 
 
//these are functions to be defined by the user: 
void setStatesParameters(int particle); //function to set the states and parameters from each particle member's values
void getStates(int particle); //function to get the new states from the model
void runModel(double inputs, double &outputs); //function to run the model

//initialise particle filter object
particleFilter pf;

//parameters
int tStates = 1;
int tParameters = 1;
int particles = 200;
int tPredictions = 1;

//number of total itterations.
int tObs = 10000;

//store to show output from the assimilation run for all timesteps.
void createStorage();
void setStorage(int count);
void outputResults(std::string);
double **parMax;
double **parMin; 
double **parMean; 
double **par95U;
double **par95L;
double **stateMean;
double **stateMax;
double **stateMin;
double **state95U;
double **state95L;
double **state; 
double *parameter; 
double *inputs;
double *outputs;
double *mStates; 
double *mParameters; 

int _tmain(int argc, _TCHAR* argv[])
{
	//generate data for particle filter application.
	generateData(0.04);

	//create storage for the model results.
	createStorage();
	
	//initialise the particle filter class
	pf.initialise(particles,tStates,tParameters,tPredictions);

	//initialise the parameters of the particle filter randomly from pre-defined prior distributions.
	pf.initialiseParameters(parmin, parmax);
	pf.initialiseStates(statemin, statemax);

	//set the frequency of sampling for the model parameters.	
	pf.sampFreq = 200;
	pf.parPert = 0.05;
	
	//run through the observations (e.g. run through the timesteps).
	int count = 1;
	double *observations; observations = new double[tPredictions];

	do {

		//loop through the particles, propagating them to the next timestep
		for (int i = 0; i < pf.tParticles; i++){
			
			setStatesParameters(i); //set the states and parameters in the model

			runModel(inputs[count],pf.predictions[i][0]); //run the model forwards
			
			getStates(i); //get the new states at the next timestep.
		}
		
		//assign the system observations from the store.
		observations[0] = outputs[count];
		
		//calculate the weights for the current timestep.
		pf.calculateWeights(observations,0.04);
	
		//update the cumulative particle weights.
		pf.updateCumWeight();
		
		//calculate statistics of the prediction before resampling.
		pf.calculateStatistics();

		//determine particle copies.
		pf.stochUniResampling(pf.weight,pf.copies);

		//copy states using the copies identified in stochastic universal resampling.
		pf.copyStates();	

		//perturb states with random noise
		pf.perturbStates(statemin,statemax);

		//updates the parameters once pf.sampCount is equal to pf.sampFreq
		pf.parameterUpdate(parmin,parmax);

		setStorage(count);

		count += 1;
	} while(count < tObs);

	outputResults("C:/Users/ch294/Desktop/testfile.txt");
			
	return 0;
}


void generateData(double std){
	
	genericFunctions gen;

	state = new double *[tStates]; for (int i = 0; i < tStates; i++){state[i] = new double[tObs];}
	parameter = new double[tParameters];
	parmin = new double[tParameters]; 
	parmax = new double[tParameters];
	statemin = new double[tStates];
	statemax = new double[tStates];
	inputs = new double[tObs];
	outputs = new double[tObs];

	mStates = new double[tStates];
	mParameters = new double[tStates];

	//specify parameters and uniform distributions to initialise the particles
	parameter[0] = 30; parmin[0] = 1; parmax[0] = 70;
	
	//specify the state initial conditions and the uniform distributions to initialise the state initial conditions
	state[0][0] = 3; statemin[0] = 0; statemax[0] = 200;
	
	for (int i = 1; i < tObs; i++){
		inputs[i] = sin(0.03*(i))+cos(0.02*(i))+2;
			
		double flux1 = state[0][i-1]/parameter[0];

		state[0][i] = state[0][i-1] + inputs[i] - flux1;
	
		outputs[i] = flux1;

		//add noise to the generated values 
		double n;
		gen.sampleNormDist(std,n);
		outputs[i] += n;
		if (outputs[i] < 0){outputs[i] = 0.0001;}

	}
			
}

void setStatesParameters(int particle){

	//set the states of the model
	for (int i = 0; i < pf.tParameters;i++){
		mParameters[i] = pf.parameters[particle][i];
	}

	//set the states of the model
	for (int i = 0; i < pf.tStates;i++){
		mStates[i] = pf.states[particle][i];
	}

}

void getStates(int particle){

	
	//get the states of the model
	for (int i = 0; i < pf.tStates;i++){
		 pf.states[particle][i] = mStates[i];
	}

}

void runModel(double inputs, double & outputs){

		double flux1 = mStates[0]/mParameters[0];

		mStates[0] += inputs - flux1;

		outputs = flux1;

}

void createStorage(){
	
	parMax = new double *[tParameters];
	parMin = new double *[tParameters];
	parMean = new double *[tParameters];
	par95U = new double *[tParameters];
	par95L = new double *[tParameters];

	stateMax = new double *[tStates];
	stateMin = new double *[tStates];
	state95U = new double *[tStates];
	state95L = new double *[tStates];
	stateMean = new double *[tStates];

	for (int i = 0; i < tParameters; i++){
		parMax[i] = new double[tObs];
		parMin[i] = new double[tObs];
		parMean[i] = new double[tObs];
		par95U[i] = new double[tObs];
		par95L[i] = new double[tObs];
	}

	for (int i = 0; i < tStates; i++){
		stateMax[i] = new double[tObs]; 
		stateMin[i] = new double[tObs]; 
		stateMean[i] = new double[tObs];
		state95U[i] = new double[tObs];
		state95L[i] = new double[tObs];
	}

}

void setStorage(int count){

	for (int i = 0; i < tParameters; i++){
		parMax[i][count-1] = pf.parRanMax[i];
		parMin[i][count-1] = pf.parRanMin[i];
		parMean[i][count-1] = pf.parMean[i];
		par95U[i][count-1] = pf.parRan95U[i];
		par95L[i][count-1] = pf.parRan95L[i];
	}

	for (int i = 0; i < tStates; i++){
		stateMax[i][count-1] = pf.stateRanMax[i];
		stateMin[i][count-1] = pf.stateRanMin[i];
		stateMean[i][count-1] = pf.stateMean[i];
		state95U[i][count-1] = pf.stateRan95U[i];
		state95L[i][count-1] =  pf.stateRan95L[i];
	}

}

void outputResults(std::string filename){

	std::ofstream outa(filename,std::ios::out);

	for (int j = 0; j < tParameters; j++){
		outa<<"p"<<(j+1)<<"_min"<<' '<<"p"<<(j+1)<<"_lower95"<<' '<<"p"<<(j+1)<<"_mean"<<' '<<"p"<<(j+1)<<"_upper95"<<' '<<"p"<<(j+1)<<"_max"<<' ';
	}
	
	for (int j = 0; j < tStates; j++){
		outa<<"s"<<(j+1)<<"_min"<<' '<<"s"<<(j+1)<<"_lower95"<<' '<<"s"<<(j+1)<<"_mean"<<' '<<"s"<<(j+1)<<"_upper95"<<' '<<"s"<<(j+1)<<"_max"<<' ';
	}

	outa<<"true_state"<<std::endl;
	
	for (int i = 0; i < tObs-1; i++){
		for (int j = 0; j < tParameters; j++){
			outa<<parMin[j][i]<<' '<<par95L[j][i]<<' '<<parMean[j][i]<<' '<<par95U[j][i]<<' '<<parMax[j][i]<<' ';
		}
		for (int j = 0; j < tStates; j++){
			outa<<stateMin[j][i]<<' '<<state95L[j][i]<<' '<<stateMean[j][i]<<' '<<state95U[j][i]<<' '<<stateMax[j][i]<<' ';
		}
		outa<<state[0][i]<<' ';
		outa<<std::endl;
	}

	outa.close();
		
}







