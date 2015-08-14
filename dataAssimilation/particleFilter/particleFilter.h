#include "stdlib.h"
#include <iostream>
#include "time.h"
#include <fstream>
#include <istream>
#include <cstdlib>
#include "math.h"
#include <string.h>

#include "genericFunctions.h"

#ifndef PARTICLEFILTER_HPP_INCLUDED
#define PARTICLEFILTER_HPP_INCLUDED

class particleFilter {
	
public:

	int tStates; //number of states in the model to update
	int tParameters; //number of parameters in the model to update
	int tParticles; //ensemble size - number of particles.
	int tPredictions; //number of predictions available at each observation timestep 
	int sampFreq; //frequency at which the parameters are updated
	int sampCount; //count used to determn
	double delta; //parameter when pertubing a distribution
	double parPert; //perturbation of the parameter values to maintain diversity and coverage of the observations
	double statePert; //perturbations of the model states to help maintain diversity and coverage of the systme states.

	double **states; //array storing the states of a model.
	double **parameters; //array storing the parameters of a model.
	double **predictions; //array storing the predictions from the model, used in assimilation.
	double *weight; //weights assigned to each particle during resampling.
	double *copies; //number of copies made of each particle during resampling.
	double *cumWeight; //this is the cumulative weight of a particle summed over a number of timesteps.

	//stores of summary information on the ensemble.
	double *parMean; //mean value of the ensemble 
	double *stateMean; //mean value of the ensemble.
	double *parRanMin; //parameter range min.
	double *parRanMax; //parameter range max.
	double *parRan95L; //95 confidence interval.
	double *parRan95U; //95 confidence interval.
	double *stateRanMin; //state range min.
	double *stateRanMax; //state range max.
	double *stateRan95L; //95 state confidence interval.
	double *stateRan95U; //95 state confidence interval.
	
	genericFunctions gen; //object of generic functions class.

	void initialise(int TPARTICLES, int TSTATES, int TPARAMETERS, int TPREDICTIONS);

	void initialiseStates(double *stateMin, double *stateMax);

	void initialiseParameters(double *parMin, double *parMax);

	//methods for particle copying
	void stochUniResampling(double *weight, double *copies);

	void copyStates();

	void copyParameters();

	void perturbStates(double *stateMin, double *stateMax);

	void perturbParameters(double *parMin, double *parMax);

	void calculateWeights(double *observations, double std); //calculate the weights of the 
	
	void updateCumWeight();

	void resetCumWeight();

	void parameterUpdate(double *parMin, double *parMax);

	void calculateStatistics();	

}; //end of class ParticleFilter

#endif //PARTICLEFILTER_HPP_INCLUDED