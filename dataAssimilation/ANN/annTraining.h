//standard libraries
#include <math.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>

#include "ANN.h"
#include "genericFunctions.h"

#ifndef ANNTRAINING
#define ANNTRAINING

class annTraining {

public: 

	int tInputs; //total number of inputs in the data
	int tOutputs; //total number of potential outputs in the data.
	int tLength; //total length of the training data sets


	double ** inputs; //potential number of inputs for ANN training
	double ** outputs; //potential number of outputs for ANN training

	double **rInputs; //randomised inputs for training
	double **rOutputs; //randomised outputs for training
	
	//conversion factors
	double * meanInput;
	double * stdevaInput;
	double * multiInput;
	double * meanOutput;
	double * stdevaOutput;
	double * multiOutput;
	
	double *optWeights;
	double ** optPredictions;

	//parameter used to control the learning procedure
	int itt; //total number of different starting points to run the calibration. 
	int loops; //total number of loops through the data set for each itt.
	double learningRate;
	double momentum;
	double wRange; //if wRange is set equal to 3 then ANN weights are initialised from the interval [-3,3].


	void initialiseData(int TINPUTS, int TOUTPUTS, int TLENGTH, double **INPUTS, double **OUTPUTS);

	void convertData(ANN &nn, double MIN, double MAX); //function to convert/normalise data for use with the neurons,

	void setParameters(double LEARNINGRATE, double MOMENTUM, int ITT, int LOOPS, double WRANGE);
	
	void trainNetwork(ANN &nn);

	void initialiseWeights(ANN &nn);

	double backPropTrain(ANN &nn, int onoff);

	void backPropagate(ANN &nn, double *output);

	void updateWeights(ANN &nn);

	void passConversions(ANN &nn);

	void storeWeights(ANN &nn);

	void setWeights(ANN &nn);

	void writePredictions(ANN &nn, std::string filename);

	
private:

	void randomise();
	
};

#endif