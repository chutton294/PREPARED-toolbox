//standard libraries
#include <math.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>

#ifndef NNETWORK
#define NNETWORK

class ANN{

public:

	//number of neurons in each layer, and number of hidden layers.
	int nInput, nHidden, nLayers, nOutput;

	//store for neurons.
	double* inputNeurons;
	double** hiddenNeurons;
	double* outputNeurons;

	//error gradients at nodes
	double** hiddenErrorGradients;
	double* outputErrorGradients;

	//network weights.
	double** wInputHidden;
	double** wHiddenOutput;
	double*** wHiddenHidden;

	//delta (change) in network weights.
	double** dInputHidden;
	double** dHiddenOutput;
	double*** dHiddenHidden;

	//delta from previous itteration for momentum learning.
	double** doldInputHidden;
	double** doldHiddenOutput;
	double*** doldHiddenHidden;

	int actuatorHidden;
	int actuatorOutput;

	//data conversion parameters, so that we can input real data and derive output data in the appropriate conversion.
	double * mInput; 
	double * sdInput;
	double * mulInput;
	double * mOutput;
	double * sdOutput;
	double * mulOutput;

	void createNetwork(int NINPUT, int NHIDDEN, int NLAYERS, int NOUTPUT, int ACTUATORHIDDEN, int ACTUATOROUTPUT);

	void feedForwardTrain(double *inputs);

	double activationFunction(int code, double value);

	double getOutputGradient(int code, double observed, double predicted);

	double getHiddenGradient(int code, double weightsum, double neuronValue);

	void writeWeightFile(std::string filename);

	void callWeightFile(std::string filename);

	void feedForwardRun(double *inputs, double *outputs);

	void batchRun(double **inputs, double **outputs, int length);
};

#endif