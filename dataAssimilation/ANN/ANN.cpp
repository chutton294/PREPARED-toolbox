//standard libraries
#include <math.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

//developed libraries
#include "ANN.h"

void ANN::createNetwork(int NINPUT, int NHIDDEN, int NLAYERS, int NOUTPUT, int ACTUATORHIDDEN, int ACTUATOROUTPUT){

	nInput = NINPUT;
	nHidden = NHIDDEN;
	nLayers = NLAYERS;
	nOutput = NOUTPUT;

	//default actuators for the model.
	actuatorHidden = ACTUATORHIDDEN;
	actuatorOutput = ACTUATOROUTPUT;
	
	//define network nodes and node error gradients----------------------------------------
	inputNeurons = new double[nInput +1];
	inputNeurons[nInput] = -1;//set bias term as extra node for first hidden units.
	

	outputNeurons = new double[nOutput];
	outputErrorGradients = new double[nOutput];

	hiddenNeurons = new double*[nLayers];
	hiddenErrorGradients = new double*[nLayers];
	for (int i = 0; i < nLayers; i++){
		hiddenNeurons[i] = new double[nHidden + 1];
		hiddenNeurons[i][nHidden] = -1; //set bias term as extra node for hidden and output
		hiddenErrorGradients[i] = new double[nHidden]; //only for the actual nodes in the network, not bias node.
	} 
	//defined network nodes -------------------------------------

	//define network weights and delta functions for bp learning-------------------------------------
	wInputHidden = new double*[nInput +1];
	dInputHidden = new double*[nInput +1];
	doldInputHidden = new double*[nInput +1];
	for (int i = 0; i <= nInput; i++){
		wInputHidden[i] = new double[nHidden];
		dInputHidden[i] = new double[nHidden];
		doldInputHidden[i] = new double[nHidden];
	}

	wHiddenOutput = new double*[nHidden+1];
	dHiddenOutput = new double*[nHidden+1];
	doldHiddenOutput = new double*[nHidden+1];
	for (int i = 0; i <= nHidden; i++){
		wHiddenOutput[i] = new double[nOutput];
		dHiddenOutput[i] = new double[nOutput];
		doldHiddenOutput[i] = new double[nOutput];
	}

	wHiddenHidden = new double**[nLayers-1];
	dHiddenHidden = new double**[nLayers-1];
	doldHiddenHidden = new double**[nLayers-1];
	for (int i = 0; i < nLayers -1; i++){

		wHiddenHidden[i] = new double*[nHidden +1];
		dHiddenHidden[i] = new double*[nHidden +1];
		doldHiddenHidden[i] = new double*[nHidden +1];

		for (int j = 0; j <= nHidden; j++){
			wHiddenHidden[i][j] = new double[nHidden];
			dHiddenHidden[i][j] = new double[nHidden];
			doldHiddenHidden[i][j] = new double[nHidden];
		}
	}

	//converters for input and output data.
	mInput = new double[nInput];
	sdInput = new double[nInput];
	mulInput = new double[nInput];
	mOutput = new double[nOutput];
	sdOutput = new double[nOutput];
	mulOutput = new double[nOutput];

	
}

void ANN::feedForwardTrain(double *inputs){

	//set input neurons to input values.
		for(int i = 0; i < nInput; i++){
			inputNeurons[i] = inputs[i];
		}
		
	//loop through hidden layers assigning values------------------------------------------------
		for (int i = 0; i < nLayers; i++){
			if (i == 0){
				for (int j = 0; j < nHidden; j++){

					hiddenNeurons[i][j] = 0;

					for (int k = 0; k <= nInput; k++){ //equals nInput as we have a bias node in the input layer for the hidden nodes.
						hiddenNeurons[i][j] += inputNeurons[k]*wInputHidden[k][j];
					}

					hiddenNeurons[i][j] = activationFunction(actuatorHidden, hiddenNeurons[i][j]);
				}
			}

			if (i > 0){
				for (int j = 0; j < nHidden; j++){

					hiddenNeurons[i][j] = 0;

					for (int k = 0; k <= nHidden; k++){
						hiddenNeurons[i][j] += hiddenNeurons[i-1][k]*wHiddenHidden[i-1][k][j];
					}

					hiddenNeurons[i][j] = activationFunction(actuatorHidden, hiddenNeurons[i][j]);
				}
			}
		} //end of loop through hidden layers-------------------------------------------------------

		//loop final hidden layer to output layers---------------------------------------------------

		for(int i=0; i < nOutput; i++){
		
			outputNeurons[i] = 0;				
			
			//get weighted sum of inputs and bias neuron
			for(int j=0; j <= nHidden; j++){
				outputNeurons[i] += hiddenNeurons[nLayers-1][j] * wHiddenOutput[j][i];
			}
			
			//set to result of sigmoid
			outputNeurons[i] = activationFunction(actuatorOutput, outputNeurons[i]);	
		}




}

double ANN::activationFunction(int code, double value){

	double retval = 0;

	if(code == 1){ retval = 1/(1+exp(-value));}

	if(code ==2){retval = tanh(value);}

	if(code ==3){retval = value;}

	if (code ==4){
		if(value > 0){retval = 1;}
		if (value <= 0){retval = 0;}
	}

	
	return retval;
}

double ANN::getOutputGradient(int code, double observed, double predicted){

	double retval = 0;

	retval = observed - predicted;

	if (code == 1){ retval *= predicted * (1 - predicted);}

	if (code == 2){
		retval *= 1/(cosh(2*(predicted))+1)/2;
	}
		
	if (code == 4){
		retval = retval;
	}

	return retval;
}

double ANN::getHiddenGradient(int code, double weightsum, double neuronValue){

	double retval = weightsum;

	
	if (code == 1){ retval *= neuronValue * (1 - neuronValue);}

	if (code == 2){
		retval *= 1/(cosh(2*(neuronValue))+1)/2;
	}
	

	return retval;
}

void ANN::writeWeightFile(std::string filename){

	std::ofstream outdata(filename, std::ios::out);

	outdata<<"input "<<nInput<<std::endl;
	outdata<<"hidden "<<nHidden<<std::endl;
	outdata<<"hidlay "<<nLayers<<std::endl;
	outdata<<"output "<<nOutput<<std::endl;
	outdata<<"hidactiation "<<actuatorHidden<<std::endl;
	outdata<<"outactivation "<<actuatorOutput<<std::endl;
	outdata<<"con_factors: "<<std::endl;

	for (int i = 0; i < nInput; i++){
		outdata<<mInput[i]<<' '<<sdInput[i]<<' '<<mulInput[i]<<std::endl;
	}

	for (int i = 0; i < nOutput; i++){
		outdata<<mOutput[i]<<' '<<sdOutput[i]<<' '<<mulOutput[i]<<std::endl;
	}

	//now, write out the weights in the network...

	for (int i = 0; i <= nInput; i++){
		for (int j = 0; j < nHidden; j++){
			outdata<<wInputHidden[i][j]<<std::endl;
		}
	}
	for (int i = 0; i < nLayers -1; i++){
		for (int j = 0; j <= nHidden; j++){
			for (int k = 0; k < nHidden; k++){
				outdata<<wHiddenHidden[i][j][k]<<std::endl;
			}
		}
	}
	for (int i = 0; i <= nHidden; i++){
		for (int j = 0; j < nOutput; j++){
			outdata<<wHiddenOutput[i][j]<<std::endl;
		}
	}

	outdata.close();

}

void ANN::callWeightFile(std::string filename){

	char junk[20];

	std::ifstream indata(filename, std::ios::in);

	//call node specification.
	indata>>junk; indata>>nInput;
	indata>>junk; indata>>nHidden;
	indata>>junk; indata>>nLayers;
	indata>>junk; indata>>nOutput;
	indata>>junk; indata>>actuatorHidden;
	indata>>junk; indata>>actuatorOutput;
	
	//call network weights.
	createNetwork(nInput, nHidden, nLayers, nOutput,actuatorHidden,actuatorOutput);

	indata>>junk;


	for (int i = 0; i < nInput; i++){
		indata>>mInput[i];
		indata>>sdInput[i];
		indata>>mulInput[i];
	}
	for (int i = 0; i < nOutput; i++){
		indata>>mOutput[i];
		indata>>sdOutput[i];
		indata>>mulOutput[i];
	}

	for (int i = 0; i <= nInput; i++){
		for (int j = 0; j < nHidden; j++){
			indata>>wInputHidden[i][j];
		}
	}
	for (int i = 0; i < nLayers -1; i++){
		for (int j = 0; j <= nHidden; j++){
			for (int k = 0; k < nHidden; k++){
				indata>>wHiddenHidden[i][j][k];
			}
		}
	}
	for (int i = 0; i <= nHidden; i++){
		for (int j = 0; j < nOutput; j++){
			indata>>wHiddenOutput[i][j];
		}
	}

	indata.close();

}

void ANN::feedForwardRun(double *input, double *output){


	//set input neurons to input values.
		for(int i = 0; i < nInput; i++){
			inputNeurons[i] = (input[i]-mInput[i])/(sdInput[i]*mulInput[i]);
		}
		
	//loop through hidden layers assigning values------------------------------------------------
		for (int i = 0; i < nLayers; i++){
			if (i == 0){
				for (int j = 0; j < nHidden; j++){

					hiddenNeurons[i][j] = 0;

					for (int k = 0; k <= nInput; k++){ //equals nInput as we have a bias node in the input layer for the hidden nodes.
						hiddenNeurons[i][j] += inputNeurons[k]*wInputHidden[k][j];
					}

					hiddenNeurons[i][j] = activationFunction(actuatorHidden, hiddenNeurons[i][j]);
				}
			}

			if (i > 0){
				for (int j = 0; j < nHidden; j++){

					hiddenNeurons[i][j] = 0;

					for (int k = 0; k <= nHidden; k++){
						hiddenNeurons[i][j] += hiddenNeurons[i-1][k]*wHiddenHidden[i-1][k][j];
					}

					hiddenNeurons[i][j] = activationFunction(actuatorHidden, hiddenNeurons[i][j]);
				}
			}
		} //end of loop through hidden layers-------------------------------------------------------

		//loop final hidden layer to output layers---------------------------------------------------

		for(int i=0; i < nOutput; i++){
		
			outputNeurons[i] = 0;				
			
			//get weighted sum of inputs and bias neuron
			for(int j=0; j <= nHidden; j++){
				outputNeurons[i] += hiddenNeurons[nLayers-1][j] * wHiddenOutput[j][i];
			}
			
			//set to result of sigmoid
			outputNeurons[i] = activationFunction(actuatorOutput, outputNeurons[i]);	

			output[i] = (outputNeurons[i]*mulOutput[i]*sdOutput[i])+mOutput[i];
		}


	
}

void ANN::batchRun(double **inputs, double **outputs, int length){

	double *input; input = new double[nInput];
	double *output; output = new double[nOutput];


	for (int i = 0; i < length; i++){

		for (int j = 0; j < nInput; j++){
			input[j] = inputs[j][i];
		}
			
		feedForwardRun(input,output);

		for (int j = 0; j < nOutput; j++){
			outputs[j][i] = output[j];
		}

	}
	
}




