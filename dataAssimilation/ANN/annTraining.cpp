//standard libraries
#include <math.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

//developed libraries
#include "annTraining.h"


void annTraining::initialiseData(int TINPUTS, int TOUTPUTS, int TLENGTH, double **INPUTS, double **OUTPUTS){

	tInputs = TINPUTS;
	tOutputs = TOUTPUTS;
	tLength = TLENGTH;

	inputs = INPUTS;
	outputs = OUTPUTS;
	
	meanInput = new double[tInputs];
	stdevaInput = new double[tInputs];
	multiInput = new double[tInputs];
	meanOutput = new double[tOutputs];
	stdevaOutput = new double[tOutputs];
	multiOutput = new double[tOutputs];


}

void annTraining::randomise(){

	double *rn; rn = new double[tLength];
	int *rank; rank = new int[tLength];

	for (int i = 0; i < tLength; i++){
		rn[i] = rand();
	}

	genericFunctions gen;

	gen.bubbleSortRank(rn,rank,tLength);

	//dimension the randomised ordering of the training data.
	rInputs = new double *[tInputs];
	rOutputs =  new double *[tOutputs];
	
	for (int i = 0; i < tInputs; i++){
		rInputs[i] = new double [tLength];
		for (int j = 0; j < tLength; j++){
			rInputs[i][j] = inputs[i][rank[j]];
		}
	}

	for (int i = 0; i < tOutputs; i++){
		rOutputs[i] = new double [tLength];
		for (int j = 0; j < tLength; j++){
			rOutputs[i][j] = outputs[i][rank[j]];
		}
	}


	delete [] rn;
	delete [] rank;

}

void annTraining::convertData(ANN &ann, double MIN, double MAX){

	//function converts and normalised all data for use in the neural network. this is achieved through Z-score normalisation, is the two function arguments both equal 0.
	//if a range if specified, the Z-score normalisation is applied, and a multiplier is used, specified differently for each data vector, to force the Z-score normalisation within the defined bounds.
	
	for (int i = 0; i < tInputs; i++){

		//calculate mean and standard deviation.
		double sum = 0;
		double dev = 0;
		for (int j = 0; j < tLength; j++){sum += inputs[i][j];}
		meanInput[i] = sum / static_cast<double>(tLength);
		for (int j = 0; j < tLength; j++){dev += (inputs[i][j]-meanInput[i])*(inputs[i][j]-meanInput[i]);}
		dev = dev/static_cast<double>(tLength); 
		stdevaInput[i] = sqrt(dev);
		 
		multiInput[i] = 1; //default value


		if (MIN !=0){
			if(MAX != 0){
				//then calculate a normalisation between the two values.
			
				//find min and max of the data
				double min = inputs[i][0];
				double max = inputs[i][0];

				for (int j = 0; j < tLength; j++){
					if (inputs[i][j] > max){ max = inputs[i][j];}
					if (inputs[i][j] < min){ min = inputs[i][j];}
				}

				//find multiplier to force data on the interval specified for use in the ANN.
				double multiplier = 0;
				double checker = 0;
				do {
					multiplier += 0.1;
					double maxi = (max - meanInput[i])/(multiplier*stdevaInput[i]);
					double mini = (min - meanInput[i])/(multiplier*stdevaInput[i]);

					if (mini > MIN ){ if (maxi < MAX){ checker = 1;}}

				} while (checker == 0);
				
				multiInput[i] = multiplier;

			}
		}

		//convert all data using the normalisation parameters.
		for (int j = 0; j < tLength; j++){
			inputs[i][j] = (inputs[i][j] - meanInput[i])/(multiInput[i]*stdevaInput[i]);
		}

	}


	for (int i = 0; i < tOutputs; i++){

		// only want to make these conversions if the output activation is not binary, as we then dont have to force any data on the interval of the activation function.
		//rather, we want to keep the original binary classification.

		//need corresponding information when calling in data for execution. 

		if (ann.actuatorOutput == 4){
			meanOutput[i] = 0;
			stdevaOutput[i] = 1;
			multiOutput[i] = 1;
		}

		if (ann.actuatorOutput != 4){

			//calculate meand and standard deviation for the outputs
			double sum = 0;
			double dev = 0;
			for (int j = 0; j < tLength; j++){sum += outputs[i][j];}
			meanOutput[i] = sum / static_cast<double>(tLength);
			for (int j = 0; j < tLength; j++){dev += (outputs[i][j]-meanOutput[i])*(outputs[i][j]-meanOutput[i]);}
			dev = dev/static_cast<double>(tLength); 
			stdevaOutput[i] = sqrt(dev);
		 
			multiOutput[i] = 1; //default value


			if (MIN !=0){
				if(MAX != 0){
					//then calculate a normalisation between the two values.
			
					//find min and max of the data
					double min = outputs[i][0];
					double max = outputs[i][0];

					for (int j = 0; j < tLength; j++){
						if (outputs[i][j] > max){ max = outputs[i][j];}
						if (outputs[i][j] < min){ min = outputs[i][j];}
					}

					//find multiplier to force data on the interval specified for use in the ANN.
					double multiplier = 0;
					double checker = 0;
					do {
						multiplier += 0.1;
						double maxi = (max - meanOutput[i])/(multiplier*stdevaOutput[i]);
						double mini = (min - meanOutput[i])/(multiplier*stdevaOutput[i]);

						if (mini > MIN ){ if (maxi < MAX){ checker = 1;}}

					} while (checker == 0);
				
					multiOutput[i] = multiplier;

				}
			}

			//convert all data using the normalisation parameters.
			for (int j = 0; j < tLength; j++){
				outputs[i][j] = (outputs[i][j] - meanOutput[i])/(multiOutput[i]*stdevaOutput[i]);
			}

		
		}
		
	}

	//before running the model the data conversions are passed to the neural network class such that they are stand alone.
	passConversions(ann);

	//produces a randomised vector of the predictions and observations.
	randomise();
	
}

void annTraining::setParameters(double LEARNINGRATE, double MOMENTUM, int ITT, int LOOPS, double WRANGE){

	learningRate = LEARNINGRATE;
	momentum = MOMENTUM;
	itt = ITT;
	wRange = WRANGE;
	loops = LOOPS;
	
}

void annTraining::trainNetwork(ANN &nn){

	//create a store for the optimal model predictions derived from the calibration.
	optPredictions = new double *[nn.nOutput];
	for (int i = 0; i < nn.nOutput; i++){
		optPredictions[i] = new double[tLength];
	}

	//create a store to store the optimal weights in the neural network model.
	int dim = ((nn.nInput+1)*nn.nHidden)+((nn.nLayers -1)*(nn.nHidden+1)*(nn.nHidden))+((nn.nHidden+1)*(nn.nOutput));
	optWeights = new double[dim];

	double min_error = 1000000000000000;

	//run through data itterations
	for (int j = 0; j < itt; j++){

		initialiseWeights(nn); //initialise weights for the training itteration.
		
		double error = 0;

		//the total number of loops through the data on which the error is calculated.
		for(int k = 0; k < loops; k++){
			error += backPropTrain(nn,1);
		}
		
		//loops back through data to calculate prediction errors
		error = 0; error = backPropTrain(nn, 0);
		
		if (j == 0){min_error = error;}

		if (error < min_error){
			min_error = error;
			//function to store optimal weights in memory for the optimal ANN network.
			storeWeights(nn);
			error = 0; error = backPropTrain(nn, 0);
			std::cout<<min_error<<std::endl;
		}


	}//end of itterations through the data.

	setWeights(nn);	

}

void annTraining::initialiseWeights(ANN &nn){

	int range = 2*wRange;
	int min = -1*wRange;

	float value = (((double)rand() / (RAND_MAX))*range) + min;


	for (int i = 0; i <= nn.nInput; i++){
		for (int j = 0; j < nn.nHidden; j++){
			nn.wInputHidden[i][j] = (((double)rand() / (RAND_MAX))*range) + min;
			nn.dInputHidden[i][j] = 0;
			nn.doldInputHidden[i][j] = 0;
		}
	}

	for (int i = 0; i < nn.nLayers -1; i++){
		for (int j = 0; j <= nn.nHidden; j++){
			for (int k = 0; k < nn.nHidden; k++){
				nn.wHiddenHidden[i][j][k] = (((double)rand() / (RAND_MAX))*range) + min;
				nn.	dHiddenHidden[i][j][k] = 0;
				nn.doldHiddenHidden[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i <= nn.nHidden; i++){
		for (int j = 0; j < nn.nOutput; j++){
			nn.wHiddenOutput[i][j] = (((double)rand() / (RAND_MAX))*range) + min;
			nn.dHiddenOutput[i][j] = 0;
			nn.doldHiddenOutput[i][j] = 0;
		}
	}

} //end of initialise weights.

double annTraining::backPropTrain(ANN &nn, int onoff){

		double *input; input = new double[nn.nInput];
		double *output; output = new double[nn.nOutput];
		double error = 0;

		if (onoff == 1){

			for (int i = 0; i < tLength; i++){ //run through BP for each MC sample.
				
						//assign input node values for the sample. need to modify to account for the differing input variables....
						for (int aa = 0; aa < nn.nInput; aa++){ 
							input[aa] = rInputs[aa][i];
						}

						//run ANN for data
						nn.feedForwardTrain(input); 
				
						//calculate real output data and calculate errors.
						for(int aa = 0; aa < nn.nOutput; aa++){
							output[aa] = rOutputs[aa][i]; 
							error += abs(nn.outputNeurons[aa]-output[aa]);
						}
				
					
						backPropagate(nn, output); 
						updateWeights(nn);
									
			}//end of run through the training length.

		}

		//so the data are generated from the optimal predictions.
		if (onoff == 0){

				for (int i = 0; i < tLength; i++){ //run through BP for each MC sample.
				
						//assign input node values for the sample. need to modify to account for the differing input variables....
						for (int aa = 0; aa < nn.nInput; aa++){ 
							input[aa] = inputs[aa][i];
						}

						//run ANN for data
						nn.feedForwardTrain(input); 
				
						//calculate real output data and calculate errors.
						for(int aa = 0; aa < nn.nOutput; aa++){
							output[aa] = outputs[aa][i]; 
							error += abs(nn.outputNeurons[aa]-output[aa]);
							optPredictions[aa][i] = nn.outputNeurons[aa];
						}			

					
				}//end of run through the training length.
								
		}

		delete[] input;
		delete[] output;


		return error/tLength;
}

void annTraining::backPropagate(ANN &nn, double *outputs){

	//check the network works first of all...
	//activation functions determine weight training, need the activation function codes in a separate class...

	//update weights between the final hidden later and the output layer-------------------------------------
	for (int i = 0; i < nn.nOutput; i++){
			nn.outputErrorGradients[i] = nn.getOutputGradient(nn.actuatorOutput, outputs[i], nn.outputNeurons[i]);

			//calculate change in weight from each hidden neuron in the final hidden neuron layer to the output in question.
			for(int j=0; j <= nn.nHidden; j++){
				double test = nn.wHiddenOutput[j][i];
				nn.dHiddenOutput[j][i] += learningRate * nn.hiddenNeurons[nn.nLayers -1][j] * nn.outputErrorGradients[i];
			}
	}
	//-------------------------------------------------------------------------------------------------------
	
	//update weights between the hidden layers in the network------------------------------------------------
	for (int i = nn.nLayers -1; i > -1; i--){
	
			for( int j = 0; j < nn.nHidden; j++){

				double weightedSum = 0;

				if (i == nn.nLayers-1){ //if the next layer is the output layer, then calculate sum.
					for( int k = 0; k < nn.nOutput; k++){
						weightedSum += nn.wHiddenOutput[j][k] * nn.outputErrorGradients[k];
					}
				}
				if (i < nn.nLayers-1){
					for( int k = 0; k < nn.nHidden; k++){ //if the next layer is another hidden layer then calculate sum.
						weightedSum += nn.wHiddenHidden[i][j][k] * nn.hiddenErrorGradients[i+1][k];
					}
				}

				nn.hiddenErrorGradients[i][j] = nn.getHiddenGradient(nn.actuatorHidden,weightedSum,nn.hiddenNeurons[i][j]);
			
				if (i > 0){ //we want the weights to the previous hidden layer updated.
					for( int k = 0; k <= nn.nHidden; k++){ 
						nn.dHiddenHidden[i-1][k][j] += learningRate * nn.hiddenNeurons[i-1][k]*nn.hiddenErrorGradients[i][j];
					}
				}

				if (i==0){
					for( int k = 0; k <= nn.nInput; k++){ 
						nn.dInputHidden[k][j] += learningRate * nn.inputNeurons[k]*nn.hiddenErrorGradients[i][j];
					}
				}

			}
	} //end of loop through hidden layers for backpropagation--------------------------------------------------------
	


}

void annTraining::updateWeights(ANN &nn){

	for (int i = 0; i <= nn.nInput; i++){
		for (int j = 0; j < nn.nHidden; j++){
			nn.wInputHidden[i][j] += nn.dInputHidden[i][j] + (momentum*nn.doldInputHidden[i][j]);
			nn.doldInputHidden[i][j] = nn.dInputHidden[i][j];
			nn.dInputHidden[i][j] = 0;
		}
	}

	for (int i = 0; i < nn.nLayers -1; i++){
		for (int j = 0; j <= nn.nHidden; j++){
			for (int k = 0; k < nn.nHidden; k++){
				nn.wHiddenHidden[i][j][k] += nn.dHiddenHidden[i][j][k] + (momentum*nn.doldHiddenHidden[i][j][k]); 
				nn.doldHiddenHidden[i][j][k] = nn.dHiddenHidden[i][j][k];
				nn.dHiddenHidden[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i <= nn.nHidden; i++){
		for (int j = 0; j < nn.nOutput; j++){
			nn.wHiddenOutput[i][j] += nn.dHiddenOutput[i][j] + (momentum*nn.doldHiddenOutput[i][j]);
			nn.doldHiddenOutput[i][j] = nn.dHiddenOutput[i][j];
			nn.dHiddenOutput[i][j] = 0;
		}
	}

}

void annTraining::passConversions(ANN &nn){

	//function that passes the conversions for the network training to the neural network itself, such that it can be used independently.

	for (int i = 0; i < nn.nInput; i++){
		nn.mInput[i] = meanInput[i];
		nn.sdInput[i] = stdevaInput[i];
		nn.mulInput[i] = multiInput[i];
	}

	for (int i = 0; i < nn.nOutput; i++){
		nn.mOutput[i] = meanOutput[i];
		nn.sdOutput[i] = stdevaOutput[i];
		nn.mulOutput[i] = multiOutput[i];
	}

	
} //end of write conversions.

void annTraining::storeWeights(ANN &nn){

	int count = 0;
	
	for (int i = 0; i <= nn.nInput; i++){
		for (int j = 0; j < nn.nHidden; j++){
			optWeights[count] = nn.wInputHidden[i][j];
			count += 1;
		}
	}
	for (int i = 0; i < nn. nLayers -1; i++){
		for (int j = 0; j <= nn.nHidden; j++){
			for (int k = 0; k < nn.nHidden; k++){
				optWeights[count] = nn.wHiddenHidden[i][j][k];
				count += 1;
			}
		}
	}
	for (int i = 0; i <= nn.nHidden; i++){
		for (int j = 0; j < nn.nOutput; j++){
			optWeights[count] = nn.wHiddenOutput[i][j];
			count += 1;
		}
	}

	
}

void annTraining::setWeights(ANN &nn){

	int count = 0;


	for (int i = 0; i <= nn.nInput; i++){
		for (int j = 0; j < nn.nHidden; j++){
			nn.wInputHidden[i][j] = optWeights[count];
			count += 1;
		}
	}
	for (int i = 0; i < nn. nLayers -1; i++){
		for (int j = 0; j <= nn.nHidden; j++){
			for (int k = 0; k < nn.nHidden; k++){
				nn.wHiddenHidden[i][j][k] = optWeights[count];
				count += 1;
			}
		}
	}
	for (int i = 0; i <= nn.nHidden; i++){
		for (int j = 0; j < nn.nOutput; j++){
			nn.wHiddenOutput[i][j] = optWeights[count];
			count += 1;
		}
	}
	

} //end of setWeights

void annTraining::writePredictions(ANN &nn, std::string filename){

	//writes out the predictions from the model for the calibration period, using the conversion factors.
	
	std::ofstream outdata(filename, std::ios::out);

	for (int i = 0; i < nn.nOutput; i++){
			outdata<<"index"<<' '<<"observation"<<' '<<"prediction"<<' ';
	}

	outdata<<std::endl;

	for (int j = 0; j < tLength; j++){
		outdata<<j<<' ';
		for (int i = 0; i < nn.nOutput; i++){
			double predictions = (optPredictions[i][j]*(multiOutput[i]*stdevaOutput[i]))+meanOutput[i];
			double observations = (outputs[i][j]*(multiOutput[i]*stdevaOutput[i]))+meanOutput[i];
			outdata<<observations<<' '<<predictions<<' ';
		}
		outdata<<std::endl;
	}
	
	

}

