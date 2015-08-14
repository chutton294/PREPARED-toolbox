#include "stdafx.h"
#include "iostream"
#include "string.h"

#include "ANN.h"
#include "annTraining.h"


//model training data generation.
double **outputs; 
double **inputs;
void generateData();
void callInputData();

//network inputs and outputs.
int inputNodes = 1;
int outputNodes = 1;

//number of training observations.
int tObs = 100;

int _tmain(int argc, _TCHAR* argv[])
{
	//call input data.
	callInputData();
		
	//initialise an object of the neural network training class.
	annTraining train;

	//initialise an object of the neural network class.
	ANN ann;

	//create the network structure and activation functions.
	ann.createNetwork(inputNodes,5,1,outputNodes,1,3);
			
	//initialise the neural network data for training.
	train.initialiseData(inputNodes,outputNodes,tObs,inputs,outputs);

	//convert the input data as necessary for use in the ANN.
	train.convertData(ann,-2,2);

	//set the parameters required for neural network training.
	train.setParameters(0.1,0.9,4000,20,4);

	//trains the neural network passed to the function, 
	//sets the weights of that ann as the optimal weights.
	train.trainNetwork(ann);

	//output the weights and data conversion factors to file.
	ann.writeWeightFile("C:/Users/ch294/Documents/papers/demandForecasting/outputs/weightfileM6.txt"); 

	//output the predictions to file.
	train.writePredictions(ann,"C:/Users/ch294/Documents/papers/demandForecasting/outputs/predictionsM6.txt");
		
	return 0;
}


void generateData(){
	
	outputs = new double *[inputNodes]; for (int i = 0; i < inputNodes; i++){outputs[i] = new double[tObs];}
	inputs = new double *[outputNodes]; for (int i = 0; i < outputNodes; i++){inputs[i] = new double[tObs];}
	
	for (int i = 0; i < 100; i++){
		inputs[0][i] = i*0.03;
		outputs[0][i] = sin(double(2*inputs[0][i]))*exp(double(-inputs[0][i]));
	}
			
}

void callInputData(){

	std::ifstream indata("C:/Users/ch294/Documents/papers/demandForecasting/inputs/inputDataM6_rand.txt", std::ios::in);

	indata>>inputNodes;
	indata>>outputNodes;
	indata>>tObs;

	outputs = new double *[outputNodes]; for (int i = 0; i < outputNodes; i++){outputs[i] = new double[tObs];}
	inputs = new double *[inputNodes]; for (int i = 0; i < inputNodes; i++){inputs[i] = new double[tObs];}

	for (int i = 0; i < tObs; i++){
		for (int j = 0; j < inputNodes; j++){ 
			indata>>inputs[j][i];
			//std::cout<<inputs[j][i]<<' ';
		}
		for (int j = 0; j < outputNodes; j++){ 
			indata>>outputs[j][i];
			//std::cout<<outputs[j][i]<<' ';
		}
	}

	indata.close();

}