#include "genericFunctions.h"


void genericFunctions::MAE(double *v1, double * v2, int & length,double & result){

	result = 0;

	for (int i = 0; i < length; i ++){
		result += abs(v1[i]-v2[i]);
	}
	result /=  static_cast<double>(length);

}

void genericFunctions::bubbleSort(double *values, int length){
	   
int swaps;
double temp;
	
do {
	swaps = 0;
	for(int i = 0; i < length-1;i++){
		if (values[i]<values[i+1]){
			temp = values[i];
			values[i] = values[i+1];
			values[i+1] = temp;
			swaps += 1;}
	}
} while (swaps > 0);


} //end bubbleSort

void genericFunctions::bubbleSortRank(double *scores, int *rank, int length){

double *values; values = new double[length];
for (int i = 0; i < length; i++){ 
	rank[i] = i;
	values[i] = scores[i];
}
   
int swaps;
double temp;
int temp1;
		
	do {
		swaps = 0;
		for(int i = 0; i < length-1;i++){
			if (values[i]<values[i+1]){
				temp = values[i];
				values[i] = values[i+1];
				values[i+1] = temp;
				temp1 = rank[i];
				rank[i] = rank[i+1];
				rank[i+1] = temp1;
				swaps += 1;
			}
		}
	} while (swaps > 0);
	// end sort

delete [] values;

} //end bubbleSort 

void genericFunctions::normalise(double *values, double *norms, int length){

double sum = 0;

for (int i = 0; i < length; i++){
	sum += values[i];
}


for (int i = 0; i < length; i++){
	norms[i] = values[i]/sum;
}


}//end normalise

void genericFunctions::probDist(double *par, double * prob, int & length, double & min, double & max, int & bins, double *PDF, double * binCent){

	double check;
	int count;

	double sum = 0;

	max += max*0.000001; //add very small value to ensure the max value is assigned into the final bin.
	min -= min*0.000001;

	//set bin width and centre value (parameter) of each bin.
	double binWidth = (max - min)/static_cast <double>(bins);
	binCent[0] = min + 0.5*binWidth;
	PDF[0] = 0;
	for (int i = 1; i < bins; i++){
		binCent[i] = binCent[i-1] + binWidth; 
		PDF[i] = 0;
	}

	//assign probabilities to each bin.
	for (int i = 0; i < length; i++){

		double bin1 = (par[i]-min)/binWidth;
		double bin2;
		double decimal = modf(bin1,&bin2);
		PDF[static_cast<int>(bin2)] += prob[i];

		//check = 0;
		//count = 0;
		//do {
		//	double lower = min + (count*binWidth);
		//	if (par[i] >= lower){
		//		if (par[i] < (lower + binWidth)){
		//			PDF[count] += prob[i];
		//			check = 1;
		//			//if (PDF[count] < 0){
		//			//	float der = 10*10;
		//			//}
		//		}
		//	}
		//	count += 1;
		//} while (check == 0);
	}
	
} //end of histogram function.

void genericFunctions::cumDist(double *PDF, int & bins, double * CDF){

 CDF[0] = PDF[0];

 for (int i = 1; i < bins; i++){
	CDF[i] = CDF[i-1] + PDF[i];
 }
 
}//end of CDF function

void genericFunctions::CDFDiff(double *CDF, double *binCent, int & bins, double & area, double &parmin, double &parmax){

	area = 0;
	double range = parmax - parmin;


	for (int i = 0; i < bins; i++){
		double uni = (binCent[i]-parmin)/range;
		area += abs(CDF[i]-uni);
	}

	area /= static_cast<double>(bins);
} //end of CDFDiff

void genericFunctions::coeffDeterm(double *px, double * py, int & length, double & R2){

	//calculate means
	double *coeff; coeff = new double [6];

	for (int i = 0; i < 6; i++){coeff[i] = 0;}

	for (int i = 0; i < length; i++){
		coeff[0] += py[i]*px[i];
		coeff[1] += py[i];
		coeff[2] += px[i];
		coeff[3] += px[i]*px[i];
		coeff[4] += py[i]*py[i];
	}

	coeff[0]*=length;
	coeff[3]*=length;
	coeff[4]*=length;

	coeff[0] = coeff[0]-(coeff[1]*coeff[2]);
	coeff[1] = coeff[4]-(coeff[1]*coeff[1]);
	coeff[3] = coeff[3]-(coeff[2]*coeff[2]);
	coeff[1] = sqrt(coeff[1]);
	coeff[3] = sqrt(coeff[3]);

	R2 = coeff[0]/(coeff[1]*coeff[3]);
	R2 = R2*R2;
	
	delete[]coeff;

} //end of coeffDeterm

void genericFunctions::max(double *p1, int & length, double & maxVal, int & rank){
	
	maxVal = p1[0];
	rank = 0;

	for (int i = 1; i < length; i++){
		if (p1[i] > maxVal){ maxVal = p1[i]; rank = i;}
	}

} //end of max

void genericFunctions::min(double *p1, int & length, double & minVal, int & rank){

	minVal = p1[0];
	rank = 0;

	for (int i = 1; i < length; i++){
		if (p1[i] < minVal){ minVal = p1[i]; rank = i;}
	}


} //end of min

void genericFunctions::confInt(double *bin, double *CDF, int & length, double & ciL, double & ciU, double & ci){

	double ciLower = (1-ci)/2;
	double ciUpper = 1- ciLower;


	for (int i = 0; i < length -1; i ++){

		if(CDF[0] > ciLower){
			ciL = bin[0];
		}
		else{
			if (CDF[i] <= ciLower){
				if (CDF[i+1] > ciLower){
					double range = CDF[i+1]-CDF[i];
					ciL = (((CDF[i+1]-ciLower)*bin[i])/range)+(((ciLower-CDF[i])*bin[i+1])/range);
				}
			}
		}

		
		if (CDF[i] <= ciUpper){
			if (CDF[i+1] > ciUpper){
				double range = CDF[i+1]-CDF[i];
				ciU = (((CDF[i+1]-ciUpper)/range)*bin[i])+(((ciUpper-CDF[i])/range)*bin[i+1]);
			}
		}

		}

	
} //end of conf interval.

void genericFunctions::sampleNormDist(double std, double & result){

	//adds noise purturbation using the mean and variance specified in the function call.
		

	double rand1 = 0;
	double rand2 = 0;

	do{
		rand1 =  rand()/static_cast<double>(RAND_MAX);
	} while (rand1 == 0);

	do{
		rand2 =  rand()/static_cast<double>(RAND_MAX);
	} while (rand2 == 0);
		
	result = std*sqrt(-2*log(rand1))*cos(2*3.14159265*(rand2));

} //end of sampleNormDist

void genericFunctions::calcMean(double *vect1, int &length, double & mean){

	mean = 0;

	for (int i = 0; i < length; i++){
		mean += vect1[i];
	}

	float der = 10*10;
	mean /= static_cast<double>(length);

}

void genericFunctions::calcStd(double *vect1, int &length, double &std){

	double mean;
	std = 0;

	calcMean(vect1,length,mean);

	for (int i = 0; i < length; i++){
		std += (vect1[i]-mean)*(vect1[i]-mean);
	}

	std/=static_cast<double>(length);

	std = sqrt(std);
}

void genericFunctions::initialiseRand(){

	srand(time(0));

}

void genericFunctions::randInt(double lower, double upper, double & result){

	result = rand();
	result /= RAND_MAX;
	result *= (upper -lower);
	result += lower;

} //end of function randInt;

void genericFunctions::perturbKernSmooth(double *vect1, int length, double dirac, double setStd){
	
			double mean = 0;
			double stdeva = 0;

			//calculate mean for the parameter.
			calcMean(vect1,length,mean);

			//calculate standard devation for the parameter.
			calcStd(vect1,length,stdeva);			
			

			if (setStd != 0){
				if (stdeva < setStd){
				stdeva = setStd;
				}
			}
			
			//now calculate perturbations of these values.
			double val = ((1/dirac)-1)*stdeva*stdeva;
			val = sqrt(val);
			double scaler = ((3*dirac)-1)/(2*dirac);
			
			for(int i = 0; i < length; i++){
				double noise = 0;
				sampleNormDist(val,noise);
				vect1[i] = (scaler*vect1[i])+((1-scaler)*mean)+noise;	
			}

			double mean1 = 0;
			double stdeva1 = 0;
			calcMean(vect1,length,mean1);
			calcStd(vect1,length,stdeva1);	

	
} //end of perturbKernSmooth

		
