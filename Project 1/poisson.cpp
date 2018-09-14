#include <iostream>
#include <vector>
using namespace std;

double[] subtract_vectors(double[], double[]);
double calc_norm(double[], double[], int);
double[] solve_jacobi();
double[] initialize_random_vector(int);

int main(int argc, char* argv[]){
	int x, y;

	cout<<"Enter x dimension: ";
	cin >> x;
	cout<<"Enter y dimension: ";
	cin >> y;

	int left = 0;
	int top = 1;
	int right = 1;
	int bottom = 0;

	double poissonMatrix[y][x];
	double initVariables[x] = initialize_random_vector(x, 0, 100);
	double gaussSeidel[y] = solve_gauss_seidel(1000, x, initVariables, poissonMatrix, ..., 1.0);

	return 0;
}

double

// TODO: Implement
double[] solve_jacobi(){
	return 0;
}

double[] solve_gauss_seidel(int iters, int numVars, double XO[], double a[][], double b[], double tolerance){
	int k = 1;
	double x[numVars];

	while(k <= iters){
		for(int i = 0; i < numVars; i++){
			double firstSum = 0.0;
			double secondSum = 0.0;
			for(int j = 0; j < i; j++){
				firstSum += a[i][j]*x[j];
			}
			for(int j = i; j < iters; j++){
				secondSum += a[i][j]*XO[j];
			}
			x[i] = (1.0/a[i][i]) * (-firstSum - secondSum + b[i]); // Finish this line
		}

		double vectorDiff[] = subtract_vectors(x, XO);
		double vecNorm = calc_norm(vectorDiff, 2); // euclidean norm
		if(vecNorm < tolerance){
			return x; // We solved it
		}
		k++;
		for(int i = 0; i < numVars; i++){
			XO[i] = x[i];
		}
	}

	return x;
}

double[] initialize_random_vector(int n, int min, int max){
	double output[n];
	for(int i = 0; i < n; i++){
		double randVal = ((rand() % 100) / 100.0) * (max - min) + min;
		output[i] = randVal;
	}

	return output;
}

// Subtract two vectors from each other
double[] subtract_vectors(double[] t1, double[] t2){
	int t1Size = t1.size();
	int t2Size = t2.size();
	if(t1Size != t2Size){
		return NULL;
	}

	double output[t1Size];
	for(int i = 0; i < t1Size; i++){
		int diff = t2[i] - t1[i];
		output[i] = diff;
	}

	return output;
}

// Calculate the p-norm between two vectors t1 and t2
double calc_norm(double vec[], double p){
	double sum = 0;
	for(i = 0; i < sizeof(vec)/sizeof(double); i++){
		sum += (p == 1) ? abs(vec[i]) : pow(vec[i], p);
	}
	double norm = pow(sum, 1.0/p);
	return norm;
}
