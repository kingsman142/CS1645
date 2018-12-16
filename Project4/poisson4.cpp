#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <iomanip>
#include "openacc.h"
#include <stdlib.h>
// NOTE: If you want to create the contour map, or .vtk file, uncomment this line when the vtk.c file is in your directory
#include "vtk.c"
using namespace std;

#define CONVERGENCE .000000000001

double calc_norm(double*, int, int, double, double);
double* initialize_arr(double*, int, int);
double* set_boundary_values(double*, int, int, double, double, double* l(double*, int, int, double), double* t(double*, int, int, double), double* r(double*, int, int, double), double* b(double*, int, int, double));
double* verif_left(double*, int, int, double);
double* verif_top(double*, int, int, double);
double* verif_right(double*, int, int, double);
double* verif_bottom(double*, int, int, double);
void print_solution(double*, int, int);

int main(int argc, char* argv[]){
	int x, y;

	cout<<"Enter x dimension: ";
	cin >> x;
	cout<<"Enter y dimension: ";
	cin >> y;

	double xConstraintMin = 0.0;
	double yConstraintMin = 0.0;
	double xConstraint = 2.0;
	double yConstraint = 1.0;
	double dx = (double) xConstraint / (x-1);
	double dy = (double) yConstraint / (y-1);
	int n = y*x;

	double* poissonMatrix = (double*) malloc((y*x) * sizeof(double));
	double* tempArr = (double*) malloc((y*x) * sizeof(double));

	initialize_arr(poissonMatrix, y, x);
	set_boundary_values(poissonMatrix, y, x, dx, dy, verif_left, verif_top, verif_right, verif_bottom);

	clock_t t1, t2; // create timers
	t1 = clock();

	double dy2 = dy*dy; // dy^2
	double dx2 = dx*dx; // dx^2

	#pragma acc data copy(poissonMatrix[0:n]) copy(tempArr[0:n])
	while(true){ // NOTE: Use iterations only for the problem solving example; otherwise, if you know the ground-truth value, use lines 106-108
		// Important for Jacobi algorithm: keep track of old values from previous iteration
		#pragma acc parallel loop gang
		for(int i = 0; i < y*x; i++) tempArr[i] = poissonMatrix[i];

		// Perform matrix value updates
		#pragma acc parallel loop gang
		for(int i = 1; i < y-1; i++){
			#pragma acc loop worker
			for(int j = 1; j < x-1; j++){
				double below = poissonMatrix[(i+1)*x + j];
				double above = poissonMatrix[(i-1)*x + j];
				double right = poissonMatrix[i*x + j+1];
				double left = poissonMatrix[i*x + j-1];
				double sourceTerm = (j*dx)*exp(i*dy);
				tempArr[i*x + j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * sourceTerm) / (2.0*dy2 + 2.0*dx2);
			}
		}

		double convergenceCheck = fabs(tempArr[x + 1] - poissonMatrix[x + 1]);

		// Calculate Norm
		#pragma acc parallel loop gang reduction(max:convergenceCheck)
		for(int i = 1; i < y-1; i++){
			#pragma acc loop worker
			for(int j = 1; j < x-1; j++){
				double diff = fabs(tempArr[i*x + j] - poissonMatrix[i*x + j]);
				convergenceCheck = max(convergenceCheck, diff);
			}
		}

		// "Move" to the next iteration
		double* temp = poissonMatrix;
		poissonMatrix = tempArr;
		tempArr = temp;

		if(convergenceCheck <= CONVERGENCE){
			break;
		}
	}

	// calculate final time difference
	t2 = clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;

	// NOTE: Uncomment this line to print out a square matrix of values containing the solution.
	// print_solution(poissonMatrix, y, x);

	// NOTE: The following three lines of code, when uncommented, allow the user to create a .vtk file of the output array to view the
	//		 contour map in any visualization software, such as Paraview.
	//		 Also required is uncommenting line 7, which is an include line, to gain access to VTK_out(...).
	// double* oneDPoisson = convertTo1D(poissonMatrix);
	VTK_out(x, y, &xConstraintMin, &xConstraint, &yConstraintMin, &yConstraint, poissonMatrix, 1);
	delete[] poissonMatrix;
	delete[] tempArr;

	// NOTE: The following two lines can only be used for the verification test, because we know the ground truth T(x, y) values.
	double error = calc_norm(poissonMatrix, y, x, dx, dy);
	cout<<"Error: "<<error<<endl;
	cout<<"Time to execute: "<<seconds<<endl;

	return 0;
}

// Print out any 2D matrix
void print_solution(double* a, int m, int n){
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			cout<<a[i*n + j]<<" ";
		}
		cout<<endl;
	}
}

// The user may pass in any four functions to set necessary boundary conditions of the matrix
double* set_boundary_values(double* a, int m, int n, double dx, double dy, double* left(double*, int, int, double), double* top(double*, int, int, double), double* right(double*, int, int, double), double* bottom(double*, int, int, double)){
	a = left(a, m, n, dy);
	a = right(a, m, n, dy);
	a = bottom(a, m, n, dx);
	a = top(a, m, n, dx);
	return a;
}

// Initialize the array to all 0s since C++ does not explicitly handle this
double* initialize_arr(double* output, int m, int n){
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[i*n + j] = 0.0;
		}
	}

	return output;
}

// Calculate the L-inf norm between the 'solved' matrix and its underlying true value
double calc_norm(double* a, int m, int n, double dx, double dy){
	double maxVal = fabs(a[n + 1] - dx*exp(dy)); // baseline max

	for(int i = 1; i < m-1; i++){
		for(int j = 1; j < n-1; j++){
			double sourceTerm = ((double) j*dx)*exp((double) i*dy); //NOTE: This is the verification source term (j*dx)*exp(i*dy); for the problem solving section, this value is 0.2
			double diff = fabs(a[i*n + j] - sourceTerm);
			maxVal = max(maxVal, diff);
		}
	}

	return maxVal;
}

// Verification problem, left boundary condition T(0, y) = 0
double* verif_left(double* a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[i*n] = 0;
	}

	return a;
}

// Verification problem, right boundary condition T(2, y) = 2e^y
double* verif_right(double* a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[i*n + n-1] = 2.0*exp(i*dy);
	}

	return a;
}

// Verification problem, bottom boundary condition T(x, 0) = x
double* verif_bottom(double* a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[i] = i*dx;
	}

	return a;
}

// Verification problem, top boundary condition T(x, 1) = xe
double* verif_top(double* a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[(m-1)*n + i] = (i*dx)*exp(1);
	}

	return a;
}
