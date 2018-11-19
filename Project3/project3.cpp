#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
// NOTE: If you want to create the contour map, or .vtk file, uncomment this line when the vtk.c file is in your directory
#include "vtk.c"
using namespace std;

#define CONVERGENCE .000000000001
#define HEAP_SIZE 8192

double calc_norm(double*, double, double, int, int);
double* solve_poisson_equation(double*, double, double, int, int, int, int);
double* initialize_arr(int, int);
double* set_boundary_values(double*, int, int, double, double, double* l(double*, int, int, double), double* t(double*, int, int, double), double* r(double*, int, int, double), double* b(double*, int, int, double));
double* verif_left(double*, int, int, double);
double* verif_top(double*, int, int, double);
double* verif_right(double*, int, int, double);
double* verif_bottom(double*, int, int, double);
void print_matrix(double*, int, int);
double* subtract_vectors(double*, double*, int, int);
double calc_normal_norm(double*, int, int);
void get_border_values(double*, double*, double*, double*, double*, int, int, int, int, int);
double* copy_arr(double*, int, int);

int main(int argc, char* argv[]){
	int x, y, npX, npY;

	cout<<"Enter grid x dimension: ";
	cin >> x;
	cout<<"Enter grid y dimension: ";
	cin >> y;
	cout<<"Threads in x dimension: ";
	cin >> npX;
	cout<<"Threads in y dimension: ";
	cin >> npY;

	if(x % npX != 0 || y % npY != 0){ // can the grid be evenly divided by the processor layout specified by the user?
		cout<<"Mismatch in the number of processors!"<<endl;
		return 1;
	} else if((x * y * sizeof(double))/1024 > HEAP_SIZE){ // can the grid even be stored on the heap?
		cout<<"Grid size too large!"<<endl;
		return 1;
	}

	double* poissonMatrix = initialize_arr(y, x);
	poissonMatrix = set_boundary_values(poissonMatrix, y, x, dx, dy, verif_left, verif_top, verif_right, verif_bottom);
	omp_set_num_threads(npX*npY);
	poissonMatrix = solve_poisson_equation(poissonMatrix, dx, dy, y, x, npX, npY);
	// NOTE: Uncomment these three line to print out a square matrix of values containing the solution.
	// cout<<"Processor "<<rank<<"'s matrix:"<<endl;
	// print_matrix(poissonMatrix, subMatrixY, subMatrixN);
	// cout<<endl;

	// NOTE: The following line of code, when uncommented, allow the user to create a .vtk file of the output array to view the
	//		 contour map in any visualization software, such as Paraview.
	//		 Also required is uncommenting line 7, which is an include line, to gain access to VTK_out(...).
	VTK_out(x, y, &xConstraintMin, &xConstraint, &yConstraintMin, &yConstraint, poissonMatrix, 1);

	double error = calc_norm(poissonMatrix, dx, dy, y, x);
	cout<<"Error: "<<error<<endl;

	delete[] poissonMatrix; // free the malloc'd memory to prevent memory leaks

	return 0;
}

// Print any given matrix given its dimensions m and n
void print_matrix(double* a, int m, int n){
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			cout<<a[i*n + j]<<" ";
		}
		cout<<endl;
	}
}

// Use a Gauss-Seidel and Jacobi hybrid multi-processor algorithm to solve our Poisson equation
double* solve_poisson_equation(double* a, double dx, double dy, int m, int n, int npX, int npY){
	double dx2 = pow(dx, 2); // dx^2
	double dy2 = pow(dy, 2); // dy^2
	int np = npX * npY; // number of processors
	//double startTime = MPI_Wtime();

	int iters = 0;
	while(true){
		iters++; // update the number of iterations
		double* output = copy_arr(a, m, n); // store a copy of the old values of the array A, so we can find the convergence between iterations and know when to stop

		#pragma omp parallel for num_threads(npY)
			for(int i = 1; i < m-1; i++){
				#pragma omp parallel for num_threads(npX)
					for(int j = 1; j < n-1; j++){
						double below = output[n*(i+1) + j];
						double here = output[i*n + j];
						double above = output[n*(i-1) + j];
						double right = output[i*n + j+1];
						double left = output[i*n + j-1];
						// NOTE: Verification source term: (j*dx)*exp(i*dy)
						double sourceTerm = (j*dx)*exp(i*dy); // ensure chunkStartX and chunkStartY add the necessary offsets
						a[n*i + j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * sourceTerm) / (2.0*dy2 + 2.0*dx2);
					}
			}

		double* iterDiff = subtract_vectors(output, a, m, n); // difference between iterations so we can take the norm next and check convergence

		// find the norm for the chunk on each processor, then take the MAX of those norms and check if the entire grid has converged
		double convergenceCheck = calc_normal_norm(iterDiff, m, n);
		if(convergenceCheck <= CONVERGENCE){ // based on what has been assessed from all other processors, has the entire grid converged yet?
			break;
		}

		// free the malloc'd memory to prevent memory leaks
		delete output;
		delete iterDiff;
	}

	cout<<"Iterations: "<<iters<<endl;

	// calculate final time difference
	/*double endTime = MPI_Wtime();
	double timeDiff = endTime - startTime;
	double timeDiffsSum = 5.0;
	//MPI_Reduce(&timeDiff, &timeDiffsSum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD); // grab the sum of all processor times so we can average it later
	double avgTimeDiff = timeDiffsSum / (double) np; // average the time taken by all processors
	cout<<"Avg. time to execute: "<<avgTimeDiff<<" seconds"<<endl<<endl;*/

	return a;
}

// Make a deep copy of an array A and return the new array
double* copy_arr(double* a, int m, int n){
	double* output = new double[m*n];
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[n*i + j] = a[n*i + j];
		}
	}

	return output;
}

// Subtract any two vector<vector<double>> 2D vectors element-wise
double* subtract_vectors(double* a, double* b, int m, int n){
	double* output = new double[m*n];
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[n*i + j] = a[i*n + j] - b[i*n + j];
		}
	}

	return output;
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
double* initialize_arr(int m, int n){
	double* output = new double[m*n];
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[n*i + j] = 0;
		}
	}

	return output;
}

// calc_norm(...) calculates the L-inf norm, comparing a vector A to the exact solution of the Poisson equation
// This function, calc_normal_norm(...) calculates the L-inf norm on a vector a (typically done after two vectors have been subtracted on an iteration-by-iteration basis)
double calc_normal_norm(double* a, int m, int n){
	double max = a[n*1 + 1];
	for(int i = 1; i < m-1; i++){
		for(int j = 1; j < n-1; j++){
			double diff = fabs(a[i*n + j]);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Calculate the L-inf norm between the 'solved' matrix and its underlying true value
double calc_norm(double* a, double dx, double dy, int m, int n){
	double max = fabs(a[n+1] - dx*exp(dy)); // baseline max
	int maxX = 0;
	int maxY = 0;
	double maxSource = 0;
	for(int i = 1; i < m-1; i++){
		for(int j = 1; j < n-1; j++){
			double sourceTerm = ((double) j*dx)*exp((double) i*dy); //NOTE: This is the verification source term (j*dx)*exp(i*dy); for the problem solving section, this value is 0.2
			double diff = fabs(a[n*i + j] - sourceTerm);
			if(diff > max) max = diff;
		}
	}

	return max;
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
		a[n*(i+1) - 1] = 2.0*exp(i*dy);
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
		a[n*(m-1) + i] = (i*dx)*exp(1);
	}

	return a;
}
