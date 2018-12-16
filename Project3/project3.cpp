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
#define ROOT 0

double* initialize_arr(int, int);
void set_boundary_values(double*, int, int, int, int, double, double, int, int, int, void l(double*, int, double, int, int, int, int), void r(double*, int, double, int, int, int, int), void t(double*, int, int, double, int, int, int, int), void b(double*, double, int, int, int, int));
double* solve_poisson_equation(double*, double, double, int, int, int, int, int, int, int);
void verif_left(double*, int, double, int, int, int, int);
void verif_right(double*, int, double, int, int, int, int);
void verif_top(double*, int, int, double, int, int, int, int);
void verif_bottom(double*, double, int, int, int, int);
void subtract_vectors(double*, double*, double*, int, int, int, int, int); // subtract two vectors
void print_matrix(double*, int, int);
double calc_norm(double*, double, double, int, int); // used to compute error once final vector/matrix is calculated
double calc_normal_norm(double*, int, int, int, int, int); // used to calculate the norm of the difference between two vectors

int main(int argc, char* argv[]){
	int x, y, npX, npY;

	cout<<"Enter grid x dimension: ";
	cin >> x;
	cout<<"Enter grid y dimension: ";
	cin >> y;
	cout<<"Processors in x dimension: ";
	cin >> npX;
	cout<<"Processors in y dimension: ";
	cin >> npY;

	if(x % npX != 0 || y % npY != 0){ // can the grid be evenly divided by the processor layout specified by the user?
		cout<<"Mismatch in the number of processors!"<<endl;
		return 1;
	} else if((x * y * sizeof(double))/1024 > HEAP_SIZE){ // can the grid even be stored on the heap?
		cout<<"Grid size too large!"<<endl;
		return 1;
	}

	int np = npX * npY;
	omp_set_num_threads(np);
	double* poissonMatrix = initialize_arr(y, x);

	double xConstraintMin = 0.0;
	double yConstraintMin = 0.0;
	double xConstraint = 2.0;
	double yConstraint = 1.0;
	double dx = (double) xConstraint / (x-1);
	double dy = (double) yConstraint / (y-1);
	double subMatrixX = x / npX;
	double subMatrixY = y / npY;

	#pragma omp parallel shared(poissonMatrix) num_threads(np)
	{
		int rank = omp_get_thread_num();
		set_boundary_values(poissonMatrix, subMatrixY, subMatrixX, y, x, dx, dy, npX, np, rank, verif_left, verif_right, verif_top, verif_bottom);
		#pragma omp barrier
		poissonMatrix = solve_poisson_equation(poissonMatrix, dx, dy, subMatrixY, subMatrixX, rank, npX, npY, y, x);
	}
	// NOTE: Uncomment these two lines to print out a square matrix of values containing the solution.
	// print_matrix(poissonMatrix, y, x);
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

// Use the Jacobi algorithm to solve our Poisson equation
double* solve_poisson_equation(double* a, double dx, double dy, int subMatrixM, int subMatrixN, int rank, int npX, int npY, int m, int n){
	double dx2 = pow(dx, 2); // dx^2
	double dy2 = pow(dy, 2); // dy^2
	int np = npX * npY; // number of threads
	clock_t t1, t2; // create timers
	t1 = clock();

	// respective indices of this thread's chunk in the larger matrix (kind of creates a 'rectangle' of coordinates around the chunk)
	int leftBoundaryStartIndex = (rank % npX) * subMatrixN;
	int rightBoundaryEndIndex = leftBoundaryStartIndex + subMatrixN;
	int topBoundaryStartIndex = (rank / npX) * subMatrixM;
	int bottomBoundaryEndIndex = topBoundaryStartIndex + subMatrixM;

	// if a chunk's boundary starts on the matrix's boundary, make sure to make the chunk a little bit smaller by moving the indices inward
	if(leftBoundaryStartIndex == 0) leftBoundaryStartIndex = 1;
	if(rightBoundaryEndIndex == n) rightBoundaryEndIndex = n-1;
	if(topBoundaryStartIndex == 0) topBoundaryStartIndex = 1;
	if(bottomBoundaryEndIndex == m) bottomBoundaryEndIndex = m-1;

	// boundary values, but include the boundary values as well this time
	int trueLeftIndex = leftBoundaryStartIndex == 1 ? 0 : leftBoundaryStartIndex;
	int trueRightIndex = rightBoundaryEndIndex == n-1 ? n : rightBoundaryEndIndex;
	int trueTopIndex = topBoundaryStartIndex == 1 ? 0 : topBoundaryStartIndex;
	int trueBottomIndex = bottomBoundaryEndIndex == m-1 ? m : bottomBoundaryEndIndex;

	double* output; // just a placeholder for memory for the weights since we're performing the Jacobi algorithm
	double* iterDiff; // current iteration's weights minus the previous iteration's weights
	int iters = 0;
	double maxConvergence; // max convergence across all threads' chunks; this essentially helps compute l-infinity norm when comparing weight iterations
	bool breakFlag = false; // has our matrix converged? should we break? OpenMP doesn't allow us to place an explicit 'break' statement

	#pragma omp parallel shared(a, iters, output, iterDiff, maxConvergence)
	{
		while(!breakFlag){
			#pragma omp single
			{
				iters++; // update the number of iterations
				maxConvergence = 0.00; // reset the maximum convergence value across all threads for this iteration
			}
			#pragma omp barrier

			#pragma omp single
			output = new double[m*n]; // only allocate one new output buffer

			// copy all the values from the array, a, to the output placeholder matrix
			for(int i = trueTopIndex; i < trueBottomIndex; i++){
				for(int j = trueLeftIndex; j < trueRightIndex; j++){
					output[i*n + j] = a[i*n + j];
				}
			}
			#pragma omp barrier

			// perform the actual weight updates
			for(int i = topBoundaryStartIndex; i < bottomBoundaryEndIndex; i++){
				for(int j = leftBoundaryStartIndex; j < rightBoundaryEndIndex; j++){
					double below = a[(i+1)*n + j];
					double here = a[i*n + j];
					double above = a[(i-1)*n + j];
					double right = a[i*n + j+1];
					double left = a[i*n + j-1];
					// NOTE: Verification source term: (j*dx)*exp(i*dy)
					double sourceTerm = (j*dx)*exp(i*dy);
					output[i*n + j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * sourceTerm) / (2.0*dy2 + 2.0*dx2);
				}
			}

			// create a new array to subtract the current iteration's weights from the previous iteration's weights
			#pragma omp single
			iterDiff = new double[m*n];

			// subtract the current iteration's weights from the previous iteration's weights and store it in iterDiff
			subtract_vectors(output, a, iterDiff, n, trueLeftIndex, trueRightIndex, trueTopIndex, trueBottomIndex); // difference between iterations so we can take the norm next and check convergence

			// copy all values from the output array (newly updated matrix) into our official 'answer' matrix, which is 'a'
			for(int i = topBoundaryStartIndex; i < bottomBoundaryEndIndex; i++){
				for(int j = leftBoundaryStartIndex; j < rightBoundaryEndIndex; j++){
					a[i*n + j] = output[i*n + j];
				}
			}
			#pragma omp barrier // make sure all threads are here before we check for convergence

			// find the norm for the chunk on each processor, then take the MAX of those norms and check if the entire grid has converged
			double convergenceCheck = calc_normal_norm(iterDiff, n, leftBoundaryStartIndex, rightBoundaryEndIndex, topBoundaryStartIndex, bottomBoundaryEndIndex);
			#pragma omp critical
			{
				if(convergenceCheck > maxConvergence) maxConvergence = convergenceCheck; // basically performs a MAX reducation across all threads
			}
			#pragma omp barrier // make sure every thread was given the opportunity to update the MAX value

			if(maxConvergence <= CONVERGENCE){ // based on what has been assessed from all other processors, has the entire grid converged yet?
				breakFlag = true;
			}

			// free the malloc'd memory to prevent memory leaks
			#pragma omp single
			delete iterDiff;
		}
	}

	delete output;

	if(rank == ROOT) cout<<"Iterations: "<<iters<<endl;

	// calculate final time difference
	t2 = clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	#pragma omp parallel reduction(+ : seconds)
	{
		// Just need something here as a placeholder; without it, OpenMP will spit out an error about avgTimeDiff not being declared
	}

	double avgTimeDiff = seconds / (double) np; // average the time taken by all processors
	if(rank == ROOT) cout<<"Avg. time to execute: "<<avgTimeDiff<<" seconds"<<endl<<endl;

	return a;
}

// Subtract any two vector<vector<double>> 2D vectors element-wise
void subtract_vectors(double* a, double* b, double* diff, int n, int leftIndex, int rightIndex, int topIndex, int bottomIndex){
	for(int i = topIndex; i < bottomIndex; i++){
		for(int j = leftIndex; j < rightIndex; j++){
			diff[n*i + j] = a[i*n + j] - b[i*n + j];
		}
	}
}

// The user may pass in any four functions to set necessary boundary conditions of the matrix
void set_boundary_values(double* a, int subMatrixM, int subMatrixN, int m, int n, double dx, double dy, int npX, int np, int rank, void left(double*, int, double, int, int, int, int), void right(double*, int, double, int, int, int, int), void top(double*, int, int, double, int, int, int, int), void bottom(double*, double, int, int, int, int)){
	int leftIndex = (rank % npX) * subMatrixN; // left index of the larger matrix for this specific thread's 'chunk'
	int rightIndex = leftIndex + subMatrixN; // right index of the larger matrix for this specific thread's 'chunk'
	int topIndex = (rank / npX) * subMatrixM; // top index of the larger matrix for this specific thread's 'chunk'
	int bottomIndex = topIndex + subMatrixM; // bottom index of the larger matrix for this specific thread's 'chunk'

	if(rank % npX == 0) left(a, n, dy, rank, npX, topIndex, bottomIndex); // this processor has left boundary conditions to set
	if((rank+1) % npX == 0) right(a, n, dy, rank, npX, topIndex, bottomIndex); // this processor has right boundary conditions to set
	if(rank < npX) top(a, m, n, dx, rank, npX, leftIndex, rightIndex); // this processor has top boundary conditions to set
	if((np - rank) <= npX) bottom(a, dx, rank, npX, leftIndex, rightIndex); // this processor has bottom boundary conditions to set
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
double calc_normal_norm(double* a, int n, int leftIndex, int rightIndex, int topIndex, int bottomIndex){
	double max = a[n*topIndex + leftIndex]; // initialize the max to some arbitrary value, such as the top left corner (outside of the boundaries)
	for(int i = topIndex; i < bottomIndex; i++){
		for(int j = leftIndex; j < rightIndex; j++){
			double diff = fabs(a[i*n + j]);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Calculate the L-inf norm between the 'solved' matrix and its underlying true value
double calc_norm(double* a, double dx, double dy, int m, int n){
	double max = fabs(a[n + 1] - (dx)*exp(dy)); // baseline max

	for(int i = 1; i < m-1; i++){
		for(int j = 1; j < n-1; j++){
			//NOTE: This is the verification source term (j*dx)*exp(i*dy)
			double sourceTerm = ((double) j*dx)*exp((double) i*dy);
			double diff = fabs(a[i*n + j] - sourceTerm);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Verification problem, left boundary condition T(0, y) = 0
void verif_left(double* a, int n, double dy, int rank, int npX, int topIndex, int bottomIndex){
	for(int i = topIndex; i < bottomIndex; i++){
		a[i*n] = 0;
	}
}

// Verification problem, right boundary condition T(2, y) = 2e^y
void verif_right(double* a, int n, double dy, int rank, int npX, int topIndex, int bottomIndex){
	for(int i = topIndex; i < bottomIndex; i++){
		a[n*(i+1) - 1] = 2.0*exp(i*dy);
	}
}

// Verification problem, top boundary condition T(x, 1) = xe
void verif_top(double* a, int m, int n, double dx, int rank, int npX, int leftIndex, int rightIndex){
	for(int i = leftIndex; i < rightIndex; i++){
		a[n*(m-1) + i] = (i*dx)*exp(1);
	}
}

// Verification problem, bottom boundary condition T(x, 0) = x
void verif_bottom(double* a, double dx, int rank, int npX, int leftIndex, int rightIndex){
	for(int i = leftIndex; i < rightIndex; i++){
		a[i] = i*dx;
	}
}
