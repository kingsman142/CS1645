#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <mpi.h>
// NOTE: If you want to create the contour map, or .vtk file, uncomment this line when the vtk.c file is in your directory
#include "vtk.c"
using namespace std;

#define CONVERGENCE .000000000001
#define HEAP_SIZE 8192

double calc_norm(double*, double, double, int, int);
void solve_jacobi();
double* solve_poisson_equation(double*, double, double, int, int, int, int, int);
double* initialize_arr(int, int);
double* set_boundary_values(double*, int, int, double, double, int, int, int, double* l(double*, int, int, double, int, int), double* t(double*, int, int, double, int, int), double* r(double*, int, int, double, int, int), double* b(double*, int, int, double, int, int));
double* verif_left(double*, int, int, double, int, int);
double* verif_top(double*, int, int, double, int, int);
double* verif_right(double*, int, int, double, int, int);
double* verif_bottom(double*, int, int, double, int, int);
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
	#pragma parallel shared(poissonMatrix){
		int rank = omp_get_thread_num();

		double xConstraintMin = 0.0;
		double yConstraintMin = 0.0;
		double xConstraint = 2.0;
		double yConstraint = 1.0;
		double dx = (double) xConstraint / (x-1);
		double dy = (double) yConstraint / (y-1);
		double subMatrixX = x / npX;
		double subMatrixY = y / npY;

		double chunkXMin = (rank % npX) * (subMatrixX) * dx; // minimum X value in the grid for this processor
		double chunkXMax = chunkXMin + (subMatrixX*dx); // maximum X value in the grid for this processor
		double chunkYMin = (rank / npX) * (subMatrixY) * dy; // minimum Y value in the grid for this processor
		double chunkYMax = chunkYMin + (subMatrixY*dy); // maximum Y value in the grid for this processor

		poissonMatrix = set_boundary_values(poissonMatrix, subMatrixY, subMatrixX, dx, dy, npX, np, rank, verif_left, verif_top, verif_right, verif_bottom);
		poissonMatrix = solve_poisson_equation(poissonMatrix, dx, dy, subMatrixY, subMatrixX, rank, npX, npY);
	}
	// NOTE: Uncomment these three line to print out a square matrix of values containing the solution.
	// cout<<"Processor "<<rank<<"'s matrix:"<<endl;
	// print_matrix(poissonMatrix, subMatrixY, subMatrixN);
	// cout<<endl;

	// NOTE: The following line of code, when uncommented, allow the user to create a .vtk file of the output array to view the
	//		 contour map in any visualization software, such as Paraview.
	//		 Also required is uncommenting line 7, which is an include line, to gain access to VTK_out(...).
	VTK_out(x, y, &xConstraintMin, &xConstraint, &yConstraintMin, &yConstraint, poissonMatrix, 1);

	double error = calc_norm(poissonMatrix, dx, dy, y, x);
	double minError = -1.0; // placeholder value, just in case the message passing fails
	cout<<"Error: "<<minError<<endl;

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
double* solve_poisson_equation(double* a, double dx, double dy, int subMatrixM, int subMatrixN, int rank, int npX, int npY){
	double dx2 = pow(dx, 2); // dx^2
	double dy2 = pow(dy, 2); // dy^2
	double chunkStartX = (rank % npX) * subMatrixN * dx; // the starting X value of this processor's chunk
	double chunkStartY = (rank / npX) * subMatrixM * dy; // the starting Y value of this processor's chunk
	int np = npX * npY; // number of processors
	//double startTime = MPI_Wtime();

	// some processors have to watch out for boundary conditions, so some update values on index 0, n-1, and m-1, but others have boundary conditions there
	int leftBoundaryStartIndex = (rank % npX) == 0 ? 1 : 0; // on the left side of the grid? start at index 1; otherwise, start at 0
	int rightBoundaryEndIndex = ((rank+1) % npX) == 0 ? subMatrixN-1 : subMatrixN; // on the right side of the grid? go up to index subMatrixN-1; otherwise, go to subMatrixN
	int topBoundaryStartIndex = (rank < npX) ? 1 : 0; // on the top side of the grid? start at index 1; otherwise, start at 0
	int bottomBoundaryEndIndex = ((np - rank) <= npX) ? subMatrixM-1 : subMatrixM; // on the bottom of the grid? go up to index subMatrixM-1; otherwise, go to subMatrixM

	int iters = 0;
	while(true){
		iters++; // update the number of iterations
		double* output = copy_arr(a, subMatrixM, subMatrixN); // store a copy of the old values of the array A, so we can find the convergence between iterations and know when to stop

		for(int i = topBoundaryStartIndex; i < bottomBoundaryEndIndex; i++){
			for(int j = leftBoundaryStartIndex; j < rightBoundaryEndIndex; j++){
				double below = (i < subMatrixM-1) ? output[subMatrixN*(i+1) + j];
				double here = output[i*subMatrixN + j];
				double above = (i > 0) ? output[subMatrixN*(i-1) + j] : top[j];
				double right = (j < subMatrixN-1) ? output[i*subMatrixN + j+1];
				double left = (j > 0) ? output[i*subMatrixN + j-1];
				// NOTE: Verification source term: (j*dx)*exp(i*dy)
				double sourceTerm = (chunkStartX + j*dx)*exp(chunkStartY + i*dy); // ensure chunkStartX and chunkStartY add the necessary offsets
				a[subMatrixN*i + j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * sourceTerm) / (2.0*dy2 + 2.0*dx2);
			}
		}
		#pragma omp barrier

		double* iterDiff = subtract_vectors(output, a, subMatrixM, subMatrixN); // difference between iterations so we can take the norm next and check convergence

		// find the norm for the chunk on each processor, then take the MAX of those norms and check if the entire grid has converged
		double processorConvergenceCheck = calc_normal_norm(iterDiff, subMatrixM, subMatrixN);
		double convergenceCheck = 1.0; // placeholder value, just in case the parallel processing screws up
		if(convergenceCheck <= CONVERGENCE){ // based on what has been assessed from all other processors, has the entire grid converged yet?
			break;
		}

		// free the malloc'd memory to prevent memory leaks
		delete output;
		delete iterDiff;
	}

	if(rank == ROOT) cout<<"Iterations: "<<iters<<endl;

	// calculate final time difference
	double endTime = MPI_Wtime();
	double timeDiff = endTime - startTime;
	double timeDiffsSum;
	MPI_Reduce(&timeDiff, &timeDiffsSum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD); // grab the sum of all processor times so we can average it later
	double avgTimeDiff = timeDiffsSum / (double) np; // average the time taken by all processors
	if(rank == ROOT) cout<<"Avg. time to execute: "<<avgTimeDiff<<" seconds"<<endl<<endl;

	return a;
}

// Make a deep copy oe an array A and return the new array
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
double* set_boundary_values(double* a, int m, int n, double dx, double dy, int npX, int np, int rank, double* left(double*, int, int, double, int, int), double* top(double*, int, int, double, int, int), double* right(double*, int, int, double, int, int), double* bottom(double*, int, int, double, int, int)){
	if(rank % npX == 0) a = left(a, m, n, dy, rank, npX); // this processor has left boundary conditions to set
	if((rank+1) % npX == 0) a = right(a, m, n, dy, rank, npX); // this processor has right boundary conditions to set
	if((np - rank) <= npX) a = top(a, m, n, dx, rank, npX); // this processor has bottom boundary conditions to set
	if(rank < npX) a = bottom(a, m, n, dx, rank, npX); // this processor has top boundary conditions to set
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

// Get this processor's border values from parallel processes; for example, if this processor has a neighbor to the right and bottom, grab and send border values for processors to the right and bottom
void get_border_values(double* a, double* top, double* bottom, double* left, double* right, int subMatrixM, int subMatrixN, int npX, int np, int rank){
	MPI_Request topRecvReq, bottomRecvReq, leftRecvReq, rightRecvReq, topSendReq, bottomSendReq, leftSendReq, rightSendReq;

	if(rank >= npX){ // this processor's chunk is not on the top row of the grid, so receive its respective top border values and send its top row to the top process for that processor's bottom border values
		MPI_Irecv(top, subMatrixN, MPI_DOUBLE, rank - npX, 0, MPI_COMM_WORLD, &topRecvReq);
		MPI_Isend(a, subMatrixN, MPI_DOUBLE, rank - npX, 0, MPI_COMM_WORLD, &bottomSendReq);
		MPI_Wait(&topRecvReq, MPI_STATUS_IGNORE);
		MPI_Wait(&bottomSendReq, MPI_STATUS_IGNORE);
	}
	if((np - rank) > npX){ // this processor's chunk is not on the bottom row of the grid, so receive its respective bottom border values and send its bottom row to the bottom process for that processor's top border values
		MPI_Irecv(bottom, subMatrixN, MPI_DOUBLE, rank + npX, 0, MPI_COMM_WORLD, &bottomRecvReq);
		MPI_Isend(a + subMatrixN*(subMatrixM-1), subMatrixN, MPI_DOUBLE, rank + npX, 0, MPI_COMM_WORLD, &topSendReq);
		MPI_Wait(&bottomRecvReq, MPI_STATUS_IGNORE);
		MPI_Wait(&topSendReq, MPI_STATUS_IGNORE);
	}
	if((rank+1) % npX != 0){ // this processor's chunk is not on the right column of the grid, so receive its respective right border values and send its right column to the right process for that processor's left border values
		MPI_Irecv(right, subMatrixM, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &rightRecvReq);
		double leftSend[subMatrixM]; // grab this chunk's left-most column
		for(int i = 0; i < subMatrixM; i++) leftSend[i] = a[subMatrixN*(i+1) - 1];
		MPI_Isend(leftSend, subMatrixM, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &leftSendReq);
		MPI_Wait(&rightRecvReq, MPI_STATUS_IGNORE);
		MPI_Wait(&leftSendReq, MPI_STATUS_IGNORE);
	}
	if(rank % npX != 0){ // this processor's chunk is not on the left column of the grid, so receive its respective left border values and send its left column to the left process for that processor's right border values
		MPI_Irecv(left, subMatrixM, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &leftRecvReq);
		double rightSend[subMatrixM]; // grab this chunk's right-most column
		for(int i = 0; i < subMatrixM; i++) rightSend[i] = a[subMatrixN*i];
		MPI_Isend(rightSend, subMatrixM, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &rightSendReq);
		MPI_Wait(&leftRecvReq, MPI_STATUS_IGNORE);
		MPI_Wait(&rightSendReq, MPI_STATUS_IGNORE);
	}
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
double calc_norm(double* a, double dx, double dy, int m, int n, int rank, int npX){
	double chunkStartX = (rank % npX) * n * dx; // the starting X value of this processor's chunk
	double chunkStartY = (rank / npX) * m * dy; // the starting Y value of this processor's chunk
	double max = fabs(a[n + 1] - (chunkStartX + dx)*exp(chunkStartY + dy)); // baseline max

	for(int i = 1; i < m-1; i++){
		for(int j = 1; j < n-1; j++){
			//NOTE: This is the verification source term (j*dx)*exp(i*dy)
			double sourceTerm = ((double) chunkStartX + j*dx)*exp((double) chunkStartY + i*dy);
			double diff = fabs(a[i*n + j] - sourceTerm);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Verification problem, left boundary condition T(0, y) = 0
double* verif_left(double* a, int m, int n, double dy, int rank, int npX){
	for(int i = 0; i < m; i++){
		a[i*n] = 0;
	}

	return a;
}

// Verification problem, right boundary condition T(2, y) = 2e^y
double* verif_right(double* a, int m, int n, double dy, int rank, int npX){
	int chunkNumY = rank / npX;
	double startY = chunkNumY * m * dy;
	for(int i = 0; i < m; i++){
		a[n*(i+1) - 1] = 2.0*exp(startY + i*dy);
	}

	return a;
}

// Verification problem, bottom boundary condition T(x, 0) = x
double* verif_bottom(double* a, int m, int n, double dx, int rank, int npX){
	int chunkNumX = rank % npX;
	double startX = chunkNumX * n * dx;
	for(int i = 0; i < n; i++){
		a[i] = startX + i*dx;
	}

	return a;
}

// Verification problem, top boundary condition T(x, 1) = xe
double* verif_top(double* a, int m, int n, double dx, int rank, int npX){
	int chunkNumX = rank % npX;
	double startX = chunkNumX * n * dx;
	for(int i = 0; i < n; i++){
		a[n*(m-1) + i] = (startX + i*dx)*exp(1);
	}

	return a;
}
