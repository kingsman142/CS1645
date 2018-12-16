#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <iomanip>
// NOTE: If you want to create the contour map, or .vtk file, uncomment this line when the vtk.c file is in your directory
// #include "vtk.c"
using namespace std;

#define vvd vector< vector<double> >
#define CONVERGENCE .000000000001

double calc_norm(vvd, double, double);
vvd solve_jacobi(vvd, int, int, double, double);
vvd solve_gauss_seidel(vvd, int, int, double, double);
vvd initialize_arr(int, int);
vvd set_boundary_values(vvd, int, int, double, double, vvd l(vvd, int, double), vvd t(vvd, int, int, double), vvd r(vvd, int, int, double), vvd b(vvd, int, double));
vvd verif_left(vvd, int, double);
vvd verif_top(vvd, int, int, double);
vvd verif_right(vvd, int, int, double);
vvd verif_bottom(vvd, int, double);
vvd examp_left(vvd, int, double);
vvd examp_top(vvd, int, int, double);
vvd examp_right(vvd, int, int, double);
vvd examp_bottom(vvd, int, double);
void print_solution(vvd);
double* convertTo1D(vvd);
vvd subtract_vectors(vvd, vvd);
double calc_normal_norm(vvd);

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

	vvd poissonMatrix = initialize_arr(y, x);
	poissonMatrix = set_boundary_values(poissonMatrix, y, x, dx, dy, verif_left, verif_top, verif_right, verif_bottom);
	vvd poissonMatrixJacobi = solve_jacobi(poissonMatrix, y, x, dx, dy);
	vvd poissonMatrixGauss = solve_gauss_seidel(poissonMatrix, y, x, dx, dy);
	// NOTE: Uncomment these two lines to print out a square matrix of values containing the solutions for Jacobi and Gauss-Seidel
	print_solution(poissonMatrixJacobi);
	// print_solution(poissonMatrixGauss);

	// NOTE: The following three lines of code, when uncommented, allow the user to create a .vtk file of the output array to view the
	//		 contour map in any visualization software, such as Paraview.
	//		 Also required is uncommenting line 7, which is an include line, to gain access to VTK_out(...).
	// double* oneDPoisson = convertTo1D(poissonMatrixJacobi);
	// VTK_out(x, y, &xConstraintMin, &xConstraint, &yConstraintMin, &yConstraint, oneDPoisson, 1);
	// delete[] oneDPoisson;

	// NOTE: The following two lines can only be used for the verification test, because we know the ground truth T(x, y) values.
	// 	 Comment these lines when not solving the verification problem.
	double error = calc_norm(poissonMatrixGauss, dx, dy);
	cout<<endl<<"Error: "<<error<<endl;

	return 0;
}

// Convert the 2D vector<vector<double>> to a 1D array
double* convertTo1D(vvd a){
	double* output = new double[a.size()*a[0].size()];
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a[0].size(); j++){
			output[a[0].size() * i + j] = a[i][j];
		}
	}

	return output;
}

// Print out any 2D vector<vector<double>>
void print_solution(vvd a){
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a[0].size(); j++){
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
}

// Use the Jacobi algorithm to solve our Poisson equation
vvd solve_jacobi(vvd a, int m, int n, double dx, double dy){
	clock_t t1, t2; // create timers
	t1 = clock();
	double dx2 = pow(dx, 2); // dx^2
	double dy2 = pow(dy, 2); // dy^2

	int iters = 0;
	while(iters < 10){ // NOTE: Use iterations only for the problem solving example; otherwise, if you know the ground-truth value, use lines 106-108
		iters++;
		vvd output(a);
		for(int i = 1; i < m-1; i++){
			for(int j = 1; j < n-1; j++){
				double below = a[i+1][j];
				double here = a[i][j];
				double above = a[i-1][j];
				double right = a[i][j+1];
				double left = a[i][j-1];
				// NOTE: Verification source term: (j*dx)*exp(i*dy) -- if solving the verification proble, swap out 0.2 for the aforementioned source term
				// NOTE: If doing the problem solving portion, swap out the (j*dx)*exp(i*dy) for 0.2
				output[i][j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * (j*dx)*exp(i*dy) ) / (2.0*dy2 + 2.0*dx2);
			}
		}
		vvd iterDiff = subtract_vectors(output, a);
		a = output;

		double convergenceCheck = calc_normal_norm(iterDiff);
		if(convergenceCheck <= CONVERGENCE){
			break;
		}
	}

	cout<<"Iterations: "<<iters<<endl;
	// calculate final time difference
	t2 = clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	cout<<"Time to execute (Jacobi): "<<seconds<<endl<<endl;

	return a;
}

// Use the Gauss-Seidel algorithm to solve our Poisson equation
vvd solve_gauss_seidel(vvd a, int m, int n, double dx, double dy){
	clock_t t1, t2; // create timers
	t1 = clock();
	double dx2 = pow(dx, 2); // dx^2
	double dy2 = pow(dy, 2); // dy^2

	while(true){ // NOTE: Use iterations only for the problem solving example; otherwise, if you know the ground-truth value, use lines 106-108
		vvd output(a);
		for(int i = 1; i < m-1; i++){
			for(int j = 1; j < n-1; j++){
				double below = a[i+1][j];
				double here = a[i][j];
				double above = a[i-1][j];
				double right = a[i][j+1];
				double left = a[i][j-1];
				// NOTE: Verification source term: (j*dx)*exp(i*dy) -- if solving the verification proble, swap out 0.2 for the aforementioned source term
				// NOTE: If doing the problem solving portion, swap out the (j*dx)*exp(i*dy) for 0.2
				a[i][j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2 * (j*dx)*exp(i*dy) ) / (2.0*dy2 + 2.0*dx2);
			}
		}
		vvd iterDiff = subtract_vectors(output, a);

		if(calc_normal_norm(iterDiff) <= CONVERGENCE){
			break;
		}
	}

	// calculate final time difference
	t2 = clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	cout<<"Time to execute (Gauss-Seidel): "<<seconds<<endl<<endl;

	return a;
}

// Subtract any two vector<vector<double>> 2D vectors element-wise
vvd subtract_vectors(vvd a, vvd b){
	int m = a.size();
	int n = a[0].size();
	vvd output(m, vector<double>(n, 0));
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a[0].size(); j++){
			output[i][j] = a[i][j] - b[i][j];
		}
	}

	return output;
}

// The user may pass in any four functions to set necessary boundary conditions of the matrix
vvd set_boundary_values(vvd a, int m, int n, double dx, double dy, vvd left(vvd, int, double), vvd top(vvd, int, int, double), vvd right(vvd, int, int, double), vvd bottom(vvd, int, double)){
	a = left(a, m, dy);
	a = right(a, m, n, dy);
	a = bottom(a, n, dx);
	a = top(a, m, n, dx);
	return a;
}

// Initialize the array to all 0s since C++ does not explicitly handle this
vvd initialize_arr(int m, int n){
	vvd output(m, vector<double>(n, 0));
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[i][j] = 0;
		}
	}

	return output;
}

// calc_norm(...) calculates the L-inf norm, comparing a vector A to the exact solution of the Poisson equation
// This function, calc_normal_norm(...) calculates the L-inf norm on a vector a (typically done after two vectors have been subtracted on an iteration-by-iteration basis)
double calc_normal_norm(vvd a){
	double max = a[1][1];
	for(int i = 1; i < a.size()-1; i++){
		for(int j = 1; j < a[0].size()-1; j++){
			double diff = fabs(a[i][j]);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Calculate the L-inf norm between the 'solved' matrix and its underlying true value
double calc_norm(vvd a, double dx, double dy){
	double max = fabs(a[1][1] - dx*exp(dy)); // baseline max
	int maxX = 0;
	int maxY = 0;
	double maxSource = 0;
	for(int i = 1; i < a.size()-1; i++){
		for(int j = 1; j < a[0].size()-1; j++){
			double sourceTerm = ((double) j*dx)*exp((double) i*dy); //NOTE: This is the verification source term (j*dx)*exp(i*dy); for the problem solving section, this value is 0.2
			double diff = fabs(a[i][j] - sourceTerm);
			if(diff > max) max = diff;
		}
	}

	return max;
}

// Verification problem, left boundary condition T(0, y) = 0
vvd verif_left(vvd a, int m, double dy){
	for(int i = 0; i < m; i++){
		a[i][0] = 0;
	}

	return a;
}

// Verification problem, right boundary condition T(2, y) = 2e^y
vvd verif_right(vvd a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[i][n-1] = 2.0*exp(i*dy);
	}

	return a;
}

// Verification problem, bottom boundary condition T(x, 0) = x
vvd verif_bottom(vvd a, int n, double dx){
	for(int i = 0; i < n; i++){
		a[0][i] = i*dx;
	}

	return a;
}

// Verification problem, top boundary condition T(x, 1) = xe
vvd verif_top(vvd a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[m-1][i] = (i*dx)*exp(1);
	}

	return a;
}

// Problem solving example, left boundary condition T(0, y) = 0
vvd examp_left(vvd a, int m, double dy){
	for(int i = 0; i < m; i++){
		a[i][0] = 0;
	}

	return a;
}

// Problem solving example, right boundary condition T(1, y) = 0
vvd examp_right(vvd a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[i][n-1] = 0;
	}

	return a;
}

// Problem solving example, bottom boundary condition T(x, 0) = 0
vvd examp_bottom(vvd a, int n, double dx){
	for(int i = 0; i < n; i++){
		a[0][i] = 0;
	}

	return a;
}

// Problem solving example, top boundary condition T(x, 1) = 0
vvd examp_top(vvd a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[m-1][i] = 0;
	}

	return a;
}
