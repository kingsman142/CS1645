#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

#define vvd vector< vector<double> >
#define CONVERGENCE .000000000001

double calc_norm(vvd, double, double);
void solve_jacobi();
vvd solve_gauss_seidel(vvd, int, int, double, double);
vvd initialize_arr(int, int);
vvd set_boundary_values(vvd, int, int, double, double, vvd l(vvd, int, double), vvd t(vvd, int, double), vvd r(vvd, int, int, double), vvd b(vvd, int, int, double));
vvd verif_left(vvd, int, double);
vvd verif_top(vvd, int, double);
vvd verif_right(vvd, int, int, double);
vvd verif_bottom(vvd, int, int, double);
void print_solution(vvd);

int main(int argc, char* argv[]){
	int x, y;

	cout<<"Enter x dimension: ";
	cin >> x;
	cout<<"Enter y dimension: ";
	cin >> y;

	double xConstraint = 2.0;
	double yConstraint = 1.0;
	double dx = (double) xConstraint / (x-1);
	double dy = (double) yConstraint / (y-1);

	vvd poissonMatrix = initialize_arr(y, x);
	poissonMatrix = set_boundary_values(poissonMatrix, y, x, dx, dy, verif_left, verif_top, verif_right, verif_bottom);
	poissonMatrix = solve_gauss_seidel(poissonMatrix, y, x, dx, dy);
	print_solution(poissonMatrix);
	double error = calc_norm(poissonMatrix, dx, dy);
	cout<<"Error: "<<error<<endl;
	return 0;
}

void print_solution(vvd a){
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a[0].size(); j++){
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
}

// TODO: Implement
void solve_jacobi(){
	return;
}

vvd solve_gauss_seidel(vvd a, int m, int n, double dx, double dy){
	double dx2 = pow(dx, 2);
	double dy2 = pow(dy, 2);
	while(true){
		vvd output(a);
		double convergenceCheck = a[1][1];
		for(int i = 1; i < m-1; i++){
			for(int j = 1; j < n-1; j++){
				double below = a[i+1][j];
				double here = a[i][j];
				double above = a[i-1][j];
				double right = a[i][j+1];
				double left = a[i][j-1];
				output[i][j] = (dy2*left + dy2*right + dx2*above + dx2*below - dx2*dy2*(j*dx)*exp(i*dy)) / (2*dy2 + 2*dx2);
			}
		}
		a = output;
		if(a[1][1] - convergenceCheck <= CONVERGENCE){
			break;
		}
	}

	return a;
}

vvd set_boundary_values(vvd a, int m, int n, double dx, double dy, vvd left(vvd, int, double), vvd top(vvd, int, double), vvd right(vvd, int, int, double), vvd bottom(vvd, int, int, double)){
	a = left(a, m, dy);
	a = right(a, m, n, dy);
	a = bottom(a, m, n, dx);
	a = top(a, n, dx);
	return a;
}

vvd initialize_arr(int m, int n){
	vvd output(m, vector<double>(n, 0));
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[i][j] = 0;
		}
	}

	return output;
}

// Calculate the p-norm between two vectors t1 and t2
double calc_norm(vvd a, double dx, double dy){
	double min = abs(a[1][1] - dx*exp(dy));
	cout<<"min: "<<min<<endl;
	for(int i = 1; i < a.size()-1; i++){
		for(int j = 1; j < a[0].size()-1; j++){
			double sourceTerm = (j*dx)*exp(i*dy);
			double diff = abs(a[i][j] - sourceTerm);
			if(diff < min) min = diff;
		}
	}

	return min;
}

vvd verif_left(vvd a, int m, double dy){
	for(int i = 0; i < m; i++){
		a[i][0] = 0;
	}

	return a;
}

vvd verif_right(vvd a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[m-i-1][n-1] = 2.0*exp(i*dy);
	}

	return a;
}

vvd verif_bottom(vvd a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[m-1][i] = i*dx;
	}

	return a;
}

vvd verif_top(vvd a, int n, double dx){
	for(int i = 0; i < n; i++){
		a[0][i] = (i*dx)*exp(1);
	}

	return a;
}
